//
// Created by fritsche on 29/12/2021.
//

#pragma once

#include <compact_map/compact_map.h>
#include <Benchmark.h>
#include <kmer_processing/SimpleKmerIterator.h>
#include <kmer_processing/SimpleKmerPutter.h>
#include <omp.h>
#include <bloom_filter.h>
#include "OptionsContainer.h"
#include "classification.h"


namespace BuildMode {
    const std::string DIVIDER = std::string(60, '-');

    /* #######################################################################
     *  Save Bucket Sizes
     * #################################################################### */
    static void SaveBucketSizes(std::string path, size_t* buckets, size_t bucket_size) {
        ofstream ofs(path, ios::binary);
        ofs.write((char *) &bucket_size, sizeof(bucket_size));
        ofs.write((char *) buckets, sizeof(buckets) * bucket_size);
        ofs.close();
    }


    /* #######################################################################
     *  Estimate Bucket Sizes
     * #################################################################### */
    static void EstimateBucketSizes(VarkitOptionsContainer &options) {
        std::cout << "Estimating bucket sizes (Bloom filter)" << std::endl;

        ProgressBar bar;

        std::string shape_str = options.shape_str;
        const size_t shape_length = options.ShapeLength();
        bool *shape = options.Shape();

        size_t total_size = options.ReferencesTotalSize();

        // Lambda for computing the space requirement in gb
        auto compute_size = [](size_t size) { return (size / 8) / ((double) 1024 * 1024 * 1024); };
        auto compute_size_raw = [](size_t size) { return (size / 8) / ((double) 1024 * 1024 * 1024); };

        // Set up bloom filter to estimate unique k-mer count and bucket sizes
        bloom_parameters bloom_params;
        bloom_params.projected_element_count = total_size * 3;
        bloom_params.false_positive_probability = 0.0000000000000001;
        bloom_params.random_seed = 0xA5A5A5A5;

        bool bloom_filter_success{false};
        size_t max_tries{20'000};
        size_t count_tries{0};
        double change_fp{1.005};
        double fp_lower{0.01};
        double change_element{0.8};

        bloom_filter *filter;

        std::cout << "Try to initialize bloom filter with limit " << options.memory_upper_gb_ << std::endl;
        std::cout << "Note: if " << options.memory_upper_gb_
                  << " is too high, specify the available RAM in GB through the options" << std::endl;

        while (!bloom_filter_success) {
            if (++count_tries == max_tries) {
                exit(1);
            }
            bloom_params.compute_optimal_parameters();

            if (compute_size(bloom_params.optimal_parameters.table_size) <= options.memory_upper_gb_) {
                bloom_filter_success = true;
            } else {
                if (bloom_params.false_positive_probability > fp_lower) {
                    bloom_params.projected_element_count *= change_element;
                } else {
                    bloom_params.false_positive_probability *= change_fp;
                }
                bloom_filter_success = false;
            }
        }

        std::cout << "Bloom filter false positive rate: " << bloom_params.false_positive_probability << std::endl;
        std::cout << "Bloom filter space requirements: " << compute_size(bloom_params.optimal_parameters.table_size)
                  << " GB" << std::endl;

        filter = new bloom_filter(bloom_params);
        filter->clear();

        // Important for initial db size
        double compensate_fp{5.f};

        const size_t block_size = (10 * 1024 * 1024);
        const size_t max_key = (1llu << options.key_bits_);
        const size_t bloom_threshold = 5;

        /* #########################################################
         *  Estimating Bucket size
         ######################################################## */

        ds::Benchmark bucket_bm("Estimating bucket size");
        bucket_bm.start();

        size_t *offsets = new std::size_t[options.bucket_count_];
        std::fill_n(offsets, options.bucket_count_, 0llu);

        size_t key;

        RollingKmerIterator iterator(options.k_, shape, shape_length);

        BufferedFastxReader reader;
        FastxRecord record;

        bar.reset(total_size);
        size_t progress = 0;

        for (auto file : options.references_) {
            ifstream is(file, ios::in);

            while (true) {
                bool ok = false;

                ok = reader.LoadBlock(is, block_size);
                if (!ok) break;

                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;

                    if (record.sequence.length() < shape_length) continue;

                    iterator.SetRecord(record);

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

                        if (key == 8344059844850) {
                            exit(9);
                        }

                        assert(key <= max_key);

                        uint64_t bucket_key = IndexedMap::ComputeOffsetKey(key, options.key_bits_);

                        assert(bucket_key < (1llu << options.offset_key_bits_));

                        if (!filter->contains(key)) {
                            ++offsets[bucket_key];
                            filter->insert(key);
                        }

                        // sometimes a k-mer is detected as seen in the bloom filter (fp)
                        // If the bucket of the FP is empty, we will run into a problem during inserting as
                        // there is at least one k-mer with that prefix and therefore the bucket size cant be 0;
                        // Account for that!
                        if (offsets[bucket_key] < bloom_threshold) {
                            ++offsets[bucket_key];
                        }
                    }
                }
                progress += reader.LastBlockSize();

                bar.Update(progress);
            }

            is.close();
        }

        std::cout << filter->size() << std::endl;
        delete filter;

        std::cout << "Estimating bucket size done." << std::endl;
        size_t total_capacity = 0;

        for (uint64_t i = 0; i < options.bucket_count_; i++) {
            offsets[i] *= compensate_fp;
            total_capacity += offsets[i];
        }

        std::cout << "Estimated size demands: " << total_capacity << std::endl;
        std::cout << "Estimated DB size: " << ((total_capacity * 8 * 1.25) / (1024 * 1024 * 1024.0)) << " GB"
                  << std::endl;

        bucket_bm.stop();
        bucket_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;

        std::cout << "Saving buckets.." << std::endl;

        // Save bucket sizes
        SaveBucketSizes(options.BucketsFile(), offsets, options.bucket_count_);

        std::cout << std::string(60, '-') << std::endl;

        delete[] offsets;
    }

    /* #######################################################################
     *  Determine Bucket Sizes
     * #################################################################### */
    static void DetermineBucketSizes(VarkitOptionsContainer &options) {
        std::cout << "Determining bucket sizes" << std::endl;

        ProgressBar bar;

        size_t k = options.k_;

        size_t *offsets = new std::size_t[options.bucket_count_];
        std::fill_n(offsets, options.bucket_count_, 0llu);


        std::string shape_str = options.shape_str;
        const size_t shape_length = options.ShapeLength();
        bool *shape = options.Shape();

        size_t total_size = options.ReferenceSize();
        size_t total_kmers = 0;


        const size_t block_size = (16 * 1024 * 1024);
        const uint64_t max_key = (1llu << options.key_bits_);

        const uint64_t key_bits_offset = options.offset_key_bits_;
        const uint64_t max_offset = (1llu << key_bits_offset);

        std::cout << "offset: " << key_bits_offset << std::endl;
        std::cout << "internal: " << options.internal_key_bits_ << std::endl;
        std::cout << "key total: " << options.key_bits_ << std::endl;
        std::cout << "value bits: " << options.value_bits_ << std::endl;

        /* #########################################################
         * # Estimating Bucket size
         ######################################################## */

        ds::Benchmark bucket_bm("Estimating bucket size (overestimate)");
        bucket_bm.start();

        omp_set_num_threads(options.threads_);

        ifstream* is = nullptr;

        size_t file_index = 0;

        bar.reset(options.ReferenceSize());
        size_t progress = 0;

#pragma omp parallel
        {
            size_t key;

            RollingKmerIterator iterator(k, shape, shape_length);
            BufferedFastxReader reader;
            FastxRecord record;

            is = new ifstream(options.ReferenceFile(), ios::in);
            while (true) {
                bool ok = false;
#pragma omp critical(reader)
                {
//                    if ((!is || is->eof()) && file_index < options.references_.size()) {
//                        if (is) is->close(); delete is;
//                        is = new ifstream(options.references_[file_index++], ios::in);
//                    }
                    if (is)
                        ok = reader.LoadBlock(*is, block_size);
                }
                if (!ok) break;

                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;
                    if (record.sequence.length() < shape_length) continue;
                    iterator.SetRecord(record);

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

                        size_t bucket_key = key >> options.internal_key_bits_;

                        assert(bucket_key < max_offset);
                        assert(key < (1llu << (options.key_bits_)));

#pragma omp atomic
                        offsets[bucket_key]++;
#pragma omp atomic
                        total_kmers++;
                    }
                }

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);
            }
        }
        delete is;

        size_t total_capacity = 0;
        for (int i = 0; i < options.bucket_count_; i++) {
            total_capacity += offsets[i];
        }

        std::cout << "Determined bucket size done." << std::endl;
        std::cout << "Determined size demands: " << total_kmers << std::endl;
        std::cout << "Estimates DB size: "
                  << ((total_capacity * 8 * (1.f / options.load_factor_)) / (1024 * 1024 * 1024.0)) << " GB"
                  << std::endl;

        bucket_bm.stop();
        bucket_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;

        // Save bucket sizes
        SaveBucketSizes(options.BucketsFile(), offsets, options.bucket_count_);

        std::cout << std::string(60, '-') << std::endl;

        delete[] offsets;
    }

    /* #######################################################################
     *  Calculate Bucket Sizes
     * #################################################################### */
    static void CalculateBucketSizes(VarkitOptionsContainer &options) {
        std::cout << "calculate bucket sizes" << std::endl;
        DetermineBucketSizes(options);
    }


    /* #######################################################################
     *  Build index Routine
     * #################################################################### */
    static void BuildIndex(VarkitOptionsContainer &options) {
        std::cout << "Start building Index" << std::endl;

        ds::Benchmark build_bm("Build process");

        ProgressBar bar;

        IndexedMap *map = new IndexedMap(options.internal_key_bits_,
                       options.value_bits_,
                       options.offset_key_bits_,
                       options.BucketsFile(),
                       options.load_factor_);

        options.SetMap(map);

        bar.reset(options.ReferenceSize());

        auto taxonomy = options.InternalTaxonomy();

        const size_t block_size = (10 * 1024 * 1024);

        ifstream *is = nullptr;
        size_t file_index = 0;

        omp_set_num_threads(options.threads_);

        size_t total_reads = 0;
        size_t skipped_reads = 0;
        size_t progress = 0;
        size_t shape_length = options.ShapeLength();

#pragma omp parallel
        {
            size_t key = 0;
            size_t internal_taxid;
            size_t marker_gene_id;
            size_t marker_gene_pos;

            FastxRecord record;
            BufferedFastxReader reader;

            RollingKmerIterator iterator(options.k_, options.Shape(), options.ShapeLength());
            SimpleMarkerKmerProcessor processor(map, taxonomy, options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);

            is = new ifstream(options.ReferenceFile(), ios::in);
            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                {
//                    if ((!is || is->eof()) && file_index < options.references_.size()) {
//                        if (is) is->close(); delete is;
//                        is = new ifstream(options.references_[file_index++], ios::in);
//                    }
                    if (is)
                        ok = reader.LoadBlock(*is, block_size);
                }
                if (!ok) break;


                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;

                    ++total_reads;

                    // Skip read if too short!
                    if (record.sequence.length() < shape_length) {
                        ++skipped_reads;
                        continue;has an initial distance of 0)
                    }

                    SimpleMarkerKmerProcessor::ParseHeader(record.header, internal_taxid, marker_gene_id);
                    if (!taxonomy->IsLeaf(internal_taxid)) {
                        std::cerr << "References are labelled incorrectly for header" << std::endl;
                        std::cerr << record.header << std::endl;
                        std::cerr << "denoted taxid " << internal_taxid << " has to be a leaf in the taxonomy as it represents a single biological genome" << std::endl;
                    }

                    // Pass new record on to iterator
                    iterator.SetRecord(record);

//                    std::cout << record.id << std::endl;

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

                        if (record.id == "2571_1") {
                            std::cerr << "2571_1\t" << iterator.GetPos() << '\t' << key << std::endl;
                        }

                        marker_gene_pos = iterator.GetPos();

                        if (marker_gene_pos > options.max_marker_gene_pos_) {
                            std::cerr << "Marker gene too long, just stop here: " << marker_gene_pos << " max: " << options.max_marker_gene_pos_ << std::endl;
                            break;
                        }

#pragma omp critical(put_operation)
                        {
                            processor.operator()(key, internal_taxid, marker_gene_id, marker_gene_pos);
                        }

//                        if (record.id.starts_with("21852_109") && iterator.GetPos() > 4186 && iterator.GetPos() < 4190) {
//                            std::cout << "record.id.starts_with(\"21852_109\") && iterator.GetPos() > 4260 && iterator.GetPos() < 4265" << std::endl;
//                            std::cout << record.to_string() << std::endl;
//                            std::cout << record.sequence.substr(iterator.GetPos(), 31) << std::endl;
//                            std::cout << iterator.GetPos() << " " << ShapeUtils::ApplyShape(record.sequence, iterator.GetPos(), options.Shape(), iterator.ShapeSize(), true) << " " << key << std::endl;
//                            std::string stop;
//                            std::cin >> stop;
//                        }

//                        if (key == 8344059844850) {
//                            std::cout << std::string(79, '-') << std::endl;
//                            std::cout << record.id << std::endl;
//                            std::cout << 8344059844850 << std::endl;
//                            std::cout << ShapeUtils::ApplyShape(record.sequence, iterator.GetPos(), options.Shape(), iterator.ShapeSize(), true) << std::endl;
//                            std::cout << iterator.GetPos() << std::endl;
//                            auto entry = map->Find(key);
//
//                            if (entry) {
//                                std::cout << "IteratorPos: " << iterator.GetPos() << std::endl;
//                                auto value = entry->value(map->value_bits_);
//                                classification::ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);
//                                std::cout << "LOOKUP: " << extractor.GetTaxid(value) << " " << extractor.GetGeneId(value) << " " << extractor.GetGenePos(value) << " " << ShapeUtils::ApplyShape(record.sequence, iterator.GetPos(), options.Shape(), options.shape_str.length(), true) << std::endl;
//                            }
//
//                            std::string stop;
//                            std::cin >> stop;
//                        }
                    }
                }

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);
            }
        };

        std::cout << "Save database ... ";
        map->Save(options.DatabaseFile());
        std::cout << "done." << std::endl;
    }

    static void Verify(VarkitOptionsContainer &options) {
        std::cout << "Verify build process .. " << std::endl;

        ds::Benchmark build_bm("Verify build process");

        ProgressBar bar;

        IndexedMap *map;
        if (options.map == nullptr) {
            map = IndexedMap::Load(options.DatabaseFile());
            options.SetMap(map);
        } else {
            map = options.map;
        }


        bar.reset(options.ReferencesTotalSize());

        auto taxonomy = options.InternalTaxonomy();

        const size_t block_size = (10 * 1024 * 1024);

        ifstream *is = nullptr;
        size_t file_index = 0;

        omp_set_num_threads(options.threads_);

        size_t total_reads = 0;
        size_t skipped_reads = 0;
        size_t progress = 0;
        size_t shape_length = options.ShapeLength();

#pragma omp parallel
        {
            size_t key = 0;
            size_t internal_taxid;
            size_t marker_gene_id;
            size_t marker_gene_pos;

            FastxRecord record;
            BufferedFastxReader reader;

            RollingKmerIterator iterator(options.k_, options.Shape(), options.ShapeLength());
            classification::ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);

            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                {
                    if ((!is || is->eof()) && file_index < options.references_.size()) {
                        if (is) is->close(); delete is;
                        is = new ifstream(options.references_[file_index++], ios::in);
                    }
                    if (is)
                        ok = reader.LoadBlock(*is, block_size);
                }
                if (!ok) break;


                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;

                    ++total_reads;

                    // Skip read if too short!
                    if (record.sequence.length() < shape_length) {
                        ++skipped_reads;
                        continue;
                    }

                    SimpleMarkerKmerProcessor::ParseHeader(record.header, internal_taxid, marker_gene_id);
                    if (!taxonomy->IsLeaf(internal_taxid)) {
                        std::cerr << "References are labelled incorrectly for header" << std::endl;
                        std::cerr << record.header << std::endl;
                        std::cerr << "denoted taxid " << internal_taxid << " has to be a leaf in the taxonomy as it represents a single biological genome" << std::endl;
                    }

                    // Pass new record on to iterator
                    iterator.SetRecord(record);

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

                        marker_gene_pos = iterator.GetPos();

                        if (marker_gene_pos > options.max_marker_gene_pos_) {
                            std::cerr << "Marker gene too long, just stop here: " << marker_gene_pos << " max: " << options.max_marker_gene_pos_ << std::endl;
                            break;
                        }

                        auto found = map->Find(key);

                        if (!found) {
                            std::cerr << "key " << key << " in " << record.header << " not found." << std::endl;
                            exit(26);
                        }
                        auto value = found->value(options.value_bits_);
                        auto found_taxid = extractor.GetTaxid(value);
                        auto found_geneid = extractor.GetGeneId(value);
                        auto found_genepos = extractor.GetGenePos(value);

                        if (!taxonomy->IsNodeAncestor(found_taxid, internal_taxid)) {
                            std::cerr << "Faulty cell." << std::endl;
                            exit(26);
                        }
                        if (taxonomy->IsLeaf(internal_taxid) && found_geneid > 0) {
                            if (found_geneid != marker_gene_id) {
                                std::cout << record.id << std::endl;
                                std::cout << found_taxid << " " << internal_taxid << std::endl;
                                std::cout << found_geneid << " " << marker_gene_id << std::endl;
                                std::cout << found_genepos << " " << marker_gene_pos << std::endl;
                                exit(30);
                            }
                        }
                    }
                }

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);
            }
        };

        std::cout << "done." << std::endl;
    }

    static void ComputeAbundanceCorrection(VarkitOptionsContainer &options) {
        IndexedMap* map = options.map;
        ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);

        size_t total_size = 0;
        for (auto file : options.references_)
            total_size += Utils::GetFileSize(file);
        const size_t block_size = (1024 * 1024);

        std::string shape_str = options.shape_str;
        const size_t shape_length = shape_str.length();
        bool * shape = ShapeUtils::GetShape(shape_str);

        size_t progress = 0;
        omp_set_num_threads(options.threads_);

        std::unordered_map<size_t, std::unordered_map<size_t, size_t>> taxon_total_hits;
        std::unordered_map<size_t, std::unordered_map<size_t, size_t>> taxon_classified_hits;
        std::ofstream correction_out(options.AbundanceCorrectionFile(), ios::out);

        ProgressBar bar;
        bar.reset(total_size);


        ifstream is(options.references_[0], ios::in);

        std::cout << "Compute abundance correction values and save to " << options.AbundanceCorrectionFile() << std::endl;

        PatternDB pattern_db(options.ShapeDatabaseFile());
        auto taxonomy = options.InternalTaxonomy();

#pragma omp parallel
        {
            RollingKmerIterator iterator(options.k_, shape, shape_length);
            BufferedFastxReader reader;
            FastxRecord record;

            ReadClassifier read_classifier(extractor, *taxonomy, shape_length);
            read_classifier.InitProcessor(&pattern_db, true);
            std::vector<size_t> values;

            uint64_t key;

            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                ok = reader.LoadBlock(is, block_size);
                if (!ok) break;

                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;

                    // continue if sequence too short.
                    if (record.sequence.length() < shape_length) continue;

                    iterator.SetRecord(record);

                    auto tokens = Utils::split(record.id, "_");
                    size_t true_taxid = stoul(tokens[0]);
                    size_t true_geneid = stoul(tokens[1]);

                    size_t leaf_hit_count = 0lu;
                    size_t total_count = 0lu;

                    values.clear();


                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);


                        if (iterator.GetPos() > options.max_marker_gene_pos_) {
                            std::cerr << "Marker gene too long, just stop here: " << iterator.GetPos() << " max: " << options.max_marker_gene_pos_ << std::endl;
                            break;
                        }

                        auto entry = map->Find(key);

                        if (entry) {
                            auto value = entry->value(map->value_bits_);
                            values.push_back(value);
                        } else {
                            std::cerr << "Oops. Something went wrong during database construction." << std::endl;
                            std::cerr << "Key " << key << " not found. Occurred at position " << iterator.GetPos() << " in the sequence " << record.header << ". Sequence context: " << record.sequence.substr(iterator.GetPos(), iterator.GetPos() + 100 < record.sequence.length() ? 100 : record.sequence.length() - iterator.GetPos()) << std::endl;
                            exit(9);
                        }
                        total_count++;
                    }

//                    std::cout << "values size " << values.size() << std::endl;

                    // Go over windows
                    int window_size = 100 - options.ShapeLength() + 1;
                    int min_overlap = 20;
                    int step_size = 1;


                    for (int i = min_overlap - window_size; i < ((int) values.size() + (window_size - min_overlap)); i += step_size) {
#pragma omp critical(taxon_hits)
                        {
                            // Commented out for DEBUG
                            if (!taxon_total_hits.contains(true_taxid)) {
                                std::unordered_map<size_t, size_t> inner;
                                taxon_total_hits.insert( { true_taxid, inner } );
                            }
                            if (!taxon_total_hits[true_taxid].contains(true_geneid)) {
                                taxon_total_hits[true_taxid][true_geneid] = 0;
                            }
                            taxon_total_hits[true_taxid][true_geneid]++;

                            // Read correctly classified
                            if (!taxon_classified_hits.contains(true_taxid)) {
                                std::unordered_map<size_t, size_t> inner;
                                taxon_classified_hits.insert( { true_taxid, inner } );
                            }
                            if (!taxon_classified_hits[true_taxid].contains(true_geneid)) {
                                taxon_classified_hits[true_taxid][true_geneid] = 0;
                            }
                        }

                        if (!taxon_classified_hits.contains(true_taxid)) {
                            std::cout << "SOMETHING IS OFF " << true_taxid << taxon_classified_hits.contains(true_taxid) << std::endl;
                        }

                        read_classifier.Reset();
                        read_classifier.SetLength(window_size);

                        // Simulated read
                        for (auto rpos = i; rpos < i + window_size; rpos++) {
                            if (rpos < 0 || rpos >= values.size())
                                continue;
                            read_classifier.AddHit(values[rpos]);
                        }

                        size_t taxid = 0, best_count = 0, leaf_hits = 0, best_total_counts = 0;
                        int gene_pos = -1, pattern_idx = -1, gene_id = -1;
                        double rank_confidence = 0, candidate_confidence = 0, classification_confidence = 0;
                        bool forward = true;

                        read_classifier.EvaluateTreeReturnBest(taxid, gene_id, gene_pos, forward, best_count, leaf_hits, best_total_counts, pattern_idx);
                        gene_pos -= !forward * record.sequence.length();
                        bool is_classified = (taxid != 0);

//                        if (!taxon_classified_hits.contains(true_taxid)) {
//                            std::cout << "ASDLASKJDSAD: " << true_taxid << std::endl;
//                        }

                        if (is_classified) {
                            // Commented out for DEBUG
                            auto x = 0;
                            if (true_taxid == taxid) {
#pragma omp critical(taxon_hits)
                                {
//                                    if (!taxon_classified_hits.contains(true_taxid)) {
//                                        std::cout << "indicator taxid: " << true_taxid << std::endl;
//                                    } else if (!taxon_classified_hits[true_taxid].contains(true_geneid)) {
//                                        std::cout << "indicator gene: " << true_geneid << std::endl;
//                                    }
                                    taxon_classified_hits[true_taxid][true_geneid]++;
                                }
                            }
                        }
                    }


                } // End Read Block

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);

            } // End Reader

        } // Close OMP parallel

    // Commented out for DEBUG
        for (auto& taxon : taxon_total_hits) {
//            std::cout << "taxon: " << taxon.first << std::endl;
            auto taxid = taxon.first;
            for (auto& gene : taxon.second) {
                auto geneid = gene.first;
                auto total_hits = gene.second;
                auto true_hits = taxon_classified_hits[taxid][geneid];
                correction_out << taxid << '\t' << geneid << '\t' << (double)true_hits/total_hits << '\n';
            }
        }
        correction_out.close();
    }

    static void EvaluateIndexRankComposition(VarkitOptionsContainer &options) {
        std::cout << "Evaluate Rank Level Composition of Database\n" << std::endl;

        ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);

        if (!options.IsMapInitialized()) {
            std::cerr << "First initialize map in options." << std::endl;
            exit(9);
        }

        auto map = options.map;
        Cell *cells = map->Map();

        tsl::robin_map<std::string, size_t> rank_counts;

        size_t total = 0;
        size_t leaf_count = 0;
        size_t unique_count = 0;
        size_t non_leaf_with_gid = 0;

        for (size_t i = 0; i < map->Capacity(); i++) {
            auto& cell = cells[i];

            if (!cell.empty()) {
                auto key = cell.key(map->value_bits_);
                auto value = cell.value(map->value_bits_);

                total++;

                size_t internal_taxid = extractor.GetTaxid(value);
                size_t marker_gene_id = extractor.GetGeneId(value);
                size_t marker_gene_pos = extractor.GetGenePos(value);

                std::string rank = options.taxonomy->Get(internal_taxid).rank;

                if (rank.empty()) {
                    std::cout << "Empty:" << std::endl;
                    std::cout << "key: " << key << std::endl;
                    std::cout << "value: " << value << std::endl;
                    std::cout << "internal_taxid: " << internal_taxid << std::endl;
                }
                if (!options.taxonomy->Get(internal_taxid).IsLeaf() && marker_gene_id > 0) {
                    non_leaf_with_gid++;
//                    printf("%zu\t%zu\t%zu\t%s\n", internal_taxid, marker_gene_id, marker_gene_pos, options.taxonomy->Get(internal_taxid).scientific_name.c_str());
                }

                assert(!rank.empty());

                if (!rank_counts.contains(rank)) {
                    rank_counts.insert({rank, 0});
                }

                if (options.taxonomy->Get(internal_taxid).IsLeaf()) {
                    leaf_count++;
                    if (marker_gene_id > 0) {
                        unique_count++;
                    }
                }

                rank_counts[rank]++;
            }
        }

        for (auto pair : rank_counts) {
            std::cout << "Rank: " << pair.first << "\t" << pair.second << '\t' << ((double) pair.second / total)
                      << std::endl;
        }
        std::cout << "------------" << std::endl;
        std::cout << "leaves: " << (double) leaf_count / total << std::endl;
        std::cout << "unique: " << (double) unique_count / total << std::endl;
        std::cout << "non leaf with gene id information: " << non_leaf_with_gid << std::endl;
        std::cout << "total k-mers: " << total << std::endl;
    }

    static void Run(VarkitOptionsContainer &options) {
        ds::Benchmark total_bm("Total process");
        total_bm.start();

        ProgressBar bar;

        // Sanity check
        if (!options.IsValidDBPreBuild()) {
            std::cerr << "Error" << std::endl;
            exit(9);
        }

        if (!options.HasBuckets() || options.force_rebuild_) {
            CalculateBucketSizes(options);
            std::cout << DIVIDER << std::endl;
        }

        if (!options.HasDatabase() || options.force_rebuild_) {
            BuildIndex(options);
            std::cout << "Save bucket sizes .. ";
            options.map->SaveBucketSizes(options.BucketsFile());
            std::cout << "done." << std::endl;
            std::cout << DIVIDER << std::endl;
        }

        if (options.HasDatabase() && !options.IsMapInitialized()) {
            options.LoadResources();
        }

//        options.map->print(options.map->offset_[42007363], 20, true);

//        options.map = nullptr;
//        std::cout << std::string(60, '-') << std::endl;
//        std::cout << "Rebuild condense" << std::endl;
//        BuildIndex(options);
//        options.map->PrintMeta();

//        EvaluateIndexRankComposition(options);
//        std::cout << DIVIDER << std::endl;
//
//        Verify(options);
//        std::cout << DIVIDER << std::endl;
//
//        options.map->PrintMeta();
//        std::cout << DIVIDER << std::endl;

//        std::cout << "Compute abundance correction for DB" << std::endl;
//        ComputeAbundanceCorrection(options);
        std::cout << DIVIDER << std::endl;

        total_bm.stop();
        total_bm.printResults();
    }
};