//
// Created by fritsche on 11/03/2021.
//

#pragma once

#include <cstdio>
#include <string>
#include <vector>
#include "bloom_filter.h"
#include <err.h>
#include <sysexits.h>
#include <FastxReader.h>
#include <fstream>
#include <SpacedKmerIterator.h>
#include <SimpleKmerIterator.h>
#include <SimpleKmerPutter.h>
#include <omp.h>
#include "Utils.h"
#include "progress_bar.h"
#include "compact_map.h"
#include "ShapeUtils.h"
#include "TaxonomyNew.h"
#include "robin_set.h"
#include "robin_map.h"
//#include "multimap.h"
#include <utils/Benchmark.h>

#define USE_OMP 0
#define ASSERT 1


namespace Build {

    enum CellType { UINT16 = 16, UINT32 = 32, UINT64 = 64 };

    inline static const std::string TAXONOMY_NODES = "nodes.dmp";
    inline static const std::string TAXONOMY_NAMES = "names.dmp";
    inline static const std::string TAXONOMY_CUSTOM = "internal_taxonomy.dmp";
    inline static const std::string DATABASE = "database.bin";
    inline static const std::string SHAPE = "shape.txt";
    inline static const std::string PATTERN_DATABASE = "shape.bin";
    inline static const std::string ENTRY_META = "entry.meta";

    struct BuildOptions {
        BuildOptions(
                const std::string db_dir,
                const std::string shape,
                const size_t value_bits,
                const size_t key_bit_total,
                const size_t threads,
                const size_t bucket_bits,
                const size_t key_bits_internal,
                const std::vector<std::string> &references,
                const CellType cell_type,
                const double load_factor) :
            value_bits_(value_bits),
            key_bits_total_(key_bit_total),
            db_dir_(db_dir),
            k_ (key_bit_total/2),
            shape_ (shape),
            shape_path_(db_dir + "/" + SHAPE),
            db_path_(db_dir + "/" + DATABASE),
            custom_taxonomy_path_(db_dir + "/" + TAXONOMY_CUSTOM),
            bucket_bits_ (bucket_bits),
            bucket_count_(1llu << bucket_bits),
            key_bits_internal_(key_bits_internal),
            threads_(threads),
            references_(references),
            cell_type_(cell_type),
            load_factor_(load_factor)
            {};


    public:
        const uint32_t threads_;
        const std::string db_dir_;
        const std::string db_path_;
        const std::string shape_path_;
        const std::string entry_meta_path_;
        const std::string custom_taxonomy_path_;

        const std::string shape_;

        const size_t key_bits_internal_;
        const size_t value_bits_;
        const size_t key_bits_total_;
        const size_t k_;

        const double load_factor_;

        const size_t bucket_bits_;
        const size_t bucket_count_;
        const std::vector<std::string> references_;

        const CellType cell_type_;


    };

    struct BuildOptionsMG {
    public:
        const uint32_t threads_;

        const std::string working_directory_;
        const std::string db_path_;
        const std::string pattern_db_path_;
        const std::string taxonomy_nodes_path_;
        const std::string taxonomy_names_path_;
        const std::string custom_taxonomy_;
        const std::string shape_path_;
        const std::string entry_meta_path_;

        const std::string target_rank_;

        const size_t key_bits_internal_;
        const size_t value_bits_;
        const size_t key_bits_total_;
        const size_t k_;

        const size_t bucket_bits_;
        const size_t bucket_count_;

        const double load_factor_;
        const size_t bloom_limit_gb_;

        const bool eval_;
        const bool update_;

        const std::vector<std::string> references_;

        BuildOptionsMG(size_t key_bits_internal,
                       size_t value_bits,
                       size_t key_bits_total,
                       size_t k,
                       uint32_t threads,
                       double load_factor,
                       size_t bloom_limit_gb,
                       size_t bucket_bits,
                       std::string working_directory,
                       bool eval,
                       std::string target_rank,
                       bool update,
                       std::vector<std::string> &references) :
                key_bits_internal_(key_bits_internal),
                key_bits_total_(key_bits_total),
                value_bits_(value_bits),
                bucket_bits_(bucket_bits),
                k_(k),
                eval_(eval),
                update_(update),
                bucket_count_(1llu << bucket_bits),
                references_(references),
                threads_(threads),
                target_rank_(target_rank),
                bloom_limit_gb_(bloom_limit_gb),
                load_factor_(load_factor),
                working_directory_(working_directory),
                custom_taxonomy_(working_directory + "/" + TAXONOMY_CUSTOM),
                taxonomy_names_path_(working_directory + "/" + TAXONOMY_NAMES),
                taxonomy_nodes_path_(working_directory + "/" + TAXONOMY_NODES),
                db_path_(working_directory + "/" + DATABASE),
                pattern_db_path_(working_directory + "/" + PATTERN_DATABASE),
                shape_path_(working_directory + "/" + SHAPE),
                entry_meta_path_(working_directory + "/" + ENTRY_META) {}

        void Print() {
            std::cout << "key_bits_total:   \t" << key_bits_total_ << std::endl;
            std::cout << "bucket_bits   :   \t" << bucket_bits_ << std::endl;
            std::cout << "key_internal_bits:\t" << key_bits_internal_ << std::endl;
            std::cout << "value_bits:       \t" << value_bits_ << std::endl;
            std::cout << "k:                \t" << k_ << std::endl;
            std::cout << "bucket_size:      \t" << bucket_count_ << std::endl;
        }

        bool Check() {
            if ((key_bits_internal_ + bucket_bits_) != key_bits_total_) {
                errx(EX_CONFIG,
                     "Internal key bit size (%llu) + bucket bit size (%llu) must be equal to the total key bit size (%llu)",
                     (unsigned long long) key_bits_internal_,
                     (unsigned long long) bucket_bits_,
                     (unsigned long long) key_bits_total_);
            }
            if ((key_bits_total_ != (2 * k_))) {
                errx(EX_CONFIG,
                     "Total key bit size (%llu) must be equal to 2 * k-mer size (%llu)",
                     (unsigned long long) key_bits_total_,
                     (unsigned long long) k_);
            }
            if ((64 - key_bits_internal_) != value_bits_) {
                errx(EX_CONFIG,
                     "Entry size 64 - internal key bit size (%llu) must be equal to value bit size (%llu)",
                     (unsigned long long) key_bits_internal_,
                     (unsigned long long) value_bits_);
            }
//            if (!Utils::exists(db_path_)) {
//                errx(EX_IOERR,
//                     "Database file %s does not exist.", db_path_.c_str());
//            }
            if (!Utils::exists(custom_taxonomy_)) {
                errx(EX_IOERR,
                     "Taxonomy names file %s does not exist.", taxonomy_names_path_.c_str());
            }
//            if (!Utils::exists(taxonomy_nodes_path_)) {
//                errx(EX_IOERR,
//                     "Taxonomy nodes file %s does not exist.", taxonomy_nodes_path_.c_str());
//            }
            for (auto file : references_) {
                if (!Utils::exists(file)) {
                    errx(EX_IOERR,
                         "Reference file %s does not exist.", file.c_str());
                }
            }

            return true;
        }


    };


    struct BuildOptionsMMG {
    public:
        const uint32_t threads_;

        const std::string working_directory_;
        const std::string db_path_;
        const std::string pattern_db_path_;
        const std::string taxonomy_nodes_path_;
        const std::string taxonomy_names_path_;
        const std::string custom_taxonomy_;
        const std::string shape_path_;
        const std::string entry_meta_path_;


        const size_t key_bits_internal_;
        const size_t value_bits_;
        const size_t key_bits_total_;
        const size_t k_;

        const size_t bucket_bits_;
        const size_t bucket_count_;

        const double load_factor_;
        const size_t bloom_limit_gb_;

        const std::vector<std::string> references_;

        BuildOptionsMMG(size_t key_bits_internal,
                       size_t value_bits,
                       size_t key_bits_total,
                       size_t k,
                       uint32_t threads,
                       double load_factor,
                       size_t bloom_limit_gb,
                       size_t bucket_bits,
                       std::string working_directory,
                       std::vector<std::string> &references) :
                key_bits_internal_(key_bits_internal),
                key_bits_total_(key_bits_total),
                value_bits_(value_bits),
                bucket_bits_(bucket_bits),
                k_(k),
                bucket_count_(1llu << bucket_bits),
                references_(references),
                threads_(threads),
                bloom_limit_gb_(bloom_limit_gb),
                load_factor_(load_factor),
                working_directory_(working_directory),
                custom_taxonomy_(working_directory + "/" + TAXONOMY_CUSTOM),
                taxonomy_names_path_(working_directory + "/" + TAXONOMY_NAMES),
                taxonomy_nodes_path_(working_directory + "/" + TAXONOMY_NODES),
                db_path_(working_directory + "/" + DATABASE),
                pattern_db_path_(working_directory + "/" + PATTERN_DATABASE),
                shape_path_(working_directory + "/" + SHAPE),
                entry_meta_path_(working_directory + "/" + ENTRY_META) {}

        void Print() {
            std::cout << "key_bits_total:   \t" << key_bits_total_ << std::endl;
            std::cout << "bucket_bits   :   \t" << bucket_bits_ << std::endl;
            std::cout << "key_internal_bits:\t" << key_bits_internal_ << std::endl;
            std::cout << "value_bits:       \t" << value_bits_ << std::endl;
            std::cout << "k:                \t" << k_ << std::endl;
            std::cout << "bucket_size:      \t" << bucket_count_ << std::endl;
        }

        bool Check() {
            if ((key_bits_internal_ + bucket_bits_) != key_bits_total_) {
                errx(EX_CONFIG,
                     "Internal key bit size (%llu) + bucket bit size (%llu) must be equal to the total key bit size (%llu)",
                     (unsigned long long) key_bits_internal_,
                     (unsigned long long) bucket_bits_,
                     (unsigned long long) key_bits_total_);
            }
            if ((key_bits_total_ != (2 * k_))) {
                errx(EX_CONFIG,
                     "Total key bit size (%llu) must be equal to 2 * k-mer size (%llu)",
                     (unsigned long long) key_bits_total_,
                     (unsigned long long) k_);
            }
            if ((62 - key_bits_internal_) != value_bits_) {
                errx(EX_CONFIG,
                     "Entry size 62 - internal key bit size (%llu) must be equal to value bit size (%llu)",
                     (unsigned long long) key_bits_internal_,
                     (unsigned long long) value_bits_);
            }
            for (auto file : references_) {
                if (!Utils::exists(file)) {
                    errx(EX_IOERR,
                         "Reference file %s does not exist.", file.c_str());
                }
            }

            return true;
        }
    };

    struct ValueExtractorMG {
    public:
        ValueExtractorMG(size_t internal_taxid_bits,
                         size_t gene_id_bits,
                         size_t gene_pos_bits) :
                internal_taxid_bits_(internal_taxid_bits),
                gene_id_bits_(gene_id_bits),
                gene_pos_bits_(gene_pos_bits) {}

        ValueExtractorMG(const ValueExtractorMG &t) {
            internal_taxid_bits_ = t.internal_taxid_bits_;
            gene_id_bits_ = t.gene_id_bits_;
            gene_pos_bits_ = t.gene_pos_bits_;
        }


        static ValueExtractorMG Load(std::string path) {
            std::ifstream is(path);
            std::string line;

            size_t internal_taxid_bits;
            size_t gene_id_bits;
            size_t gene_pos_bits;

            if (getline(is, line)) {
                internal_taxid_bits = stoull(line);
            }
            if (getline(is, line)) {
                gene_id_bits = stoull(line);
            }
            if (getline(is, line)) {
                gene_pos_bits = stoull(line);
            }
            is.close();

            return ValueExtractorMG(internal_taxid_bits, gene_id_bits, gene_pos_bits);
        }

        void Save(std::string path) {
            std::ofstream ofs(path);

            ofs << internal_taxid_bits_ << '\n';
            ofs << gene_id_bits_ << '\n';
            ofs << gene_pos_bits_;

            ofs.close();
        }

        size_t GetTaxid(size_t value) {
            return value >> (gene_id_bits_ + gene_pos_bits_);
        }

        size_t GetGeneId(size_t value) {
            return value >> (gene_pos_bits_) & ((1 << gene_id_bits_) - 1);
        }

        size_t GetGenePos(size_t value) {
            return value & ((1 << gene_pos_bits_) - 1);
        }


    public:
        size_t internal_taxid_bits_;
        size_t gene_id_bits_;
        size_t gene_pos_bits_;
    };



    static void SaveBucketSizes(std::string path, size_t* buckets, size_t bucket_size) {
        ofstream ofs(path, ios::binary);
        ofs.write((char *) &bucket_size, sizeof(bucket_size));
        ofs.write((char *) buckets, sizeof(buckets) * bucket_size);
        ofs.close();
    }
//
//    static void ChangeBucketSizes(std::string path, std::vector<size_t> targets, double factor) {
//        auto result = IndexedMap::LoadBucketSizes(path);
//
//        for (auto target : targets) {
//            result.first[target] *= factor;
//        }
//
//        SaveBucketSizes(path, result.first, result.second);
//        delete[] result.first;
//    }

    size_t *EstimateBucketSizes(BuildOptionsMG &options) {
        std::cout << "Estimating bucket sizes (Bloom filter)" << std::endl;

        ProgressBar bar;

        std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        const size_t shape_length = shape_str.length();
        bool *shape = ShapeUtils::GetShape(shape_str);

        size_t total_size = 0;

        for (auto file : options.references_) {
            total_size += Utils::GetFileSize(file);
        }



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

        std::cout << "Try to initialize bloom filter with limit " << options.bloom_limit_gb_ << std::endl;
        std::cout << "Note: if " << options.bloom_limit_gb_
                  << " is too much, specify the available RAM in GB through the options" << std::endl;

        while (!bloom_filter_success) {
            if (++count_tries == max_tries) {
                exit(1);
            }
            bloom_params.compute_optimal_parameters();

            if (compute_size(bloom_params.optimal_parameters.table_size) <= options.bloom_limit_gb_) {
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
        const uint64_t max_key = (1llu << options.key_bits_total_);

//        std::cout << "key_bits_total: " << options.key_bits_total_ << std::endl;
//        std::cout << "max_key: " << max_key << std::endl;


        /* #########################################################
         * # Estimating Bucket size
         ######################################################## */

        ds::Benchmark bucket_bm2("Hello");
        ds::Benchmark bucket_bm("Estimating bucket size");
        bucket_bm.start();


        size_t *offsets = new std::size_t[options.bucket_count_];
        std::fill_n(offsets, options.bucket_count_, 0llu);

        omp_set_num_threads(3);

#if USE_OMP
#pragma omp parallel
        {

            const size_t buffer_size = 1000;
            uint64_t buffer[buffer_size];
            size_t buffer_pos = 0;
#endif
        uint64_t key;

//        SimpleKmerIterator iterator(options.k_, shape, shape_length);
        RollingKmerIterator iterator(options.k_, shape, shape_length);

        BufferedFastxReader reader;
        FastxRecord record;

        bar.reset(total_size);
        size_t progress = 0;

        for (auto file : options.references_) {
            ifstream is(file, ios::in);

            while (true) {
                bool ok = false;
#if USE_OMP
#pragma omp critical(reader)
#endif
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

                        if (key == UINT64_MAX) {
                            continue;
                        }

#if ASSERT
                        assert(key <= max_key);
#endif

#if USE_OMP
                        if (buffer_pos >= buffer_size) {
                                for (uint64_t i = 0; i < buffer_size; i++) {
                                    uint64_t tmp_key = buffer[i];

#if ASSERT
                                    assert(tmp_key <= max_key);
#endif

                                    uint64_t bucket_key = IndexedMap::ComputeOffsetKey(tmp_key, options.key_bits_internal_);
#pragma omp critical(filter)
                                    {
                                        if (!filter->contains(tmp_key)) {
                                            ++offsets[bucket_key];
                                            filter->insert(tmp_key);
                                        }
                                        // sometimes a k-mer is detected as seen in the bloom filter (fp)
                                        // If the bucket of the FP is empty, we will run into a problem during inserting as
                                        // there is at least one k-mer with that prefix and therefore the bucket size cant be 0;
                                        // Account for that!

                                        if (bucket_key == 1450) {
                                            std::cerr << "size: " << offsets[bucket_key] << std::endl;
                                        }

                                        if (offsets[bucket_key] < 10) {
                                            ++offsets[bucket_key];
                                        }
                                    }
                                }
                                buffer_pos = 0;
                            }

                            buffer[buffer_pos] = key;
                            buffer_pos++;
#else
                        uint64_t bucket_key = IndexedMap::ComputeOffsetKey(key, options.key_bits_internal_);

                        assert(bucket_key < (1llu << options.bucket_bits_));

                        if (!filter->contains(key)) {
                            ++offsets[bucket_key];
                            filter->insert(key);
                        }

                        // sometimes a k-mer is detected as seen in the bloom filter (fp)
                        // If the bucket of the FP is empty, we will run into a problem during inserting as
                        // there is at least one k-mer with that prefix and therefore the bucket size cant be 0;
                        // Account for that!
                        if (offsets[bucket_key] < 5) {
                            ++offsets[bucket_key];
                        }
#endif
                    }
                }

#if USE_OMP
#pragma omp atomic
#endif
                progress += reader.LastBlockSize();
#if USE_OMP
#pragma omp critical(updatebar)
#endif
                bar.Update(progress);
            }

            is.close();
        }
#if USE_OMP
        }
#endif
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
        ds::Benchmark savebuckets_bm("Saving buckets");
        savebuckets_bm.start();

        // Save bucket sizes
        SaveBucketSizes(options.db_path_ + ".buckets", offsets, options.bucket_count_);

        std::cout << "Saving buckets done" << std::endl;
        savebuckets_bm.stop();
        savebuckets_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;

        return offsets;
    }

    static pair<size_t*, size_t> LoadBucketSizes(std::string file) {
        ifstream ifs(file, ios::binary);
        if (!ifs.is_open()) {
            std::cerr << "replace by exception" << std::endl;
        }

        size_t bucket_size;
        ifs.read((char *) &bucket_size, sizeof(bucket_size));

        size_t* buckets = new size_t[bucket_size];

        ifs.read((char *) buckets, sizeof(buckets) * bucket_size);
        ifs.close();

        return { buckets, bucket_size };
    }

    size_t *DetermineBucketSizes(BuildOptions &options) {
        std::cout << "Determining bucket sizes" << std::endl;

        ProgressBar bar;

        std::string shape_str = options.shape_;//ShapeUtils::LoadShape(options.shape_path_);
        const size_t shape_length = shape_str.length();
        bool *shape = ShapeUtils::GetShape(shape_str);

        size_t total_size = 0;

        for (auto file : options.references_) {
            total_size += Utils::GetFileSize(file);
        }

        // Lambda for computing the space requirement in gb
        auto compute_size = [](size_t size) { return (size / 8) / ((double) 1024 * 1024 * 1024); };

        // power with base two with size_t
        auto two_pow = [](size_t exp) {
            assert(exp >= 0);
            size_t res = 1;
            while (exp-- > 0) res *= 2;
            return res;
        };

        size_t element_count = two_pow(options.key_bits_total_);
        size_t filter_size = element_count / 32 + ((bool) (element_count % 32)) * 1;
        uint32_t *filter = new uint32_t[filter_size];

        // insert k-mer key into filter
        auto insert_filter = [filter, element_count](size_t kmer) {
            assert(kmer < element_count);
#pragma omp atomic
            filter[kmer / 32] |= 1 << (kmer % 32);
        };
        auto retrieve_filter = [filter, element_count](size_t kmer) {
            assert(kmer < element_count);
            return filter[kmer / 32] & (1 << (kmer % 32));
        };

//        const size_t block_size = (10 * 1024 * 1024);
        const size_t block_size = (16 * 1024 * 1024);
        const uint64_t max_key = (1llu << options.key_bits_total_);

        /* #########################################################
         * # Estimating Bucket size
         ######################################################## */

        ds::Benchmark bucket_bm("Estimating bucket size");
        bucket_bm.start();

        omp_set_num_threads(options.threads_);

        ifstream* is = nullptr;

        size_t file_index = 0;

        bar.reset(total_size);
        size_t progress = 0;

//        for (auto file : options.references_) {
//            is = new ifstream(file, ios::in);

#pragma omp parallel
            {
//                const size_t buffer_size = 1000;
//                uint64_t buffer[buffer_size];
//                size_t buffer_pos = 0;

                uint64_t key;
//                 SimpleKmerIterator iterator(options.k_, shape, shape_length);
                RollingKmerIterator iterator(options.k_, shape, shape_length);
                BufferedFastxReader reader;
                FastxRecord record;


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
//                    ok = reader.LoadBlock(*is, block_size);
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

//                            if (key == 3761522364695) {
//                                std::cout << record.header << std::endl;
//                            }

                            if (key == UINT64_MAX) {
                                continue;
                            }
                            insert_filter(key);
                            assert(retrieve_filter(key));
                        }
                    }

#pragma omp atomic
                    progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                    bar.Update(progress);
                }
            }

//        }

        size_t *offsets = new std::size_t[options.bucket_count_];
        std::fill_n(offsets, options.bucket_count_, 0llu);

        size_t total_capacity = 0;
        for (size_t i = 0; i < element_count; i++) {
            if (retrieve_filter(i)) {
                uint64_t bucket_key = IndexedMap::ComputeOffsetKey(i, options.key_bits_internal_);
                offsets[bucket_key]++;
                total_capacity++;
            };
        }
        delete[] filter;

        std::cout << "Determined bucket size done." << std::endl;
        std::cout << "Determined size demands: " << total_capacity << std::endl;
        std::cout << "Estimated DB size: "
                  << ((total_capacity * (options.cell_type_/8) * (1.f / options.load_factor_)) / (1024 * 1024 * 1024.0)) << " GB"
                  << std::endl;

        bucket_bm.stop();
        bucket_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;


        std::cout << "Saving buckets.." << std::endl;
        ds::Benchmark savebuckets_bm("Saving buckets");
        savebuckets_bm.start();

        // Save bucket sizes
        SaveBucketSizes(options.db_path_ + ".buckets", offsets, options.bucket_count_);

        std::cout << "Saving buckets done" << std::endl;
        savebuckets_bm.stop();
        savebuckets_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;

        return offsets;
    }



    size_t *DetermineBucketSizes2(std::string db_path, size_t bucket_count, double load_factor, size_t threads, size_t key_bits_total, size_t key_bits_internal, std::string shape_str, const std::vector<std::string> references) {
//    size_t *DetermineBucketSizes2(BuildOptionsMG &options) {
        std::cout << "Determining bucket sizes (2)" << std::endl;

        ProgressBar bar;

        size_t k = key_bits_total/2;

        size_t *offsets = new std::size_t[bucket_count];
        std::fill_n(offsets, bucket_count, 0llu);

        const size_t shape_length = shape_str.length();
        bool *shape = ShapeUtils::GetShape(shape_str);

        size_t total_size = 0;
        size_t total_kmers = 0;

        for (auto file : references) {
            total_size += Utils::GetFileSize(file);
        }


//        const size_t block_size = (10 * 1024 * 1024);
        const size_t block_size = (16 * 1024 * 1024);
        const uint64_t max_key = (1llu << key_bits_total);
        const uint64_t key_bits_offset = key_bits_total - key_bits_internal;
        const uint64_t max_offset = (1llu << key_bits_offset);

        /* #########################################################
         * # Estimating Bucket size
         ######################################################## */

        ds::Benchmark bucket_bm("Estimating bucket size (overestimate)");
        bucket_bm.start();

        omp_set_num_threads(threads);

        ifstream* is = nullptr;

        size_t file_index = 0;

        bar.reset(total_size);
        size_t progress = 0;

//        for (auto file : options.references_) {
//            is = new ifstream(file, ios::in);

//        std::cout << "k: " << k << std::endl;
//        std::cout << "shape: " << shape << std::endl;
//        std::cout << "internal: " << key_bits_internal << std::endl;
//        std::cout << "total: " << key_bits_total << std::endl;
//        std::cout << "buckets: " << bucket_count << std::endl;
//        std::string stop;
//        std::cin >> stop;

#pragma omp parallel
        {
//                const size_t buffer_size = 1000;
//                uint64_t buffer[buffer_size];
//                size_t buffer_pos = 0;

            uint64_t key;
            uint64_t key2;
//            std::cout << "kmer iterator: " << k << ", " << shape << ", " << shape_length << std::endl;

            SimpleKmerIterator iterator(k, shape, shape_length);
            RollingKmerIterator iterator2(k, shape, shape_length);
            BufferedFastxReader reader;
            FastxRecord record;


            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                {
                    if ((!is || is->eof()) && file_index < references.size()) {
                        if (is) is->close(); delete is;
                        is = new ifstream(references[file_index++], ios::in);
                    }
                    if (is)
                        ok = reader.LoadBlock(*is, block_size);
                }
//                    ok = reader.LoadBlock(*is, block_size);
                if (!ok) break;

                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;
                    if (record.sequence.length() < shape_length) continue;
                    iterator.SetRecord(record);
                    iterator2.SetRecord(record);
                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);
                        iterator2.operator()(key2);

                        if (key2 != key) {
                            std::cout << record.header << std::endl;
                            std::cout << record.sequence.size() << std::endl;
                            std::cout << key << " != " << key2 << std::endl;

                            std::cout << iterator.GetString() << std::endl;
                            std::cout << "key simple: " << KmerUtils::ExpandShape(KmerUtils::ToString(key, k*2), shape, shape_length) << std::endl;
                            std::cout << "key roll:   " << KmerUtils::ExpandShape(KmerUtils::ToString(key2, k*2), shape, shape_length) << std::endl;


                            exit(8);
                        }

                        if (key == UINT64_MAX) {
                            continue;
                        }
                        assert(key < (1llu << (key_bits_total)));

                        uint64_t bucket_key = key >> key_bits_internal;

//                        if (key == 29978872988162) {
//                            std::cout << key << std::endl;
//                            std::cout << bucket_key << std::endl;
//                            std::cout << iterator.GetString() << std::endl;
//                            std::cout << iterator.GetString() << std::endl;
//                            std::cout << KmerUtils::ExpandShape(KmerUtils::ToString(key, k*2), shape, shape_length) << " simple" << std::endl;
//                            std::cout << KmerUtils::ExpandShape(KmerUtils::ToString(key2, k*2), shape, shape_length) << " roll" << std::endl;
//
//
//                            exit(12);
//                        }

                        assert(bucket_key < max_offset);
//                        std::cout << iterator.GetString() << " " << bucket_key << " " << key << " " << key_bits_internal << std::endl;
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

//        }
        size_t total_capacity = 0;
        for (int i = 0; i < bucket_count; i++) {
            total_capacity += offsets[i];
        }

        std::cout << "Determined bucket size done." << std::endl;
        std::cout << "Determined size demands: " << total_kmers << std::endl;
        std::cout << "Estimates DB size: "
                  << ((total_capacity * 8 * (1.f / load_factor)) / (1024 * 1024 * 1024.0)) << " GB"
                  << std::endl;

        bucket_bm.stop();
        bucket_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;


        std::cout << "Saving buckets.." << std::endl;
        ds::Benchmark savebuckets_bm("Saving buckets");
        savebuckets_bm.start();

        // Save bucket sizes
        SaveBucketSizes(db_path + ".buckets", offsets, bucket_count);

        std::cout << "Saving buckets done" << std::endl;
        savebuckets_bm.stop();
        savebuckets_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;

        delete shape;

        return offsets;
    }

    size_t *DetermineBucketSizes2(BuildOptionsMG options) {
        const std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        return DetermineBucketSizes2(options.db_path_, options.bucket_count_, options.load_factor_, options.threads_, options.key_bits_total_, options.key_bits_internal_, shape_str, options.references_);
    }

    size_t *DetermineBucketSizes2(BuildOptionsMMG options) {
        const std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        return DetermineBucketSizes2(options.db_path_, options.bucket_count_, options.load_factor_, options.threads_, options.key_bits_total_, options.key_bits_internal_, shape_str, options.references_);
    }

    size_t *DetermineBucketSizes2(BuildOptions options) {
//        const std::string shape_str = ShapeUtils::LoadShape(options.shape);
        return DetermineBucketSizes2(options.db_path_, options.bucket_count_, options.load_factor_, options.threads_, options.key_bits_total_, options.key_bits_internal_, options.shape_, options.references_);
    }


    static void RankLevelComposition(std::string path);



    void ComputeAbundanceCorrection(BuildOptionsMG &options) {
        IndexedMap* map = IndexedMap::Load(options.db_path_);
        Build::ValueExtractorMG extractor = Build::ValueExtractorMG::Load(options.entry_meta_path_);

        size_t total_size = 0;
        for (auto file : options.references_)
            total_size += Utils::GetFileSize(file);
        const size_t block_size = (1024 * 1024);

        std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        const size_t shape_length = shape_str.length();
        bool * shape = ShapeUtils::GetShape(shape_str);

        size_t progress = 0;
        omp_set_num_threads(options.threads_);

        std::unordered_map<size_t, std::unordered_map<size_t, size_t>> taxon_total_hits;
        std::unordered_map<size_t, std::unordered_map<size_t, size_t>> taxon_classified_hits;
        std::ofstream correction_out(options.working_directory_ + "correction.tsv", ios::out);

        ProgressBar bar;
        bar.reset(total_size);

        Taxonomy::IntTaxonomy taxonomy(options.working_directory_ + "/internal_taxonomy.dmp");

        ifstream is(options.references_[0], ios::in);

        std::cout << "Compute abundance correction values. " << std::endl;

#pragma omp parallel
        {
            RollingKmerIterator iterator(options.k_, shape, shape_length);
            BufferedFastxReader reader;
            FastxRecord record;

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

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

                        auto entry = map->Find(key);

                        if (entry) {

                            auto value = entry->value(map->value_bits_);

                            size_t taxid = extractor.GetTaxid(value);
                            size_t gene_id = extractor.GetGeneId(value);
                            size_t gene_pos = extractor.GetGenePos(value);

                            if (gene_id > 0 && taxonomy.Get(taxid).IsLeaf()) {
                                leaf_hit_count++;
                            }

                        } else {
                            std::cerr << "Oops. Something went wrong during database construction." << std::endl;
                            exit(9);
                        }
                        total_count++;
                    }

#pragma omp critical(taxon_hits)
                    {
                        if (!taxon_total_hits.contains(true_taxid)) {
                            std::unordered_map<size_t, size_t> inner;
                            taxon_total_hits.insert({true_taxid, inner});
                        }
                        if (!taxon_total_hits[true_taxid].contains(true_geneid)) {
                            taxon_total_hits[true_taxid][true_geneid] = 0;
                        }
                        taxon_total_hits[true_taxid][true_geneid] += total_count;

                        // Read correctly classified
                        if (!taxon_classified_hits.contains(true_taxid)) {
                            std::unordered_map<size_t, size_t> inner;
                            taxon_classified_hits.insert({true_taxid, inner});
                        }
                        if (!taxon_classified_hits[true_taxid].contains(true_geneid)) {
                            taxon_classified_hits[true_taxid][true_geneid] = 0;
                        }
                        taxon_classified_hits[true_taxid][true_geneid] += leaf_hit_count;
                    }

                } // End Read Block

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);

            } // End Reader

        } // Close OMP parallel



        for (auto& taxon : taxon_total_hits) {
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


    static void BuildMarkerGenes(BuildOptionsMG &options) {
        ds::Benchmark total_bm("Total process");
        total_bm.start();


        bool use_precomputed_buckets = Utils::exists(options.db_path_ + ".buckets");

        size_t *offsets = nullptr;
        if (!use_precomputed_buckets) {
//            size_t divisor = 1024llu * 1024 * 1024 * 8;
//            if ((1llu << options.key_bits_total_) / divisor < options.bloom_limit_gb_) {
//                offsets = DetermineBucketSizes(options);
//            } else {
////                offsets = EstimateBucketSizes(options);s
//                offsets = DetermineBucketSizes2(options);
//            }
            // This counts the number of k-mers disregarding of uniqueness
            offsets = DetermineBucketSizes2(options);
        }


        std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        const size_t shape_length = shape_str.length();
        bool *shape = ShapeUtils::GetShape(shape_str);

        ProgressBar bar;

        // First estimate bucket size
        size_t total_size = 0;

        for (auto file : options.references_) {
            total_size += Utils::GetFileSize(file);
        }

        /* #########################################################
         * # START BUILDING DATABASE
         ######################################################## */

        std::cout << "Build database" << std::endl;
        ds::Benchmark build_bm("Building database");
        build_bm.start();

        double load_factor = options.load_factor_;

        uint64_t key;
        const uint64_t max_key = (1llu << options.key_bits_total_);

        BufferedFastxReader reader;
        const size_t block_size = (10 * 1024 * 1024);

        FastxRecord record;

        IndexedMap *map;
        if (use_precomputed_buckets) {
            map = new IndexedMap(options.key_bits_internal_,
                                 options.value_bits_,
                                 options.bucket_bits_,
                                 options.db_path_ + ".buckets",
                                 load_factor);
        } else {
            map = new IndexedMap(options.key_bits_internal_,
                                 options.value_bits_,
                                 options.bucket_bits_,
                                 offsets,
                                 load_factor);
        }
        delete[] offsets;


        bar.reset(total_size);
        size_t progress = 0;
        size_t total_reads = 0;
        size_t skipped_reads = 0;

        // Compute that based on taxonomy etc....


        size_t internal_taxid_bits = 13;
        size_t marker_gene_id_bits = 9;
        size_t marker_gene_pos_bits = 14; // Used to be 22

        // Debug
        ValueExtractorMG extractor = ValueExtractorMG(internal_taxid_bits, marker_gene_id_bits, marker_gene_pos_bits);

        if (Utils::exists(options.working_directory_ + "/entry.meta")) {
            std::cout << "Use existing meta file: " << (options.working_directory_ + "/entry.meta") << std::endl;
            extractor = ValueExtractorMG::Load(options.working_directory_ + "/entry.meta");
        }


        assert(extractor.internal_taxid_bits_ + extractor.gene_id_bits_ + extractor.gene_pos_bits_ <= options.value_bits_);


        size_t max_marker_gene_pos = (1 << extractor.gene_pos_bits_) - 1;



        Taxonomy::IntTaxonomy taxonomy(options.custom_taxonomy_);

        ifstream *is = nullptr;
        size_t file_index = 0;


        omp_set_num_threads(options.threads_);
//        omp_set_num_threads(1);

#pragma omp parallel
        {
            uint64_t key = 0;
//            SimpleKmerIterator iterator(options.k_, shape, shape_length);
            RollingKmerIterator iterator(options.k_, shape, shape_length);
            BufferedFastxReader reader;
            FastxRecord record;

            size_t internal_taxid;
            size_t marker_gene_id;
            size_t marker_gene_pos;

            SimpleMarkerKmerProcessor processor(map, &taxonomy, extractor.internal_taxid_bits_, extractor.gene_id_bits_,
                                                extractor.gene_pos_bits_);

//            CountKmerProcessor processor(map, &taxonomy, extractor.internal_taxid_bits_, extractor.gene_id_bits_,
//                                                extractor.gene_pos_bits_);



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
//                    ok = reader.LoadBlock(*is, block_size);
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

                    iterator.SetRecord(record);
                    SimpleMarkerKmerProcessor::ParseHeader(record.header, internal_taxid, marker_gene_id);

//                    // If target rank is specified map internal_taxid onto that rank
//                    if (!options.target_rank_.empty() &&
//                        strcmp(options.target_rank_.c_str(), taxonomy.Get(internal_taxid).rank.c_str()) != 0) {
//                        Taxonomy::IntNode node = taxonomy.Get(taxonomy.Get(internal_taxid).parent_id);
//                        while (strcmp(node.rank.c_str(), options.target_rank_.c_str()) != 0) {
//                            node = taxonomy.Get(taxonomy.Get(internal_taxid).parent_id);
//                        }
//                        internal_taxid = node.id;
//                    }


                    if (!taxonomy.IsLeaf(internal_taxid)) {
                        std::cout << record.id << std::endl;
                        std::cout << "iid: " << internal_taxid << std::endl;
                        std::cout << "geneid: " << marker_gene_id << std::endl;
                        std::cout << "genepos: " << marker_gene_pos << std::endl;
                    }

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

//                        if (key == UINT64_MAX) {
//                            continue;
//                        }
                        marker_gene_pos = iterator.GetPos();

                        if (marker_gene_pos > max_marker_gene_pos) {
                            std::cout << "Marker gene too long, just stop here: " << marker_gene_pos << " max: " << max_marker_gene_pos << std::endl;
                            break;
                        }

#pragma omp critical(put_operation)
                        {
                            processor.operator()(key, internal_taxid, marker_gene_id, marker_gene_pos);
                        }
                    }
                }

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);
            }
        }
        delete[] shape;

        std::cout << "Building database done." << std::endl;
        build_bm.stop();
        build_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;

        ds::Benchmark save_bm("Saving db");
        save_bm.start();
        map->Save(options.db_path_);
        map->SaveBucketSizes(options.db_path_ + ".buckets");
        save_bm.stop();

        // Save values for internal_taxid_bits etc...
        extractor.Save(options.entry_meta_path_);

        std::cout << "Saving database done." << std::endl;
        save_bm.printResults();

        std::cout << std::string(60, '-') << std::endl;

        if (options.eval_) {

//            std::cout << "Look at first 100 buckets." << std::endl;
//            map->print(0, 100, false);
//            std::cout << std::string(60, '-') << std::endl;

            std::cout << "Check buckets." << std::endl;
            map->CheckBuckets();
            std::cout << std::string(60, '-') << std::endl;

            map->PrintMeta();

            std::cout << std::string(60, '-') << std::endl;

            // Get rank-level distribution in DB
            std::cout << "Rank level composition of DB\n" << std::endl;

            RankLevelComposition(options.working_directory_);
            std::cout << std::string(60, '-') << std::endl;
        }



//        bool count_occurrences = true;
//        size_t value_mask = ((1lu << (extractor.gene_id_bits_ + extractor.gene_pos_bits_)) - 1);
//        if (count_occurrences) {
//            for (auto& cell : *map) {
//                auto value = cell.value(map->ValueBits());
//                if (cell.empty() || (value & value_mask) < 50) {
//                    continue;
//                }
//                size_t taxid = extractor.GetTaxid(value);
//                std::cout << (cell.value(map->ValueBits()) & value_mask) << "  #children: " << taxonomy.GetSubtreeSize(taxid) << " " << taxonomy.Get(taxid).scientific_name << std::endl;
//            }
//        }
        ComputeAbundanceCorrection(options);

        total_bm.stop();
        total_bm.printResults();

        delete map;
    }


//#######################################################################################################
// asdasdasdasdasdasd
//#######################################################################################################

    template <typename T>
    static void BuildPreMap(BuildOptions &options) {
        ds::Benchmark total_bm("Total process");
        total_bm.start();


        bool use_precomputed_buckets = Utils::exists(options.db_path_ + ".buckets");

        bool multimap = false;

        size_t *offsets = nullptr;
        if (!use_precomputed_buckets) {
//            offsets = DetermineBucketSizes2(options);
            offsets = DetermineBucketSizes(options);
        } else {
            auto buckets = LoadBucketSizes(options.db_path_ + ".buckets");
            assert(buckets.second == options.bucket_count_);
            offsets = buckets.first;
        }


//        std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        std::string shape_str = options.shape_;
        const size_t shape_length = shape_str.length();
        bool *shape = ShapeUtils::GetShape(shape_str);

        ProgressBar bar;

        // First estimate bucket size
        size_t total_size = 0;

        for (auto file : options.references_) {
            total_size += Utils::GetFileSize(file);
        }

        /* #########################################################
         * # START BUILDING DATABASE
         ######################################################## */

        std::cout << "Build pre-database" << std::endl;
        ds::Benchmark build_bm("Building pre-database");
        build_bm.start();

        double load_factor = options.load_factor_;

        uint64_t key;
        const uint64_t max_key = (1llu << options.key_bits_total_);

        BufferedFastxReader reader;
        const size_t block_size = (10 * 1024 * 1024);

        FastxRecord record;

//        IndexedMap *map;
//        if (use_precomputed_buckets) {
//            map = new MiniMap::MiniMap<T>(options.key_bits_internal_,
//                                 options.value_bits_,
//                                 options.bucket_bits_,
//                                 options.db_path_ + ".buckets",
//                                 load_factor);
//        }
//        else {
//            map = new MiniMap::MiniMap<T>(options.key_bits_internal_,
//                                 options.value_bits_,
//                                 options.bucket_bits_,
//                                 offsets,
//                                 load_factor);
//        }
        auto map = new MiniMap::MiniMap<T>(
                options.key_bits_internal_,
                options.value_bits_,
                options.bucket_bits_,
                offsets,
                load_factor);
        delete[] offsets;


        bar.reset(total_size);
        size_t progress = 0;
        size_t total_reads = 0;
        size_t skipped_reads = 0;

        // Compute that based on taxonomy etc....
        size_t internal_taxid_bits = 13;


        Taxonomy::IntTaxonomy taxonomy(options.custom_taxonomy_path_);

        ifstream *is = nullptr;
        size_t file_index = 0;

        omp_set_num_threads(options.threads_);

        size_t count_n = 0;

#pragma omp parallel
        {
            uint64_t key = 0;
            RollingKmerIterator iterator(options.k_, shape, shape_length);
            BufferedFastxReader reader;
            FastxRecord record;

            size_t internal_taxid;
            size_t internal_taxid_bits = options.value_bits_;


            SimpleKmerProcessor<T> processor(map, &taxonomy, internal_taxid_bits);


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

//                    std::cout << record.id << std::endl;

                    ++total_reads;

                    // Skip read if too short!
                    if (record.sequence.length() < shape_length) {
                        ++skipped_reads;
                        continue;
                    }
                    size_t marker_gene_id = 0;
                    iterator.SetRecord(record);
                    SimpleMarkerKmerProcessor::ParseHeader(record.header, internal_taxid, marker_gene_id);

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);




                        if (key == UINT64_MAX) {
#pragma omp atomic
                            count_n++;
                            continue;
                        }
                        if (key >= (1llu << (options.bucket_bits_ + options.key_bits_internal_))) {
                            std::cout << std::endl;
                            std::cout << iterator.GetString() << std::endl;
                            std::cout << key << std::endl;
                            std::cout << iterator.K() << std::endl;
                            std::cout << iterator.ShapeSize() << std::endl;
                        }

                        assert(key < (1llu << (options.bucket_bits_ + options.key_bits_internal_)));


#pragma omp critical(put_operation)
                        {
                            processor.operator()(key, internal_taxid);
                        }
                    }
                }

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);
            }
        }
        delete[] shape;

        std::cout << "Building database done." << std::endl;
        build_bm.stop();
        build_bm.printResults();
        std::cout << std::string(60, '-') << std::endl;

        ds::Benchmark save_bm("Saving db");
        save_bm.start();
        map->Save(options.db_path_);
//        map->SaveBucketSizes(options.db_path_ + ".buckets");
        save_bm.stop();


        std::cout << "Saving database done." << std::endl;
        save_bm.printResults();

        std::cout << std::string(60, '-') << std::endl;

        if (true) {
            std::cout << "Check buckets." << std::endl;
            map->CheckBuckets();
            std::cout << std::string(60, '-') << std::endl;

            map->PrintMeta();

            std::cout << std::string(60, '-') << std::endl;

            // Get rank-level distribution in DB
            std::cout << "Rank level composition of DB\n" << std::endl;

            RankLevelCompositionMiniMap<T>(options.db_dir_);
            std::cout << std::string(60, '-') << std::endl;
        }

        std::cout << "kmers containing Ns: " << count_n << std::endl;
        std::cout << "kmers containing Ns(%): " << (float)count_n/map->Size() << std::endl;

        total_bm.stop();
        total_bm.printResults();

        delete map;
    }


    /* New Building Index Routine
     * For Multimap that implements binary search and includes pack sizes */
    static void BuildIndex(BuildOptionsMMG &options) {
        ds::Benchmark total_bm("Total process");
        total_bm.start();

        // If buckets have once been computed use existing buckets from file.
        bool use_precomputed_buckets = Utils::exists(options.db_path_ + ".buckets");

        // Object holding the bucket sizes
        size_t *offsets = nullptr;
//        if (!use_precomputed_buckets) {
//            // This counts the number of k-mers disregarding of uniqueness
//
//        }
        offsets = DetermineBucketSizes2(options);


        std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        const size_t shape_length = shape_str.length();
        bool *shape = ShapeUtils::GetShape(shape_str);

        ProgressBar bar;

        // First estimate bucket size
        size_t total_size = 0;

        for (auto file : options.references_)
            total_size += Utils::GetFileSize(file);

        /* #########################################################
         * # START BUILDING DATABASE
         ######################################################## */

        std::cout << "Build database" << std::endl;
        ds::Benchmark build_bm("Building database");
        build_bm.start();

        // Load factor to determine size of HashMap
        // Normal is around 70-80%
        // Does not matter for small buckets
        double load_factor = options.load_factor_;

        uint64_t key;
        const uint64_t max_key = (1llu << options.key_bits_total_);

        // FastaReader:
        // -Block size amount of bytes to read and keep in buffer
        const size_t block_size = (10 * 1024 * 1024);

        // MultiMap object
        MultiMap::MultiMap *map;
//        if (use_precomputed_buckets) {
////            map = new MultiMap::MultiMap(options.key_bits_internal_,
////                                 options.value_bits_,
////                                 options.bucket_bits_,
////                                 options.db_path_ + ".buckets",
////                                 load_factor,
////                                 multimap);
//        } else {
//            map = new MultiMap::MultiMap(options.key_bits_internal_,
//                                 options.value_bits_,
//                                 options.bucket_bits_,
//                                 offsets,
//                                 load_factor);
//        }
        map = new MultiMap::MultiMap(options.key_bits_internal_,
                                     options.value_bits_,
                                     options.bucket_bits_,
                                     offsets,
                                     load_factor);
        delete[] offsets;


        // some counts to keep track of
        size_t progress = 0;
        size_t total_reads = 0;
        size_t skipped_reads = 0;

        // Compute that based on taxonomy etc....
        size_t internal_taxid_bits = 13;
        size_t marker_gene_id_bits = 9;
        size_t marker_gene_pos_bits = 14; // Used to be 22

        // Debug
        ValueExtractorMG extractor = ValueExtractorMG(internal_taxid_bits, marker_gene_id_bits, marker_gene_pos_bits);

        std::cout << "bits for key: " << (64 - map->ValueBits() - 2 + map->OffsetBits()) << std::endl;
        std::cout << "bits for key needed: " << options.k_ * 2 << std::endl;

        // Sanity check
        if ( (64 - map->ValueBits() - 2 + map->OffsetBits()) > options.k_ * 2 ) {
            std::cout << "Abort. " << (64 - map->ValueBits() - 2 + map->OffsetBits()) << " > " << options.k_ * 2 << std::endl;
            exit(9);
        }

        // If meta data exists use this (currently default..)
        if (Utils::exists(options.working_directory_ + "/entry.meta")) {
            std::cout << "Use existing meta file: " << (options.working_directory_ + "/entry.meta") << std::endl;
            extractor = ValueExtractorMG::Load(options.working_directory_ + "/entry.meta");
        }

        assert(extractor.internal_taxid_bits_ + extractor.gene_id_bits_ + extractor.gene_pos_bits_ <= options.value_bits_);

        // Determine max marker gene size according to bits reserved for marker gene position
        // Cut off above that size.
        size_t max_marker_gene_pos = (1 << extractor.gene_pos_bits_) - 1;

        // Load taxonomy for this
        // not necessary for multimap
        Taxonomy::IntTaxonomy taxonomy(options.custom_taxonomy_);

        ifstream *is = nullptr;
        size_t file_index = 0;

        // Threads for multithreaded opteration
        omp_set_num_threads(options.threads_);

        // Set progress bar to combined size of files
        bar.reset(total_size);

#pragma omp parallel
        {
            uint64_t key = 0;

            // Fast kmer iterator for shapes
            // Iterates all k-mers for given sequences
            RollingKmerIterator iterator(options.k_, shape, shape_length);



            /* Reader and object holding the next sequence item
             * >Header
             * ACGTCGTCGC
             */
            BufferedFastxReader reader;
            FastxRecord record;

            // Variables holding the taxid, marker gene id and pos
            // Will be passed to functions as reference
            size_t internal_taxid;
            size_t marker_gene_id;
            size_t marker_gene_pos;

            // Processor puts kmers into map.
            MultiMapMarkerKmerProcessor processor(map, extractor.internal_taxid_bits_,
                                                  extractor.gene_id_bits_,
                                                  extractor.gene_pos_bits_);

            /* ##########################################################################
             * Iterate all reference sequences (Marker genes)
             * ########################################################################## */
            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                {
                    // Switch to next marker gene reference file if old done processing
                    // This section is need if many multiple input files exist
                    // that not all threads have to wait for one reference to finish,
                    // but that idleing threads can already move on to the next reference file
                    // While the old ones still working on the last bit of the last file
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

                    // pass on record object to iterator.
                    iterator.SetRecord(record);

                    // Extract Information from header
                    SimpleMarkerKmerProcessor::ParseHeader(record.header, internal_taxid, marker_gene_id);

                    // Sanity check. Internal Taxid has to be leaf
                    // Should be omitted in final production code
                    if (!taxonomy.IsLeaf(internal_taxid)) {
                        errx(EX_DATAERR,
                             "Internal taxid (%llu) must be, but is not a leaf in the given taxonomy.",
                             (unsigned long long) internal_taxid);
                    }

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

                        if (key == UINT64_MAX) {
                            // Add counter for keys...
                            continue;
                        }

                        marker_gene_pos = iterator.GetPos();

                        // Stop if position is larger than max marker gene pos
                        if (marker_gene_pos > max_marker_gene_pos) {
                            // Add counter for marker genes that are too short
                            break;
                        }

#pragma omp critical(put_operation)
                        {
                            processor.operator()(key, internal_taxid, marker_gene_id, marker_gene_pos);
                        }

                    }
                }

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);
            }
        }
        delete[] shape;

        std::cout << "Inserting into database done." << std::endl;
        build_bm.stop();
        build_bm.printResults();


        std::cout << std::string(60, '-') << std::endl;
        std::cout << "Annotate entries, calculate block sizes and sort multi-key blocks." << std::endl;
        map->GeneratePackSizes();

        // Save DB dialog
        std::cout << std::string(60, '-') << std::endl;
        ds::Benchmark save_bm("Saving db");
        save_bm.start();
        map->Save(options.db_path_);
        //map->SaveBucketSizes(options.db_path_ + ".buckets");
        save_bm.stop();

        // Save values for internal_taxid_bits etc...
        extractor.Save(options.entry_meta_path_);

        std::cout << "Saving database done." << std::endl;
        save_bm.printResults();

        std::cout << std::string(60, '-') << std::endl;

        // Check buckets (What?)
        std::cout << "Check buckets." << std::endl;
        map->CheckBuckets();

        // Map meta information
        std::cout << std::string(60, '-') << std::endl;
        map->PrintMeta();

        // Get rank-level distribution in DB
        std::cout << std::string(60, '-') << std::endl;
        std::cout << "Rank level composition of DB\n" << std::endl;
        RankLevelComposition(options.working_directory_);

        // Display total build time
        std::cout << std::string(60, '-') << std::endl;
        total_bm.stop();
        total_bm.printResults();

        delete map;
    }


    static void QueryDB(std::string path) {
        IndexedMap *map = IndexedMap::Load(path + "/database.bin");

        Taxonomy::IntTaxonomy taxonomy(path + "/internal_taxonomy.dmp");


        size_t size = map->Capacity();
        Cell *cells = map->Map();

        for (size_t i = 0; i < 1000; i++) {
            auto cell = cells[i];

            if (!cell.empty()) {
                auto key = cell.key(map->value_bits_);
                auto value = cell.value(map->value_bits_);

            }
        }

        delete map;
    }

    template <typename T>
    static void RankLevelCompositionMiniMap(std::string path) {
        MiniMap::MiniMap<T> *map = MiniMap::MiniMap<T>::Load(path + "/database.bin");
        Taxonomy::IntTaxonomy taxonomy(path + "/internal_taxonomy.dmp");


    /* New Building Index Routine
     * For Multimap that implements binary search and includes pack sizes */
    static void BuildIndex(BuildOptionsMMG &options) {
        ds::Benchmark total_bm("Total process");
        total_bm.start();

        // If buckets have once been computed use existing buckets from file.
        bool use_precomputed_buckets = Utils::exists(options.db_path_ + ".buckets");

        // Object holding the bucket sizes
        size_t *offsets = nullptr;
//        if (!use_precomputed_buckets) {
//            // This counts the number of k-mers disregarding of uniqueness
//
//        }
        offsets = DetermineBucketSizes2(options);


        std::string shape_str = ShapeUtils::LoadShape(options.shape_path_);
        const size_t shape_length = shape_str.length();
        bool *shape = ShapeUtils::GetShape(shape_str);

        ProgressBar bar;

        // First estimate bucket size
        size_t total_size = 0;

        for (auto file : options.references_)
            total_size += Utils::GetFileSize(file);

        /* #########################################################
         * # START BUILDING DATABASE
         ######################################################## */

        std::cout << "Build database" << std::endl;
        ds::Benchmark build_bm("Building database");
        build_bm.start();

        // Load factor to determine size of HashMap
        // Normal is around 70-80%
        // Does not matter for small buckets
        double load_factor = options.load_factor_;

        uint64_t key;
        const uint64_t max_key = (1llu << options.key_bits_total_);

        // FastaReader:
        // -Block size amount of bytes to read and keep in buffer
        const size_t block_size = (10 * 1024 * 1024);

        // MultiMap object
        MultiMap::MultiMap *map;
//        if (use_precomputed_buckets) {
////            map = new MultiMap::MultiMap(options.key_bits_internal_,
////                                 options.value_bits_,
////                                 options.bucket_bits_,
////                                 options.db_path_ + ".buckets",
////                                 load_factor,
////                                 multimap);
//        } else {
//            map = new MultiMap::MultiMap(options.key_bits_internal_,
//                                 options.value_bits_,
//                                 options.bucket_bits_,
//                                 offsets,
//                                 load_factor);
//        }
        map = new MultiMap::MultiMap(options.key_bits_internal_,
                                     options.value_bits_,
                                     options.bucket_bits_,
                                     offsets,
                                     load_factor);
        delete[] offsets;


        // some counts to keep track of
        size_t progress = 0;
        size_t total_reads = 0;
        size_t skipped_reads = 0;

        // Compute that based on taxonomy etc....
        size_t internal_taxid_bits = 13;
        size_t marker_gene_id_bits = 9;
        size_t marker_gene_pos_bits = 14; // Used to be 22

        // Debug
        ValueExtractorMG extractor = ValueExtractorMG(internal_taxid_bits, marker_gene_id_bits, marker_gene_pos_bits);

        std::cout << "bits for key: " << (64 - map->ValueBits() - 2 + map->OffsetBits()) << std::endl;
        std::cout << "bits for key needed: " << options.k_ * 2 << std::endl;

        // Sanity check
        if ( (64 - map->ValueBits() - 2 + map->OffsetBits()) > options.k_ * 2 ) {
            std::cout << "Abort. " << (64 - map->ValueBits() - 2 + map->OffsetBits()) << " > " << options.k_ * 2 << std::endl;
            exit(9);
        }

        // If meta data exists use this (currently default..)
        if (Utils::exists(options.working_directory_ + "/entry.meta")) {
            std::cout << "Use existing meta file: " << (options.working_directory_ + "/entry.meta") << std::endl;
            extractor = ValueExtractorMG::Load(options.working_directory_ + "/entry.meta");
        }

        assert(extractor.internal_taxid_bits_ + extractor.gene_id_bits_ + extractor.gene_pos_bits_ <= options.value_bits_);

        // Determine max marker gene size according to bits reserved for marker gene position
        // Cut off above that size.
        size_t max_marker_gene_pos = (1 << extractor.gene_pos_bits_) - 1;

        // Load taxonomy for this
        // not necessary for multimap
        Taxonomy::IntTaxonomy taxonomy(options.custom_taxonomy_);

        ifstream *is = nullptr;
        size_t file_index = 0;

        // Threads for multithreaded opteration
        omp_set_num_threads(options.threads_);

        // Set progress bar to combined size of files
        bar.reset(total_size);

#pragma omp parallel
        {
            uint64_t key = 0;

            // Fast kmer iterator for shapes
            // Iterates all k-mers for given sequences
            RollingKmerIterator iterator(options.k_, shape, shape_length);



            /* Reader and object holding the next sequence item
             * >Header
             * ACGTCGTCGC
             */
            BufferedFastxReader reader;
            FastxRecord record;

            // Variables holding the taxid, marker gene id and pos
            // Will be passed to functions as reference
            size_t internal_taxid;
            size_t marker_gene_id;
            size_t marker_gene_pos;

            // Processor puts kmers into map.
            MultiMapMarkerKmerProcessor processor(map, extractor.internal_taxid_bits_,
                                                  extractor.gene_id_bits_,
                                                  extractor.gene_pos_bits_);

            /* ##########################################################################
             * Iterate all reference sequences (Marker genes)
             * ########################################################################## */
            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                {
                    // Switch to next marker gene reference file if old done processing
                    // This section is need if many multiple input files exist
                    // that not all threads have to wait for one reference to finish,
                    // but that idleing threads can already move on to the next reference file
                    // While the old ones still working on the last bit of the last file
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

                    // pass on record object to iterator.
                    iterator.SetRecord(record);

                    // Extract Information from header
                    SimpleMarkerKmerProcessor::ParseHeader(record.header, internal_taxid, marker_gene_id);

                    // Sanity check. Internal Taxid has to be leaf
                    // Should be omitted in final production code
                    if (!taxonomy.IsLeaf(internal_taxid)) {
                        errx(EX_DATAERR,
                             "Internal taxid (%llu) must be, but is not a leaf in the given taxonomy.",
                             (unsigned long long) internal_taxid);
                    }

                    // Iterate over all k-mers
                    while (iterator.HasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);

                        if (key == UINT64_MAX) {
                            // Add counter for keys...
                            continue;
                        }

                        marker_gene_pos = iterator.GetPos();

                        // Stop if position is larger than max marker gene pos
                        if (marker_gene_pos > max_marker_gene_pos) {
                            // Add counter for marker genes that are too short
                            break;
                        }

#pragma omp critical(put_operation)
                        {
                            processor.operator()(key, internal_taxid, marker_gene_id, marker_gene_pos);
                        }

                    }
                }

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);
            }
        }
        delete[] shape;

        std::cout << "Inserting into database done." << std::endl;
        build_bm.stop();
        build_bm.printResults();


        std::cout << std::string(60, '-') << std::endl;
        std::cout << "Annotate entries, calculate block sizes and sort multi-key blocks." << std::endl;
        map->GeneratePackSizes();

        // Save DB dialog
        std::cout << std::string(60, '-') << std::endl;
        ds::Benchmark save_bm("Saving db");
        save_bm.start();
        map->Save(options.db_path_);
        //map->SaveBucketSizes(options.db_path_ + ".buckets");
        save_bm.stop();

        // Save values for internal_taxid_bits etc...
        extractor.Save(options.entry_meta_path_);

        std::cout << "Saving database done." << std::endl;
        save_bm.printResults();

        std::cout << std::string(60, '-') << std::endl;

        // Check buckets (What?)
        std::cout << "Check buckets." << std::endl;
        map->CheckBuckets();

        // Map meta information
        std::cout << std::string(60, '-') << std::endl;
        map->PrintMeta();

        // Get rank-level distribution in DB
        std::cout << std::string(60, '-') << std::endl;
        std::cout << "Rank level composition of DB\n" << std::endl;
        RankLevelComposition(options.working_directory_);

        // Display total build time
        std::cout << std::string(60, '-') << std::endl;
        total_bm.stop();
        total_bm.printResults();

        delete map;
    }


    static void QueryDB(std::string path) {
        IndexedMap *map = IndexedMap::Load(path + "/database.bin");

        Taxonomy::IntTaxonomy taxonomy(path + "/internal_taxonomy.dmp");


        size_t size = map->Capacity();
        Cell *cells = map->Map();

        for (size_t i = 0; i < 1000; i++) {
            auto cell = cells[i];

            if (!cell.empty()) {
                auto key = cell.key(map->value_bits_);
                auto value = cell.value(map->value_bits_);

            }
        }

        delete map;
    }

    template <typename T>
    static void RankLevelCompositionMiniMap(std::string path) {
        MiniMap::MiniMap<T> *map = MiniMap::MiniMap<T>::Load(path + "/database.bin");
        Taxonomy::IntTaxonomy taxonomy(path + "/internal_taxonomy.dmp");


        size_t size = map->Capacity();
        size_t size = map->Capacity();
        auto *cells = map->Map();

        tsl::robin_map<std::string, size_t> rank_counts;


        size_t total = 0;
        size_t leaf_count = 0;

        for (size_t i = 0; i < size; i++) {
            auto cell = cells[i];

            if (!cell.empty()) {
                auto key = cell.key(map->ValueBits());
                auto value = cell.value(map->ValueBits());

                total++;


                std::string rank = taxonomy.Get(value).rank;

                assert(!rank.empty());

                if (!rank_counts.contains(rank)) {
                    rank_counts.insert({rank, 0});
                }

                if (taxonomy.Get(value).IsLeaf()) {
                    leaf_count++;
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
        std::cout << "total k-mers: " << total << std::endl;


        delete map;
    }

    static void DBCompositionMM(std::string path) {
        MultiMap::MultiMap *map = MultiMap::MultiMap::Load(path + "/database.bin");
        ValueExtractorMG extractor = ValueExtractorMG::Load(path + "/entry.meta");

        size_t size = map->Capacity();
        auto *cells = map->Map();
        auto offset_ = map->Offset();

        size_t offset_index = 0;
        size_t bucket_start = 0;
        size_t bucket_end = 0;
        size_t bucket_size = 0;

        size_t total = 0;
        size_t leaf_count = 0;

        std::vector<size_t> kmer_pack_sizes(1000);


        for (size_t i = 0; i < size; i++) {

            while (bucket_end <= i) {
                offset_index++;
                bucket_start = bucket_end;
                bucket_end = offset_[offset_index];
                bucket_size = bucket_end - bucket_start;
            }

            auto cell = cells[i];

            if (cell.get_flag() == MultiMap::UNIQUE) {
                kmer_pack_sizes[1] += 1;
            } else if (cell.get_flag() == MultiMap::MULTI_KEY_START) {
                size_t next_idx = i + 1 == bucket_end ? bucket_start : i + 1;
                auto next = cells + next_idx;
                size_t pack_size = next->key(map->ValueBits(), map->KeyBits()) + 1;

                if (pack_size >= kmer_pack_sizes.size()) {
                    kmer_pack_sizes.resize(kmer_pack_sizes.size() * 2);
                    std::cout << "resize kmer_pack_sizes: " << kmer_pack_sizes.size() << std::endl;
                }

                kmer_pack_sizes[pack_size]++;
            }


        }

        for (auto i = 0; i < kmer_pack_sizes.size(); i++) {
            std::cout << i << '\t' << kmer_pack_sizes[i] << std::endl;
        }

        delete map;
    }


    static void RankLevelComposition(std::string path) {
        IndexedMap *map = IndexedMap::Load(path + "/database.bin");
        Taxonomy::IntTaxonomy taxonomy(path + "/internal_taxonomy.dmp");
        ValueExtractorMG extractor = ValueExtractorMG::Load(path + "/entry.meta");


        size_t size = map->Capacity();
        Cell *cells = map->Map();

        tsl::robin_map<std::string, size_t> rank_counts;


        size_t total = 0;
        size_t leaf_count = 0;
        size_t non_leaf_with_gid = 0;

        for (size_t i = 0; i < size; i++) {
            auto cell = cells[i];

            if (!cell.empty()) {
                auto key = cell.key(map->value_bits_);
                auto value = cell.value(map->value_bits_);

                total++;

                size_t internal_taxid = extractor.GetTaxid(value);
                size_t marker_gene_id = extractor.GetGeneId(value);
                size_t marker_gene_pos = extractor.GetGenePos(value);

                std::string rank = taxonomy.Get(internal_taxid).rank;

                if (rank.empty()) {
                    std::cout << "Empty:" << std::endl;
                    std::cout << "key: " << key << std::endl;
                    std::cout << "value: " << value << std::endl;
                    std::cout << "internal_taxid: " << internal_taxid << std::endl;
                }
                if (!taxonomy.Get(internal_taxid).IsLeaf() && marker_gene_id > 0) {
                    non_leaf_with_gid++;
                    printf("%lu\t%lu\t%lu\t%s\n", internal_taxid, marker_gene_id, marker_gene_pos, taxonomy.Get(internal_taxid).scientific_name.c_str());

                }

                assert(!rank.empty());

                if (!rank_counts.contains(rank)) {
                    rank_counts.insert({rank, 0});
                }

                if (taxonomy.Get(internal_taxid).IsLeaf()) {
                    leaf_count++;
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
        std::cout << "non leaf with gene id information: " << non_leaf_with_gid << std::endl;
        std::cout << "total k-mers: " << total << std::endl;


        delete map;
    }
}