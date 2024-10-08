////
//// Created by fritsche on 26/05/2021.
////
//
//#pragma once
//
//#include <cstddef>
//#include <string>
//#include <Benchmark.h>
//#include <compact_map.h>
//#include <SimpleKmerIterator.h>
//#include <FastxReader.h>
//#include <ShapeUtils.h>
//#include <omp.h>
//#include <TaxonomyNew.h>
//#include <build_scripts.h>
//#include <ReadClassifier.h>
//#include <IOHandler.h>
//#include <BufferedOutput.h>
//
//#define CLASSIFY_OMP 1
//
//#define KNOWNREF 1
//#define KNOWNREF_DEBUG 0
//#define WRITEOUT_FP 0
//
//namespace Classify {
//    static const size_t LEAF_HIT_THRESHOLD = 2;
//    static const size_t BEST_TOTAL_THRESHOLD = 2;
//
//    struct ClassifyOptions {
//    public:
//        const size_t threads_;
//
//        const std::string output_dir_;
//        const std::string db_dir_;
//        const std::string meta_path_;
//        const std::string db_path_;
//        const std::string custom_taxonomy_path_;
//        const bool sensitive_snps_;
//
//        const std::vector<std::string> reads_;
//
//        const bool output_unclassified_ = false;
//
//        ClassifyOptions(size_t threads,
//                        std::string db_dir,
//                        std::string output_dir,
//                        bool output_unclassified,
//                        bool sensitive_snps,
//                        std::vector<std::string> &reads) :
//                threads_(threads),
//                output_dir_(output_dir),
//                db_dir_(db_dir),
//                output_unclassified_(output_unclassified),
//                sensitive_snps_(sensitive_snps),
//                db_path_(db_dir + "/" + Build::DATABASE),
//                custom_taxonomy_path_(db_dir + "/" + Build::TAXONOMY_CUSTOM),
//                meta_path_(db_dir + "/" + Build::ENTRY_META),
//                reads_(reads) {};
//
//        bool Check() {
//            if (!Utils::exists(db_path_)) {
//                errx(EX_IOERR,
//                     "Database file %s does not exist.", db_path_.c_str());
//            }
//            if (!Utils::exists(custom_taxonomy_path_)) {
//                errx(EX_IOERR,
//                     "Taxonomy names file %s does not exist.", custom_taxonomy_path_.c_str());
//            }
//            for (auto file : reads_) {
//                if (!Utils::exists(file)) {
//                    errx(EX_IOERR,
//                         "Reference file %s does not exist.", file.c_str());
//                }
//            }
//
//            return true;
//        }
//    };
//
//
//    struct ClassifyOptionsMM {
//    public:
//        const size_t threads_;
//
//        const std::string output_dir_;
//        const std::string db_dir_;
//        const std::string meta_path_;
//        const std::string db_path_;
//        const std::string custom_taxonomy_path_;
//
//        const std::vector<std::string> reads_;
//
//        const bool output_unclassified_ = false;
//
//        ClassifyOptionsMM(size_t threads,
//                        std::string db_dir,
//                        std::string output_dir,
//                        bool output_unclassified,
//                        std::vector<std::string> &reads) :
//                threads_(threads),
//                output_dir_(output_dir),
//                db_dir_(db_dir),
//                output_unclassified_(output_unclassified),
//                db_path_(db_dir + "/" + Build::DATABASE),
//                custom_taxonomy_path_(db_dir + "/" + Build::TAXONOMY_CUSTOM),
//                meta_path_(db_dir + "/" + Build::ENTRY_META),
//                reads_(reads) {};
//
//        bool Check() {
//            if (!Utils::exists(db_path_)) {
//                errx(EX_IOERR,
//                     "Database file %s does not exist.", db_path_.c_str());
//            }
//            if (!Utils::exists(custom_taxonomy_path_)) {
//                errx(EX_IOERR,
//                     "Taxonomy names file %s does not exist.", custom_taxonomy_path_.c_str());
//            }
//            for (auto file : reads_) {
//                if (!Utils::exists(file)) {
//                    errx(EX_IOERR,
//                         "Reference file %s does not exist.", file.c_str());
//                }
//            }
//
//            return true;
//        }
//    };
//
//    static void QueryDB(std::string db, std::vector<std::string> queries, bool rank_composition) {
//        IndexedMap* map = IndexedMap::Load(db + "/database.bin");
//        Taxonomy::IntTaxonomy taxonomy(db + "/internal_taxonomy.dmp");
//        Build::ValueExtractorMG extractor = Build::ValueExtractorMG::Load(db + "/entry.meta");
//
//        size_t k { (map->key_bits_ + map->offset_bits_) / 2 };
//        std::string shape_str = ShapeUtils::LoadShape(db + "/shape.txt");
//        const size_t shape_length = shape_str.length();
//        bool * shape = ShapeUtils::GetShape(shape_str);
//
//        RollingKmerIterator iterator(k, shape, shape_length);
//        BufferedFastxReader reader;
//        FastxRecord record;
//
//        const size_t block_size = (1024 * 1024);
//
//
//
//        class Ranks {
//        public:
//            std::unordered_map<std::string, uint64_t> ranks;
//            std::string name_;
//
//            Ranks(std::string name) : name_(name) {
//                ranks.insert( { "no rank", 0 } );
//                ranks.insert( { "domain", 0 } );
//                ranks.insert( { "phylum", 0 } );
//                ranks.insert( { "class", 0 } );
//                ranks.insert( { "order", 0 } );
//                ranks.insert( { "family", 0 } );
//                ranks.insert( { "genus", 0 } );
//                ranks.insert( { "species", 0 } );
//            }
//
//            Ranks() {
//                ranks.insert( { "no rank", 0 } );
//                ranks.insert( { "domain", 0 } );
//                ranks.insert( { "phylum", 0 } );
//                ranks.insert( { "class", 0 } );
//                ranks.insert( { "order", 0 } );
//                ranks.insert( { "family", 0 } );
//                ranks.insert( { "genus", 0 } );
//                ranks.insert( { "species", 0 } );
//            }
//
//            std::string GetString() {
//                std::string line = "";
//
//                line += to_string(ranks["species"]);
//                line += '\t';
//                line += to_string(ranks["genus"]);
//                line += '\t';
//                line += to_string(ranks["family"]);
//                line += '\t';
//                line += to_string(ranks["order"]);
//                line += '\t';
//                line += to_string(ranks["class"]);
//                line += '\t';
//                line += to_string(ranks["phylum"]);
//                line += '\t';
//                line += to_string(ranks["domain"]);
//                line += '\t';
//                line += to_string(ranks["no rank"]);
//
//                return line;
//            }
//        };
//
//        std::unordered_map<std::string, Ranks> rankmap;
//
//
//        uint64_t key;
//        for (auto& path : queries) {
//            ifstream is(path);
//
//            while (true) {
//                bool ok = false;
//
//                ok = reader.LoadBlock(is, block_size);
//                if (!ok) break;
//
//
//                // Read records from datablock
//                while (true) {
//                    auto valid_fragment = reader.NextSequence(record);
//                    if (!valid_fragment) break;
//
//                    iterator.SetRecord(record);
//
//                    // Iterate over all k-mers
//                    int pos = 0;
//                    while (iterator.HasNext()) {
//                        iterator.operator()(key);
//
//                        if (key == UINT64_MAX) {
//                            continue;
//                        }
//
//                        auto entry = map->Find(key);
//
//                        if (entry) {
//
//                            auto value = entry->value(map->value_bits_);
//
//                            size_t taxid = extractor.GetTaxid(value);
//                            size_t gene_id = extractor.GetGeneId(value);
//                            size_t gene_pos = extractor.GetGenePos(value);
//
//                            auto& node = taxonomy.Get(taxid);
//
//                            if (!rank_composition) {
//                                std::cout << pos << '\t' << record.id << '\t' << taxid << '\t' << node.scientific_name
//                                          << '\t' << node.rank << '\t' << gene_id << '\t' << gene_pos << '\t' << value
//                                          << std::endl;
//                            } else {
//                                auto tokens = Utils::split(record.id, "_");
//                                auto id = tokens[0];
//
//                                if (!rankmap.contains(id)) {
//                                    auto tokens2 = Utils::split(path, "/");
//                                    auto filename = tokens2[tokens2.size()-1];
//                                    tokens2 = Utils::split(filename, ".");
//                                    std::string name = tokens2[0];
//
//                                    rankmap.insert ({ id, Ranks(name) });
//                                }
//                                rankmap[id].ranks[node.rank] += 1;
//                            }
//                            pos++;
//                        }
//                    }
//                }
//            }
//
//
//
//            is.close();
//        }
//        for (auto& rid : rankmap) {
//            std::cout << rid.first << '\t';
//            std::cout << rid.second.name_ << '\t';
//            std::cout << rid.second.GetString() << std::endl;
//        }
//    }
//
//    static void ClassifyMarkerGenes(ClassifyOptions &options) {
//        Benchmark total_bm("Total process");
//        Benchmark preload_bm("Preloading database");
//        Benchmark classify_bm("Classifying reads");
//
//        total_bm.start();
//
//
//
//        std::cout << "DB: " << options.db_dir_ + "/database.bin" << std::endl;
//
//        preload_bm.start();
//        IndexedMap* map = IndexedMap::Load(options.db_dir_ + "/database.bin");
//        preload_bm.stop();
//        preload_bm.printResults();
//
//        Taxonomy::IntTaxonomy taxonomy(options.db_dir_ + "/internal_taxonomy.dmp");
//
//
//        Build::ValueExtractorMG extractor = Build::ValueExtractorMG::Load(options.meta_path_);
//
//        size_t k { (map->key_bits_ + map->offset_bits_) / 2 };
//
//        size_t total_size = 0;
//
//        for (auto file : options.reads_)
//            total_size += Utils::GetFileSize(file);
//
////        const size_t block_size = (10 * 1024 * 1024);
//        const size_t block_size = (1024 * 1024);
//
//        // count records
//        size_t records = 0;
//
//        std::string shape_str = ShapeUtils::LoadShape(options.db_dir_ + "/shape.txt");
//        const size_t shape_length = shape_str.length();
//        bool * shape = ShapeUtils::GetShape(shape_str);
//
//        size_t progress = 0;
//
//        // Maybe get that into the loop at some point, for now only work with one file
//        ifstream is(options.reads_[0], ios::in);
//
//        omp_set_num_threads(options.threads_);
//
//        size_t histo[150] = { 0 };
//
//        std::string snp_detector_path = options.db_dir_ + "/shape.db";
//
//        std::ofstream classification_ofs(options.output_dir_ + ".tsv", ios::binary);
//        std::ofstream snp_ofs(options.output_dir_ + ".snps", ios::binary);
//        const size_t snp_ofs_buffer_size = 1024 * 1024  * 128;
//        const size_t classification_ofs_buffer_size = 1024 * 1024;
////        const size_t snp_ofs_buffer_size = 1;
//
//#if WRITEOUT_FP == 1
//        std::ofstream fp_ofs(options.output_dir_ + ".fq", ios::out);
//#endif
//
//#if KNOWNREF == 1
////        std::unordered_set<int> geneids;
//        std::unordered_map<size_t, std::unordered_map<size_t, size_t>> taxon_total_hits;
//        std::unordered_map<size_t, std::unordered_map<size_t, size_t>> taxon_classified_hits;
//        std::ofstream correction_out(options.output_dir_ + "_correction.tsv", ios::out);
//#endif
//
//        std::cout << std::endl;
//        std::cout << std::string(60, '-') << std::endl;
//        std::cout << "Start classification" << std::endl;
//        std::cout << std::string(60, '-') << std::endl;
//
//        ProgressBar bar;
//        bar.reset(total_size);
//
//        classify_bm.start();
//
//        size_t classified = 0;
//
//        size_t false_positives = 0;
//        size_t no_leaf_hit = 0;
//
//#pragma omp parallel
//        {
//
////            SimpleKmerIterator iterator(k, shape, shape_length);
//            RollingKmerIterator iterator(k, shape, shape_length);
//            BufferedFastxReader reader;
//            FastxRecord record;
//
//            IO::SNP snp;
//            IO::ClassificationLine line;
//            BufferedOutput<IO::SNP> snp_output(snp_ofs_buffer_size);
//            BufferedOutput<IO::ClassificationLine> line_output(classification_ofs_buffer_size);
//
////            HitEvaluator eval(extractor, &taxonomy);
//            ReadClassifier read_classifier(extractor, taxonomy, shape_length);
//            read_classifier.InitProcessor(snp_detector_path, options.sensitive_snps_);
//
//#pragma omp master
//{
//            read_classifier.GetProcessor()->snp_detector_.PrintSpecs();
//            if (shape_str != read_classifier.GetProcessor()->snp_detector_.GetShapeStr()) {
//                std::cout << "shape str is different for shape.db snp_detector and shape.txt" << std::endl;
//                std::cout << shape_str << " " << read_classifier.GetProcessor()->snp_detector_.GetShapeStr() << std::endl;
//                exit(27);
//            }
//            std::cout << std::string(60, '-') << std::endl;
//}
////            VariantInferer vi;
////            std::cout << "Init variant infererer" << std::endl;
//
//            uint64_t key;
//
//            uint64_t record_id = omp_get_thread_num();
//
//
//            while (true) {
//                bool ok = false;
//
//#pragma omp critical(reader)
//                ok = reader.LoadBlock(is, block_size);
//                if (!ok) break;
//
//
//                // Read records from datablock
//                while (true) {
//                    auto valid_fragment = reader.NextSequence(record);
//                    if (!valid_fragment) break;
//#pragma omp atomic
//                    records++;
//
//                    // continue if sequence too short.
//                    if (record.sequence.length() < shape_length) continue;
//
//                    iterator.SetRecord(record);
//
//#if KNOWNREF == 1
//                    auto tokens = Utils::split(record.id, "_");
//                    size_t true_taxid = stoul(tokens[0]);
//                    size_t true_geneid = stoul(tokens[1]);
//
//#pragma omp critical(taxon_hits)
//                    {
//                        if (!taxon_total_hits.contains(true_taxid)) {
//                            std::unordered_map<size_t, size_t> inner;
//                            taxon_total_hits.insert( { true_taxid, inner } );
//                        }
//                        if (!taxon_total_hits[true_taxid].contains(true_geneid)) {
//                            taxon_total_hits[true_taxid][true_geneid] = 0;
//                        }
//                        taxon_total_hits[true_taxid][true_geneid]++;
//
//                        // Read correctly classified
//                        if (!taxon_classified_hits.contains(true_taxid)) {
//                            std::unordered_map<size_t, size_t> inner;
//                            taxon_classified_hits.insert( { true_taxid, inner } );
//                        }
//                        if (!taxon_classified_hits[true_taxid].contains(true_geneid)) {
//                            taxon_classified_hits[true_taxid][true_geneid] = 0;
//                        }
//                    }
//#endif
//
//                    // Debug print header
////                    std::cout << record.header << std::endl;
//
//                    bool call_snps = true;
//
//                    // Temporary
//                    int offset = INT32_MIN;
//
//                    // Reset
//                    read_classifier.Reset();
//                    read_classifier.SetLength(record.sequence.size());
//
//                    int hit_count = 0;
//                    int total_count = 0;
//                    int min_hits = 0;
//
//                    std::string lookups = "";
//
//                    // Iterate over all k-mers
////                    std::cout << std::string(50, '-') << std::endl;
//                    while (iterator.HasNext()) {
//                        // Extract key from k-mer
//                        iterator.operator()(key);
//
//                        if (key == UINT64_MAX) {
//                            call_snps = false;
//                            continue;
//                        }
////                        std::cout << "key: " << key << std::endl;
//
//                        auto entry = map->Find(key);
//
////                        std::cout << KmerUtils::ExpandShape(KmerUtils::ToString(key, k*2), shape, shape_length) << " " << key << " ";
//
//                        if (entry) {
//                            hit_count++;
//                            auto value = entry->value(map->value_bits_);
//
//                            size_t taxid = extractor.GetTaxid(value);
//                            size_t gene_id = extractor.GetGeneId(value);
//                            size_t gene_pos = extractor.GetGenePos(value);
//
//
////                            geneids.insert(gene_id);
//#if KNOWNREF == 1
//
//#if KNOWNREF_DEBUG == 1
//                            lookups += to_string(taxid) + '\t' + to_string(gene_id) + '\t' + to_string(gene_pos) + '\t' + to_string(read_classifier.pos_) + '\t' + taxonomy.Get(taxid).scientific_name + '\n';
////                            if (record.id.starts_with("15581_42"))
////                                printf("%llu\t%llu\t%llu\t%llu\t%s\n", taxid, gene_id, gene_pos, read_classifier.pos_, taxonomy.Get(taxid).scientific_name.c_str());
//#endif
//#endif
//
////                             If gene id is present, taxon has to be a leaf. (not true anymore)
////                            assert(!gene_id or (gene_id > 0 and taxonomy.Get(taxid).IsLeaf()));
//
//                            read_classifier.AddHit(value);
//
//                        } else {
////                            std::cout << "miss" << std::endl;
//                            read_classifier.AddMiss();
//                        }
//                        total_count++;
//                    }
//
////                    size_t taxid = 0, gene_id = -1, best_count = 0, best_total_counts = 0;
//                    size_t taxid = 0, best_count = 0, leaf_hits = 0, best_total_counts = 0;
//                    int gene_pos = -1, pattern_idx = -1, gene_id = -1;
//                    double rank_confidence = 0, candidate_confidence = 0, classification_confidence = 0;
//                    bool forward = true;
//
//                    // Debug, get hits etc.
//#pragma omp atomic
//                    histo[hit_count]++;
//
////                    read_classifier.PrintLeafs();
//                    read_classifier.EvaluateTreeReturnBest(taxid, gene_id, gene_pos, forward, best_count, leaf_hits, best_total_counts, pattern_idx);
//
//
//
//                    gene_pos -= !forward * record.sequence.length();
//
////                    std::cout << "pattern index: " << pattern_idx << " taxid: " << taxid << std::endl;
//                    if (best_total_counts < 4) continue;
//
//#pragma omp atomic
//                    classified++;
//
//                    int mutation_count = read_classifier.EvaluatePatterns(pattern_idx);
//
////                    std::cout << "mutation count " << mutation_count << std::endl;
//
//                    auto& mutations = read_classifier.GetProcessor()->GetMutations();
//
//                    auto cov = read_classifier.GetProcessor()->GetHitCov();
//
////                    std::cout << "cov: " << cov << std::endl;
//
////                    std::cout << "patternidx: " << pattern_idx << std::endl;
////                    if (pattern_idx >= 0) {
////                        read_classifier.PrintPatterns(pattern_idx, total_count);
////                    }
//
//                    bool is_classified = (taxid != 0);
//
//                    // Number of hits supporting the best candidates rank downstream divided by all hits supporting that candidate
//                    rank_confidence = best_total_counts ? best_count / (double) best_total_counts : 0;
//                    // Total counts of the best candidate with respect to all lookup-hits
//                    candidate_confidence = hit_count ? best_total_counts / (double) hit_count : 0;
//                    // Total counts of the best candidate with respect to all lookups
//                    classification_confidence = total_count ? best_total_counts / (double) total_count : 0;
//
////                    if (gene_pos < 0 && gene_pos != -1 && best_count > 1) {
////                        std::cout << std::string(40, '-') << std::endl;
////                        std::cout << record.id << std::endl;
////                        std::cout << "gene_pos: " << gene_pos << std::endl;
////                        std::cout << forward << std::endl;
////                        std::cout << record.sequence << std::endl;
////                        read_classifier.PrintLeafs();
////                        std::cout << "hitcount: " << hit_count << std::endl;
////                        std::cout << "best supporting: " << best_total_counts << std::endl;
////                        std::cout << "best: " << best_count << std::endl;
////                        std::cout << "taxon: " << taxid << " > " << taxonomy.Get(taxid).scientific_name << std::endl;
////                        std::string stop;
////                        std::cin >> stop;
////                    }
//
////                    /*
////                     * VARKIT READ CLASSIFICATION OUTPUT
////                     */
////
//////                    // Parse mutations in header DEBUG
////                    int start = record.header.find("_SNPS=[") + 7;
////                    int n = record.header.length() - start - 1;
////                    auto s = record.header.substr(start, n);
////                    auto tokens = Utils::split(s, ",");
////                    std::unordered_set<int> ref_muts;
////                    for (auto& e : tokens) {
////                        ref_muts.insert(stoi(Utils::split(e, ":")[0]));
////                    }
////                    //DEBUG END
//////                    for (auto& rm : ref_muts) {
//////                        std::cout << rm << ",";
//////                    }
//////                    std::cout << std::endl;
//
//#define BUFFERED_OUTPUT 1
//
//#ifdef BUFFERED_OUTPUT
//                    if (is_classified) {
//#if KNOWNREF == 1
//
//                        if (true_taxid != taxid) {
//
////                            std::cout << "leaf_hits: " << leaf_hits << std::endl;
//                            if (leaf_hits == 0) {
//#pragma omp atomic
//                                no_leaf_hit++;
////                            } else {
////                                std::string stop;
////                                std::cin >> stop;
//                            } else {
//
////                                std::cout << lookups << std::endl;
////                                std::cout << record.id << std::endl;
////                                read_classifier.PrintLeafs();
////                                std::string stop;
////                                std::cin >> stop;
//                            }
//
//
//#pragma omp critical(taxon_hits)
//                            false_positives++;
//#if WRITEOUT_FP == 1
//#pragma omp critical(fpofs)
//{
//                            fp_ofs << record.to_string();
//}
//#endif
//                        } else {
//
//#pragma omp atomic
//                            taxon_classified_hits[true_taxid][true_geneid]++;
//                        }
//#endif
//
//                        line.classified = 'C';
//                        line.read_length = record.sequence.length();
//                        line.record_id = record_id;
//                        line.mutation_count = mutation_count;
//
//                        line.taxid = taxid;
//                        line.geneid = gene_id;
//                        line.genepos = gene_pos;
//
//                        line.total_hits_best = best_total_counts;
//                        line.total_hits = hit_count;
//                        line.leaf_hits_best = best_count;
//                        line.rank_confidence = rank_confidence;
//                        line.candidate_confidence = candidate_confidence;
//                        line.classification_confidence = classification_confidence;
//                        line.SetHeader(record.id);
//
//
////                        std::cout << line.ToString() << std::endl;
//
//                        if (!line_output.Write(line)) {
//#pragma omp critical(line_ofs)
//                            {
//                                line_output.Write(classification_ofs);
//                            }
//                        }
//
//                        // SNP OUTPUT
//                        for (auto& pos : mutations) {
//                            if (!(pos.gene_pos_ >= gene_pos)) {
//                                continue;
//                            }
//                            if (pos.read_pos_ > record.sequence.length() || pos.read_pos_ < 0) {
//                                continue;
//                            }
//
//                            snp.read_id = record_id;
//                            snp.first_read = true;
//                            snp.snp_position = pos.gene_pos_;
//                            snp.snp_quality = record.quality[pos.read_pos_];
//                            snp.snp_base = record.sequence[pos.read_pos_];
//
//                            if (!snp_output.Write(snp)) {
//#pragma omp critical(snp_ofs)
//                                {
//                                    snp_output.Write(snp_ofs);
//                                }
//                            }
//                        }
//                    }
//#else
//
//#pragma omp critical(output_std)
//                    {
//
//                        if (is_classified) {
////// Debug how fast
//                            classification_ofs << "C\t"; //0
//                            classification_ofs << taxonomy.Get(taxid).scientific_name << '\t'; //1
//                            classification_ofs << taxid << '\t'; //2
//                            classification_ofs << best_total_counts << '\t'; //3
//                            classification_ofs << best_count << '\t'; //4
//                            classification_ofs << hit_count << '\t'; //5
//                            classification_ofs << rank_confidence << '\t'; //6
//                            classification_ofs << candidate_confidence << '\t'; //7
//                            classification_ofs << classification_confidence << '\t'; //8
//                            classification_ofs << gene_id << '\t'; //9
//                            classification_ofs << record.id << '\t'; //10
//                            classification_ofs << gene_pos << '\t'; //11
//                            classification_ofs << mutation_count << '\t'; //12
//                            classification_ofs << cov << '\t'; //13
//                            classification_ofs << record.sequence.length() << '\t'; //14
//                            classification_ofs << forward; //15
//                            classification_ofs << '\n'; //
//                        }
//
////                        std::vector<int> v;
////                        v.assign(ref_muts.begin(), ref_muts.end());
////                        std::sort(v.begin(), v.end());
//
////                        std::cout << "mutation_count: " << mutations.size() << std::endl;
////                        if (mutations.size()) {
////                            std::cout << "Ref mutations: ";
////                            for (auto& rm : v) {
////                                std::cout << rm << ",";
////                            }
////                            std::cout << std::endl;
////                            for (auto& mut : mutations) {
////                                std::cout << "(rpos: " << mut.read_pos_ << ", gpos: " << gene_pos << "), ";
////                                std::cout << std::endl;
////                            }
////                        }
//
//
//
//                        for (auto& pos : mutations) {
//////                            std::cout << pos.gene_pos_ << ",";
////                            // DEBUG
////                            if (!ref_muts.contains(pos.gene_pos_)) {
////                                std::cout << "pos: " << pos.gene_pos_ << " false positive______________________" << std::endl;
////                                std::cout << record.id << std::endl;
////                                std::cout << record.sequence << std::endl;
////
////                                auto& hitmiss = read_classifier.GetProcessor()->hitmiss_;
////                                int rpos = 0;
////                                int total_hits = 0;
////                                for (auto& lookup : hitmiss) {
////                                    std::cout << (rpos++) << ":" << lookup.taxid << "(" << lookup.IsGeneSet() << "),";
////                                    if (lookup.taxid != -1) total_hits++;
////                                }
////                                std::cout << std::endl;
////                                std::cout << "total hits: " << total_hits << std::endl;
////
////                                std::cout << pos.gene_pos_ << " rel: " << pos.read_pos_ << std::endl;
////                                std::cout << "forward: " << forward << std::endl;
////                                std::cout << "first: " << read_classifier.GetProcessor()->first_hit_ << std::endl;
////                                std::cout << "last: " << read_classifier.GetProcessor()->last_hit_ << std::endl;
////                                std::cout << "predicted muts: " << mutations.size() << std::endl;
////
//////                                std::string stop;
//////                                std::cin >> stop;
////                            }
////                            //DEBUG END
//////                            std::cout << std::endl;
//
//
//
////                            assert(pos.gene_pos_ >= gene_pos);
//                            if (!(pos.gene_pos_ >= gene_pos)) {
//                                continue;
//                            }
//                            if (pos.read_pos_ > record.sequence.length() || pos.read_pos_ < 0) {
////                                std::cout << "skip" << std::endl;
//                                continue;
//                            }
////                            if (pos.gene_pos_ < gene_pos) {
////                                std::cout << "Thats bad" << std::endl;
////                                std::cout << record.id << std::endl;
////                                std::cout << "taxid: " << taxid << std::endl;
////                                std::cout << pos.gene_pos_ << " read_pos: " << gene_pos << std::endl;
////
//////                                exit(7);
////                            }
//
//                            snp.read_id = record_id;
//                            snp.first_read = true;
//                            snp.snp_position = pos.gene_pos_;
//                            snp.snp_quality = record.quality[pos.read_pos_];
//                            snp.snp_base = record.sequence[pos.read_pos_];
//
////                            char tp = ref_muts.contains(pos.gene_pos_) ? 'O' : 'X';
////                            std::cout << snp.ToString() << " rpos:" << pos.read_pos_ << "   " << tp << std::endl;
//
//                            if (!snp_output.Write(snp)) {
//#pragma omp critical(snp_ofs)
//                                {
//                                    snp_output.Write(snp_ofs);
//                                }
//                            }
//                        }
//                    } // End omp critical output
//#endif
//
//                    record_id += options.threads_;
//
////                    if (record.id.starts_with("6336_19_SNPS"))
////                        exit(8);
//
////                    std::cout << std::string(60, '-') << std::endl;
//                } // End Read Block
//
//#pragma omp atomic
//                progress += reader.LastBlockSize();
//
//#pragma omp critical(updatebar)
//                bar.Update(progress);
//
////#pragma omp critical(updatebar)
////                bar.Update(progress);
//            } // End Reader
//
//            snp_output.Write(snp_ofs);
//            line_output.Write(classification_ofs);
//
//        } // Close OMP parallel
//
//        // close input and output stream
//        snp_ofs.close();
//        classification_ofs.close();
//        is.close();
//
//
////        // Print historgram of hits
////        for (int i = 0; i < 71; i++) {
////            std::cout << i << '\t' << histo[i] << std::endl;
////        }
//
//        delete shape;
//        delete map;
//
//        classify_bm.stop();
//        total_bm.stop();
//
//        auto seconds = classify_bm.GetDuration(Time::seconds);
//
//
//        std::cout << "Processed " << records << fixed << setprecision(2) << " reads (" << ((double)progress/1000000000) << " Gb)" << std::endl;
//        std::cout << "Records per minute: " << ((uint64_t) ((double)records/seconds * 60)) << std::endl;
//        std::cout << "classified: " << classified << std::endl;
//#if KNOWNREF == 1
//        std::cout << "_________" << std::endl;
//        for (auto& taxon : taxon_total_hits) {
//            auto taxid = taxon.first;
//            for (auto& gene : taxon.second) {
//                auto geneid = gene.first;
//                auto total_hits = gene.second;
//                auto true_hits = taxon_classified_hits[taxid][geneid];
//                correction_out << taxid << '\t' << geneid << '\t' << (double)true_hits/total_hits << '\n';
//            }
//        }
//
//        correction_out.close();
//        std::cout << "fp: " << (double)false_positives/classified << std::endl;
//        std::cout << "noleaf fp: "  << (double)no_leaf_hit/false_positives << std::endl;
////        std::cout << "geneids: ";
////        for (auto& gid : geneids) {
////            std::cout << to_string(gid) << ",";
////        }
////        std::cout << std::endl;
//#endif
//
//        preload_bm.printResults();
//        classify_bm.printResults();
//        total_bm.printResults();
//    }
//
//    static void ClassifyMarkerGenesMM(ClassifyOptionsMM &options) {
//        Benchmark total_bm("Total process");
//        Benchmark preload_bm("Preloading database");
//        Benchmark classify_bm("Classifying reads");
//
//        total_bm.start();
//
//        preload_bm.start();
//        MultiMap::MultiMap* map = MultiMap::MultiMap::Load(options.db_dir_ + "/database.bin");
//
//        preload_bm.stop();
//        preload_bm.printResults();
//
//        Taxonomy::IntTaxonomy taxonomy(options.db_dir_ + "/internal_taxonomy.dmp");
//
//
//        Build::ValueExtractorMG extractor = Build::ValueExtractorMG::Load(options.meta_path_);
//
//        size_t k { (map->KeyBits() + map->OffsetBits()) / 2 };
//
//        size_t total_size = 0;
//
//        for (auto file : options.reads_) {
//            total_size += Utils::GetFileSize(file);
//        }
//
//        const size_t block_size = (1024 * 1024);
//
//        size_t records = 0;
//
//        std::string shape_str = ShapeUtils::LoadShape(options.db_dir_ + "/shape.txt");
//        const size_t shape_length = shape_str.length();
//        bool * shape = ShapeUtils::GetShape(shape_str);
//
//        size_t progress = 0;
//
//        // Maybe get that into the loop at some point, for now only work with one file
//        ifstream is(options.reads_[0], ios::in);
//
//        omp_set_num_threads(options.threads_);
//
//        size_t dummy = 0;
//
//        size_t histo[150] = { 0 };
//        size_t count_fp_snps = 0;
//
//        std::string snp_detector_path = options.db_dir_ + "/shape.db";
//
//        std::ofstream classification_ofs(options.output_dir_ + ".tsv", ios::out);
//        std::ofstream snp_ofs(options.output_dir_ + ".snps", ios::out);
//
//
//        std::cout << std::endl;
//        std::cout << std::string(60, '-') << std::endl;
//        std::cout << "Start classification" << std::endl;
//        std::cout << std::string(60, '-') << std::endl;
//
//        ProgressBar bar;
//        bar.reset(total_size);
//
//        classify_bm.start();
//
//        auto lookups = 0;
//        auto nonunique_lookups = 0;
//
//        // record DEBUG multimap unfolds
//        std::vector<uint32_t> lookups_per_read(1000);
//        std::vector<uint32_t> uniques_per_read(1000);
////        std::fill(lookups_per_read, 0, lookups_per_read_size * sizeof(uint32_t));
//
//
//#pragma omp parallel
//        {
//            RollingKmerIterator iterator(k, shape, shape_length);
//            BufferedFastxReader reader;
//            FastxRecord record;
//
//            MultiMap::MapIterator map_iterator(map->ValueBits(), map->KeyBits());
//            map_iterator.map = map->Map();
//
//            ReadClassifierMM read_classifier(extractor, shape_length);
//
//            uint64_t key;
//
//            while (true) {
//                bool ok = false;
//
//#pragma omp critical(reader)
//                ok = reader.LoadBlock(is, block_size);
//                if (!ok) break;
//
//                // Read records from datablock
//                while (true) {
//                    auto valid_fragment = reader.NextSequence(record);
//                    if (!valid_fragment) break;
//#pragma omp atomic
//                    records++;
//
//                    // continue if sequence too short.
//                    if (record.sequence.length() < shape_length) continue;
//
//                    iterator.SetRecord(record);
//
//                    bool call_snps = true;
//
//                    // Temporary
//                    int offset = INT32_MIN;
//
//                    // Reset
//                    read_classifier.Reset();
//                    read_classifier.SetLength(record.sequence.size());
//
//                    int hit_count = 0;
//                    int total_count = 0;
//                    int min_hits = 0;
//
//                    // multimap debug lookups per read
//                    size_t total_lookups = 0;
//                    size_t total_uniques = 0;
//
//                    // Iterate over all k-mers
//                    while (iterator.HasNext()) {
//                        // Extract key from k-mer
//                        iterator.operator()(key);
//
////                        if (key == UINT64_MAX) {
////                            call_snps = false;
////                            continue;
////                        }
//
//                        map->Find(key, map_iterator);
//
////                        std::cout << KmerUtils::ExpandShape(KmerUtils::ToString(key, k*2), shape, shape_length) << " " << key << " ";
//
//
//                        uint64_t main_key = key & ((1llu << map->KeyBits()) - 1);
//                        size_t value;
//                        if (map_iterator.Found()) {
//                            auto it = 0;
//                            hit_count++;
//                            bool unique = true;
//
//                            lookups++;
//
////                            map->print(map_iterator.map_pos, map_iterator.pack_size);
////                            std::cout << main_key << ": " << map_iterator.pack_size << std::endl;
//                            total_lookups += map_iterator.pack_size;
//                            if (map_iterator.pack_size == 1)
//                                total_uniques += 1;
//
////                            while (map_iterator.Next(value)) {
////                                size_t taxid = extractor.GetTaxid(value);
////                                size_t gene_id = extractor.GetGeneId(value);
////                                size_t gene_pos = extractor.GetGenePos(value);
////
////                                // If gene id is present, taxon has to be a leaf.
////                                assert(taxonomy.Get(taxid).IsLeaf());
////
////                                read_classifier.AddHit(value);
////                                it++;
////                                total_lookups++;
////                            }
//
//
////                            std::cout << std::string(10, '-') << std::endl;
//
////                            if (map_iterator.Find(value, insert_value)) {
////                                std::cout << value << " " << insert_value << std::endl;
////                                found = true;
////                                break;
//////                                while ( iterator.Next(key, value, insert_value) ) {
//////                                    if (value == insert_value) {
//////                                        found = true;
//////                                        break;
//////                                    }
//////                                }
////                            }
//
////                            while (!entry->empty() && main_key == entry->key(map->ValueBits())) {
////                                if (it > 0) {
////                                    unique = false;
////                                }
////
////
//////                                std::cout << "before lookup: " << iterator.GetPos() << std::endl;
////                                auto value = entry->value(map->ValueBits());
//////                                std::cout << "after lookup: " << value << std::endl;
////
////
////                                size_t taxid = extractor.GetTaxid(value);
////                                size_t gene_id = extractor.GetGeneId(value);
////                                size_t gene_pos = extractor.GetGenePos(value);
////
//////                                if (it > 300) { // debug
//////                                    std::cout << "it: " << it << " mainkey: " << main_key << " tid: " << taxid << ":" << gene_id << ":" << gene_pos << " " << iterator.GetString() << std::endl;
//////                                }
////
////                                // If gene id is present, taxon has to be a leaf.
////                                assert(taxonomy.Get(taxid).IsLeaf());
////
////                                read_classifier.AddHit(value);
//////                                std::cout << "after add hit " << taxid << ", " << gene_id << ", " << gene_pos << std::endl;
////                                entry++;
////                                it++;
////                                total_lookups++;
////                            }
////                            if (it > 1)
////                                std::cout << "keycount: " << it << std::endl;
//
////                            if (!unique)
////                                nonunique_lookups++;
//
//
//                        } else {
////                            read_classifier.AddMiss();
//                        }
//                        // move read_pos by one
//                        read_classifier.Move();
//                        total_count++;
//                    }
//
//#pragma omp critical(lps)
//                    {
//                    if (total_lookups >= lookups_per_read.size()) {
//                        lookups_per_read.resize(lookups_per_read.size() * 2);
//                        std::cout << "resize lookups per read: " << lookups_per_read.size() << std::endl;
//                    }
//                    if (total_uniques >= uniques_per_read.size()) {
//                        uniques_per_read.resize(lookups_per_read.size() * 2);
//                            std::cout << "resize uniques_per_read per read: " << uniques_per_read.size() << std::endl;
//                    }
//                    }
//#pragma omp atomic
//                    lookups_per_read[total_lookups]++;
//
//
////                    size_t taxid = 0, gene_id = -1, best_count = 0, best_total_counts = 0;
//                    size_t taxid = 0, best_count = 0, best_total_counts = 0;
//                    int gene_pos = -1, pattern_idx = -1, gene_id = -1;
//                    double rank_confidence = 0, candidate_confidence = 0, classification_confidence = 0;
//                    bool forward = true;
//
//                    // Debug, get hits etc.
//#pragma omp atomic
//                    histo[hit_count]++;
//
////                    if (read_classifier.GetLeavesSize() == 0 and read_classifier.GetMaxCount() < 5){
////                        continue;
////                    }
////                    std::cout << "max: " << read_classifier.GetMaxCount() << std::string(20, '-') << std::endl;
////                    std::cout << record.id << std::endl;
////                    read_classifier.PrintLeafs();
////                    read_classifier.Evaluate();
//
//                    continue;
//
//                    gene_pos -= !forward * record.sequence.length();
//
////                    int mutation_count = read_classifier.EvaluatePatterns(pattern_idx);
//
////                    auto& mutations = read_classifier.GetProcessor()->GetMutations();
////                    auto cov = read_classifier.GetProcessor()->GetHitCov();
//
//                    bool is_classified = (taxid != 0);
//
//                    // Number of hits supporting the best candidates rank downstream divided by all hits supporting that candidate
//                    rank_confidence = best_total_counts ? best_count / (double) best_total_counts : 0;
//                    // Total counts of the best candidate with respect to all lookup-hits
//                    candidate_confidence = hit_count ? best_total_counts / (double) hit_count : 0;
//                    // Total counts of the best candidate with respect to all lookups
//                    classification_confidence = total_count ? best_total_counts / (double) total_count : 0;
//
////
////#pragma omp critical(output_std)
////                    {
////
////
////                        if (is_classified) {
////                            classification_ofs << "C\t"; //0
////                            classification_ofs << taxonomy.Get(taxid).scientific_name << '\t'; //1
////                            classification_ofs << taxid << '\t'; //2
////                            classification_ofs << best_total_counts << '\t'; //3
////                            classification_ofs << best_count << '\t'; //4
////                            classification_ofs << hit_count << '\t'; //5
////                            classification_ofs << rank_confidence << '\t'; //6
////                            classification_ofs << candidate_confidence << '\t'; //7
////                            classification_ofs << classification_confidence << '\t'; //8
////                            classification_ofs << gene_id << '\t'; //9
////                            classification_ofs << record.id << '\t'; //10
////                            classification_ofs << gene_pos << '\t'; //11
////                            classification_ofs << mutation_count << '\t'; //12
////                            classification_ofs << cov << '\t'; //13
////                            classification_ofs << record.sequence.length() << '\t'; //14
////                            classification_ofs << forward; //15
////                            classification_ofs << '\n'; //
////
////
////                        } else {
////                            if (options.output_unclassified_) {
////                                classification_ofs << "U\t";//0
////                                classification_ofs << taxonomy.Get(taxid).scientific_name << '\t'; //1
////                                classification_ofs << taxid << '\t'; //2
////                                classification_ofs << best_total_counts << '\t'; //3
////                                classification_ofs << best_count << '\t'; //4
////                                classification_ofs << hit_count << '\t'; //5
////                                classification_ofs << rank_confidence << '\t'; //6
////                                classification_ofs << candidate_confidence << '\t'; //7
////                                classification_ofs << classification_confidence << '\t'; //8
////                                classification_ofs << gene_id << '\t'; //9
////                                classification_ofs << record.id << '\t'; //10
////                                classification_ofs << gene_pos << '\t'; //11
////                                classification_ofs << -1 << '\t'; //12
////                                classification_ofs << cov; //13
////                                classification_ofs << '\n'; //
////                            }
////                        }
////
////
////                        for (auto& pos : mutations) {
////                            assert(pos.gene_pos_ >= gene_pos);
////                            if (pos.gene_pos_ < gene_pos) {
////                                std::cout << record.id << std::endl;
////                                std::cout << pos.gene_pos_ << " read_pos: " << gene_pos << std::endl;
////                                exit(0);
////                            }
////                            snp_ofs << taxid;
////                            snp_ofs << '\t';
////                            snp_ofs << gene_id;
////                            snp_ofs << '\t';
////                            snp_ofs << pos.gene_pos_;
////                            snp_ofs << '\t';
////                            snp_ofs << record.id;
////                            snp_ofs << '\n';
////                        }
////                    } // End output
//
//
//                } // End Read Block
//
//#pragma omp atomic
//                progress += reader.LastBlockSize();
//
//#pragma omp critical(updatebar)
//                bar.Update(progress);
//
////#pragma omp critical(updatebar)
////                bar.Update(progress);
//
//            } // End Reader
//        } // Close OMP parallel
//
//
//        std::cout << "total lookups: " << lookups << std::endl;
//        std::cout << "non unique lookups: " << nonunique_lookups << std::endl;
//
//        // close input and output stream
//        snp_ofs.close();
//        classification_ofs.close();
//        is.close();
//
//
//        delete shape;
//        delete map;
//
//        classify_bm.stop();
//        total_bm.stop();
//
//        auto seconds = classify_bm.GetDuration(Time::seconds);
//
//
//        std::cout << "Processed " << records << fixed << setprecision(2) << " reads (" << ((double)progress/1000000000) << " Gb)" << std::endl;
//        std::cout << "Records per minute: " << ((uint64_t) ((double)records/seconds * 60)) << std::endl;
//
//#define WRITE_COR 1
//#if WRITE_COR == 1
//        // multimap output debug
//        std::ofstream hits_per_read_ofs(options.output_dir_ + "_hpr.tsv", ios::out);
//        bool start = false;
//        size_t succ_lookups = 0;
//        for (int i = lookups_per_read.size()-1; i >= 0; i--) {
//            if (lookups_per_read[i]) {
//                start = true;
//            }
//            if (!start) continue;
//            hits_per_read_ofs << i << "\t" << lookups_per_read[i] << std::endl;
//            succ_lookups += i * lookups_per_read[i];
//        }
//        hits_per_read_ofs.close();
//#endif
//        std::cout << "total succ lookups: " << succ_lookups << std::endl;
//
//        preload_bm.printResults();
//        classify_bm.printResults();
//        total_bm.printResults();
//    }
//}