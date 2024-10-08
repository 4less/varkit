//
// Created by fritsche on 03/01/2022.
//


#pragma once

#include <compact_map/compact_map.h>
#include <Benchmark.h>
#include <kmer_processing/SimpleKmerIterator.h>
#include <kmer_processing/SimpleKmerPutter.h>
#include <omp.h>
#include <bloom_filter.h>
//#include <build/build_scripts.h>
//#include <taxonomy/HitsClassifier.h>
#include <classify/classification.h>
#include "OptionsContainer.h"
//#include "IOHandler.h"
#include "BufferedOutput.h"
#include <utils/BinaryClassifierEvaluator.h>

#define KNOWNREF 1

namespace ClassifyMode {
    const std::string DIVIDER = std::string(60, '-');

    std::string g_read_orientation_str[7] = {
            "WithinGene","ForwardOverhangEnd","ForwardOverhangBegin","ForwardOverhangBoth","ReverseOverhangEnd","ReverseOverhangBegin","ReverseOverhangBoth"
    };
    enum ReadOrientation {
        WithinGene,
        ForwardOverhangEnd,
        ForwardOverhangBegin,
        ForwardOverhangBoth,
        ReverseOverhangEnd,
        ReverseOverhangBegin,
        ReverseOverhangBoth
    };

    static const ReadOrientation GetReadOrientation(bool forward, int gene_pos, size_t read_length, size_t gene_length) {
        if (gene_pos > 0 && gene_pos + read_length < gene_length) {
            return WithinGene;
        }
        if (gene_pos < 0) return forward ? ForwardOverhangBegin : ReverseOverhangBegin;
        if ((gene_pos + read_length) > gene_length) return forward ? ForwardOverhangEnd : ReverseOverhangEnd;
        return WithinGene;
    }

    static const std::pair<uint32_t, uint32_t> OverlappingKmers(bool forward, int gene_pos, size_t read_length, size_t shape_length, size_t gene_length) {
        uint32_t n = std::min(gene_pos + read_length, gene_length) - std::max(0, gene_pos) - shape_length + 1;
        uint32_t start;

        if (gene_pos < 0) {
            if (forward) {
                start = -gene_pos;
            } else {
                start = 0;
            }
        } else {
            if (forward) {
                start = 0;
            } else {
                start = read_length - n - shape_length + 1;
            }
        }
        return { start, n };
    }



    static std::string SummaryString(size_t classified_reads, size_t total_reads, size_t milliseconds) {
        std::string summary_str = "";
        summary_str += "classified " + std::to_string(classified_reads);
        summary_str += "/" + std::to_string(total_reads);
        summary_str += " in " + Utils::FormatMilliseconds(milliseconds);
        return summary_str;
    }

//
//    static void Classify(VarkitOptionsContainer &options, Sample& m_sample) {
//
//        IndexedMap* map = options.map;
//
//        std::ofstream classification_ofs(options.ClassificationOutputFile(m_sample), std::ios::app | std::ios::binary);
//        std::ofstream snp_ofs(options.SNPOutputFile(), std::ios::app | std::ios::binary);
//
//
//        const size_t snp_ofs_buffer_size = 1024 * 1024  * 128;
//        const size_t classification_ofs_buffer_size = 1024 * 1024;
//        const size_t block_size = (1024 * 1024);
//
//        ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);
//
//        ProgressBar bar;
//        bar.reset(options.ReadsTotalSize());
//
//
//        if (!Utils::exists(m_sample.GetSingleEndPath())) {
//            std::cerr << m_sample.GetSingleEndPath() << " does not exist." << std::endl;
//        }
//        ifstream is(m_sample.GetSingleEndPath(), ios::in);
//
//
//        size_t records = 0;
//
//        omp_set_num_threads(options.threads_);
//
//        size_t classified = 0;
//        size_t unclassified = 0;
//        size_t false_positives = 0;
//        size_t no_leaf_hit = 0;
//        size_t progress = 0;
//
//        uint32_t shape_size = options.ShapeLength();
//
//
////        PatternDB* pattern_db = options.HasPatternDB() ? options.pattern_db : new PatternDB(options.ShapeDatabaseFile());
//        PatternDB* pattern_db = new PatternDB(options.ShapeDatabaseFile());
//        ds::PatternMap* pattern_map = new ds::PatternMap(options.ShapeDatabaseFile());
//
////        pattern_db->Print();
//
//        ds::Benchmark classify_bm("Raw classification");
//
//        size_t dummy = 0;
//
//        Benchmark classifier_test_bm("ClassifierTest: ");
//
//#pragma omp parallel
//        {
//            RollingKmerIterator iterator(options.k_, options.Shape(), options.ShapeLength());
//            BufferedFastxReader reader;
//            FastxRecord record;
//
//
//
//            IO::IOSNP snp;
//            IO::ClassificationLine line;
//            BufferedOutput<IO::IOSNP> snp_output(snp_ofs_buffer_size);
//            BufferedOutput<IO::ClassificationLine> line_output(classification_ofs_buffer_size);
//
//
////            ReadClassifier read_classifier(extractor, *options.taxonomy, options.ShapeLength());
//            ReadClassifier read_classifier(extractor, *options.taxonomy, options.ShapeLength());
//
//            // Test
////            HitsClassifier classifier(*options.taxonomy);
//
//            // Performance test
//            Benchmark init_classifier_bm("Init classifiertest");
//            init_classifier_bm.start();
//            classification::Classifier classifier_test(*options.taxonomy);
//            classification::Classifier2 classifier_test2(*options.taxonomy);
//            classifier_test2.InitializePatternProcessor(pattern_map);
////            classification::Classifier2 classifier_test2(*options.taxonomy);
//            init_classifier_bm.stop();
//            init_classifier_bm.printResults();
//
//            read_classifier.InitProcessor(pattern_db, true);
//            classifier_test.InitializePatternProcessor(pattern_map);
//
//            classification::ClassificationResult cresult;
//
////            read_classifier.InitProcessor(options.ShapeDatabaseFile(), true);
//
//#pragma omp master
//            {
////                read_classifier.GetProcessor()->snp_detector_.PrintSpecs();
////                if (options.shape_str != read_classifier.GetProcessor()->pattern_db->shape_str) {
//                if (options.shape_str != read_classifier.GetProcessor()->pattern_db->shape_str) {
//                    std::cout << "shape str is different for shape.db snp_detector and shape.txt" << std::endl;
//                    std::cout << options.shape_str << " " << read_classifier.GetProcessor()->pattern_db->shape_str << std::endl;
////                    if (options.shape_str != read_classifier.GetProcessor()->snp_detector_.GetShapeStr()) {
////                    std::cout << options.shape_str << " " << read_classifier.GetProcessor()->snp_detector_.GetShapeStr() << std::endl;
//                    exit(27);
//                }
//                std::cout << DIVIDER << std::endl;
//
//            }
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
//                // Read records from datablock
//                while (true) {
//                    auto valid_fragment = reader.NextSequence(record);
//                    if (!valid_fragment) break;
//#pragma omp atomic
//                    records++;
//
//                    // continue if sequence too short.
//                    if (record.sequence.length() < options.ShapeLength()) continue;
//
//                    iterator.SetRecord(record);
//
////#if KNOWNREF == 1
////                    auto tokens = Utils::split(record.id, "_");
////                    size_t true_taxid = stoul(tokens[0]);
////                    size_t true_geneid = stoul(tokens[1]);
////
////#pragma omp critical(taxon_hits)
////                    {
////                        if (!taxon_total_hits.contains(true_taxid)) {
////                            std::unordered_map<size_t, size_t> inner;
////                            taxon_total_hits.insert( { true_taxid, inner } );
////                        }
////                        if (!taxon_total_hits[true_taxid].contains(true_geneid)) {
////                            taxon_total_hits[true_taxid][true_geneid] = 0;
////                        }
////                        taxon_total_hits[true_taxid][true_geneid]++;
////
////                        // Read correctly classified
////                        if (!taxon_classified_hits.contains(true_taxid)) {
////                            std::unordered_map<size_t, size_t> inner;
////                            taxon_classified_hits.insert( { true_taxid, inner } );
////                        }
////                        if (!taxon_classified_hits[true_taxid].contains(true_geneid)) {
////                            taxon_classified_hits[true_taxid][true_geneid] = 0;
////                        }
////                    }
////#endif
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
////                    classifier.Reset();
////                    classifier_test.Clear();
////                    classifier_test2.Clear();
//
//                    uint32_t read_length = record.sequence.size();
//                    read_classifier.SetLength(read_length);
//
//                    int hit_count = 0;
//                    int total_count = 0;
//                    int min_hits = 0;
//
//                    std::string lookups = "";
//
//// CLASSIFIER TEST DEBUG
////                    std::string target_header = "21852_18_SNPS=[22:A->T,23:T->C,66:C->G,99:G->A,142:A->T,207:C->G,214:C->G,224:G->T,237:G->A,240:A->T,289:";
////                    if (record.id.starts_with(target_header)) {
////                        std::cout << "\n\n";
////                    }
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
//
//                            hit_count++;
//                            auto value = entry->value(map->value_bits_);
//
//                            uint32_t taxid = extractor.GetTaxid(value);
//                            uint32_t gene_id = extractor.GetGeneId(value);
//                            uint32_t gene_pos = extractor.GetGenePos(value);
//                            uint32_t read_pos = iterator.GetPos();
//
//                            classifier_test_bm.start();
//
//                            read_classifier.AddHit(value);
////                            classifier_test.AddHit(taxid, gene_id, gene_pos, read_pos);
////                            classifier_test2.AddHit(taxid, gene_id, gene_pos, read_pos);
////                            classifier_test_bm.stop();
//
//                            if (taxid == 1) continue;
////                                std::cout << "_____________________________" << taxid << " " <<  gene_id << " " << gene_pos << std::endl;
//
//// CLASSIFIER TEST DEBUG
////                            if (record.id.starts_with(target_header)) {
////                                std::cout << "taxid:   " << taxid;
////                                std::cout << "\tgene_id: " << gene_id;
////                                std::cout << "\tgene_pos: " << gene_pos;
////                                std::cout << "\tread_pos: " << read_pos << "\t" << gene_pos - read_pos << "\t" << gene_pos + read_pos << " " << !classifier_test.m_consistent_strains << "   " << ShapeUtils::ApplyShape(record.sequence, iterator.GetPos(), options.shape, shape_size) << std::endl;
////                            }
//
//
//                        } else {
//                            read_classifier.AddMiss();
////                            classifier_test.AddMiss();
////                            classifier_test2.AddMiss();
//                        }
//                        total_count++;
//                    }
//
////                    std::cout << record.id << std::endl;
////                    classifier_test_bm.start();
////                    classifier_test2.Process();
////                    classifier_test_bm.stop();
//
////
////// DEBUG DEBUG DEBUG DEBUG
////                    if (classifier_test2.GetLookupResults().size() > 10) {
//////                        for (auto& lr : classifier_test2.GetLookupResults()) {
//////                            std::cout << lr.ToString() << std::endl;
//////                        }
//////                        std::cout << "Tree size: " << classifier_test2.GetTree().Data().size() << std::endl;
//////                        std::cout << "_____________________________________" << std::endl;
//////                        std::cout << "HasRoot: " << classifier_test2.HasNode(1) << std::endl;
//////                        classifier_test.Tree().PrettyPrint(classification::PrintCNode);
//////                        std::cout << "____" << std::endl;
//////                        classifier_test2.Process();
//////                        std::cout << "____" << std::endl;
//////
//////                        std::cout << "____" << std::endl;
//////                        classifier_test2.GetTree().PrettyPrint(classification::PrintCNode);
//////                        std::cout << "____" << std::endl;
////
////                        classifier_test2.ProcessHits(cresult);
////
////
////                        classifier_test.Tree().PrettyPrint(classification::PrintC1Node);
////                        std::cout << "____________" << std::endl;
////                        classifier_test2.GetTree().PrettyPrint(classification::PrintCNode);
////
//////                        for (auto& lr : classifier_test2.GetLookupResults()) {
//////                            if (lr.GetLineage().GetTaxonomicId() == 50034) {
//////                                auto& hit = lr.GetHit();
////////                                std::cout << hit.ToString() << std::endl;
//////                                std::cout << hit.m_total_count << std::endl;
//////                                std::cout << hit.m_offset.ToString2() << std::endl;
//////                            }
//////                            std::cout << lr.ToString();
//////                            std::cout << "          IsLEAF " << classifier_test2.IsLeaf(lr.GetLineage().GetTaxonomicId(), lr.GeneId(), 1) << std::endl;
//////                        }
////
////                        std::string stop;
////                        std::cin >> stop;
////
////                    }
//
//
////                    if (!classifier_test.m_consistent_strains) {
////                        std::cout << "Id:" << std::endl;
////                        std::cout << record.id << std::endl;
////                        std::cout << record.sequence << std::endl;
////                        exit(9);
////                    }
//
////                    if (classifier_test.Tree().Data().size() > 20) {
//
////                    if (classifier_test.IsAmbiguous()) {
////                        std::cout << "Record.id: " << record.id << std::endl;
////                        std::cout << classifier_test.Tree().Data().size() << std::endl;
////                        std::cout << classifier_test.Newick() << std::endl;
////                        classifier_test.Tree().PrettyPrint(classification::PrintCNode);
////                        classifier_test.PrintAmbiguous();
////                        std::string stop;
////                        std::cin >> stop;
////                    }
//
//                    size_t taxid = 0, best_count = 0, leaf_hits = 0, best_total_counts = 0;
//                    int gene_pos = -1, pattern_idx = -1, gene_id = -1;
//                    double rank_confidence = 0, candidate_confidence = 0, classification_confidence = 0;
//                    bool forward = true;
//
////                    // Debug, get hits etc.
////#pragma omp atomic
////                    histo[hit_count]++;
//
//
//                    read_classifier.EvaluateTreeReturnBest(taxid, gene_id, gene_pos, forward, best_count, leaf_hits, best_total_counts, pattern_idx);
//
//                    dummy += classifier_test.Dummy();
//                    dummy += classifier_test2.Dummy();
////                    dummy += read_classifier.pos_;
//
////                    std::cout << "test" << std::endl;
//
////                    classifier.Evaluate();
////                    std::cout << "Checkpoint9" << std::endl;
////                    dummy += classifier.dummy_var;
////                    std::string stop;
////                    std::cin >> stop;
//
//                    if (best_total_counts < 1) continue;
////                    if (record.id == "S0R2690790/1") {
////                        read_classifier.PrintLeafs();
////
////                        std::cout << "genepos : " << gene_pos << std::endl;
////                        std::cout << "forward : " << forward << std::endl;
////                        std::cout << " record.sequence.length() : " <<  record.sequence.length() << std::endl;
////                        std::cout << "!forward * record.sequence.length(); : " << !forward * record.sequence.length() << std::endl;
////                    }
//                    gene_pos -= !forward * record.sequence.length();
//
//#pragma omp atomic
//                    classified += 1;
//
//
//                    int mutation_count = read_classifier.EvaluatePatterns(pattern_idx);
////                    int mutation_count = 0;
//
//                    auto& mutations = read_classifier.GetProcessor()->GetMutations();
//
//                    auto cov = read_classifier.GetProcessor()->GetHitCov();
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
//                    dummy += is_classified;
//
////                    if (is_classified) {
////#if KNOWNREF == 1
////                        if (true_taxid != taxid) {
////
//////                            std::cout << "leaf_hits: " << leaf_hits << std::endl;
////                            if (leaf_hits == 0) {
////#pragma omp atomic
////                                no_leaf_hit++;
//////                            } else {
//////                                std::string stop;
//////                                std::cin >> stop;
////                            } else {
////
//////                                std::cout << lookups << std::endl;
//////                                std::cout << record.id << std::endl;
//////                                read_classifier.PrintLeafs();
//////                                std::string stop;
//////                                std::cin >> stop;
////                            }
////
////
////#pragma omp critical(taxon_hits)
////                            false_positives++;
////#if WRITEOUT_FP == 1
////#pragma omp critical(fpofs)
////{
////                            fp_ofs << record.to_string();
////}
////#endif
////                        } else {
////
////#pragma omp atomic
////                            taxon_classified_hits[true_taxid][true_geneid]++;
////                        }
////#endif
////
////                        line.classified = 'C';
////                        line.read_length = record.sequence.length();
////                        line.record_id = record_id;
////                        line.mutation_count = mutation_count;
////
////                        line.taxid = taxid;
////                        line.geneid = gene_id;
////                        line.genepos = gene_pos;
////
////                        line.total_hits_best = best_total_counts;
////                        line.total_hits = hit_count;
////                        line.leaf_hits_best = best_count;
////                        line.rank_confidence = rank_confidence;
////                        line.candidate_confidence = candidate_confidence;
////                        line.classification_confidence = classification_confidence;
////                        line.SetHeader(record.id);
////
//////                        std::cout << "header:  " << std::string(line.buffer) << std::endl;
////
//////                        std::cout << line.ToString() << std::endl;
////
////                        if (!line_output.Write(line)) {
////#pragma omp critical(line_ofs)
////                            {
////                                line_output.Write(classification_ofs);
////                            }
////                        }
////
////                        // SNP OUTPUT
////                        for (auto& pos : mutations) {
////                            if (!(pos.gene_pos_ >= gene_pos)) {
////                                continue;
////                            }
////                            if (pos.read_pos_ > record.sequence.length() || pos.read_pos_ < 0) {
////                                continue;
////                            }
////
////                            snp.read_id = record_id;
////                            snp.first_read = true;
////                            snp.snp_position = pos.gene_pos_;
////                            snp.snp_quality = record.quality[pos.read_pos_];
////                            snp.snp_base = record.sequence[pos.read_pos_];
////
////                            if (!snp_output.Write(snp)) {
////#pragma omp critical(snp_ofs)
////                                {
////                                    snp_output.Write(snp_ofs);
////                                }
////                            }
////                        }
////                    }
//
//                    record_id += options.threads_;
//                } // End Read Block
//
//#pragma omp atomic
//                progress += reader.LastBlockSize();
//
//#pragma omp critical(updatebar)
//                bar.Update(progress);
//
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
//        delete pattern_db;
//
//        auto seconds = classify_bm.GetDuration(static_cast<ds::Time>(Time::seconds));
//        auto milliseconds = classify_bm.GetDuration(static_cast<ds::Time>(Time::milliseconds));
//
//        // classifier debug and test
//        std::cout << "Dummyvar: " << dummy << std::endl;
//        classifier_test_bm.printResults();
//
//
//        std::cout << SummaryString(classified, records, milliseconds) << std::endl;
////        std::cout << "classified: " << classified << std::endl;
////        std::cout << "Processed " << records << fixed << setprecision(2) << " reads (" << ((double)progress/1000000000) << " Gb)" << std::endl;
////        std::cout << "Records per minute: " << ((uint64_t) ((double)records/seconds * 60)) << std::endl;
//    }

//
//    template<bool TestPerfectReads=false>
//    static void ClassifyPEReadsTestFrame(VarkitOptionsContainer &options, Sample& m_sample) {
//        IndexedMap* map = options.map;
//
//        std::ofstream classification_ofs(options.ClassificationOutputFile(m_sample), std::ios::binary | std::ios::app);
//        std::ofstream snp_ofs(options.SNPOutputFile(m_sample), std::ios::binary | std::ios::app);
//
//        const size_t snp_ofs_buffer_size = 1024 * 1024  * 128;
//        const size_t classification_ofs_buffer_size = 1024 * 1024;
//        const size_t block_size = (1024 * 1024);
//
//        classification::ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);
//
//        ProgressBar bar;
//
//        bar.reset(Utils::GetFileSize(m_sample.GetPairedEndPaths().first));
//
//        ifstream is(m_sample.GetPairedEndPaths().first, ios::in);
//        ifstream is2(m_sample.GetPairedEndPaths().second, ios::in);
//
//        size_t records = 0;
//
//        omp_set_num_threads(options.threads_);
//
//        size_t classified = 0;
//        size_t false_positives = 0;
//        size_t no_leaf_hit = 0;
//        size_t progress = 0;
//        size_t agreeing_pair_count = 0;
//        size_t total_kmer_hits = 0;
//
////        PatternDB pattern_db = options.HasPatternDB() ? *options.pattern_db : PatternDB(options.ShapeDatabaseFile());
//        ds::PatternMap* pattern_db = new ds::PatternMap(options.ShapeDatabaseFile());
//
//
//        ds::Benchmark classify_bm("Classification");
//
//#pragma omp parallel
//        {
//            RollingKmerIterator iterator(options.k_, options.Shape(), options.ShapeLength());
//
//            //DEBUG
//            SimpleKmerIterator iterator2(options.k_, options.Shape(), options.ShapeLength());
//
//            BufferedFastxReader reader;
//            BufferedFastxReader reader2;
//
//            FastxRecord record;
//            FastxRecord record2;
//
//            classification::ClassificationResult cresult1;
//            classification::ClassificationResult cresult2;
//
//            IO::IOSNP snp;
//            IO::ClassificationLine line;
//
//            BufferedOutput<IO::IOSNP> snp_output(snp_ofs_buffer_size);
//
//            BufferedOutput<IO::ClassificationLine> line_output(classification_ofs_buffer_size);
//
////            ReadClassifier read_classifier(extractor, *options.taxonomy, options.ShapeLength());
////            read_classifier.InitProcessor(pattern_db, true);
////
////            ReadClassifier read_classifier2(extractor, *options.taxonomy, options.ShapeLength());
////            read_classifier2.InitProcessor(pattern_db, true);
//
//            classification::Classifier2 classifier_test(*options.taxonomy);
//            classifier_test.InitializePatternProcessor(pattern_db);
//
//            classification::Classifier2 classifier_test1(*options.taxonomy);
//            classifier_test1.InitializePatternProcessor(pattern_db);
//            classification::Classifier2 classifier_test2(*options.taxonomy);
//            classifier_test2.InitializePatternProcessor(pattern_db);
//
////            if (read_classifier.GetProcessor()->pattern_db == nullptr) exit(77);
//
//#pragma omp master
//            {
////                if (options.shape_str != read_classifier.GetProcessor()->pattern_db->shape_str) {
////                    std::cout << "shape str is different for shape.db snp_detector and shape.txt" << std::endl;
////                    std::cout << options.shape_str << " " << read_classifier.GetProcessor()->snp_detector_.GetShapeStr() << std::endl;
////                    exit(27);
////                }
//                std::cout << DIVIDER << std::endl;
//            }
//
//            uint64_t key;
//
//            uint64_t record_id = omp_get_thread_num();
//
//            while (true) {
//                bool ok = false;
//                bool ok2 = false;
//
//                // READ FIRST READ
//#pragma omp critical(reader)
//                ok = reader.LoadBatch(is, 2048);
//                if (!ok) break;
//
////                // Read records from datablock
////                while (true) {
////                    auto valid_fragment = reader.NextSequence(record);
////                    if (!valid_fragment) break;
//
//                // READ SECOND READ
//#pragma omp critical(reader2)
//                ok2 = reader2.LoadBatch(is2, 2048);
//                if (!ok2) break;
//
//                // Read records from datablock
//                while (true) {
//                    // Read from read1
//                    auto valid_fragment = reader.NextSequence(record);
//                    // Read from read2
//                    auto valid_fragment2 = reader2.NextSequence(record2);
//                    if (!valid_fragment2) break;
//
//                    assert(FastxUtils::IsValidPair(record, record2));
//
//#pragma omp atomic
//                    records++;
//
////#########################################################################
//                    // FIRST READ
//                    // continue if sequence too short.
//                    if (record.sequence.length() < options.ShapeLength()) continue;
//
//                    bool call_snps = true;
//                    iterator.SetRecord(record);
//                    iterator2.SetRecord(record);
//
//                    //DEBUG
//                    std::string id = Utils::split(record.id, "_").at(0);
//                    std::size_t truth_taxid = std::stoll(id);
//                    //DEBUG end
//
//                    std::cout << "Truth: " << truth_taxid << std::endl;
//
//                    // Reset
//                    classifier_test.Clear();
//                    classifier_test.SetReadLength(record.sequence.size(), false);
//                    classifier_test1.Clear();
//                    classifier_test1.SetReadLength(record.sequence.size(), false);
//
//                    int hit_count = 0;
//                    int total_count = 0;;
//                    int miss_count = 0;
//
//                    //DEBUG
//                    size_t dummy;
//
//                    // Iterate over all k-mers
//                    while (iterator.HasNext()) {
//                        // Extract key from k-mer
//                        iterator.operator()(key);
//
//                        if (key == UINT64_MAX) {
//                            call_snps = false;
//                            continue;
//                        }
//
//                        auto entry = map->Find(key);
//
//                        if (entry) {
//                            hit_count++;
//                            total_kmer_hits++;
//                            auto value = entry->value(map->value_bits_);
//
//                            auto taxid = extractor.GetTaxid(value);
//                            if (!options.taxonomy->IsNodeAncestor(taxid, truth_taxid))
//                                std::cout << extractor.GetTaxid(value) << " " << extractor.GetGeneId(value) << " " << extractor.GetGenePos(value) << " " << ShapeUtils::ApplyShape(record.sequence, iterator.GetPos(), options.Shape(), options.shape_str.length(), true) << " " << key << std::endl;
//
//                            if (!options.taxonomy->IsNodeAncestor(taxid, truth_taxid)) {
//                                std::cout << record.id << " Pos : "  << iterator.GetPos() << std::endl;
//                                std::cout << extractor.GetTaxid(value) << " " << extractor.GetGeneId(value) << " " << extractor.GetGenePos(value) << " " << ShapeUtils::ApplyShape(record.sequence, iterator.GetPos(), options.Shape(), options.shape_str.length(), true) << " " << key << std::endl;
//                                std::cout << "KmerToString: " << KmerUtils::ToString(key, 46) << std::endl;
//                                std::cout << record.sequence << std::endl;
//                                std::cout << KmerUtils::reverseComplement(record.sequence) << std::endl;
//                                std::cout << record.sequence.substr(iterator.GetPos(), 31) << std::endl;
//                                std::cout << record.sequence.substr(iterator.GetPos(), 31) << std::endl;
//                            }
//                            assert(options.taxonomy->IsNodeAncestor(taxid, truth_taxid));
//
//                            classifier_test1.AddHit(
//                                    extractor.GetTaxid(value),
//                                    extractor.GetGeneId(value),
//                                    extractor.GetGenePos(value),
//                                    iterator.GetPos());
//                            classifier_test.AddHit(
//                                    extractor.GetTaxid(value),
//                                    extractor.GetGeneId(value),
//                                    extractor.GetGenePos(value),
//                                    iterator.GetPos());
//                        } else {
//                            classifier_test.AddMiss();
//                            classifier_test1.AddMiss();
//                            miss_count++;
//                        }
//                        total_count++;
//                    }
//
//                    if (classifier_test1.GetLookupResults().size() == 0)
//                        continue;
//
////#########################################################################
//
//                    // SECOND READ
//                    // continue if sequence too short.
//                    if (record.sequence.length() < options.ShapeLength()) continue;
//
//                    bool call_snps2 = true;
//                    iterator.SetRecord(record2);
//                    iterator2.SetRecord(record2);
//
//                    // Reset
//                    classifier_test2.Clear();
//                    classifier_test2.SetReadLength(record.sequence.size(), true);
//
//                    int hit_count2 = 0;
//                    int total_count2 = 0;
//                    int miss_count2 = 0;
//
//                    // Iterate over all k-mers
//                    while (iterator.HasNext()) {
//                        // Extract key from k-mer
//                        iterator.operator()(key);
//
//                        if (key == UINT64_MAX) {
//                            call_snps2 = false;
//                            continue;
//                        }
//
//                        auto entry = map->Find(key);
//
//                        if (entry) {
//                            hit_count2++;
//                            total_kmer_hits++;
//                            auto value = entry->value(map->value_bits_);
//
//                            auto taxid = extractor.GetTaxid(value);
//
//                            if (!options.taxonomy->IsNodeAncestor(taxid, truth_taxid)) {
//                                std::cout << record2.id << " Pos : "  << iterator.GetPos() << std::endl;
//                                std::cout << extractor.GetTaxid(value) << " " << extractor.GetGeneId(value) << " " << extractor.GetGenePos(value) << " " << ShapeUtils::ApplyShape(record2.sequence, iterator.GetPos(), options.Shape(), options.shape_str.length(), true) << " " << key << std::endl;
//                                std::cout << "KmerToString: " << KmerUtils::ToString(key, 46) << std::endl;
//                                std::cout << record2.sequence << std::endl;
//                                std::cout << KmerUtils::reverseComplement(record2.sequence) << std::endl;
//                                std::cout << record2.sequence.substr(iterator.GetPos(), 31) << std::endl;
//                                std::cout << KmerUtils::reverseComplement(record2.sequence.substr(iterator.GetPos(), 31)) << std::endl;
//                            }
//                            assert(options.taxonomy->IsNodeAncestor(taxid, truth_taxid));
//
//                            classifier_test.AddHit(
//                                    extractor.GetTaxid(value),
//                                    extractor.GetGeneId(value),
//                                    extractor.GetGenePos(value),
//                                    iterator.GetPos(),
//                                    true);
//                            classifier_test2.AddHit(
//                                    extractor.GetTaxid(value),
//                                    extractor.GetGeneId(value),
//                                    extractor.GetGenePos(value),
//                                    iterator.GetPos());
//                        } else {
//                            classifier_test.AddMiss(true);
//                            classifier_test2.AddMiss(true);
//                            miss_count2++;
//                        }
//                        total_count2++;
//                    }
//
////#########################################################################
//
//                    classifier_test.ProcessHits(cresult1);
//
//
////                    if (classifier_test1.GetLookupResults().size() > 3) {
////                        std::cout << "First read" << std::endl;
////
////                        std::cout << "Miss count:  " << miss_count << std::endl;
////                        std::cout << "Hit count:   " << hit_count << std::endl;
////                        std::cout << "Total count: " << total_count << std::endl;
////
////                        std::cout << "Merged" << std::endl;
////                        classifier_test.GetTree().PrettyPrint(classification::PrintCNode);
////
////                        auto& best = *cresult1.m_hits[0];
////                        auto& best_hit = best.Data();
////
////                        auto& first_read = best_hit.m_read[0];
////                        auto& second_read = best_hit.m_read[1];
////
////                        std::cout << "Best: " << classification::NodeId::ToString(best.GetNodeId()) << "   " << first_read.ToString() << "  " << second_read.ToString() << std::endl;
////
////                        std::cout << "Hit Type: " << best_hit.GetHitType() << std::endl;
////
////
////                        using classification::HitPE;
////                        if (best_hit.GetHitType() == HitPE::HIT_TYPE::BOTH_STRAIN) {
////                            std::string first = best_hit.m_read[0].m_offset.Forward() ? record.sequence : KmerUtils::reverseComplement(record.sequence);
////                            std::string second = best_hit.m_read[1].m_offset.Forward() ? record2.sequence : KmerUtils::reverseComplement(record2.sequence);
////
////                            std::cout << record.id << std::endl;
////                            int64_t offset_first = best_hit.m_read[0].m_offset.Forward() ? std::get<0>(best_hit.m_read[0].m_offset.GetHits()) : std::get<0>(best_hit.m_read[0].m_offset.GetHits()) + options.ShapeLength() - record.sequence.length();
////                            int64_t offset_second = best_hit.m_read[1].m_offset.Forward() ? std::get<0>(best_hit.m_read[1].m_offset.GetHits()) : std::get<0>(best_hit.m_read[1].m_offset.GetHits()) + options.shape_str.length() - record2.sequence.length();
////
////                            auto min = std::min(offset_first, offset_second);
////                            offset_first -= min;
////                            offset_second -=  min;
////
////                            assert(offset_first >= 0);
////                            assert(offset_second >= 0);
////                            assert(offset_first < 500);
////                            assert(offset_second < 500);
////
////                            std::cout << std::string(offset_first, ' ') << first << " reversed: " << (!first_read.m_offset.Forward()) << std::endl;
////                            std::cout << std::string(offset_second, ' ') << second << " reversed: " << (!second_read.m_offset.Forward()) << std::endl;
////                        }
////
////
//////
//////                        std::cout << "Single1" << std::endl;
//////                        classifier_test1.GetTree().PrettyPrint(classification::PrintCNode);
//////                        std::cout << "\n\nSingle2" << std::endl;
//////                        classifier_test2.GetTree().PrettyPrint(classification::PrintCNode);
//////                        std::cout << "Second read" << std::endl;
////////                        classifier_test2.GetTree().PrettyPrint(classification::PrintCNode);
//////                        std::cout << "Miss count:  " << miss_count2 << std::endl;
//////                        std::cout << "Hit count:   " << hit_count2 << std::endl;
//////                        std::cout << "Total count: " << total_count2 << std::endl;
////
////                        std::string stop;
////                        std::cin >> stop;
////                    }
//
//                    continue;
////
////                    if (rcr.best_total_counts < 1) continue;
////
////                    rcr.gene_pos -= !rcr.forward * record.sequence.length();
////                    rcr2.gene_pos -= !rcr2.forward * record2.sequence.length();
////
////#pragma omp atomic
////                    classified += 1;
////
////                    if(rcr.taxid == rcr2.taxid)
////#pragma omp atomic
////                        agreeing_pair_count++;
////
////
////                    int mutation_count = read_classifier.EvaluatePatterns(rcr.pattern_idx);
////
////                    auto& mutations = read_classifier.GetProcessor()->GetMutations();
////
////                    auto cov = read_classifier.GetProcessor()->GetHitCov();
////
////
////                    bool is_classified = (rcr.taxid != 0);
////
////
////
////                    // Number of hits supporting the best candidates rank downstream divided by all hits supporting that candidate
////                    double rank_confidence = rcr.best_total_counts ? rcr.best_count / (double) rcr.best_total_counts : 0;
////                    // Total counts of the best candidate with respect to all lookup-hits
////                    double candidate_confidence = hit_count ? rcr.best_total_counts / (double) hit_count : 0;
////                    // Total counts of the best candidate with respect to all lookups
////                    double classification_confidence = total_count ? rcr.best_total_counts / (double) total_count : 0;
////
////
////                    if (is_classified) {
////
////                        line.classified = 'C';
////                        line.read_length = record.sequence.length();
////                        line.record_id = record_id;
////                        line.mutation_count = mutation_count;
////
////                        line.taxid = rcr.taxid;
////                        line.geneid = rcr.gene_id;
////                        line.genepos = rcr.gene_pos;
////
////                        line.total_hits_best = rcr.best_total_counts;
////                        line.total_hits = hit_count;
////                        line.leaf_hits_best = rcr.best_count;
////                        line.rank_confidence = rank_confidence;
////                        line.candidate_confidence = candidate_confidence;
////                        line.classification_confidence = classification_confidence;
////                        line.SetHeader(record.id);
////
////
////
////                        if (!line_output.Write(line)) {
////#pragma omp critical(line_ofs)
////                            {
////                                line_output.Write(classification_ofs);
////                            }
////                        }
////
////                        // SNP OUTPUT
////                        for (auto& pos : mutations) {
////                            if (!(pos.gene_pos_ >= rcr.gene_pos)) {
////                                continue;
////                            }
////                            if (pos.read_pos_ > record.sequence.length() || pos.read_pos_ < 0) {
////                                continue;
////                            }
////
////                            snp.read_id = record_id;
////                            snp.first_read = true;
////                            snp.snp_position = pos.gene_pos_;
////                            snp.snp_quality = record.quality[pos.read_pos_];
////                            snp.snp_base = record.sequence[pos.read_pos_];
////
////                            if (!snp_output.Write(snp)) {
////#pragma omp critical(snp_ofs)
////                                {
////                                    snp_output.Write(snp_ofs);
////                                }
////                            }
////                        }
////                    }
//
//                    record_id += options.threads_;
//                } // End Read Block
//
//#pragma omp atomic
//                progress += reader.LastBlockSize();
//
//#pragma omp critical(updatebar)
//                bar.Update(progress);
//
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
//        delete pattern_db;
//
//        auto milliseconds = classify_bm.GetDuration(static_cast<ds::Time>(Time::milliseconds));
//
//        std::cout << "total_kmer_hits: " << total_kmer_hits << std::endl;
//        std::cout << SummaryString(classified, records, milliseconds) << std::endl;
////        std::cout << "classified: " << classified << std::endl;
////        std::cout << "Both reads agree: " << agreeing_pair_count << std::endl;
////        std::cout << "Processed " << records << fixed << setprecision(2) << " reads (" << ((double)progress/1000000000) << " Gb)" << std::endl;
////        std::cout << "Records per minute: " << ((uint64_t) ((double)records/seconds * 60)) << std::endl;
//    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class RunStatistics {
    private:
        size_t m_total_reads = 0;
        size_t m_total_hits = 0;
        size_t m_total_kmers = 0;

    public:
        inline void AddHits(const size_t hits) {
            m_total_hits += hits;
        }
        inline void AddReads(const size_t reads) {
            m_total_reads += reads;
        }
        inline void AddTotalKmers(const size_t kmers) {
            m_total_kmers += kmers;
        }
        inline void IncrementTotalKmers() {
            m_total_kmers++;
        }
        inline void IncrementKmerHits() {
            m_total_hits++;
        }
        inline size_t GetKmerHits() {
            return m_total_hits;
        }
        inline void Clear() {
            m_total_hits = 0;
            m_total_kmers = 0;
            m_total_reads = 0;
        }
    };

    class BenchmarkClassification {
        using TID = classification::TID;
        std::string m_id;

        Taxonomy::IntTaxonomy& m_taxonomy;

        std::unordered_set<TID> found_ids;
        std::unordered_set<TID> true_ids;
        const size_t m_predictors = 0;
        bool m_actual_positive = false;

    public:
        std::vector<BinaryClassifierEvaluator> m_eval_identity;
        std::vector<BinaryClassifierEvaluator> m_eval_ancestry;;

        void WriteTrueIds(std::ostream& os = std::cout) {
            os << "Ids present in reads: " << std::endl;
            for (auto& id : found_ids) {
                os << id << ", ";
            }
            os << std::endl;
        }

        void AddPrediction(size_t predictor, TID prediction) {
            if (predictor >= m_predictors) {
                exit(9);
            }

            if (true_ids.empty()) {
                if (prediction == m_truth) {
                    m_eval_identity[predictor].tp++;
                } else if (!m_taxonomy.IsNodeAncestor(prediction, m_truth)) {
                    m_eval_identity[predictor].fp++;
                }

                if (m_taxonomy.IsNodeAncestor(prediction, m_truth)) {
                    m_eval_ancestry[predictor].tp++;
                } else {
                    m_eval_ancestry[predictor].fp++;
                }
            } else {
                // IF INFORMATION ON TRUE READS IS AVAILABLE (AS IN WHAT IS CONTAINED WITHIN THE DATABASE)
                if (m_actual_positive) {
                    if (prediction == m_truth) {
                        m_eval_identity[predictor].tp++;
                    } else if (!m_taxonomy.IsNodeAncestor(prediction, m_truth)) {
                        m_eval_identity[predictor].fp++;
                    }

                    if (m_taxonomy.IsNodeAncestor(prediction, m_truth)) {
                        m_eval_ancestry[predictor].tp++;
                    } else {
                        m_eval_ancestry[predictor].fp++;
                    }
                } else {
                    m_eval_identity[predictor].fp++;
                    m_eval_ancestry[predictor].fp++;
                }
            }
        }

        void LoadTrueIds(std::string path) {
            std::ifstream ifs(path, std::ios::in);
            std::string line;
            while (std::getline(ifs, line)) {
                true_ids.insert(std::stoi(line));
            }
            ifs.close();
        }



        void AddNoPrediction(size_t predictor) {
            if (true_ids.empty()) {
                m_eval_identity[predictor].fn++;
                m_eval_ancestry[predictor].fn++;
            } else {
                if (m_actual_positive) {
                    m_eval_identity[predictor].fn++;
                    m_eval_ancestry[predictor].fn++;
                } else {
                    m_eval_identity[predictor].tn++;
                    m_eval_ancestry[predictor].tn++;
                }
            }
        }

        BenchmarkClassification(Taxonomy::IntTaxonomy& taxonomy, int predictors) : m_predictors(predictors), m_taxonomy(taxonomy) {
            for (auto i = 0; i < m_predictors; i++) {
                m_eval_identity.emplace_back(BinaryClassifierEvaluator{});
                m_eval_ancestry.emplace_back(BinaryClassifierEvaluator{});
            }
        }

        void SetHeader(std::string id) {
            m_id = id;
            std::size_t truth_taxid = std::stoll(Utils::split(id, "_").at(0));
            m_truth = truth_taxid;
            found_ids.insert(m_truth);
            if (!true_ids.empty() && true_ids.contains(m_truth)) {
                m_actual_positive = true;
            } else {
                m_actual_positive = false;
            }
        }

        TID m_truth;
    };

    struct DebugStruct {
        std::ofstream os;
        size_t truth = 0;
        size_t best = 0;

        DebugStruct(std::string output_file) :
                os(std::ofstream( output_file, std::ios::out )) {
        }

        void DebugWithEvenness(classification::Classifier2 &classifier, Taxonomy::IntTaxonomy &taxonomy, FastxRecord record, classification::ClassificationResult &result) {
            if (!result.HasHits()) return;

            std::string vector_leaf;
            std::string vector_total;

            std::string id = Utils::split(record.id, "_").at(0);
            std::size_t truth_taxid = std::stoll(id);
            truth = truth_taxid;


            auto& best_node = *result.m_hits.at(0);
            auto& best_hit = best_node.Data();

            vector_leaf += std::to_string(
                    best_hit.NodeHitsOnly());
            vector_total += std::to_string(
                    best_hit.TotalHits());


            auto last_hits = best_hit.TotalHits();
            auto all_equal = true;

            for (auto i = 1; i < result.m_hits.size(); i++) {
                vector_leaf += ',' + std::to_string(
                        result.m_hits.at(i)->Data().NodeHitsOnly());
                vector_total += ',' + std::to_string(
                        result.m_hits.at(i)->Data().TotalHits());

                if (last_hits != result.m_hits.at(i)->Data().TotalHits()) all_equal = false;
                last_hits = result.m_hits.at(i)->Data().TotalHits();
            }

            best = best_hit.m_taxonomic_id;
            auto correct = taxonomy.IsNodeAncestor(best_hit.m_taxonomic_id, truth_taxid);

//            auto vec = Utils::split(vector_leaf, ",");
////            if (result.m_hits.size() == 2 && result.m_hits.at(0)->Data().TotalHits() == result.m_hits.at(1)->Data().TotalHits() && correct) {
////            if (vector_leaf == vector_total && all_equal && best_hit.TotalHits() > 1 && result.m_hits.size() > 2) {
////            if (vec.size() == 2 && vec.at(0) == vec.at(1) && correct) {
//            if (result.m_hits.size() > 2 && !correct) {
//                std::cout << "sus.." << std::endl;
//                std::cout << "leaf: " << vector_leaf << "   total: " << vector_total << std::endl;
//                std::cout << record.id << std::endl;
//                for (auto& res : result.m_hits) {
//                    std::cout << res->Data().m_taxonomic_id << '\t';
//                    std::cout << res->Data().m_read[0].m_offset.ToString() << '\t';
//                    std::cout << res->Data().m_read[1].m_offset.ToString() << std::endl;
//                }
//                classifier.GetTree().PrettyPrint(classification::PrintCNode);
//                std::cout << "Prediction: " << best << " correct? " << correct << std::endl;
//                std::cout << "Spotted? " << result.MultipleNodesSameScore() << std::endl;
//
//                std::string stop;
//                std::cin >> stop;
//            }

            os << correct << '\t' << truth_taxid << '\t' << best_hit.m_taxonomic_id << '\t';
            os << vector_leaf << '\t' << vector_total << '\t' << best_hit.m_rank_confidence << '\t' << best_hit.m_candidate_confidence << '\t';
            os << best_hit.TotalHits() << '\t' << best_hit.NodeHitsOnly() << std::endl;



//            std::cout << (truth_taxid == best_hit.m_taxonomic_id) << '\t' << truth_taxid << '\t' << best_hit.m_taxonomic_id << '\t';
//            std::cout << vector_leaf << '\t' << vector_total << '\t' << best_hit.m_rank_confidence << std::endl;
        }
    };


    template<
            typename KmerIterator = RollingKmerIterator,
            typename Map = IndexedMap,
            typename Classifier = classification::Classifier2,
            typename Extractor = ValueExtractorMG,
            typename Statistics = RunStatistics>
    static inline void ProcessRead(
            Map &map,
            FastxRecord &record,
            Classifier &classifier,
            KmerIterator &iterator,
            Extractor &extractor,
            Statistics &statistics,
            bool second) {

        if (record.sequence.length() < iterator.ShapeSize()) return;

        uint64_t key;
        bool call_snps = true;

        iterator.SetRecord(record);

        // Reset
        classifier.SetReadLength(record.sequence.size(), second);

        // Iterate over all k-mers
        while (iterator.HasNext()) {
            // Extract key from k-mer
            iterator.operator()(key);

            if (key == UINT64_MAX) {
                call_snps = false;
                continue;
            }

            auto entry = map.Find(key);

            if (entry) {
                auto value = entry->value(map.value_bits_);
                statistics.IncrementKmerHits();
                classifier.AddHit(
                        extractor.GetTaxid(value),
                        extractor.GetGeneId(value),
                        extractor.GetGenePos(value),
                        iterator.GetPos(),
                        second);
            } else {
                classifier.AddMiss(second);
            }
            statistics.IncrementTotalKmers();
        }
    }

    template<typename Classifier = classification::Classifier2>
    static void ClassifierDebugPE(Classifier &classifier, FastxRecord &read1, FastxRecord &read2, size_t shape_size, classification::ClassificationResult result) {
        //                    if (classifier_test1.GetLookupResults().size() > 3) {

        std::cout << std::string(79, '-') << std::endl;

        std::cout << read1.to_string() << std::endl;
        std::cout << read2.to_string() << std::endl;

        std::cout << "Merged" << std::endl;
        classifier.GetTree().PrettyPrint(classification::PrintCNode);

        auto& best = *result.m_hits[0];
        auto& best_hit = best.Data();

        auto& first_read_hit = best_hit.m_read[0];
        auto& second_read_hit = best_hit.m_read[1];

        std::cout << "Best: " << classification::NodeId::ToString(best.GetNodeId()) << "   " << first_read_hit.ToString() << "  " << second_read_hit.ToString() << std::endl;

        std::cout << "Hit Type: " << best_hit.GetHitType() << std::endl;


        using classification::HitPE;
        if (best_hit.GetHitType() == HitPE::HIT_TYPE::BOTH_STRAIN) {
            std::string first = best_hit.m_read[0].m_offset.Forward() ? read1.sequence : KmerUtils::reverseComplement(read1.sequence);
            std::string second = best_hit.m_read[1].m_offset.Forward() ? read2.sequence : KmerUtils::reverseComplement(read2.sequence);

            std::cout << read1.id << std::endl;
            int64_t offset_first = best_hit.m_read[0].m_offset.Forward() ? std::get<0>(best_hit.m_read[0].m_offset.GetHits()) : std::get<0>(best_hit.m_read[0].m_offset.GetHits()) + shape_size - read1.sequence.length();
            int64_t offset_second = best_hit.m_read[1].m_offset.Forward() ? std::get<0>(best_hit.m_read[1].m_offset.GetHits()) : std::get<0>(best_hit.m_read[1].m_offset.GetHits()) + shape_size - read2.sequence.length();

            auto min = std::min(offset_first, offset_second);


            offset_first -= min;
            offset_second -=  min;

            assert(offset_first >= 0);
            assert(offset_second >= 0);
            assert(offset_first < 500);
            assert(offset_second < 500);

            std::cout << std::string(offset_first, ' ') << first << " reversed: " << (!first_read_hit.m_offset.Forward()) << std::endl;
            std::cout << std::string(offset_second, ' ') << second << " reversed: " << (!second_read_hit.m_offset.Forward()) << std::endl;
        }


        std::string stop;
        std::cin >> stop;
    }



    static void ProcessAndOutput(
            classification::ClassificationResult &result,
            classification::PatternHandler &pattern_handler,
            classification::GeneLengths &gene_lengths,
            IO::IOSNP& output_snp,
            IO::ClassificationLine &line,
            FastxRecord &record,
            size_t record_id,
            size_t shape_length,
            uint8_t read_num,
            BufferedOutput<IO::IOSNP> &snp_output,
            BufferedOutput<IO::ClassificationLine> &line_output,
            std::ofstream &snp_ofs,
            std::ofstream &classification_ofs,
            classification::LineageTIDMap lineages) {

        constexpr bool debug = false;

        size_t debug_id = result.m_taxonomic_id;

        if (result.Success()) {
            if (result.IsLeafHit()) {
                // LEAF HIT
                size_t gene_length = gene_lengths.Get(result.m_taxonomic_id, result.m_gene_id);
                pattern_handler.Evaluate3(result, record.sequence.length(), shape_length, gene_length, lineages);



                // Write SNPs out
                for (auto& snp_pos : pattern_handler.mutations_) {
                    if (snp_pos.read_pos_ >= record.sequence.length()) {
                        std::cerr << snp_pos.ToString() << std::endl;
                        std::cerr << "record length: " << record.sequence.length() << std::endl;
                        std::cerr << result.m_hit->ToString() << std::endl;
                        exit(9);
                    }
                    classification::IoSnpFromMutation(snp_pos, result, output_snp, record, record_id, read_num);
                    if constexpr(debug) {
                        std::cout << snp_pos.ToString() << " -> " << output_snp.ToString() << std::endl;
                    }

//                    // DEBUGDEBUGDEBUG
//                    if (output_snp.snp_position == 223 || output_snp.snp_position == 224) {
//                        std::cout << std::string(79, '-') << std::endl;
//                        std::cout << record.header << std::endl;
//                        std::cout << result.m_hit->ToString() << std::endl;
//                        IO::ClassificationLine templine;
//                        classification::LineFromClassificationResult(templine, result, record, record_id);
//                        std::cout << templine.ToString() << std::endl;
//                        Utils::Input();
//                    }

                    if (output_snp.tax_id != result.m_taxonomic_id || output_snp.gene_id != result.m_gene_id) {
                        std::cout << result.m_hit->ToString() << std::endl;
                        std::cout << output_snp.ToString() << std::endl;
                        exit(23);
                    }

                    if (!snp_output.Write(output_snp)) {
#pragma omp critical(snp_ofs)
                        snp_output.Write(snp_ofs);
                    }
                }
//                if (pattern_handler.mutations_.size()) Utils::Input();

            } else {
                // NO LEAF HIT
                result.SetTotalKmers(record.sequence.length() - shape_length + 1);
            }
            if (result.m_total_kmers == 0 || result.m_result_kmers > result.m_total_kmers) {
                std::cerr << "ERROR1XX: " << result.m_total_kmers << "," << result.m_result_kmers << "," << result.IsLeafHit() << " " << result.m_hit->ToString() << std::endl;
            }
            classification::LineFromClassificationResult(line, result, record, record_id, read_num);

            if (line.taxid != result.m_taxonomic_id || line.geneid != result.m_gene_id) {
                std::cout << result.m_hit->ToString() << std::endl;
                std::cout << line.ToString() << std::endl;
                exit(23);
            }

            if (!line_output.Write(line)) {
#pragma omp critical(line_ofs)
                line_output.Write(classification_ofs);
            }
        }
    }

    static void ClassifyPEReads(VarkitOptionsContainer &options, Sample& sample) {
        constexpr bool debug = false;

        IndexedMap* map = options.map;

        std::ofstream classification_ofs(options.ClassificationOutputFile(sample), std::ios::binary | std::ios::app);
        std::ofstream snp_ofs(options.SNPOutputFile(sample), std::ios::binary | std::ios::app);

        const size_t snp_ofs_buffer_size = 1024 * 1024  * 128;
        const size_t classification_ofs_buffer_size = 1024 * 1024;
        const size_t block_size = (1024 * 1024);

        classification::ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_, options.genepos_bits_);

        ProgressBar bar;

        std::cout << "Filesize:  " <<  Utils::GetFileSize(sample.GetPairedEndPaths().first) << std::endl;
        bar.reset(Utils::GetFileSize(sample.GetPairedEndPaths().first));

        ifstream is(sample.GetPairedEndPaths().first, ios::in);
        ifstream is2(sample.GetPairedEndPaths().second, ios::in);

        size_t records = 0;

        omp_set_num_threads(options.threads_);

        size_t classified = 0;
        size_t false_positives = 0;
        size_t no_leaf_hit = 0;
        size_t progress = 0;
        size_t agreeing_pair_count = 0;
        size_t dummy = 0;

//        PatternDB pattern_db = options.HasPatternDB() ? *options.pattern_db : PatternDB(options.ShapeDatabaseFile());
        ds::PatternMap* pattern_db = new ds::PatternMap(options.ShapeDatabaseFile());

        classification::GeneLengths gene_lengths(options.MarkerGeneLengthsFile(), options.InternalTaxonomy()->MaxLeafId()+1, 121);

        // DEBUG
        std::string debug_struct_out = options.OutputPrefix() + "_debug.tsv";
        DebugStruct debug_struct(debug_struct_out);


        auto min_hits_threshold = 2;

#if KNOWNREF == 1
        // ////////////////////////////////////////////////////////////////////////////////////////
        // ////////////////////////////////////////////////////////////////////////////////////////
        //
        BenchmarkClassification bc(*options.InternalTaxonomy(), 5*5 + 1);
        auto min_hits_threshold_str = " " + std::to_string(min_hits_threshold);

        std::vector<double> anis;

        std::string true_ids = "/usr/users/QIB_fr017/fritsche/Projects/Results/varkit/DBs/gtdb_mg_k23b/library/true_ids.txt";

        bc.LoadTrueIds(true_ids);

        std::string row_header_name = "Name\tMethod\tConf_Threshold\tMin_Hits";

        bc.m_eval_identity[0].SetName("Identity\told\t0\t0");
        bc.m_eval_ancestry[0].SetName("Ancestry\told\t0\t0");

        auto pred_num = 1;
        for (double conf : { 0.5, 0.6, 0.7, 0.8, 0.9 } ) {
            for (int hits : { 1, 2, 3, 4, 5 } ) {
                std::string name_suffix = "";
                name_suffix += '\t';
                name_suffix += std::to_string(conf);
                name_suffix += '\t';
                name_suffix += std::to_string(hits);

                bc.m_eval_identity[pred_num].SetName("Identity\tnew" + name_suffix);
                bc.m_eval_ancestry[pred_num].SetName("Ancestry\tnew" + name_suffix);
                pred_num++;
            }
        }


        //
        // ////////////////////////////////////////////////////////////////////////////////////////
        // ////////////////////////////////////////////////////////////////////////////////////////
#endif

        std::cout << "Classification to " << options.ClassificationOutputFile(sample) << std::endl;
        std::cout << "SNPs to           " << options.SNPOutputFile(sample) << std::endl;
        std::cout << "Debug to          " << debug_struct_out << std::endl;

        int debug_info_counter = 0;
        // DEBUG

        ds::Benchmark classify_bm("Classification");

#pragma omp parallel
        {
            RollingKmerIterator iterator(options.k_, options.Shape(), options.ShapeLength());

            BufferedFastxReader reader;
            BufferedFastxReader reader2;

            FastxRecord record;
            FastxRecord record2;

            classification::ClassificationResult result;
            classification::ClassificationResult result2;
            classification::ClassificationResult result2b;
            classification::ClassificationResult result3;
            classification::ClassificationResult result3b;

            RunStatistics read_statistics;

            IO::IOSNP snp;
            IO::ClassificationLine line;

            BufferedOutput<IO::IOSNP> snp_output(snp_ofs_buffer_size);

            BufferedOutput<IO::ClassificationLine> line_output(classification_ofs_buffer_size);

            classification::Classifier2 classifier(*options.taxonomy);

            classifier.InitializePatternProcessor(pattern_db);

            if (classifier.GetProcessor()->pattern_db == nullptr) exit(77);

#pragma omp master
            {
//                if (options.shape_str != read_classifier.GetProcessor()->pattern_db->shape_str) {
//                    std::cout << "shape str is different for shape.db snp_detector and shape.txt" << std::endl;
//                    std::cout << options.shape_str << " " << read_classifier.GetProcessor()->snp_detector_.GetShapeStr() << std::endl;
//                    exit(27);
//                }
                std::cout << DIVIDER << std::endl;
            }

            uint64_t key;

            uint64_t record_id = omp_get_thread_num();

            while (true) {
                bool ok = false;
                bool ok2 = false;

                // READ FIRST READ
#pragma omp critical(reader)
                ok = reader.LoadBatch(is, 2048);
                if (!ok) break;

                // READ SECOND READ
#pragma omp critical(reader2)
                ok2 = reader2.LoadBatch(is2, 2048);
                if (!ok2) break;

                // Read records from datablock
                while (true) {
                    // Read from read1
                    auto valid_fragment = reader.NextSequence(record);
                    // Read from read2
                    auto valid_fragment2 = reader2.NextSequence(record2);
                    if (!valid_fragment) break;
                    if (!valid_fragment2) break;

                    //assert(FastxUtils::IsValidPair(record, record2));

#pragma omp atomic
                    records++;

                    read_statistics.Clear();
                    classifier.Clear();

                    // FIRST READ
                    ProcessRead(*map, record, classifier, iterator, extractor, read_statistics, false);
                    // SECOND READ
                    ProcessRead(*map, record2, classifier, iterator, extractor, read_statistics, true);

                    if (read_statistics.GetKmerHits() < min_hits_threshold) {
                        continue;
                    }

                    if constexpr(debug) {
                        std::cout << "\nNEXT " << std::string(99,'#') << std::endl;
                        std::cout << record.id << std::endl;
                    }

                    result.Clear();
                    result2.Clear();
                    result2b.Clear();

                    classifier.ProcessHits(result);
                    classifier.ProcessHits3(result2, result2b, 0.7, 3);

                    ProcessAndOutput(result2, *classifier.GetProcessor(0), gene_lengths, snp, line, record, record_id,
                                     iterator.ShapeSize(), true, snp_output, line_output, snp_ofs, classification_ofs,
                                     classifier.GetLineages());
                    ProcessAndOutput(result2b, *classifier.GetProcessor(1), gene_lengths, snp, line, record2, record_id,
                                     iterator.ShapeSize(), true, snp_output, line_output, snp_ofs, classification_ofs,
                                     classifier.GetLineages());

                    if constexpr(false) {
                        if (result2.Success()) {
                            if (result2.IsLeafHit()) {
                                // LEAF HIT
                                size_t gene_length = gene_lengths.Get(result2.m_taxonomic_id, result2.m_gene_id);
                                auto &processor = classifier.GetProcessor(0);
                                processor->Evaluate3(result2, record.sequence.length(), iterator.ShapeSize(),
                                                     gene_length, classifier.GetLineages());

                                // Write SNPs out
                                for (auto &snp_pos: processor->mutations_) {
                                    classification::IoSnpFromMutation(snp_pos, result, snp, record, record_id, true);
                                    if (!snp_output.Write(snp)) {
#pragma omp critical(snp_ofs)
                                        snp_output.Write(snp_ofs);
                                    }
                                }
                            } else {
                                // NO LEAF HIT
                                result2.SetTotalKmers(record.sequence.length() - iterator.ShapeSize() + 1);
                            }

                            classification::LineFromClassificationResult(line, result2, record, record_id);
                            if (!line_output.Write(line)) {
#pragma omp critical(line_ofs)
                                line_output.Write(classification_ofs);
                            }

                            if constexpr(debug) {
                                std::cout << "\n________First read classified as " << result2.m_taxonomic_id
                                          << std::endl;
                                result2.WriteInfo(iterator.ShapeSize(), record.sequence.length());
                            }

                        }
                        if (result2b.Success()) {
                            if (result2b.IsLeafHit()) {
                                size_t gene_length = gene_lengths.Get(result2b.m_taxonomic_id, result2b.m_gene_id);
                                auto &processor = classifier.GetProcessor(1);
                                processor->Evaluate3(result2b, record2.sequence.length(), iterator.ShapeSize(),
                                                     gene_length, classifier.GetLineages());

                                // Write SNPs out
                                for (auto &snp_pos: processor->mutations_) {
                                    classification::IoSnpFromMutation(snp_pos, result, snp, record2, record_id, false);
                                    if (!snp_output.Write(snp)) {
#pragma omp critical(snp_ofs)
                                        snp_output.Write(snp_ofs);
                                    }
                                }
                            } else {
                                result2.SetTotalKmers(record2.sequence.length() - iterator.ShapeSize() + 1);
                            }

                            classification::LineFromClassificationResult(line, result2b, record, record_id);
                            if (!line_output.Write(line)) {
#pragma omp critical(line_ofs)
                                line_output.Write(classification_ofs);
                            }

                            if constexpr(debug) {
                                std::cout << "\n________Second read classified as " << result2.m_taxonomic_id
                                          << std::endl;
                                result2b.WriteInfo(iterator.ShapeSize(), record2.sequence.length());
                            }
                        }
                    }

                    // ####################################################################################################################################################################################


//                    if (is_classified) {
//
//                        if (auto hit1 = result2.GetHitRead1()) {
//                            auto &&[offset1, hits1, forward1] = hit1->GetOffset().GetHits(iterator.ShapeSize(),
//                                                                                          record.sequence.length());
//                            auto &processor1 = classifier.GetProcessor(0);
////                            auto read_orientation1 = GetReadOrientation(forward1, offset1, record.sequence.length(),
////                                                                        gene_length);
//                            auto [start1, n1] = OverlappingKmers(forward1, offset1, record.sequence.length(),
//                                                                 iterator.ShapeSize(), gene_length);
//                            if (hit1->IsLeaf()) {
//                                processor1->Evaluate2(result2.m_taxonomic_id, result2.m_gene_id, offset1, forward1,
//                                                      iterator.ShapeSize(), record.sequence.length(),
//                                                      start1, n1);
//                            }
//
//                        }
//                        if (auto hit2 = result2.GetHitRead2()) {
//
//                        }
//                    }




#if KNOWNREF == 1
                    bc.SetHeader(record.id);

                    if (result.HasBothReadsClassified()) {
                        auto first = result.GetBothClassified().m_taxonomic_id;

                        if constexpr(debug) {
                            if (bc.m_truth != first) {
                                if (!options.InternalTaxonomy()->IsNodeAncestor(first, bc.m_truth)) {
                                    std::cout << "FALSE POSITIVE IDENTITY " << std::endl;
                                } else {
                                    std::cout << "FALSE POSITIVE ANCESTRAL" << std::endl;
                                }
                            }
                        }

//                        std::cout << (first == result2.taxid) << " " << first << " " << result2.taxid << std::endl;
                        bc.AddPrediction(0, first);
                    } else {
                        bc.AddNoPrediction(0);
                    }

                    int pred = 1;
//                    for (double conf : { 0.5, 0.6, 0.7, 0.8, 0.9 } ) {
//                        for (int hits : { 1, 2, 3, 4, 5 } ) {

                    for (double conf : { 0.7 } ) {
                        for (int hits : { 3 } ) {
                            std::cout << "______________________________________________________BENCHMARK" << std::endl;
                            result3.Clear();
                            result3b.Clear();
                            classifier.ProcessHits3(result3, result3b, conf, hits);
                            std::cout << "\tRead1: " <<  result3.Success() << " and Leaf? " << result3.IsLeafHit() << "    " << (result3.m_hit ? result3.m_hit->ToString() : "") << std::endl;
                            std::cout << "\tRead2: " <<  result3b.Success() << " and Leaf? " << result3b.IsLeafHit() << "    " << (result3b.m_hit ? result3b.m_hit->ToString() : "") << std::endl;

                            auto classified = result3.Success() || result3b.Success();
                            if (classified) {

                                result3.SetTotalKmers(record.sequence.length() - iterator.ShapeSize() + 1);
                                result3b.SetTotalKmers(record2.sequence.length() - iterator.ShapeSize() + 1);
                                auto taxid = result3.Success() ? result3.m_taxonomic_id : result3b.m_taxonomic_id;

                                double ratio1 = (double)result3.m_result_kmers/result3.m_total_kmers;
                                double ratio2 = (double)result3b.m_result_kmers/result3b.m_total_kmers;
                                if (result3.Success() && result3.IsLeafHit()) {
                                    std::cout << "\n________Evaluate read 1 (LEAF) " << result3.m_taxonomic_id << " " << result3.m_hit->ToString() << std::endl;
                                    std::cout << record.id << std::endl;
                                    auto gene_length = gene_lengths.Get(taxid, result3.m_gene_id);
                                    classifier.GetProcessor(0)->Evaluate3(result3, record.sequence.length(), iterator.ShapeSize(), gene_length, classifier.GetLineages());
                                    ratio1 = (double)result3.m_result_kmers/result3.m_total_kmers;
                                    auto& snps = classifier.GetProcessor(0)->GetMutations();
                                    if (!snps.empty()) {
                                        auto &&[offset, hits, forward] = result3.m_hit->GetOffset().GetHits(iterator.ShapeSize(), record.sequence.length());
                                        std::cout << "____List snps\t\tRead 1: " << offset << "-" << offset+record.sequence.length() << std::endl;
                                        for (auto snp: snps) {
                                            std::cout << snp.ToString() << std::endl;
                                        }
//                                        Utils::Input();
                                    }
                                }
                                if (result3b.Success() && result3b.IsLeafHit()) {
                                    std::cout << "\n________Evaluate read 2 (LEAF)" << result3b.m_taxonomic_id << " " << result3b.m_hit->ToString() << std::endl;
                                    std::cout << record2.id << std::endl;
                                    auto gene_length = gene_lengths.Get(taxid, result3b.m_gene_id);
                                    classifier.GetProcessor(1)->Evaluate3(result3b, record2.sequence.length(), iterator.ShapeSize(), gene_length, classifier.GetLineages());
                                    ratio2 = (double)result3b.m_result_kmers/result3b.m_total_kmers;
                                    auto& snps = classifier.GetProcessor(1)->GetMutations();
                                    if (!snps.empty()) {
                                        auto &&[offset, hits, forward] = result3b.m_hit->GetOffset().GetHits(iterator.ShapeSize(), record2.sequence.length());
                                        std::cout << "____List snps\t\tRead 2: " << offset << "-" << offset+record2.sequence.length() << std::endl;
                                        for (auto snp: snps) {
                                            std::cout << snp.ToString() << std::endl;
                                        }
//                                        Utils::Input();
                                    }
                                }




                                if constexpr(debug) {
                                    if (bc.m_truth != result3.m_taxonomic_id) {
                                        if (options.InternalTaxonomy()->IsNodeAncestor(taxid, bc.m_truth)) {
                                            std::cout << "FALSE POSITIVE IDENTITY (SECOND METHOD)" << bc.m_truth << " "
                                                      << taxid << "   ";
                                            if (result3.Success()) {
                                                std::cout << "   1: " << std::to_string(ratio1);
                                                if (result3.IsLeafHit()) {
                                                    auto ani = AbundanceEstimator::EstimateAni(result3.m_result_kmers, result3.m_total_kmers, options.k_);
                                                    std::cout << " ani: " << ani;
                                                    anis.emplace_back(ani);                                               }
                                            }
                                            if (result3b.Success()) {
                                                std::cout << "   2: " << std::to_string(ratio2);
                                                if (result3b.IsLeafHit()) {
                                                    auto ani = AbundanceEstimator::EstimateAni(result3b.m_result_kmers, result3b.m_total_kmers, options.k_);
                                                    std::cout << " ani: " << ani;
                                                    anis.emplace_back(ani);
                                                }
                                            }
                                            std::cout << std::endl;
                                        } else {
                                            std::cout << "FALSE POSITIVE ANCESTRAL (SECOND METHOD)" << bc.m_truth << " "
                                                      << taxid << "   ";
                                            if (result3.Success()) {
                                                std::cout << "   1: " << std::to_string(ratio1);
                                                if (result3.IsLeafHit()) {
                                                    auto ani = AbundanceEstimator::EstimateAni(result3.m_result_kmers, result3.m_total_kmers, options.k_);
                                                    std::cout << " ani: " << ani;
                                                    anis.emplace_back(ani);                                               }
                                            }
                                            if (result3b.Success()) {
                                                std::cout << "   2: " << std::to_string(ratio2);
                                                if (result3b.IsLeafHit()) {
                                                    auto ani = AbundanceEstimator::EstimateAni(result3b.m_result_kmers, result3b.m_total_kmers, options.k_);
                                                    std::cout << " ani: " << ani;
                                                    anis.emplace_back(ani);
                                               }
                                            }
                                            std::cout << std::endl;
                                        }
                                    } else {
                                        std::cout << "TRUE POSITIVE  " << bc.m_truth << " "
                                                  << taxid << "   ";
                                        if (result3.Success()) {
                                            std::cout << "   1: " << std::to_string(ratio1);
                                            if (result3.IsLeafHit()) {
                                                auto ani = AbundanceEstimator::EstimateAni(result3.m_result_kmers, result3.m_total_kmers, options.k_);
                                                std::cout << " ani: " << ani;
                                                anis.emplace_back(ani);
                                            }
                                        }
                                        if (result3b.Success()) {
                                            std::cout << "   2: " << std::to_string(ratio2);
                                            if (result3b.IsLeafHit()) {
                                                auto ani = AbundanceEstimator::EstimateAni(result3b.m_result_kmers, result3b.m_total_kmers, options.k_);
                                                std::cout << " ani: " << ani;
                                                anis.emplace_back(ani);

                                            }
                                        }
                                        std::cout << std::endl;
                                    }
                                }
                                bc.AddPrediction(pred, taxid);
                            } else {
                                bc.AddNoPrediction(pred);
                            }
                            pred++;
                        }
                    }
#endif



                    if constexpr (debug) {
                        std::cout << std::endl << std::string(79, '#') << std::endl;
                        classifier.GetTree().PrettyPrint(classification::PrintCNode);
                        std::cout << "Header: " << record.id << std::endl;

                        if (result.HasBothReadsClassified()) {
                            std::cout << std::string(10, '#') << " First classifier" << std::endl;
                            std::cout << "result1:  " << result.GetBothClassified().ToVerboseString() << "\n"
                                      << std::endl;
                        }
//                        if (true) {
//                            std::cout << std::string(10, '#') << " Second classifier" << std::endl;
//                            std::cout << "result2:  " << result2.GetBothClassified().ToVerboseString() << std::endl;
//                            std::cout << "HasHit(0):  " << result2.GetBothClassified().HasHit(0) << std::endl;
//                            std::cout << "HasHit(1):  " << result2.GetBothClassified().HasHit(1) << std::endl;
//                            auto first = result2.GetBothClassified().GetHit(false);
//                            auto second = result2.GetBothClassified().GetHit(true);
//                            bool first_maps = first.m_total_count >= min_hits_threshold;
//                            bool second_maps = second.m_total_count >= min_hits_threshold;
//                            bool only_one_maps = first.m_total_count < min_hits_threshold ||
//                                                 second.m_total_count < min_hits_threshold;
//                            std::cout << only_one_maps << " " << first.m_total_count << ", " << second.m_total_count
//                                      << std::endl;
//
//                            auto taxonomic_id = result2.GetBothClassified().m_taxonomic_id;
//                            auto gene_id = result2.GetBothClassified().m_gene_id;
//                            auto gene_length = gene_id > 0 ? gene_lengths.Get(taxonomic_id, gene_id) : 0;
//
//                            if (first_maps && first.GetOffset().IsLeaf()) {
//                                auto &&[offset1, hits1, forward1] = first.GetOffset().GetHits(iterator.ShapeSize(),
//                                                                                              record.sequence.length());
//
//                                auto &processor1 = classifier.GetProcessor(0);
//                                auto read_orientation1 = GetReadOrientation(forward1, offset1, record.sequence.length(),
//                                                                            gene_length);
//
//                                auto [start1, n1] = OverlappingKmers(forward1, offset1, record.sequence.length(),
//                                                                     iterator.ShapeSize(), gene_length);
//
//                                std::cout << "Record1 length: " << record.sequence.length() << std::endl;
//                                std::cout << "Read1: " << g_read_orientation_str[read_orientation1] << "  start: "
//                                          << start1 << " n: " << n1 << std::endl;
//                                std::cout << "Range1: " << offset1 << " - " << (offset1 + record.sequence.length())
//                                          << std::endl;
//                                std::cout << "Evaluate: \n"
//                                          << processor1->Evaluate2(taxonomic_id, gene_id, offset1, forward1,
//                                                                   iterator.ShapeSize(), record.sequence.length(),
//                                                                   start1, n1) << std::endl;
//                            }
//
//                            if (second_maps && second.GetOffset().IsLeaf()) {
//                                auto &&[offset2, hits2, forward2] = second.GetOffset().GetHits(iterator.ShapeSize(),
//                                                                                               record2.sequence.length());
//
//                                auto &processor2 = classifier.GetProcessor(1);
//                                auto read_orientation2 = GetReadOrientation(forward2, offset2,
//                                                                            record2.sequence.length(), gene_length);
//
//                                auto [start2, n2] = OverlappingKmers(forward2, offset2, record2.sequence.length(),
//                                                                     iterator.ShapeSize(), gene_length);
//
//                                std::cout << "Record1 length: " << record.sequence.length() << std::endl;
//                                std::cout << "Read1: " << g_read_orientation_str[read_orientation2] << "  start: "
//                                          << start2 << " n: " << n2 << std::endl;
//                                std::cout << "Range1: " << offset2 << " - " << (offset2 + record.sequence.length())
//                                          << std::endl;
//                                std::cout << "Evaluate: \n"
//                                          << processor2->Evaluate2(taxonomic_id, gene_id, offset2, forward2,
//                                                                   iterator.ShapeSize(), record.sequence.length(),
//                                                                   start2, n2) << std::endl;
//                            }
//
//
//                            if (only_one_maps) {
//                                std::cout << "Only one maps." << std::endl;
//
////                            Utils::Input();
//
//
//                            } else {
//                                std::cout << "Both map";
//                                if (first.GetOffset().IsLeaf() && second.GetOffset().IsLeaf()) {
//                                    std::cout << ", both are leafs" << std::endl;
//                                    auto &&[offset1, hits1, forward1] = first.GetOffset().GetHits(iterator.ShapeSize(),
//                                                                                                  record.sequence.length());
//                                    auto &&[offset2, hits2, forward2] = second.GetOffset().GetHits(iterator.ShapeSize(),
//                                                                                                   record2.sequence.length());
//                                    auto min = std::min(offset1, offset2);
//                                    int insert_size = (offset1 < offset2 ? offset2 -
//                                                                           (offset1 + record.sequence.length()) :
//                                                       offset1 - (offset2 + record2.sequence.length()));
//                                    std::cout << "offset1: " << offset1 << "  offset2: " << offset2 << std::endl;
//                                    std::cout << "offset1: " << offset1 - min << "  offset2: " << offset2 - min
//                                              << std::endl;
//                                    std::cout << std::string(offset1 - min, ' ')
//                                              << (!forward1 ? KmerUtils::reverseComplement(record.sequence)
//                                                            : record.sequence) << std::endl;
//                                    std::cout << std::string(offset2 - min, ' ')
//                                              << (!forward2 ? KmerUtils::reverseComplement(record2.sequence)
//                                                            : record2.sequence) << std::endl;
//                                    std::cout << "INSERT SIZE: " << insert_size << std::endl;
//
//                                } else {
//                                    std::cout << ", at least one is no leaf" << std::endl;
//                                }
//                            }
//                        }
                    }



                    if constexpr(false) {
//                    result2.Clear();
//                    auto succ = classifier.ProcessHits2(result2, 0.9, 5);
//                    if (succ && !options.InternalTaxonomy()->IsNodeAncestor(result2.taxid, bc.m_truth)) {
//                        classifier.GetTree().PrettyPrint(classification::PrintCNode);
//                        classifier.PrintLookupResults();
//                        std::cout << "TRUTH: " << bc.m_truth << std::endl;
//                        std::cout << "PREDI: " << result2.taxid << std::endl;
//                        std::string stop;
//                        std::cin >> stop;
//                    }

//
//
//                    auto& hitpe = result.GetBothClassified();
////                    if (classifier.GetLookupResults().size() > 10) {
//                    if (hitpe.TotalHits() > 2 && result.HasBothReadsClassified() && hitpe.m_taxonomic_id != result2.taxid) {
//                        classifier.PrintLookupResults();
//                        classifier.GetTree().PrettyPrint(classification::PrintCNode);
//
//                        std::cout << "TRUTH:  " << record.id << std::endl;
//
//                        if (result.HasBothReadsClassified())
//                            std::cout << "pred1: " << result.GetBothClassified().m_taxonomic_id << std::endl;
//                        std::cout << "pred2: " << result2.taxid << std::endl;F
//
//
//                        std::string stop;
//                        std::cin >> stop;
//                    } else {
//                        if (result.HasBothReadsClassified() && hitpe.m_taxonomic_id == result2.taxid) {
//                            std::cout << "yay :) " << std::endl;
//                        }
//                    }


//                        // Change to LCA later
//                        if (result.MultipleNodesSameScore()) {
//                            continue;
//                        }

                        // DEBUG
//#pragma omp critical(debug_struct)
//                    debug_struct.DebugWithEvenness(classifier, *options.taxonomy, record, result);

//                    if (debug_struct.truth)
//                    if (debug_struct.best == 9097 && debug_struct.truth == 12270) {
//                        //12270,9097
//                        std::cout << record.id << std::endl;
//                        std::cout << classifier.Newick() << std::endl;
//                        std::string stop;
//                        std::cin >> stop;
//                    }
                        // DEBUG END

//
//                    std::cout << result.m_hits.size() << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
//                    for (auto& node : result.m_hits) {
//                        auto& data = node->Data();
//                        std::cout << node->Data().m_taxonomic_id << " " << node->Data().TotalHits() << std::endl;
//                    }
//                    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

//                    ClassifierDebugPE(classifier, record, record2, iterator.ShapeSize(), result);
                    }


                    if (result.HasBothReadsClassified()) {
                        // Both reads support same evidence
                        auto& hit = result.GetBothClassified();
                        classified += 1;
                        dummy += hit.TotalHits();
                        classification::LineFromHitPE(line, hit, record_id, record, record2, iterator.ShapeSize(), read_statistics.GetKmerHits());

                        auto taxid = hit.m_taxonomic_id;
                        auto gene_id = hit.m_gene_id;


                        if constexpr (debug) {
                            if (hit.GetHit().IsLeaf()) {
                                std::cout << "TAXID: " << taxid << " GENEID: " << gene_id << ":  "
                                          << gene_lengths.Get(taxid, gene_id) << std::endl;
                            }
                        }
//                        auto input = Utils::Input();

//                        if (line.genepos < -20'000 && line.genepos > INT32_MIN) {
//                            std::cout << line.ToString() << std::endl;
//                            classifier.PrintLookupResults();
//                            std::string stop;
//                            std::cin >> stop;
//                        }

                        if (debug_info_counter > 0) {
                            std::cout << std::string(79, '-') << std::endl;
                            classifier.PrintLookupResults();
                            std::cout << std::string(39, '-') << std::endl;
                            classifier.GetTree().PrettyPrint(classification::PrintCNode);
                            std::cout << hit.ToVerboseString() << std::endl;
                            std::cout << line.ToString() << std::endl;
                            debug_info_counter--;
                        }

//                        classifier.PrintLookupResults();
//                        ClassifierDebugPE(classifier, record, record2, iterator.ShapeSize(), result);

//                        std::cout << line.ToString() << std::endl;

                        if (!line_output.Write(line)) {
#pragma omp critical(line_ofs)
                            {
                                line_output.Write(classification_ofs);
                            }
                        }

                    } else {
                        // One read is does not have hits
                        if (result.HasFirstReadClassified()) {
                            auto& hit = result.GetFirst();
//                            line.SetHit(hit, record, record_id, false);
                        }

                        if (result.HasSecondReadClassified()) {
                            auto& hit = result.GetSecond();
//                            line.SetHit(hit, record, record_id, true);
                        }
                    }

                    auto& best = *result.m_hits[0];
                    auto& best_hit = best.Data();
                    auto taxonomic_id = classification::NodeId::GetTaxonomicId(best.GetNodeId());
                    auto gene_id = classification::NodeId::GetTaxonomicId(best.GetNodeId());

                    dummy += taxonomic_id;
                    dummy += gene_id;


                    if (best_hit.m_total_count < min_hits_threshold) continue;

                    // Required output
                    auto type = best_hit.GetHitType();

//                    //processor->Evaluate(leaf.taxid, leaf.geneid, rel_pos, forward, shape_size_);
//
//                    int mutation_count = classifier.EvaluatePatterns(rcr.pattern_idx);
//
//                    int mutation_count = read_classifier.EvaluatePatterns(rcr.pattern_idx);
//                    auto& mutations = read_classifier.GetProcessor()->GetMutations();
//                    auto cov = read_classifier.GetProcessor()->GetHitCov();
//
//
//                    bool is_classified = (rcr.taxid != 0);
//
//
//
//                    // Number of hits supporting the best candidates rank downstream divided by all hits supporting that candidate
//                    double rank_confidence = rcr.best_total_counts ? rcr.best_count / (double) rcr.best_total_counts : 0;
//                    // Total counts of the best candidate with respect to all lookup-hits
//                    double candidate_confidence = hit_count ? rcr.best_total_counts / (double) hit_count : 0;
//                    // Total counts of the best candidate with respect to all lookups
//                    double classification_confidence = total_count ? rcr.best_total_counts / (double) total_count : 0;
//
//
//                    if (is_classified) {
//
//                        line.classified = 'C';
//                        line.read_length = record.sequence.length();
//                        line.record_id = record_id;
//                        line.mutation_count = mutation_count;
//
//                        line.taxid = rcr.taxid;
//                        line.geneid = rcr.gene_id;
//                        line.genepos = rcr.gene_pos;
//
//                        line.total_hits_best = rcr.best_total_counts;
//                        line.total_hits = hit_count;
//                        line.leaf_hits_best = rcr.best_count;
//                        line.rank_confidence = rank_confidence;
//                        line.candidate_confidence = candidate_confidence;
//                        line.classification_confidence = classification_confidence;
//                        line.SetHeader(record.id);
//
//
//
//                        if (!line_output.Write(line)) {
//#pragma omp critical(line_ofs)
//                            {
//                                line_output.Write(classification_ofs);
//                            }
//                        }
//
////                        // SNP OUTPUT
////                        for (auto& pos : mutations) {
////                            if (!(pos.gene_pos_ >= rcr.gene_pos)) {
////                                continue;
////                            }
////                            if (pos.read_pos_ > record.sequence.length() || pos.read_pos_ < 0) {
////                                continue;
////                            }
////
////                            snp.read_id = record_id;
////                            snp.first_read = true;
////                            snp.snp_position = pos.gene_pos_;
////                            snp.snp_quality = record.quality[pos.read_pos_];
////                            snp.snp_base = record.sequence[pos.read_pos_];
////
////                            if (!snp_output.Write(snp)) {
////#pragma omp critical(snp_ofs)
////                                {
////                                    snp_output.Write(snp_ofs);
////                                }
////                            }
////                        }
//                    }

                    record_id += options.threads_;
                } // End Read Block

#pragma omp atomic
                progress += reader.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);

            } // End Reader

            snp_output.Write(snp_ofs);
            line_output.Write(classification_ofs);

        } // Close OMP parallel

        // close input and output stream
        snp_ofs.close();
        classification_ofs.close();
        is.close();
        delete pattern_db;

        auto milliseconds = classify_bm.GetDuration(static_cast<ds::Time>(ds::Time::milliseconds));

#if KNOWNREF == 1
        bc.WriteTrueIds();
        bc.m_eval_identity[0].WriteRowHeader(std::cout, row_header_name);
        for (auto i = 0; i < bc.m_eval_identity.size(); i++) {
            bc.m_eval_identity[i].WriteRowStats();
        }
        for (auto i = 0; i < bc.m_eval_identity.size(); i++) {
            bc.m_eval_ancestry[i].WriteRowStats();
        }
        std::cout << "Mean ani: " << Statistics::Mean(anis);
        std::cout << "Median ani: " << Statistics::Median(anis);
#endif


        std::cout << dummy << std::endl;
        std::cout << SummaryString(classified, records, milliseconds) << std::endl;
//        std::cout << "classified: " << classified << std::endl;
//        std::cout << "Both reads agree: " << agreeing_pair_count << std::endl;
//        std::cout << "Processed " << records << fixed << setprecision(2) << " reads (" << ((double)progress/1000000000) << " Gb)" << std::endl;
//        std::cout << "Records per minute: " << ((uint64_t) ((double)records/seconds * 60)) << std::endl;
    }


    static void ClassifyPEReadsNew(VarkitOptionsContainer &options, Sample& sample) {
        constexpr bool debug = true;

        IndexedMap *map = options.map;

//        std::ofstream classification_ofs(options.ClassificationOutputFile(sample), std::ios::binary | std::ios::app);
//        std::ofstream snp_ofs(options.SNPOutputFile(sample), std::ios::binary | std::ios::app);
        std::ofstream classification_ofs(options.ClassificationOutputFile(sample), std::ios::binary);
        std::ofstream snp_ofs(options.SNPOutputFile(sample), std::ios::binary);


        const size_t snp_ofs_buffer_size = 1024 * 1024 * 128;
        const size_t snp_ofs_buffer_size_mb = 256;

        const size_t classification_ofs_buffer_size = 1024 * 1024;
        const size_t classification_ofs_buffer_size_mb = 256;

        const size_t block_size = (1024 * 1024);


        classification::ValueExtractorMG extractor(options.internal_taxid_bits_, options.geneid_bits_,
                                                   options.genepos_bits_);

        ProgressBar bar;

        std::cout << "Filesize:  " << Utils::GetFileSize(sample.GetPairedEndPaths().first) << std::endl;
        bar.reset(Utils::GetFileSize(sample.GetPairedEndPaths().first));

        ifstream is(sample.GetPairedEndPaths().first, ios::in);
        ifstream is2(sample.GetPairedEndPaths().second, ios::in);

        size_t records = 0;

        omp_set_num_threads(options.threads_);

        size_t classified = 0;
        size_t progress = 0;
//        size_t false_positives = 0;
//        size_t no_leaf_hit = 0;
//        size_t agreeing_pair_count = 0;

        ds::PatternMap *pattern_db = new ds::PatternMap(options.ShapeDatabaseFile());

        classification::GeneLengths gene_lengths(options.MarkerGeneLengthsFile(), options.InternalTaxonomy()->MaxLeafId()+1, 121);


        auto min_hits_threshold = 2;

        std::cout << "Classification to " << options.ClassificationOutputFile(sample) << std::endl;
        std::cout << "SNPs to           " << options.SNPOutputFile(sample) << std::endl;
        std::cout << "options.MarkerGeneLengthsFile(): " << options.MarkerGeneLengthsFile() << std::endl;

        // DEBUG
        ds::Benchmark classify_bm("Classification");

        if (!options.taxonomy) {
            exit(8);
        }

#pragma omp parallel
        {
            RollingKmerIterator iterator(options.k_, options.Shape(), options.ShapeLength());

            BufferedFastxReader reader_r1;
            BufferedFastxReader reader_r2;

            FastxRecord record_r1;
            FastxRecord record_r2;

            classification::ClassificationResult result_r1;
            classification::ClassificationResult result_r2;

            RunStatistics read_statistics;

            IO::IOSNP snp;
            IO::ClassificationLine line;

            BufferedOutput<IO::IOSNP> snp_output(snp_ofs_buffer_size_mb, BufferedOutputUnit::MB);
//            BufferedOutput<IO::IOSNP> snp_output(snp_ofs_buffer_size);
            BufferedOutput<IO::ClassificationLine> line_output(classification_ofs_buffer_size_mb, BufferedOutputUnit::MB);
//            BufferedOutput<IO::ClassificationLine> line_output(classification_ofs_buffer_size);

            classification::Classifier2 classifier(*options.taxonomy);
            classifier.InitializePatternProcessor(pattern_db);

            if (classifier.GetProcessor(true)->pattern_db == nullptr) exit(77);
            if (classifier.GetProcessor(false)->pattern_db == nullptr) exit(77);

#pragma omp master
            {
//                if (options.shape_str != read_classifier.GetProcessor()->pattern_db->shape_str) {
//                    std::cout << "shape str is different for shape.db snp_detector and shape.txt" << std::endl;
//                    std::cout << options.shape_str << " " << read_classifier.GetProcessor()->snp_detector_.GetShapeStr() << std::endl;
//                    exit(27);
//                }
                std::cout << DIVIDER << " Master: " << omp_get_thread_num() << std::endl;
            }

            uint64_t key;

            uint64_t record_id = omp_get_thread_num();

            while (true) {
                bool ok = false;
                bool ok2 = false;

                // READ FIRST READ
#pragma omp critical(reader)
                ok = reader_r1.LoadBatch(is, 2048);
                if (!ok) break;

                // READ SECOND READ
#pragma omp critical(reader2)
                ok2 = reader_r2.LoadBatch(is2, 2048);
                if (!ok2) break;

                // Read records from datablock
                while (true) {
                    // Read from read1
                    auto valid_fragment = reader_r1.NextSequence(record_r1);
                    // Read from read2
                    auto valid_fragment2 = reader_r2.NextSequence(record_r2);
                    if (!valid_fragment) break;
                    if (!valid_fragment2) break;

//                    assert(FastxUtils::IsValidPair(record_r1, record_r2));

#pragma omp atomic
                    records++;

                    read_statistics.Clear();
                    classifier.Clear();

                    // FIRST READ
                    ProcessRead(*map, record_r1, classifier, iterator, extractor, read_statistics, false);
                    // SECOND READ
                    ProcessRead(*map, record_r2, classifier, iterator, extractor, read_statistics, true);

                    if (read_statistics.GetKmerHits() < min_hits_threshold) {
                        continue;
                    }

                    result_r1.Clear();
                    result_r2.Clear();

                    classifier.ProcessHits3(result_r1, result_r2, 0.7, 3);

#pragma omp atomic
                    classified += result_r1.Success() || result_r2.Success();

                    ProcessAndOutput(result_r1, *classifier.GetProcessor(0), gene_lengths, snp, line, record_r1, record_id,
                                     iterator.ShapeSize(), 0, snp_output, line_output, snp_ofs, classification_ofs,
                                     classifier.GetLineages());
                    ProcessAndOutput(result_r2, *classifier.GetProcessor(1), gene_lengths, snp, line, record_r2, record_id,
                                     iterator.ShapeSize(), 1, snp_output, line_output, snp_ofs, classification_ofs,
                                     classifier.GetLineages());

                    record_id += options.threads_;

                } // End Read Block

#pragma omp atomic
                progress += reader_r1.LastBlockSize();

#pragma omp critical(updatebar)
                bar.Update(progress);

            } // End Reader

#pragma omp critical(snp_ofs)
            snp_output.Write(snp_ofs);
#pragma omp critical(line_ofs)
            line_output.Write(classification_ofs);

        } // Close OMP parallel

        // close input and output stream
        snp_ofs.close();
        classification_ofs.close();
        is.close();
        delete pattern_db;

        auto milliseconds = classify_bm.GetDuration(static_cast<ds::Time>(ds::Time::milliseconds));
        std::cout << SummaryString(classified, records, milliseconds) << std::endl;
    }


    static void Run(VarkitOptionsContainer &options, Sample& sample) {
        std::cout << "CLASSIFY::RUN" << std::endl;
        auto class_file = options.ClassificationOutputFile(sample);
        auto snp_file = options.SNPOutputFile(sample);

        if (Utils::exists(class_file)) {
            std::remove(class_file.c_str());
        }
        if (Utils::exists(snp_file)) {
            std::remove(snp_file.c_str());
        }


//        m_sample.SetSingleEnd(m_sample.GetPairedEndPaths().first);
        if (sample.HasSingleEnd()) {
//            Classify(options, m_sample);
            exit(8);
        }
        if (sample.HasPairedEnd()) {
//            ClassifyPEReads(options, m_sample);
            ClassifyPEReadsNew(options, sample);
//            ClassifyPEReadsTestFrame(options, m_sample);
        }
    }
}