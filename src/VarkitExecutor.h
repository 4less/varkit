//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by joachim on 01/07/2020.
//

#pragma once


#include <iostream>
#include <FastaReader.h>
#include <FastxReader.h>
#include <thread>
#include <SpacedKmerIterator.h>
#include "io/MetaDataDB.h"
#include "utils/ShapeUtils.h"
#include "kmer_processing/KmerProcessor.h"
#include "OptionsContainer.h"
#include <omp.h>
#include <KmerBuffer.h>
#include <set>
#include <VariantInferer.h>

using namespace std;

namespace VarkitExecutor {
    
    static const string PATTERN_MAP_PATH = "/shape/pattern_map.csv";
    static const string SHAPE_PATH = "/shape/shape";
    static const string TAXONOMY_NODES_PATH = "/taxonomy/nodes.dmp";
    static const string TAXONOMY_NAMES_PATH = "/taxonomy/names.dmp";
    static const string DB_PATH = "/kmer.db";
    static const string META_PATH = "/meta.txt";
    
    
    static bool* getShape(string s) {
        bool * shape = new bool[s.length()];
        const char * c = s.c_str();
        for (int i = 0; i < s.length(); i++) {
            shape[i] = c[i] == '_';
        }
        return shape;
    }
    
    static inline void getKey(uint8_t* key, size_t k, size_t key_bytes, const char * sequence, bool* shape, size_t shape_size, size_t pos) {
        memset(key, 0, key_bytes);
    
        int offset = 0;
        for (int i = 0; i < k; i++) {
            while (shape[i+offset]) {
                offset++;
            }
            key[i/4] = key[i/4] | KmerUtils::getBitFromBase(*(sequence + (pos + i + offset))) << (2*(3-(i%4)));
        }
    }
    
    static inline void getKeyRC(uint8_t* key, size_t k, size_t key_bytes, const char * sequence, bool* shape, size_t shape_size, size_t pos) {
        memset(key, 0, key_bytes);
        
        int offset = 0;
        for (int i = 0; i < k; i++) {
            while (shape[shape_size-1-(i+offset)]) {
                offset++;
            }
            key[i/4] = key[i/4] | KmerUtils::getBitFromBaseC(*(sequence + (pos + (shape_size-1) - (i + offset)))) << (2*(3-(i%4)));
        }
    }
    
    static inline void getKeyTest(uint8_t* key_f, uint8_t* key_r, size_t key_bytes, const char * sequence, bool* shape, size_t shape_size) {
        memset(key_f, 0, key_bytes);
        int k_offset = 0;
        for (int s_pos = 0; s_pos < shape_size; s_pos++) {
            k_offset += !shape[s_pos];
            key_f[(s_pos - k_offset) / 4] |= shape[s_pos] * KmerUtils::getBitFromBase(*(sequence + s_pos)) << (2 * (3 - ((s_pos - k_offset) % 4)));
            key_r[(s_pos - k_offset) / 4] |= shape[s_pos] * KmerUtils::getBitFromBaseC(*(sequence + (shape_size - s_pos - 1))) << (2 * (3 - ((s_pos - k_offset) % 4)));
        }
    }
    
    static inline bool isSmallerThan(uint8_t * a, uint8_t * b, size_t key_bytes) {
        for (int i = 0; i < key_bytes; i++) {
            if (a[i] != b[i])
                return (a[i] < b[i]);
        }
        return true;
    }
    
    
    
    template<int KeyBits, int ValueBits, class T>
    static void runVerify(BHashMap<KeyBits, ValueBits, T> &map, BuildOptionsContainer &options) {
        string shape_str = MetaDataDB::loadShape(options.db + MetaDataDB::SHAPE_FILE);
        const size_t shape_length = shape_str.length();
        bool * shape = getShape(shape_str);
        
        int k = KeyBits/2;
        const int key_bytes = (2*k + 8 - 1) /  8;
        const int value_bytes = (ValueBits + 8 - 1) /  8;
        
        
        static size_t block_size = 200000;
        
        // Input stream
        std::istream* is = nullptr;
        
        // Iterate over all input files
        // set up file
        is = new std::ifstream(options.reference[0]);
        if (!is) {
            cerr << options.reference[0] << " does not exist." << endl;
            exit(0);
        }
        
        omp_set_num_threads(options.threads);
        
        //Multithreading
#pragma omp parallel
        {
#pragma omp critical(print_thread_start)
            std::cout << "start: " << omp_get_thread_num() << std::endl;
            
            BufferedFastxReader reader = BufferedFastxReader();
            FastxRecord record = FastxRecord();
            
            uint8_t *key = new uint8_t[key_bytes];
            uint64_t value = 1;
            memset(key, 0, key_bytes);
            
            // Iterator
            SpacedKmerIterator iterator(k, shape, shape_length);
            
            unsigned int read_id = 0;
            
            int unclassified = 0;
            int classified = 0;
            
            // Outer loop to load data in blocks (blocksize specified in number of records)
            while (true) {
                bool ok = false;

#pragma omp critical(seqread)
                ok = reader.LoadBlock(*is, block_size);
                if (!ok) break;
                
                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;
                    iterator.setRecord(record);
                    
                    ++read_id;
                    
                    // count processed records
                    if ((read_id % 1000) == 0)
                        cout << "processed record: " << read_id << endl;
                    
                    
                    // Iterate over all k-mers
                    while (iterator.hasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);
                        value = map.search(key);
                        if (!value) {
                            cout << record.to_string();
                            cout << record.header << endl;
                            cout << iterator.getSubstring(iterator.getPos()-1,45);
                            exit(0);
                        }
                    }
                }
            }
            delete[] key;
        } // parallel end
        
        cout << "Validation process finished. All k-mers were found." << endl;
    }

    
    template<int KeyBits, int ValueBits, class T>
    static void runBuildWorker(BuildOptionsContainer &options) {
        ds::Benchmark build_bm("Build process");
        build_bm.start();
        
        // Shape information
        cout << "shape_path: " << options.db + MetaDataDB::SHAPE_FILE << endl;
        string shape_str = MetaDataDB::loadShape(options.db + MetaDataDB::SHAPE_FILE);
        const size_t shape_length = shape_str.length();
        bool * shape = getShape(shape_str);
        
        const int k = KeyBits/2;
        const int key_bytes = (2*k + 8 - 1) /  8;
        const int value_bytes = (ValueBits + 8 - 1) /  8;
    
        static size_t block_size = 200000;
    
        // Load Taxonomy
        TaxonomyInterface *taxonomy;
        if (options.taxonomy.compare("NCBI") == 0) {
            taxonomy = new NCBITaxonomy();
            taxonomy->loadCustomNodes(options.db + MetaDataDB::TAX_NODES_FILE);
            taxonomy->loadCustomNames(options.db + MetaDataDB::TAX_NAMES_FILE);
        } else if (options.taxonomy.compare("GTDB") == 0) {
            taxonomy = new GTDBTaxonomy();
            taxonomy->loadCustomNodes(options.db + MetaDataDB::TAX_NODES_FILE);
            taxonomy->loadCustomNames(options.db + MetaDataDB::TAX_NAMES_FILE);
        } else {
            std::cerr << options.taxonomy << " is no valid taxonomy. Valid taxonomies are NCBI and GTDB." << endl;
            exit(0);
        }
        
        // Initialize Map
        cout << "init capa: " << options.initial_capacity << endl;
        BHashMap<KeyBits,ValueBits,T> map( options.initial_capacity , 0.8 , 1.5);
        cout << "Initial size in MB: " << map.getMemoryInMB() << endl;
        cout << "files: " << options.reference.size() << endl;
        cout << "shape: " << shape_str << endl;
        cout << "shape_length: " << shape_length << endl;
        
        // Input stream
        std::istream* is = nullptr;
        // FIX set IS BEFORE PARALLEL
        is = new std::ifstream(options.reference[0]);
        if (!is) {
            cerr << "There is no such file: " << options.reference[0] << endl;
            exit(0);
        }
    
        //stats
        int record_count = 0;
        int skipped_count = 0;
        
        
        std::cout << "Number of threads: " << omp_get_thread_num() << std::endl;
        omp_set_num_threads(1);
        
        uint8_t* key_temp = nullptr;
        uint8_t* value_temp = nullptr;
        
        //Multithreading
#pragma omp parallel
        {
            #pragma omp critical(print_thread_start)
            std::cout << "start: " << omp_get_thread_num() << std::endl;
            
            BufferedFastxReader reader = BufferedFastxReader();
            FastxRecord record = FastxRecord();
            
            static thread_local uint8_t *key = new uint8_t[key_bytes];
            static thread_local uint8_t *key_rc = new uint8_t[key_bytes];
            static thread_local uint64_t value = 1;
            static thread_local uint64_t tax_id = 0;
        
            memset(key, 0, key_bytes);
            memset(key_rc, 0, key_bytes);
        
            // Iterator and putter
            SpacedKmerIterator iterator(k, shape, shape_length);
            KmerPutter<KeyBits, ValueBits, T> putter = KmerPutter<KeyBits, ValueBits, T>();

            // Initialize Putter
            putter.setMap(&map);
            putter.setTaxonomy(taxonomy);
        
            //Init KmerBuffer
            //KmerBuffer buffer(3125000, key_bytes, value_bytes);
            
            // stats
            uint64_t kmer_count = 0;
            uint64_t rand_calc = 0;
        
            // Iterate over all input files
            for (const auto &file : options.reference) {
            
                // Outer loop to load data in blocks (blocksize specified in number of records)
                while (true) {
                    bool ok = false;
                    
                    #pragma omp critical(reader)
                    ok = reader.LoadBlock(*is, block_size);
                    if (!ok) break;
                
                    // Read records from datablock
                    while (true) {
                        auto valid_fragment = reader.NextSequence(record);
                        if (!valid_fragment) break;
                        iterator.setRecord(record);
                        cerr << record.sequence.length() << endl;
                        
                        // count processed records ecoli k12
                        if ((record_count % 10) == 0)
                            cout << "processed record: " << record_count << "(of which " << skipped_count
                                 << " were skipped)" << endl;
                        ++record_count;
                        
//                        if (record_count >= 5570) {
//                            cout << record.header << endl;
//                            map.save(options.db + "/index_size");
//                        }

                        // Extract taxonomic identifier from sequence header (has to be a ncbi identifier e.g.: >813)
                        value = stoll(record.header.substr(1));
                        
//                        tax_id = options.taxonomy.compare("GTDB") == 0 ? taxonomy->getCustom(value) : taxonomy->getCustom(record.header.substr(1));

                        // depends on header prep of references.
                        //tax_id = taxonomy->getCustom(value);
                        tax_id = value;
                        if (tax_id == -1) {
                            cerr << "unknown taxid: " << tax_id << " for value " << value << endl;
                            exit(0);
                            continue;
                        }
                        
                        // Skip if ncbi id is unknown to local taxonomy subset
                        if (!taxonomy->getNode(tax_id)) {
                            cerr << "local taxonomy subset has no node (entry) for NCBI identifier: " << tax_id << endl;
                            cerr << "Skip read (please refer to the option -w/--write_unknown to extract skipped reads)." << endl;
                            skipped_count++;
                            exit(0);
                            continue;
                        }
                        
                        // Iterate over all k-mers
                        while (iterator.hasNext()) {
                            // Extract key from k-mer
                            iterator.operator()(key);
                            kmer_count++;
                            
                            //buffer.pushKV(key, ((uint8_t*)&value));
                            
                            //#pragma omp critical(putter)
                            putter.operator()(key, tax_id);
                            
                            if (!map.search(key)) {
                                cout << "fail with: " << endl;
                                cout << record.to_string() << endl;
                                exit(0);
                            }
                        }
                    }
//                    while (buffer.next()) {
//                        //(key_temp = buffer.popKey();
//                        //rand_calc += (int) key_temp[0] + (int) key_temp[1] + (int) key_temp[2];
//                        #pragma omp critical(putter)
//                        putter.operator()(buffer.popKey(), tax_id);
//                    }
//                    buffer.reset();
                }
                cout << "kmercount: " << kmer_count << endl;
                cout << "randcalc: " << rand_calc << endl;
                break; //FIX CODE
            }
    
            delete[] key;
            delete[] key_rc;
        }
        
        build_bm.stop();
        build_bm.printResults();
        
        // exclude while testing key extracting methods etc.
        map.printVars();
        map.printStats();
        cout << "save db to " << (options.db + "/index.db") << endl;
        map.save(options.db + "/index");
        cout << "done saving." << endl;
        
        delete is;
        delete[] shape;
        delete taxonomy;
        
        cout << "validate: " << options.validate << endl;
        if (options.validate) {
            cout << "validate" << endl;
            runVerify(map, options);
        }
    }
    
    
    static void runBuild(BuildOptionsContainer &options) {
        int key_bits = 44;
        int value_bits = 20;
        if (key_bits == 44 && value_bits == 20) {
            cout << options.db << endl;
            runBuildWorker<44, 20, CustomHash>(options);
        }
    }

    
    template<int KeyBits, int ValueBits, class T>
    static void runClassifyThreadWorker(BHashMap<KeyBits, ValueBits, T> &map, FastaReader &reader, ClassifyOptionsContainer &options) {
    
        string shape_str = MetaDataDB::loadShape(options.meta_db.path + MetaDataDB::SHAPE_FILE);
        const size_t shape_length = shape_str.length();
        bool * shape = getShape(shape_str);
        int k = KeyBits/2;
        SpacedKmerIterator iterator(k, shape, shape_length);
        
        // Load Taxonomy
        NCBITaxonomy taxonomy = NCBITaxonomy();
        taxonomy.loadCustomNodes(options.meta_db.path + MetaDataDB::TAX_NODES_FILE);
        taxonomy.loadCustomNames(options.meta_db.path + MetaDataDB::TAX_NAMES_FILE);
    }

    
    template<int KeyBits, int ValueBits, class T>
    static void runClassifyWorkerOMP(ClassifyOptionsContainer &options, BHashMap<KeyBits, ValueBits, T> &map, std::string file) {
        cout << " start classify" << endl;
        ds::Benchmark classify_bm("Classify process");
        classify_bm.start();
        
        // Shape information
        cout << "shape_path: " << options.meta_db.path + MetaDataDB::SHAPE_FILE << endl;
        string shape_str = MetaDataDB::loadShape(options.meta_db.path + MetaDataDB::SHAPE_FILE);
        const size_t shape_length = shape_str.length();
        bool * shape = getShape(shape_str);
    
        cout << "shape: " << shape_str << endl;
        cout << "shape_length: " << shape_length << endl;
        
        const int k = KeyBits/2;
        const int key_bytes = (2*k + 8 - 1) /  8;
        const int value_bytes = (ValueBits + 8 - 1) /  8;
        
        static size_t block_size = 200000;
//
//        // Load Taxonomy
//        NCBITaxonomy taxonomy = NCBITaxonomy();
//        taxonomy.loadCustomNodes(options.meta_db.path + MetaDataDB::TAX_NODES_FILE);
//        taxonomy.loadCustomNames(options.meta_db.path + MetaDataDB::TAX_NAMES_FILE);
//
        // Load Taxonomy
        TaxonomyInterface *taxonomy;
        if (options.meta_db.taxonomy.compare("NCBI") == 0) {
            taxonomy = new NCBITaxonomy();
            taxonomy->loadCustomNodes(options.meta_db.path + MetaDataDB::TAX_NODES_FILE);
            taxonomy->loadCustomNames(options.meta_db.path + MetaDataDB::TAX_NAMES_FILE);
        } else if (options.meta_db.taxonomy.compare("GTDB") == 0) {
            taxonomy = new GTDBTaxonomy();
            taxonomy->loadCustomNodes(options.meta_db.path + MetaDataDB::TAX_NODES_FILE);
            taxonomy->loadCustomNames(options.meta_db.path + MetaDataDB::TAX_NAMES_FILE);
        } else {
            std::cerr << options.meta_db.taxonomy << " is no valid taxonomy. Valid taxonomies are NCBI and GTDB." << endl;
            exit(0);
        }
        
        
        //Initialize patternMap
        BHashMap<32,64,SuperFastHash> pattern_map;
        pattern_map.load(options.meta_db.path + "/shape/shape.db", options.meta_db.path + "/shape/shape.meta");
        
        // Input stream
        // Iterate over all input files
        // set up file
        std::istream* is = new std::ifstream(file);
        if (!is) {
            cerr << file << " does not exist." << endl;
            exit(0);
        } else {
            cout << "open: " << file << endl;
        }
        
        std::string output_basename = Utils::stripExtension(file);
        std::string classified_file = options.output + "/" + output_basename + ".c";
        std::string unclassified_file = options.output + "/" + output_basename + ".u";
        std::ostream* os_classified = new std::ofstream(classified_file, ios::out);
        std::ostream* os_unclassified = new std::ofstream(unclassified_file, ios::out);
        if (!os_classified) {
            cerr << "unable to open " << classified_file << endl;
            exit(0);
        } else {
            cout << "open: " << classified_file << endl;
        }
        if (!os_unclassified) {
            cerr << "unable to open " << unclassified_file << endl;
            exit(0);
        } else {
            cout << "open: " << unclassified_file << endl;
        }
        
        
        
        omp_set_num_threads(options.threads);
        
        unsigned int global_read_num = 0;
        unsigned int global_classified = 0;
        unsigned int global_unclassified = 0;
        

        
        //Multithreading
#pragma omp parallel
        {
#pragma omp critical(print_thread_start)
            std::cout << "start: " << omp_get_thread_num() << std::endl;
            
            BufferedFastxReader reader = BufferedFastxReader();
            FastxRecord record = FastxRecord();
            
            uint8_t *key = new uint8_t[key_bytes];
            uint8_t *value_cache = new uint8_t[8];
            uint64_t value = 1;
            
            memset(key, 0, key_bytes);
            
            // Iterator
            SpacedKmerIterator iterator(k, shape, shape_length);
            
            // variant inferer
//            VariantInferer vi = VariantInferer(options.meta_db.path + options.meta_db.SNP_PATTERN_FILE, 1000);
            VariantInferer vi;
            vi.setMap(&pattern_map);
//            vi.load(options.meta_db.path + "/shape/shape");
            uint8_t hit_pattern[4];
            static const size_t vmer_length = options.vmer_length;
            
            unsigned int read_id = 0;
            unsigned int thread_id = omp_get_thread_num();
            
            uint64_t found_count = 0;
            uint64_t unfound_count = 0;
            
            int unclassified = 0;
            int classified = 0;
            
            // taxon counter
            TaxonomyNodeComparer t_cmp(taxonomy);
            set<int, TaxonomyNodeComparer> taxa(t_cmp);
            unordered_map<int, int> taxon_counter;
            
            ostringstream oss_classified;
            ostringstream oss_unclassified;
        
            // Outer loop to load data in blocks (blocksize specified in number of records)
            while (true) {
                bool ok = false;

#pragma omp critical(seqread)
                ok = reader.LoadBlock(*is, block_size);
                if (!ok) break;
                
                oss_classified.str(""); // reset
                oss_unclassified.str("");
                
                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;
                    iterator.setRecord(record);
                    
                    //reset counter for taxonomic classification
                    taxon_counter.clear();
                    taxa.clear();
                    
                    // reset for variant infering
                    memset(hit_pattern,0,4);
                    vi.clear();
                    vi.setLength(record.sequence.length());
    
    
                    int hit_count = 0;
                    int read_kmer_count = 0;
                    
                    ++read_id;
    
                    // count processed records
                    if ((read_id % 1000) == 0)
                        cout << "processed record: " << read_id << endl;

                
                    bool keep = false;
                    int kmer_count = 0;
                    
                    // Iterate over all k-mers
                    while (iterator.hasNext()) {
                        // Extract key from k-mer
                        iterator.operator()(key);
//                        cout << *((uint64_t*)key) << ": " << iterator.getSubstring(iterator.getPos(), 45);
                        // Do some classification magic
                        value = map.search(key);
                        
                        Utils::swapEndianess((uint8_t*)&value, 3);
                        int t = value;
                        
                        // Shift hit pattern by one
                        Utils::shift(hit_pattern, 4, 1);
                        
                        if (value) {
                            taxa.insert(t);
                            taxon_counter[t]++;
                            hit_count++;
                            //value = taxonomy->getNCBI(value);
                            keep = true;
                            found_count++;
                            //Set rightmost bit to one if value was found
                            hit_pattern[3] |= 1u;
                        } else {
                            unfound_count++;
                        }
                        
                        ++read_kmer_count;
                        ++kmer_count;
                        // For variant infering
                        if (kmer_count >= 32) {
                            vi.add(*((uint32_t *) hit_pattern), kmer_count-31);
                        }
                    }
    
                    // REsolve taxon id
                    int max = 0;
                    int tax_id = 0;
                    for (auto t : taxa) {
                        auto node = taxonomy->getNode(t);
                        //cout << node->id << ", ncbi: " << taxonomy.getNCBI(node->id) << ": " << node->level << " count: " << taxon_counter[t] << " -> ";
                        
                        // node has no parent, so no parent id;
                        while (node->id != 1 && (node = node->parent)->id != 1)
                            taxon_counter[t] += taxon_counter[node->id];
                        if (taxon_counter[t] > max) {
                            max = taxon_counter[t];
                            tax_id = t;
                        }
                        if (taxon_counter[t] == max) tax_id = 0;
                        //cout << taxon_counter[t] << endl;
                    }
                    if (!tax_id) {
                        for (auto t : taxa) {
                            if (taxon_counter[t] == max) {
                                if (!tax_id) tax_id = t;
                                else tax_id = taxonomy->lca(tax_id, t);
                            }
                        }
                    }
                    
//                    cout << "taxid: " << tax_id << endl;
//                    cout << "ncbi: " << taxonomy.getNCBI(tax_id) << endl;
//
//                    string stop;
//                    cin >> stop;

                    if (keep) {
                        std::string ratio = Utils::to_string_precision((double)hit_count / read_kmer_count, 3);
#pragma omp critical(output_classified)
                        {
                            classified++;
                            if (!vi.getVariants().empty()) {
                                for (auto var : vi.getVariants()) {
                                    if (var < 0) {
                                        cout << var << endl;
                                        exit(0);
                                    }
                                    oss_classified << thread_id;
                                    oss_classified << '_';
                                    oss_classified << read_id;
                                    oss_classified << '\t';
                                    oss_classified << (tax_id == 1 ? "root" : taxonomy->getNode(tax_id)->parent->name) << ";" << taxonomy->getOriginal(tax_id);
                                    oss_classified << '\t';
                                    oss_classified << iterator.getSubstring(var - 1 - (vmer_length / 2), vmer_length);
                                    oss_classified << '\t' << ratio << '\n'; // << '\t' << read_kmer_count << '\n';
                                }
                            } else {
                                oss_classified << thread_id;
                                oss_classified << '_';
                                oss_classified << read_id;
                                oss_classified << '\t';
                                oss_classified << (tax_id == 1 ? "root" : taxonomy->getNode(tax_id)->parent->name) << ";" << taxonomy->getOriginal(tax_id);
                                oss_classified << '\t' << '\t' << ratio << '\n'; // << hit_count << '\t' << read_kmer_count << '\n';
                            }
                        }
                    } else {
                        unclassified++;
#pragma omp critical(output_unclassified)
                        oss_unclassified << record.to_string();
                    }
                }
    
                (*os_classified) << oss_classified.str();
                (*os_unclassified) << oss_unclassified.str();
    
            }
            global_read_num += read_id;
            global_classified += classified;
            global_unclassified += unclassified;
            
            cout << "found count: " << found_count << endl;
            cout << "unfound count: " << unfound_count << endl;
        
            delete[] key;
            delete[] value_cache;
        } // parallel end
        delete is;
        delete[] shape;
        delete os_classified;
        delete os_unclassified;
        delete taxonomy;
    
        cout << "processed " << global_read_num << " reads" << endl;
        cout << "classified " << global_classified << endl;
        cout << "unclassified " << global_unclassified << endl;
        
        classify_bm.stop();
        classify_bm.printResults();
    }
    
    template<int KeyBits, int ValueBits, class T>
    static void runClassifyWorkerOMP(ClassifyOptionsContainer &options) {
        // Initialize Map
        BHashMap<KeyBits,ValueBits,T> map;
        map.load(options.meta_db.path + "/index.db", options.meta_db.path + "/index.meta");
        
        for (auto file : options.reads) {
            cout << "runClassify for: " << file << endl;
            runClassifyWorkerOMP(options, map, file);
        }
    }
    
    template<int KeyBits, int ValueBits, class T>
    static void runClassifyWorker(ClassifyOptionsContainer &options) {
        
        BHashMap<KeyBits, ValueBits, T> map (options.meta_db.capacity, options.meta_db.load_factor, options.meta_db.growth_factor);
        map.load(options.meta_db.path + "/kmer.db", options.meta_db.path + "/meta.txt");
        
        FastaReader reader(options.reads[0]);
        
        std::vector<std::thread> threads;
        
        if (options.threads == 1) {
            runClassifyThreadWorker<KeyBits, ValueBits, T>(map, reader, options);
        } else if (options.threads > 1) {
            for (int t = 0; t < options.threads; t++) {
                string name = std::to_string(t);
                threads.emplace_back([&]() {
                    runClassifyThreadWorker<KeyBits, ValueBits, T>(std::ref(map), std::ref(reader), std::ref(options));
                });
            }
            
            for (auto &t : threads)
                t.join();
        }
    }


    static void runBuildMarkerGenes(BuildOptionsContainer &options) {

    }

    template<int KeyBits, int ValueBits, class T>
    static void runBuildMarkerGenesWorkerOMP(BuildOptionsContainer &options) {
        ds::Benchmark build_bm("Build process");
        build_bm.start();

        // Shape information
        cout << "shape_path: " << options.db + MetaDataDB::SHAPE_FILE << endl;
        string shape_str = MetaDataDB::loadShape(options.db + MetaDataDB::SHAPE_FILE);
        const size_t shape_length = shape_str.length();
        bool * shape = getShape(shape_str);

        const int k = KeyBits/2;
        const int key_bytes = (2*k + 8 - 1) /  8;
        const int value_bytes = (ValueBits + 8 - 1) /  8;

        static size_t block_size = 200000;

        // Load Taxonomy
        TaxonomyInterface *taxonomy;
        if (options.taxonomy.compare("NCBI") == 0) {
            taxonomy = new NCBITaxonomy();
            taxonomy->loadCustomNodes(options.db + MetaDataDB::TAX_NODES_FILE);
            taxonomy->loadCustomNames(options.db + MetaDataDB::TAX_NAMES_FILE);
        } else if (options.taxonomy.compare("GTDB") == 0) {
            taxonomy = new GTDBTaxonomy();
            taxonomy->loadCustomNodes(options.db + MetaDataDB::TAX_NODES_FILE);
            taxonomy->loadCustomNames(options.db + MetaDataDB::TAX_NAMES_FILE);
        } else {
            std::cerr << options.taxonomy << " is no valid taxonomy. Valid taxonomies are NCBI and GTDB." << endl;
            exit(0);
        }

        // Initialize Map
        cout << "init capa: " << options.initial_capacity << endl;
        BHashMap<KeyBits,ValueBits,T> map( options.initial_capacity , 0.8 , 1.5);
        cout << "Initial size in MB: " << map.getMemoryInMB() << endl;
        cout << "files: " << options.reference.size() << endl;
        cout << "shape: " << shape_str << endl;
        cout << "shape_length: " << shape_length << endl;

        // Input stream
        std::istream* is = nullptr;
        // FIX set IS BEFORE PARALLEL
        is = new std::ifstream(options.reference[0]);
        if (!is) {
            cerr << "There is no such file: " << options.reference[0] << endl;
            exit(0);
        }

        //stats
        int record_count = 0;
        int skipped_count = 0;


        std::cout << "Number of threads: " << omp_get_thread_num() << std::endl;
        static int num_threads = 2;
        omp_set_num_threads(num_threads);

        uint8_t* key_temp = nullptr;
        uint8_t* value_temp = nullptr;

        //Multithreading
#pragma omp parallel
        {
#pragma omp critical(print_thread_start)
            std::cout << "start: " << omp_get_thread_num() << std::endl;

            BufferedFastxReader reader = BufferedFastxReader();
            FastxRecord record = FastxRecord();

            static thread_local uint8_t *key = new uint8_t[key_bytes];
            static thread_local uint8_t *key_rc = new uint8_t[key_bytes];
            static thread_local uint64_t value = 1;
            static thread_local uint64_t tax_id = 0;

            memset(key, 0, key_bytes);
            memset(key_rc, 0, key_bytes);

            // Iterator and putter
            static thread_local SpacedKmerIterator iterator(k, shape, shape_length);
            static thread_local KmerPutter<KeyBits, ValueBits, T> putter = KmerPutter<KeyBits, ValueBits, T>();

            // Initialize Putter
            putter.setMap(&map);
            putter.setTaxonomy(taxonomy);

            //Init KmerBuffer
            //KmerBuffer buffer(3125000, key_bytes, value_bytes);

            // stats
            uint64_t kmer_count = 0;
            uint64_t rand_calc = 0;

            // Iterate over all input files
            for (const auto &file : options.reference) {

                // Outer loop to load data in blocks (blocksize specified in number of records)
                while (true) {
                    bool ok = false;

#pragma omp critical(reader)
                    ok = reader.LoadBlock(*is, block_size);
                    if (!ok) break;

                    // Read records from datablock
                    while (true) {
                        auto valid_fragment = reader.NextSequence(record);
                        if (!valid_fragment) break;
                        iterator.setRecord(record);
                        cerr << record.sequence.length() << endl;

                        // count processed records
                        if ((record_count % 10) == 0)
                            cout << "processed record: " << record_count << "(of which " << skipped_count
                                 << " were skipped)" << endl;

                        #pragma omp atomic update
                        ++record_count;

                        // Extract taxonomic identifier from sequence header (has to be a ncbi identifier e.g.: >813)
                        value = stoll(record.header.substr(1));

//                        tax_id = options.taxonomy.compare("GTDB") == 0 ? taxonomy->getCustom(value) : taxonomy->getCustom(record.header.substr(1));

                        // depends on header prep of references.
                        //tax_id = taxonomy->getCustom(value);
                        tax_id = value;
                        if (tax_id == -1) {
                            cerr << "unknown taxid: " << tax_id << " for value " << value << endl;
                            exit(0);
                            continue;
                        }

                        // Skip if ncbi id is unknown to local taxonomy subset
                        if (!taxonomy->getNode(tax_id)) {
                            cerr << "local taxonomy subset has no node (entry) for NCBI identifier: " << tax_id << endl;
                            cerr << "Skip read (please refer to the option -w/--write_unknown to extract skipped reads)." << endl;
                            skipped_count++;
                            exit(0);
                            continue;
                        }

                        // Iterate over all k-mers
                        while (iterator.hasNext()) {
                            // Extract key from k-mer
                            iterator.operator()(key);
                            kmer_count++;

                            //buffer.pushKV(key, ((uint8_t*)&value));

                            //#pragma omp critical(putter)
                            putter.operator()(key, tax_id);

                            if (!map.search(key)) {
                                cout << "fail with: " << endl;
                                cout << record.to_string() << endl;
                                exit(0);
                            }
                        }
                    }
//                    while (buffer.next()) {
//                        //(key_temp = buffer.popKey();
//                        //rand_calc += (int) key_temp[0] + (int) key_temp[1] + (int) key_temp[2];
//                        #pragma omp critical(putter)
//                        putter.operator()(buffer.popKey(), tax_id);
//                    }
//                    buffer.reset();
                }
                cout << "kmercount: " << kmer_count << endl;
                cout << "randcalc: " << rand_calc << endl;
                break; //FIX CODE
            }

            delete[] key;
            delete[] key_rc;
        }

        build_bm.stop();
        build_bm.printResults();

        // exclude while testing key extracting methods etc.
        map.printVars();
        map.printStats();
        cout << "save db to " << (options.db + "/index.db") << endl;
        map.save(options.db + "/index");
        cout << "done saving." << endl;

        delete is;
        delete[] shape;
        delete taxonomy;

        cout << "validate: " << options.validate << endl;
        if (options.validate) {
            cout << "validate" << endl;
            runVerify(map, options);
        }
    }

    static void runClassify(ClassifyOptionsContainer &options) {
        
        if (options.meta_db.key_bits == 44 && options.meta_db.value_bits == 20)
            runClassifyWorker<44,20, CustomHash>(options);
        if (options.meta_db.key_bits == 46 && options.meta_db.value_bits == 18)
            runClassifyWorker<46,18, CustomHash>(options);
        if (options.meta_db.key_bits == 48 && options.meta_db.value_bits == 24)
            runClassifyWorker<48,24, CustomHash>(options);
        if (options.meta_db.key_bits == 48 && options.meta_db.value_bits == 24)
            runClassifyWorker<48,24, CustomHash>(options);
        
    }
    
    static void runClassifyOMP(ClassifyOptionsContainer &options) {
        cout << "path: " << options.meta_db.path << endl;
        cout << "key_bits: " << options.meta_db.key_bits << endl;
        cout << "value_bits: " << options.meta_db.value_bits << endl;
        
        if (options.meta_db.key_bits == 44 && options.meta_db.value_bits == 20)
            runClassifyWorkerOMP<44,20, CustomHash>(options);
    }
    
    
    
    static void runNewShape() {
    
    }
    
    static void runDBStats() {
    
    }
};


//#pragma clang diagnostic pop