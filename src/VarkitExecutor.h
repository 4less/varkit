//
// Created by joachim on 01/07/2020.
//

#ifndef VARKIT_VARKITEXECUTOR_H
#define VARKIT_VARKITEXECUTOR_H

#include <iostream>
#include <FastaReader.h>
#include <thread>
#include <SpacedKmerIterator.h>
#include "io/MetaDataDB.h"
#include "utils/ShapeUtils.h"
#include "kmer_processing/KmerProcessor.h"
#include "OptionsContainer.h"

using namespace std;

namespace VarkitExecutor {
    
    static const string PATTERN_MAP_PATH = "/shape/pattern_map.csv";
    static const string SHAPE_PATH = "/shape/shape";
    static const string TAXONOMY_NODES_PATH = "/taxonomy/nodes.dmp";
    static const string TAXONOMY_NAMES_PATH = "/taxonomy/names.dmp";
    static const string DB_PATH = "/kmer.db";
    static const string META_PATH = "/meta.txt";
    
    
    static bool * getShape(string s) {
        bool * shape = new bool[s.length()];
        const char * c = s.c_str();
        for (int i = 0; i < s.length(); i++) {
            shape[i] = c[i] == '_';
        }
        return shape;
    }
    
    template<int KeyBits, int ValueBits, class T>
    static void runBuildWorker(BuildOptionsContainer &options) {
        Benchmark build_bm("Build process");
        build_bm.start();
        
        // Shape information
        cout << "shape_path: " << options.db + MetaDataDB::SHAPE_FILE << endl;
        string shape_str = MetaDataDB::loadShape(options.db + MetaDataDB::SHAPE_FILE);
        const size_t shape_length = shape_str.length();
        bool * shape = getShape(shape_str);
        
    
        const int k = KeyBits/2;
        const int key_bytes = (2*k + 8 - 1) /  8;
        const int value_bytes = (ValueBits + 8 - 1) /  8;
        
        static thread_local uint8_t * key = new uint8_t[key_bytes];
        static thread_local uint64_t value;
        
        memset(key, 0, key_bytes);
    
        SpacedKmerIterator iterator(k, shape, shape_length);
        KmerPutter<KeyBits, ValueBits, T> putter = KmerPutter<KeyBits, ValueBits, T>();
    
        FastaRecord record;
        
        BHashMap<KeyBits,ValueBits,T> map( options.initial_capacity , 0.8 , 1.5);
        
        NCBITaxonomy taxonomy = NCBITaxonomy();
        taxonomy.loadCustomNodes(options.db + MetaDataDB::TAX_NODES_FILE);
        taxonomy.loadCustomNames(options.db + MetaDataDB::TAX_NAMES_FILE);
        
        putter.setMap(&map);
        putter.setTaxonomy(&taxonomy);
        
        cout << "files: " << options.reference.size() << endl;
        cout << "shape: " << shape_str << endl;
        cout << "shape_length: " << shape_length << endl;
        
        int record_count = 0;
        for (auto file : options.reference) {
            FastaReader reader(file);
            
            
            while (reader.hasNext()) {
                if ((++record_count % 10) == 0)
                    cout << "processed record: " << record_count << endl;
                record = reader.next();
                iterator.setFastaRecord(record);
                
                cout << record.header << endl;
                
                cout << taxonomy.getNode(taxonomy.getCustom(stoll(record.header.substr(1))))->id << endl;
                
                while (iterator.hasNext()) {
                    iterator.operator()(key);
                    value = stoll(record.header.substr(1));
                    putter.operator()(key, value);
                }
            }
        }
        
        build_bm.stop();
        build_bm.printResults();
        putter.put.printResults();
        putter.swap.printResults();
        putter.routine.printResults();
        map.printVars();
        map.printStats();
        
        
        map.save(options.db + "index.db");
        
        
    }
    
    static void runBuild(BuildOptionsContainer &options) {
        int key_bits = 44;
        int value_bits = 20;
        
        if (key_bits == 44 && value_bits == 20)
            runBuildWorker<44,20, CustomHash>(options);
    }
    
    
    template<int KeyBits, int ValueBits, class T>
    static void runClassifyThreadWorker(BHashMap<KeyBits, ValueBits, T> &map, FastaReader &reader, ClassifyOptionsContainer &options) {
        bool shape_bool[options.meta_db.shape.length()];
        memset(shape_bool, 0, options.meta_db.shape.length());
        shape_utils::getBool(options.meta_db.shape, shape_bool);
        SpacedKmerIterator iterator(options.meta_db.key_bits/2, shape_bool, options.meta_db.shape.length());
        
        
    }
    
    template<int KeyBits, int ValueBits, class T>
    static void runClassifyWorker(ClassifyOptionsContainer &options) {
        
        BHashMap<KeyBits, ValueBits, T> map (options.meta_db.capacity, options.meta_db.load_factor, options.meta_db.growth_factor);
        map.load(options.meta_db.path + "/kmer.db", options.meta_db.path + "/meta.txt");
    

        FastaReader reader(options.read);
    
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
    
    static void runNewShape() {
    
    }
    
    static void runDBStats() {
    
    }
};


#endif //VARKIT_VARKITEXECUTOR_H
