//
// Created by joachim on 08/06/2020.
//
#pragma once

#include <Taxonomy.h>

class KmerProcessor {
public:
    virtual uint64_t operator () (uint8_t * key, uint8_t * value) = 0;
    virtual uint64_t operator () (uint8_t * key, uint64_t &value) = 0;
};


template <int KeyBits, int ValueBits, class T>
class KmerPutter : public KmerProcessor {
private:
    BHashMap<KeyBits,ValueBits,T> * map_;
    TaxonomyInterface * taxonomy_ = nullptr;
    const int value_bytes = (ValueBits + 8 - 1) /  8;
    
public:
    Benchmark put = Benchmark("put");
    Benchmark swap = Benchmark("swap");
    Benchmark routine = Benchmark("routine");
    
    void setMap(BHashMap<KeyBits, ValueBits, T> * map) {
        map_ = map;
    }

    void setTaxonomy(TaxonomyInterface * taxonomy) {
        taxonomy_ = taxonomy;
    }
    
    void Debug() {
        uint8_t key[6];
        memset(key, 0, 6);
        string stop;
        cin >> stop;
        while(stop != "stop") {
            cout << "lookup: " << stop << endl;
            KmerUtils::dnaToBytes(stop, key);
            cout << *((uint32_t*) key) << endl;
            cout << KmerUtils::bytesToDNA(key, 22) << endl;
            cout << KmerUtils::dnaToBitString(stop) << endl;
            cout << map_->search(key) << endl;
            cin >> stop;
        }
    }
    
    void Debug(uint8_t * k) {
        cout << "test: " << KmerUtils::bytesToDNA(k, KeyBits/2) << endl;
        cout << *((uint64_t*) k) << endl;
        Debug();
    }
    
    uint64_t operator () (uint8_t * key, uint8_t * value) override {
        int id = 0;
        
        if ((id = map_->search(key)) != 0 && taxonomy_) {
            map_->put(key, value);
        }
        else {
            map_->put(key, value);
        }
    };
    
    uint64_t operator () (uint8_t * key, uint64_t &value) override {
        routine.start();
        static thread_local int id = 0;
        static thread_local uint8_t value_c[8];
        static thread_local uint8_t value_c2[8];
        static thread_local int count = 0;
        static thread_local uint8_t fixed[] = { 0x00, 0x00, 0x01 };
        
        if ((++count % 5000000) == 0) {
            cout << "size: " << map_->getElements() << endl;
            cout << "processed: " << count << endl;
        }
        
        //cout << "search: " << map_->search(key) << " " << KmerUtils::bytesToDNA(key, KeyBits/2) << endl;
        
        if ((id = map_->search(key)) != 0 && taxonomy_) {
            swap.start();
            *((uint64_t*)value_c2) = id;
            Utils::swapEndianess(value_c2, value_bytes);
            *((uint64_t*)value_c) = taxonomy_->lca(value, *((uint64_t*)value_c2));
            cout << "lca: " << *((uint64_t*)value_c) << endl;
            
            Utils::swapEndianess(value_c, value_bytes);
            swap.stop();
            
            put.start();
            map_->put(key, value_c+5);
            put.stop();
        }
        else {
            swap.start();
            *((uint64_t*)value_c) = value;
            Utils::swapEndianess(value_c, value_bytes);
            swap.stop();
            
            put.start();
            map_->put(key, value_c+5);
            put.stop();
        }
        routine.stop();
    };
};

