//
// Created by joachim on 08/06/2020.
//
#pragma once

#include <Taxonomy.h>

class KmerProcessor {
public:
    virtual uint64_t operator () (uint8_t * key, uint8_t * value) = 0;
    virtual uint64_t operator () (uint8_t * key, uint64_t &value) = 0;
    virtual uint64_t operator () () = 0;
};


template <int KeyBits, int ValueBits, class T>
class KmerPutter : public KmerProcessor {
private:
    BHashMap<KeyBits,ValueBits,T>* map_ = nullptr;
    TaxonomyInterface * taxonomy_ = nullptr;
    const int value_bytes = (ValueBits + 8 - 1) /  8;
    int count = 0;
    
public:
    // remove after debugging
    uint8_t test_key[6];
    
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
        return 0L;
    };
    
    uint64_t operator () () override {
        cout << "bollocks" << endl;
        return 0L;
    }
    
    static bool equals(uint8_t* one, uint8_t* two, size_t byte) {
        for (int i = 0; i < byte; i++) {
            if (one[i] != two[i]) return false;
        }
        return true;
    }
    
    uint64_t operator () (uint8_t* key, uint64_t &value) override {
        static thread_local uint64_t id_cached = 0;
        static thread_local uint64_t value_cached = 0;
        
        static thread_local uint64_t id = 0;
        static thread_local uint8_t value_c[8];
        static thread_local uint8_t value_c2[8];
        static uint64_t count = 0;
        static int offset = 8 - value_bytes;
        
        count += 1;
    
        
        if ((count % 5000000) == 0) {
            cout << "processed: " << count << endl;
        }
        
        id = 1;
    
//        //DEBUGGING
//        if (equals(key, test_key, 6)) {
//            cout << "#find before" << endl;
//            cout << KmerUtils::bytesToDNA(key, 22) << endl;
//            cout << "value: " << value << endl;
//            uint64_t map_value = map_->search(key);
//            *((uint64_t*)value_c2) = map_value;
//            Utils::swapEndianess(value_c2, value_bytes);
//            cout << "map_value: " << *((uint64_t*)value_c2) << endl;
//        }
        
        if ((id = map_->search(key)) != 0) {
            // If the last value and the last retrieved value did not change then just put the key in
            if (value == value_cached && id == id_cached) {
                //cout << "cached... " << endl;
//                if (equals(key, test_key, 6)) {
//                    cout << "#find before put cached area" << endl;
//                    cout << KmerUtils::bytesToDNA(key, 22) << endl;
//                    cout << "value: " << value << endl;
//                    cout << "value_cached: " << value_cached << endl;
//                    uint64_t map_value = map_->search(key);
//                    *((uint64_t*)value_c2) = map_value;
//                    Utils::swapEndianess(value_c2, value_bytes);
//                    cout << "map_value: " << *((uint64_t*)value_c2) << endl;
//                }
                
                map_->put(key, (value_c+offset));
    
//                if (equals(key, test_key, 6)) {
//                    cout << "#find cached area" << endl;
//                    cout << KmerUtils::bytesToDNA(key, 22) << endl;
//                    cout << "value: " << value << endl;
//                    cout << "value_cached: " << value_cached << endl;
//                    uint64_t map_value = map_->search(key);
//                    *((uint64_t*)value_c2) = map_value;
//                    Utils::swapEndianess(value_c2, value_bytes);
//                    cout << "map_value: " << *((uint64_t*)value_c2) << endl;
//                }
                
                return 0L;
            }
            
            *((uint64_t*)value_c2) = id;
            Utils::swapEndianess(value_c2, value_bytes);
            if (value == *((uint64_t*)value_c2) || *((uint64_t*)value_c2) == 1) {
                //cout << "new_value: " << value << "; map_value: " <<  *((uint64_t*)value_c2) << endl;
    
//                if (equals(key, test_key, 6)) {
//                    cout << "#find mapval = 1 or mapval same area" << endl;
//                    cout << KmerUtils::bytesToDNA(key, 22) << endl;
//                    cout << "value: " << value << endl;
//                    uint64_t map_value = map_->search(key);
//                    *((uint64_t*)value_c2) = map_value;
//                    Utils::swapEndianess(value_c2, value_bytes);
//                    cout << "map_value: " << *((uint64_t*)value_c2) << endl;
//                }
                
                return 0L;
            }
            
            *((uint64_t*)value_c) = taxonomy_->lca(value, *((uint64_t*)value_c2));
//            if (*((uint64_t*)value_c) == 1) {
//                cout << "lca(" << value << "," << *((uint64_t*)value_c2) << "): " << *((uint64_t*)value_c) << endl;
//            }
            
//            if (lca_c < 1000 || value == 92 || *((uint64_t*)value_c2) == 92)
//                cout << "lca(" << value << "," << *((uint64_t*)value_c2) << "): " << *((uint64_t*)value_c) << endl;
            
            Utils::swapEndianess(value_c, value_bytes);
            map_->put(key, (value_c+offset));
            //map_->printMapPSL(map_->Hash(key), 5);
            //string stop;
            //std::cin >> stop;
            id_cached = id;
            value_cached = *((uint64_t*)value_c);
        }
        else {
            *((uint64_t*)value_c) = value;
            Utils::swapEndianess(value_c, value_bytes);
            int offset = 8-value_bytes;
            map_->put(key, (value_c+offset));
        }
    
//        //DEBUGGING
//        if (equals(key, test_key, 6)) {
//            cout << "#find after" << endl;
//            cout << KmerUtils::bytesToDNA(key, 22) << endl;
//            cout << "value: " << value << endl;
//            uint64_t map_value = map_->search(key);
//            *((uint64_t*)value_c2) = map_value;
//            Utils::swapEndianess(value_c2, value_bytes);
//            cout << "map_value: " << *((uint64_t*)value_c2) << endl;
//        }
        
        return 0L;
    };
};

