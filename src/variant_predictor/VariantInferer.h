//
// Created by joachim on 15/06/2020.
//

#pragma once

#include "robin_map.h"
#include <string>
#include <vector>
#include <set>
#include <Benchmark.h>
//#include "BHashMap.h"

using namespace std;
using namespace tsl;

class VariantInferer {
    BHashMap<32,64,SuperFastHash> *pattern_map_2_;
    bool del = false;
    uint32_t read_length_;
    
    int current_length_;
    set<uint32_t> variants;
    
    uint8_t byte_key[4];
    uint8_t value_cache[8];
    uint8_t key_cache[4];
    
    uint32_t added_count_;
    
    vector<uint32_t> split(const string &s);
    void loadDB(string db_path);
    uint32_t keyToInt(string key);
public:
    
    ds::Benchmark bm = ds::Benchmark("VariantI");
    
    VariantInferer () {};
    VariantInferer (string db_path, uint32_t read_length);
    ~VariantInferer();
    void clear();
    void setLength(uint32_t length);
    set<uint32_t> getVariants();
    void add(uint32_t key, uint32_t pos);
    uint32_t get(uint32_t key);
    void writeVarPos(string path);
    
    void load(std::string path);
    
    void setMap(BHashMap<32, 64, SuperFastHash> *map);
};


