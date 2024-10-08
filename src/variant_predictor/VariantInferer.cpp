//
// Created by joachim on 15/06/2020.
//

#include "VariantInferer.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include "Utils.h"
#include "BitUtils.h"

VariantInferer::VariantInferer(string db_path, uint32_t read_length) : read_length_(read_length) {
    ifstream db_in;
    int line_count = 0;
    db_in.open(db_path);
    if (!db_in) {
        cerr << "Could not find pattern map at " << db_path << endl;
        exit(123);
    } else {
        std::string line;
        while (getline(db_in, line)) line_count++;
    }
    db_in.close();
    
    cout << "linecount: " << line_count << endl;
    
    double load_factor = 0.3;
    uint64_t capacity = line_count / load_factor + 1;
    cout << "capacity: " << ++capacity << endl;
    
    this->pattern_map_2_ = new BHashMap<32,64,SuperFastHash>(capacity, load_factor, 1.5);
    cout << "capamap: " << pattern_map_2_->getCapacity() << endl;
    clear();
    loadDB(db_path);
    
}

VariantInferer::~VariantInferer() {
//    cout << omp_get_thread_num() << " delete VariantInferer" << endl;
//    if (pattern_map_2_) {
//        cout << "delete VariantInferer";
//        delete pattern_map_2_;
//    }
}

void VariantInferer::loadDB(string db_path) {
    
    string key;
    uint32_t keyInt;
    string value;
    string delimiter = " ";
    size_t del_pos;
    string line;
    
    cout << "load db: " << db_path << endl;
    
    // read db file
    ifstream db_in;
    db_in.open(db_path);
    if (db_in.is_open()) {
        
        while (getline(db_in, line)) {
            del_pos = line.find(delimiter);
            key = line.substr(0, del_pos);
            value = line.substr(del_pos + delimiter.length(), line.length());
            
            if (value == "")
                continue;
            
            memset(value_cache, 255, 8);//4);
            int cache_pos = 0;
            for (auto rel_pos : split(value))
                value_cache[cache_pos++] = (uint8_t) rel_pos;
        
            keyInt = keyToInt(key);
            *((uint32_t *) key_cache) = keyInt;
        
            pattern_map_2_->put(key_cache, value_cache);
        }
    }
    pattern_map_2_->printStats();
    cout << "size: " << pattern_map_2_->getElements() << endl;
    
    cout << "save out" << endl;
    pattern_map_2_->save("/home/joachim/CLionProjects/varkit/data/bactsub_db/shape/shape");
    
}

uint32_t VariantInferer::keyToInt(string key) {
    memset(byte_key,0,4);
    const char * key_char = key.c_str();
    
    for (int i = 0; i <  key.length(); i++) {
        if (key_char[i] == '_')
            byte_key[i / 8] |= (uint8_t) 1u << (7 - (i % 8));
    }
    //cout << "key: " << key << "   " << BitUtils::ToBitString(byte_key, 4) << endl;
    return *((uint32_t *) byte_key);
}

void VariantInferer::add(uint32_t key, uint32_t pos) {
//    static int counter = 0;
//    counter++;
//    if ((counter % 5000000) == 0)
//        cout << counter << ": " << key << "," << variants.size() << endl;
    uint64_t value = pattern_map_2_->search((uint8_t*)&key);
    
    if (value == 0) return;

    *(uint64_t*)value_cache = value;

    // changed from 4 to 8
    uint8_t vpos = 0;
    uint32_t rpos = 0;
    uint32_t qualifies = 0;
    for (int i = 0; i < 8 && (vpos = value_cache[i]) != 0xFF && (rpos = vpos+pos) < current_length_; i++) {
//        cout << "cuckoooo" << endl;
        variants.insert(rpos);
    }
    
}

vector<uint32_t> VariantInferer::split(const string &s) {
    vector<uint32_t> values;
    uint32_t token;
    
    string line = s;
    int lastpos = 0;
    int pos = 0;
    
    while ((pos = line.find(",", lastpos)) != string::npos) {
        std::cout << "stoi incoming..." << std::endl;
        token = stoi(line.substr(lastpos, pos-lastpos));
        std::cout << "after stoi..." << std::endl;
        values.push_back(token);
        lastpos = pos+1;
    }
    std::cout << "stoi incoming...2" << std::endl;
    token = stoi(line.substr(lastpos, line.length()-lastpos));
    std::cout << "after stoi...2" << std::endl;
    values.push_back(token);
    return values;
}

void VariantInferer::clear() {
    added_count_ = 0;
    variants.clear();
}



set<uint32_t> VariantInferer::getVariants() {
    return variants;
}

uint32_t VariantInferer::get(uint32_t key) {
    uint64_t value = pattern_map_2_->search(((uint8_t*)&key));
    cout << value << endl;
    return (uint32_t)(value << 32u);
}

void VariantInferer::setLength(uint32_t length) {
    current_length_ = length;
}

void VariantInferer::writeVarPos(string path) {
    ofstream varout;
    varout.open(path);
    if (!varout) {
        cout << "Cannot open file!" << endl;
        exit(1);
    }
    
    for (auto var : variants) {
        varout << var << "\n";
    }
    
    varout.close();
}

void VariantInferer::load(std::string path) {
    pattern_map_2_ = new BHashMap<32,64,SuperFastHash>();
    pattern_map_2_->load(path + ".db", path + ".meta");
}

void VariantInferer::setMap(BHashMap<32,64,SuperFastHash>* map) {
    pattern_map_2_ = map;
}