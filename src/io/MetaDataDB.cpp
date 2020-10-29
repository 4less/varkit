//
// Created by joachim on 01/07/2020.
//

#include "MetaDataDB.h"
#include <iostream>

using namespace std;

MetaDataDB::MetaDataDB(string path, int key_bits, int value_bits, int capacity, int size, double load_factor, double growth_factor, string shape, std::string taxonomy)  : path(path), key_bits(key_bits), value_bits(value_bits), capacity(capacity), size_(size), load_factor(load_factor), growth_factor(growth_factor), shape(shape), taxonomy(taxonomy) {}

MetaDataDB loadMetaDataDB(const string& path) {
    const string delimiter = "=";
    string line;
    
    int key_bits, value_bits, capacity, size;
    double load_factor,growth_factor;
    string shape, taxonomy;
    
    // read meta file_
    ifstream metain;
    metain.open(path + "/index.meta");
    
    if (!metain) {
        cerr << "Given database folder \"" << path << "\" has no file called index.meta" << endl;
        exit(0);
    }
    if (metain.is_open()) {
        while (getline (metain, line)) {
            if (line.substr(0,8) == "key_bits") {
                string token = line.substr(line.find(delimiter) + delimiter.length(), line.length());
                key_bits = stoi(token);
            }
            else if (line.substr(0,10) == "value_bits") {
                string token = line.substr(line.find(delimiter) + delimiter.length(), line.length());
                value_bits = stoi(token);
            }
            else if (line.substr(0,8) == "capacity") {
                string token = line.substr(line.find(delimiter) + delimiter.length(), line.length());
                capacity = stoll(token);
            }
            else if (line.substr(0,4) == "size") {
                string token = line.substr(line.find(delimiter) + delimiter.length(), line.length());
                size = stoll(token);
            }
            else if (line.substr(0,11) == "load_factor") {
                string token = line.substr(line.find(delimiter) + delimiter.length(), line.length());
                load_factor = stod(token);
            }
            else if (line.substr(0,13) == "growth_factor") {
                string token = line.substr(line.find(delimiter) + delimiter.length(), line.length());
                growth_factor = stod(token);
            }
            else if (line.substr(0,8) == "taxonomy") {
                string token = line.substr(line.find(delimiter) + delimiter.length(), line.length());
                taxonomy = token;
            }
        }
    }
    metain.close();
    
    return MetaDataDB(path, key_bits, value_bits, capacity, size, load_factor, growth_factor, shape, taxonomy);
}