//
// Created by joachim on 01/07/2020.
//

#pragma once

#include <iostream>
#include <fstream>
#include <iostream>
#include <regex>

using namespace std;

struct MetaDataDB {
public:
    const string path = "";
    const int key_bits = -1;
    const int value_bits = -1;
    int key_bytes_ = -1;
    int value_bytes_ = -1;
    const int size_ = -1;
    const int capacity = -1;
    const int load_factor = -1;
    const int growth_factor = -1;
    const string shape = "";
    const string taxonomy = "";
    
    inline const static string SHAPE_FILE = "/shape/shape";;
    inline const static string SNP_PATTERN_FILE = "/shape/map.csv";
    inline const static string TAX_NODES_FILE = "/taxonomy/nodes.dmp";
    inline const static string TAX_NAMES_FILE = "/taxonomy/names.dmp";
    
    MetaDataDB(string path, int key_bits, int value_bits, int capacity, int size, double load_factor, double growth_factor, string shape, std::string taxonomy);
    
    static string loadShape(string path) {
        ifstream shapein;
        string line;
        string shape;
        
        shapein.open(path);
        if (shapein.is_open()) {
            getline(shapein, line);
            if (std::regex_match (line, std::regex("[X|_]+") ))
                shape=line;
        }
        shapein.close();
        return shape;
    }
    
    static bool* getShape(string s) {
        bool * shape = new bool[s.length()];
        const char * c = s.c_str();
        for (int i = 0; i < s.length(); i++) {
            shape[i] = c[i] == '_';
        }
        return shape;
    }
};

MetaDataDB loadMetaDataDB(const string& filename);

