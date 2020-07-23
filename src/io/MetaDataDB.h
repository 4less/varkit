//
// Created by joachim on 01/07/2020.
//

#ifndef VARKIT_METADATADB_H
#define VARKIT_METADATADB_H

#include <iostream>
#include <fstream>
#include <iostream>
#include <regex>

using namespace std;

struct MetaDataDB {
public:
    const string path;
    const int key_bits;
    const int value_bits;
    int key_bytes_;
    int value_bytes_;
    const int size_;
    const int capacity;
    const int load_factor;
    const int growth_factor;
    const string shape;
    
    inline const static string SHAPE_FILE = "/shape/shape";;
    inline const static string SNP_PATTERN_FILE = "/shape/map.csv";
    inline const static string TAX_NODES_FILE = "/taxonomy/nodes.dmp";
    inline const static string TAX_NAMES_FILE = "/taxonomy/names.dmp";
    
    MetaDataDB(string path, int key_bits, int value_bits, int capacity, int size, double load_factor, double growth_factor, string shape);
    
    static string loadShape(string path) {
        ifstream shapein;
        string line;
        string shape;
        
        shapein.open(path );
        if (shapein.is_open()) {
            getline(shapein, line);
            if (std::regex_match (line, std::regex("[X|_]+") ))
                shape=line;
        }
        shapein.close();
        return shape;
    }
};

MetaDataDB loadMetaDataDB(const string& filename);


#endif //VARKIT_METADATADB_H
