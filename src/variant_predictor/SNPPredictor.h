//
// Created by fritsche on 11/05/2021.
//

#pragma once
#include "sparse_map.h"
#include "ShapeBuilder.h"
#include "data_structures.h"
#include <fstream>

#include <string>



class SNPPredictor {
public:
    bool* shape;

    size_t shape_size;
    size_t pattern_size;
    size_t limit;
    size_t lower_limit;

    tsl::sparse_map<uint32_t,ds::Mutations> *map;

    uint32_t read_length_;

    uint32_t pattern = 0;
    uint32_t cur_pos = 0;

    int current_length_;
    set<uint32_t> variants;

    ~SNPPredictor() {
        delete[] shape;
    }

    void Load(std::string path);
    void SetMap(tsl::sparse_map<uint32_t,ds::Mutations> &map);
    void Add(bool state);
};
