//
// Created by fritsche on 11/05/2021.
//

#include "SNPPredictor.h"

void SNPPredictor::Load(std::string path) {
    ifstream ifs(path, ifstream::in);

    std::string line;


    ifs.read((char *) &shape_size, sizeof(shape_size));
    ifs.read((char *) &pattern_size, sizeof(pattern_size));
    ifs.read((char *) &limit, sizeof(limit));
    ifs.read((char *) &lower_limit, sizeof(lower_limit));

    delete[] shape;
    shape = new bool[shape_size];
    ifs.read((char *) shape, shape_size * sizeof(bool));

    map = new tsl::sparse_map<uint32_t,ds::Mutations>();

    while (ifs.peek() != EOF) {
        int32_t pattern;
//        size_t value;
        ds::Mutations value;
        ifs.read((char*) &pattern, 4);
        ifs.read((char*) &value, sizeof(ds::Mutations));
        map->insert({ pattern, value });
    }

    ifs.close();
}

void SNPPredictor::Add(bool state) {
    pattern <<= 1;
    pattern |= ((int) state) << (32 - pattern_size);
    cur_pos++;

    if (cur_pos >= pattern_size) {
        auto res = map->find(pattern);


    }
}

void SNPPredictor::SetMap(tsl::sparse_map<uint32_t, ds::Mutations> &map) {
    this->map = &map;
}
