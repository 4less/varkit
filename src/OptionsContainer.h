//
// Created by joachim on 06/07/2020.
//
#pragma once

#include "io/MetaDataDB.h"

class ClassifyOptionsContainer {
public:
    const MetaDataDB meta_db;
    const int threads;
    const string read;
    
    ClassifyOptionsContainer(MetaDataDB &meta_db, string read, int threads) :
            meta_db(meta_db),
            read(read),
            threads(threads) {};
};

class BuildOptionsContainer {
public:
    const int threads;
    const string db;
    const vector<string> reference;
    const int initial_capacity;
    
    BuildOptionsContainer(int threads, string db, vector<string> reference, int initial_capacity) :
            reference(reference),
            db(db),
            threads(threads),
            initial_capacity(initial_capacity) {};
};
