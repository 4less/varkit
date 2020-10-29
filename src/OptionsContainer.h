//
// Created by joachim on 06/07/2020.
//
#pragma once

#include "io/MetaDataDB.h"

class ClassifyOptionsContainer {
public:
    const MetaDataDB meta_db;
    const int threads;
//    const string read;
    const std::vector<std::string> reads;
    const string output;
    const size_t vmer_length;
    
    ClassifyOptionsContainer(MetaDataDB &meta_db, std::vector<std::string> reads, int threads, string output, size_t vmer_length) :
            meta_db(meta_db),
            reads(reads),
            threads(threads),
            output(output),
            vmer_length(vmer_length) {};
};

class BuildOptionsContainer {
public:
    const int threads;
    const string db;
    const vector<string> reference;
    const int64_t initial_capacity;
    const string taxonomy;
    const bool validate;
    
    BuildOptionsContainer(int threads, string db, vector<string> reference, int64_t initial_capacity, string taxonomy, bool validate) :
            reference(reference),
            db(db),
            threads(threads),
            initial_capacity(initial_capacity),
            taxonomy(taxonomy),
            validate(validate) {};
};