//
// Created by joachim on 22/09/2020.
//
#pragma once

#include "KmerBuffer.h"
//#include <zstd.h>
#include "zstr.hpp"
#include "strict_fstream.hpp"
#include <zlib.h>

static inline uint8_t getBitFromBaseTest(char base) {
    switch (base) {
        case 'A':
            return (0);
        case 'C':
            return (1);
        case 'G':
            return (2);
        case 'T':
            return (3);
    }
    return (-1);
}

void test_key() {
    char a = 'A';
    const char ca = 'C';
    string test = "ACGT";
    
    cout << "1: " << (int) getBitFromBaseTest('G') << endl;
    cout << "2: "  << (int) getBitFromBaseTest(a) << endl;
    cout << "3: "  << (int) getBitFromBaseTest(ca) << endl;
    cout << "3: "  << (int) getBitFromBaseTest(*(test.c_str())) << endl;
    cout << "3: "  << (int) getBitFromBaseTest(*(test.c_str() + 1)) << endl;
    cout << "3: "  << (int) getBitFromBaseTest(*(test.c_str() + 2)) << endl;
    cout << "3: "  << (int) getBitFromBaseTest(*(test.c_str() + 3)) << endl;
}

static inline void test_kmer_buffer() {
    KmerBuffer buffer(100, 8, 8);
    
    for (int i = 1; i <= 105; i++) {
        uint64_t key = i;
        uint64_t value = i + 1;

        buffer.pushKV((uint8_t*) &key, (uint8_t*) &value);
    }

    cout << "retrieve" << endl;


    uint64_t key, value;
    while (buffer.next()) {
        key = *((uint64_t *) buffer.popKey());
        value = *((uint64_t *) buffer.popValue());
        cout << key << " : " << value << endl;
    }
    
    buffer.reset();
    
    for (int i = 1; i <= 105; i++) {
        uint64_t key = i;
        uint64_t value = i + 1;
        
        buffer.pushKV((uint8_t*) &key, (uint8_t*) &value);
    }
    
    cout << "retrieve" << endl;
    
    while (buffer.next()) {
        key = *((uint64_t *) buffer.popKey());
        value = *((uint64_t *) buffer.popValue());
        cout << key << " : " << value << endl;
    }
}

static void gzip_test(const char* in, const char* out) {
    
    std::istream* is = new ifstream(in);
    std::ostream* os = new zstr::ofstream(out, ios::app);

    const std::streamsize buff_size = 1 << 16;
    char * buff = new char[buff_size];

    while (true) {
        is->read(buff, buff_size);
        std::streamsize cnt = is->gcount();
        if (cnt == 0) break;
        
        os->write(buff, cnt);
    }

    delete[] buff;
    delete is;
    delete os;
}

static void read_auto(const char* in) {
    
    std::istream* is = new zstr::ifstream(in);
//
//    const std::streamsize buff_size = 1 << 16;
//    const std::streamsize buff_size = 1 << 5;
//    cout << "buffsize: " << buff_size << endl;
//    char * buff = new char[buff_size];
//    memset(buff, 0 , buff_size);
    string line;
    
    while (getline(*is, line)) {
        cout << line << "\n";
    }
    delete is;
}

static void read_unzipped(const char* in) {
    
    std::istream* is = new ifstream(in);
    
    const std::streamsize buff_size = 1 << 16;
    char * buff = new char[buff_size];
    memset(buff, 0 , buff_size);
    
    while (true) {
        is->read(buff, buff_size);
        std::streamsize cnt = is->gcount();
        if (cnt == 0) break;
    
        cout << buff;
    }
    
    delete[] buff;
    delete is;
}

static void read_getline(const char* in) {
    
    std::istream* is = new ifstream(in);
    
    std::string line = "";
    
    while (getline(*is, line)) {
        cout << line << endl;
    }
    
    delete is;
}

static void gunzip_test(const char* in, const char* out) {
    
    std::istream* is = new zstr::ifstream(in);
    std::ostream* os = new ofstream(out);
    
    const std::streamsize buff_size = 1 << 16;
    char * buff = new char[buff_size];
    
    while (true) {
        is->read(buff, buff_size);
        std::streamsize cnt = is->gcount();
        if (cnt == 0) break;
        os->write(buff, cnt);
    }
    
    delete[] buff;
    delete is;
    delete os;
}