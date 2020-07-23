//
// Created by joachim on 08/07/2020.
//

#ifndef VARKIT_IOTEST_H
#define VARKIT_IOTEST_H


#include <iostream>
#include <cstring>
#include <fstream>
#include "Benchmark.h"
#include <sstream>

using namespace std;


namespace iotest {
    
    const unsigned long long int bufSize = 1<<24;
    const unsigned long long int smallBufSize = 1<<12;
    const unsigned int headBufSize = 1024;
    
    static void fopenTest(string &path) {
        uintmax_t char_counter = 0;
        FILE* f;
        string line;
        char* buf = new char[bufSize];
        
        if ((f = fopen(path.c_str(), "r"))) {
            while (fgets(buf, bufSize, f) != NULL) {
                buf[strlen(buf)-1] = '\0';
                line = string(buf);
                char_counter = char_counter + line.size() - 3;
            }
        }
        cout << "chars counted: " << char_counter << endl;
        delete[] buf;
    }
    
    static void fopenFastaTest(string &path) {
        uintmax_t char_counter = 0;
        FILE* f;
        int pos;
        char head[headBufSize];
        char* seq = new char[bufSize];
        string seq_s;
        uint32_t seq_pos = -1;
        
        if ((f = fopen(path.c_str(), "r"))) {
            while (fgets(head, sizeof(head), f) != NULL) {
                while (!feof(f) && (seq[++seq_pos] = fgetc(f)) != '>')
                    if (seq[seq_pos] == '\n') --seq_pos;
                seq[seq_pos-1] = '\0';
                char_counter += strlen(seq);
                seq_pos = -1;
            }
        }
        cout << "chars counted: " << char_counter << endl;
        delete[] seq;
    }
    
    static void inline doSth(int& counter, string& sequence) {
        for (int i = 0; i < sequence.length(); i++) {
            counter += sequence.c_str()[i];
        }
    }
    
    static void ifstreamFastaTest(string &path) {
        int counter = 0;
        
        const int reserve=200;
        ifstream f;
        //char* seq = new char[bufSize];
        char seq[smallBufSize];
        f.rdbuf()->pubsetbuf(seq, smallBufSize);
        f.open(path);
        
        string line;
        string sequence;
        sequence.reserve(reserve);
        uintmax_t char_counter = 0;
        
        
        getline(f, line);
        while (!f.eof()) {
            sequence = "";
            sequence.reserve(reserve);
            while (getline(f, line) && line[0] != '>') {
                sequence += line;
            }
            char_counter += sequence.length();
            
            doSth(counter, sequence);
        }

        cout << "chars counted: " << char_counter << endl;
        cout << "counter: " << counter << endl;
        //delete[] seq;
    }
    

    
    static void ifstreamTest(string &path) {
        ifstream f;
        f.open(path);
        string line;
        uintmax_t char_counter = 0;
        
        while (getline(f, line)) {
            char_counter = char_counter + line.size() - 3;
        }
        cout << "chars counted: " << char_counter << endl;
    }
    
    static void compareIO(string path) {
        Benchmark if_bm("ifstream");
        Benchmark fo_bm("fopen");
        

        fo_bm.start();
        fopenTest(path);
        fo_bm.stop();
    
    
        if_bm.start();
        ifstreamTest(path);
        if_bm.stop();
    
    
        if_bm.printResults();
        fo_bm.printResults();
    }
    
    static void compareFastaIO(string path) {
        Benchmark if_bm("ifstream");
        Benchmark fo_bm("fopen");
    
        fo_bm.start();
        fopenFastaTest(path);
        fo_bm.stop();
        
        if_bm.start();
        ifstreamFastaTest(path);
        if_bm.stop();
        

        
        
        if_bm.printResults();
        fo_bm.printResults();
    }
};


#endif //VARKIT_IOTEST_H
