//
// Created by joachim on 02/07/2020.
//
#pragma once

#include <iostream>
#include <thread>
#include <mutex>
#include <vector>
#include "FastaReader.h"

//#include <seqan3/alphabet/nucleotide/dna4.hpp>
//#include <seqan3/io/sequence_file/input.hpp>
//#include <seqan3/io/sequence_file/output.hpp>
//#include <seqan3/std/filesystem>

double computation(double x) {
    return x*3.14;
}

double cached_computation(double x)
{
    std::cout << "x: " << x << std::endl;
    
    static double cached_x = 0.0;                       // (1)
    static double cached_result = 0.0;  // (2)
    
    std::cout << "cached_x:      " << cached_x << std::endl;
    std::cout << "cached_result: " << cached_result << std::endl;
    
    double result;
    
    if (cached_x == x)                                  // (1)
        return cached_result;                           // (2)
    std::cout << "compute from scratch" << endl;
    result = computation(x);
    cached_x = x;                                       // (1)
    cached_result = result;                             // (2)
    return result;
}

void doWork(string name) {
    for (int i = 0; i < 100; i++) {
        cout << name << ": " << i << endl;
        std::this_thread::sleep_for(1s);
    }
}

double threadTest() {
    std::thread worker_a(doWork, "A");
    std::thread worker_b(doWork, "B");
    
    worker_a.join();
    cout << "a finished" << endl;
    worker_b.join();
    cout << "b finished" << endl;
    return 0.0;
}

std::mutex m_read;
void fastaWorker(FastaReader& reader) {
    thread_local static int processed = 0;
    thread_local static long length = 0;
    const uintmax_t hundredth = reader.getFileSize()/100;
    static uintmax_t percent_counter = hundredth;
    
    //print process
    cout << "Process: " << (100*((double)reader.getProcessedBytes()/reader.getFileSize())) << "%            ";
    
    thread_local static FastaRecord record;
    while (true) {
        // print process
        if (reader.getProcessedBytes() > percent_counter) {
            cout << "\rProcess: " << (100*((double)reader.getProcessedBytes()/reader.getFileSize())) << "%            ";
            percent_counter += hundredth;
        }
        {
            std::lock_guard<std::mutex> lck(m_read);
            if (reader.hasNext()) {
                record = reader.next();
                length = length + record.sequence.length() + record.header.length();
                for (int i = 0; i < record.sequence.length(); i++) {
                    processed += record.sequence.c_str()[i];
                }
                
            } else {
                cout << endl;
                cout << std::this_thread::get_id() << " processed " << processed << " sum. " << endl;
                cout << std::this_thread::get_id() << " length " << length << " sum. " << endl;
                return;
            }
        }
    }
}


//static void next()


void fastaTest(int thread_count) {
    FastaReader reader("/home/joachim/CLionProjects/varkit/data/test/ultra_test.fa");
    
    std::vector<std::thread> threads;

    for (int t = 0; t < thread_count; t++) {
        string name = std::to_string(t);
        threads.emplace_back([&]() {fastaWorker(std::ref(reader));});
    }
    
    for (auto& t : threads)
        t.join();
    
    cout << "read_num: " << reader.getReadNum() << endl;
    cout << "processed_bytes: " << reader.getProcessedBytes() << endl;
}