//
// Created by joachim on 27/05/2020.
//

#ifndef KALAMITY_FASTAREADER_H
#define KALAMITY_FASTAREADER_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
using namespace std;

struct FastaRecord {
    string header;
    string sequence;

    FastaRecord(string const& header, string const& sequence) {
        this->header = header;
        this->sequence = sequence;
    }

    FastaRecord() = default;

    void print() {
        cout << header << endl;
        cout << sequence << endl;
    }
};

class FastaReader {
private:
    string file_;
    uintmax_t file_size_;
    uintmax_t processed_bytes_ = 0;
    string line_;
    string header_;
    string sequence_;
    ifstream inFile;
    uint32_t read_num_ = 0;
    
    int buffer_size_;
    char* sequence_buffer_;
    long processed = 0;
    
public:
    FastaReader(string file, int buffer_size=1<<12);
    ~FastaReader();
    
    inline FastaRecord next() {
        header_ = line_;
        sequence_ = "";
        
        while (getline(inFile, line_) && line_[0] != '>') {
            processed_bytes_ = processed_bytes_ + line_.size() + 1;
            sequence_ += line_;
        }
        processed_bytes_ = processed_bytes_ + line_.size() + 1;
        read_num_++;
        
        return FastaRecord(header_, sequence_);
    }
    
    inline bool hasNext() {
        if (inFile.eof()) {
            processed_bytes_-1;
            inFile.close();
            return false;
        }
        return true;
    }
    
    int processedCount();
    uint32_t getReadNum();
    uintmax_t getProcessedBytes();
    string getFile();
    
    uintmax_t getFileSize();
};





#endif //KALAMITY_FASTAREADER_H
