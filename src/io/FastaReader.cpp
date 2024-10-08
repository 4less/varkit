//
// Created by joachim on 27/05/2020.
//

#include <cstring>
#include "FastaReader.h"
#include <experimental/filesystem>


FastaReader::FastaReader(string file, int buffer_size) : file_ (std::move(file)), buffer_size_(buffer_size), file_size_(std::experimental::filesystem::file_size(file_)) {
    sequence_buffer_ = new char[buffer_size_];
    inFile.rdbuf()->pubsetbuf(sequence_buffer_, buffer_size_);
    
    cout << "filesize: " << this->file_size_ << endl;
    
    inFile.open(this->file_);

    if (!inFile) {
        cerr << "FastaReader::Unable to open file_.\n" << this->file_;
        exit(1);
    }

    // read first line_
    getline(inFile, line_);

    if (line_[0] != '>') {
        cerr << "First line in fasta file_ has to start with '>' (header)";
        exit(1);
    }
    processed_bytes_ += line_.size()+1;
}

FastaReader::~FastaReader() {
    delete[] sequence_buffer_;
}


inline int FastaReader::processedCount() {
    return 0;
}

string FastaReader::getFile() {
    return file_;
}



uint32_t FastaReader::getReadNum() {
    return read_num_;
}

uintmax_t FastaReader::getProcessedBytes() {
    return processed_bytes_;
}

uintmax_t FastaReader::getFileSize() {
    return std::experimental::filesystem::file_size(file_);
}



