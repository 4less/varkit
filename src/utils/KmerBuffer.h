//
// Created by joachim on 23/09/2020.
//

#pragma once

#include <cstring>
#include "cstdint"
#include <iostream>
#include <iostream>

//using namespace std;

class KmerBuffer {
    uint8_t* buffer_;
    size_t buffer_size_ = 8;
    size_t fill_ = 0;
    size_t pos_ = -1;
    const size_t key_bytes_ = 0;
    const size_t value_bytes_ = 0;
    
    
    void resize(double by) {
        buffer_size_ *= by;
        buffer_ = (uint8_t*) realloc(buffer_, buffer_size_ * 16);
    }

    
public:
    KmerBuffer(size_t buffer_size, size_t key_bytes, size_t value_bytes) : buffer_size_(buffer_size), key_bytes_(key_bytes), value_bytes_(value_bytes) {
        cout << "malloc: " << buffer_size_ * 16 << endl;
        this->buffer_ = (uint8_t *) malloc(buffer_size_ * 16);
        memset(buffer_, 0, buffer_size_ * 16);
    }
    
    ~KmerBuffer() {
        free(buffer_);
    }
    
    void pushKV(uint8_t* key, uint8_t* value) {
        if (fill_ >= buffer_size_) resize(1.5);
        memcpy((buffer_ + 16 * fill_), key, key_bytes_);
        memcpy((buffer_ + 16 * fill_ + 8), value, value_bytes_);
        fill_++;
    }
    
    uint8_t* popKey() {
        return buffer_ + 16 * pos_;
    }
    
    uint8_t* popValue() {
        return buffer_ + 16 * pos_ + 8;
    }
    
    bool next() {
        pos_++;
        return pos_ < fill_;
    }
    
    
    void reset() {
        memset(buffer_, 0, buffer_size_ * 16);
        pos_ = -1;
        fill_ = 0;
    }
};
