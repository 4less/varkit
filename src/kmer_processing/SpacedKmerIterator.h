//
// Created by joachim on 13/06/2020.
//

#ifndef VARKIT_SPACEDKMERITERATOR_H
#define VARKIT_SPACEDKMERITERATOR_H

#include <FastxReader.h>
#include "KmerIterator.h"

class SpacedKmerIterator : public KmerIterator {
    uint32_t key_bytes_;
    uint8_t * key_;
    uint8_t * key_rc_;
    uint32_t k_;
    uint32_t pos_= 0;
    int buffer_;
    bool * shape_;
    uint32_t shape_size_;
    
    FastaRecord record_;
    FastxRecord recordx_;
    size_t sequence_length_ = 0;
    const char * sequence_ = nullptr;
    
    void getKey(uint8_t * key) {
        memset(key, 0, key_bytes_);
        int offset = 0;
        for (int i = 0; i < k_; i++) {
            while (shape_[i+offset]) {
                offset++;
            }
            key[i/4] = key[i/4] | KmerUtils::getBitFromBase(sequence_[pos_ + i + offset]) << (2*(3-(i%4)));
            
            if (KmerUtils::getBitFromBase(sequence_[pos_ + i + offset]) == -1) {
                memset(key, 1, key_bytes_);
                return;
            }
        }
    }
    
    void getKeyRC(uint8_t * key) {
        memset(key, 0, key_bytes_);
        
        int offset = 0;
        for (int i = 0; i < k_; i++) {
            while (shape_[shape_size_-1-(i+offset)]) {
                offset++;
            }
            key[i/4] = key[i/4] | KmerUtils::getBitFromBaseC(sequence_[pos_ + (shape_size_-1) - (i + offset)]) << (2*(3-(i%4)));
        }
    }
    
    void inline getKeyBoth(uint8_t * key, uint8_t * key_rc) {
        memset(key, 0, key_bytes_);
        memset(key_rc, 0, key_bytes_);
        
        int offset = 0;
        for (int i = 0; i < k_; i++) {
            while (shape_[i+offset])
                offset++;
            key[i/4] |= KmerUtils::getBitFromBase(sequence_[pos_ + i + offset]) << (2*(3-(i%4)));
            key_rc[i/4] |= KmerUtils::getBitFromBaseC(sequence_[pos_ + (shape_size_-1) - (i + offset)]) << (2*(3-(i%4)));
        }
    }
    
public:
    void setFastaRecord(FastaRecord record) {
        record_ = record;
        /*
        if (buffer_ == 0 && sequence_) {
            //cout << "free" << endl;
            free(sequence_);
        }
        
        if (buffer_ == 0)
            sequence_ = (char*)malloc(record.sequence.length() * sizeof(char));
        
        memcpy(sequence_, record.sequence.c_str(), record.sequence.length() * sizeof(char));
        */
        pos_ = 0;
    }
    
    void setRecord(FastxRecord &record) {
        sequence_ = record.sequence.c_str();
        sequence_length_ = record.sequence.size();
        recordx_ = record;
        pos_ = 0;
    }
    
    string getSubstring(int pos, int length) {
        int offset = 0;
        if (pos < 0)
            offset = pos;
        else if (pos+length >= sequence_length_) {
            offset = pos + length - sequence_length_;
        }
        return recordx_.sequence.substr(pos - offset, length);
    }
    
    uint32_t currentLength() {
        return sequence_length_;
    }
    
    inline bool hasNext() {
        return (pos_ < sequence_length_ - shape_size_+1);
    }
    
    bool inline isSmallerThan(uint8_t * a, uint8_t * b) {
        for (int i = 0; i < key_bytes_; i++) {
            if (a[i] != b[i])
                return (a[i] < b[i]);
        }
        //cout << "equal" << endl;
        return true;
    }
    
    explicit SpacedKmerIterator(uint32_t k, bool * shape, uint32_t shape_size, int buffer = 0) : k_(k), shape_(shape), shape_size_(shape_size), buffer_(buffer) {
        assert(k % 2 == 0);
        
        key_bytes_ = (2*k + 8 - 1) /  8;
        
        key_ = new uint8_t[key_bytes_];
        key_rc_ = new uint8_t[key_bytes_];

        if (buffer > 0) {
            sequence_ = (char*)malloc(buffer * sizeof(char));
        }
    }
    
    ~SpacedKmerIterator() {
        delete[] key_;
        delete[] key_rc_;
    }
    
    void inline operator () (uint8_t * key) {
        //getKey(key_);
        //getKeyRC(key_rc_);
        
        getKeyBoth(key_, key_rc_);
        
        if (isSmallerThan(key_, key_rc_))
            memcpy(key, key_, key_bytes_);
        else
            memcpy(key, key_rc_, key_bytes_);
        pos_++;
    }
    
    uint32_t getPos() {
        return pos_-1;
    }
};


#endif //KALAMITY_SPACEDKMERITERATOR_H
