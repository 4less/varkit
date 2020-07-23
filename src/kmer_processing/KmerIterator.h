//
// Created by joachim on 08/06/2020.
//

#ifndef VARKIT_KMERITERATOR_H
#define VARKIT_KMERITERATOR_H


#include "FastaReader.h"
#include "KmerUtils.h"
//#include "KmerProcessor.h"
#include <assert.h>
#include <cstring>

class KmerIterator {
public:
    virtual void operator() (uint8_t * key) = 0;
};

class FastKmerIterator : KmerIterator {
private:
    FastaRecord record_;
    uint32_t key_bytes_;
    uint32_t last_pos_shift_;
    uint32_t k_;
    uint32_t pos_ = 1;
    char * sequence_ = nullptr;
    bool first_ = true;
    
    void getKey(uint8_t * key) {
        for (int i = 0; i < k_; i++) {
            if (i % 4 == 0) key[i/4] = 0;
            key[i/4] = key[i/4] | KmerUtils::getBitFromBase(sequence_[i]) << (2*(3-(i%4)));
        }
    }
    
    void rollKey(uint8_t * key) {
        for (int i = 0; i < key_bytes_-1; i++) {
            key[i] = (key[i] << 2) | ((key[i+1]) >> 6);
        }
        key[key_bytes_-1] = (key[key_bytes_-1] << 2) | (KmerUtils::getBitFromBase((sequence_ + pos_)[k_-1]) << last_pos_shift_);
    }

public:
    void setFastaRecord(FastaRecord record) {
        record_ = record;
        if (!sequence_) free(sequence_);
        sequence_ = (char*)malloc(record.sequence.length() * sizeof(char));
        memcpy(sequence_, record.sequence.c_str(), record.sequence.length() * sizeof(char));
        pos_ = 1;
        first_ = true;
    }
    
    bool hasNext() {
        return (pos_ < record_.sequence.length()-k_+1);
    }
    
    explicit FastKmerIterator(uint32_t k) : k_(k) {
        key_bytes_ = (2*k + 8 - 1) /  8;
        last_pos_shift_ = (8 - ((2*k) % 8)) % 8;
    }
    
    void operator () (uint8_t * key) {
        if (first_) {
            getKey(key);
            first_ = false;
        } else {
            rollKey(key);
            pos_++;
        }
    }
};

class CanonicalKmerIterator : KmerIterator {
private:
    FastaRecord record_;
    uint32_t key_bytes_;
    uint32_t last_pos_shift_;
    uint32_t k_;
    uint32_t pos_ = 1;
    char * sequence_ = nullptr;
    bool first_ = true;
    
    uint8_t * key_;
    uint8_t * key_can_;
    
    void getKey(uint8_t * key) {
        for (int i = 0; i < k_; i++) {
            if (i % 4 == 0) key[i/4] = 0;
            key[i/4] = key[i/4] | KmerUtils::getBitFromBase(sequence_[k_]) << (2*(3-(i%4)));
        }
    }
    
    void getKeyC(uint8_t * key) {
        for (int i = 0; i < k_; i++) {
            if (i % 4 == 0) key[i/4] = 0;
            key[i/4] = key[i/4] | KmerUtils::getBitFromBaseC(sequence_[k_ - 1 - i]) << (2*(3-(i%4)));
        }
    }
    
    void rollKey(uint8_t * key) {
        for (int i = 0; i < key_bytes_-1; i++) {
            key[i] = (key[i] << 2) | ((key[i+1]) >> 6);
        }
        key[key_bytes_-1] = (key[key_bytes_-1] << 2) | (KmerUtils::getBitFromBase((sequence_ + pos_)[k_-1]) << last_pos_shift_);
    }
    
    void rollKeyC(uint8_t * key) {
        for (int i = 1; i < key_bytes_; i++) {
            key[i] = (key[i] >> 2) | ((key[i-1]) << 6);
        }
        key[0] = (key[0] >> 2) | (KmerUtils::getBitFromBaseC((sequence_ + pos_)[k_-1]) << 6);
        key[key_bytes_-1] &= ((1 << last_pos_shift_)-1);
    }

public:
    void setFastaRecord(FastaRecord record) {
        record_ = record;
        if (!sequence_) free(sequence_);
        sequence_ = (char*)malloc(record.sequence.length() * sizeof(char));
        memcpy(sequence_, record.sequence.c_str(), record.sequence.length() * sizeof(char));
        pos_ = 1;
        first_ = true;
    }

    bool hasNext() {
        return (pos_ < record_.sequence.length()-k_+1);
    }
    
    explicit CanonicalKmerIterator(uint32_t k) : k_(k) {
        key_bytes_ = (2*k + 8 - 1) /  8;
        last_pos_shift_ = (8 - ((2*k) % 8)) % 8;
    
        key_ = new uint8_t[key_bytes_];
        key_can_ = new uint8_t[key_bytes_];
    }
    ~CanonicalKmerIterator() {
        delete[] key_;
        delete[] key_can_;
    }
    
    bool isSmallerThan(uint8_t * a, uint8_t * b) {
        for (int i = 0; i < key_bytes_; i++) {
            if (a[i] != b[i])
                return (a[i] < b[i]);
        }
        cout << "equal" << endl;
        return true;
    }
    
    void operator () (uint8_t * key) {
        if (first_) {
            getKey(key_);
            getKeyC(key_can_);
            first_ = false;
        } else {
            rollKey(key_);
            rollKeyC(key_can_);
        }
        
        //cout << "normal, canonical " << BitUtils::ToBitString(key_, key_bytes_) << "   " << BitUtils::ToBitString(key_can_, key_bytes_) << endl;
        
        if (isSmallerThan(key_, key_can_)) {
            //cout << "norm" << endl;
            memcpy(key, key_, key_bytes_);
        }
        else {
            //cout << "can" << endl;
            memcpy(key, key_can_, key_bytes_);
        }
        pos_++;
        //if (pos_ > 10) exit(0);
    }
};


class CanonicalSingleSpacedKmerIterator : KmerIterator {
private:
    FastaRecord record_;
    uint32_t key_bytes_;
    uint32_t k_;
    uint32_t pos_ = 0;
    char * sequence_ = nullptr;
    bool first_ = true;
    int skip_;
    int buffer_;
    uint8_t * key_;
    uint8_t * key_can_;
    
    void getKey(uint8_t * key) {
        //cout << "getKey:  " << endl;
        memset(key, 0, key_bytes_);
        int offset = 0;
        for (int i = 0; i < k_; i++) {
            if (i == skip_) {
                //cout << ".";
                offset++;
            }
            //cout << sequence_[pos_+i+offset];
            //cout << KmerUtils::getBaseFromBit(KmerUtils::getBitFromBase(sequence_[pos_ + offset + i])) ;
            key[i/4] = key[i/4] | KmerUtils::getBitFromBase(sequence_[pos_ + i + offset]) << (2*(3-(i%4)));
        }
        //cout << endl << KmerUtils::bytesToDNA(key, k_) << endl;
    }
    
    void getKeyC(uint8_t * key) {
        
        //TODO: print out skip part to see were the middle position actually is ! :)
        //cout << "getKeyRC: " << endl;
        memset(key, 0, key_bytes_);
        
        int offset = 0;
        for (int i = 0; i < k_; i++) {
            if (i == skip_) {
                //cout << ".";
                offset++;
            }
            //cout << KmerUtils::getBaseFromBit(KmerUtils::getBitFromBaseC(sequence_[pos_ + k_ - offset - i]));
            key[i/4] = key[i/4] | KmerUtils::getBitFromBaseC(sequence_[pos_ + k_ - i - offset]) << (2*(3-(i%4)));
        }
        //cout << endl << KmerUtils::bytesToDNA(key, k_) << endl;
    }

public:
    void setFastaRecord(FastaRecord record) {
        record_ = record;
        
        if (buffer_ == 0 && sequence_) {
            cout << "free" << endl;
            free(sequence_);
        }

        if (buffer_ == 0)
            sequence_ = (char*)malloc(record.sequence.length() * sizeof(char));
        
        memcpy(sequence_, record.sequence.c_str(), record.sequence.length() * sizeof(char));
        
        pos_ = 0;
        first_ = true;
        
        //for (int i = 0; i < 30; i++)
        //    cout << sequence_[i];
        //cout << endl << record.sequence.substr(0, 30) << endl;
        
    }
    
    bool hasNext() {
        return (pos_ < record_.sequence.length()-k_);
    }

    
    explicit CanonicalSingleSpacedKmerIterator(uint32_t k, int buffer = 0) : k_(k), buffer_(buffer) {
        assert(k % 2 == 0);
        
        key_bytes_ = (2*k + 8 - 1) /  8;
        //last_pos_shift_ = (8 - ((2*k) % 8)) % 8;
        
        key_ = new uint8_t[key_bytes_];
        key_can_ = new uint8_t[key_bytes_];
        
        skip_ = k/2;
        
        if (buffer > 0) {
            sequence_ = (char*)malloc(buffer * sizeof(char));
        }
    }
    
    ~CanonicalSingleSpacedKmerIterator() {
        delete[] key_;
        delete[] key_can_;
        free(sequence_);
    }
    
    bool isSmallerThan(uint8_t * a, uint8_t * b) {
        for (int i = 0; i < key_bytes_; i++) {
            if (a[i] != b[i])
                return (a[i] < b[i]);
        }
        //cout << "equal" << endl;
        return true;
    }
    
    string getLastKmer() {
        //cout << "pos_" << pos_ << " k_" << k_ << "  begin: " << (pos_-2) << " end: " <<  (pos_ - 2 + k_ + 1) << ", length: " << record_.sequence.length() << endl;
        return record_.sequence.substr(pos_-2, k_+1);
    }
    
    string getLastKmer(int offset) {
        return record_.sequence.substr(pos_+offset-1, k_+1);
    }
    
    string getLastKmer2() {
        string kmer;
        for (int i = pos_; i < pos_+k_+1; i++) {
            kmer += sequence_[i];
        }
        return kmer;
    }
    
    void operator () (uint8_t * key) {
        getKey(key_);
        getKeyC(key_can_);
        
        
        
        if (isSmallerThan(key_, key_can_)) {
            //cout << "norm" << endl;
            memcpy(key, key_, key_bytes_);
        }
        else {
            //cout << "can" << endl;
            memcpy(key, key_can_, key_bytes_);
        }
        pos_++;
        //if (pos_ > 10) exit(0);
    }
};

#endif //KALAMITY_KMERITERATOR_H
