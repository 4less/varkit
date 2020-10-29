//
// Created by joachim on 08/06/2020.
//

#ifndef VARKIT_KMERUTILS_H
#define VARKIT_KMERUTILS_H


#include <cstdint>
#include <string>

using namespace std;

namespace KmerUtils {
    static inline uint8_t getBitFromBase(char base) {
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
    
    static inline uint8_t getBitFromBase(const char* base) {
        switch((char)(*base)) {
            case 'A':
                return (0);
            case 'C':
                return (1);
            case 'G':
                return (2);
            case 'T':
                return (3);
        }
        return(-1);
    }
    
    static inline uint8_t getBitFromBaseC(char base) {
        switch(base) {
            case 'A':
                return (3);
            case 'C':
                return (2);
            case 'G':
                return (1);
            case 'T':
                return (0);
        }
        return(-1);
    }
    
    static inline uint8_t getBitFromBaseC(const char* base) {
        switch((char)(*base)) {
            case 'A':
                return (3);
            case 'C':
                return (2);
            case 'G':
                return (1);
            case 'T':
                return (0);
        }
        return(-1);
    }
    
    static inline char getBaseFromBit(uint8_t bits) {
        switch(bits) {
            case 0:
                return ('A');
            case 1:
                return ('C');
            case 2:
                return ('G');
            case 3:
                return ('T');
        }
        return('#');
    }
    
    static inline char complement(char c) {
        switch(c) {
            case 'A': return ('T');
            case 'C': return ('G');
            case 'G': return ('C');
            case 'T': return ('A');
        }
        return ('#');
    }
    
    static string bytesToDNA(uint8_t *kmer, int len) {
        string s;
        for (int i = 0; i < len; i++) {
            s += getBaseFromBit( (kmer[i/4] >> (2*(3-(i%4)))) & 0x03);
        }
        return s;
    }
    
    static string reverseComplement(string forward) {
        string reverse = "";
        const char * seq = forward.c_str();
        for (int i = forward.length()-1; i >= 0; i--) {
            reverse += complement(seq[i]);
        }
        return reverse;
    }
    
    static string dnaToBitString(string kmer) {
        const char* cseq = kmer.c_str();
        string int_string = "";
        for (int i = 0; i < kmer.length(); i++) {
            int_string += to_string(getBitFromBase(cseq[i]));
        }
        return int_string;
    }
    
    static void dnaToBytes(string kmer, uint8_t * key) {
        const char* cstr = kmer.c_str();
        
        for (int i = 0; i < kmer.length(); i++) {
            if (i % 4 == 0) key[i/4] = 0;
            key[i/4] = key[i/4] | KmerUtils::getBitFromBase(cstr[i]) << (2*(3-(i%4)));
        }
    }
}


#endif //KALAMITY_KMERUTILS_H
