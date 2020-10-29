//
// Created by joachim on 20/07/2020.
//

#pragma once

#include <vector>
#include <string>
#include <c++/7/bits/locale_facets.tcc>
#include <iomanip>

namespace Utils {
    static void shift(uint8_t * array, size_t length, size_t by) {
        for (int i = 0; i < length-1; i++) {
            array[i] <<= by;
            array[i] |= array[i + 1] >> (8-by);
        }
        array[length-1] <<= by;
    }
    
    static std::vector<std::string> split(const std::string& s, std::string delimiter)
    {
        std::vector<std::string> tokens;
        std::string token;
        
        std::string line = s;
        int lastpos = 0;
        int pos = 0;
        
        while ((pos = line.find(delimiter, lastpos)) != std::string::npos) {
            token = line.substr(lastpos, pos-lastpos);
            tokens.push_back(token);
            lastpos = pos+delimiter.length();
        }
        token = line.substr(lastpos, line.length()-lastpos);
        tokens.push_back(token);
        return tokens;
    }
    
    static std::vector<std::string> split(std::vector<std::string> &tokens, const std::string& s, std::string delimiter)
    {
        std::string token;
        
        std::string line = s;
        int lastpos = 0;
        int pos = 0;
        
        while ((pos = line.find(delimiter, lastpos)) != std::string::npos) {
            token = line.substr(lastpos, pos-lastpos);
            tokens.push_back(token);
            lastpos = pos+delimiter.length();
        }
        token = line.substr(lastpos, line.length()-lastpos);
        tokens.push_back(token);
        return tokens;
    }
    
    template <class T>
    static bool vectorContains(vector<T> v, T element) {
        return (v.find(element) == v.end());
    }
    
    static inline bool exists (const std::string& name) {
        if (FILE *file = fopen(name.c_str(), "r")) {
            fclose(file);
            return true;
        } else {
            return false;
        }
    }
    
    static std::string to_string_precision(double l, int precision) {
        std::ostringstream strs;
        strs << std::fixed << std::setprecision(precision)  << l;
        return strs.str();
    }
    
    template <class T, class U>
    static bool mapContains(unordered_map<T, U> v, T element) {
        return (v.find(element) == v.end());
    }
    
    static uint64_t swapEndianess(uint64_t &value, uint32_t bytes) {
        uint8_t swap[8];
        memset(swap, 0 ,8);
        for (int i = 7; i > (7-bytes); i--) {
            swap[7 - i] = ((uint8_t *) &value)[i];
        }
        return *((uint64_t*)swap);
    }
    
    static void swapEndianess(uint8_t* value, int32_t bytes) {
        uint8_t cache;
        for (int i = 7; i > (7-bytes) && (i >= 4); i--) {
            cache = value[i];
            value[i] = value[7-i];
            value[7-i] = cache;
        }
    }
    
    static std::string stripExtension(std::string file) {
        int index = file.rfind('.');
        int slash_i = file.rfind('/');
        if (index > 0 && slash_i > 0 && index > slash_i) {
            return file.substr(slash_i+1, index-slash_i-1);
        }
        return file;
    }
}