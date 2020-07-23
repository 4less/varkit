//
// Created by joachim on 20/07/2020.
//

#pragma once

#include <vector>
#include <string>
#include <c++/7/bits/locale_facets.tcc>

namespace Utils {
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
    
    static uint64_t swapEndianess(uint64_t &value, uint32_t bytes) {
        uint8_t swap[8];
        memset(swap, 0 ,8);
        for (int i = 7; i > (7-bytes); i--) {
            swap[7 - i] = ((uint8_t *) &value)[i];
        }
        return *((uint64_t*)swap);
    }
    
    static void swapEndianess(uint8_t* value, uint32_t bytes) {
        uint8_t cache;
        for (int i = 7; i > (7-bytes); i--) {
            cache = value[i];
            value[i] = value[7-i];
            value[7-i] = cache;
        }
    }
}