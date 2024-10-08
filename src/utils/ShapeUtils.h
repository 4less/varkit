//
// Created by joachim on 01/07/2020.
//

#pragma once

#include <fstream>
#include <regex>
#include <sysexits.h>
#include "err.h"
#include "assert.h"

namespace ShapeUtils {
    static void getBool(std::string shape, bool *target) {
        const char *shape_array = shape.c_str();
        
        for (int i = 0; i < shape.length(); i++) {
            target[i] = shape_array[i] == '_';
        }
    }
    
    static std::string GetString(const bool *shape, int len) {
        std::string result = "";
        for (int i = 0; i < len; i++) {
            if (shape[i]) result += "_";
            else result += "X";
        }
        return result;
    }

    static std::string LoadShape(std::string path) {
        std::ifstream shapein(path, std::ios::in);
        std::string line;
        getline(shapein, line);
        if (!std::regex_match (line, std::regex("[X|_]+") )) {
            errx(EX_IOERR,
                 "Shape file may only contain the characters X and _. Instead the shape string looks like %s",
                 line.c_str());
        }
        shapein.close();
        return line;
    }

    static size_t LongestGap(const bool* shape, size_t shape_size) {
        size_t longest_gap = 0;
        size_t gap_size = 0;
        for (auto i = 0; i < shape_size; i++) {
            gap_size = shape[i] * (gap_size + 1);
            longest_gap = gap_size > longest_gap ? gap_size : longest_gap;
        }
        return longest_gap;
    }

    static size_t GetK(std::string shape) {
        int count = 0;
        for (int i = 0; i < shape.size(); i++)
            if (shape[i] == 'X') count++;

        return count;
    }

    static std::string ApplyShape(std::string seq, size_t pos, const bool* shape, size_t shape_len, bool include_gaps=false) {
        assert(seq.size() >= pos + shape_len);

        std::string result = "";

        for (int i = 0; i < shape_len; i++) {
            auto s_pos = i + pos;
            if (!shape[i])
                result += seq[s_pos];
            else if (include_gaps) {
                result += '_';
            }
        }
        return result;
    }

    static bool* GetShape(std::string s) {
        bool * shape = new bool[s.length()];
        const char * c = s.c_str();
        for (int i = 0; i < s.length(); i++) {
            shape[i] = c[i] == '_';
        }
        return shape;
    }
};
