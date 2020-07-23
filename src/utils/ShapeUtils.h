//
// Created by joachim on 01/07/2020.
//

#ifndef VARKIT_SHAPEUTILS_H
#define VARKIT_SHAPEUTILS_H

#include <iostream>

using namespace std;

namespace shape_utils {
    static void getBool(string shape, bool *target) {
        const char *shape_array = shape.c_str();
        
        for (int i = 0; i < shape.length(); i++) {
            target[i] = shape_array[i] == '_';
        }
    }
    
    static string getString(bool shape[], int len) {
        string result = "";
        for (int i = 0; i < len; i++) {
            if (shape[i]) result += "_";
            else result += "X";
        }
        return result;
    }
};

#endif //VARKIT_SHAPEUTILS_H
