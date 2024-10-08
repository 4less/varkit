//
// Created by fritsche on 12/04/2021.
//

#pragma once


#include <cstddef>
#include <numeric>
#include <iostream>
#include <cstring>

class CombinationIterator {
    size_t n = 0; // size of areas
    size_t k = 0; // number of distances

    size_t last_sum = 0;
    bool init = true;

public:
    CombinationIterator(size_t n, size_t k) : n(n), k(k) {};

    void Init(size_t* distances) {
        memset(distances, 0, k * sizeof(size_t));
        init = true;
    }

    bool Next(size_t* distances) {
        if (init) {
            init = false;
            return true;
        }
//        for (int i = 0; i < k; i++) {
//            std::cout << i << ": " << distances[i] << std::endl;
//        }

        size_t index = k-1;
        distances[index]++;
        last_sum++;
//        std::cout << "next sum: " << sum << " index: " << index << std::endl;
//        std::cout << "next k: " << k << " n: " << n << std::endl;
        while ((last_sum + k) >= n) {
            if (index == 0)
                return false;
            distances[index--] = 0;
            distances[index]++;
            last_sum = std::accumulate(distances, distances+k, 0);
        }
        return true;
    }
};