//
// Created by fritsche on 12/04/2021.
//

#pragma once

#include <progress_bar.h>
#include "Progress.h"
#include <bits/unordered_set.h>
#include <set>
//#include <condition_variable>
#include "omp.h"
#include <iostream>
#include <cstring>
#include <numeric>
#include <fstream>
//#include <thread>
#include "time.h"

#include "robin_map.h"
#include "robin_set.h"
#include "sparse_map.h"
#include "sparse_set.h"

#include "CombinationIterator.h"
#include "Benchmark.h"
#include "KmerUtils.h"
#include "data_structures.h"
#include "csv.h"
#include "ShapeUtils.h"
#include "ShapeHistory.h"
//#include "BitUtils.h"
#include "Utils.h"
//#include <ThreadPool.h>
#include <thread_pool.h>

/**
 * Shapes are encoded as bools with
 * GAP = TRUE
 * TAKE = FALSE
 *
 * In a read a mutation is a TRUE
 * if GAP !XOR MUTATION we will get
 */


class ShapeResult;
class ShapeBuilder {
public:
    // External interface
    ShapeBuilder() {};
    ShapeBuilder(std::string shape_str);

    void Load(std::string path);
    void LoadArray(std::string path);
    void Save(std::string path);

    void TrainByPatternMT(size_t threads_num);
    void TrainByPatternTP(size_t threads_num);

    bool* StringToShape(std::string);
    std::string ShapeToString(bool* shape, size_t size);
    static std::string ReadToString(bool* read, size_t size);
    double Test(size_t read_length, size_t iterations, size_t mutations);

    std::string GetShapeStr();
    void PrintSpecs();

    ProgressInterface *progress = nullptr;

//    tsl::sparse_map<uint32_t,varkit::Mutations> map2;
    ds::Mutations* map_array;
    tsl::sparse_map<uint64_t,ds::Mutations> map2;

    size_t pattern_size = 22;
    size_t limit = 16;
    size_t lower_limit = 0;
    size_t maximum_chunk = 0;

    ~ShapeBuilder();

    static void ShapeSpaceExplorer(size_t take, size_t space, string output_folder, double random_prob, size_t num_threads, size_t init_pattern_size, size_t init_limit, size_t upper_limit, size_t iterations);
    static void DeepTrain(std::string output_folder, std::string output_stats, std::string output_pattern_distr, std::string shape, size_t num_threads, size_t p, size_t l_end);

    static double Sum(double* array, size_t size) {
        double result = 0;
        for (auto i = 0; i < size; i++) {
            result += array[i];
        }
        return result;
    }
    void TestTP(size_t read_length, size_t iterations, std::vector<size_t> mutations, size_t num_threads, ShapeResult &result);


    void TestTP(size_t read_length, size_t iterations, size_t *mutations, double *sensitivities, size_t mutations_size,
                size_t num_threads);

    size_t pot_mut_counter[100];

    static size_t median(size_t* ls, size_t size) {
        size_t median;
        std::vector<size_t> v;
        for (auto i = 0; i < size; i++) {
            if (ls[i])
                v.emplace_back(ls[i]);
        }
        std::sort (v.begin(), v.end());
        return v.at(v.size()/2);
    }

    static size_t GenerateReadRandom(bool *&read, size_t &read_length, double mutation_rate, double synonymous_prob);
    static size_t GenerateReadRandom(bool* &read, size_t &read_length, size_t mutation_num);
    void GenerateHitMissPattern(bool* read, size_t read_length, bool *pattern);

    static std::string RandomShape(size_t take, size_t space) {
        size_t shape_length = take + space;
        size_t partial_k = take / 2;
        size_t partial_space = space / 2;
        size_t range = partial_k - 1 +  partial_space;

//        srand(time(0));

        std::string shape(take + space , 'X');
        size_t s_counter = 0;
        while (s_counter < partial_space) {
            int random = (rand() % range) + 1;
            if (shape[random] != '_') {
                shape[random] = '_';
                shape[shape_length - 1 - random] = '_';
                s_counter++;
            }
        }
        return shape;
    }

    static std::string AlterShape(std::string shape, double rand_prob, size_t &take, size_t &space, std::string &operation);

    void Test(bool* read, size_t read_length, bool* hm_pattern, size_t hm_size, std::unordered_set<size_t> &mutations, size_t &tp, size_t &fp, size_t &fn, size_t* fn_pos, size_t &total_hits);

    static void Intersect(std::unordered_set<size_t> s1, std::unordered_set<size_t> s2, std::vector<size_t> &target) {
        if (s1.size() < s2.size()) {
            for (auto e : s1) {
                if (s2.contains(e)) target.emplace_back(e);
            }
        } else {
            for (auto e : s2) {
                if (s1.contains(e)) target.emplace_back(e);
            }
        }
    }

    size_t shape_size = 0;
    bool* shape = nullptr;

    typedef tsl::sparse_map<size_t, ds::Mutations> PatternMap;
    void Train(size_t threads);
    void Train(PatternMap& pattern_map, size_t from, size_t to);

private:

    std::mutex pot_mut_mutex;

    inline static char take = 'X';
    inline static char space = '_';

    //debug chars
    size_t count_total = 0;
    size_t count_below = 0;
    size_t count_no_valid_patterns = 0;

    // Internal methods

    void TestPattern(size_t mutations, size_t iterations, size_t &tp, size_t &fp, size_t &fn);

    void Update(size_t mutation_count, size_t iterations);
    void Update(size_t mutation_count, size_t read_length, size_t iterations);
    void GenerateRead(bool* &read, size_t &read_length, size_t* distances, size_t mutation_count);
    size_t GetMutations(bool* read, size_t read_length, std::unordered_set<size_t> &mutations);
    size_t GetMutations(bool* read, size_t read_length, size_t* mutations, size_t mutations_size);
    size_t GenerateReadRandom(bool* &read, size_t &read_length, size_t mutation_num, size_t* distances, size_t max_dist);
    ds::Mutations TrainByPatternOMP(uint32_t pattern, bool* read, bool* potential_mutations, bool *hm_pattern, size_t *mutations, size_t &total_count, size_t &valid_count, size_t* mut_counter);

//    void TrainByPatternMT(uint32_t* p, size_t size, varkit::Mutations* tmp);//, size_t *tmp);
    void TrainByPatternMT(uint64_t* p, size_t size, ds::Mutations* tmp);//, size_t *tmp);

    // Utility
    void SetBit(size_t pos, bool to, uint32_t& pattern);
    void SetBit(size_t pos, bool to, uint64_t& pattern);

    static void Intersect(tsl::robin_set<uint8_t> &s1, tsl::robin_set<uint8_t> &s2, tsl::robin_set<uint8_t> &target) {
        for (auto e : s1) {
            if (s2.contains(e))
                target.insert(e);
        }
    }

    static bool Contains(size_t &s, uint8_t &e) {
        uint8_t* sa = (uint8_t*) &s;
        for (auto i = 0; i < 8 && sa[i] > 0; i++) {
            if (sa[i] == e) return true;
        }
        return false;
    }

    static bool Contains(size_t *s, size_t s_size, size_t &e) {
        for (auto i = 0; i < s_size; i++) {
            if (s[i] == e) return true;
        }
        return false;
    }

    static void Intersect(size_t &s1, size_t &s2, size_t &s3) {
        s3 = 0;
        uint8_t* s1a = (uint8_t*) &s1;
        uint8_t* s2a = (uint8_t*) &s2;
        uint8_t* s3a = (uint8_t*) &s3;
        size_t s3i = 0;
        for (auto i = 0; i < 8 && s1a[i] > 0; i++) {
            for (auto j = 0; j < 8 && s2a[j] > 0; j++) {
                if (s1a[i] == s2a[j] && !Contains(s3, s1a[i]))
                    s3a[s3i++] = s1a[i];
            }
        }
        if (s3a[0] && s3a[0] == s3a[1]) {
            std::cout << ToString(s3) << std::endl;
            exit(0);
        }
    }

    static void Intersect(size_t &s1, size_t* s2, size_t s2_size, size_t &s3) {
        s3 = 0;
        uint8_t* s1a = (uint8_t*) &s1;
        uint8_t* s3a = (uint8_t*) &s3;
        size_t s3i = 0;
        for (auto i = 0; i < 8 && s1a[i] > 0; i++) {
            for (auto j = 0; j < s2_size; j++) {
                if (s1a[i] == s2[j] && !Contains(s3, s1a[i]))
                    s3a[s3i++] = s1a[i];
            }
        }
        if (s3a[0] && s3a[0] == s3a[1]) {
            std::cout << ToString(s3) << std::endl;
            exit(0);
        }
    }

    static size_t Intersect(size_t* s1, size_t s1_size, size_t* s2, size_t s2_size, size_t* s3, size_t s3_size) {
        memset(s3, 0, s3_size * sizeof(size_t));

        size_t s3i = 0;
        for (auto i = 0; i < s1_size; i++) {
            for (auto j = 0; j < s2_size; j++) {
                if (s1[i] == s2[j] && !Contains(s3, s3_size, s1[i])) {
                    s3[s3i++] = s1[i];
                }
            }
        }
        return s3i;
    }

    static std::string ToBitString(uint32_t num) {
        std::string bits(32, '0');
        for (auto i = 0; i < 32; i++) {
//            std::cout << "bit " << (num & (1 << (31 - i))) << std::endl;
//            std::cout << "bit " << num << " " << (1 << (31 - i)) << std::endl;
            if (num & (1 << (31 - i)))
                bits[i] = '1';
        }
        return bits;
    }

    static std::string ToBitString(uint64_t num) {
        std::string bits(64, '0');
        for (uint64_t i = 0; i < 64; i++) {
            if (num & (1llu << (64 - i))) {
                bits[i] = '1';
            }
        }
        return bits;
    }

    static std::string ToBitString(uint32_t num, size_t size) {
        std::string bits(size, '0');
        for (auto i = 0; i < size; i++) {
            if (num & (1 << (size-1 - i)))
                bits[i] = '1';
        }
        return bits;
    }

    static size_t ToULL(std::string shape_str) {
        uint32_t shape = 0;
        for (size_t pos = 0; pos < shape_str.size(); pos++) {
            if (shape_str[pos] == '_') continue;
            shape |= 1 << (shape_str.size() - 1 - pos);
        }
        return shape;
    }

    static double RandomProb() {
//        srand(time(NULL));
        size_t precision = 1000;
        return ((double) (rand() % precision) / precision);
    }

    static bool BiasedCoinFlip(double true_rate) {
        return RandomProb() < true_rate;
    }

    static void MoveSpace(std::string &shape, int dir, size_t pos, size_t take, size_t space, std::string &operation);

    static void Add(size_t &s, size_t e) {
        uint8_t* sa = (uint8_t*) &s;
        for (auto i = 0; i < 8; i++) {
            if (sa[i] == e) return;
            if (sa[i] == 0xFF) {
                sa[i] = e;
                return;
            }
        }
    }

    // This size_t is treated as uint8_t[8] and printed as an array
    static std::string ToString(size_t &s) {
        uint8_t* sa = (uint8_t*) &s;
        std::string res = "[";
        for (auto si = 0;; si++) {
            if (si >= 8) {// || sa[si] == 0xFF) {
                res += "]";
                break;
            } else {
                if (si != 0) res += ", ";
                res += std::to_string((uint32_t) sa[si]);
            }
        }
        return res;
    }

    static std::string ToString(std::unordered_set<size_t> &s) {
        std::string res = "{";
        bool isfirst = true;
        for (auto e : s) {
            if (isfirst) {
                isfirst = false;
            } else {
                res += ", ";
            }
            res += std::to_string(e);
        }

        res += "}";
        return res;
    }

    static std::string ToString(std::vector<size_t> &s) {
        std::string res = "[";
        bool isfirst = true;
        for (auto e : s) {
            if (isfirst) {
                isfirst = false;
            } else {
                res += ", ";
            }
            res += std::to_string(e);
        }

        res += "]";
        return res;
    }

    static size_t ChooseN(size_t k, size_t n_start, size_t threshold);

    static size_t FirstFromLeft(uint32_t pattern, size_t size, bool shift) {
        if (shift)
            pattern >>= 32 - size;
        size_t ppos = (bool) pattern * 1;
        while ((pattern /= 2))
            ppos++;

        return size - ppos;
    }

    static size_t InformationContent(uint32_t pattern, size_t size) {
        uint32_t mask = 1 << 31;
        size_t res = 0;
        bool last = (pattern & mask);
        for (auto pos = 1; pos < size; pos++) {
            mask = 1 << (31 - pos);
            bool current = pattern & mask;
            if ( current != last)
                res++;
            last = current;
        }
        return res;
    }

    static unsigned int CountGapNum(std::string shape) {
        uint32_t gap_num = 0;
        char old_char = 'X';
        for (auto i = 0; i < shape.size(); i++) {
            if (shape[i] == 'X' && shape[i] != old_char)
                gap_num++;
            if (gap_num > 100) {
                std::cout << "major suspicious : "<< i << std::endl;
                exit(8);
            }
            old_char = shape[i];
        }
        return gap_num;
    }

    static unsigned int CountSetBits(unsigned int n)
    {
        unsigned int count = 0;
        while (n) {
            count += n & 1;
            n >>= 1;
        }
        return count;
    }

    static void dummyTask() {
        std::cout << "dummy" << std::endl;
        size_t dummy;
        for (size_t t = 0; t < 100000000; t++)
            dummy += t;
        std::cout << dummy << std::endl;
    }

    static void AwithoutB(std::unordered_set<size_t> a, std::unordered_set<size_t> b, std::vector<size_t> &target) {
        for (auto e : a) {
            if (!b.contains(e)) target.emplace_back(e);
        }
    }

    static size_t DifferentSpacesLengths(bool* shape, size_t shape_size) {
        std::unordered_set<size_t> lengths;
        size_t spaces = 0;
        for (auto i = 0; i < shape_size; i++) {
            if (shape[i]) spaces++;
            if (!shape[i] && spaces != 0) {
                lengths.insert(spaces);
                spaces = 0;
            }
        }
        return lengths.size();
    }

    static size_t LongestConsecutive(bool* shape, size_t shape_size) {
        size_t max = 0;
        size_t consecutive = 0;
        for (auto i = 0; i < shape_size; i++) {
            if (!shape[i]) consecutive++;
            if (shape[i] && consecutive != 0) {
                if (consecutive > max) max = consecutive;
                consecutive = 0;
            }
        }
        return max;
    }


    // Tests
    void TrainPatternBased();
    void TrainByPattern();

    static size_t RandomBiasedOption(double probabilities[3], size_t size);

    void TestTP(size_t read_length, size_t iterations, size_t mutations, double *sensitivity);
    void TestTP(size_t read_length, size_t iterations, double ani, ShapeResult &result);

    static bool RandomCoinFlip(double true_prob);

    void ShapeSpaceExplorerSynonymous(size_t take, size_t space, string outfile, double random_prob, size_t num_threads,
                                      size_t init_pattern_size, size_t init_limit, size_t upper_limit,
                                      size_t iterations);


    void Print();

};