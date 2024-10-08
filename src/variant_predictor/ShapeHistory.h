//
// Created by fritsche on 11/08/2021.
//

#ifndef VARKIT_SHAPEHISTORY_H
#define VARKIT_SHAPEHISTORY_H


#include <string>
#include <unordered_set>
#include <vector>
#include <iostream>
#include "csv.h"


#include "ShapeBuilder.h"

class ShapeBuilder;

class StatsContainer {
public:
    size_t tp, fp, tn, fn;

    double Sensitivity() {
        // Also recall, hitrate or true positive rate
        return (double) tp / (tp + fn);
    }

    double Precision() {
        // Also positive predictive value
        return (double) tp / (tp + fp);
    }

    double Specificity() {
        return (double) tn / (tn + fp);
    }
};

class ShapeResult;
//typedef int& (* Comparator)(ShapeResult &a, ShapeResult &b);
//typedef int (* Comparator)(::ShapeResult &,::ShapeResult &);
typedef std::function<int(ShapeResult&, ShapeResult&)> Comparator;

class ShapeResult {
public:
    std::string shape_;

    std::unordered_map<std::string, StatsContainer> statistics;
    std::unordered_map<std::string, size_t> hits;
    std::unordered_map<std::string, size_t> kmers;
    std::unordered_map<std::string, double> sensitivities;
    std::vector<size_t> positional;

    size_t total_hits = 0;
    size_t total_kmers = 0;
    size_t total_reads = 0;

    ShapeResult(std::string shape);

//    static std::vector<ShapeResult> LoadResults(std::string path);
//    static ShapeResult LineToResult(std::string line);
    int Compare(ShapeResult &other, Comparator function);
    StatsContainer& GetStatistic(std::string key);
    void AddPosition(size_t pos, size_t counts);
    void Print(ostream &os);
    void CalculateStats();

    void PrintMean(ostream &os, vector<double> anis);
};


class ShapeHistory {
    std::vector<ShapeResult> shapes;
    std::unordered_set<std::string> tested_shapes;

    ShapeHistory(std::string path);


    static double StringToAni(std::string ani_str);

public:
    std::vector<std::string> sens_colnames;
    std::vector<std::string> hit_colnames;
    std::vector<double> anis;

    static inline std::string prefix = "ANI_";

    ShapeHistory();

    bool Empty();


    ShapeResult& GetBest(Comparator fun);
//    ShapeResult& GetBest();

    void Add(ShapeResult &res);

    bool HasShape(std::string shape_str);

    void LoadHistory(std::string history);

    static std::string AniToString(double ani);

    void Write(std::string file);
};

class FinderOptions {
public:
    size_t take;
    size_t space;
    std::string output;
    int threads;
    size_t pattern_size;
    size_t training_max;
    size_t training_init;
    size_t iterations;
    std::string checked_options;
    std::string previous_findings = "";

};

class  ShapeFinder {
    FinderOptions options_;
    ShapeHistory history_;
    ProgressSingleton *progress_ = ProgressSingleton::GetInstance();

    bool ShapeIsValid(ShapeResult &result);
    std::pair<std::string, std::string> AlterShape(std::string &shape_str);
    void Train(ShapeBuilder &builder);
    void Benchmark(ShapeBuilder &builder, ShapeResult &result, vector<size_t> mutations);



public:
    static void Test(ShapeBuilder &shape, size_t read_length, size_t iterations, std::vector<size_t> mutations, size_t num_threads, ShapeResult &result);
    static void Test(ShapeBuilder &shape, size_t read_length, size_t iterations, double ani, ShapeResult &result);

    ShapeFinder(FinderOptions options);

    void Find();

    static void Predict(ShapeBuilder &builder, bool *read, size_t read_length, bool *hm_pattern, size_t hm_size,
                 unordered_set<size_t> &mutations, unordered_set<size_t> &predictions);

    static void SetBit(size_t pos, bool to, uint32_t &pattern);

    static void GenerateHitMissPattern(bool *read, size_t read_length, bool *pattern, ShapeBuilder &builder);
};

#endif //VARKIT_SHAPEHISTORY_H
