//
// Created by fritsche on 11/08/2021.
//

#include <fstream>
//#include <semaphore.h>
#include "ShapeHistory.h"

ShapeResult::ShapeResult(std::string shape) : shape_(shape) {}

StatsContainer &ShapeResult::GetStatistic(std::string key) {
//    std::cout << "key: " << key << std::endl;
    if (not statistics.contains(key)) {
        statistics[key] = StatsContainer();
    }
    return statistics[key];
}

void ShapeResult::AddPosition(size_t pos, size_t hit_count) {
    if (pos >= positional.size()) {
        positional.resize(pos + 1);
    }
    positional[pos] += hit_count;
}


int ShapeResult::Compare(ShapeResult &other, Comparator function) {
    return function(*this, other);
}

void ShapeResult::Print(ostream &os) {
    os << shape_;

    for (auto sens : sensitivities) {
        os << '\t' << sens.first << ":";
        os << sens.second;
    }

    for (auto hit : hits) {
        os << '\t' << hit.first << ":";
        os << (double) hit.second / kmers[hit.first];
    }
    os << std::endl;
}

void ShapeResult::PrintMean(ostream &os, std::vector<double> selected_anis) {
    os << shape_;

    double ani_sens = 0;
    double ani_hits = 0;

    for (auto ani : selected_anis) {
        std::string sens_key = ShapeHistory::AniToString(ani) + "_sens";
        ani_sens += sensitivities[sens_key];
        ani_hits += hits[ShapeHistory::AniToString(ani)];


    }
    os << '\t' << "sensmean:";
    os << ani_sens/selected_anis.size();
    os << '\t' << "hitsmean:";
    os << ani_hits/selected_anis.size();


//    for (auto hit : hits) {
//        os << '\t' << hit.first << ":";
//        os << (double) hit.second / kmers[hit.first];
//    }
    os << std::endl;
}

void ShapeResult::CalculateStats() {
    for (auto entry : statistics) {
        if (entry.second.fp > 0) {
            exit(1);
        }
        sensitivities[entry.first + "_sens"] = entry.second.Sensitivity();
    }
}

void ShapeHistory::LoadHistory(std::string history) {
    if (Utils::exists(history)) {
        auto str = csv::load_file(history.c_str());
        auto parser = csv::make_parser(str, '\t');
        int col = 0;

        bool header = true;

        for (auto &&row : parser) {
            auto it = row.begin();

            if (header) {
                header = false;
//                while (it != row.end()) {
//                    std::string name = (*it).to_string();
//                    colnames.push_back(name);
//                    ++it;
//                }
                continue;
            }

            std::string shape = (*it).to_string();
            ShapeResult result{shape};

            double sensitivity;
            for (int i = 0; i < 10; i++) {
                sensitivity = (*(++it)).to_double();
            }
            col++;
        }
    }
}


void ShapeHistory::Write(std::string file) {
    bool write_header = not Utils::exists(file);

    std::ofstream ofs(file, std::ofstream::app);

    if (write_header) {
        ofs << "shape";
        for (auto& ani : anis) {
            std::string sens_key = AniToString(ani) + "_sens";
            ofs << '\t' << sens_key;
        }
        for (auto& ani : anis) {
            std::string hit_key = AniToString(ani) + "_hits";
            ofs << '\t' << hit_key;
        }
        ofs << std::endl;
    }
//    for (auto& ani : anis) {
//        ofs << '\t' << AniToString(ani) << "_sens";
//    }
//
//    for (auto& ani : anis) {
//        ofs << '\t' << AniToString(ani) << "_hits";
//    }
        auto& last = shapes[shapes.size() - 1];
        ofs << last.shape_;

        for (auto& ani : anis) {
//            std::cout << "ani: " << ani;
            std::string sens_key = AniToString(ani) + "_sens";
            ofs << '\t' << last.sensitivities[sens_key];
        }
        for (auto& ani : anis) {
            std::string hits_key = AniToString(ani);// + "_hits";
            ofs << '\t' << (double) last.hits[hits_key] / last.kmers[hits_key];
        }
    ofs << std::endl;

    ofs.close();
}


std::string ShapeHistory::AniToString(double ani) {
    return prefix + to_string((int) (ani * 100));
}

double ShapeHistory::StringToAni(std::string ani_str) {
    return stod(ani_str.substr(prefix.length())) / 100;
}

bool ShapeHistory::Empty() {
    return shapes.empty();
}

void ShapeHistory::Add(ShapeResult &res) {
    shapes.push_back(res);
    tested_shapes.insert(res.shape_);
}

bool ShapeHistory::HasShape(std::string shape_str) {
    return tested_shapes.contains(shape_str);
}

ShapeHistory::ShapeHistory() {}

//ShapeResult &ShapeHistory::GetBest() {
//    return shapes.at(0);
//}

ShapeResult &ShapeHistory::GetBest(Comparator fun) {
    ShapeResult& best = shapes[0];
//    std::cout << "compare_to: " << best.shape_ << std::endl;
//    for (auto test : best.sensitivities) {
//        std::cout << test.first << ", " << test.second << std::endl;
//    }

    for (auto i = 1; i < shapes.size(); i++) {
//        std::cout << "compare_to: " << shapes[i].shape_ << std::endl;
        if (fun(best, shapes[i]) < 0)
            best = shapes[i];
    }
    return best;
}


ShapeFinder::ShapeFinder(FinderOptions options) : options_(options) {
    history_.LoadHistory(options.previous_findings);
    progress_ = ProgressSingleton::GetInstance();
}

void ShapeFinder::Test(ShapeBuilder &shape, size_t read_length, size_t iterations, std::vector<size_t> mutations, size_t num_threads,
                       ShapeResult &result) {
    srand(time(0));

    size_t block_size = 500;
    size_t iteration_blocks = iterations / block_size;

//    if (progress_)
//        progress->Reset(iteration_blocks * mutations.size());

    size_t prog_chunk = 1;


    {
        thread_pool pool(num_threads);
        for (auto it_block = 0; it_block < iteration_blocks; it_block++) {
            for (auto mi = 0; mi < mutations.size(); mi++) {
                size_t m = mutations[mi];
//                std::cout << "mi: " << mi << " : " << mutations[mi] << " -> ani: " << (1 - ((double) m / read_length)) << std::endl;
                pool.enqueue_work([&shape, read_length, block_size, m, &result ]() {
                    Test(shape, read_length, block_size, 1 - ((double) m / read_length), result);
                });
            }
        }
    }
}

void ShapeFinder::SetBit(size_t pos, bool to, uint32_t &pattern) {
    if (to) {
        pattern |= 1 << (32 - 1 - pos);
    }
    else {
        pattern &= ~(1 << (32 - 1 - pos));
    }
}

void ShapeFinder::GenerateHitMissPattern(bool* read, size_t read_length, bool* pattern, ShapeBuilder& builder) {
    for (auto rpos = 0; rpos < read_length - builder.shape_size + 1; rpos++) {
        bool hit = true;
        for (auto spos = 0; spos < builder.shape_size; spos++) {
            if (read[rpos + spos] && !builder.shape[spos]) {
                hit = false;
                break;
            }
        }
        if (rpos < 0 || rpos > read_length - builder.shape_size) {
            std::cout << read_length << std::endl;
            std::cout << builder.pattern_size << std::endl;
            exit(3);
        }

        pattern[rpos] = hit;
    }
}

void ShapeFinder::Predict(ShapeBuilder &builder, bool* read, size_t read_length, bool* hm_pattern, size_t hm_size, std::unordered_set<size_t> &mutations, std::unordered_set<size_t> &predictions) {
    GenerateHitMissPattern(read, read_length, hm_pattern, builder);


    mutations.clear();
    for (int i = 0; i < read_length; i++)
        if (read[i])
            mutations.insert(i);

    uint32_t pattern = 0;


    // (Initialize) generate pattern (first digits)
    for (auto ppos = 0; ppos < builder.pattern_size; ppos++) {
        if (hm_pattern[ppos]) {
            SetBit(ppos, true, pattern);
        }
    }


    for (auto pattern_pos = 0; pattern_pos < hm_size - builder.pattern_size + 1; pattern_pos++) {
        // Update pattern
        SetBit(builder.pattern_size - 1, hm_pattern[builder.pattern_size - 1 + pattern_pos], pattern);

        if (pattern == 0 || pattern == ((-1u) << (32 - builder.pattern_size)))  {
            pattern <<= 1;
            continue;
        }

        size_t tmp = 0;

        if (builder.map2.contains(pattern)) {
            auto predicted = builder.map2[pattern];
            for (auto pi = 0; pi < predicted.Capacity() && predicted.mutations[pi] != 0xFF; pi++)
                predictions.insert(pattern_pos + predicted.mutations[pi]);
        }

        pattern <<= 1;
    }
}


std::mutex test_tp_mtx2;
void ShapeFinder::Test(ShapeBuilder &shape, size_t read_length, size_t iterations, double ani, ShapeResult &result) {
//    std::cout << "ani: " << ani << std::endl;
    bool *read = new bool[read_length];
    size_t *counts = new size_t[read_length];


    size_t hm_size = read_length - shape.shape_size + 1;
    bool *hm_pattern = new bool[hm_size];

    // init?
    memset(read, 0, read_length * sizeof(bool));
    memset(hm_pattern, 0, hm_size * sizeof(bool));
    memset(counts, 0, read_length * sizeof(size_t));

    std::unordered_set<size_t> mutation_set;

    size_t fp = 0;
    size_t tp = 0;
    size_t fn = 0;

    size_t total_hits = 0;
    size_t total_kmers = iterations * (read_length - shape.shape_size + 1);

    size_t prog_counter = 0;
    for (auto iteration = 0; iteration < iterations; iteration++) {
        double synonymous_probability = 0.8;
        size_t mutations = 0;
        while (!mutations)
//        while (mutations != 1)
            mutations = ShapeBuilder::GenerateReadRandom(read, read_length, 1-ani, synonymous_probability);
//        ShapeBuilder::GenerateReadRandom(read, read_length, 5);

        shape.Test(read, read_length, hm_pattern, hm_size, mutation_set, tp, fp, fn, counts, total_hits);

//        std::cout << "hits(" << mutations << "," << ani << "): " << total_hits << "/" << (iteration + 1) * (read_length - shape.shape_size + 1) << " = " << (double)total_hits/((iteration + 1) * (read_length - shape.shape_size + 1)) << std::endl;

        bool debug = false;
        if (debug || fp != 0) {
            std::cout << ShapeBuilder::ReadToString(read, read_length) << std::endl;
            for (auto m : mutation_set) {
                std::cout << (int) m << ",";
            }
            std::cout << std::endl;

            std::unordered_set<size_t> predictions;
            Predict(shape, read, read_length, hm_pattern, hm_size, mutation_set,predictions);
            std::cout << "predictions: " << std::endl;
            for (auto m : predictions) {
                std::cout << (int) m << ",";
            }
//            exit(9);
            std::string stop;
            std::cin >> stop;

        } else {
//            std::cout << "fp: " << fp << ", tp: " << tp << std::endl;
        }
    }

    assert (fp == 0);

    {
        std::lock_guard<std::mutex> sens_guard(test_tp_mtx2);
        auto& statistics = result.GetStatistic(ShapeHistory::AniToString(ani));
        statistics.fp += fp;
        statistics.tp += tp;
        statistics.fn += fn;



        result.hits[ShapeHistory::AniToString(ani)] += total_hits;
        result.kmers[ShapeHistory::AniToString(ani)] += total_kmers;

        for (int i = 0; i < read_length; i++) {
            result.AddPosition(i, counts[i]);

        }
    }

//    auto progress = ProgressSingleton::GetInstance();
//    if (progress)
//        progress->UpdateAdd(1);

    delete[] read;
    delete[] counts;
    delete[] hm_pattern;
}

std::pair<std::string, std::string> ShapeFinder::AlterShape(string &shape_str) {
    std::string operation = "";
    std::string new_shape_str = ShapeBuilder::AlterShape(shape_str, 0.5, options_.take, options_.space, operation);
    while (history_.HasShape(new_shape_str)) {
        new_shape_str = ShapeBuilder::AlterShape(shape_str, 0.5, options_.take, options_.space, operation);
    }

    return { new_shape_str, operation };
}

bool ShapeFinder::ShapeIsValid(ShapeResult &result) {
    size_t sc = 0, tc = 0; //space and take counter to assure validity of shape
    for (auto ci = 0; ci < result.shape_.size(); ci++) {
        char c = result.shape_[ci];
        c == 'X' ? tc++ : sc++;
    }

    return (tc == options_.take and sc == options_.space);
}

void ShapeFinder::Train(ShapeBuilder &builder) {
//        builder.progress = progress_;
    builder.pattern_size = options_.pattern_size;
    builder.limit = options_.training_max;
    builder.TrainByPatternTP(options_.threads);
}

void ShapeFinder::Benchmark(ShapeBuilder &builder, ShapeResult &result, vector<size_t> mutations) {
    Test(builder, 100, options_.iterations, mutations, options_.threads, result);
}

void ShapeFinder::Find() {
    size_t read_length = 100;
    vector<size_t> mutations = {1, 2, 3, 4, 5, 6, 7, 8, 9 };

    for (auto mut : mutations) {
        double ani = (1 - (double) mut/read_length);
        history_.anis.push_back(ani);
    }

    std::vector<std::string> compare_cols;
    for (auto mutation : { 1, 2, 3, 4, 5 }) {
//            std::cout << "ani: " << (1 - ((double)mutation/read_length)) << " col: " << history_.AniToString(1 - ((double)mutation/read_length)) << std::endl;
        compare_cols.push_back(history_.AniToString(1 - ((double)mutation/read_length)));
    }
    auto attributes_to_value = [compare_cols] (ShapeResult &result) {
        double sensitivity_sum = 0;
        for (auto v : compare_cols) {
//                std::cout << "header: " << v + "_sens" << std::endl;
//                std::cout << "... "<< result.sensitivities[v + "_sens"] << std::endl;
            sensitivity_sum += result.sensitivities[v + "_sens"];
        }
//            std::cout << sensitivity_sum << std::endl;
        return sensitivity_sum/compare_cols.size();
    };

    std::function<int(ShapeResult&, ShapeResult&)> comparator = [attributes_to_value] (ShapeResult &result1, ShapeResult &result2) {
        auto a = attributes_to_value(result1);
        auto b = attributes_to_value(result2);

//            std::cout << "a: " << a << ", b: " << b << std::endl;
        return a > b ? 1 : a == b ? 0 : -1;
    };

    ShapeResult shape {
            history_.Empty() ?
            ShapeResult(ShapeBuilder::RandomShape(options_.take, options_.space)) :
            history_.GetBest(comparator) };

    size_t iterations = 0;


    for (auto mutation : mutations) {
        history_.sens_colnames.push_back(history_.AniToString((int) (1 - (double)mutation/read_length)) + "_sens");
        history_.hit_colnames.push_back(history_.AniToString((int) (1 - (double)mutation/read_length)) + "_hits");
    }



    auto continue_search = [&] {
        return iterations < options_.iterations;
    };




    while (continue_search()) {
        assert(ShapeIsValid(shape));

        ShapeBuilder builder(shape.shape_);
        Train(builder);

        Benchmark(builder, shape, mutations);

        shape.CalculateStats();
//            shape.Print(std::cout);
        shape.PrintMean(std::cout, { 0.95, 0.96, 0.97, 0.98, 0.99 });

        history_.Add(shape);

        auto& best = history_.GetBest(comparator);

        std::cout << "Best (" << attributes_to_value(best) << ")" << std::endl;
//            best.Print(std::cout);
        best.PrintMean(std::cout, { 0.95, 0.96, 0.97, 0.98, 0.99 });

        std::cout << "____" << std::endl;

        auto [new_shape_c, operation_c] = AlterShape(shape.shape_);
        std::string new_shape = new_shape_c;
        std::string operation = operation_c;
        shape = ShapeResult(new_shape);


        history_.Write(options_.previous_findings);
    }
}
