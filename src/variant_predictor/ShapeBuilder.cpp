//
// Created by fritsche on 12/04/2021.
//

#include "ShapeBuilder.h"
#include "varkit_config.h"
#include "constants.h"
#include "version.h"

using namespace ds;

ShapeBuilder::ShapeBuilder(std::string shape_str) {
    this->shape = StringToShape(shape_str);
    this->shape_size = shape_str.length();
}

ShapeBuilder::~ShapeBuilder() {
    map2.clear();
    delete[] shape;
}

void ShapeBuilder::LoadArray(std::string path) {
    ifstream ifs(path, ifstream::in);

    std::string line;


    ifs.read((char *) &shape_size, sizeof(shape_size));
    ifs.read((char *) &pattern_size, sizeof(pattern_size));
    ifs.read((char *) &limit, sizeof(limit));
    ifs.read((char *) &lower_limit, sizeof(lower_limit));

    delete[] shape;
    shape = new bool[shape_size];
    ifs.read((char *) shape, shape_size * sizeof(bool));

    auto map_array_length = 1llu << pattern_size;

    std::cout << "map_array_length: " << map_array_length << std::endl;
    std::cout << "MB requirement for map (pattern size: " << pattern_size << "): " << ((double)(map_array_length * 4)/(1024*1024)) << std::endl;

    map_array = new Mutations[map_array_length];

    while (ifs.peek() != EOF) {
        uint64_t pattern;
//        ifs.read((char*) &pattern, 8);
//        uint32_t pattern;
//        size_t value;
        Mutations value;
        ifs.read((char*) &pattern, 8);
        ifs.read((char*) &value, sizeof(Mutations));

        uint32_t trunc_pattern = pattern >> (64llu - pattern_size);

        map_array[trunc_pattern] = value;
    }
//    exit(9);

    ifs.close();
}


void ShapeBuilder::Load(std::string path) {
    ifstream ifs(path, ifstream::in);

    std::string line;


    ifs.read((char *) &shape_size, sizeof(shape_size));
    ifs.read((char *) &pattern_size, sizeof(pattern_size));
    ifs.read((char *) &limit, sizeof(limit));
    ifs.read((char *) &lower_limit, sizeof(lower_limit));

    delete[] shape;
    shape = new bool[shape_size];
    ifs.read((char *) shape, shape_size * sizeof(bool));


    while (ifs.peek() != EOF) {
        uint64_t pattern;
//        ifs.read((char*) &pattern, 8);
//        uint32_t pattern;
//        size_t value;
        Mutations value;
        ifs.read((char*) &pattern, 8);
        ifs.read((char*) &value, sizeof(Mutations));

        uint32_t trunc_pattern = pattern >> 32llu;

//        std::cout << BitUtils::ToBitString(trunc_pattern, pattern_size) << " -> " << value.ToString() << std::endl;
        map2.insert({ trunc_pattern, value });
    }
//    exit(9);

    ifs.close();
}

void ShapeBuilder::Save(std::string path) {
    using namespace constants;

    ofstream ofs(path, ofstream::out);

    Version version(varkit_VERSION_MAJOR, varkit_VERSION_MINOR, varkit_VERSION_PATCH);
    ofs.write((char *) &version, sizeof(version));

    size_t pattern_sizeof = sizeof(map2.begin()->first);
    size_t mutations_sizeof = sizeof(map2.begin()->second);

    ofs.write((char *) &pattern_sizeof, sizeof(pattern_sizeof));
    ofs.write((char *) &mutations_sizeof, sizeof(mutations_sizeof));

    ofs.write((char *) &shape_size, sizeof(shape_size));
    ofs.write((char *) &pattern_size, sizeof(pattern_size));
    ofs.write((char *) &limit, sizeof(limit));
    ofs.write((char *) &lower_limit, sizeof(lower_limit));

    ofs.write((char *) shape, shape_size * sizeof(bool));

    for (auto pair : map2) {
//        std::cout << ShapeBuilder::ToBitString(pair.first) << '\t' << pair.second.ToString() << std::endl;
        ofs.write((char *) &pair.first, sizeof(pair.first));
        ofs.write((char *) &pair.second, sizeof(pair.second));
    }

    ofs.close();
}





//typedef tsl::sparse_map<size_t, Mutations> PatternMap;
void ShapeBuilder::Train(const size_t threads) {
    const size_t max_pattern = 1llu << pattern_size;
    const size_t block_size = max_pattern / (threads*20);
    const size_t blocks = max_pattern / block_size;
    PatternMap pattern_maps[blocks];
    PatternMap final_map;

    if (progress) progress->Reset(blocks);


    {
        thread_pool pool(threads);
        for (auto i = 0; i < blocks; i++) {
            pool.enqueue_work([this, &pattern_maps, i, block_size, blocks, max_pattern]() {
                size_t lower = i * block_size;
                // check if last block
                size_t upper = i == blocks - 1 ? max_pattern : (i+1) * block_size;
                Train(pattern_maps[i], lower, upper);
            });
        }
    }


    // Insert pattern from tmp storage into map
    size_t pattern_count = 0;
    for (auto map : pattern_maps) {
        for (auto pair : map) {
            pattern_count++;

            uint64_t key = pair.first << (64 - pattern_size);
            map2.insert( { key, pair.second } );
        }
    }
}

void ShapeBuilder::Train(PatternMap& pattern_map, size_t from, size_t to) {
    size_t read_length = pattern_size + shape_size - 1;

    bool read_dummy[read_length];
    bool read[read_length];
    bool hm_pattern[pattern_size];

    size_t pot_mut_counter_loc[100];
    memset(pot_mut_counter_loc, 0, 100* sizeof(size_t));

    size_t mutations[100];

    size_t total_count = 0;
    size_t valid_count = 0;



    for (size_t pattern = from; pattern < to; pattern++) {
        Mutations value = TrainByPatternOMP(pattern, read, read_dummy, hm_pattern, mutations, total_count, valid_count, pot_mut_counter);

        if (!value.Empty()) {
            pattern_map.insert({pattern, value});
        }

    }


    {
        std::lock_guard<std::mutex> lock(pot_mut_mutex);
        for (auto i = 0; i < 100; i++) {
            pot_mut_counter[i] += pot_mut_counter_loc[i];
        }
    }

    if (progress) progress->UpdateAdd(1);
}


std::mutex error_log_mtx;
void ShapeBuilder::Test(bool* read, size_t read_length, bool* hm_pattern, size_t hm_size, std::unordered_set<size_t> &mutations, size_t &tp, size_t &fp, size_t &fn, size_t* fn_pos, size_t &total_hits) {
    GenerateHitMissPattern(read, read_length, hm_pattern);

    for (auto i = 0; i < hm_size; i++)
        total_hits += hm_pattern[i];

    GetMutations(read, read_length, mutations);
    uint64_t pattern = 0;

    std::unordered_set<size_t> predictions;


    // (Initialize) generate pattern (first digits)
    for (auto ppos = 0; ppos < pattern_size; ppos++) {
        if (hm_pattern[ppos]) {
            SetBit(ppos, true, pattern);
        }
    }


    for (auto pattern_pos = 0; pattern_pos < hm_size - pattern_size + 1; pattern_pos++) {
        // Update pattern
        SetBit(pattern_size - 1, hm_pattern[pattern_size - 1 + pattern_pos], pattern);

        if (pattern == 0 || pattern == ((-1llu) << (64 - pattern_size))) {
            pattern <<= 1;
            continue;
        }

        size_t tmp = 0;


        if (map2.contains(pattern)) {
            auto predicted = map2[pattern];
            for (auto pi = 0; pi < predicted.Capacity() && predicted.mutations[pi] != 0xFF; pi++) {
                assert(pattern_pos + predicted.mutations[pi] < read_length);
                predictions.insert(pattern_pos + predicted.mutations[pi]);
            }
        }

        pattern <<= 1;
    }

    std::vector<size_t> intersection;
    ShapeBuilder::Intersect(predictions, mutations, intersection);

    tp += intersection.size();
    fp += predictions.size() - intersection.size();
    fn += mutations.size() - intersection.size();

    bool print_ani = true;
    if (print_ani) {
        double ani = (double) mutations.size() / (double) read_length;
        size_t hm_size = read_length - shape_size + 1;
        size_t count_hits = 0;
        for (int i = 0; i < hm_size; i++) {
            count_hits += hm_pattern[i];
        }
    }

    if (fp != 0) {
        std::cout << "\ndebug_______________" << map2.size() << std::endl;
        std::cout << ShapeBuilder::ReadToString(read, read_length) << std::endl;
        std::cout << "predictions: " << std::endl;
        for (auto m : predictions) {
            std::cout << (int) m << ",";
        }
        std::cout << "\nreal: " << std::endl;
        for (auto m : mutations) {
            std::cout << (int) m << ",";
        }
        std::cout << "\t____debugend" << std::endl;
    }


    if (fn_pos) {
        std::vector<size_t> a_without_b;
        ShapeBuilder::AwithoutB(mutations, predictions, a_without_b);
        for (auto e : a_without_b)
            fn_pos[e]++;
    }
}



void ShapeBuilder::TestPattern(size_t mutation_num, size_t iterations, size_t &tp, size_t &fp, size_t &fn) {
    size_t max_distance = shape_size + pattern_size - 1;
//    size_t* mutations = new size_t[mutation_num];

    size_t read_length = (mutation_num+1) * max_distance;

    size_t* distances = new size_t[mutation_num];
    bool* read = new bool[read_length];
    bool* hm_pattern = new bool[read_length];


    memset(distances, 0, mutation_num * sizeof(size_t));
    memset(read, 0, read_length * sizeof(bool));
    memset(hm_pattern, 0, read_length * sizeof(bool));

    std::unordered_set<size_t> mutations;

    size_t rlen = 0;

    srand(time(0));

    ProgressBar bar(iterations);
    size_t prog = 0;
    while (iterations-- > 0) {
        std::cout << "it" << iterations << std::endl;
        if (prog++ % 1000 == 0) {
            bar.Update(prog);
        }
        rlen = GenerateReadRandom(read, read_length, mutation_num, distances, max_distance);
        size_t hm_size = rlen - shape_size + 1;
//        std::cout << "random read: " << std::endl;
//        std::cout << ShapeBuilder::ToString(mutations) << std::endl;
//        std::cout << ReadToString(read, rlen) << std::endl;
        size_t total_hits = 0;
        Test(read, rlen, hm_pattern, hm_size, mutations, tp, fp, fn, nullptr, total_hits);
    }

    delete[] read;
    delete[] distances;
    delete[] hm_pattern;
}


void ShapeBuilder::GenerateHitMissPattern(bool* read, size_t read_length, bool* pattern) {
    for (auto rpos = 0; rpos < read_length - shape_size + 1; rpos++) {
        bool hit = true;
        for (auto spos = 0; spos < shape_size; spos++) {
            if (read[rpos + spos] && !shape[spos]) {
                hit = false;
                break;
            }
        }
        if (rpos < 0 || rpos > read_length - shape_size) {
            std::cout << read_length << std::endl;
            std::cout << pattern_size << std::endl;
            exit(3);
        }

        pattern[rpos] = hit;
    }
}

void ShapeBuilder::GenerateRead(bool* &read, size_t &read_length, size_t* distances, size_t mutation_count) {
    size_t margin = shape_size;

    size_t begin_margin = shape_size + pattern_size - 2;
    size_t end_margin = shape_size;

    size_t required_length = begin_margin + end_margin;

    // If more than one mutation then update required length
    if (mutation_count > 1) {
        size_t sum = std::accumulate(distances, distances + mutation_count - 1, mutation_count);
        required_length += sum;
    }

    if (read_length < required_length) {
        std::cerr << "read length < required length in Generate Read" << std::endl;
        exit(3);
//        std::cout << (read_length < required_length) << std::endl;
//        std::cout << "rl: " << read_length << std::endl;
//        std::cout << "reql: " << required_length << std::endl;
        delete[] read;

        if (read_length == 0)
            read_length = 1 << 10;

        while (read_length < required_length) {
            read_length *= 2;
        }
        read = new bool[read_length];
    }

    // reset read
    memset(read, false, read_length * sizeof(bool));

    size_t pos = begin_margin;
    read[pos] = true;

    for (auto i = 0; i < mutation_count - 1; i++) {
        read[pos + distances[i] + 1] = true;
        pos += distances[i] + 1;
    }
}


size_t ShapeBuilder::GenerateReadRandom(bool* &read, size_t &read_length, size_t mutation_num, size_t* distances, size_t max_dist) {
    size_t distance_size = mutation_num;

    srand(time(NULL));

    for (auto i = 0; i < distance_size; i++) {
        distances[i] = rand() % max_dist;
    }

    memset(read, false, read_length * sizeof(bool));

    int64_t pos = -1;
    for (auto i = 0; i < distance_size; i++) {
        read[pos + distances[i] + 1] = true;
        pos += distances[i] + 1;
    }

    return pos + max_dist;
}


size_t ShapeBuilder::GenerateReadRandom(bool *&read, size_t &read_length, size_t mutation_num) {

    memset(read, false, read_length * sizeof(bool));

    size_t mutation_counter = 0;
    while (mutation_counter != mutation_num) {
        auto pos = rand() % read_length;
        if (read[pos]) continue;
        read[pos]= true;
        mutation_counter++;
    }
    return 0;
}

size_t ShapeBuilder::GenerateReadRandom(bool *&read, size_t &read_length, double mutation_rate, double synonymous_prob) {
    int offset = rand() % 3;

    double non_synonymous_prob { (1 - synonymous_prob) / 2 };

    memset(read, false, read_length * sizeof(bool));

    int mutation_count = 0;

    for (int i = 0; i < read_length; i++) {
        double prob = mutation_rate * 3 * (((i % 3) == offset) * synonymous_prob + ((i % 3) != offset) * non_synonymous_prob);

//        std::cout << "prob:" << prob << std::endl;
        read[i] = ShapeBuilder::RandomCoinFlip(prob);
        mutation_count += read[i];
    }

    return mutation_count;
}


// Utility functions

void ShapeBuilder::SetBit(size_t pos, bool to, uint32_t &pattern) {
    if (to) {
        pattern |= 1 << (32 - 1 - pos);
    }
    else {
        pattern &= ~(1 << (32 - 1 - pos));
    }
}

void ShapeBuilder::SetBit(size_t pos, bool to, uint64_t &pattern) {
    if (to) {
        pattern |= 1llu << (64 - 1 - pos);
    }
    else {
        pattern &= ~(1llu << (64 - 1 - pos));
    }
}

bool *ShapeBuilder::StringToShape(std::string shape_str) {

    bool* shape = new bool[shape_str.length()];

    for (auto i = 0; i < shape_str.length(); i++) {
        if (shape_str[i] == take) {
            // take is false
            shape[i] = false;
        } else if (shape_str[i] == space) {
            // space is true
            shape[i] = true;
        } else {
            std::cerr << "Shape string may only contain " << take << " and " << space << std::endl;
            exit(0);
        }
    }

    return shape;
}

std::string ShapeBuilder::ShapeToString(bool *shape, size_t size) {
    std::string shape_str(size, 'X');

    for (auto i = 0; i < size; i++)
        shape_str[i] = take * !shape[i] + space * shape[i];

    return shape_str;
}

std::string ShapeBuilder::ReadToString(bool *read, size_t size) {
    std::string read_str = std::string(size, '_');

    for (int i = 0; i < size; i++)
        read_str[i] = take * read[i] + space * !read[i];

    return read_str;
}


size_t ShapeBuilder::GetMutations(bool *read, size_t read_length, std::unordered_set<size_t> &pos) {
    pos.clear();
    size_t mutation_num = 0;

    for (int i = 0; i < read_length; i++) {
        if (read[i]) {
            pos.insert(i);
            mutation_num++;
        }
    }

    return mutation_num;
}



//void ShapeBuilder::TrainByPatternMT(uint32_t* patterns, size_t size, Mutations* tmp) {
void ShapeBuilder::TrainByPatternMT(uint64_t* patterns, size_t size, Mutations* tmp) {
    size_t read_length = pattern_size + shape_size - 1;

    bool read_dummy[read_length];
    bool read[read_length];
    bool hm_pattern[pattern_size];


    size_t pot_mut_counter_loc[100];
    memset(pot_mut_counter_loc, 0, 100* sizeof(size_t));

    size_t mutations[100];

    size_t total_count = 0;
    size_t valid_count = 0;


    for (size_t pi = 0; pi < size; pi++) {
        size_t p = patterns[pi];

        Mutations value = TrainByPatternOMP(p, read, read_dummy, hm_pattern, mutations, total_count, valid_count, pot_mut_counter);


        if (!value.Empty()) {
            tmp[p] = value;
        }
    }


    {
        std::lock_guard<std::mutex> lock(pot_mut_mutex);
        for (auto i = 0; i < 100; i++) {
            pot_mut_counter[i] += pot_mut_counter_loc[i];
        }
    }

    if (progress) progress->UpdateAdd(1);
}

void ShapeBuilder::TrainByPatternMT(size_t threads_num) {
    size_t pattern_perm = 1 << pattern_size;

//    size_t *tmp_storage = new size_t[pattern_perm];
//    memset(tmp_storage, 0xFF, pattern_perm * sizeof(size_t));
    Mutations* tmp_storage = new Mutations[pattern_perm];


    std::vector<std::thread> threads;


    size_t lower = 0;
    size_t upper = 1;
    size_t chunk = pattern_perm / threads_num;


    size_t bucket_size = chunk + 1;

//    uint32_t* pattern_sizes = new uint32_t[threads_num];
//    uint32_t* patterns = new uint32_t[bucket_size * threads_num];
//    memset(pattern_sizes, 0, threads_num * sizeof(uint32_t));
//    memset(patterns, 0, bucket_size * threads_num * sizeof(uint32_t));


    uint64_t* pattern_sizes = new uint64_t[threads_num];
    uint64_t* patterns = new uint64_t[bucket_size * threads_num];
    memset(pattern_sizes, 0, threads_num * sizeof(uint64_t));
    memset(patterns, 0, bucket_size * threads_num * sizeof(uint64_t));

    std::cout << "bucket_size: " << bucket_size << std::endl;

    size_t offset = 0;
    for (int p = 1; p < pattern_perm; p++) {
        size_t bucket = (p + offset) % threads_num;
        size_t bucket_offset = bucket * bucket_size;
        patterns[bucket_offset + pattern_sizes[bucket]++] = p;
        if ((p % threads_num) == 0) bucket_offset++;
    }


    for (auto t = 0; t < threads_num; t++) {
        lower = t * chunk;
        if (!lower) lower = 1;
        size_t upper = t * chunk + chunk;

        if (t == threads_num - 1) upper = pattern_perm;

        // Give each thread patterns distributed across the pattern space
        threads.emplace_back([this, t, patterns, bucket_size, pattern_sizes, &tmp_storage]() {
            TrainByPatternMT(patterns + (t * bucket_size), pattern_sizes[t], tmp_storage);
        });
    }


    for (auto &t : threads)
        t.join();

    std::cout << map2.size() << std::endl;


    // Insert from tmp storage into map
    size_t pattern_count = 0;
    for (auto p = 0; p < pattern_perm; p++) {
        if (!tmp_storage[p].Empty()) {
//            std::cout << ToString(tmp_storage[p]) << std::endl;
            pattern_count++;
            uint32_t key = p << (32 - pattern_size);
            map2.insert( { key, tmp_storage[p] });
        }
    }
    delete[] tmp_storage;
    std::cout << "pattern_count: " << pattern_count << std::endl;
}


void ShapeBuilder::TrainByPattern() {
//    size_t pattern_max = 1 << pattern_size;
//
////    ProgressBar bar(pattern_max);
////    size_t prog_counter = 0;
//
//
//    bool verbose = false;
//
//    std::cout << "size: " << map2.size() << std::endl;
//
//    size_t min_1 = 0;
//
//
//    Benchmark bm("total");
//    bm.start();
//
//    omp_set_num_threads(1);
//
//#pragma omp parallel for
//    for (uint32_t p = 1; p < pattern_max; p++) {
//
//        size_t read_length = pattern_size + shape_size - 1;
//
//        bool read_dummy[read_length];
//        bool read[read_length];
//
////        size_t value = TrainByPatternOMP(p, read, read_dummy);
//
//        size_t value = 0;
//        if (value) {
//#pragma omp critical(insert)
//            {
//                uint32_t key = p << (32 - pattern_size);
//                map2.insert( { key, value } );
//            }
//        }
//
////        double prog = ((double) p / pattern_max);
////        if ((++prog_counter % 250) == 0) {
////            bar.Update(prog_counter);
////        }
//    }
//
//    std::cout << "no valid patterns: " << count_no_valid_patterns << std::endl;
//    std::cout << "no valid ration: " << ((double) count_no_valid_patterns / count_total) << std::endl;
//    std::cout << "ratio: " << ((double) count_below / count_total) << std::endl;
//    std::cout << "valid patterns: " << ((double) (count_below-count_no_valid_patterns) / count_total) << std::endl;
//
//    bm.stop();
//    bm.printResults();
//
//
//    std::cout << "min 4 bits set: " << min_1 << std::endl;
//    std::cout << "size: " << map2.size() << std::endl;

}

void ShapeBuilder::TrainPatternBased() {
//    size_t pattern_max = 1 << pattern_size;
//
//    ProgressBar bar(pattern_max);
//    size_t prog_counter = 0;
//
//    size_t read_length = pattern_size + shape_size - 1;
//    bool* read_dummy = new bool[read_length];
//
//    bool verbose = false;
//
//    std::cout << "size: " << map2.size() << std::endl;
//
//    size_t min_1 = 0;
//
//    bool read[read_length];
//
//    Benchmark bm("total");
//    bm.start();
//
//
//    for (uint32_t p = 1; p < pattern_max; p++) {
//        TrainByPattern(p, read, read_dummy);
//
//        double prog = ((double) p/pattern_max);
//
//        if ((++prog_counter % 250) == 0)  {
//            bar.Update(prog_counter);
//        }
//
//        continue;
//
//        if (CountSetBits(p) < 6) continue;
//
//
//        uint32_t key = p << (32-pattern_size);
//        if (ShapeBuilder::InformationContent(key, pattern_size) < 9)  {
//            continue;
//        }
//
//        min_1++;
//
//
//
//        memset(read_dummy, true, read_length * sizeof(bool));
//
//        for (size_t pattern_pos = 0; pattern_pos < pattern_size; pattern_pos++) {
//            if (!((p >> (pattern_size - 1 - pattern_pos)) & 1)) continue;
//
//            for (size_t shape_pos = 0; shape_pos < shape_size; shape_pos++) {
//                read_dummy[pattern_pos + shape_pos] = shape[shape_pos] && read_dummy[pattern_pos + shape_pos];
//            }
//        }
//
//        std::set<size_t> muts;
//
//        uint32_t pc = p;
//        size_t ppos = (bool) pc * 1;
//        while ((pc /= 2))
//            ppos++;
//
//        auto ffl = (pattern_size - ppos);
//
//        if (verbose) {
//            std::cout << "first from left: " << ffl << std::endl;
//            std::cout << ReadToString(read_dummy, read_length) << " dummy read" << std::endl;
//        }
//
//        size_t m_margin = 5;
//        if (ppos > 1) {
//            for (int rpos = std::min((pattern_size - ppos), m_margin); rpos < read_length - m_margin; rpos++) {
//                if (read_dummy[rpos]) muts.insert(rpos);
//            }
//        }
//
//
//        size_t mut_int = 0;
//        size_t mut_cnt = 0;
//        for (auto m : muts) {
//            if (m < ffl) continue;
//            if (mut_cnt == 8) break;
//            Add(mut_int, m);
//            mut_cnt++;
//        }
////        std::cout << ToBitString(key) << " pattern key, information content: " << ShapeBuilder::InformationContent(key, pattern_size) << std::endl;
////        std::cout << "insert: " << ToBitString(p, pattern_size) << " " << ToString(mut_int) << std::endl;
//        map2.insert({ key, mut_int });
//
//        if (mut_int == 0) continue;
//
//        // Test instantly
//        bool test[150];
//        memset(test, 0, 150);
//        size_t offset = 0;
//        for (auto m : muts)
//            test[offset + m] = true;
//
//        bool hm_pattern[150];
//        size_t fn_pos[150];
//        size_t fp, tp, fn;
//        std::unordered_set<size_t> set;
//
//
//        Test(test, 150, hm_pattern, 150, set, tp, fp, fn, fn_pos);
//
//        if (verbose) {
//            std::cout << ShapeToString(shape, shape_size) << " shape" << std::endl;
//            std::cout << ToBitString(p, pattern_size) << " pattern" << std::endl;
//            std::cout << ReadToString(read_dummy, read_length) << " read" << std::endl;
//            if (muts.size() > 8) continue;
//            for (auto m : muts) {
//                std::cout << m << ", ";
//            }
//            std::cout << std::endl;
//
//            std::cout << "________________________________________MAP" << std::endl;
//            for (auto pair : map2) {
////                std::cout << ToBitString(pair.first) << " " << ToString(pair.second) << std::endl;
//                std::cout << ToBitString(pair.first) << " " << pair.second.ToString() << std::endl;
//            }
//
//
//            std::string stop;
//            std::cin >> stop;
//        } else {
//            bar.Update(++prog_counter);
//        }
//
//    }
//
//
//    std::cout << "no valid patterns: " << count_no_valid_patterns << std::endl;
//    std::cout << "no valid ration: " << ((double) count_no_valid_patterns / count_total) << std::endl;
//    std::cout << "ratio: " << ((double) count_below / count_total) << std::endl;
//    std::cout << "valid patterns: " << ((double) (count_below-count_no_valid_patterns) / count_total) << std::endl;
//
//    bm.stop();
//    bm.printResults();
//
//
//    std::cout << "min 4 bits set: " << min_1 << std::endl;
//    std::cout << "size: " << map2.size() << std::endl;
//
//    delete[] read_dummy;
}

size_t ShapeBuilder::GetMutations(bool *read, size_t read_length, size_t* mutations, size_t mutations_size) {
//    *((size_t*)mutations) = 0;
    memset(mutations, 0, mutations_size * sizeof(size_t));
    size_t mutation_num = 0;

    for (int i = 0; i < read_length; i++) {
        if (read[i]) {
            mutations[mutation_num++] = i;
        }
    }

    return mutation_num;
}


size_t ShapeBuilder::ChooseN(size_t k, size_t n, size_t threshold) {
    n++;
    while (Utils::nChoosek(--n, k) > threshold);
    return n;
}



Mutations ShapeBuilder::TrainByPatternOMP(uint32_t pattern, bool *read, bool *potential_mutations, bool *hm_pattern, size_t *mutations, size_t &total_count, size_t &valid_count, size_t* mut_counter) {
    size_t dummy_read_size = pattern_size + shape_size - 1;

    // reset
    memset(potential_mutations, true, dummy_read_size);
    memset(read, false, dummy_read_size);
    memset(hm_pattern, false, pattern_size);
    memset(mutations, 0, 100 * sizeof(size_t));

    size_t core_mutations_count[dummy_read_size];
    memset(core_mutations_count, 0, dummy_read_size * sizeof(size_t));

    uint32_t original_p = pattern << (32 - pattern_size);

//    std::cout << "Pattern:" << pattern << std::endl;
//    std::cout << BitUtils::ToBitString(pattern, 16) << std::endl;
//    std::cout << BitUtils::ToBitString(original_p, 16) << std::endl;

    for (size_t pattern_pos = 0; pattern_pos < pattern_size; pattern_pos++) {
        // if pattern at pos is 0 we do nothing
        if (!((pattern >> (pattern_size - 1 - pattern_pos)) & 1)) continue;

        for (size_t shape_pos = 0; shape_pos < shape_size; shape_pos++) {
            potential_mutations[pattern_pos + shape_pos] = shape[shape_pos] && potential_mutations[pattern_pos + shape_pos];
        }
    }

//    std::cout << ShapeBuilder::ReadToString(potential_mutations, dummy_read_size) << std::endl;

    size_t potential_mutation_count = 0;
    for (auto i = 0; i < dummy_read_size; i++) {
        potential_mutation_count += 1 * potential_mutations[i];
    }

    total_count++;

    // IMPORTANT
    size_t perm_threshold = limit;

    mut_counter[potential_mutation_count]++;

    Mutations value;

//    std::cout << "proceed: " << potential_mutation_count << std::endl;
    if (potential_mutation_count > perm_threshold || potential_mutation_count <= lower_limit || potential_mutation_count == 0)
        return value;
//    std::cout << "YES" << std::endl;

    size_t permutations = (size_t) pow(2, potential_mutation_count);

    GetMutations(potential_mutations, dummy_read_size, mutations, 100);


    for (auto m = 0; m < potential_mutation_count; m++) {
        read[mutations[m]] = true;
    }

    size_t equal_pattern_count = 0;

    // count mutations



    for (auto i = 0; i < permutations; i++) {
        read[mutations[potential_mutation_count - 1]] = !read[mutations[potential_mutation_count - 1]];
        for (int m = potential_mutation_count-2; m >= 0; m--) {
            int interval = pow(2, m + 1);
            auto pos = mutations[m];

            if ((i % interval) == 0)
                read[pos] = !read[pos];
        }

        if (i == 0) continue;

        //__________________________________________________________________________________
        // DO STH WITH READ
        //__________________________________________________________________________________

        GenerateHitMissPattern(read, dummy_read_size, hm_pattern);
        uint32_t p = 0;
        for (auto ppos = 0; ppos < pattern_size; ppos++) {
            if (hm_pattern[ppos]) {
                SetBit(ppos, true, p);
            }
        }

        if (p == original_p) {
            for (auto m = 0; m < potential_mutation_count; m++) {
                auto mpos = mutations[m];
                core_mutations_count[mpos] += read[mpos];
            }

            equal_pattern_count++;
        }
    }
    if (equal_pattern_count == 0) {
//        return (-1llu);
        return value;
    }
    valid_count++;

    for (auto m = 0; m < potential_mutation_count; m++) {
        auto mpos = mutations[m];
        if (core_mutations_count[mpos] == equal_pattern_count) {
            value.Add(mpos);
        }
    }


    return value;
}

void ShapeBuilder::DeepTrain(std::string output_db, std::string output_stats, std::string output_pattern_distr, std::string shape, size_t num_threads, size_t p, size_t l_end) {
    std::cout << "Deep train" << std::endl;

    ofstream ofs(output_stats, ios::out);

    std::cout << "Allocate builder... ";
    ShapeBuilder builder(shape);
    std::cout << "done" << std::endl;

    ProgressSingleton *progress = ProgressSingleton::GetInstance();
    builder.progress = progress;
    builder.pattern_size = p;

    size_t max_l = 0;

    bool auto_reg = l_end == 0;

    for (int i = 1; i <= l_end; i++) {
        std::cout << "limit: " << i << " lend: " << l_end << std::endl;
        builder.limit = i;

        builder.TrainByPatternTP(num_threads * 3);
//        std::cout << "Train by pattern single threaded" << std::endl;
//        builder.TrainByPatternTP(1);

//        for (auto entry : builder.map2) {
//            std::cout << entry.first << " . " << entry.second.ToString() << std::endl;
//        }

        if (max_l == 0) {
            size_t max_pot_mut = 0;

            ofstream distr_ofs(output_pattern_distr, ios::out);


            for (auto i = 0; i < 100; i++) {
                if (builder.pot_mut_counter[i] != 0)
                    max_pot_mut = i;
            }

            for (auto i = 0; i <= max_pot_mut; i++) {
                distr_ofs << i << "," << builder.pot_mut_counter[i] << std::endl;
            }
            distr_ofs.close();
            if (auto_reg)
                l_end = max_pot_mut;
        }


        size_t muts[10] = { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
        double sens[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


        builder.TestTP(100, 10000, muts, sens, 10, num_threads * 3);

        ofs << i;
        for (auto j = 0; j < 10; j++) {
            ofs << "," << sens[j];
        }
        ofs << "," << builder.pot_mut_counter[i];
        ofs << std::endl;

        builder.Save(output_db);
    }

    builder.Print();

    ofs.close();
}


void ShapeBuilder::ShapeSpaceExplorer(size_t take, size_t space, std::string outfile, double random_prob, size_t num_threads, size_t init_pattern_size, size_t init_limit, size_t upper_limit, size_t iterations) {

    std::unordered_set<size_t> checked_shapes;
    std::unordered_map<std::string,double> shape_history;

    checked_shapes.insert(0);

    ProgressSingleton *progress = ProgressSingleton::GetInstance();
//    ProgressBar *pb = new ProgressBar();

    std::string last_best = "";
    double last_mean = 0;
    size_t last_median = 0;

//    std::string outfile = output_folder + "/" + to_string(take) + "_" + to_string(space) + "_sumtest.csv";
    if (Utils::exists(outfile)) {
        std::cout << "exists.." << std::endl;
        auto str = csv::load_file(outfile.c_str());
        auto parser = csv::make_parser(str);
        for (auto&& row : parser ) {
            auto it = row.begin();
            std::string shape = (*it).to_string();
            double mean = (*(++++++++++it)).to_double();

            std::cout << shape << "," << mean << std::endl;

            size_t shape_num = ShapeBuilder::ToULL(shape);
            checked_shapes.insert(shape_num);
            shape_history.insert( { shape, mean } );

            if (mean < 1 && mean > last_mean) {
                last_best = shape;
                last_mean= mean;
            }
        }

    }

    std::cout << "start with " << last_best << std::endl;
    std::cout << "last_mean: " << last_mean << std::endl;


    ofstream ofs(outfile, std::ios_base::app);

    // Generate Random Shape to Start with
    size_t shape_length = take + space;
    size_t partial_k = take / 2;
    size_t partial_space = space / 2;
    size_t range = partial_k - 1 +  partial_space;

    srand(time(0));

    std::string init_shape = "";

    if (init_shape.length() == 0) {
        if (last_best.size() == 0) {
            std::cout << "random shape: ";
            init_shape = ShapeBuilder::RandomShape(take, space);
        } else {
            init_shape = last_best;
        }
    }

    std::string new_shape = init_shape;
    size_t shape_num = ShapeBuilder::ToULL(new_shape);

    std::cout << init_shape << std::endl;

    std::string operation = "INIT";
    size_t steps = iterations;


    size_t random_op_counter = 0;

    double sensitivity_ani[10];
    bool take_shape = false;

    while (steps-- > 0) {
        std::cout << "____________________________________________________________" << std::endl;
        take_shape = false;

        // reset sensitivity array
        memset(sensitivity_ani, 0, 10 * sizeof(double));

        if (random_op_counter > 10) {
            std::cout << "too many randoms.... reset to .. " << std::endl;

            double max = last_mean;
            double second_max = 0;
            std::string second_shape = "";

            for (auto p : shape_history) {
                if (p.second > second_max && p.second < max) {
                    second_max = p.second;
                    second_shape = p.first;
                }
            }

            init_shape = second_shape;
            last_mean = second_max;

            std::cout << " .." << init_shape << " - " << last_mean << std::endl;

            random_op_counter = 0;
        }

        while (checked_shapes.contains(shape_num)) {
            new_shape = ShapeBuilder::AlterShape(init_shape, 0.5, take, space, operation);
            std::cout << "newshape: " << new_shape << std::endl;
            shape_num = ShapeBuilder::ToULL(new_shape);
        }

        checked_shapes.insert(shape_num);

        // avoid too many randoms in a row
        if (operation.compare("RANDOM") == 0) {
            random_op_counter++;
        } else {
            random_op_counter = 0;
        }
        std::cout << "random counter: " << random_op_counter << std::endl;

        size_t sc = 0, tc = 0; //space and take counter to assure validity of shape
        for (auto ci = 0; ci < new_shape.size(); ci++) {
            char c = new_shape[ci];
            c == 'X' ? tc++ : sc++;
        }

        if (tc != take || sc != space){
            std::cout << tc << std::endl;
            std::cout << sc << std::endl;
            std::cout << take << std::endl;
            std::cout << space << std::endl;
            std::cout << init_shape << std::endl;
            std::cout << new_shape << std::endl;
            std::cout << operation << std::endl;
        }
//        init_shape = new_shape;

        assert(tc == take);
        assert(sc == space);



        ShapeBuilder builder(new_shape);
        builder.progress = progress;
//        builder.progress = pb;

        builder.pattern_size = init_pattern_size;
        builder.limit = init_limit;
//        builder.pattern_size = 12;
//        builder.limit = 8;

        builder.PrintSpecs();

        builder.TrainByPatternTP(num_threads * 3);
        std::cout << "done training patterns" << std::endl;

        size_t iterations = 10000;

//        double sensitivity90 = builder.Test(100, iterations, 10);

        size_t muts[10] = { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
        double sens[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//        Benchmark st_test("single_threaded");
//        st_test.start();
//        double sensitivity95 = builder.Test(100, iterations, 5);

        std::cout << "TestTP" << std::endl;
        builder.TestTP(100, 1000, muts, sens, 10, num_threads * 3);

        double mean = ShapeBuilder::Sum(sens, 10) / 10;
        std::cout << "mean: " << mean << std::endl;



        ofs << new_shape << "," << operation << "," << ShapeBuilder::CountGapNum(new_shape);// << std::endl;
        ofs << "," << DifferentSpacesLengths(builder.shape, builder.shape_size);
        ofs  << "," << LongestConsecutive(builder.shape, builder.shape_size);
        ofs << "," << mean;

        ShapeBuilder median_builder(new_shape);
        median_builder.pattern_size = 28;
        median_builder.limit = 1;

        std::cout << "Train BY PAtternTP" << std::endl;
        median_builder.TrainByPatternTP(num_threads * 3);

        size_t max_pot_mut = 0;
        for (auto i = 0; i < 100; i++) {
            if (median_builder.pot_mut_counter[i] != 0)
                max_pot_mut = i;
        }
        auto median = ShapeBuilder::median(median_builder.pot_mut_counter, max_pot_mut);

// TESTTEST
        std::cout << "TestTP one______________________________________________________" << std::endl;
        memset(sens, 0, 10 * sizeof(double));
        builder.TestTP(100, 10000, muts, sens, 10, num_threads * 3);

        std::cout << "TestTP two____________________________________________________________" << std::endl;
        auto result = ShapeResult(new_shape);
        ShapeFinder::Test(builder, 100, 10000, { 1, 2, 3, 4, 5, 6 ,7 ,8 ,9 }, num_threads * 3, result);
// TESTTEST

        if (median > last_median) {
//        if (mean > last_mean) {

            init_shape = new_shape;
            last_mean = mean;
            last_median = median;

            std::cout << "_________________________________________________retrain" << std::endl;
            builder.limit += upper_limit;

            builder.maximum_chunk = 100;
            builder.TrainByPatternTP(num_threads * 3);

            memset(sens, 0, 10 * sizeof(double));
            builder.TestTP(100, 10000, muts, sens, 10, num_threads * 3);

            mean = ShapeBuilder::Sum(sens, 10) / 10;


            shape_history.insert( { new_shape, mean } );

            take_shape = true;
        } else {
            steps++;
        }


        for (int i = 0; i < 10; i++) {
            ofs << "," << sens[i];
            if (i != 0)
                std::cout << ", ";
            std::cout << sens[i];
        }
        std::cout << std::endl;

        ofs << "," << mean;
        ofs << "," << (take_shape ? "TAKE" : "OMIT");
        ofs << "," << median;

        ofs << std::endl;
    }
}



void ShapeBuilder::ShapeSpaceExplorerSynonymous(size_t take, size_t space, std::string outfile, double random_prob, size_t num_threads, size_t init_pattern_size, size_t init_limit, size_t upper_limit, size_t iterations) {

    std::unordered_set<size_t> checked_shapes;
    std::unordered_map<std::string,double> shape_history;

    checked_shapes.insert(0);

    ProgressSingleton *progress = ProgressSingleton::GetInstance();
//    ProgressBar *pb = new ProgressBar();

    std::string last_best = "";
    double last_mean = 0;
    size_t last_median = 0;

//    std::string outfile = output_folder + "/" + to_string(take) + "_" + to_string(space) + "_sumtest.csv";
    if (Utils::exists(outfile)) {
        std::cout << "exists.." << std::endl;
        auto str = csv::load_file(outfile.c_str());
        auto parser = csv::make_parser(str);
        for (auto&& row : parser ) {
            auto it = row.begin();
            std::string shape = (*it).to_string();
            double mean = (*(++++++++++it)).to_double();

            std::cout << shape << "," << mean << std::endl;

            size_t shape_num = ShapeBuilder::ToULL(shape);
            checked_shapes.insert(shape_num);
            shape_history.insert( { shape, mean } );

            if (mean < 1 && mean > last_mean) {
                last_best = shape;
                last_mean= mean;
            }
        }

    }

    std::cout << "start with " << last_best << std::endl;
    std::cout << "last_mean: " << last_mean << std::endl;


    ofstream ofs(outfile, std::ios_base::app);

    // Generate Random Shape to Start with
    size_t shape_length = take + space;
    size_t partial_k = take / 2;
    size_t partial_space = space / 2;
    size_t range = partial_k - 1 +  partial_space;

    srand(time(0));

    std::string init_shape = "";

    if (init_shape.length() == 0) {
        if (last_best.size() == 0) {
            std::cout << "random shape: ";
            init_shape = ShapeBuilder::RandomShape(take, space);
        } else {
            init_shape = last_best;
        }
    }

    std::string new_shape = init_shape;
    size_t shape_num = ShapeBuilder::ToULL(new_shape);

    std::cout << init_shape << std::endl;

    std::string operation = "INIT";
    size_t steps = iterations;


    size_t random_op_counter = 0;

    double sensitivity_ani[10];
    bool take_shape = false;

    while (steps-- > 0) {
        std::cout << "____________________________________________________________" << std::endl;
        take_shape = false;

        // reset sensitivity array
        memset(sensitivity_ani, 0, 10 * sizeof(double));

        if (random_op_counter > 10) {
            std::cout << "too many randoms.... reset to .. " << std::endl;

            double max = last_mean;
            double second_max = 0;
            std::string second_shape = "";

            for (auto p : shape_history) {
                if (p.second > second_max && p.second < max) {
                    second_max = p.second;
                    second_shape = p.first;
                }
            }

            init_shape = second_shape;
            last_mean = second_max;

            std::cout << " .." << init_shape << " - " << last_mean << std::endl;

            random_op_counter = 0;
        }

        while (checked_shapes.contains(shape_num)) {
            new_shape = ShapeBuilder::AlterShape(init_shape, 0.5, take, space, operation);
            std::cout << "newshape: " << new_shape << std::endl;
            shape_num = ShapeBuilder::ToULL(new_shape);
        }

        checked_shapes.insert(shape_num);

        // avoid too many randoms in a row
        if (operation.compare("RANDOM") == 0) {
            random_op_counter++;
        } else {
            random_op_counter = 0;
        }
        std::cout << "random counter: " << random_op_counter << std::endl;

        size_t sc = 0, tc = 0; //space and take counter to assure validity of shape
        for (auto ci = 0; ci < new_shape.size(); ci++) {
            char c = new_shape[ci];
            c == 'X' ? tc++ : sc++;
        }

        if (tc != take || sc != space){
            std::cout << tc << std::endl;
            std::cout << sc << std::endl;
            std::cout << take << std::endl;
            std::cout << space << std::endl;
            std::cout << init_shape << std::endl;
            std::cout << new_shape << std::endl;
            std::cout << operation << std::endl;
        }
//        init_shape = new_shape;

        assert(tc == take);
        assert(sc == space);



        ShapeBuilder builder(new_shape);
        builder.progress = progress;
//        builder.progress = pb;

        builder.pattern_size = init_pattern_size;
        builder.limit = init_limit;
//        builder.pattern_size = 12;
//        builder.limit = 8;


        builder.TrainByPatternTP(num_threads * 3);
        std::cout << "done training patterns" << std::endl;

        size_t iterations = 10000;

//        double sensitivity90 = builder.Test(100, iterations, 10);

        size_t muts[10] = { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
        double sens[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//        Benchmark st_test("single_threaded");
//        st_test.start();
//        double sensitivity95 = builder.Test(100, iterations, 5);

        builder.TestTP(100, 1000, muts, sens, 10, num_threads * 3);
        double mean = ShapeBuilder::Sum(sens, 10) / 10;
        std::cout << "mean: " << mean << std::endl;



        ofs << new_shape << "," << operation << "," << ShapeBuilder::CountGapNum(new_shape);// << std::endl;
        ofs << "," << DifferentSpacesLengths(builder.shape, builder.shape_size);
        ofs  << "," << LongestConsecutive(builder.shape, builder.shape_size);
        ofs << "," << mean;

        ShapeBuilder median_builder(new_shape);
        median_builder.pattern_size = 28;
        median_builder.limit = 1;

        median_builder.TrainByPatternTP(num_threads * 3);

        size_t max_pot_mut = 0;
        for (auto i = 0; i < 100; i++) {
            if (median_builder.pot_mut_counter[i] != 0)
                max_pot_mut = i;
        }
        auto median = ShapeBuilder::median(median_builder.pot_mut_counter, max_pot_mut);

        if (median > last_median) {
//        if (mean > last_mean) {

            init_shape = new_shape;
            last_mean = mean;
            last_median = median;

            std::cout << "_________________________________________________retrain" << std::endl;
            builder.limit += upper_limit;

            builder.maximum_chunk = 100;
            builder.TrainByPatternTP(num_threads * 3);


            memset(sens, 0, 10 * sizeof(double));
            builder.TestTP(100, 10000, muts, sens, 10, num_threads * 3);
            mean = ShapeBuilder::Sum(sens, 10) / 10;

            shape_history.insert( { new_shape, mean } );

            take_shape = true;
        } else {
            steps++;
        }


        for (int i = 0; i < 10; i++) {
            ofs << "," << sens[i];
            std::cout << ", " << sens[i];
        }
        ofs << "," << mean;
        ofs << "," << (take_shape ? "TAKE" : "OMIT");
        ofs << "," << median;

        ofs << std::endl;


    }
}


void ShapeBuilder::TestTP(size_t read_length, size_t iterations, size_t* mutations, double* sensitivities, size_t mutations_size, size_t num_threads) {
    srand(time(0));

    size_t block_size = 500;
    size_t iteration_blocks = iterations / block_size;

    auto progress = ProgressSingleton::GetInstance();
    if (progress)
        progress->Reset(iteration_blocks * mutations_size);

    size_t prog_chunk = 1;

    {
        thread_pool pool(num_threads);
        for (auto it_block = 0; it_block < iteration_blocks; it_block++) {
            for (auto mi = 0; mi < mutations_size; mi++) {
                size_t m = mutations[mi];
                pool.enqueue_work([this, block_size, m, &sensitivities, mi]() {

                    TestTP(100, block_size, m, &sensitivities[mi]);
                });
            }
        }
    }

    for (size_t mi = 0; mi < mutations_size; mi++) {
        sensitivities[mi] /= iteration_blocks;
    }
}
void ShapeBuilder::TestTP(size_t read_length, size_t iterations, std::vector<size_t> mutations,
                          size_t num_threads, ShapeResult &result) {
    srand(time(0));

    size_t block_size = 500;
    size_t iteration_blocks = iterations / block_size;

    if (progress)
        progress->Reset(iteration_blocks * mutations.size());

    size_t prog_chunk = 1;

    {
        thread_pool pool(num_threads);
        for (auto it_block = 0; it_block < iteration_blocks; it_block++) {
            for (auto mi = 0; mi < mutations.size(); mi++) {
                size_t m = mutations[mi];
                pool.enqueue_work([this, block_size, m, read_length, &result]() {
                    TestTP(100, block_size, 1-((double)m/read_length), result);
                });
            }
        }
    }
}

std::mutex test_tp_mtx;
void ShapeBuilder::TestTP(size_t read_length, size_t iterations, size_t mutations, double* sensitivity) {
    bool *read = new bool[read_length];
    size_t *counts = new size_t[read_length];

    size_t hm_size = read_length - shape_size + 1;
    bool *hm_pattern = new bool[hm_size];

    // init?
    memset(read, 0, read_length * sizeof(bool));
    memset(hm_pattern, 0, hm_size * sizeof(bool));
    memset(counts, 0, read_length * sizeof(size_t));

    std::unordered_set<size_t> mutation_set;

    size_t fp = 0;
    size_t tp = 0;
    size_t fn = 0;

    double mutation_rate = (double) mutations / read_length;

    size_t total_hits = 0;

    size_t prog_counter = 0;
    while (iterations-- > 0) {
        double synonymous_probability = 0.8;
        GenerateReadRandom(read, read_length, mutation_rate, synonymous_probability);

        Test(read, read_length, hm_pattern, hm_size, mutation_set, tp, fp, fn, counts, total_hits);
    }
    if (fp != 0) {
        std::cout << ShapeBuilder::ReadToString(read, read_length) << std::endl;
        for (auto m : mutation_set) {
            std::cout << (int) m << ",";
        }
        std::cout << std::endl;
        std::string stop;
        std::cin >> stop;
    }
    assert (fp == 0);


    {
        std::lock_guard<std::mutex> sens_guard(test_tp_mtx);
        *sensitivity += ((double) tp / (tp + fn));
        if (*sensitivity == 0) {
            std::cout << mutations << " zero" << std::endl;
            exit(9);
        }
    }

    if (progress)
        progress->UpdateAdd(1);

    delete[] read;
    delete[] counts;
    delete[] hm_pattern;
}



void ShapeBuilder::TestTP(size_t read_length, size_t iterations, double ani, ShapeResult &result) {
    bool *read = new bool[read_length];
    size_t *counts = new size_t[read_length];

    size_t hm_size = read_length - shape_size + 1;
    bool *hm_pattern = new bool[hm_size];

    // init?
    memset(read, 0, read_length * sizeof(bool));
    memset(hm_pattern, 0, hm_size * sizeof(bool));
    memset(counts, 0, read_length * sizeof(size_t));

    std::unordered_set<size_t> mutation_set;

    size_t fp = 0;
    size_t tp = 0;
    size_t fn = 0;

    size_t prog_counter = 0;
    size_t total_hits = 0;
    size_t total_length = iterations * read_length;
    while (iterations-- > 0) {
        double synonymous_probability = 0.8;
        GenerateReadRandom(read, read_length, 1-ani, synonymous_probability);

        Test(read, read_length, hm_pattern, hm_size, mutation_set, tp, fp, fn, counts, total_hits);
    }
    assert (fp == 0);

    {
        std::lock_guard<std::mutex> sens_guard(test_tp_mtx);
//        auto& statistics = result.GetStatistic(to_string(ani));
//        statistics.fp += fp;
//        statistics.tp += tp;
//        statistics.fn += fn;
//        for (int i = 0; i < read_length; i++)
//            result.AddPosition(counts[i]);
    }

    if (progress)
        progress->UpdateAdd(1);

    delete[] read;
    delete[] counts;
    delete[] hm_pattern;
}

double ShapeBuilder::Test(size_t read_length, size_t iterations, size_t mutations) {
    srand(time(0));

    bool* read = new bool[read_length];
    size_t* counts = new size_t[read_length];

    size_t hm_size = read_length - shape_size + 1;
    bool* hm_pattern = new bool[hm_size];

    // init?
    memset(read, 0, read_length * sizeof(bool));
    memset(hm_pattern, 0, hm_size * sizeof(bool));
    memset(counts, 0, read_length * sizeof(size_t));

    std::unordered_set<size_t> mutation_set;

    size_t fp = 0;
    size_t tp = 0;
    size_t fn = 0;

    size_t total_hits = 0;

    size_t prog_counter = 0;
    while (iterations-- > 0) {
        GenerateReadRandom(read, read_length, mutations);
        Test(read, read_length, hm_pattern, hm_size, mutation_set, tp, fp, fn, counts, total_hits);
    }

    double sensitivity = ((double) tp / (tp + fn));
//    std::cout << "sensitivity: " << sensitivity << std::endl;

    delete[] read;
    delete[] counts;
    delete[] hm_pattern;

    return sensitivity;
}

void ShapeBuilder::TrainByPatternTP(size_t threads_num) {
    if (lower_limit == limit) {
        std::cerr << "shape has already been trained to a limit of " << limit << std::endl;
        std::cerr << "Please increase the limit" << std::endl;
    }

//    std::cout << "patternsize: " << pattern_size << std::endl;
    size_t pattern_perm = 1llu << pattern_size;
    size_t blocks = threads_num * 5;
    size_t chunk = pattern_perm / blocks;

    size_t max_chunk_size = maximum_chunk == 0 ? 2000 : maximum_chunk;

//    std::cout << "pattern_perm: " << pattern_perm << std::endl;
//    std::cout << "maximum_chunk: " << maximum_chunk << std::endl;
//    std::cout << "blocks: " << blocks << std::endl;
//    std::cout << "chunk: " << chunk << std::endl;
//    std::cout << "max_chunk_size: " << max_chunk_size << std::endl;


    while (chunk > max_chunk_size)
        chunk = pattern_perm / ++blocks;


//        std::cout << "maximum_chunk: " << maximum_chunk << std::endl;
//        std::cout << "blocks: " << blocks << std::endl;
//        std::cout << "chunk: " << chunk << std::endl;
//        std::cout << "max_chunk_size: " << max_chunk_size << std::endl;
//        std::cout << "pattern_perm: " << pattern_perm << std::endl;


    Mutations* tmp_storage = new Mutations[pattern_perm];

    if (progress) progress->Reset(blocks);

    size_t bucket_size = chunk + 1;

//    std::cout << "checkpoint 2" << std::endl;
//    uint32_t* pattern_sizes = new uint32_t[blocks];
//    uint32_t* patterns = new uint32_t[bucket_size * blocks];
//    memset(pattern_sizes, 0, blocks * sizeof(uint_t));
//    memset(patterns, 0, bucket_size * blocks * sizeof(uint32_t));

    uint64_t* pattern_sizes = new uint64_t[blocks];
    uint64_t* patterns = new uint64_t[bucket_size * blocks];
    memset(pattern_sizes, 0, blocks * sizeof(uint64_t));
    memset(patterns, 0, bucket_size * blocks * sizeof(uint64_t));

//    std::cout << "checkpoint 3" << std::endl;
    memset(pot_mut_counter, 0, 100 * sizeof(size_t));

    size_t offset = 0;

//    std::cout << "checkpoint 4" << std::endl;
    for (int p = 1; p < pattern_perm; p++) {
        size_t bucket = (p + offset) % blocks;
        size_t bucket_offset = bucket * bucket_size;
        patterns[bucket_offset + pattern_sizes[bucket]++] = p;
        if ((p % blocks) == 0) bucket_offset++;
    }

//    std::cout << "TrainByPattern:  " << threads_num << std::endl;
    {
        thread_pool pool(threads_num);
//        thread_pool pool2;
        for (auto i = 0; i < blocks; i++) {
//            std::cout << "start block : " << i << std::endl;
            pool.enqueue_work([this, i, patterns, bucket_size, pattern_sizes, &tmp_storage]() {
//            pool.enqueue_work([&]() {
                TrainByPatternMT(patterns + (i * bucket_size), pattern_sizes[i], tmp_storage);
            });
        }
    }

    // Insert pattern from tmp storage into map
    size_t pattern_count = 0;
//    for (auto p = 0; p < pattern_perm; p++) {

//    std::cout << "Insert: " << pattern_perm << std::endl;
    int empty_pattern_count = 0;
    for (uint64_t p = 0; p < pattern_perm; p++) {
        if (!tmp_storage[p].Empty()) {
//            std::cout << p << " " << tmp_storage[p].ToString() << std::endl;
            pattern_count++;
            uint64_t key = p << (64 - pattern_size);

            map2.insert( { key, tmp_storage[p] });

        } else {
            empty_pattern_count++;
        }
    }


//    for (auto entry : map2) {
//        std::cout << entry.first << " ? " << entry.second.ToString() << std::endl;
//    }
//
//    std::string stop;
//    std::cin >> stop;

    delete[] tmp_storage;
    delete[] pattern_sizes;
    delete[] patterns;

//    std::cout << "done training with pattern_count: " << pattern_count << std::endl;
    lower_limit = limit;
}

void ShapeBuilder::MoveSpace(std::string &shape, int dir, size_t pos, size_t take, size_t space, std::string &operation) {
    if (shape[pos] != '_') return;

    size_t shape_length = take + space;
    size_t partial_k = take / 2;
    size_t partial_space = space / 2;
    size_t range = partial_k - 1 +  partial_space;

//    std::cout << "starting_pos: " << pos << " dir: " << dir << " partial_k+1:" << partial_k << " : " << partial_k + partial_space << std::endl;

    const size_t SWAP_RIM = 0;
    const size_t SWAP = 1;
    const size_t MOVE = 2;

    double SWAP_RIM_P = 0.3;
    double SWAP_P = 0.3;
    double MOVE_P = 1 - SWAP_RIM_P - SWAP_P;

    double operation_probs[3];
    operation_probs[SWAP_RIM] = SWAP_RIM_P;
    operation_probs[SWAP] = SWAP_P;
    operation_probs[MOVE] = MOVE_P;

    size_t swap_with = 0;
    size_t swap_with_counter = 3;
    while (swap_with == 0 || swap_with > range) {
        if (swap_with_counter-- == 0) exit(3);
        if (dir < 0) {
            swap_with = pos - 1;
            while (swap_with != 0 && shape[swap_with] == '_') swap_with--;
        } else {
            swap_with = pos + 1;
            while (swap_with <= range && shape[swap_with] == '_') swap_with++;
        }
        dir *= -1;
    }

    size_t option = ShapeBuilder::RandomBiasedOption(operation_probs, 2);

    operation = "SWAP";

    switch (option) {
        case SWAP_RIM:
            pos = swap_with + dir;
            operation = "RIM_SWAP";
            break;
        case MOVE:
            size_t rim_pos = swap_with + dir;
            while (shape[rim_pos] == '_') rim_pos += dir;
            if (rim_pos > 0 && rim_pos < range) {
                pos = rim_pos - dir;
                operation = "MOVE";
            }
            break;
    };

    assert(pos > 0 && swap_with > 0 && pos <= range && swap_with <= range && pos != swap_with);

    shape[pos] = 'X';
    shape[shape.size() - 1 - pos] = 'X';
    shape[swap_with] = '_';
    shape[shape.size() - 1 - swap_with] = '_';
}


size_t ShapeBuilder::RandomBiasedOption(double *probabilities, size_t size) {
    double prob = RandomProb();
    int i;
    double sum = 0;
    for (i = 0; i < size; i++) {
        sum += probabilities[i];
        if (prob < sum) return i;
    }
    return i;
}

bool ShapeBuilder::RandomCoinFlip(double true_prob) {
    double prob = RandomProb();
    return prob < true_prob;
}

std::string ShapeBuilder::GetShapeStr() {
    return ShapeUtils::GetString(shape, shape_size);
}

void ShapeBuilder::PrintSpecs() {
    std::cout << "ShapeBuilder specs:\n" << std::endl;
    std::cout << "Shape length:\t" << shape_size << std::endl;
    std::cout << "Shape:       \t" << ShapeUtils::GetString(shape, shape_size) << std::endl;
    std::cout << "limit:       \t" << limit << std::endl;
    std::cout << "lower limit: \t" << lower_limit << std::endl;
    std::cout << "pattern size:\t" << pattern_size << std::endl;
}

void ShapeBuilder::Print() {
}

std::string
ShapeBuilder::AlterShape(std::string shape, double rand_prob, size_t &take, size_t &space, string &operation) {
    std::string new_shape = shape;

    int rand_dir = rand() % 2 ? 1 : -1;
    size_t partial_k = take / 2;
    size_t partial_space = space / 2;
    size_t range = partial_k - 1 +  partial_space;

    double total_rand_prob = rand_prob;

    bool random_read = BiasedCoinFlip(total_rand_prob);
    if (random_read) {
        std::cout << "random" << std::endl;
        operation = "RANDOM";
        new_shape = RandomShape(take, space);
    } else {
        int random = (rand() % range) + 1;
        while (shape[random] != '_')
            random = (rand() % range) + 1;
        ShapeBuilder::MoveSpace(new_shape, rand_dir, random, take, space, operation);
    }
    return new_shape;
}





// Tests
