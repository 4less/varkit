//
// Created by fritsche on 08/03/2021.
//

#pragma once

#include <cstdint>
#include <cstddef>
#include <iostream>
#include <err.h>
#include <sysexits.h>
#include <unordered_set>
#include <tgmath.h>
#include <fstream>
#include <sys/stat.h>
#include <assert.h>
#include <functional>


namespace MultiMap {
    static inline uint32_t reduce(uint32_t x, uint32_t N) {
        return x % N;
//        return ((uint64_t) x * (uint64_t) N) >> 32 ;
    }


    static inline uint32_t MurmurHash32(uint32_t k) {
        k *= 0xff51afd7;
        k ^= k >> 17;
        k *= 0xc4ceb9fe;
        k ^= k >> 17;
        return k;
    }

    static inline uint32_t Hash(uint64_t key, uint32_t capacity) {
        return reduce(MurmurHash32(key), capacity);
    }

    static inline uint64_t MurmurHash3(uint64_t k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccd;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53;
//    std::cout << k << std::endl;
        //k ^= k >> 33;
        return k;
    }

    static inline uint64_t Hash(uint64_t key, size_t value_bits, size_t max) {
        auto hash = MurmurHash3(key);
        if (hash > max)
            return hash >> value_bits;
        return hash;
    }

    static long GetFileSize(std::string filename)
    {
        struct stat stat_buf;
        int rc = stat(filename.c_str(), &stat_buf);
        return rc == 0 ? stat_buf.st_size : -1;
    }
    enum flag {
        UNSET = 0, UNIQUE = 1, MULTI_KEY_START = 2, MULTI_KEY = 3
    };
    static const std::string flag_str[] = { "Unset", "Unique", "Start of multi key", "Member of multi key" };

    struct Cell {
        uint64_t data = 0;

    public:
        inline bool empty() {
            return data == 0;
        }

        inline uint64_t key(const size_t value_bits, const size_t key_bits) {
            return (data >> value_bits) & ((1llu << key_bits) - 1);
        }

        inline uint64_t value(const size_t value_bits) {
            return (data & ((1llu << value_bits) - 1));
        }

//        inline bool unique(const size_t value_bits, const size_t key_bits) {
//            assert((data >> (key_bits + value_bits)) == 1llu || (data >> (key_bits + value_bits)) == 0llu);
//            return !(data >> (key_bits + value_bits));
//        }

        inline void inc_key(const size_t value_bits, const size_t key_bits) {
            set_key(value_bits, key(value_bits, key_bits)+1,key_bits);
        }

        inline std::string ToString(const size_t value_bits, const size_t key_bits) {
            std::string result = "";
            result += "key: " + std::to_string(key(value_bits, key_bits)) + " ";
            result += "value: " + std::to_string(value(value_bits)) + " ";
            result += "flag: " + flag_str[get_flag()];
            return result;
        }

//        inline void mark_nonunique(const size_t value_bits, const size_t key_bits) {
//            assert(key_bits + value_bits < 64);
//            data |= 1llu << (key_bits + value_bits);
//        }

        inline void set_flag(const size_t value_bits, const size_t key_bits, flag flag) {
            data &= ~(3llu << 62llu);
            data |= (uint64_t) flag << 62llu;
        }

        inline flag get_flag() {
            if ((data >> 62llu) > 3) {
                std::cout << "data: " << data << std::endl;
                std::cout << "flag: " << (data >> 62llu) << std::endl;
                exit(9);
            }
            return static_cast<flag>(data >> 62llu);
        }

        inline void paste(const size_t key, const size_t value, const size_t value_bits) {
            data = 0;
            data = key << value_bits;
            data |= value;
        }

        inline void set_key(const size_t value_bits, const size_t key, const size_t key_bits) {

            size_t value_before = value(value_bits);
//            std::cout << value_bits << ": " << ~(((1llu << key_bits) - 1) << value_bits) << std::endl;
            data &= ~(((1llu << key_bits) - 1) << value_bits);
            assert(value_before == value(value_bits));
            data |= key << value_bits;
            assert(value_before == value(value_bits));
        }
    };


    class MapIterator {

        size_t main_key;

        size_t bucket_begin;
        size_t bucket_end;
        size_t bucket_size;

        size_t iteration;
        bool hit = false;

    public:
        Cell* map;
        Cell* current_cell;
        size_t pack_size;
        size_t map_pos;
        size_t bucket_pos;
        const size_t value_bits_;
        const size_t key_bits_;
        void SetHit(size_t main_key, size_t bucket_begin, Cell* cell, size_t bucket_end, size_t pack_size, size_t bucket_pos) {
            this->bucket_pos = bucket_pos;
            this->main_key = main_key;
            this->bucket_begin = bucket_begin;
            this->bucket_end = bucket_end;
            this->bucket_size = bucket_end - bucket_begin;
            this->pack_size = pack_size;
            this->current_cell = cell;
            this->iteration = 0;
            this->hit = true;
            this->map_pos = bucket_begin + bucket_pos;
        }


        inline size_t PrevCell(size_t bucket_pos) {
            if (bucket_pos == 0)
                return bucket_size - 1;
            return bucket_pos - 1;
        }

        inline size_t NextCell(size_t bucket_pos) {
            if (bucket_pos == bucket_size - 1)
                return 0;
            return bucket_pos + 1;
        }

        bool Next(size_t &value) {
//            std::cout << "Next: " << bucket_pos << " bucket_begin: " << bucket_begin << " bucket_end: " << bucket_end << " bucket_size: " << bucket_size << " iteration: " << iteration << " pack_size: " << pack_size << std::endl;
            if (iteration == pack_size) {
//                std::cout << "iteration: " << iteration << " packsize: " << pack_size << std::endl;
                return false;
            }

//            std::cout << "current empty? " << current_cell->empty() << std::endl;


            value = current_cell->value(value_bits_);

            bucket_pos++;
            iteration++;

            if (bucket_begin + bucket_pos == bucket_end) {
                current_cell = map + bucket_begin;
                bucket_pos = 0;
            } else {
                current_cell = map + bucket_begin + bucket_pos;
            }

            return true;
        }

        bool Next(uint64_t &key, uint64_t &value, uint64_t &expected_value) {
            if (current_cell->key(value_bits_, key_bits_) != key || current_cell->value(value_bits_) != expected_value)
                return false;

            value = current_cell->value(value_bits_);

            bucket_pos++;
            if (bucket_begin + bucket_pos == bucket_end)
                current_cell = map + bucket_begin;
            else
                current_cell = map + bucket_begin + bucket_pos;
            iteration++;

            return true;
        }


        bool Find(size_t &value, size_t expected_value) {
            size_t start = bucket_pos;
            size_t end = bucket_pos + pack_size;



//            if (pack_size > 1) {
//                std::cout << "STARTFIND:  " << std::endl;
//                std::cout << "start: " << start << " end: " << end << std::endl;
//                std::cout << "pos: " << bucket_pos << " packsize: " << pack_size << std::endl;
////                exit(9);
//            }


//            std::cout << "find: " << expected_value << std::endl;
            auto iterations = 0;
            while (end >= start) {
                iterations++;
//                if (iterations > 1000) {
//                    exit(9);
//                }


                size_t mid = start + ((end - start)/2);

                auto cell = map + bucket_begin + (mid % bucket_size);

                if (cell->value(value_bits_) == expected_value) {
                    while ((map + bucket_begin + PrevCell(mid))->value(value_bits_) == expected_value)
                        mid = PrevCell(mid);
                    this->bucket_pos = mid % bucket_size;
                    value = cell->value(value_bits_);
                    return true;
                }
                if (cell->value(value_bits_) < expected_value) {
                    start = mid + 1;
                } else {
                    end = mid - 1;
                }
            }
            return false;
        }

        bool Found() {
            return hit;
        }

        void SetMiss() {
            hit = false;
        }

        MapIterator(size_t value_bits, size_t key_bits) : value_bits_(value_bits), key_bits_(key_bits) {}

    };

    class MultiMap {
        size_t capacity_;
        size_t size_;

        const size_t key_bits_;
        const size_t value_bits_;
        const size_t pack_size_bits_;

        const double load_factor_;

        const size_t offset_bits_;
        const size_t offset_size_;
        uint64_t *offset_;

        Cell *map_;

        // debug
        size_t *psl_;

    public:


        ~MultiMap() {
            // Destructor
            delete[] offset_;
            delete[] map_;
        }


        static inline uint64_t ComputeOffsetKey(const uint64_t key, size_t const key_bits) {
            return key >> key_bits;
        }

        static inline uint64_t ComputeInternalKey(const uint64_t key, size_t const key_bits) {
            return key & ((1llu << key_bits) - 1);
        }

        void FeedOffset(const uint64_t key) {
            offset_[ComputeOffsetKey(key, key_bits_)]++;
        }

        size_t Capacity() const {
            return capacity_;
        }

        size_t K() {
            return (key_bits_ + offset_bits_) / 2;
        }

        size_t OffsetBits() const {
            return offset_bits_;
        }

        size_t KeyBits() const {
            return key_bits_;
        }

        size_t ValueBits() const {
            return value_bits_;
        }

        Cell* Map() {
            return map_;
        }

        size_t* Offset() {
            return offset_;
        }

        size_t Size() const {
            return size_;
        }

        bool AllBucketsEmpty() {
            auto begin = map_;
            auto end = map_ + size_;

            for (auto cell = begin; cell != end; cell++) {
                if (!cell->empty())
                    return false;
            }
            return true;
        }

        size_t CountNonEmptyCells() {
            auto begin = map_;
            auto end = map_ + capacity_;

            size_t non_empty_count = 0;
            for (auto cell = begin; cell != end; cell++) {
//                std::cout << cell->data << std::endl;
                if (!cell->empty())
                    non_empty_count++;
            }
            return non_empty_count;
        }

        inline uint32_t psl(const uint64_t key, const uint64_t pos, const uint32_t bucket_size) {
            uint32_t hash = Hash(key, bucket_size);
            return pos + ((hash > pos) * bucket_size) - hash;
        }

        void CheckBuckets() {

            size_t min_size = std::numeric_limits<size_t>::max();
            size_t max_size = 0;
            double min_load = 1.0;
            double max_load = 0.0;

            std::unordered_set<size_t> key_set;
            size_t key_count = 0;
            for (int i = 0; i < offset_size_; i++) {
                key_set.clear();
                key_count = 0;


                auto bucket_index = offset_[i];
                auto bucket_size = offset_[i + 1] - bucket_index;

                if (bucket_size > max_size) max_size = bucket_size;
                if (bucket_size < min_size) min_size = bucket_size;

                for (int j = 0; j < bucket_size; j++) {
                    auto cell = map_[bucket_index + j];
                    if (!cell.empty()) {
                        key_set.insert(cell.key(value_bits_, key_bits_));
                        key_count++;
                    }
                }

                auto load_factor = (double) key_count / bucket_size;
                if (load_factor > max_load) max_load = load_factor;
                if (load_factor < min_load) min_load = load_factor;

                std::cerr << "bucket " << i;
                if (!key_count) {
                    std::cerr << " empty" << std::endl;
                } else {
                    std::cerr << " uniqueness " << ((double) key_set.size() / key_count) << std::endl;
                }
//                if (key_count != key_set.size()) {
//                    std::cerr << "wait a minute... " << std::endl;
//                    std::cerr << bucket_index << std::endl;
//                    std::cerr << key_set.size() << std::endl;
//                    std::cerr << key_count << std::endl;
//                    exit(10);
//                }
            }

            std::cout << std::endl;
            std::cout << "min_size: \t" << min_size << std::endl;
            std::cout << "max_size: \t" << max_size << std::endl;
            std::cout << "min_load: \t" << min_load << std::endl;
            std::cout << "max_load: \t" << max_load << std::endl;
        }

        void PrintMeta() {
            std::cout << std::endl;
            std::cout << "### Meta Data for Map:\n\n";
            std::cout << "Capacity:        \t" << capacity_ << std::endl;
            std::cout << "Size:            \t" << size_ << std::endl;
            std::cout << "Offset Size:     \t" << offset_size_ << std::endl;

            std::cout << "\nKey Bits:        \t" << key_bits_ << std::endl;
            std::cout << "Value Bits:      \t" << value_bits_ << std::endl;
            std::cout << "Offset Bits:     \t" << offset_bits_ << std::endl;
            std::cout << "Key Total Bits:  \t" << (offset_bits_ + key_bits_) << std::endl;
            std::cout << "\nLoad Factor:     \t" << load_factor_ << std::endl;
            std::cout << "Real Load Factor:\t" << ((double) size_ / capacity_) << std::endl;

            auto mem_map = (capacity_ * sizeof(Cell)) / (1024 * 1024);
            auto mem_offset = (offset_size_ * 8) / (1024 * 1024);

            std::cout << "\nMemory:" << std::endl;
            std::cout << "Map (MB):        \t" << mem_map << std::endl;
            std::cout << "Map Offset (MB): \t" << mem_offset << std::endl;
            std::cout << "Total (MB):      \t" << (mem_offset + mem_map) << std::endl;
            std::cout << std::endl;
        }

        /* #####################################################################
         * # Debug
         * #####################################################################
         */

        void printPSL() {
            for (int i = 0; i < 100; i++) {
                std::cout << i << ": " << psl_[i] << std::endl;
                if (!psl_[i]) break;
            }
        }

        static inline char GetBaseFromBit(uint64_t bits) {
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

        static std::string ToString(uint64_t kmer, size_t bits) {
            std::string result (bits/2, '_');
            for (int i = bits/2 - 1; i >= 0; i--) {
                result[i] = GetBaseFromBit(kmer & 3llu);
                kmer = kmer >> 2;
            }
            return result;
        }

        void print() {
            size_t offset_index = 0;
            size_t offset_value_1 = 0;
            size_t offset_value_2 = 0;
            size_t size = offset_[offset_index] - offset_[offset_index - 1];
            for (int i = 0; i < capacity_; i++) {
                while (offset_value_2 == i) {
                    std::cout << "--- " << offset_index
                              << " --------------------------------------------------------------- size: ";
                    offset_index++;
                    offset_value_1 = offset_value_2;
                    offset_value_2 = offset_[offset_index];
                    size = offset_value_2 - offset_value_1;
                    std::cout << size << std::endl;
                }
                auto c = map_[i];
                auto key = c.key(value_bits_, key_bits_);
                auto value = c.value(value_bits_);
                if (c.empty()) std::cout << (i - offset_value_1) << std::endl;
                else {
                    auto hash = Hash(key, (uint32_t) size);
                    auto pos = (i - offset_value_1);
                    std::cout << pos << ": " << ToString(offset_index - 1, offset_bits_) << " "
                              << ToString(key, key_bits_) << ": " << value << "\t\tHasholasho: " << hash
                              << "\t\tpsl: " << psl(key, pos, size) << std::endl;
                }

            }
        }

        using KeyExtractor = size_t(size_t value);

        size_t SelectPivot(size_t low, size_t high, size_t bucket_start, size_t bucket_size) {
            return high - 1;
        }

        long int Partition(long int low, long int high, size_t bucket_start, size_t bucket_size, size_t level) {
            auto get = [=](size_t bucket_pos) {
                return map_ + bucket_start + (bucket_pos % bucket_size);
            };

            auto key = [=](Cell* cell) {
                return cell->value(value_bits_);
            };


            auto pivot_pos = high;
            auto pivot_value = key(get(pivot_pos));

            // Pivot index
            auto i = low - 1ll;

            for (auto j = low; j < high; j++) {
                auto cell = get(j);

                if (key(cell) <= pivot_value)
                    Swap(cell, get(++i));
            }
            Swap(get(i + 1), get(pivot_pos));

            return i + 1;
        }

        bool IsSorted(size_t start, size_t pack_size, size_t bucket_start, size_t bucket_end) {
            auto cell = map_ + bucket_start + start;
            auto bucket_begin_cell { map_ + bucket_start };
            auto bucket_end_cell { map_ + bucket_end };

            uint64_t value = 0llu;
            for (auto i = 0; i < pack_size; i++) {
                if (value > cell->value(value_bits_)) {
                    std::cout << value << " > " << cell->value(value_bits_) << std::endl;
                    return false;
                }
                value = cell->value(value_bits_);
                cell = Next(cell, bucket_begin_cell, bucket_end_cell);
            }
            return true;
        }

        void Quicksort(long int low, long int high, size_t bucket_start, size_t bucket_size, size_t level) {

            if (low >= high || high < 0)
                return;

            auto p = Partition(low, high, bucket_start, bucket_size, level);

            // Left of pivot
            Quicksort(low, p - 1, bucket_start, bucket_size, level + 1);

            // Right of pivot
            Quicksort(p + 1, high, bucket_start, bucket_size, level + 1);
        }



        void GeneratePackSizes() {
            size_t offset_index = 0;
            size_t bucket_start = 0;
            size_t bucket_end = 0;
            size_t bucket_size = 0;

            size_t block = 0;

            auto cell_begin = map_;
            auto cell_end = map_;
            size_t progress_steps = capacity_/100;
            for (size_t i = 0llu; i < capacity_; i++) {
                if ((i % progress_steps) == 0)
                    std::cout << (i* 100/capacity_) << "%" << "\t\r" << std::flush;


                while (bucket_end <= i) {
                    offset_index++;
                    bucket_start = bucket_end;
                    bucket_end = offset_[offset_index];
                    bucket_size = bucket_end - bucket_start;

                    cell_begin = cell_end;
                    cell_end = map_ + bucket_end;
                }

                assert(bucket_start <= i && i < bucket_end);


                auto cell = map_ + i;
                if (cell->get_flag() != UNSET || cell->empty())
                    continue;

                // Debug
                if (i == 8) {
                    std::cout << "\n________________________Sus: " << i << std::endl;
                    std::cout << cell->ToString(value_bits_, key_bits_) << std::endl;
                    std::cout << "________________________" << std::endl;
                    print(i, 30);
                    std::cout << "_____Multikey start not followed by multikey member___________________" << std::endl;
                }

                auto prev = Prev(cell, cell_begin, cell_end);
                size_t pack_size = 0;

                assert(Next(Prev(cell, cell_begin, cell_end),cell_begin, cell_end) == cell);

                if (cell->get_flag() == UNSET && // Necessary condition
                    ((prev->get_flag() != UNSET || // Either previous has been set (meaning this key cant be a part of a multikey, otherwise it would have been set
                      prev->empty() || // Or the previous key is empty, then it will never be set
                      cell->key(value_bits_, key_bits_) != prev->key(value_bits_, key_bits_)))) // or this cell is unset, non empty and its content is different from the last cell
                {
                    // Key must either be unique or start of multikey

                    // Count size of
                    auto next = Next(cell, cell_begin, cell_end);

                    while (next->key(value_bits_, key_bits_) == cell->key(value_bits_, key_bits_) && next->get_flag() == UNSET && !next->empty()) {
                        pack_size++;
                        next = Next(next, cell_begin, cell_end);
                    }

                    if (pack_size > 0) {
                        auto bucket_pos = i - bucket_start;

                        Quicksort(bucket_pos, bucket_pos + pack_size, bucket_start, bucket_size, 0);


                        if (!(IsSorted(bucket_pos, pack_size, bucket_start, bucket_end))) {
                            print(i, pack_size + 2);

                            std::cout << "IsSorted(" << bucket_pos << ", " << pack_size << ", " << bucket_start << ", " << bucket_end << ");" << std::endl;
                        }

                        assert(IsSorted(bucket_pos, pack_size, bucket_start, bucket_end));


                        // Key is not uniquq
                        cell->set_flag(value_bits_, key_bits_, MULTI_KEY_START);
                        next = Next(cell, cell_begin, cell_end);
                        while (pack_size > 0 && next->get_flag() == UNSET) {

                            next->set_flag(value_bits_, key_bits_, MULTI_KEY);
                            if (pack_size >= (1llu << key_bits_)) {
                                std::cout << "packsize: " << pack_size << std::endl;
                                std::cout << "(1llu << key_bits_): " << (1llu << key_bits_) << std::endl;
                                exit(87);
                            }

                            if (!((Prev(next, cell_begin, cell_end)->get_flag() == MULTI_KEY_START && next->get_flag() == MULTI_KEY) || Prev(next, cell_begin, cell_end)->get_flag() != MULTI_KEY_START )) {
                                std::cout << "\n________________________" << std::endl;
                                print(i, 30);
                                std::cout << "\n_____Multikey start not followed by multikey member___________________" << std::endl;
                            }

                            assert(pack_size < (1llu << key_bits_));
                            assert((Prev(next, cell_begin, cell_end)->get_flag() == MULTI_KEY_START && next->get_flag() == MULTI_KEY) || Prev(next, cell_begin, cell_end)->get_flag() != MULTI_KEY_START );

                            next->set_key(value_bits_, pack_size, key_bits_);
                            pack_size--;
                            next = Next(next, cell_begin, cell_end);
                        }
                    } else {
                        // Key is unique
                        cell->set_flag(value_bits_, key_bits_, UNIQUE);
                    }
                }
                auto next = Next(cell, cell_begin, cell_end);

                if (!((Prev(next, cell_begin, cell_end)->get_flag() == MULTI_KEY_START && next->get_flag() == MULTI_KEY) || prev->get_flag() != MULTI_KEY_START )) {
                    std::cout << "\n________________________: " << i << std::endl;
                    std::cout << cell->ToString(value_bits_, key_bits_) << std::endl;
                    std::cout << "________________________" << std::endl;
                    print(i, 30);
                    std::cout << "_____Multikey start not followed by multikey member___________________" << std::endl;
                }

                assert((Prev(next, cell_begin, cell_end)->get_flag() == MULTI_KEY_START && next->get_flag() == MULTI_KEY) || prev->get_flag() != MULTI_KEY_START );
                // DEBUG FRO HERE THIS ASSERT SHOULD GO THROUGH
            }

            for (size_t i = 0llu; i < capacity_; i++) {
                if (map_[i].get_flag() == UNSET) {
                    if (!map_[i].empty()) {

                        size_t offset_key = 0;
                        while (offset_[offset_key] < i) offset_key++;
                        auto offset_start = map_ + offset_[offset_key - 1];
                        auto offset_end = map_ + offset_[offset_key];

                        std::cout << '\n' << offset_[offset_key] << " < " << i << " < " << offset_[offset_key + 1] << std::endl;

                        auto cell = Prev(map_ + i, offset_start, offset_end);

                        size_t backsteps = 0;
                        while (cell->key(value_bits_, key_bits_) == map_[i].key(value_bits_, key_bits_)) {
                            backsteps++;
                            std::cout << "back..." << std::endl;
                            cell = Prev(cell, offset_start, offset_end);
                        }
                        std::cout << "backsteps: " << backsteps << std::endl;

                        cell = Next(cell, offset_start, offset_end);

                        std::cout << "OG: \n" << (map_ + i)->ToString(value_bits_, key_bits_) << std::endl;

                        std::cout << cell->ToString(value_bits_, key_bits_) << std::endl;
                        std::cout << Next(cell, offset_start, offset_end)->ToString(value_bits_, key_bits_) << std::endl;

                        std::cout << std::string(50, '-') << std::endl;

                        std::cout << "!map_[i].empty()" << std::endl;
                        std::cout << i << std::endl;
                        print(i > 5 ? i-5 : 0, 20);
                        std::cout << "----" << std::endl;
                        std::cout << offset_[offset_key + 1] - 15 << std::endl;
                        print(offset_[offset_key + 1] - 15, 20);
                        std::string stop;
                        std::cin >> stop;
                    }
                    assert(map_[i].empty());
                }
            }
        }

        void printCell(size_t key, size_t value) {
            size_t offset_key = key >> key_bits_;
            size_t internal_key = key & ((1llu << key_bits_) - 1);

            size_t bucket_start { offset_[offset_key] };
            size_t bucket_size { offset_[offset_key + 1] - bucket_start };

            size_t bucket_pos = Hash(internal_key, bucket_size);

            std::cout << bucket_pos << ": " << offset_key << " " << ToString(offset_key, offset_bits_) << " "
                      << ToString(key, key_bits_) << " " << internal_key << " (" << key << "): " << value << "\t\tHash: " << bucket_pos
                      << "\t\tMOCK" << std::endl;
        }

        size_t GetIdeal(size_t key) {
            size_t offset_key = key >> key_bits_;
            size_t internal_key = key & ((1llu << key_bits_) - 1);

            size_t bucket_start { offset_[offset_key] };
            size_t bucket_size { offset_[offset_key + 1] - bucket_start };
            size_t bucket_pos = Hash(internal_key, bucket_size);
            std::cout << "IDEAL" << std::endl;
            std::cout << "bucket_start: " << bucket_start << std::endl;
            std::cout << "offset_key: " << offset_key << std::endl;
            std::cout << "offset_size: " << bucket_size << std::endl;
            return bucket_start + bucket_pos;
        }

        void print(size_t pos, size_t length) {
            size_t offset_index = 0;
            while (true) {
                if (offset_[offset_index] > pos) {
                    offset_index--;
                    break;
                }
                offset_index++;
            }
            size_t offset_value_1 = offset_[offset_index];
            size_t offset_value_2 = offset_[offset_index + 1];
            size_t size = offset_value_2 - offset_value_1;

            for (int i = pos; i < capacity_ && length > 0; i++, length--) {
                while (offset_value_2 == i) {
                    std::cout << "--- " << offset_index
                              << " --------------------------------------------------------------- size: ";
                    offset_index++;
                    offset_value_1 = offset_value_2;
                    offset_value_2 = offset_[offset_index + 1];
                    size = offset_value_2 - offset_value_1;
                    std::cout << size << std::endl;
                }
                auto c = map_[i];
                auto key = c.key(value_bits_, key_bits_);
                auto value = c.value(value_bits_);

                size_t full_key = (offset_index << key_bits_) | key;

                if (c.empty()) std::cout << (i - offset_value_1) << " empty" << std::endl;
                else {
                    auto hash = Hash(key, (uint32_t) size);
                    auto pos = (i - offset_value_1);
                    std::cout << pos << ": " << offset_index << " " << ToString(offset_index, offset_bits_) << " "
                              << ToString(key, key_bits_) << " " << key << " (" << full_key << ") -> " << value << "\t\tHash: " << hash
                              << "\t\tpsl: " << psl(key, pos, size) << "\t\tflag: " << flag_str[c.get_flag()]  << "\t\t\tdata: " << c.data << std::endl;
                }

            }
        }

        MultiMap(size_t key_bits, size_t value_bits, size_t offset_bits, const size_t *offset_sizes, double load_factor) :
                size_(0),
                key_bits_(key_bits),
                value_bits_(value_bits),
                pack_size_bits_(sizeof(Cell)*8 - value_bits - 1),
                offset_bits_(offset_bits),
                offset_size_(1llu << offset_bits_),
                load_factor_(load_factor) {

            offset_ = new uint64_t[offset_size_ + 1];

            size_t total = 0;
            for (auto i = 0; i < offset_size_; i++) {
                offset_[i] = total;
                total += ceil((1 / load_factor_) * (double) offset_sizes[i]);
            }

            offset_[offset_size_] = total;
            capacity_ = total;

            map_ = new Cell[capacity_];

            // analyse psls
            psl_ = new size_t[100];
            std::fill_n(psl_, 100, 0llu);
        }

        MultiMap(size_t capacity, size_t size, const size_t key_bits,
                               const size_t value_bits, const double load_factor, const size_t offset_bits,
                               const size_t offset_size, uint64_t *offset, Cell *map) : capacity_(capacity),
                                                                                        size_(size),
                                                                                        key_bits_(key_bits),
                                                                                        value_bits_(value_bits),
                                                                                        load_factor_(load_factor),
                                                                                        offset_bits_(offset_bits),
                                                                                        offset_size_(offset_size),
                                                                                        offset_(offset), map_(map),
                                                                                        pack_size_bits_(64 - value_bits - 1) {}


        inline Cell* Next(Cell* cell, Cell* bucket_begin, Cell* bucket_end) {
            cell++;
            if (cell == bucket_end)
                cell = bucket_begin;
            return cell;
        }

        inline Cell* Next(Cell* cell, Cell* bucket_begin, Cell* bucket_end, size_t& bucket_pos) {
            cell++;
            bucket_pos++;
            if (cell == bucket_end) {
                cell = bucket_begin;
                bucket_pos = 0;
            }
            assert(cell == (bucket_begin + bucket_pos));
            return cell;
        }

        inline Cell* Next(Cell* cell_begin, size_t &bucket_pos, size_t jump_by, size_t bucket_size) {
            bucket_pos = (bucket_pos + jump_by) % bucket_size;
            return (cell_begin + bucket_pos);
        }

        inline Cell* Prev(Cell* cell, Cell* bucket_begin, Cell* bucket_end) {
            if (cell == bucket_begin)
                cell = bucket_end - 1;
            else {
                cell--;
            }
            return cell;
        }

        void Find(uint64_t &key, MapIterator &iterator) {
            // Overlap 0 < key < pre_key_bits
            // Key that is being implicitly stored in offset_
            uint64_t offset_key = key >> key_bits_;

            // Key that is being explicitly stored in map_
            uint64_t main_key = key & ((1llu << key_bits_) - 1);

            // Bucket boundaries
            size_t bucket_start { offset_[offset_key] };
            size_t bucket_size { offset_[offset_key + 1] - bucket_start };
            size_t bucket_end { offset_[offset_key] + bucket_size };

            // new entry psl and resident entry psl
            size_t bucket_pos = Hash(main_key, bucket_size);
            size_t residual_psl = 0;

            // navigate to bucket in map_
            auto cell_begin { map_ + bucket_start };
            auto cell_end { map_ + bucket_start + bucket_size };
            auto cell = cell_begin + bucket_pos;

            assert(cell == (map_ + bucket_start + bucket_pos));


            // probing until found
            auto total_steps = 0;

            if (cell->get_flag() == MULTI_KEY) {
//                std::cout << "query multikey start" << std::endl;
                // If the starting pos is in the middle of a set of multikeys, the key cannot possible part of that.
                // Jump to the next multikey start (encoded in cells key)!
                cell = Next(cell_begin, bucket_pos, cell->key(value_bits_, key_bits_), bucket_size);
                assert(cell == (map_ + bucket_start + bucket_pos));
                assert(Prev(cell, cell_begin, cell_end)->get_flag() == MULTI_KEY);
            }
            //total_steps <= residual_psl && !cell->empty()
            while (true) {
//                if (cell != (map_ + bucket_start + bucket_pos)) {
//                    std::cout << "cell: " << cell->ToString(value_bits_, key_bits_) << std::endl;
//                    std::cout << "pos:  " << (map_ + bucket_start + bucket_pos)->ToString(value_bits_, key_bits_) << std::endl;
//                }
                assert(cell == (map_ + bucket_start + bucket_pos));
                assert(total_steps <= bucket_size);

                size_t jump_by;
                switch(cell->get_flag()) {
                    case UNSET:
                        // EMPTY BUCKET - SEARCH FAILED
                        iterator.bucket_pos = bucket_start + Hash(main_key, bucket_size);
                        iterator.SetMiss();
//                        std::cout << "empty bucket" << std::endl;
                        return;

                    case MULTI_KEY:
                        // MULTI_KEY -> IF MULTI_KEY START WAS NO HIT, USE THIS TO JUMP TO THE NEXT START/UNIQUE
                        jump_by = cell->key(value_bits_, key_bits_);
                        total_steps += jump_by;
//                        print(bucket_start + bucket_pos, 10);
//                        std::cout << "jump " << jump_by << " from " << cell->ToString(value_bits_, key_bits_) << std::endl;
//                        std::cout << "jump " << jump_by << std::endl;
                        assert(cell == (map_ + bucket_start + bucket_pos));
                        cell = Next(cell_begin, bucket_pos, jump_by, bucket_size);
                        assert(cell == (map_ + bucket_start + bucket_pos));
                        assert(Prev(cell, cell_begin, cell_end)->get_flag() == MULTI_KEY);
//                        std::cout << "to " << cell->ToString(value_bits_, key_bits_) << std::endl;
//                        print(bucket_start + bucket_pos, 10);
                        continue;

                    case UNIQUE:
                    case MULTI_KEY_START:
//                        std::cout << "multikey start" << std::endl;
                        // CELL EITHER UNIQUE OR START OF MULTIKEY
                        // CHECK IF PROBE SEQUENCE LENGTH (PSL) OF CURRENT CELL IS LOWER THAN TOTAL LOOKUPS
                        residual_psl = psl(cell->key(value_bits_, key_bits_), bucket_pos, bucket_size);
                        if (total_steps > residual_psl) {
//                            std::cout << "not found... total_steps: " << total_steps << " residual_psl: " << residual_psl << std::endl;
                            iterator.bucket_pos = bucket_start + Hash(main_key, bucket_size);
                            iterator.SetMiss();
                            return;
                        }
                        // CHECK IF QUERY KEY MATCHES CELL KEY
                        if (cell->key(value_bits_, key_bits_) == main_key) {
//                            std::cout << "found" << residual_psl << std::endl;
                            iterator.SetHit(main_key, bucket_start, cell, bucket_end, cell->get_flag() == UNIQUE ? 1 : Next(cell, cell_begin, cell_end)->key(value_bits_, key_bits_) + 1, bucket_pos);

//                            long int lookuppos = bucket_pos + bucket_start - 5;
//                            print(lookuppos >= 0 ? lookuppos : 0, 10);
//                            std::cout << "lookup :  " << key << std::endl;
//                            std::cout << "bucketpos :  " << bucket_pos << std::endl;
                            return;
                        }
                        // Next moves on to the next cell within the offset bucket.
                        // If next cell would leave the bucket, set it to the begin of the bucket
                        // and update bucket_pos to 0
                        assert(cell == (map_ + bucket_start + bucket_pos));
                        cell = Next(cell, cell_begin, cell_end, bucket_pos);
                        assert(cell == (map_ + bucket_start + bucket_pos));
                        total_steps++;
                        break;
                }
            }
        }

        inline void Swap(Cell* cell, uint64_t &main_key, uint64_t &value) {
            auto tmp_key = cell->key(value_bits_, key_bits_);
            auto tmp_value = cell->value(value_bits_);

            // update cell values with current key and value
            cell->paste(main_key, value, value_bits_);

            // continue with swapped key
            main_key = tmp_key;
            value = tmp_value;
        }

        inline void Swap(Cell* cell_a, Cell* cell_b) {
            auto data_b = cell_b->data;
            cell_b->data = cell_a->data;
            cell_a->data = data_b;
        }

        void Insert(uint64_t key, uint64_t value) {

            // Overlap 0 < key < pre_key_bits
            // Key that is being implicitly stored in offset_
            uint64_t offset_key = key >> key_bits_;

            // Key that is being explicitly stored in map_
            uint64_t main_key = key & ((1llu << key_bits_) - 1);
//
            // Bucket boundaries
            size_t bucket_start { offset_[offset_key] };
            size_t bucket_end { offset_[offset_key + 1] };
            size_t bucket_size { offset_[offset_key + 1] - bucket_start };

            // new entry psl and resident entry psl
            size_t bucket_pos = Hash(main_key, bucket_size);
            size_t insert_psl = 0;
            size_t residual_psl = 0;

            // navigate to bucket in map_
            auto cell_begin { map_ + bucket_start };
            auto cell_end { map_ + bucket_start + bucket_size };
            auto cell = cell_begin + bucket_pos;

            Cell* previous_cell = nullptr;
            size_t previous_key = -1;

            // probing until found
            auto total_steps = 0;

//            if (key == 90 and value == 11) {
//                std::cout << key << std::endl;
//                std::cout << value << std::endl;
//                std::cout << offset_key << std::endl;
//                std::cout << main_key << std::endl;
////                exit(9);
//            }

            while (!cell->empty()) {

                assert(total_steps < bucket_size);

                residual_psl = psl(cell->key(value_bits_, key_bits_), bucket_pos, bucket_size);

                auto swap = false;
//                if (value == 11) {
//                    std::cout << "\n___________________" << std::endl;
//                    std::cout << "Bucket pos: " << bucket_pos << std::endl;
//                    print(bucket_pos - 1, 5);
//                    std::cout << "psl: " << insert_psl << " > " << residual_psl << " then swap" << std::endl;
//                    std::cout << "prev cell: " << (previous_cell != nullptr) << " " << previous_key << " current key: " << main_key << std::endl;
//                    std::cout << "swap: " << ((insert_psl == residual_psl && previous_cell && previous_key == main_key)) << std::endl;
//                }

                if (insert_psl > residual_psl || (insert_psl == residual_psl && previous_cell && previous_key == main_key)) {;
                    swap = true;
                    Swap(cell, main_key, value);
                    insert_psl = residual_psl;
                }

//                if (swap) {
//                    print(bucket_pos, 5);
//                    std::cout << "mainkey: " << main_key << " value: " << value << std::endl;
//                }

                previous_cell = cell->empty() ? nullptr : cell;
                previous_key = cell->key(value_bits_, key_bits_);

                bucket_pos++;
                insert_psl++;
                cell++;

                // wrap over
                if (cell == cell_end) {
                    cell = cell_begin;
                    bucket_pos = 0;
                }
                total_steps++;


//                auto roll = 0;
//                while (cell->key(value_bits_, key_bits_) == main_key && !cell->empty()) {
////                    std::cout << "cellkey: " << cell->key(value_bits_, key_bits_) << " mainkey: " << main_key
////                              << " value: " << value << " value: " << cell->value(value_bits_) << std::endl;
//                    total_steps++;
//                    insert_psl++;
//                    bucket_pos++;
//                    if (bucket_pos == bucket_size) bucket_pos = 0;
//                    cell = cell_begin + bucket_pos;
//                    roll++;
//                }
//
//                if (roll) {
//                    std::cout << "roll: " << roll << std::endl;
//                    exit(8);
//                } else
//                    std::cout << "noll" << std::endl;


            }

            size_++;
            cell->paste(main_key, value, value_bits_);

//f
        }

        void Save(std::string file) {
            std::ofstream ofs(file, std::ofstream::binary);
            ofs.write((char *) &capacity_, sizeof(capacity_));
            ofs.write((char *) &size_, sizeof(size_));
            ofs.write((char *) &load_factor_, sizeof(load_factor_));

            ofs.write((char *) &key_bits_, sizeof(key_bits_));
            ofs.write((char *) &value_bits_, sizeof(value_bits_));

            ofs.write((char *) &offset_bits_, sizeof(offset_bits_));
            ofs.write((char *) &offset_size_, sizeof(offset_size_));

            ofs.write((char *) offset_, sizeof(offset_) * (offset_size_ + 1));
            std::cout << "save... entry byte: " << sizeof(Cell) << std::endl;
            ofs.write((char *) map_, sizeof(Cell) * capacity_);
            ofs.close();
        }


        static MultiMap *Load(std::string file) {
            size_t capacity, size, offset_bits, offset_size, key_bits, value_bits;
            double load_factor;

            std::ifstream ifs(file, std::ifstream::binary);

            ifs.read((char *) &capacity, 8);
            ifs.read((char *) &size, 8);
            ifs.read((char *) &load_factor, 8);

            ifs.read((char *) &key_bits, 8);
            ifs.read((char *) &value_bits, 8);

            ifs.read((char *) &offset_bits, 8);
            ifs.read((char *) &offset_size, 8);

            std::cout << "capacity:   \t" << capacity << std::endl;
            std::cout << "size:       \t" << size << std::endl;
            std::cout << "load_factor:\t" << load_factor << std::endl;
            std::cout << "key_bits:   \t" << key_bits << std::endl;
            std::cout << "value_bits: \t" << value_bits << std::endl;
            std::cout << "offset_bits:\t" << offset_bits << std::endl;
            std::cout << "offset_size:\t" << offset_size << std::endl;


            auto *offsets = new uint64_t[offset_size + 1];
            ifs.read((char *) offsets, sizeof(offsets) * (offset_size + 1));


            Cell *map = new Cell[capacity];
            ifs.read((char *) map, sizeof(Cell) * capacity);

            if (!ifs.get())
                errx(EX_IOERR, "Failed to load database file %s", file.c_str());

            auto *imap = new MultiMap(capacity, size, key_bits, value_bits, load_factor, offset_bits, offset_size,
                                           offsets, map);
            return imap;
        }
    };
};