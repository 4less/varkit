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
#include "KmerUtils.h"
#include <tgmath.h>
#include <fstream>
#include <sys/stat.h>
#include <assert.h>


namespace MiniMapUtils {
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
}


namespace MiniMap {

    template<typename Type>
    struct Cell {
        static_assert(std::is_same<Type, uint32_t>() or std::is_same<Type, uint64_t>());
        Type data = 0;

    public:
        inline bool empty() {
            return data == 0;
        }

        inline uint64_t key(const size_t value_bits) {
            return data >> value_bits;
        }

        inline uint64_t value(const size_t value_bits) {
            return (data & ((1llu << value_bits) - 1));
        }

        inline void paste(const size_t key, const size_t value, const size_t value_bits) {
            data = key << value_bits;
            data |= value;
        }
    };


    template<typename Type>
    class MiniMap {
        size_t capacity_;
        size_t size_;

        const size_t key_bits_;
        const size_t value_bits_;

        const double load_factor_;

        const size_t offset_bits_;
        const size_t offset_size_;
        uint64_t *offset_;

        Cell<Type> *map_;

        // debug
        size_t *psl_;

    public:


        ~MiniMap() {
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

        Cell<Type>* Map() {
            return map_;
        }

        size_t Size() const {
            return size_;
        }

        inline uint32_t psl(const uint64_t key, const uint64_t pos, const uint32_t bucket_size) {
            uint32_t hash = MiniMapUtils::Hash(key, bucket_size);
            return pos + ((hash > pos) * bucket_size) - hash;
        }

        void CheckBuckets() {

            size_t min_size = numeric_limits<size_t>::max();
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
                        key_set.insert(cell.key(value_bits_));
                        key_count++;
                    }
                }

                auto load_factor = (double) key_count / bucket_size;
                if (load_factor > max_load) max_load = load_factor;
                if (load_factor < min_load) min_load = load_factor;

                if (key_count != key_set.size()) {
                    std::cerr << "wait a minute... " << std::endl;
                    std::cerr << bucket_index << std::endl;
                    std::cerr << key_set.size() << std::endl;
                    std::cerr << key_count << std::endl;
                    exit(10);
                }
            }

            std::cout << " " << std::endl;
            std::cout << "min_size: \t" << min_size << std::endl;
            std::cout << "max_size: \t" << max_size << std::endl;
            std::cout << "min_load: \t" << min_load << std::endl;
            std::cout << "max_load: \t" << max_load << std::endl;
        }

        void PrintMeta() {
            std::cout << " " << std::endl;
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

            auto mem_map = (capacity_ * sizeof(Type)) / (1024 * 1024);
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
                auto key = c.key(value_bits_);
                auto value = c.value(value_bits_);

                if (c.empty(value_bits_)) std::cout << (i - offset_value_1) << std::endl;
                else {
                    auto hash = MiniMapUtils::Hash(key, (uint32_t) size);
                    auto pos = (i - offset_value_1);
                    std::cout << pos << ": " << KmerUtils::ToString(offset_index - 1, offset_bits_) << " "
                              << KmerUtils::ToString(key, key_bits_) << ": " << value << "\t\tHasholasho: " << hash
                              << "\t\tpsl: " << psl(key, pos, size) << std::endl;
                }

            }
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
                auto key = c.key(value_bits_);
                auto value = c.value(value_bits_);

                if (c.empty(value_bits_)) std::cout << (i - offset_value_1) << std::endl;
                else {
                    auto hash = MiniMapUtils::Hash(key, (uint32_t) size);
                    auto pos = (i - offset_value_1);
                    std::cout << pos << ": " << KmerUtils::ToString(offset_index - 1, offset_bits_) << " " << key << " "
                              << KmerUtils::ToString(key, key_bits_) << ": " << value << "\t\tHash: " << hash
                              << "\t\tpsl: " << psl(key, pos, size) << std::endl;
                }

            }
        }


        size_t Size2() const {
            size_t total = 0;
            for (int i = 0; i < capacity_; i++) {
                if (!map_[i].empty(value_bits_)) total++;
            }
            return total;
        }

        MiniMap(size_t key_bits, size_t value_bits, size_t offset_bits, const size_t *offset_sizes,
                      double load_factor) :
                size_(0),
                key_bits_(key_bits),
                value_bits_(value_bits),
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

            map_ = new Cell<Type>[capacity_];

            // analyse psls
            psl_ = new size_t[100];
            std::fill_n(psl_, 100, 0llu);
        }

        MiniMap(size_t capacity, size_t size, const size_t key_bits,
                               const size_t value_bits, const double load_factor, const size_t offset_bits,
                               const size_t offset_size, uint64_t *offset, Cell<Type> *map) : capacity_(capacity),
                                                                                        size_(size),
                                                                                        key_bits_(key_bits),
                                                                                        value_bits_(value_bits),
                                                                                        load_factor_(load_factor),
                                                                                        offset_bits_(offset_bits),
                                                                                        offset_size_(offset_size),
                                                                                        offset_(offset), map_(map) {}


        Cell<Type> *Find(uint64_t &key) {
            uint64_t pre_key = key >> key_bits_;
            uint64_t post_key = key & ((1llu << key_bits_) - 1);

            assert(pre_key < (1llu << offset_bits_));
            assert(key < (1llu << (offset_bits_ + key_bits_)));

            size_t bucket_index = offset_[pre_key];
            uint32_t bucket_size = offset_[pre_key + 1] - bucket_index;
            size_t cur_pos = MiniMapUtils::Hash(post_key, bucket_size);
            Cell<Type> *cell = map_ + bucket_index + cur_pos;
            auto probes = 0;

            for (int i = 0; i < 1; i++) {
                if (cell->key(value_bits_) == post_key && !cell->empty()) {
                    return cell;
                }

                ++cell;
                if (++cur_pos == bucket_size) {
                    cell = map_ + bucket_index;
                    cur_pos = 0;
                }
                ++probes;
            }

            size_t c_psl = psl(cell->key(value_bits_), cur_pos, bucket_size);

            while (probes <= c_psl && !cell->empty()) {
                if (cell->key(value_bits_) == post_key) {
//            psl_[c_psl]++;
                    return cell;
                }

                ++cell;
                ++probes;
                if (++cur_pos == bucket_size) {
                    cell = map_ + bucket_index;
                    cur_pos = 0;
                }
                c_psl = psl(cell->key(value_bits_), cur_pos, bucket_size);
            }

            return nullptr;
        }

        void Insert(uint64_t key, uint64_t value) {
            // Overlap 0 < key < pre_key_bits
            // Key that is being implicitly stored in offset_
            uint64_t offset_key = key >> key_bits_;

            // Key that is being explicitly stored in map_
            uint64_t main_key = key & ((1llu << key_bits_) - 1);

            size_t bucket_index = offset_[offset_key];
            uint32_t bucket_size = offset_[offset_key + 1] - bucket_index;

            // new entry psl and resident entry psl
            size_t cur_pos = MiniMapUtils::Hash(main_key, bucket_size);
            size_t ins_psl = 0;
            size_t res_psl = 0;

            // navigate to bucket in map_
            Cell<Type> *cell = map_ + bucket_index + cur_pos;

            // probing until found

            auto stuck_counter = 0;

            while (!cell->empty()) {
                stuck_counter++;

                if (cell->key(value_bits_) == main_key) {
                    size_--;
                    break;
                }

                if (stuck_counter > bucket_size) {
                    std::cerr << "Bucket size too small: " << bucket_size << std::endl;
                    exit(9);
                }

                res_psl = psl(cell->key(value_bits_), cur_pos, bucket_size);

                if (ins_psl > res_psl) {
                    // find
                    // Swap entries by copying cell data into tmp key and value
                    auto tmp_key = cell->key(value_bits_);
                    auto tmp_value = cell->value(value_bits_);

                    // update cell values with current key and value
                    cell->paste(main_key, value, value_bits_);

                    // continue with swapped key
                    main_key = tmp_key;
                    value = tmp_value;

                    ins_psl = res_psl;
                }

                cur_pos++;
                ins_psl++;
                cell++;

                if (cur_pos == bucket_size) {
                    cell = map_ + bucket_index;
                    cur_pos = 0;
                }
            }

            size_++;

            cell->paste(main_key, value, value_bits_);
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
            std::cout << "save... entry byte: " << sizeof(Cell<Type>) << std::endl;
            std::cout << "save... entry byte: " << sizeof(Type) << std::endl;
            ofs.write((char *) map_, sizeof(Cell<Type>) * capacity_);
            ofs.close();
        }


        static MiniMap<Type> *Load(std::string file) {
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


            uint64_t *offsets = new uint64_t[offset_size + 1];
            ifs.read((char *) offsets, sizeof(offsets) * (offset_size + 1));

//            for (int i = 0; i < offset_size + 1; i++) {
//                std::cout << i << ": " << offsets[i] << std::endl;
//            }


            Cell<Type> *map = new Cell<Type>[capacity];
            ifs.read((char *) map, sizeof(Cell<Type>) * capacity);

//            for (int i = 0; i < capacity; i++) {
//                std::cout << map[i].data << std::endl;
//            }

            if (!std::char_traits<char>::eof() == ifs.get())
                errx(EX_IOERR, "Failed to load database file %s", file.c_str());

            auto *imap = new MiniMap<Type>(capacity, size, key_bits, value_bits, load_factor, offset_bits, offset_size,
                                           offsets, map);
            return imap;
        }

        bool IsSane() {
            return (this->key_bits_ + this->value_bits_) <= sizeof(Type) * 8;
        }
    };
};