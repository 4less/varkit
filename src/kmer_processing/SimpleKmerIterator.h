//
// Created by fritsche on 11/03/2021.
//

#pragma once

#include <iostream>
#include "FastxReader.h"
#include "KmerUtils.h"
#include "ShapeUtils.h"

class KmerIterator64 {
public:
    virtual void operator() (uint64_t& key) = 0;
};

class SimpleKmerIterator : KmerIterator64 {
    const bool * shape_;
    const uint32_t shape_size_;
    const uint32_t k_;
    FastxRecord record_;
    size_t pos_;

    uint64_t forward_key_;
    uint64_t reverse_key_;

    bool skip_N = false;

    inline void ExtractKeys(uint64_t &forward_key, uint64_t &reverse_key) {
        int offset = 0;
        forward_key = 0;
        reverse_key = 0;

        for (auto si = 0; si < shape_size_; si++) {
            if (skip_N && KmerUtils::BaseToInt(record_.sequence[pos_ + si]) > 3) {
                forward_key = UINT64_MAX;
                reverse_key = UINT64_MAX;
                return;
            }
            forward_key |= (!shape_[si]) * KmerUtils::BaseToInt(record_.sequence[pos_ + si], 0) << (2llu * (k_ + offset - si - 1llu ));
            reverse_key |= (!shape_[si]) * KmerUtils::BaseToIntC(record_.sequence[pos_ + si], 3) << (2llu * (si - offset));
            offset += shape_[si];
        }
    }

public:
    SimpleKmerIterator(uint32_t k, bool * shape, uint32_t shape_size) :
        shape_(shape),
        shape_size_(shape_size),
        k_(k) {}

    void SetRecord(FastxRecord &record) {
        assert(record.sequence.length() >= shape_size_);
        record_ = record;
        pos_ = 0;
    }

    inline bool HasNext() {
        return (pos_ < (record_.sequence.length() - shape_size_ + 1));
    }

    void inline operator () (uint64_t &key) {
        ExtractKeys(forward_key_, reverse_key_);

//        std::string forward_ex = ShapeUtils::ApplyShape(record_.sequence, pos_, shape_, shape_size_);
//        std::string reverse_ex = KmerUtils::reverseComplement(forward_ex);

//        if (forward_key_ != UINT64_MAX && (forward_key_ & (((1llu << (64 - k_*2)) - 1) << k_*2) ||
//            strcmp(KmerUtils::ToString(forward_key_, k_*2).c_str(), forward_ex.c_str()) ||
//            strcmp(KmerUtils::ToString(reverse_key_, k_*2).c_str(), reverse_ex.c_str())))
//        {
////            std::string forward_ex = ShapeUtils::ApplyShape(record_.sequence, pos_, shape_, shape_size_);
////            std::string reverse_ex = KmerUtils::reverseComplement(forward_ex);
//            std::cout << "sequence:        \t" << record_.sequence.substr(pos_, shape_size_) << std::endl;
//            std::cout << "shape:           \t" << ShapeUtils::GetString(shape_, shape_size_) << std::endl;
//            std::cout << "Applied forward: \t" << forward_ex << std::endl;
//            std::cout << "forward from key:\t" << KmerUtils::ToString(forward_key_, k_*2) << std::endl;
//            std::cout << "Applied reverse: \t" << reverse_ex << std::endl;
//            std::cout << "reverse from key:\t" << KmerUtils::ToString(reverse_key_, k_*2) << std::endl;
//        }

//        assert(forward_key_ == UINT64_MAX || !(forward_key_ & (((1llu << (64 - k_*2)) - 1) << k_*2)));
//        assert(forward_key_ == UINT64_MAX || strcmp(KmerUtils::ToString(forward_key_, k_*2).c_str(), forward_ex.c_str()) == 0);
//        assert(forward_key_ == UINT64_MAX || strcmp(KmerUtils::ToString(reverse_key_, k_*2).c_str(), reverse_ex.c_str()) == 0);

        key = (forward_key_ < reverse_key_) * forward_key_ + (forward_key_ >= reverse_key_) * reverse_key_;
        pos_++;
    }

    std::string GetString() {
        if ((pos_ - 2 + shape_size_) >= record_.sequence.length()) {
            std::cout << record_.to_string() << std::endl;
            std::cout << record_.sequence.length() << std::endl;
            std::cout << pos_ << std::endl;
        }
        assert(pos_ - 2 + shape_size_ < record_.sequence.length());
        return record_.sequence.substr(pos_ - 1 , shape_size_);
    }

    size_t GetPos() {
        return pos_ - 1;
    }
};


class RollingKmerIterator : KmerIterator64 {
    const bool * shape_;
    const uint32_t shape_size_;
    const uint32_t k_;
    FastxRecord record_;
    size_t pos_;

    size_t kmer_mask;

    size_t past_kmer_size;

    size_t* past_kmers_fwd;
    size_t* past_kmers_rev;
    size_t* past_kmer_masks;
    size_t* past_kmer_masks_rev;

    size_t skip_ = 0;
    bool skip_N_ = false;

    uint64_t forward_key_;
    uint64_t reverse_key_;

    static constexpr uint64_t N_REPLACEMENT = 0;
    static constexpr uint64_t N_REPLACEMENT_C = 0;

    void Init() {
        auto longest_gap = ShapeUtils::LongestGap(shape_, shape_size_);
        past_kmer_size = longest_gap + 1;

        past_kmer_masks = new size_t[past_kmer_size];
        past_kmer_masks_rev = new size_t[past_kmer_size];
        past_kmers_fwd = new size_t[past_kmer_size];
        past_kmers_rev = new size_t[past_kmer_size];

        std::fill_n(past_kmers_fwd, past_kmer_size, 0);
        std::fill_n(past_kmers_rev, past_kmer_size, 0);
        std::fill_n(past_kmer_masks, past_kmer_size, 0);
        std::fill_n(past_kmer_masks_rev, past_kmer_size, 0);

        size_t x_count = std::count(shape_, shape_ + shape_size_, false);

        kmer_mask = (1llu << (x_count * 2)) - 1;

        size_t rows = longest_gap + 2;
        size_t cols = x_count + 1;
        uint32_t matrix[rows][cols];
        std::fill_n(&(matrix[0][0]), cols * rows, 0);

        size_t last_row = rows - 1;
        for (int i = rows-1; i >= 0; i--) {
            size_t col = i == last_row ? 1 : 0;
            for (int j = 0; j < shape_size_; j++) {
                if (!shape_[j]) {
                    assert(col < cols);
                    uint32_t value = i + j;

                    if (i != last_row) {
                        if (matrix[last_row][col] != value)
                            value = 0;
                    }

                    matrix[i][col] = value;
                    col++;
                }
            }
        }

        for (auto i = 0; i < past_kmer_size; i++) {
            for (auto kpos = 0; kpos < cols; kpos++) {
                if (matrix[i][kpos]) {
                    past_kmer_masks[i] |= 3llu << (2*(cols - kpos - 1));
                    past_kmer_masks_rev[i] |= 3llu << (2*k_ - (2*(cols - kpos - 1)));
                }
            }
            past_kmer_masks_rev[i] >>= 2;
        }
    }

    inline void ExtractKeys(uint64_t &forward_key, uint64_t &reverse_key) {
        int offset = 0;
        forward_key = 0;
        reverse_key = 0;


        for (auto si = 0; si < shape_size_; si++) {
            // Not A C G T chars will be replaced by 0 and 3 respectively
            auto fwd_char = KmerUtils::BaseToInt(record_.sequence[pos_ + si]);
            auto rev_char = KmerUtils::BaseToIntC(record_.sequence[pos_ + si]);
            if (fwd_char > 3) {
                skip_ = shape_size_ - si;
//                forward_key = UINT64_MAX;
//                reverse_key = UINT64_MAX;
                fwd_char = N_REPLACEMENT;
                rev_char = N_REPLACEMENT_C;
//                rev_char = 3; // DEBUG ATTEMPT
            }
            forward_key |= (!shape_[si]) * fwd_char << (2llu * (k_ + offset- si - 1llu ));
            reverse_key |= (!shape_[si]) * rev_char << (2llu * (si - offset));
            offset += shape_[si];
        }

    }

    inline void RollNext(size_t &forward_key, size_t &reverse_key) {
        forward_key = 0;
        reverse_key = 0;

        for (auto i = 0; i < past_kmer_size; i++) {
            forward_key |= (past_kmers_fwd[i] & past_kmer_masks[i]);
            reverse_key |= (past_kmers_rev[i] & past_kmer_masks_rev[i]);
        }

        auto fwd_char = KmerUtils::BaseToInt(record_.sequence[pos_ + shape_size_ - 1], N_REPLACEMENT);

        skip_ = (fwd_char > 3) * shape_size_ + !(fwd_char > 3) * skip_;

        auto rev_char = KmerUtils::BaseToIntC(record_.sequence[pos_ + shape_size_ - 1], N_REPLACEMENT_C);

        forward_key |= fwd_char;
        reverse_key |= rev_char << ((k_ - 1) * 2);

        for (auto i = 0; i < past_kmer_size-1; i++) {
            past_kmers_fwd[i] = past_kmers_fwd[i+1];
            past_kmers_rev[i] = past_kmers_rev[i+1];
        }
        past_kmers_fwd[past_kmer_size-1] = (forward_key << 2) & kmer_mask;
        past_kmers_rev[past_kmer_size-1] = (reverse_key >> 2) & kmer_mask;
    }


public:
    RollingKmerIterator(uint32_t k, bool * shape, uint32_t shape_size) :
            shape_(shape),
            shape_size_(shape_size),
            k_(k) {
        Init();
    }

    ~RollingKmerIterator() {
        delete[] past_kmers_fwd;
        delete[] past_kmers_rev;
        delete[] past_kmer_masks_rev;
        delete[] past_kmer_masks;
    }

    size_t K() {
        return k_;
    }

    size_t ShapeSize() {
        return shape_size_;
    }

    void SetRecord(FastxRecord &record) {
        assert(record.sequence.length() >= shape_size_);
        record_ = record;
        pos_ = 0;
        skip_ = 0;
    }

    inline bool HasNext() {
        return (pos_ < (record_.sequence.length() - shape_size_ + 1));
    }

    void inline operator () (uint64_t &key) {
        if (pos_ < past_kmer_size) {
            ExtractKeys(forward_key_, reverse_key_);

            past_kmers_fwd[pos_] = (forward_key_ << 2) & kmer_mask;
            past_kmers_rev[pos_] = (reverse_key_ >> 2) & kmer_mask;

        } else {
            RollNext(forward_key_, reverse_key_);
        }

        key = skip_ && skip_N_ ? UINT64_MAX : (forward_key_ < reverse_key_) * forward_key_ + (forward_key_ >= reverse_key_) * reverse_key_;
        skip_ -= (skip_ > 0);
        pos_++;
    }

    std::string GetString() {
        if ((pos_ - 2 + shape_size_) >= record_.sequence.length()) {
            std::cout << record_.to_string() << std::endl;
            std::cout << record_.sequence.length() << std::endl;
            std::cout << pos_ << std::endl;
        }
        assert(pos_ - 2 + shape_size_ < record_.sequence.length());
        return record_.sequence.substr(pos_ - 1 , shape_size_);
    }

    size_t GetPos() {
        return pos_ - 1;
    }
};