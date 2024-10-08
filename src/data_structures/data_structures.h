//
// Created by fritsche on 07/06/22.
//

#ifndef VARKIT_DATA_STRUCTURES_H
#define VARKIT_DATA_STRUCTURES_H

#include <cstddef>
#include <string>
#include <fstream>
#include <iostream>
//#include "constants.h"
#include "version.h"
#include "Utils.h"

#include "ShapeUtils.h"

namespace ds {
    struct SNP {
    public:
        int read_pos_;
        int gene_pos_;

        friend bool operator<(const SNP& l, const SNP& r) {
            return l.gene_pos_ < r.gene_pos_;
        }
        friend bool operator>(const SNP& l, const SNP& r) {
            return l.gene_pos_ > r.gene_pos_;
        }
        friend bool operator >=(const SNP& l, const SNP& r) {
            return l.gene_pos_ >= r.gene_pos_;
        }
        friend bool operator <=(const SNP& l, const SNP& r) {
            return l.gene_pos_ <= r.gene_pos_;
        }
        friend bool operator ==(const SNP& l, const SNP& r) {
            return l.gene_pos_ == r.gene_pos_;
        }
        friend bool operator !=(const SNP& l, const SNP& r) {
            return l.gene_pos_ == r.gene_pos_;
        }

        SNP();
        SNP(int read_pos, int gene_pos);

        bool IsValid() const {
            return read_pos_ >= 0 && gene_pos_ >= -1;
        }

        std::string ToString() const {
            return "read_pos: " + std::to_string(read_pos_) + ", gene_pos: " + std::to_string(gene_pos_);
        }
    };

    struct HashSNP {
    public:
        size_t operator()(const SNP & snp) const {
            return std::hash<int>()(snp.gene_pos_);
        }
    };

    struct CompSNP {
    public:
        bool operator()(const SNP & snp1, const SNP & snp2) const {
            return snp1.gene_pos_ == snp2.gene_pos_;
        }
    };

    struct Mutations {
        static constexpr size_t m_size = 16;
        static constexpr size_t size_field = 15;
        static constexpr size_t max_size = m_size - 1;
        uint8_t mutations[m_size];
        const inline size_t Size() {
            assert(mutations[size_field] < size_field - 1);
            return mutations[size_field];
        };

//        static constexpr size_t m_size = 12;
//        uint8_t mutations[m_size];
//        uint16_t m_size = 0;
//        const size_t Size();

        Mutations();


        size_t Capacity();

        bool Empty();

        std::string ToString(int pos = 0);

        void Add(size_t mutation);
    };

    struct PatternMap {
        PatternMap(std::string path) {
            constexpr bool debug = false;
            /*
             * Created in void
             * ShapeBuilder::Save(std::string path);
             */
            std::ifstream ifs(path, std::ios::in);

            constexpr bool versioned = true;

            std::string line;


            if constexpr (versioned) {
                std::cout << "Load versioned.. " << std::endl;
                VersionMajorType major;
                VersionMinorType minor;
                VersionPatchType patch;
                VersionTweakType tweak;

                size_t pattern_sizeof;
                size_t mutations_sizeof;

                Version version;
                ifs.read((char *) &version, sizeof(version));
                std::cout << version.ToString() << std::endl;

                ifs.read((char *) &pattern_sizeof, sizeof(pattern_sizeof));
                ifs.read((char *) &mutations_sizeof, sizeof(mutations_sizeof));


                ifs.read((char *) &shape_size, sizeof(shape_size));
                ifs.read((char *) &pattern_size, sizeof(pattern_size));
                ifs.read((char *) &limit, sizeof(limit));
                ifs.read((char *) &lower_limit, sizeof(lower_limit));

                delete[] shape;
                shape = new bool[shape_size];
                ifs.read((char *) shape, shape_size * sizeof(bool));

                shape_str = ShapeUtils::GetString(shape, shape_size);

                auto map_array_length = 1llu << pattern_size;

                map_array = new ds::Mutations[map_array_length];
                while (ifs.peek() != EOF) {
                    uint64_t pattern;
                    ds::Mutations value;
                    ifs.read((char*) &pattern, 8);
                    ifs.read((char*) &value, sizeof(ds::Mutations));

                    uint32_t trunc_pattern = pattern >> (64llu - pattern_size);

                    map_array[trunc_pattern] = value;

                    if constexpr(debug) {
                        std::cout << trunc_pattern << " -> " << value.ToString() << std::endl;
                    }

                    if (value.Size() > 0 && value.mutations[0] == 0xFF) {
                        std::cout << value.Size() << " " << value.ToString() << std::endl;
                        exit(9);
                    }
                }
                ifs.close();
            } else {
                ifs.read((char *) &shape_size, sizeof(shape_size));
                ifs.read((char *) &pattern_size, sizeof(pattern_size));
                ifs.read((char *) &limit, sizeof(limit));
                ifs.read((char *) &lower_limit, sizeof(lower_limit));

                delete[] shape;
                shape = new bool[shape_size];
                ifs.read((char *) shape, shape_size * sizeof(bool));

                shape_str = ShapeUtils::GetString(shape, shape_size);

                auto map_array_length = 1llu << pattern_size;

                map_array = new ds::Mutations[map_array_length];
                while (ifs.peek() != EOF) {
                    uint64_t pattern;
                    ds::Mutations value;
                    ifs.read((char*) &pattern, 8);
                    ifs.read((char*) &value, sizeof(ds::Mutations));

                    uint32_t trunc_pattern = pattern >> (64llu - pattern_size);

                    map_array[trunc_pattern] = value;

                    if (value.Size()>0 && value.mutations[0] == 0xFF) {
                        std::cout << value.Size() << " " << value.ToString() << std::endl;
                        exit(9);
                    }
                }
                ifs.close();
            }

//            ifs.read((char *) &shape_size, sizeof(shape_size));
//            ifs.read((char *) &pattern_size, sizeof(pattern_size));
//            ifs.read((char *) &limit, sizeof(limit));
//            ifs.read((char *) &lower_limit, sizeof(lower_limit));
//
//            delete[] shape;
//            shape = new bool[shape_size];
//            ifs.read((char *) shape, shape_size * sizeof(bool));
//
//            shape_str = ShapeUtils::GetString(shape, shape_size);
//
//            auto map_array_length = 1llu << pattern_size;
//
//            map_array = new ds::Mutations[map_array_length];
//            while (ifs.peek() != EOF) {
//                uint64_t pattern;
//                ds::Mutations value;
//                ifs.read((char*) &pattern, 8);
//                ifs.read((char*) &value, sizeof(ds::Mutations));
//
//                uint32_t trunc_pattern = pattern >> (64llu - pattern_size);
//
//                map_array[trunc_pattern] = value;
//
//                if (value.Size()>0 && value.mutations[0] == 0xFF) {
//                    std::cout << value.Size() << " " << value.ToString() << std::endl;
//                    exit(9);
//                }
//            }
//            ifs.close();
        }

        ~PatternMap() {
            delete[] map_array;
            delete[] shape;
        }

        ds::Mutations* map_array = nullptr;
        size_t shape_size = 0;
        size_t pattern_size = 22;
        size_t limit = 16;
        size_t lower_limit = 0;
        std::string shape_str;
        bool* shape = nullptr;

        void Print() {
            std::cout << "shape_size: " << shape_size << std::endl;
            std::cout << "pattern_size: " << pattern_size << std::endl;
            std::cout << "limit: " << limit << std::endl;
            std::cout << "lower_limit: " << lower_limit << std::endl;
            std::cout << "shape_str: " << shape_str << std::endl;
            std::cout << "shape: " << (shape != nullptr) << std::endl;
        }
    };

    struct KmerHit {
        int taxid = -1;
        int geneid = -1;
        int genepos = -1;

        void Reset() {
            taxid = -1;
            geneid = -1;
            genepos = -1;
        }

        std::string ToString() {
            return (taxid > 0 ? "Hit \tTaxid: " : "Miss \tTaxid: ") + std::to_string(taxid) + "  Geneid: " + std::to_string(geneid) + "  Genepos: "
                   + std::to_string(genepos);
        }

        KmerHit() {}

        KmerHit(const KmerHit &other) {
            this->taxid = other.taxid;
            this->geneid = other.geneid;
            this->genepos = other.genepos;
        }

        KmerHit& operator=(const KmerHit& other) {
            this->taxid = other.taxid;
            this->geneid = other.geneid;
            this->genepos = other.genepos;
            return *this;
        }

        inline bool IsMiss() {
            return taxid == -1;
        }

        bool IsGeneSet() {
            return geneid != -1;
        }
    };
}


#endif //VARKIT_DATA_STRUCTURES_H
