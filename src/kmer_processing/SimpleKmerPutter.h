//
// Created by fritsche on 12/03/2021.
//

#pragma once

#include <cstdint>
#include <TaxonomyNew.h>
#include "compact_map.h"
#include "minimap.h"
#include "multimap.h"
//#include "multimap.h"

#define DEBUG 0

class MarkerKmerProcessor {
public:
    virtual void operator () (uint64_t &key, uint64_t &internal_taxid, uint64_t &marker_gene_id, uint64_t &marker_gene_position) = 0;
};

class KmerProcessor {
public:
    virtual void operator () (uint64_t &key, uint64_t &internal_taxid) = 0;
};



class CountKmerProcessor : MarkerKmerProcessor {
    using InsertKeyExists = std::function<size_t (const size_t, const size_t)>;

    Taxonomy::IntTaxonomy *taxonomy_;
    IndexedMap *map_;

    const size_t internal_taxid_bits_;
    const size_t marker_gene_id_bits_;
    const size_t marker_gene_position_bits_;

    const size_t max_internal_taxid_;
    const size_t max_marker_gene_id_;
    const size_t max_marker_gene_position_;

public:
    CountKmerProcessor(IndexedMap *map, Taxonomy::IntTaxonomy *taxonomy, size_t internal_taxid_bits, size_t marker_gene_id_bits, size_t marker_gene_position_bits) :
            internal_taxid_bits_(internal_taxid_bits),
            marker_gene_id_bits_(marker_gene_id_bits),
            marker_gene_position_bits_(marker_gene_position_bits),
            max_internal_taxid_((1llu << internal_taxid_bits) - 1),
            max_marker_gene_id_((1llu << marker_gene_id_bits) - 1),
            max_marker_gene_position_((1llu << marker_gene_position_bits) - 1),
            map_(map),
            taxonomy_(taxonomy) {
        if ((marker_gene_position_bits_ + marker_gene_id_bits_ + internal_taxid_bits_) > map_->value_bits_) {
            errx(EX_SOFTWARE,
                 "Bit sizes for internal taxid (%llu), marker gene id (%llu) and marker gene position (%llu) must add up (or be less) to the value size declared in the map (%llu).",
                 (unsigned long long) internal_taxid_bits_,
                 (unsigned long long) marker_gene_id_bits_,
                 (unsigned long long) marker_gene_position_bits_,
                 (unsigned long long) map_->value_bits_);

        }
    }

    void operator () (uint64_t &key, uint64_t &internal_taxid, uint64_t &marker_gene_id, uint64_t &marker_gene_position) override {

        uint64_t value = 1;
        uint64_t main_key = key & ((1llu << map_->KeyBits()) - 1);

        // comment in if not multimap
        auto found = map_->Find(key);

        uint64_t taxid = internal_taxid;

        if (found) {
            auto old_internal_taxid = found->value(map_->value_bits_) >> (marker_gene_id_bits_ + marker_gene_position_bits_);
            value = found->value(map_->value_bits_) & ((1lu << (marker_gene_id_bits_ + marker_gene_position_bits_)) - 1);
            value++;

            if (old_internal_taxid != internal_taxid) {
                taxid = taxonomy_->LCA(old_internal_taxid, internal_taxid);
            }

            value |= taxid << (marker_gene_id_bits_ + marker_gene_position_bits_);

            found->paste(key, value, map_->value_bits_);
        } else {
            value |= internal_taxid << (marker_gene_id_bits_ + marker_gene_position_bits_);
            assert(!taxonomy_->Get(internal_taxid).rank.empty());
            map_->Insert(key, value);
        }
    };
};



class SimpleMarkerKmerProcessor : MarkerKmerProcessor {
    using InsertKeyExists = std::function<size_t (const size_t, const size_t)>;

    Taxonomy::IntTaxonomy *taxonomy_;
    IndexedMap *map_;

    const size_t internal_taxid_bits_;
    const size_t marker_gene_id_bits_;
    const size_t marker_gene_position_bits_;

    const size_t max_internal_taxid_;
    const size_t max_marker_gene_id_;
    const size_t max_marker_gene_position_;

    const size_t taxid_mask;

    const bool update_;

    InsertKeyExists function;

public:
    SimpleMarkerKmerProcessor(IndexedMap *map, Taxonomy::IntTaxonomy *taxonomy, size_t internal_taxid_bits, size_t marker_gene_id_bits, size_t marker_gene_position_bits, bool update=false) :
            internal_taxid_bits_(internal_taxid_bits),
            marker_gene_id_bits_(marker_gene_id_bits),
            marker_gene_position_bits_(marker_gene_position_bits),
            max_internal_taxid_((1llu << internal_taxid_bits) - 1),
            max_marker_gene_id_((1llu << marker_gene_id_bits) - 1),
            max_marker_gene_position_((1llu << marker_gene_position_bits) - 1),
            map_(map),
            update_(update),
            taxid_mask(((1llu << internal_taxid_bits) - 1) << (marker_gene_id_bits + marker_gene_position_bits)),
            taxonomy_(taxonomy) {
        if ((marker_gene_position_bits_ + marker_gene_id_bits_ + internal_taxid_bits_) > map_->value_bits_) {
            errx(EX_SOFTWARE,
                 "Bit sizes for internal taxid (%llu), marker gene id (%llu) and marker gene position (%llu) must add up (or be less) to the value size declared in the map (%llu).",
                 (unsigned long long) internal_taxid_bits_,
                 (unsigned long long) marker_gene_id_bits_,
                 (unsigned long long) marker_gene_position_bits_,
                 (unsigned long long) map_->value_bits_);

        }
        function = [this] (size_t v1, size_t v2) {
            v1 = v1 >> (marker_gene_id_bits_ + marker_gene_position_bits_);
            v2 = v2 >> (marker_gene_id_bits_ + marker_gene_position_bits_);
            return taxonomy_->LCA(v1,v2) << (marker_gene_id_bits_ + marker_gene_position_bits_);
        };
    }

    size_t count_finds = 0;

    void operator () (uint64_t &key, uint64_t &internal_taxid, uint64_t &marker_gene_id, uint64_t &marker_gene_position) override {
#if DEBUG
        if (internal_taxid > max_internal_taxid_) {
            errx(EX_SOFTWARE,
                 "Internal tax id (%llu) is bigger than max internal taxid (%llu).",
                 (unsigned long long) internal_taxid,
                 (unsigned long long) max_internal_taxid_);
        }
        if (marker_gene_id > max_marker_gene_id_) {
            errx(EX_SOFTWARE,
                 "Marker gene id (%llu) is bigger than max marker gene id (%llu).",
                 (unsigned long long) marker_gene_id,
                 (unsigned long long) max_marker_gene_id_);
        }
        if (marker_gene_position > max_marker_gene_position_) {
            errx(EX_SOFTWARE,
                 "Marker gene position (%llu) is bigger than max marker gene position (%llu).",
                 (unsigned long long) marker_gene_position,
                 (unsigned long long) max_marker_gene_position_);
        }
#endif
        uint64_t value = 0;
        uint64_t bucket_key = key >> map_->key_bits_;
        uint64_t main_key = key & ((1llu << map_->KeyBits()) - 1);

        // comment in if not multimap
        auto found = map_->Find(key);

        uint64_t taxid = internal_taxid;

        if (found) {

            auto old_internal_taxid = found->value(map_->value_bits_) >> (marker_gene_id_bits_ + marker_gene_position_bits_);
            auto old_gene_id = (found->value(map_->value_bits_) >> marker_gene_position_bits_) & ((1lu << marker_gene_id_bits_) - 1);
            auto old_gene_pos = found->value(map_->value_bits_) & ((1lu << marker_gene_position_bits_) - 1);

            if (old_internal_taxid != internal_taxid) {
                taxid = taxonomy_->LCA(old_internal_taxid, internal_taxid);
//                taxid = 1;
            } else if (update_) {
                return;
            }
            assert(!taxonomy_->Get(taxid).rank.empty());

//            if (old_gene_id == marker_gene_id) {
//                printf("it is happening.. tid %llu  oldtid %llu  mgid %llu  mgpos %llu  oldpos %llu\n", taxid, old_internal_taxid, marker_gene_id, marker_gene_position, old_gene_pos);
//            }

            if (old_gene_id == marker_gene_id && old_gene_pos == marker_gene_position) {
//                printf("it is happening.. tid %llu  oldtid %llu  mgid %llu  mgpos %llu\n", taxid, old_internal_taxid, marker_gene_id, marker_gene_position);

                value = found->value(map_->value_bits_) & ~taxid_mask;
//                std::cout << "before: " << found->value(map_->value_bits_) << std::endl;
//                std::cout << "mask: " << ~taxid_mask << std::endl;
//                std::cout << "after: " << value << std::endl;
            }
//            auto flag = false;
//            if (key == 31241926069373){
//                printf("it is happening.. tid %llu  oldtirrrd %llu  oldmgid %llu  mgid %llu  mgpos %llu  oldpos %llu\n", taxid, old_internal_taxid, old_gene_id, marker_gene_id, marker_gene_position, old_gene_pos);
//                flag = true;
//            }

            value |= taxid << (marker_gene_id_bits_ + marker_gene_position_bits_);
            found->paste(key, value, map_->value_bits_);

//            if (flag){
//                std::cout << "after2: " << found->value(map_->value_bits_) << std::endl;
//                printf("it is happening.. tid %llu  oldtid %llu  mgid %llu  mgpos %llu  oldpos %llu\n", taxid, old_internal_taxid, marker_gene_id, marker_gene_position, old_gene_pos);
//                exit(9);
//            }
        } else {
            if (main_key == 70516668 && bucket_key == 42007363) {
                std::cout << "Not Found" << std::endl;
//                map_->print(map_->offset_[42007363], 10, true);

                exit(9);
            }


            value |= marker_gene_position;
            value |= marker_gene_id << (marker_gene_position_bits_);
            value |= internal_taxid << (marker_gene_id_bits_ + marker_gene_position_bits_);
            assert(!taxonomy_->Get(internal_taxid).rank.empty());
            map_->Insert(key, value);
        }
    };

    static inline void ParseHeader(std::string &header, size_t &internal_taxid, size_t &marker_gene_id) {
        auto tokens = Utils::split(header.substr(1, string::npos), "_");

//        assert(tokens.size() >= 2);
        assert(tokens.size() > 0);

        internal_taxid = stoull(tokens[0]);
        marker_gene_id = (tokens.size() > 1) * stoull(tokens[(tokens.size() > 1)]);
    }
};



template<typename T>
class SimpleKmerProcessor : KmerProcessor {

    Taxonomy::IntTaxonomy *taxonomy_;
    MiniMap::MiniMap<T> *map_;

    const size_t internal_taxid_bits_;
    const size_t max_internal_taxid_;

public:
    SimpleKmerProcessor(MiniMap::MiniMap<T>* map, Taxonomy::IntTaxonomy *taxonomy, size_t& internal_taxid_bits) :
            internal_taxid_bits_(internal_taxid_bits),
            max_internal_taxid_((1llu << internal_taxid_bits) - 1),
            map_(map),
            taxonomy_(taxonomy) {
        if (internal_taxid_bits_ > map_->ValueBits()) {
            errx(EX_SOFTWARE,
                 "Bit size for internal taxid (%llu) must be less to the value size declared in the map (%llu).",
                 (unsigned long long) internal_taxid_bits_,
                 (unsigned long long) map_->ValueBits());
        }
    }

    void operator () (uint64_t &key, uint64_t &internal_taxid) override {
        uint64_t value = 0;
        uint64_t main_key = key & ((1llu << map_->KeyBits()) - 1);

//        std::cout << "find" << std::endl;
        // comment in if not multimap
        auto found = map_->Find(key);
//        std::cout << "found? " << found << std::endl;

        uint64_t taxid = internal_taxid;

        if (found) {
//            std::cout << "found" << std::endl;
            auto old_internal_taxid = found->value(map_->ValueBits());
            if (old_internal_taxid != internal_taxid) {
                taxid = taxonomy_->LCA(old_internal_taxid, internal_taxid);
            }

//            std::cout << "lca: " << taxid << taxonomy_->Get(taxid).scientific_name << std::endl;

            assert(!taxonomy_->Get(taxid).rank.empty());
            value |= taxid;
            found->paste(key, value, map_->ValueBits());
        } else {
//            std::cout << "new" << std::endl;
            value |= internal_taxid;
//            assert(!taxonomy_->Get(internal_taxid).rank.empty());
            map_->Insert(key, value);
        }
    };

    static inline void ParseHeader(std::string &header, size_t &internal_taxid) {
        auto tokens = Utils::split(header.substr(1, string::npos), "_");
        assert(tokens.size() > 0);
        internal_taxid = stoull(tokens[0]);
    }
};



class MultiMapMarkerKmerProcessor : MarkerKmerProcessor {
    using InsertKeyExists = std::function<size_t (const size_t, const size_t)>;

    MultiMap::MultiMap *map_;

    const size_t internal_taxid_bits_;
    const size_t marker_gene_id_bits_;
    const size_t marker_gene_position_bits_;

    const size_t max_internal_taxid_;
    const size_t max_marker_gene_id_;
    const size_t max_marker_gene_position_;

public:
    MultiMapMarkerKmerProcessor(MultiMap::MultiMap *map, size_t internal_taxid_bits, size_t marker_gene_id_bits, size_t marker_gene_position_bits) :
            internal_taxid_bits_(internal_taxid_bits),
            marker_gene_id_bits_(marker_gene_id_bits),
            marker_gene_position_bits_(marker_gene_position_bits),
            max_internal_taxid_((1llu << internal_taxid_bits) - 1),
            max_marker_gene_id_((1llu << marker_gene_id_bits) - 1),
            max_marker_gene_position_((1llu << marker_gene_position_bits) - 1),
            map_(map) {
        if ((marker_gene_position_bits_ + marker_gene_id_bits_ + internal_taxid_bits_) > map_->ValueBits()) {
            errx(EX_SOFTWARE,
                 "Bit sizes for internal taxid (%llu), marker gene id (%llu) and marker gene position (%llu) must add up (or be less) to the value size declared in the map (%llu).",
                 (unsigned long long) internal_taxid_bits_,
                 (unsigned long long) marker_gene_id_bits_,
                 (unsigned long long) marker_gene_position_bits_,
                 (unsigned long long) map_->ValueBits());

        }
    }

    void operator () (uint64_t &key, uint64_t &internal_taxid, uint64_t &marker_gene_id, uint64_t &marker_gene_position) override {
        uint64_t value = 0;

        value |= marker_gene_position;
        value |= marker_gene_id << (marker_gene_position_bits_);
        value |= internal_taxid << (marker_gene_id_bits_ + marker_gene_position_bits_);

        map_->Insert(key, value);
    };

    static inline void ParseHeader(std::string &header, size_t &internal_taxid, size_t &marker_gene_id) {
        auto tokens = Utils::split(header.substr(1, string::npos), "_");

        assert(tokens.size() > 0);

        internal_taxid = stoull(tokens[0]);
        marker_gene_id = (tokens.size() > 1) * stoull(tokens[(tokens.size() > 1)]);
    }
};