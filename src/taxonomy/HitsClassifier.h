//
// Created by fritsche on 04/04/2022.
//

#pragma once


#include <cstddef>
#include <cstdint>
#include <robin_map.h>
#include <sparse_map.h>
#include "CustomTaxonomy.h"
#include "TaxonomyNew.h"

using Taxid = uint32_t;

using IntTaxonomy = Taxonomy::IntTaxonomy;



struct Hit {
    Hit() {};
    Hit(Taxid node_id, Taxid parent_id, uint32_t gene_id) :
            m_node_id(node_id),
            m_parent_id(parent_id),
            m_gene_id(gene_id) {};

    Hit(Taxid const& node_id, Taxid const& parent_id,  uint32_t const& gene_id, int read_pos_fwd, int read_pos_rev) :
            m_node_id(node_id),
            m_parent_id(parent_id),
            m_gene_id(gene_id),
            m_read_pos_fwd(read_pos_fwd),
            m_read_pos_rev(read_pos_rev) {};

    Taxid m_node_id = 0;
    Taxid m_parent_id = 0;

    uint32_t m_gene_id = 0;

    int m_read_pos_fwd = 0;
    int m_read_pos_rev = 0;

    uint32_t hits_fwd = 1;
    uint32_t hits_rev = 1;

    uint32_t m_supportive_counts = 0;

    bool m_is_leaf = false;
    bool m_was_visited = false;


    inline size_t Hits() const {
        return std::max(hits_rev, hits_fwd) + m_supportive_counts;
    }

    std::string ToString() const {
        return "Hit: "
                + std::to_string(m_is_leaf) + " "
                + std::to_string(m_node_id)
                + " (" + std::to_string(m_parent_id) + ")\t"
                + std::to_string(m_gene_id) + "\t"
                + std::to_string(m_read_pos_fwd)
                + " (" + std::to_string(hits_fwd) + ")\t"
                + std::to_string(m_read_pos_rev)
                + " (" + std::to_string(hits_rev) + ")\t";
    }
};

class HitsClassifier {
    using HitMap = tsl::sparse_map<size_t, Hit>;
private:
    HitMap m_hits;
    IntTaxonomy &m_taxonomy;

public:
    size_t dummy_var = 0;

    HitsClassifier(IntTaxonomy &taxonomy) : m_taxonomy(taxonomy) {};

    inline void AddHit(
            const Taxid &id,
            const uint32_t &gene_id,
            const uint32_t &gene_pos,
            const uint32_t &read_pos,
            const uint32_t &read_length,
            const uint32_t &shape_length) {

        auto key = GetKey(id, gene_id);


//        if (gene_id > 0)
//            std::cout << key << '\t' << OffsetFwd(gene_pos, read_pos) << "\t" << OffsetRev(gene_pos, read_pos, read_length, shape_length) << '\t' << gene_pos << "\t" << read_pos << std::endl;


        if (m_hits.contains(key)) {
            auto& hit = m_hits[key];

//            std::cout << "hasit: ";// << (hit.m_read_pos_fwd == OffsetFwd(gene_pos, read_pos)) << " " << (hit.m_read_pos_rev == OffsetRev(gene_pos, read_pos, read_length, shape_length)) << std::endl;
//            std::cout << hit.m_read_pos_fwd << ":" << OffsetFwd(gene_pos, read_pos) << "\t";
//            std::cout << hit.m_read_pos_rev << ":" << OffsetRev(gene_pos, read_pos, read_length, shape_length) << std::endl;
            bool has_hit = false;
            if (hit.m_read_pos_fwd == OffsetFwd(gene_pos, read_pos)) {
                has_hit = true;
                hit.hits_fwd++;
            }
            if (hit.m_read_pos_rev == OffsetRev(gene_pos, read_pos, read_length, shape_length)) {
                has_hit = true;
                hit.hits_rev++;
            }
            if (hit.m_is_leaf && !has_hit) {
                std::cerr << "PROBLEMO" << std::endl;
                std::cout << "hasit: " << (hit.m_read_pos_fwd == OffsetFwd(gene_pos, read_pos)) << " " << (hit.m_read_pos_rev == OffsetRev(gene_pos, read_pos, read_length, shape_length)) << std::endl;
                std::cout << hit.m_read_pos_fwd << ":" << OffsetFwd(gene_pos, read_pos) << "\t";
                std::cout << hit.m_read_pos_rev << ":" << OffsetRev(gene_pos, read_pos, read_length, shape_length) << std::endl;
                std::cout << hit.ToString() << std::endl;

                exit(8);
            }

        } else {
            if (gene_id > 0) {
                m_hits.insert({key, { id, id, gene_id, OffsetFwd(gene_pos, read_pos),
                        OffsetRev(gene_pos, read_pos, read_length, shape_length) }});
            } else {
                m_hits.insert({key, { id, id, gene_id }});
                //                std::cout << "m_hits.size(): " << m_hits.size() << std::endl;
//                std::cout << m_hits[key].ToString() << std::endl;
            }
        }
    }

    inline void Reset() {
        m_hits.clear();
    }

    void Evaluate();
    static inline const size_t GetKey(const Taxid &id, const uint32_t &gene_id) {
        static constexpr size_t gene_id_bits = 20;
        static constexpr size_t tax_id_bits = 20;

        size_t key = (static_cast<uint64_t>(id) << gene_id_bits) | gene_id;
        return key;
    }

    inline void EvaluateTree() {
        auto total_taxa = m_hits.size();
        auto processed_taxa = 0;
        auto processed_hits = 0;

        struct linmem {
            Taxid id = 0;
            uint32_t count = 0;
            linmem() {};
            linmem(Taxid id, uint32_t count) : id(id), count(count) {}
        };

        std::vector<linmem> lineage;


        size_t best = 0;
        std::cout << "Checkpoint a" << std::endl;
        for (auto& [ key, hit ] : m_hits) {
            auto parent_tid = hit.m_node_id;
            auto parent_key = GetKey(parent_tid, 0);
            auto count = m_hits[key].Hits();

            lineage.emplace_back( linmem{ parent_tid, static_cast<uint32_t>(hit.Hits()) } );

            while (parent_tid != 1) {
                parent_tid = m_taxonomy.Get(parent_tid).parent_id;
                parent_key = GetKey(parent_tid, 0);

                if (!m_hits.contains(parent_key) || m_hits[parent_key].m_was_visited)
                    continue;

                lineage.emplace_back( linmem{ parent_tid, static_cast<uint32_t>(hit.Hits()) });
                m_hits[parent_tid].m_was_visited = true;
            }

//            for (auto& lm : lineage) {
//                std::cout << lm.id << "(" << lm.count << ")\t";
//            }
//            std::cout << std::endl;

            size_t last = 0;
            size_t best_count = 0;
            for (auto it = lineage.rbegin(); it != lineage.rend(); it++) {
                auto hit = m_hits[it->id];
                hit.m_supportive_counts += last;
                hit.m_was_visited = true;
                if (hit.Hits() > best_count) {
                    best_count = hit.Hits();
                    best = it->id;
                }
                last = hit.Hits();
            }
            std::cout << "Checkpoint b" << std::endl;

//            for (auto& id : lineage) {
//                std::cout << id << "(" << m_hits[id].Hits() << ")\t";// << m_hits[id].ToString() << "\t";
//            }
//            std::cout << std::endl;
//            std::cout << "--" << std::endl;

            lineage.clear();

            if (processed_taxa == total_taxa) {
//                std::cout << "early break" << std::endl;
                break;
            }
        }
        dummy_var = best;
    }

    static inline const int OffsetFwd(const uint32_t &gene_pos, const uint32_t &read_pos) {
        return gene_pos - read_pos;
    };

    static inline const int OffsetRev(const uint32_t &gene_pos, const uint32_t &read_pos, const uint32_t &read_length, const uint32_t &shape_length) {
        return gene_pos - (read_length - shape_length - read_pos);
    };


};

