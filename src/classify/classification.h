//
// Created by fritsche on 25/04/2022.
//

#pragma once

#include <sparse_map.h>
#include "TaxonomyNew.h"
#include <cstring>
#include <iostream>

#include "FastxReader.h"
#include "Utils.h"
#include "data_structures.h"
#include "ShapeBuilder.h"
#include "ShapeHistory.h"
#include "IOHandler.h"

namespace classification {
//    namespace IO {
//        class ClassificationLine;
//    }
    using TID = Taxonomy::TaxId;
    // Unique Node Identifier. Incorporate this into
    using GID = int32_t;
    using GPOS = int32_t;


    struct NodeId {
        // Both Taxid and Geneid get int32_t respectively
        using NodeID_t = uint64_t;
        const NodeID_t data;

        static inline NodeID_t ToNodeId(int32_t taxonomic_id, int32_t gene_id) {
            NodeID_t data = 0;
            static_assert((sizeof(taxonomic_id) + sizeof(gene_id)) == sizeof(data));
            memcpy((char*)&data+sizeof(gene_id), (char*) &taxonomic_id, sizeof(taxonomic_id));
            memcpy(((char*)&data), (char*) &gene_id, sizeof(gene_id));
            return data;
        }

        static inline NodeID_t ToNodeId(int32_t taxonomic_id) {
            NodeID_t data = 0;
            memcpy((char*)&data + sizeof(int32_t), (char*) &taxonomic_id, sizeof(taxonomic_id));
            return data;
        }

        NodeId(int32_t taxonomic_id, int32_t gene_id=0) :
                data(ToNodeId(taxonomic_id, gene_id)) {
            static_assert((sizeof(taxonomic_id) + sizeof(gene_id)) == sizeof(data));
        };

        int32_t& GetGeneId() const {
            return *((int32_t*)(data+4));
        }
        int32_t& GetTaxonomicId() const {
            return *((int32_t*)data);
        }

        static int32_t GetGeneId(NodeID_t data) {
            return *((int32_t*)(&data));
        }

        static int32_t GetTaxonomicId(NodeID_t data) {
            return data >> 32;
        }

        static std::string ToString(NodeID_t data) {
            std::string str;
            str += std::to_string(GetTaxonomicId(data));
            str += ",";
            str += std::to_string(GetGeneId(data));
            return str;
        }
    };
    using UNID = NodeId::NodeID_t;

    /**
     * Prepare lineages for easier access
     */
    class Lineage {
    public:
        constexpr static size_t ROOT_NODE_ID = 1;
        constexpr static size_t LINEAGE_SIZE = 10;
        TID m_lineage[LINEAGE_SIZE];
    private:
        int32_t m_rank_id = -1;
    public:
//        +---+---+---+---+---+---+---+
//        |   | 1 | 2 | 3 | 4 | 5 |   |
//        +---+---+---+---+---+---+---+
//          ↑   ↑               ↑   ↑
//          |   |               |   |
//       rend() |         rbegin()  end()
//          |
//        begin()
        Lineage() {
            std::fill_n(m_lineage, LINEAGE_SIZE, -1);
        }

//        const TID& operator[](int idx) {
//            return m_lineage[idx];
//        }

        TID* begin();
        TID* end();

        TID* rbegin();
        TID* rend();
        void Set(int rank_id, TID rank_taxonomic_id);
        void SetRankId(int rank_id);
        std::string ToString() const;
        std::pair<size_t, TID> GetRankIdAndTaxId() const;
        std::pair<size_t, TID> GetLCARankIdAndTaxId(const Lineage& other) const;
        auto GetRankId() const;
        TID GetTaxonomicId() const;
        bool IsRoot() const;
        const TID GetNextSetTaxonomicId(size_t rank_id) const;
        const TID GetPrevSetTaxonomicId(size_t rank_id) const;
        TID ParentId() const {
            return m_lineage[(m_rank_id - 1) * (m_rank_id != 0)];
        }

        bool NotInit();

        bool IsIdAncestor(TID taxid);
    };
    using LineageTIDMap = Lineage*;

    struct ValueExtractorMG {
    public:
        ValueExtractorMG(size_t internal_taxid_bits,
                         size_t gene_id_bits,
                         size_t gene_pos_bits) :
                internal_taxid_bits_(internal_taxid_bits),
                gene_id_bits_(gene_id_bits),
                gene_pos_bits_(gene_pos_bits) {}

        ValueExtractorMG(const ValueExtractorMG &t) {
            internal_taxid_bits_ = t.internal_taxid_bits_;
            gene_id_bits_ = t.gene_id_bits_;
            gene_pos_bits_ = t.gene_pos_bits_;
        }


        static ValueExtractorMG Load(std::string path) {
            std::ifstream is(path);
            std::string line;

            size_t internal_taxid_bits{};
            size_t gene_id_bits{};
            size_t gene_pos_bits{};

            if (getline(is, line)) {
                internal_taxid_bits = stoull(line);
            }
            if (getline(is, line)) {
                gene_id_bits = stoull(line);
            }
            if (getline(is, line)) {
                gene_pos_bits = stoull(line);
            }
            is.close();

            return ValueExtractorMG(internal_taxid_bits, gene_id_bits, gene_pos_bits);
        }

        void Save(std::string path) {
            std::ofstream ofs(path);

            ofs << internal_taxid_bits_ << '\n';
            ofs << gene_id_bits_ << '\n';
            ofs << gene_pos_bits_;

            ofs.close();
        }

        size_t GetTaxid(size_t value) {
            return value >> (gene_id_bits_ + gene_pos_bits_);
        }

        size_t GetGeneId(size_t value) {
            return value >> (gene_pos_bits_) & ((1 << gene_id_bits_) - 1);
        }

        size_t GetGenePos(size_t value) {
            return value & ((1 << gene_pos_bits_) - 1);
        }

        std::string ToString(size_t value) {
            std::string str;
            str += std::to_string(GetTaxid(value)) + ',';
            str += std::to_string(GetGeneId(value)) + ',';
            str += std::to_string(GetGenePos(value));
            return str;
        }


    public:
        size_t internal_taxid_bits_;
        size_t gene_id_bits_;
        size_t gene_pos_bits_;
    };

//    using TID = int32_t;


    using Taxonomy = Taxonomy::IntTaxonomy;


    /**
     * Class that represents a node in the taxonomic tree.
     */
    template<typename T>
    class ClassificationNode {
        UNID m_node_id = 1;
        UNID m_parent_node_id = 1;
        std::vector<UNID> m_children_node_ids;
        T m_data{};

    public:
        ClassificationNode(UNID node_id, UNID parent_node_id) :
                m_node_id(node_id),
                m_parent_node_id(parent_node_id) {}

        ClassificationNode(UNID node_id, UNID parent_node_id, T&& data) :
                m_node_id(node_id),
                m_parent_node_id(parent_node_id),
                m_data(data) {}

        ClassificationNode(UNID node_id, UNID parent_node_id, T& data) :
                m_node_id(node_id),
                m_parent_node_id(parent_node_id),
                m_data(data) {}

        ClassificationNode() {};

        constexpr UNID GetParentNodeId() const {
            return m_parent_node_id;
        }

        constexpr UNID GetNodeId() const {
            return m_node_id;
        }

        void AddChild(UNID child_node_id) {
            m_children_node_ids.emplace_back(child_node_id);
        }

        void SetParent(UNID parent_node_id) {
            m_parent_node_id = parent_node_id;
        }

        void RemoveChild(UNID child_node_id) {
            m_children_node_ids.erase(std::remove(m_children_node_ids.begin(), m_children_node_ids.end(), child_node_id));
        }

        const std::vector<UNID>& Children() const {
            return m_children_node_ids;
        }

        std::string ToString() {
            constexpr bool has_ToString = requires(T& t) {
                t.ToString();
            };

            if constexpr (has_ToString) {
                return m_data.ToString();
            } else {
                std::cout << "NO TO STRING METHOD" << std::endl;
            }
            return std::string();
        }

        T& Data() {
            return m_data;
        };

        const T& Data() const {
            return m_data;
        };

        std::string ToString() const {
            std::string str;
            str += std::to_string(m_node_id);
            str += "\tParent: ";
            str += std::to_string(m_parent_node_id);
            if (m_children_node_ids.empty()) {
                str += "\tLeaf";
            } else {
                str += "\tChildren(" + std::to_string(m_children_node_ids.at(0));
                for (auto i = 1; i < m_children_node_ids.size(); i++) {
                    str += ',' + std::to_string(m_children_node_ids.at(i));
                }
                str += ")";
            }
            return str;
        }
    };

    struct ReadOffset {
        static constexpr int16_t UNINIT_OFFSET = INT16_MIN;
        int16_t m_offset_forward = UNINIT_OFFSET; //
        int16_t m_offset_reverse = UNINIT_OFFSET;
        int16_t m_hits_forward = 0;
        int16_t m_hits_reverse = 0;
        uint16_t m_initial_hit_pos = -1;
        uint16_t m_distance_first_last_hit = 0;

        // if two k-mers span the same mutation, they might be dependant
        // Check how many independant k-mers there are
        // Note: This is oblivious of the shape and thus underestimates
        // the number of independant k-mers as shapes might interlock
        size_t IndependantKmerEvidence(const size_t shape_length) const {
            return (size_t) m_distance_first_last_hit / shape_length;
        }

        void Set(const GPOS& gene_pos, const GPOS& read_pos) {
            m_offset_forward = gene_pos - read_pos;
            m_offset_reverse = gene_pos + read_pos;
            m_hits_forward++;
            m_hits_reverse++;
            m_initial_hit_pos = read_pos;
        }


        ReadOffset(const GPOS& gene_pos, const GPOS& read_pos) {
            Set(gene_pos, read_pos);
        }

        ReadOffset() {};

        const bool IsLeaf() const {
            return m_offset_forward != UNINIT_OFFSET && m_hits_forward > 1 || m_offset_reverse != UNINIT_OFFSET && m_hits_reverse > 1 ;
        }

        const bool IsInit() const {
            return m_offset_forward != UNINIT_OFFSET || m_offset_reverse != UNINIT_OFFSET;
        }

        // Returns tuple of < offset, #hits, forward/reverse >
        const std::tuple<int16_t, int16_t, bool> GetHits() const {
            if (m_hits_forward > m_hits_reverse) {
                return { m_offset_forward, m_hits_forward, true };
            } else {
                return { m_offset_reverse, m_hits_reverse, false };
            }
        }

        // Returns tuple of < offset, #hits, forward/reverse >
        const std::tuple<int16_t, int16_t, bool> GetHits(const size_t shape_size, const size_t read_length) const {
            if (m_hits_forward > m_hits_reverse) {
                return { m_offset_forward, m_hits_forward, true };
            } else {
                return { m_offset_reverse + shape_size - read_length, m_hits_reverse, false };
            }
        }

        const bool Forward() const {
            return m_hits_forward > m_hits_reverse;
        }

        const std::string ToString() const {
            auto&& [ offset, hits, forward ] = GetHits();
            std::string str;

            str += "offset: " + std::to_string(offset);
            str += "\thits: " + std::to_string(hits) + '\t';
            str += (forward ? "forward" : "reverse");
            return str;
        }

        const std::string ToString2() const {
            std::string str;

            str += "offset_fwd: " + std::to_string(m_offset_forward);
            str += "\toffset_rev: " + std::to_string(m_offset_reverse);
            str += "\thits_fwd: " + std::to_string(m_hits_forward);
            str += "\thits_rev: " + std::to_string(m_hits_reverse);
            str += "\tfirst_hit: " + std::to_string(m_initial_hit_pos);

            return str;
        }

        void Add(GPOS gene_pos, GPOS read_pos) {
            if (m_offset_forward == gene_pos - read_pos) {
                m_hits_forward++;
            }
            if (m_offset_reverse == gene_pos + read_pos) {
                m_hits_reverse++;
            }
        }

        inline void Increment() {
            m_hits_forward++;
        }

        inline bool IncrementIfMatch(GPOS gene_pos, GPOS read_pos) {
            bool match = false;
            if (m_offset_forward == gene_pos - read_pos) {
                m_hits_forward++;
                match = true;
            }
            if (m_offset_reverse == gene_pos + read_pos) {
                m_hits_reverse++;
                match = true;
            }
            if (match) m_distance_first_last_hit = read_pos - m_initial_hit_pos;
            return match;
        }

        inline bool IsMatch(GPOS gene_pos, GPOS read_pos) const {
            auto match =
                    // Observed offset either matches the forward position
                    (m_offset_forward == gene_pos - read_pos) ||
                    // or the reverse position
                    (m_offset_reverse == gene_pos + read_pos);
            return match;
        }

        uint16_t GetMaxHits() const {
            return std::max(m_hits_forward, m_hits_reverse);
        }
    };

    struct Hit {
        size_t m_total_count = 0;
        // Keep read_offset information close
        ReadOffset m_offset;
        // If there are ambiguous hits on species level store them in vector
        std::vector<ReadOffset> m_offsets;

    public:

        std::string ToString() const {
            auto&& [ offset, hits, forward ] = m_offset.GetHits();
            std::string str = "hits." + std::to_string(hits) + "_offset." + std::to_string(offset) + "_" + (forward ? "forward" : "reverse") + "_ambiguous." + std::to_string(!m_offsets.empty());
            return str;
        }

        uint16_t GetMaxHits() const {
            return m_offset.GetMaxHits();
        }

        const ReadOffset& GetOffset() const {
            return m_offset;
        }

        const std::vector<ReadOffset>& GetOffsets() const {
            return m_offsets;
        }

        void Set(GPOS gene_pos, GPOS read_pos) {
            m_offset.Set(gene_pos, read_pos);
        }

        void AddAlternativeOffset(GPOS gene_pos, GPOS read_pos) {
            m_offsets.emplace_back( ReadOffset(gene_pos, read_pos) );
        }

        bool HasAlternativeOffsets() const {
            return !m_offsets.empty();
        }

        bool IsMatch(const GPOS& gene_pos, const GPOS& read_pos) const {
            return m_offset.IsMatch(gene_pos, read_pos);
        }

        bool IncrementIfMatch(const GPOS& gene_pos, const GPOS& read_pos) {
            return m_offset.IncrementIfMatch(gene_pos, read_pos);
        }

        void SwapIfBetter(ReadOffset& offset) {
            if (offset.GetMaxHits() > m_offset.GetMaxHits()) {
                std::swap(offset, m_offset);
            }
        }

        void Add(const GPOS& gene_pos, const GPOS& read_pos) {
            // Add function takes care of everything
            // Check if there is a match in offset
            // if not, add alternative offset
            if (!m_offset.IncrementIfMatch(gene_pos, read_pos)) {
                for (auto& offset : m_offsets) {
                    if (offset.IncrementIfMatch(gene_pos, read_pos)) {
                        SwapIfBetter(offset);
                        return;
                    }
                }
                AddAlternativeOffset(gene_pos, read_pos);
            }
        }

        bool IsLeaf() {
            return m_offset.IsLeaf();
        }

        bool IsInit() {
            return m_offset.IsInit();
        }

        void AddToTotal() {
            m_total_count = m_offset.GetMaxHits();
        }

        void AddToTotalFromOther(const Hit& other) {
            m_total_count = m_offset.GetMaxHits();
            m_total_count += other.m_total_count;
        }

        size_t Total() {
            return m_total_count;
        }

        void Increment() {
            m_offset.Increment();
            assert(m_offsets.empty());
        }
    };

    struct HitPE {
        enum HIT_TYPE { NO_STRAIN, FIRST_STRAIN, SECOND_STRAIN, BOTH_STRAIN };

        TID m_taxonomic_id;
        GID m_gene_id;
        size_t m_total_count = 0;
        double m_rank_confidence = 0;
        double m_candidate_confidence = 0;
        Hit m_read[2];

        HitPE(TID taxonomic_id) :
                m_taxonomic_id(taxonomic_id) {
        }

        HitPE(TID taxonomic_id, GID gene_id) :
                m_taxonomic_id(taxonomic_id),
                m_gene_id(gene_id) {
        }

        const Hit& GetHit(bool second = false) const {
            return m_read[second];
        }

        const bool HasHit(bool second = false) const {
            return m_read[second].m_total_count > 0;
        }

        Hit& GetHit(bool second = false) {
            return m_read[second];
        }

        const size_t TotalHits() const {
            return m_total_count;
        }

        const std::string ToString() const {
            std::string str;
            str += std::to_string(m_total_count);
            return str;
        }

        const std::string ToVerboseString() const {
            std::string str;
            str += std::to_string(m_taxonomic_id) + '\t';
            str += std::to_string(m_gene_id) + '\t';
            str += std::to_string(TotalHits()) + '\t';
            str += m_read[0].ToString() + '\t';
            str += m_read[0].ToString();
            return str;
        }

        const size_t NodeHitsOnly() const {
            return GetHit(0).GetMaxHits() + GetHit(1).GetMaxHits();
        }

        HIT_TYPE GetHitType() {
            auto first_leaf = m_read[0].IsLeaf();
            auto second_leaf = m_read[1].IsLeaf();

            if (first_leaf && !second_leaf) return FIRST_STRAIN;
            if (!first_leaf && second_leaf) return SECOND_STRAIN;
            if (first_leaf && second_leaf) return BOTH_STRAIN;
            return NO_STRAIN;
        }
    };

    using C1Node = ClassificationNode<Hit>;
    using CNode = ClassificationNode<HitPE>;

    struct ClassificationResult {
        static constexpr int GENEID_UNSET = IO::ClassificationLine::GENEID_UNSET;
        static constexpr uint32_t TAXID_UNSET = IO::ClassificationLine::TAXID_UNSET;
        static constexpr int GENEPOS_UNSET = IO::ClassificationLine::GENEPOS_UNSET;

        size_t m_taxonomic_id = TAXID_UNSET;
        int m_gene_id = GENEID_UNSET;
        int m_gene_pos = GENEPOS_UNSET;

        size_t m_total_kmers = 0;
        size_t m_result_kmers = 0;
        size_t m_hit_kmers = 0;
        size_t m_overlap = 0;

        int pattern_idx = -1;
        bool m_forward = true;

        HitPE* m_best_first = nullptr;
        HitPE* m_best_second = nullptr;
        HitPE* m_best = nullptr;

        Hit* m_hit = nullptr;
        Hit* m_read1_hit = nullptr;
        Hit* m_read2_hit = nullptr;

        bool m_success = false;



        std::vector<CNode*> m_hits;
        void SortHits() {
            std::sort(m_hits.begin(), m_hits.end(),
                      [](CNode *a, CNode *b) -> bool {
                          return a->Data().TotalHits() > b->Data().TotalHits();
                      });
        };

        std::vector<CNode*> Hits() {
            return m_hits;
        }

        bool HasHits() const {
            return !m_hits.empty();
        }

        void Clear() {
            m_best_first = nullptr;
            m_best_second = nullptr;
            m_read1_hit = nullptr;
            m_read2_hit = nullptr;
            m_best = nullptr;
            m_taxonomic_id = TAXID_UNSET;
            m_gene_id = GENEID_UNSET;
            m_gene_pos = GENEPOS_UNSET;
            m_hits.clear();
            m_success = false;
            m_hit = nullptr;
            m_total_kmers = 0;
            m_result_kmers = 0;
            m_hit_kmers = 0;
        }


        void SetNodeHit(TID taxonomic_id, size_t kmer_hits_result);
        void SetTotalKmers(size_t total_kmers);
        void SetKmerHits(size_t kmer_hits);
        void SetLeafHit(Hit* hit, GID gene_id);


        bool Success() const {
            return m_success;
        }

        bool IsLeafHit() const {
            return m_hit;
        }

        bool MultipleNodesSameScore() {
            if (m_hits.size() < 2) return false;

            return (m_hits.at(0)->Data().TotalHits() == m_hits.at(1)->Data().TotalHits());
        }

        CNode& GetBestNode() {
            if (m_hits.size() == 1) return *m_hits.at(0);

            auto& first_hit = *m_hits.at(0);
            TID lca_taxid = first_hit.Data().m_taxonomic_id;

//            for (auto i = 0) {
//
//            }
            return first_hit;
        }

        bool HasBothReadsClassified() {
            return m_best;
        }

        HitPE& GetBothClassified() {
            return *m_best;
        }

        bool HasFirstReadClassified() {
            return m_best_first;
        }

        HitPE& GetFirst() {
            return *m_best_first;
        }

        bool HasSecondReadClassified() {
            return m_best_second;
        }

        HitPE& GetSecond() {
            return *m_best_second;
        }

        Hit* GetHitRead1() {
            return m_read1_hit;
        }

        Hit* GetHitRead2() {
            return m_read2_hit;
        }

        void WriteInfo(size_t shape_size, size_t read_length, std::ostream& os=std::cout) const {
            os << "taxid:      \t" << m_taxonomic_id << std::endl;
            if (IsLeafHit()) {
                auto &&[offset, hits, forward] = m_hit->GetOffset().GetHits(shape_size, read_length);
                os << "geneid:           \t" << m_gene_id << std::endl;
                os << "genepos:          \t" << offset << std::endl;
                os << "hits:             \t" << hits << std::endl;
                os << "orientation:      \t" << forward << std::endl;
            }
            os << "taxid:            \t" << m_taxonomic_id << std::endl;
            os << "result kmers:     \t" << m_result_kmers << std::endl;
            os << "hittable kmers:   \t" << m_total_kmers << std::endl;
            os << "total hit kmers:  \t" << m_result_kmers << std::endl;
        }
    };



    class PatternHandler {
        int current_pos_ = 0;
        bool sensitive_mode_ = false;

        static constexpr size_t INIT_HITMISS_SIZE = 100;
        static constexpr size_t INIT_PATTERNS_SIZE = 100;

//    std::vector<uint32_t> patterns_;
        std::vector<uint64_t> patterns_ { INIT_PATTERNS_SIZE } ;

        size_t pattern_len_;

        Taxonomy& taxonomy_;

        std::vector<uint32_t> cov;


    public:
//        using MutationSet = std::set<ds::SNP>;

        using MutationSet = tsl::sparse_set<ds::SNP, ds::HashSNP, ds::CompSNP>;
        MutationSet mutations_;
        int first_hit_ = INT_MAX;
        int last_hit_ = 0;
        int read_length_ = 0;
        int cov_bucket_size = -1;

        size_t last_size = 0;

        std::vector<ds::KmerHit> hitmiss_ { INIT_HITMISS_SIZE };

        const MutationSet &GetMutations() const {
            return mutations_;
        }

        PatternHandler(Taxonomy &taxonomy, std::string shape_path, bool sensitive_mode);
        PatternHandler(Taxonomy &taxonomy, ds::PatternMap* db, bool sensitive_mode);

        void Reset();

        inline void AddNode(int taxid);
        inline void AddLeaf(int taxid, int geneid, int genepos);
        inline void AddMiss();
        inline void Add(size_t pos, bool hit);
        inline double GetHitCov();
        int Evaluate3(ClassificationResult& result, const size_t sequence_length, const size_t shape_size, const size_t gene_length, LineageTIDMap map);
        int Evaluate2(int taxid, int geneid, int abs_pos, bool forward, size_t shape_length, size_t read_length, int start, int n);
        int Evaluate(int taxid, int geneid, int rel_pos, bool forward, int shape_length);

        static const std::pair<uint32_t, uint32_t> OverlappingKmers(bool forward, int gene_pos, size_t read_length, size_t shape_length, size_t gene_length) {
            uint32_t n = std::min(gene_pos + read_length, gene_length) - std::max(0, gene_pos) - shape_length + 1;
            uint32_t start;

            if (gene_pos < 0) {
                if (forward) {
                    start = -gene_pos;
                } else {
                    start = 0;
                }
            } else {
                if (forward) {
                    start = 0;
                } else {
                    start = read_length - n - shape_length + 1;
                }
            }
            return { start, n };
        }

        ds::PatternMap* pattern_db = nullptr;
        ShapeBuilder snp_detector_;
    };

    static void LineFromClassificationResult(IO::ClassificationLine &line, ClassificationResult const& result, FastxRecord const& record, size_t record_id, uint8_t read_num=0) {
        line.taxid = result.m_taxonomic_id;
        line.geneid = result.m_gene_id;
        line.genepos = result.m_gene_pos;
        line.record_id = record_id;
        line.total_hits = result.m_total_kmers;
        line.total_hits_best = result.m_result_kmers;
        line.read_length = record.sequence.length();
        line.forward = result.m_forward;
        line.read_num = read_num;

        // set overlap here
        line.overlap_with_gene = result.m_overlap;
    }

    static void IoSnpFromMutation(const ds::SNP& snp, ClassificationResult const& result, IO::IOSNP& output_snp, FastxRecord &record, size_t record_id, uint8_t read_num=0) {
        if (!snp.IsValid() || !std::unordered_set<char>{'A','C','G','T'}.contains(record.sequence[snp.read_pos_])) {
            std::cout << "\n\nProblemo: " << std::endl;
            std::cout << snp.ToString() << std::endl;
            std::cout << result.m_hit->ToString() << std::endl;
            std::cout << "Gene pos: " << result.m_gene_pos << std::endl;
            std::cout << "Base: " << record.sequence[snp.read_pos_] << "(" << snp.read_pos_ << ")" << std::endl;
            std::cout << "Quality: " << record.quality[snp.read_pos_] << "(" << snp.read_pos_ << ")" << std::endl;
            std::cout << "\n" << std::endl;
        }

        output_snp.read_id = record_id;
        output_snp.tax_id = result.m_taxonomic_id;
        output_snp.gene_id = result.m_gene_id;
        output_snp.read_num = read_num;
        output_snp.snp_position = snp.gene_pos_;
        output_snp.snp_quality = record.quality[snp.read_pos_];
        output_snp.snp_base = result.m_forward ? record.sequence[snp.read_pos_] : KmerUtils::Complement(
                record.sequence[snp.read_pos_]);
    }

    static void LineFromHitPE(IO::ClassificationLine &line, const HitPE& hit_pe, size_t record_id, FastxRecord &record1, FastxRecord &record2, size_t shape_size, size_t total_hits) {
        constexpr bool debug = false;

        auto& first = hit_pe.GetHit(0);
        auto& second = hit_pe.GetHit(1);

        line.Reset();

        assert(hit_pe.TotalHits() > 0);



        int64_t offset_first = first.GetOffset().Forward() ? std::get<0>(first.GetOffset().GetHits()) : std::get<0>(first.GetOffset().GetHits()) + shape_size - record1.sequence.length();
        int64_t offset_second = second.GetOffset().Forward() ? std::get<0>(second.GetOffset().GetHits()) : std::get<0>(second.GetOffset().GetHits()) + shape_size - record2.sequence.length();

//        std::cout << offset_first << "," << offset_second << std::endl;

        if (first.GetOffset().IsLeaf() && second.GetOffset().IsLeaf()) {
            line.genepos = std::min(offset_first, offset_second);
            line.read_length = record1.sequence.length() + record2.sequence.length();
            line.geneid = hit_pe.m_gene_id;
        } else if (first.GetOffset().IsLeaf()) {
            line.genepos = offset_first;
            line.read_length = record1.sequence.length();
            line.geneid = hit_pe.m_gene_id;
        } else if (second.GetOffset().IsLeaf()) {
            line.genepos = offset_second;
            line.read_length = record2.sequence.length();
            line.geneid = hit_pe.m_gene_id;
        }


//        if(!first.GetOffset().IsLeaf() && !second.GetOffset().IsLeaf()) {
//            std::cerr << "LINE FROM HITPE REQUIRES AT LEAST ONE LEAF HIT" << std::endl;
//            exit(12);
//        }
        line.classified = 'C';
        line.record_id = record_id;
        line.mutation_count = 0;

        line.taxid = hit_pe.m_taxonomic_id;

        // Total hits for best candidate
        line.total_hits_best = hit_pe.TotalHits();
        line.total_hits = total_hits;
        line.leaf_hits_best = first.GetMaxHits() + second.GetMaxHits();
        line.rank_confidence = hit_pe.m_rank_confidence;
        line.candidate_confidence = hit_pe.m_candidate_confidence;
//        line.classification_confidence = classification_confidence;

        line.SetHeader(record1.id);

        if (!first.GetOffset().IsInit() && !second.GetOffset().IsInit() && hit_pe.m_gene_id > 0) {
            std::cout << hit_pe.ToString() << std::endl;
            std::cout << line.ToString() << std::endl;
            std::cout << "NOT LEAF1: " << std::endl;
            std::cout << first.GetOffset().IsLeaf() << "," << second.GetOffset().IsLeaf() << std::endl;
//            exit(123);
        }
    };

    static std::string PrintCNode(const ClassificationNode<HitPE>& node) {
        std::string str;
        auto taxonomic_id = NodeId::GetTaxonomicId(node.GetNodeId());
        auto gene_id = NodeId::GetGeneId(node.GetNodeId());
        auto data = node.Data();

        str += std::to_string(taxonomic_id) + '\t';
        str += std::to_string(node.GetNodeId()) + '\t';
        str += std::to_string(data.m_taxonomic_id);

        auto& first_read_offset = node.Data().GetHit(0).GetOffset();
        auto& second_read_offset = node.Data().GetHit(1).GetOffset();

        auto& [ offset, hits, forward ] = first_read_offset.GetHits();
        auto& [ offset2, hits2, forward2 ] = second_read_offset.GetHits();

        if (first_read_offset.IsLeaf() || second_read_offset.IsLeaf()) {
            str += "     gene: ";
            str += std::to_string(gene_id);
        }

        if (first_read_offset.IsLeaf()) {
            str += "  1(";
            str += (forward ? ">> " : "<< ");
            str += std::to_string(offset);
            str += ",";
            str += std::to_string(first_read_offset.m_distance_first_last_hit);
            str += ") -" + std::to_string(hits) + "-  ";
        }

        if (second_read_offset.IsLeaf()) {
            str += "  2(";
            str += (forward2 ? ">> " : "<< ");
            str += std::to_string(offset2);
            str += ",";
            str += std::to_string(second_read_offset.m_distance_first_last_hit);
            str += ") -" + std::to_string(hits2) + "-  ";
//            str += " fl(";
//            str += std::to_string(first_read_offset.m_initial_hit_pos) + ",";
//            str += std::to_string(first_read_offset.m_initial_hit_pos + first_read_offset.m_distance_first_last_hit) + "=";
//            str += std::to_string(first_read_offset.m_distance_first_last_hit);
//            str += ")  ";
        }
        if (!(first_read_offset.IsLeaf() || second_read_offset.IsLeaf())) {
            str += " -" + std::to_string(hits) + "," + std::to_string(hits2) + "-" ;
        }

//        str += " -" + std::to_string(hits) + "-" ;
        if (node.Data().GetHit(0).HasAlternativeOffsets()) {
            str += " (alternatives: " + std::to_string(node.Data().GetHit(0).GetOffsets().size()) + ")";
        }
        if (node.Data().m_total_count) {
            str += "  +" + std::to_string(node.Data().m_total_count) + "+";
        }

        return str;
    }

    static std::string PrintC1Node(const ClassificationNode<Hit>& node) {
        std::string str;
        auto taxonomic_id = NodeId::GetTaxonomicId(node.GetNodeId());
        auto gene_id = NodeId::GetGeneId(node.GetNodeId());
        auto& data = node.Data();
        auto& [ offset, hits, forward ] = node.Data().m_offset.GetHits();

        str += std::to_string(taxonomic_id);
        if (node.Data().m_offset.IsLeaf()) {
            str += "     gene: ";
            str += std::to_string(gene_id) + " offset: ";
            str += std::to_string(offset) + (forward ? " >> " : " << ");
            str += " fl(";
            str += std::to_string(node.Data().m_offset.m_initial_hit_pos) + ",";
            str += std::to_string(node.Data().m_offset.m_initial_hit_pos + node.Data().m_offset.m_distance_first_last_hit) + "=";
            str += std::to_string(node.Data().m_offset.m_distance_first_last_hit);
            str += ")  ";
        }
        str += " -" + std::to_string(hits) + "-" ;
        if (node.Data().HasAlternativeOffsets()) {
            str += " (alternatives: " + std::to_string(node.Data().m_offsets.size()) + ")";
        }
        if (node.Data().m_total_count) {
            str += "  +" + std::to_string(node.Data().m_total_count) + "+";
        }

        return str;
    }

    /**
     * Data structure holding the tree.
     * @tparam T
     */
    template<typename T>
    class ClassificationTree {
//        using Map = tsl::sparse_map<UNID, T>;

        using Node = ClassificationNode<T>;
        using Map = tsl::sparse_map<UNID, Node>;

        Map m_map;
        TID m_root_id = 1;
        TID m_max_leaf_id = -1;
        UNID m_root_node_id;

    public:
        ClassificationTree(TID max_leaf_id);
        void Insert(UNID node_id, Node&& node);
        void Insert(UNID node_id, Node& node);
        bool HasNode(UNID);
        void Clear();

        bool HasRoot();
        TID& RootId();
        UNID& RootNodeId();
        void InsertRoot(T&& data);

        std::string ToString();

        bool HasLeaf(UNID taxonomic_id) {
            return m_map.contains(taxonomic_id);
        }

        Map& Data() {
            return m_map;
        }

        Node& GetNode(UNID node_id) {
            if (!m_map.contains(node_id)) {
                std::cout << "Node does not exist " << NodeId::ToString(node_id) << std::endl;
            }
            assert(m_map.contains(node_id));
            return m_map.at(node_id);
        }

        Node& Root() {
            if (!HasRoot()) {
                std::cout << "ListNodes: " << m_map.size() << std::endl;
                ListNodes();
                std::cout << "______" << std::endl;
            }
            assert(HasRoot() && "Has no root.");
            return m_map.at(m_root_node_id);
        }

        bool IsLeaf(TID taxonomic_id, GID gene_id, GPOS gene_position) {
            return taxonomic_id > m_root_id && taxonomic_id <= m_max_leaf_id &&
                gene_id > 0 && gene_position > 0;
        }

        Node& GetLeaf(UNID leaf_id) {
            if (!m_map.contains(leaf_id)) {
                std::cout << "Leaf does not exist " << leaf_id << std::endl;
                exit(9);
            }
            assert(m_map.contains(leaf_id));
            return m_map.at(leaf_id);
        }

        void PrettyPrint(const std::string& prefix, const ClassificationNode<T>& node, bool isLast, std::function<std::string(const  ClassificationNode<T>& node)> nodeToString) {
            std::cout << prefix;
            std::cout << (!isLast ? "├──" : "└──" );

            // print the value of the node
            std::cout << nodeToString(node) << std::endl;

            // enter the next tree level - left and right branch
            if (!node.Children().empty()) {
                for (auto i = 0; i < node.Children().size() - 1; i++) {
                    auto&& child = GetNode(node.Children().at(i));
                    assert(child.GetNodeId() != node.GetNodeId());
//                    auto new_prefix = prefix + (!isLast ? "│   " : "    ");
//                    PrettyPrint(new_prefix, child, false);
                    PrettyPrint( prefix + (!isLast ? "│   " : "    "), child, false, nodeToString);
                }
                auto& child_node_id = node.Children().at(node.Children().size()-1);
                assert(child_node_id != node.GetNodeId());
                PrettyPrint( prefix + (!isLast ? "│   " : "    "), GetNode(child_node_id), true, nodeToString);
            }
        }

        void PrettyPrint(std::function<std::string(const ClassificationNode<T>& node)> nodeToString) {
            PrettyPrint("", Root(), true, nodeToString);
        }

        void ListNodes() {
            for (auto& node : m_map) {
                std::cout << NodeId::ToString(node.first) << "\tChildren(";
                for (auto& child : node.second.Children())
                    std::cout << NodeId::ToString(child) << ",";
                std::cout << ")" << std::endl;
            }
        }
    };

    template<typename T>
    void ClassificationTree<T>::Insert(UNID node_id, Node&& node) {
        m_map.insert( { node_id, node } );
    }

    template<typename T>
    void ClassificationTree<T>::Insert(UNID node_id, Node& node) {
        m_map.insert( { node_id, node } );
    }

    template<typename T>
    void ClassificationTree<T>::Clear() {
        m_map.clear();
    }

    template<typename T>
    inline bool ClassificationTree<T>::HasNode(UNID node_id) {
        return m_map.contains(node_id);
    }

    template<typename T>
    std::string ClassificationTree<T>::ToString() {
        std::string str;
        for (auto& [id, node] : m_map) {
            str += NodeId::ToString(id);
            str += ": ";
//            str += node.ToString();
            str += "\n";
        }
        return str;
    }

    template<typename T>
    ClassificationTree<T>::ClassificationTree(TID max_leaf_id) :
            m_root_id(1),
            m_root_node_id(NodeId::ToNodeId(m_root_id)),
            m_max_leaf_id(max_leaf_id) {
    }

    template<typename T>
    TID &ClassificationTree<T>::RootId() {
        return m_root_id;
    }

    template<typename T>
    bool ClassificationTree<T>::HasRoot() {
        return m_map.contains(m_root_node_id);
    }

    template<typename T>
    UNID &ClassificationTree<T>::RootNodeId() {
        return m_root_node_id;
    }

    template<typename T>
    void ClassificationTree<T>::InsertRoot(T&& data) {
        Insert(RootNodeId(), Node(RootNodeId(), RootNodeId(), data));
    }



    using C1Tree = ClassificationTree<Hit>;



    /**
     * Main class handling the recording of k-mer hits from the database for each read
     */
    class Classifier {
    public:
        Classifier(Taxonomy &taxonomy);
        ~Classifier();
        void AddHit(TID, GID gene_id, GPOS gene_position, GPOS read_position);
        void AddMiss();
        bool HasNode(UNID);
        void LoadLineage();
        Lineage* GetLineages();
        void ListNodes();
        void PrintTree(C1Node& node, unsigned int level);
        void PrintTree();
        C1Tree &Tree();
        std::string Newick();
        std::string Newick(C1Node& node, unsigned int level);
        size_t Dummy();
        void Clear();
        bool IsAmbiguous();
        void PrintAmbiguous();
        void InitializePatternProcessor(ds::PatternMap* pattern_db, bool sensitive_mode=true);
        ReadOffset& Evaluate();

        bool m_consistent_strains = true;

    private:
        Taxonomy& m_taxonomy;
        C1Tree m_tree;
        LineageTIDMap m_lineages;
        std::unique_ptr<PatternHandler> processor = nullptr;

        void AddLeaf(UNID id, TID taxonomic_id, GID gene_id, GPOS gene_pos, GPOS read_pos);
        void AddNode(UNID id, TID taxonomic_id, GID gene_id);
        void AddNode(UNID id, TID taxonomic_id, GID gene_id, Lineage& lineage, size_t init_hit=0);
    };


    class LookupResult {
    private:
        Lineage m_lineage;
        TID m_taxonomic_id;
        GID m_gene_id;

        HitPE m_hit;

    public:
        LookupResult(const Lineage& lineage, TID taxonomic_id, GID gene_id, GPOS gene_position, GPOS read_position, bool second=false);
        LookupResult(const Lineage& lineage, TID taxonomic_id, GID gene_id, bool second);
        LookupResult(const Lineage& lineage, TID taxonomic_id, GID gene_id);

        Hit& GetHit(bool second=false);
        HitPE& GetHitPE();
        const HitPE& GetHitPE() const;
        const Lineage &GetLineage() const;
        const GID& GeneId() const;
        const TID& GetTaxonomicId() const;
        const std::string ToString() const;
    };


    using CTree = ClassificationTree<HitPE>;

    using GENE_LEN = uint16_t;
    class Gene {
        static constexpr GENE_LEN UNSET_VALUE = UINT16_MAX;
        std::vector<GENE_LEN> m_lengths;
    public:
        Gene(const size_t gene_count) {
            Init(gene_count);
        }

        Gene() {};

        void Init(const size_t gene_count) {
            m_lengths.resize(gene_count);
            std::fill(m_lengths.begin(), m_lengths.end(), UNSET_VALUE);
        }

        void Set(GID gene_id, GENE_LEN length) {
            m_lengths[gene_id] = length;
        }
        GENE_LEN GetLength(GID gene_id) {
            return m_lengths[gene_id];
        }
        bool HasGene(GID gene_id) {
            return m_lengths[gene_id] != UNSET_VALUE;
        }
    };

    class GeneLengths {
        using GENE_LENGTHMAP = std::vector<Gene>;
        GENE_LENGTHMAP map;
        size_t m_gene_count;

        GeneLengths(const size_t genome_count, const size_t gene_count) : m_gene_count(gene_count) {
            map.resize(genome_count, Gene(gene_count));
        }

        void Load(const std::string path) {
            std::ifstream is(path, std::ios::in);

            std::vector<std::string> tokens;
            for (std::string line; std::getline(is, line); ) {
                Utils::split(tokens, line, "\t");
                if (tokens.size() != 3) {
                    std::cout << "Exception" << std::endl;
                    exit(2);
                }
                TID taxonomic_id = std::stoi(tokens[0]);
                GID gene_id = std::stoul(tokens[1]);
                GENE_LEN gene_length = std::stoul(tokens[2]);

                Set(taxonomic_id, gene_id, gene_length);
            }

            is.close();
        }

        void Set(TID taxonomic_id, GID gene_id, GENE_LEN length) {
            assert(gene_id < m_gene_count);
            if (taxonomic_id > map.size()) {
                map.resize(taxonomic_id+1, Gene(m_gene_count));
            }
            map[taxonomic_id].Set(gene_id, length);
        }

    public:
        GENE_LEN Get(TID taxonomic_id, GID gene_id) {
            return map[taxonomic_id].GetLength(gene_id);
        }

        GeneLengths(const std::string path, const size_t genome_count, const size_t gene_count) : m_gene_count(gene_count) {
            map.resize(genome_count+1, Gene(gene_count));
            Load(path);
        }
    };

    class Classifier2 {
        using LookupList = std::vector<LookupResult>;
        using LookupIdxMap = tsl::sparse_map<UNID, size_t>;


    private:
        Taxonomy& m_taxonomy;
        CTree m_tree;
        LineageTIDMap m_lineages = nullptr;

        size_t m_lineages_size = 0;

        LookupIdxMap m_lookup_idx_map;
        LookupList m_lookup_results;

        size_t m_max_leaf_id;
        size_t m_read_length;
        size_t m_cov_bucket_size;

        unique_ptr<PatternHandler> m_processor[2];

    public:
        Classifier2(Taxonomy &taxonomy);

        void LoadLineage();
        LineageTIDMap GetLineages();
        void SetLineage();
        // Record successful k-mer lookup (Hit)
        void AddHit(TID, GID gene_id, GPOS gene_position, GPOS read_position, bool second=false);
        // Record unsuccessful k-mer lookup (Miss)
        void AddMiss(bool reverse=false);

        // Is Leaf or Node
        bool IsLeaf(TID taxonomic_id, GID gene_id, GPOS gene_position);

        void AddLeaf(TID taxonomic_id, GID gene_id, GPOS gene_position, GPOS read_position, bool second=false);
        void AddNode(TID taxonomic_id, GID gene_id, bool second=false);
        bool HasNode(TID taxonomic_id);
        CTree &GetTree();

        void ProcessWorker();
        void BuildTree();
        void AddUpTree(ClassificationResult &result);
        void AddUpTree(ClassificationResult& result, CNode& parent, CNode& target);
        void ProcessHits(ClassificationResult& result);
        bool ProcessHits2(ClassificationResult& result, double min_confidence=0.7, size_t min_hits=3);
        bool ProcessHits3(ClassificationResult& result1, ClassificationResult& result2, double min_confidence=0.7, size_t min_hits=3);


        std::string Newick();
        std::string Newick(CNode& node, unsigned int level);

        void InitializePatternProcessor(ds::PatternMap* pattern_db, bool sensitive_mode=true);
        void SetReadLength(size_t l, bool second);

        const unique_ptr<PatternHandler> &GetProcessor(bool second=false) const;

        std::vector<LookupResult>& GetLookupResults();
        void PrintLookupResults() const;


        ~Classifier2() {
            delete[] m_lineages;
        };

        size_t Dummy() {
            return m_tree.Data().size();
        };

        void Clear();

        bool m_consistent_strains = true;

        pair<size_t, size_t> GetSumOfHitsForTaxon(TID taxid, GID gene_id);
    };
};

