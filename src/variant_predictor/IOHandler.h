//
// Created by fritsche on 30/11/2021.
//

#ifndef VARKIT_IOHANDLER_H
#define VARKIT_IOHANDLER_H

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <ostream>
#include <fstream>
#include <cstring>

namespace IO {
    struct IOSNP {
        size_t read_id = 0;
        uint32_t tax_id = 0;
        uint32_t gene_id = 0;
        uint32_t snp_position = 0;
        uint8_t snp_quality = 0;
        char snp_base = 'X';
        uint8_t read_num = false;

        IOSNP() {};

        IOSNP (size_t read_id, uint32_t tax_id, uint32_t gene_id, uint32_t snp_position, char snp_quality, char snp_base, uint8_t read_num):
                read_id(read_id),
                tax_id(tax_id),
                gene_id(gene_id),
                snp_position(snp_position),
                snp_quality(snp_quality),
                snp_base(snp_base),
                read_num(read_num) {}

        std::string Header() {
            return "read_id\ttaxonomic_id\tgene_id\tsnp_position\tsnp_quality\tsnp_base\tfirst_read";
        }

        static inline IOSNP FromLine(std::string &line) {
            auto snp = IOSNP();
            int start = 0;
            int stop = 0;

            while (++stop != '\t');
            snp.read_id = std::stoull(line.substr(start, stop - start));
            start = stop + 1;
            while (++stop != '\t');
            snp.tax_id = std::stoul(line.substr(start, stop - start));
            start = stop + 1;
            while (++stop != '\t');
            snp.gene_id = std::stoul(line.substr(start, stop - start));
            start = stop + 1;
            while (++stop != '\t');
            snp.snp_position = std::stoul(line.substr(start, stop - start));
            start = stop + 1;
            while (++stop != '\t');
            snp.snp_quality = std::stoul(line.substr(start, stop - start));
            start = stop + 1;
            snp.snp_base = std::stoul(line.substr(start, stop - start));
            start = stop + 1;
            snp.read_num = (uint8_t) std::stoul(line.substr(start, stop - start));

            return snp;
        }

        std::string ToString() const {
            return std::to_string(read_id) + "\t" + std::to_string(tax_id) + "\t" + std::to_string(gene_id) + "\t" + std::to_string(snp_position) + '\t' + std::to_string((int) snp_quality) + '\t' + snp_base + '\t' + std::to_string(int(read_num));
        }
    };

    class SNPReader {
    public:
        inline static bool Read(std::istream& is, IOSNP& snp) {
            is.read((char*) &snp, sizeof(IO::IOSNP));
            return !is.eof();
        }
    };

    struct ClassificationLine {
        static constexpr int GENEID_UNSET = 0;
        static constexpr uint32_t TAXID_UNSET = 0;
        static constexpr int GENEPOS_UNSET = INT32_MIN;
        uint32_t record_id = 0;
        uint32_t mutation_count = 0;

        uint32_t taxid = TAXID_UNSET;
        int geneid = GENEID_UNSET;
        int genepos = GENEPOS_UNSET;

        uint32_t total_hits_best = 0;
        uint32_t total_hits = 0;
        uint32_t leaf_hits_best = 0;
        uint32_t read_length = 0;
        uint32_t overlap_with_gene = 0;

        double rank_confidence = 0;
        double candidate_confidence = 0;
        double classification_confidence = 0;

        char classified = 'U';
        bool forward = false;
        uint8_t read_num = 0;

        double predicted_ani = -1;

        static const size_t buffer_size = 20;
        char buffer[buffer_size];

        void SetHeader(std::string &header) {
            size_t len = buffer_size - 1;
            if (header.length() < buffer_size - 1)
                len = header.length();
            memcpy(buffer, header.c_str(), len);
            buffer[len] = '\0';
        }

        void Reset() {
            record_id = 0;
            mutation_count = 0;
            taxid = TAXID_UNSET;
            geneid = GENEID_UNSET;
            genepos = GENEPOS_UNSET;
            total_hits_best = 0;
            total_hits = 0;
            leaf_hits_best = 0;
            read_length = 0;
            overlap_with_gene = 0;
            rank_confidence = 0.0;
            candidate_confidence = 0.0;
            classification_confidence = 0.0;
            predicted_ani = 0;
            classified = 'U';
            read_num = 0;
        }

        bool HasGene() const {
            return geneid != GENEID_UNSET;
        }
        bool HasPos() const {
            return genepos != GENEPOS_UNSET;
        }

        std::string ToString() {
            return (std::to_string(record_id) + " mutations: "
                    + std::to_string(mutation_count) + " taxid: "
                    + std::to_string(taxid) + " geneid: "
                    + std::to_string(geneid) + " genepos: "
                    + std::to_string(genepos) + " totalhitsbest: "
                    + std::to_string(total_hits_best) + " leafhitsbest: "
                    + std::to_string(leaf_hits_best) + " totalhits: "
                    + std::to_string(total_hits) + " rl: "
                    + std::to_string(read_length) + " ov: "
                    + std::to_string(overlap_with_gene) + " header: "
                    + std::string(buffer));



        }
//
//        void SetPairedHit(size_t taxonomic_id, size_t gene_id, classification::HitPE &first, classification::Hit &second, FastxRecord& record1, FastxRecord& record2) {
////            classified = 'C';
////            read_length = record1.sequence.length();
////            record_id = record_id;
////            mutation_count = mutation_count;
////
////            taxid = classification::NodeId::GetTaxonomicId(first.GetNodeId());
////            geneid = rcr.gene_id;
////            genepos = rcr.gene_pos;
////
////            total_hits_best = rcr.best_total_counts;
////            total_hits = hit_count;
////            leaf_hits_best = rcr.best_count;
////            rank_confidence = rank_confidence;
////            candidate_confidence = candidate_confidence;
////            classification_confidence = classification_confidence;
////            SetHeader(record.id);
//        }
//
//        inline void SetHit(classification::HitPE &hit_pe, FastxRecord& record, size_t record_id, bool second) {
//            auto& hit = hit_pe.GetHit(second);
//            auto [ offset, hits, forward ] = hit.GetOffset().GetHits();
//
//            classified = 'C';
//            read_length = record.sequence.length();
//            this->record_id = record_id;
//            mutation_count = 0;
//
//
//            taxid = hit_pe.m_taxonomic_id;
//            if (hit.GetOffset().IsLeaf() && hits > 2) {
//                geneid = hit_pe.m_gene_id;
//                genepos = offset;
//                total_hits_best = hit.Total();
//                total_hits = hit.Total();
//                leaf_hits_best = hits;
//                rank_confidence = 1;
//                candidate_confidence = 1;
//                classification_confidence = 1;
//            }
//
//            SetHeader(record.id);
//        }
    };

    enum class VcfColumns {
        CHROM = 0,
        POS = 1,
        ID = 2,
        REF = 3,
        ALT = 4,
        QUAL = 5,
        FILTER = 6,
        INFO = 7,
        FORMAT = 8
    };

    struct VcfHeader {
        template <typename A, typename B>
        using kv_pair = std::pair<A, B>;
        using kv_pair_str = kv_pair<std::string, std::string>;
        std::vector<kv_pair_str> header_items;

        void Add(std::string &key, std::string &value) {
            header_items.push_back({ key, value });
        }

        void Write(std::ofstream output) {
            for (auto& kv : header_items)
                output << "##" << kv.first << "=" << kv.second << std::endl;
        }
    };

    struct VcfRecord {
        // Numeric id of the reference sequence.
        int32_t rID;
        // Position on the reference.
        int32_t beginPos;
        // Textual identifier of the variant.
        std::string id;
        // Bases in the reference.
        std::string ref;
        // Bases in the alternatives, comma-separated.
        std::string alt;
        // Quality
        float qual;
        // Value of FILTER field.
        std::string filter;
        // Value of INFO field.
        std::string info;
        // Value of FORMAT field.
        std::string format;
        // The genotype infos.
        std::vector<std::string> genotypeInfos;
    };
};


#endif //VARKIT_IOHANDLER_H
