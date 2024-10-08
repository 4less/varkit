//
// Created by fritsche on 03/07/2022.
//

#pragma once

#include <string>
#include <vector>
#include "OptionsContainer.h"
#include "IOHandler.h"
#include <iostream>
#include "GenomeLoader.h"
#include "TaxonomyNew.h"


class SimpleSNPFilter;
namespace varkit {
    class SNP;
    class SNPAllele {
    private:
        uint32_t m_observations = 0;
        uint32_t m_quality_sum = 0;
        char m_base;
    public:
//        SNPAllele(char base, uint32_t quality) :
//                m_observations(1),
//                m_base(base),
//                m_quality_sum(quality) {}

        SNPAllele(char base) :
                m_base(base) {}

        void AddObservation(uint32_t quality) {
//            std::cout << "Add observation: " << ToString() << std::endl;
            m_quality_sum += quality;
            m_observations++;
        }

        std::string ToString() const {
            return std::string(1, m_base) + "(" + std::to_string(m_observations) + ")";
        }

        char Base() const {
            return m_base;
        }
        uint32_t Count() const {
            return m_observations;
        }
        double AvgQual() {
            return (double)m_quality_sum/m_observations;
        }
        friend class SNP;
    };

    class SNP {
//        std::vector<char> m_base = {};
//        std::vector<uint16_t> m_observations = {};
//        uint16_t m_observations = 0;
//        std::vector<uint8_t> m_qualities;

        using Alleles = std::vector<SNPAllele>;

        char m_reference = 'X';
        int m_position = -1;
        // TODO figure out more clever way to do this
        mutable Alleles m_alleles;
        size_t m_vertical_coverage = 0;

        SNPAllele& GetOrAdd(char allele_base) {
            for (auto& allele : m_alleles) {
                if (allele.Base() == allele_base) return allele;
            }
            m_alleles.emplace_back(SNPAllele(allele_base));
            return m_alleles.back();
        }

    public:
        int Position() const {
            return m_position;
        }

        bool operator < (const SNP& snp) const
        {
            return (m_position < snp.Position());
        }
        bool operator > (const SNP& snp) const {
            return (m_position < snp.Position());
        }
        bool operator >= (const SNP& snp) const {
            return (m_position >= snp.Position());
        }
        bool operator <= (const SNP& snp) const {
            return (m_position <= snp.Position());
        }
        bool operator == (const SNP& snp) const {
            return (m_position == snp.Position());
        }
        bool operator != (const SNP& snp) const {
            return (m_position != snp.Position());
        }

        const Alleles& GetAlleles() const {
            return m_alleles;
        }

        size_t CountAlleles() const {
            return m_alleles.size();
        }
        size_t CountTotalObservations() const {
            return std::accumulate(m_alleles.begin(), m_alleles.end(), 0, [](int count, SNPAllele const& allele) { return count + allele.Count(); });
        }

        size_t CountAllelesNoRef() const {
            return std::count_if(m_alleles.begin(), m_alleles.end(), [this](SNPAllele const& a){ return a.Base() != 'R' && a.Base() != m_reference; });
        }

        SNPAllele& MajorAlleleNoRef() const {
            auto iter = std::max_element(m_alleles.begin(), m_alleles.end(), [](SNPAllele const& a, SNPAllele const& b) { return a.m_observations < b.m_observations && a.Base() != 'R'; });
            return *iter;
        }

        std::string ToString() const {
//            size_t allele_sum = std::accumulate(m_alleles.begin(), m_alleles.end(),; 0, [](int acc, SNPAllele const& a) { return acc + a.Count(); });
            return std::to_string(m_position) + " (" + std::to_string(m_vertical_coverage) + ") [" + Utils::Join(m_alleles, ",", [](SNPAllele const& allele){ return allele.ToString(); }) + " Ref: " + m_reference + " All: " + std::to_string(CountAllelesNoRef()) + " ]";
//            return std::to_string(m_position) + ": " + m_base.at(0) + "(" + std::to_string(m_observations.at(0)) + ")";
//            return std::to_string(m_position) + ": " + m_base.at(0) + "(" + std::to_string(m_observations.at(0)) + ")";
        }

        void SetReference(char ref) {
            m_reference = ref;
        }

        void Add(IO::IOSNP const& iosnp) {
//            if (m_position == 1691) {
//                std::cout << iosnp.ToString() << std::endl;
////                std::cout << ToString() << " " << (iosnp.snp_position == m_position) << " " << iosnp.tax_id << std::endl;
//            }
            auto& allele = GetOrAdd(iosnp.snp_base);
            allele.AddObservation(iosnp.snp_quality);
            m_position = iosnp.snp_position;

//            m_base = iosnp.snp_base;
//            m_observations++;
//            m_qualities.emplace_back(iosnp.snp_quality);
        };

        void AddReference(size_t observations) {
            auto iter = std::find_if(m_alleles.begin(), m_alleles.end(), [this](SNPAllele const& allele) {
                return m_reference == allele.Base();
            });
            if (iter != m_alleles.end()) {
                std::cout << ToString() << std::endl;
                iter->m_base = 'R';
                Utils::Input();
            } else {
                m_alleles.emplace_back(SNPAllele('R'));
                m_alleles.back().m_observations = observations;
            }
        }

        void SetVerticalCoverage(size_t vertical_coverage) {
            m_vertical_coverage = vertical_coverage;
        }

        size_t GetVerticalCoverage() const {
            return m_vertical_coverage;
        }

        Alleles& GetAlleles() {
            return m_alleles;
        }

        void AddTotalCoverage(size_t vcov) {
            auto sum = std::accumulate(m_alleles.begin(), m_alleles.end(), 0, [](size_t acc, SNPAllele const& allele){ return acc + allele.Count(); });
            if (vcov < sum) {
                std::cout << "____________________Sum: " << sum << std::endl;
                std::cout << this->ToString() << std::endl;
                std::cout << "vcov: " <<  vcov << " snps: " << sum << std::endl;
                std::cout << "Vcov from reads must be larger or equal to number of SNPs at that position" << std::endl;
                //TODO this is a temporary solution
                vcov = sum;
            }

            auto ref_count = vcov - sum;
            SetVerticalCoverage(vcov);
            if (ref_count > 0) {
                AddReference(ref_count);
            }
        }
    };

    struct ReadInfo {
        size_t read_id;
        size_t start;
        size_t length;
        bool forward;

        std::string ToString() const {
            return "read_id: " + std::to_string(read_id) + " start: " + std::to_string(start) + " length: " + std::to_string(length) + " forward: " + std::to_string(forward);
        }

        bool Contains(size_t pos) const {
//            std::cout << "pos: " << pos << " contains: " << (pos >= start && pos <=(start + length)) << " " << start << " len: " << length << std::endl;
            return pos >= start && pos < (start + length);
        }
    };

    class SequenceRangeHandler;
    class SequenceRange {
        size_t m_start = 0;
        size_t m_end = 0;

        // Collect read evidence
        mutable std::vector<ReadInfo> m_read_info;

        bool Overlap(size_t start, size_t end) const;

        bool Overlap(SequenceRange const& range) const {
            return Overlap(range.m_start, range.m_end);
        }


        void Union(SequenceRange const& range) {
            assert(*this == range);
            m_start = std::min(m_start, range.m_start);
            m_end = std::max(m_end, range.m_end);

            // Copy elements from src to dest.
            // (Could move, but don't want to leave old object
            // in undefined state)
            m_read_info.insert(
                    m_read_info.end(),
                    range.m_read_info.begin(),
                    range.m_read_info.end()
            );
        }
        void Intersect(SequenceRange const& range) {
            assert(*this == range);
            m_start = std::max(m_start, range.m_start);
            m_end = std::min(m_end, range.m_end);

            auto end = std::remove_if(m_read_info.begin(),
                                      m_read_info.end(),
                                      [this](ReadInfo const& read_info) {
                                          return !this->Contains(read_info.start);    // remove odd numbers
                                      });

            m_read_info.erase(end, m_read_info.end());

        }

    public:
        SequenceRange(size_t start, size_t end) : m_start(start), m_end(end) {};

        void AddReadInfo(ReadInfo const&& read_info) {
            m_read_info.emplace_back(read_info);
        }

        void AddReadInfo(ReadInfo const& read_info) {
            m_read_info.emplace_back(read_info);
        }

        const vector<ReadInfo> GetReadInfo() const {
            return m_read_info;
        }

        void SortReadInfo() const {
            std::sort(m_read_info.begin(), m_read_info.end(), [](ReadInfo const& a, ReadInfo const& b) {
                return a.start < b.start;
            });
        }

        using CoverageVec = std::vector<uint16_t>;
        CoverageVec CoverageVector() const {
            CoverageVec coverage;

            std::vector<size_t> ends;

            for (auto& read_info : GetReadInfo()) {
                ends.emplace_back(read_info.start + read_info.length);
            }
            std::sort(ends.begin(), ends.end());

            size_t current_coverage = 0;
            size_t read_info_index = 0;
            size_t ends_index = 0;

            coverage.resize(m_end - m_start);

            for (auto i = 0; i < m_end - m_start; i++) {

                auto start = m_read_info[read_info_index].start - m_start;
                while (read_info_index < m_read_info.size() && start <= i) {
                    current_coverage++;
                    read_info_index++;
                    start = m_read_info[read_info_index].start - m_start;
                }

                auto end = ends[ends_index] - m_start;
                while (ends_index < ends.size() && end <= i) {
                    current_coverage--;
                    ends_index++;
                    end = ends[ends_index] - m_start;
                }
//                coverage.emplace_back(current_coverage);
                coverage[i] = current_coverage;
            }

            return coverage;
        }

        bool operator < (const SequenceRange& range) const {
            return (m_end < range.m_start);
        }
        bool operator > (const SequenceRange& range) const {
            return (m_start > range.m_end);
        }
        bool operator == (const SequenceRange& range) const {
            return Overlap(range);
        }
        bool operator != (const SequenceRange& range) const {
            return !(*this == range);
        }
        bool operator <= (const SequenceRange& range) const {
            return *this < range || *this == range;
        }
        bool operator >= (const SequenceRange& range) const {
            return *this > range || *this == range;
        }

        size_t Length() const {
            return m_end - m_start;
        }

        std::string ToString() const {
            return "[" + std::to_string(m_start) + ',' + std::to_string(m_end) + "](" + std::to_string(m_read_info.size()) + ")";
        }

        std::string ToVerboseString() const {
            std::string output = "[" + std::to_string(m_start) + ',' + std::to_string(m_end) + "] {\n";
            for (auto& read_info : m_read_info) {
                output += read_info.ToString() + "\n";
            }
            return output;
        }

        size_t Start() const {
            return m_start;
        }

        size_t End() const {
            return m_end;
        }

        friend class SequenceRangeHandler;

        bool Contains(size_t pos) const {
            return pos >= m_start && pos < m_end;
        }

        bool Contains(size_t start, size_t end) const {
            return start >= m_start && end < m_end;
        }

        size_t Coverage(size_t pos) const {
            auto count = std::count_if(m_read_info.begin(), m_read_info.end(), [&pos](ReadInfo const& read_info) {
                return read_info.Contains(pos);
            });
//            std::cout << "count: " << count << " total: " << m_read_info.size() << std::endl;
            return count;
        }
    };


    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator() (const std::pair<T1, T2> &pair) const {
            return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
        }
    };

    using SNPList = std::vector<SNP>;
//    using SNPKey = std::pair<size_t, char>;
    using SNPKey = size_t;
//    using SNPMap = tsl::sparse_map<SNPKey, SNP, pair_hash>;
    using SNPMap = tsl::sparse_map<SNPKey, SNP>;
    using SequenceRangeList = std::vector<SequenceRange>;

    /**
     * Handles sequence ranges.
     */
    class SequenceRangeHandler {
    private:
        using SRIterator = const std::vector<SequenceRange>::iterator;

        // Variables
        mutable SequenceRangeList m_ranges;

        // Methods
        SRIterator GetSequenceRangeIterator(const size_t start, const size_t end);
        SRIterator GetOrAddSequenceRangeIterator(const size_t start, const size_t end);

    public:
        // Iterator wrapper
        SRIterator begin() noexcept;
        SRIterator end() noexcept;

        // Navigation
        SRIterator FindSequenceRange(size_t start, size_t end);
        const SequenceRangeList& GetRanges() const;
        SequenceRangeList& GetRanges();

        SequenceRangeHandler Intersect(SequenceRangeHandler const& handler) const;
        SequenceRangeHandler Unify(SequenceRangeHandler const& handler) const;

        const SequenceRange& GetRange(size_t pos) const;
        SequenceRange& GetRange(size_t pos);
        bool HasRange(size_t pos) const;

        size_t Size() const;
        size_t SequenceLength() const;
        void Add(size_t start, size_t end);
        void Add(const IO::ClassificationLine &line);
        void Add(SequenceRange& range);
        void Add(SequenceRange&& range);

        std::string ToString() const;
        std::string ToVerboseString() const;

    };


    // Gene holds the SNPs for a Gene
    class Gene {
    private:
        // Variables
        std::string m_name{};
        int m_id = -1;
        size_t m_gene_length = 0;

        // SNPs in gene
        mutable SNPList m_snps;
        mutable SNPMap m_snps_map;
        // Sequence overlap
        SequenceRangeHandler m_ranges;

    public:
        Gene() {};
        Gene(size_t gene_length) : m_gene_length(gene_length) {};

        // Store SNPs in vector after loading into memory
        void SNPsToVector() const;

        // Modifier
        void AddRange(IO::ClassificationLine const& line);
        void AddSNP(IO::IOSNP const& iosnp);
        const void CompareSNPsToReference(std::string const& reference) const;
        SNPList GetSNPsInRange(SequenceRangeList const& ranges) const {
            SNPList snps;
            std::copy_if(m_snps.begin(), m_snps.end(), std::back_inserter(snps),
                         [&ranges](SNP const& snp){
                return std::any_of(ranges.begin(), ranges.end(), [&snp](SequenceRange const& range) {
                    return range.Contains(snp.Position());
                });
            });
            return snps;
        }

        // Getter & Setter
        SNP& GetSNP(SNPKey key);
        const SequenceRangeHandler& GetRangeHandler() const;
        SequenceRangeHandler GetRangeHandlerForCoverage(size_t minimum_coverage, size_t minimum_length) const;
        const SNPList& GetSNPList() const;
        SNPList& GetSNPList();
        SequenceRangeHandler& GetRangeHandler();
        void SetGeneLength(size_t length);
        size_t GetGeneLength() const;

        void AddReferenceCountToSNPs() const {
            for (auto& snp : m_snps) {
                auto& range = m_ranges.GetRange(snp.Position());
                auto total = range.Coverage(snp.Position());

                snp.AddTotalCoverage(total);
            }
        }
    };

    using GeneID = size_t;
    using GeneMap = tsl::sparse_map<GeneID, Gene>;


    /**
     * Genes class
     */
    class Genes {
        GeneMap m_genes;

    public:
        Genes() {};

        bool HasGene(GeneID id) const;

        void AddGene(GeneID id, size_t gene_length);
        Gene& GetOrAddGene(GeneID id) noexcept;

        // Getter
        Gene& GetGene(GeneID id);
        const Gene& GetGene(GeneID id) const;
        GeneMap& GetGeneMap();
        const GeneMap& GetGeneMap() const;


    };

    class Strain {
        std::string m_id{};
        int m_taxonomic_id;
        Genes m_genes;

    public:
        Strain() {};

        Genes& GetGenes();
        const Genes& GetGenes() const;

        void ToVerbose(std::ostream &os = std::cout) const {
            std::string result;
            std::cout << "Number of genes: " << m_genes.GetGeneMap().size() << std::endl;
            for (auto& gene : m_genes.GetGeneMap()) {
                os << gene.first << '\t';
            }
            os << std::endl;
            for (auto& gene : m_genes.GetGeneMap()) {
                os << gene.second.GetRangeHandler().SequenceLength() << '\t';
            }
            os << std::endl;
            for (auto& gene : m_genes.GetGeneMap()) {
                os << gene.second.GetRangeHandler().GetRanges().size() << '\t';
            }
            os << std::endl;
            for (auto& gene : m_genes.GetGeneMap()) {
                os << gene.second.GetSNPList().size() << '\t';
            }
            os << std::endl;
        }
    };

    struct GeneDistanceReturn {
        double distance = 0.0;
        double shared_sequence_percent = 0.0;

        size_t shared_snps = 0;
        size_t snps_a = 0;
        size_t snps_b = 0;

        size_t length_a = 0;
        size_t length_b = 0;
        size_t shared_length = 0;

        size_t shared_genes = 0;

        GeneDistanceReturn(
                double distance,
                double shared_sequence_percent,
                size_t snps_a,
                size_t snps_b,
                size_t shared_snps,
                size_t length_a,
                size_t length_b,
                size_t shared_length,
                size_t shared_genes) :
            distance(distance),
            shared_sequence_percent(shared_sequence_percent),
            snps_a(snps_a),
            snps_b(snps_b),
            shared_snps(shared_snps),
            length_a(length_a),
            length_b(length_b),
            shared_length(shared_length),
            shared_genes(shared_genes) {};

        std::string ToString() const {
            return "Distance: " + std::to_string(distance) +
                    "\tShared sequence percentage: " + std::to_string(shared_sequence_percent) +
                    "\tsnps_a: " + std::to_string(snps_a) +
                    "\tsnps_b: " + std::to_string(snps_b) +
                    "\tshared_snps: " + std::to_string(shared_snps) +
                    "\tlength_a: " + std::to_string(length_a) +
                    "\tlength_b: " + std::to_string(length_b) +
                    "\tTotal shared sequence: " + std::to_string(shared_length) +
                    "\tShared genes: " + std::to_string(shared_genes);
        }
    };


    class SimpleSNPFilter {
        using Alleles = std::vector<SNPAllele>;
        double m_threshold = 0.2;
        size_t m_min_observations = 2;
    public:
        SimpleSNPFilter() {}

        /**
         * Filter SNPAllele by assessing if SNP is likely to be an error.
         * Simple ratio
         * @param snp
         * @return
         */
        Alleles operator()(SNP const& snp) const {
            Alleles alleles;
            size_t vcov = snp.GetVerticalCoverage();
            for (auto& allele : snp.GetAlleles()) {
                if ((double)allele.Count()/vcov >= m_threshold && allele.Count() >= m_min_observations && allele.Base() != 'R') {
                    alleles.emplace_back(allele);
                }
            }
//            if (alleles.size() < snp.GetAlleles().size()) {
//                std::cout << "\nUnfiltered, vcov: " << vcov << std::endl;
//                std::cout << snp.ToString() << std::endl;
//
//                std::cout << "\nFiltered" << std::endl;
//                for (auto& a : alleles) {
//                    std::cout << a.ToString() << std::endl;
//                }
//
//                Utils::Input();
//            }

            return alleles;
        }
    };


    class SimplestGeneDistance {
        SimpleSNPFilter filter;
    public:
        GeneDistanceReturn operator()(Gene const& a, Gene const& b) const {

            auto a_cov = a.GetRangeHandlerForCoverage(2, 30);
            auto b_cov = b.GetRangeHandlerForCoverage(2, 30);

//            auto range_intersection = a_cov.Intersect(b_cov);
//            auto range_union = a_cov.Unify(b_cov);

            auto range_intersection = a.GetRangeHandler().Intersect(b.GetRangeHandler());
            auto range_union = a.GetRangeHandler().Unify(b.GetRangeHandler());

            auto snps_a = a.GetSNPsInRange(range_intersection.GetRanges());
            auto snps_b = b.GetSNPsInRange(range_intersection.GetRanges());

            size_t unique_snps_a = 0;
            size_t unique_snps_b = 0;
            size_t shared_snps = 0;

            auto a_begin = snps_a.begin(), a_end = snps_a.end();
            auto b_begin = snps_b.begin(), b_end = snps_b.end();

            static auto process_unique = [&](SNP const& snp, size_t &uniques) {
                auto alleles = filter(snp);
                if (alleles.size() > 0) uniques++;
            };

            while (a_begin != a_end || b_begin != b_end) {
                if (a_begin == a_end) {
                    // Do B Stuff

                    // Filter SNPs. Only return alleles that are not likely to be result of sequencing errors.
                    process_unique(*b_begin, unique_snps_b);

                    std::advance(b_begin, 1);
                    continue;
                }
                if (b_begin == b_end) {
                    // Do A Stuff

                    // Filter SNPs. Only return alleles that are not likely to be result of sequencing errors.
                    process_unique(*a_begin, unique_snps_a);

                    std::advance(a_begin, 1);
                    continue;
                }
                if (*a_begin == *b_begin) {
//                    std::cout << (*a_begin).ToString() << " == " << (*b_begin).ToString() << std::endl;


                    // Filter SNPs. Only return alleles that are not likely to be result of sequencing errors.
                    auto alleles_a = filter(*a_begin);
                    auto alleles_b = filter(*b_begin);

                    // Both alleles are not empty
                    if (!alleles_a.empty() && !alleles_b.empty()) {

                        auto& major_a = alleles_a.size() == 1 ? alleles_a[0] :
                                  *(std::max_element(alleles_a.begin(), alleles_a.end(), [](SNPAllele const& a, SNPAllele const& b) {
                                      return a.Count() < b.Count();
                                  }));
                        auto& major_b = alleles_a.size() == 1 ? alleles_a[0] :
                                  *(std::max_element(alleles_a.begin(), alleles_a.end(), [](SNPAllele const& a, SNPAllele const& b) {
                                      return a.Count() < b.Count();
                                  }));

                        if (major_a.Base() ==  major_b.Base()) {
                            shared_snps++;
                        } else {
                            unique_snps_a++;
                            unique_snps_b++;
                        }

                    } else if (alleles_a.empty()) {
                        unique_snps_b++;
                    } else if (alleles_b.empty()) {
                        unique_snps_a++;
                    } else {
                        std::cout << "Both are errors" << std::endl;
                    }

                    //shared_snps++;
                    std::advance(a_begin, 1);
                    std::advance(b_begin, 1);
                } else if (*a_begin < *b_begin) {
                    // Filter SNPs. Only return alleles that are not likely to be result of sequencing errors.
                    process_unique(*a_begin, unique_snps_a);
                    std::advance(a_begin, 1);
                } else {
                    // Filter SNPs. Only return alleles that are not likely to be result of sequencing errors.
                    process_unique(*b_begin, unique_snps_b);
                    std::advance(b_begin, 1);
                }
            }

            size_t difference_snps = unique_snps_a + unique_snps_b;
            double distance = (double)(difference_snps)/ range_intersection.SequenceLength();

            auto return_val = GeneDistanceReturn(
                    distance,
                    (double) range_intersection.SequenceLength() / range_union.SequenceLength(),
                    unique_snps_a,
                    unique_snps_b,
                    shared_snps,
                    a.GetRangeHandler().SequenceLength(),
                    b.GetRangeHandler().SequenceLength(),
                    range_intersection.SequenceLength(),
                    1);

//            std::cout << return_val.ToString() << std::endl;
//            std::cout << Utils::Join(a.GetSNPList(), "; ", [](SNP const& snp) { return snp.ToString(); }) << std::endl;
//            std::cout << Utils::Join(b.GetSNPList(), "; ", [](SNP const& snp) { return snp.ToString(); }) << std::endl;
//            Utils::Input();

            return return_val;
        }
    };


    class SimpleGeneDistance {
        SimpleSNPFilter filter;
    public:
        static bool Move(Gene const& a, Gene const& b, SequenceRangeHandler const& intersection, size_t& index_a, size_t& index_b, size_t& index_range) {
            constexpr bool debug = false;

            // Move a_i and b_i to the right to find a SNP within the range
            auto &a_snps = a.GetSNPList();
            auto &b_snps = b.GetSNPList();
            auto &ranges = intersection.GetRanges();

            while (index_range < ranges.size() &&
                    !(index_a < a_snps.size() && ranges[index_range].Contains(a_snps[index_a].Position()) && ((index_b < b_snps.size() && b_snps[index_b].Position() >= ranges[index_range].Start()) || (index_b == b_snps.size()))) &&
                    !(index_b < b_snps.size() && ranges[index_range].Contains(b_snps[index_b].Position()) && ((index_a < a_snps.size() && a_snps[index_a].Position() >= ranges[index_range].Start()) || (index_a == a_snps.size())))) {

                if constexpr(debug) {
                    std::cout << "First expression:  "
                              << (index_a < a_snps.size() && ranges[index_range].Contains(a_snps[index_a].Position()) &&
                                  ((index_b < b_snps.size() &&
                                    b_snps[index_b].Position() >= ranges[index_range].Start()) ||
                                   (index_b == b_snps.size()))) << std::endl;
                    std::cout << "Second expression: "
                              << (index_b < b_snps.size() && ranges[index_range].Contains(b_snps[index_b].Position()) &&
                                  ((index_a < a_snps.size() &&
                                    a_snps[index_a].Position() >= ranges[index_range].Start()) ||
                                   (index_a == a_snps.size()))) << std::endl;

                    std::cout << "Infinite move " << "range(" << index_range << "," << ranges.size() << ","
                              << ranges[index_range].ToString() << ") a(" << index_a << "," << a_snps.size() << ","
                              << (index_a < a_snps.size() ? a_snps[index_a].Position() : 0) << ") b(" << index_b << ","
                              << b_snps.size() << "," << (index_b < b_snps.size() ? b_snps[index_b].Position() : 0)
                              << ")" << std::endl;
                }
                // Iterate while index_range is still valid
                // And both SNPs at index_a and index_b are invalid
                // If one SNP is valid it must be processed before moving
                // This happens if e.g. index range catch up (B) overshot and is both larger
                // than SNPs at index_a and index_b

                auto &tmp_range = ranges[index_range];
                // (A) Try moving into range
                while (index_a < a_snps.size() && a_snps[index_a].Position() < tmp_range.Start()) index_a++;
                while (index_b < b_snps.size() && b_snps[index_b].Position() < tmp_range.Start()) index_b++;
                // (B) Catch up with range in case both SNPs moved past the range
                while (index_range < ranges.size() &&
                        (index_a == a_snps.size() || ranges[index_range].End() <= a_snps[index_a].Position()) &&
                        (index_b == b_snps.size() || ranges[index_range].End() <= b_snps[index_b].Position()) ) index_range++;
            }
            // If at this point, either one of the snps is in a range or its invalid because both snps
            // are
            bool success = (index_a < a_snps.size() ||
                           index_b < b_snps.size()) &&
                           index_range < ranges.size();

            if constexpr(debug) {
                std::cout << "Infinite move at end " << "range(" << index_range << "," << ranges.size() << ","
                          << ranges[index_range].ToString() << ") a(" << index_a << "," << a_snps.size() << ","
                          << a_snps[index_a].Position() << ") b(" << index_b << "," << b_snps.size() << ","
                          << b_snps[index_b].Position() << ")" << std::endl;
            }
            return success;
        }

        GeneDistanceReturn operator()(Gene const& a, Gene const& b) const {
            constexpr bool debug = false;


            auto range_intersection = a.GetRangeHandler().Intersect(b.GetRangeHandler());
            auto range_union = a.GetRangeHandler().Unify(b.GetRangeHandler());

            if constexpr(debug) {
                std::cout << range_intersection.SequenceLength() << "/" << range_union.SequenceLength() << " = "
                          << (double) range_intersection.SequenceLength() / range_union.SequenceLength() << "\t\t";
            }
            size_t a_i = 0;
            size_t b_i = 0;
            size_t range_i = 0;

            size_t range_snps_count_a = 0;
            size_t range_snps_count_b = 0;
            size_t shared_snps_count = 0;

            auto& a_snps = a.GetSNPList();
            auto& b_snps = b.GetSNPList();
            auto& ranges = range_intersection.GetRanges();

            // Move into position.

            while (Move(a, b, range_intersection, a_i, b_i, range_i)) {
//                auto& a_snp = a_snps[a_i];
//                auto& b_snp = b_snps[b_i];
                auto& range = ranges[range_i];

                bool a_allowed = a_i < a_snps.size();
                bool b_allowed = b_i < b_snps.size();
                bool a_within = a_allowed && range.Contains(a_snps[a_i].Position());
                bool b_within = b_allowed && range.Contains(b_snps[b_i].Position());

                bool a_smaller = a_allowed && (!b_allowed || a_snps[a_i].Position() < b_snps[b_i].Position());
                bool b_smaller = b_allowed && (!a_allowed || b_snps[b_i].Position() < a_snps[a_i].Position());

                if constexpr(debug) {
                    std::cout << "a_within: " << a_within << std::endl;
                    std::cout << "b_within: " << b_within << std::endl;
                    std::cout << "a_smaller: " << a_smaller << std::endl;
                    std::cout << "b_smaller: " << b_smaller << std::endl;
                }
                if (a_within && b_within && a_snps[a_i].Position() == b_snps[b_i].Position()) {
                    // Process SNP A and B
                    if constexpr(debug) {
                        std::cout << "shared snp pos" << std::endl;
                    }
                    shared_snps_count++;
                    a_i++;
                    b_i++;
                } else if (a_within && a_smaller) {
                    // Process SNP A
                    if constexpr(debug) {
                        std::cout << "only a" << std::endl;
                    }
                    range_snps_count_a++;
                    a_i++;
                } else if (b_within && b_smaller){
                    // Process SNP B
                    if constexpr(debug) {
                        std::cout << "only b" << std::endl;
                    }
                    range_snps_count_b++;
                    b_i++;
                }
            }

            // Difference SNPs are SNPs that are unique to either strain a or strain b
            // and thus without the snps shared between strain a and strain b.
            auto difference_snps = range_snps_count_a + range_snps_count_b;

            if constexpr(debug) {
                std::cout << "Shared SNPs: " << shared_snps_count << std::endl;
                std::cout << "A excl SNPs: " << range_snps_count_a << std::endl;
                std::cout << "B excl SNPs: " << range_snps_count_b << std::endl;
                std::cout << "A SNPS: " << a.GetSNPList().size() << " " << (shared_snps_count + range_snps_count_a)
                          << std::endl;
                std::cout << "B SNPS: " << b.GetSNPList().size() << " " << (shared_snps_count + range_snps_count_b)
                          << std::endl;
                std::cout << "A to Ref\t-> " << (double)(range_snps_count_a + shared_snps_count)/
                                                range_intersection.SequenceLength() << std::endl;
                std::cout << "B to Ref\t-> " << (double)(range_snps_count_b + shared_snps_count)/
                                                range_intersection.SequenceLength() << std::endl;
                std::cout << "A to B\t-> " << (double)(difference_snps)/ range_intersection.SequenceLength() << std::endl;
            }

            double distance = (double)(difference_snps)/ range_intersection.SequenceLength();

            return GeneDistanceReturn(
                    distance,
                    (double) range_intersection.SequenceLength() / range_union.SequenceLength(),
                    range_snps_count_a,
                    range_snps_count_b,
                    shared_snps_count,
                    a.GetRangeHandler().SequenceLength(),
                    b.GetRangeHandler().SequenceLength(),
                    range_intersection.SequenceLength(),
                    0);
        }
    };


//    class SimpleRefactoredGeneDistance {
//    private:
//
//    public:
//        GeneDistanceReturn operator()(Gene const& a, Gene const& b) const {
//            constexpr bool debug = false;
//        }
//    };


    // Std distance set as SimpleGeneDistance
    using StdGeneDistance = SimpleGeneDistance;
    using StrainKey = size_t;
    using StrainSet = tsl::sparse_set<StrainKey>;


    template<typename GeneDistance=StdGeneDistance>
    class SimpleStrainDistance {
        using GeneDistanceReturnList = std::vector<GeneDistanceReturn>;
    private:
        GeneDistanceReturn Summarize(GeneDistanceReturnList pairwise_distances) const {
            size_t total_shared_length = 0;
            double normalized_distance_sum = 0;
            double normalized_overlap_sums = 0;

            size_t total_length_a = 0;
            size_t total_length_b = 0;
            size_t total_length_shared = 0;

            size_t total_snps_a = 0;
            size_t total_snps_b = 0;
            size_t total_shared_snps = 0;

            size_t shared_genes = pairwise_distances.size();

            for (auto& distance : pairwise_distances) {
                normalized_distance_sum += distance.shared_length * distance.distance;
                normalized_overlap_sums += distance.shared_length * distance.shared_sequence_percent;

                total_snps_a += distance.snps_a;
                total_snps_b += distance.snps_b;
                total_shared_snps += distance.shared_snps;

                total_length_a += distance.length_a;
                total_length_b += distance.length_b;
                total_shared_length += distance.shared_length;
            }
//            double strain_distance = normalized_distance_sum / total_shared_length;
//            double shared_percentage = normalized_overlap_sums / total_shared_length;
            double strain_distance = (double)(total_snps_a + total_snps_b) / total_shared_length;
            double shared_percentage = normalized_overlap_sums / total_shared_length;

//            std::cout << "total_length_a: " << total_length_a << std::endl;
//            std::cout << "total_length_b: " << total_length_b << std::endl;
//            std::cout << "total_shared_length: " << total_shared_length << std::endl;

            // TODO update
            return GeneDistanceReturn(
                    strain_distance,
                    shared_percentage,
                    total_snps_a,
                    total_snps_b,
                    total_shared_snps,
                    total_length_a,
                    total_length_b,
                    total_shared_length,
                    shared_genes);
        }
    public:
        using GeneSet = tsl::sparse_set<size_t>;
        GeneDistanceReturn operator()(Strain const& a, Strain const& b) const {
//            GeneSet gene_union{};
            GeneDistanceReturnList gene_distances;
            for (const auto& [key, gene] : a.GetGenes().GetGeneMap()) {
//                std::cout << gene.GetRangeHandler().ToVerboseString();
//                Utils::Input();
                if (b.GetGenes().GetGeneMap().contains(key)) {
//                    std::cout << "GENE DISTANCE " << key << std::endl;
                    gene_distances.emplace_back(GeneDistance()(gene, b.GetGenes().GetGene(key)));
                }
            }

            auto gene_distance_return_summary = Summarize(gene_distances);

            return gene_distance_return_summary;
        }
    };

    template<typename GeneDistance=StdGeneDistance>
    class JukesCantorDistance {
        using GeneDistanceReturnList = std::vector<GeneDistanceReturn>;
    private:
        GeneDistanceReturn Summarize(GeneDistanceReturnList pairwise_distances) const {
            size_t total_shared_length = 0;
            double normalized_distance_sum = 0;
            double normalized_overlap_sums = 0;

            size_t total_length_a = 0;
            size_t total_length_b = 0;
            size_t total_length_shared = 0;

            size_t unique_snps_a = 0;
            size_t unique_snps_b = 0;
            size_t total_shared_snps = 0;

            for (auto& distance : pairwise_distances) {
                normalized_distance_sum += distance.shared_length * distance.distance;
                normalized_overlap_sums += distance.shared_length * distance.shared_sequence_percent;

                unique_snps_a += distance.snps_a;
                unique_snps_b += distance.snps_b;
                total_shared_snps += distance.shared_snps;

                total_length_a += distance.length_a;
                total_length_b += distance.length_b;
                total_shared_length += distance.shared_length;
            }

            size_t total_different_snps = unique_snps_a + unique_snps_b;

            double strain_distance = (double)total_different_snps / total_shared_length;
            double shared_percentage = normalized_overlap_sums / total_shared_length;

            // Jukes Cantor Distance (1969)
            double jc69 = -3.0f/4 * log(1 - 4.0f/3 * strain_distance);

//            std::cout << "Jukes Cantor " << jc69 << " from " << strain_distance << std::endl;

            // TODO update
            return GeneDistanceReturn(
                    jc69,
                    shared_percentage,
                    unique_snps_a,
                    unique_snps_b,
                    total_shared_snps,
                    total_length_a,
                    total_length_b,
                    total_shared_length,
                    pairwise_distances.size());
        }
    public:
        using GeneSet = tsl::sparse_set<size_t>;
        GeneDistanceReturn operator()(Strain const& a, Strain const& b) const {
            GeneDistanceReturnList gene_distances;
            for (const auto& [key, gene] : a.GetGenes().GetGeneMap()) {
                if (b.GetGenes().GetGeneMap().contains(key)) {
                    gene_distances.emplace_back(GeneDistance()(gene, b.GetGenes().GetGene(key)));
                }
            }

            auto gene_distance_return_summary = Summarize(gene_distances);

            return gene_distance_return_summary;
        }
    };

    struct DoubleMatrix {
        using DoubleVectorT = std::vector<double>;
        using DoubleMatrixT = std::vector<DoubleVectorT>;

        using Names = std::vector<std::string>;

        constexpr static double DEFAULT_VALUE = 0.0;

        // matrix[row][col]
        DoubleMatrixT m;
        Names row_names;
        Names col_names;

        DoubleMatrix(size_t rows, size_t cols) :
                m(rows, DoubleVectorT(cols, DEFAULT_VALUE)),
                row_names(rows),
                col_names(cols) {};

        DoubleMatrix(size_t n) :
                m(n, DoubleVectorT(n, DEFAULT_VALUE)),
                row_names(n),
                col_names(n) {};



        std::vector<std::string> ToStringVector(std::string delimiter="\t", bool print_rownames=true) {
            std::vector<std::string> output;
            for (int row = 0; row < m.size(); row++) {
                std::string line = print_rownames ? row_names[row] + delimiter : "";
                line += Utils::Join(m[row], delimiter.c_str());
                output.emplace_back(line);
            }
            return output;
        }

        void Fill(double value) {
            for (int i = 0; i < m.size(); i++) {
                for (int j = 0; j < m.size(); j++) {
                    m[i][j] = value;
                }
            }
        }

        void Convert() {
            for (int i = 0; i < m.size(); i++) {
                for (int j = 0; j < m.size(); j++) {
                    m[i][j] = 1 - m[i][j];
                }
            }
        }

        std::string Header(std::string delimiter="\t") {
            return "samples" + delimiter + Utils::Join(col_names, delimiter.c_str());
        }

        void Write(std::ostream& os, bool print_names=true) {
            std::cout << "Write()" << std::endl;
            std::string delimiter = "\t";
            if (print_names) {
                os << Header(delimiter) << '\n';
            }
            for (auto& line : ToStringVector(delimiter, print_names)) {
                os << line << '\n';
            }
        }

        void Write(std::string& filepath, bool print_names=true) {
            std::ofstream os(filepath, std::ios::out);
            Write(os, print_names);
            os.close();
        }

        size_t Size() {
            return m.size();
        }

        size_t Rows() {
            return m.size();
        }

        size_t Columns() {
            return m.at(0).size();
        }
    };

    // Std distance set as SimpleGeneDistance
    using StdStrainDistance = SimpleStrainDistance<StdGeneDistance>;


    using StrainMap = tsl::sparse_map<StrainKey, Strain>;

    /**
     * Sample Strains class
     * Represents a m_sample
     * Contains a set of strains.
     */
    class SampleStrains {
    private:
        int m_sample_id = -1;
        size_t m_sample_index = 0;
        Sample m_sample;
        StrainMap m_strains;
        StrainSet m_predicted_strains;

    public:
        SampleStrains(Sample& sample) : m_sample(sample) {};

        const tsl::sparse_map<size_t, Strain>::const_iterator begin() const noexcept;
        const tsl::sparse_map<size_t, Strain>::const_iterator end() const noexcept;

        StrainMap& GetStrains();

        const StrainMap& GetStrains() const;
        const StrainSet& GetPredictedStrains() const;

        void AddPredictedStrain(StrainKey strain) {
            if (!m_strains.contains(strain)) {
                auto strain_string = std::to_string(strain);
                errx(EX_SOFTWARE, "Tried to add strain prediciton that is not in list of strains (%s)", strain_string.c_str());
            }
            m_predicted_strains.insert(strain);
        }

        Strain& GetStrain(size_t internal_id);
        Sample& GetSample() {
            return m_sample;
        }
        const Sample& GetSample() const {
            return m_sample;
        }
        Strain& GetOrAddStrain(size_t internal_id) noexcept;
        void LoadSNPs(std::string const path, tsl::sparse_set<size_t>& target_species, tsl::sparse_set<std::pair<size_t, uint8_t>, Utils::pair_hash>& valid_read_ids);
//        void LoadSNPs(std::string const path, tsl::sparse_set<size_t>& target_species, tsl::sparse_set<size_t>& valid_read_ids);
    };

    using SamplesMap = tsl::robin_map<std::string, SampleStrains>;
    class SampleHandler {
        using GenomeLoaderPtr = std::shared_ptr<gsl::GenomeLoader>;
        using InternalTaxonomyPtr = Taxonomy::IntTaxonomy*;

        SamplesMap m_samples;
        GenomeLoaderPtr m_loader;
//        VarkitOptionsContainer& m_options;
        InternalTaxonomyPtr m_taxonomy = nullptr;


        SampleStrains& GetOrAddSample(std::string sample_name);
        StrainSet GetStrainList() const;
        StrainSet GetPredictedStrainList() const;

    public:
//        SampleHandler(VarkitOptionsContainer& options);
        SampleHandler() {};

        void AddSample(std::string sample_name);
        void AddSample(Sample& sample);

        void SetGenomeLoader(GenomeLoaderPtr genome_loader);
        void SetInternalTaxonomy(InternalTaxonomyPtr taxonomy);

        SampleStrains& GetSample(std::string sample_name);
        GenomeLoaderPtr GetGenomeLoader();
        void CompareSNPsToReference();

        template<typename GeneDistance=StdGeneDistance, template<class> class StrainDistance=SimpleStrainDistance>
        void CalculateDistanceMatrices(std::string output_folder) {
            std::cout << "#Samples: " << m_samples.size() << " " << typeid(GeneDistance).name() << std::endl;

            std::cout << "Predicted Strains" << std::endl;
            auto pred_strains = GetPredictedStrainList();

            for (auto strain_key : GetPredictedStrainList()) {
                std::cout << strain_key << std::endl;
            }

            for (auto strain_key : GetPredictedStrainList()) {
                std::cout << "\n__________________________Strain: " << strain_key << std::endl;
                auto distance_matrix = CalculateDistanceMatrix<GeneDistance, StrainDistance>(strain_key);

                if (m_taxonomy) {
                    std::cout << m_taxonomy->LineageStr(strain_key) << std::endl;
                }



                distance_matrix.Convert();
                distance_matrix.Write(std::cout);

                std::string strain_name = std::to_string(strain_key);
                if (m_taxonomy) {
                    strain_name = m_taxonomy->Get(strain_key).scientific_name;
                    std::replace(strain_name.begin(), strain_name.end(), ' ', '_');
                }

                std::string output_file = output_folder + '/' + strain_name + ".tsv";
                distance_matrix.Write(output_file);
            }
        }

        template<typename GeneDistance=StdGeneDistance, template<class> class StrainDistance=SimpleStrainDistance>
        DoubleMatrix CalculateDistanceMatrix(size_t strain_key, bool collapse_samples=false) {
            if (!collapse_samples) {
                DoubleMatrix matrix = DoubleMatrix(m_samples.size());
                matrix.Fill(-1.0);

                for (const auto &[sample_name1, sample1]: m_samples) {
                    auto row_idx = sample1.GetSample().GetIndex();

                    // Set Name in matrix
                    matrix.row_names[row_idx] = sample1.GetSample().GetSampleName();

                    // Exception
                    if (row_idx >= matrix.Rows()) {
                        errx(EX_SOFTWARE, "Sample index is not within matrix row range");
                    }

                    for (const auto &[sample_name2, sample2]: m_samples) {
                        auto col_idx = sample2.GetSample().GetIndex();

                        // Set Name
                        matrix.col_names[col_idx] = sample2.GetSample().GetSampleName();

                        // Exception
                        if (col_idx >= matrix.Columns()) {
                            errx(EX_SOFTWARE, "Sample index is not within matrix column range");
                        }

                        if (sample_name1 == sample_name2) continue;

                        if (sample1.GetStrains().contains(strain_key) && sample2.GetStrains().contains(strain_key)) {
                            std::cout << "SAMPLES " << sample_name1 << " - " << sample_name2 << std::endl;
                            std::cout << "Strain key: " << strain_key << std::endl;
                            auto& strain_a = sample1.GetStrains().at(strain_key);
                            auto& strain_b = sample2.GetStrains().at(strain_key);

//                            strain_a.ToVerbose();
//                            strain_b.ToVerbose();


                            size_t max_genes = std::max(strain_a.GetGenes().GetGeneMap().size(), strain_b.GetGenes().GetGeneMap().size());
                            GeneDistanceReturn return_value = StrainDistance<GeneDistance>()(sample1.GetStrains().at(strain_key),
                                                                           sample2.GetStrains().at(strain_key));


                            std::cout << return_value.ToString() << std::endl;
                            if (return_value.distance != return_value.distance) {
                                std::cout << "stop " << std::endl;
//                                Utils::Input();
                            }

                            // TODO implement this more flexibly, for now this fix has to do it
                            // Require at least 3 shared genes to calculate distance
                            if (return_value.shared_genes < 3) {
                                std::cout << "t00 little shared genes" << std::endl;
                                continue;
                            };

                            matrix.m[row_idx][col_idx] = return_value.distance;
                        } else {
                            //matrix.m[row_idx][col_idx] = 1.0;
                        }
                    }
                }
                return matrix;
            } else {
                return DoubleMatrix(0);
            }
        }
    };

    using StdSampleHandler = SampleHandler;
};


