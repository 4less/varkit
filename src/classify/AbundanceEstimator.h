//
// Created by fritsche on 07/12/2021.
//

#ifndef VARKIT_ABUNDANCEESTIMATOR_H
#define VARKIT_ABUNDANCEESTIMATOR_H

#include <IOHandler.h>
#include <cmath>
#include <unordered_map>
#include <TaxonomyNew.h>
#include <utils/BinaryClassifierEvaluator.h>
#include "strain_handler.h"


class MarkerGenome {
    int marker_gene_count = 0;
    int max_marker_gene_id = 0;
    int marker_gene_lengths_length = 0;
    std::vector<int> marker_gene_lengths;
    size_t marker_genome_length = 0;

public:
    MarkerGenome() {
//        std::cout << "default constructor: " << this << std::endl;
    };

    MarkerGenome(size_t max_marker_gene_id) :
            max_marker_gene_id (max_marker_gene_id),
            marker_gene_lengths_length (max_marker_gene_id + 1) {
        this->marker_gene_lengths.resize(marker_gene_lengths_length);
        std::fill_n(this->marker_gene_lengths.begin(), this->marker_gene_lengths_length, -1);

//        std::cout << "Constructor (this): " << this << std::endl;
//        std::cout << "marker_gene_lengths(size): " << marker_gene_lengths.size() << std::endl;
    }

    ~MarkerGenome() {
//        std::cout << "destroy MarkerGenome " << (this) << std::endl;
    }

    void Set(size_t geneid, int gene_length) {
        if (marker_gene_lengths.empty()) {
            std::cerr << "invalid write to marker_gene_lengths" << std::endl;
            exit(78);
        }
        marker_gene_lengths[geneid] = gene_length;
    }

    int Get(uint32_t gene_id) {
        if (marker_gene_lengths.empty()) {
            return -1;
        }
        return marker_gene_lengths[gene_id];
    }

    std::string ToString() {
        std::string str_representation =
                "marker_gene_count: " + std::to_string(marker_gene_count)
                + "  max_marker_gene_id: " + std::to_string(max_marker_gene_id)
                + "  marker_gene_lengths_length: " + std::to_string(marker_gene_lengths_length)
                + "  marker_genome_length: " + std::to_string(marker_genome_length) + " ";
        for (auto i = 0; i < 10; i++) {
//            if (Get(i) == -1) {
//                exit(89);
//            }
            str_representation += " [" + std::to_string(i) + "]:" + std::to_string(Get(i));
        }
        return str_representation;
    }

    std::vector<int>& MarkerGeneLengths();

    size_t MarkerGenomeLength() {
        return marker_genome_length;
    }

    size_t MarkerGeneCount() {
        return marker_gene_count;
    }

    void CalculateMGSize() {
        for (auto i = 0; i < marker_gene_lengths_length; i++) {
//            std::cout << "marker_gene_lengths[" << i << "]: " << marker_gene_lengths[i] << " is null? " << (marker_gene_lengths == nullptr) << " Get: " << Get(i) << std::endl;
            if (marker_gene_lengths[i] > 0)
                marker_gene_count++;
            marker_genome_length += marker_gene_lengths[i];
        }
    }

    const int MaxMarkerGeneId() {
        return max_marker_gene_id;
    }
};

class ReadFilter {
private:
    std::string m_name;
    double m_min_ani;
    size_t m_shape_length;
    size_t m_k;

    double m_error_margin = 0.36; // has no tn and no fn for 95 ani
    size_t m_pass_counter = 0;
    size_t m_decline_counter = 0;

    static double EstimatedHits(double ani, size_t read_length, size_t shape_size, size_t k) {
        return (read_length - shape_size + 1) * pow(ani, k);
    }

public:
    ReadFilter(double min_ani, size_t shape_length, size_t k) :
            m_min_ani(min_ani),
            m_shape_length(shape_length),
            m_k(k) {};

    bool Pass(IO::ClassificationLine &line, bool force_pass) {
        // total_hits_best = best candidate number of hits
        // must be larger than Estimated hits * m_error_margin
//        bool pass = line.total_hits_best > (EstimatedHits(m_min_ani, line.read_length, m_shape_length, k_) * (1 - m_error_margin) );
        bool pass = force_pass || line.total_hits_best > (EstimatedHits(m_min_ani, line.overlap_with_gene, m_shape_length, m_k) * (1 - m_error_margin) );

        m_pass_counter += pass;
        m_decline_counter += !pass;

        return pass;
    }

    size_t GetPassCounter() {
        return m_pass_counter;
    }

    size_t GetDeclineCounter() {
        return m_decline_counter;
    }
};

struct Gene {
    std::vector<size_t> buckets;
    size_t geneid = 0;
    size_t gene_length = 0;
    size_t bucket_count = 0;
    size_t bucket_length = 0;
    size_t mapped_bases = 0;
    size_t total_reads = 0;
    double correction_factor = 1.0;

    Gene() {};
    ~Gene() {
//        std::cout << "Delete gene... " << this << std::endl;
    }
    Gene(MarkerGenome &marker_genome, size_t geneid);
    double VerticalCoverage();
    size_t ActiveBuckets();
    size_t TotalBuckets();
    double ActiveRate();
    double ExpectedActiveBuckets();

    void Add(IO::ClassificationLine &line);
    bool HasGene();

    double VerticalCoverage2();

    bool IsPresent();
};

class Statistics {
public:
    static double Median(std::vector<double> v) {
        if (v.empty()) return 0.0;
        auto m = v.begin() + v.size()/2;
        std::nth_element(v.begin(), m, v.end());
        double median = v.size()%2 == 0 ? (v[v.size()/2 - 1] + v[v.size()/2])/2 : v[v.size()/2];
        return median;
    }
    static double Mean(std::vector<double> v) {
        if (v.empty()) return 0.0;
        double sum = 0;
        for (auto& value : v) sum += value;
        return sum/v.size();
    }
};

struct Taxon {
    size_t id = UINT64_MAX;
    size_t mapped_reads = 0;
    size_t mapped_length = 0;

    // Test to calculate ani
    std::vector<uint32_t> hits;
    std::vector<uint32_t> overlap;
    std::vector<double> anis;

    MarkerGenome *marker_genome = nullptr;

    std::vector<Gene> marker_gene_buckets;

    double m_vertical_coverage = -1;

    bool leaf = false;

    Taxon() {};
    Taxon(size_t taxid) : id(taxid)  {
//        std::cout << "Taxon default constructor " << std::endl;
    };
    Taxon(size_t taxid, MarkerGenome &marker_genome) : id(taxid), marker_genome(&marker_genome), leaf(true) {
        marker_gene_buckets.resize(marker_genome.MaxMarkerGeneId()+1);
        for (int i = 0; i <= marker_genome.MaxMarkerGeneId(); i++) {
            if (marker_genome.Get(i) > 0) {
//                std::cout << "mglen: " << marker_genome.Get(i) << std::endl;
                if (i > marker_gene_buckets.size()) {
                    std::cout << " whoopu" << std::endl;
                    exit(44);
                }
                marker_gene_buckets[i] = Gene(marker_genome, i);
            }
        }
//        std::fill(marker_gene_buckets.begin(), marker_gene_buckets.end(), 0);
    };
    size_t MarkerGeneCount();
    size_t PresentGenes();
    double NaiveVerticalCoverage();
    std::string VerboseString();
    std::vector<double> GetBucketAbundances();
    void Add(IO::ClassificationLine &line);

    std::vector<double> GetBucketAbundances2();

    std::vector<double> GetBucketAbundances3(std::unordered_map<size_t, double> &map, double threshold=0.0);

    double VerticalCoverage();
};

class TaxonFilter {
    double allowed_deviation = 0;//0.4;//0.25 is good

    size_t mapped_reads_threshold = 31;
    double expected_genes_threshold = 0.8;
    double median_ani_threshold = 0.9;

public:
    TaxonFilter(size_t mapped_reads_threshold, double expected_genes_threshold, double median_ani_threshold) :
            mapped_reads_threshold(mapped_reads_threshold),
            expected_genes_threshold(expected_genes_threshold),
            median_ani_threshold(median_ani_threshold) {

    }

    bool Pass(Taxon &taxon) {
        if (!taxon.leaf)
            return false;
        if (taxon.MarkerGeneCount() == 0)
            return false;

        auto ratio = (double) taxon.PresentGenes() / ExpectedGenes(taxon.MarkerGeneCount(), taxon.mapped_reads);


        return ratio > expected_genes_threshold && taxon.mapped_reads > mapped_reads_threshold && Statistics::Median(taxon.anis) > median_ani_threshold;
//        return taxon.PresentGenes() >= (ExpectedGenes(taxon.MarkerGeneCount(), taxon.mapped_reads) * (1 - allowed_deviation));
    }

    static double ExpectedGenes(size_t marker_genes_count, size_t mapped_reads) {
//        std::cout << "Expected Calc: " << (1.f - pow((marker_genes_count - 1)/marker_genes_count, mapped_reads));
//        std::cout << "  Expected pow(): " << pow((double) (marker_genes_count - 1)/marker_genes_count, mapped_reads) << " mg count: " << marker_genes_count << " mapped reads: " << mapped_reads << "  " << (1 - pow((double) (marker_genes_count - 1)/marker_genes_count, mapped_reads)) * marker_genes_count << std::endl;
        return (1 - pow((double) (marker_genes_count - 1)/marker_genes_count, mapped_reads)) * marker_genes_count;
//        return (1 - pow( (marker_genes_count - 1)/marker_genes_count, mapped_reads)) * marker_genes_count;
    }

    static double ExpectedGenes(Taxon &taxon) {
//        std::cout << "ExpectedGenes: " << taxon.id << '\t' << taxon.MarkerGeneCount() << "(MG)\t" << taxon.mapped_reads << "(READS)\t" << ExpectedGenes(taxon.MarkerGeneCount(), taxon.mapped_reads) << " (Expected)" << std::endl;
        return ExpectedGenes(taxon.MarkerGeneCount(), taxon.mapped_reads);
    }
};


using TaxidSet = tsl::sparse_set<size_t>;
class AbundanceEstimator {
private:
    using MarkerGenomeMap = std::unordered_map<size_t, MarkerGenome>;

    MarkerGenomeMap m_marker_genomes;

    std::unordered_map<size_t, std::unordered_map<size_t, double>> taxon2abundance_factor;
    std::unordered_map<uint32_t, Taxon> taxa; // Mapping taxid to Taxon object
    std::unordered_map<size_t, double> taxon2depth; // Mapping taxid to Taxon object


    std::unordered_map<size_t, double> m_abundances;
    std::unordered_map<std::string, double> m_total_depths;

    TaxidSet m_predicted_strains;

    double total_depth = 0;

    // Filters for reads and taxa
    ReadFilter& read_filter;
    TaxonFilter& taxon_filter;

    // Internal taxonomy to translate classifications into real-world taxonomic ids
    Taxonomy::IntTaxonomy& taxonomy;

    // Add strain level functionality
    size_t m_current_sample_idx = 0;
    varkit::StdSampleHandler m_strain_handler;
    varkit::SampleStrains* m_active_sample;

//    tsl::sparse_set<size_t> m_valid_readids;


    using ReadUID = std::pair<size_t, uint8_t>;
    using ReadUIDSet = tsl::sparse_set<ReadUID, Utils::pair_hash>;
    ReadUIDSet m_valid_readuids;



    size_t passed_reads = 0;
    size_t declined_reads = 0;

    // Current active sample_name
    std::string sample_name = "";
    size_t sample_idx = 0;


    // Variables needed for formulas etc
    size_t k;
    size_t shape_length;

    bool correct_abundances = false;
    bool verbose = false;

    // If ani is determined at an earlier step and supposed to be more
    // accurate, then set this to false to not recompute ani.
    bool m_estimate_ani = true;


    // Verbose output stream
    std::ofstream* verbose_output = nullptr;
    // Raw taxa output outputs all possible taxa with their attributes
    std::ofstream* raw_taxa = nullptr;

    // This is for benchmarking if true information for each read is available
    BinaryClassifierEvaluator bce;




    // Cami output format functions
    std::string GetBioboxHeader(std::string sample_name) const;
    std::string BioboxString(size_t id, double abundance, bool as_leaf= false, std::vector<std::string> ranks= {
            "superkingdom", "phylum", "class", "order", "family", "genus", "species"}) const;

    // Include more output formats
    //TODO: Include more output formats for abundance profiles other than CAMI output format.


public:
    AbundanceEstimator(size_t k, size_t shape_length, ReadFilter &read_filter, TaxonFilter &taxon_filter, Taxonomy::IntTaxonomy &taxonomy);
    ~AbundanceEstimator();

    /**
     * Function to estimate ani from
     * @param hits number of k-mer hits
     * @param target_length length of sequence where hits occured
     * @param shape_size shape size used to obtain k-mer hits
     * @param k weight of k-mer shape (number of take positions or shape-size minus number of gaps)
     * @return double value serving as ANI estimate
     */
    static double EstimateAni(double hits, double target_length, size_t shape_size, size_t k) {
        if ((target_length - shape_size + 1) == 0) {
            std::cerr << "target_length: " << target_length <<  "  shape size: " << shape_size << std::endl;
            std::cerr << "(mean_read_length - shape_size + 1) is 0" << std::endl;
            exit(51);
        }
//        std::cout << ((double)hits / (target_length - shape_size + 1)) << "^(" << ((double)1/k) << ")" << std::endl;
        assert(pow((double) hits / (target_length - shape_size + 1), (double)1/k) >= 0);
        return pow((double) hits / (target_length - shape_size + 1), (double)1/k);
    }

    static double EstimateAni(double hits, double total_hits, size_t k) {
        return pow((double) hits / total_hits, (double)1/k);
    }

//    tsl::sparse_set<size_t>& GetValidReadIds() {
//        return m_valid_readids;
//    }

    ReadUIDSet& GetValidReadUids() {
        return m_valid_readuids;
    }

    bool AddLine(IO::ClassificationLine &line, bool force_pass=false, int true_internal_id=-1);
    void Profile(std::string output_file);
    void SetVerbose(bool verbose);
    void SetRawOutput(std::string path);
    void ReadClassificationFile(const char* file_path);
    void LoadMarkerGenome(const char* file_path);
    void Clear();
    ReadFilter& GetReadFilter();

    void WriteBiobox(const std::string path) const;

    void LoadAbundanceCorrection(const char *file_path);


    varkit::StdSampleHandler& GetSampleHandler() {
        return m_strain_handler;
    }

    TaxidSet& GetPredictedStrains() {
        return m_predicted_strains;
    }

    void DebugFunction();

    void SetSample(std::string sample_name);
    void SetSample(Sample &sample);

    void SetEstimateAni(bool value) {
        m_estimate_ani = value;
    }
};


#endif //VARKIT_ABUNDANCEESTIMATOR_H
