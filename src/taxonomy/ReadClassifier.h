//
// Created by fritsche on 06/07/2021.
//

#pragma once

#include <iostream>
#include "TaxonomyNew.h"
#include <fstream>
#include "ShapeUtils.h"
#include "data_structures.h"
#include "ShapeBuilder.h"


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

        size_t internal_taxid_bits;
        size_t gene_id_bits;
        size_t gene_pos_bits;

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


public:
    size_t internal_taxid_bits_;
    size_t gene_id_bits_;
    size_t gene_pos_bits_;
};


struct PatternDB {
    PatternDB(std::string path) {
        ifstream ifs(path, ifstream::in);

        std::string line;

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
        }
        ifs.close();
    }

    ~PatternDB() {
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

struct ReadClassifierResult {
    size_t taxid = 0;
    size_t best_count = 0;
    size_t leaf_hits = 0;
    size_t best_total_counts = 0;
    int gene_pos = -1;
    int pattern_idx = -1;
    int gene_id = -1;
    double rank_confidence = 0;
    double candidate_confidence = 0;
    double classification_confidence = 0;
    bool forward = true;
};


struct Mutation {
public:
    int read_pos_;
    int gene_pos_;

    friend bool operator<(const Mutation& l, const Mutation& r) {
        return l.gene_pos_ < r.gene_pos_;
    }

    friend bool operator>(const Mutation& l, const Mutation& r) {
        return l.gene_pos_ > r.gene_pos_;
    }

    Mutation(int read_pos, int gene_pos) : read_pos_(read_pos), gene_pos_(gene_pos) {}
};

class PatternProcessor2 {
    int current_pos_ = 0;
    bool sensitive_mode_ = false;

    static constexpr size_t INIT_HITMISS_SIZE = 100;
    static constexpr size_t INIT_PATTERNS_SIZE = 100;

//    std::vector<uint32_t> patterns_;
    std::vector<uint64_t> patterns_ { INIT_PATTERNS_SIZE } ;

    size_t pattern_len_;

    Taxonomy::IntTaxonomy& taxonomy_;
    std::set<Mutation> mutations;

    std::vector<uint32_t> cov;


public:
    int first_hit_ = INT_MAX;
    int last_hit_ = 0;
    int read_length_ = 0;
    int cov_bucket_size = -1;

    size_t last_size = 0;

    std::vector<ds::KmerHit> hitmiss_ {INIT_HITMISS_SIZE };

    const set<Mutation> &GetMutations() const {
        return mutations;
    }

    PatternProcessor2(Taxonomy::IntTaxonomy &taxonomy, std::string shape_path, bool sensitive_mode) :
            taxonomy_(taxonomy),
            sensitive_mode_(sensitive_mode) {

        hitmiss_.resize(100);
        patterns_.resize(100);
        cov.resize(10);

        if (!Utils::exists(shape_path)) {
            std::cerr << "File " << shape_path << " does not exist. Abort." << std::endl;
            exit(7);
        }
        snp_detector_.LoadArray(shape_path);

        pattern_len_ = snp_detector_.pattern_size;
    }

    PatternProcessor2(Taxonomy::IntTaxonomy &taxonomy, PatternDB* db, bool sensitive_mode) :
            taxonomy_(taxonomy),
            sensitive_mode_(sensitive_mode) {

        hitmiss_.resize(100);
        patterns_.resize(100);
        cov.resize(10);

        pattern_db = db;
        pattern_len_ = db->pattern_size;
    }

    void Reset() {
        patterns_[0] = 0;
        current_pos_ = 0;
        mutations.clear();
        first_hit_ = INT_MAX;
        last_hit_ = 0;
        std::fill(cov.begin(), cov.end(), 0);
    }

    inline void AddNode(int taxid) {
        last_size = hitmiss_.size();

        hitmiss_[current_pos_].Reset();
        hitmiss_[current_pos_++].taxid = taxid;


        if (current_pos_ == hitmiss_.size()) {
            hitmiss_.resize(current_pos_ * 2);
//            std::cout << "len: " << hitmiss_.size() << std::endl;
        }
    }

    inline void AddLeaf(int taxid, int geneid, int genepos) {
        last_size = hitmiss_.size();

        std::cout << "current_pos_: " << current_pos_ << std::endl;
        if (current_pos_ >= hitmiss_.size()) {
            std::cerr << "current_pos_ >= hitmiss_.size()" << std::endl;
        }

        hitmiss_[current_pos_].taxid = taxid;
        hitmiss_[current_pos_].geneid = geneid;
        hitmiss_[current_pos_++].genepos = genepos;


        if (current_pos_ == hitmiss_.size()) {
            std::cout << "RESIZE current_pos_: " << current_pos_ << std::endl;
            hitmiss_.resize(current_pos_ * 2);
        }
    }

    inline void AddMiss() {
        hitmiss_[current_pos_].Reset();
        current_pos_++;


        if (current_pos_ == hitmiss_.size()) {
            hitmiss_.resize(current_pos_ * 2);
        }
    }

    inline void Add(size_t pos, bool hit) {
        if (pos/cov_bucket_size > cov.size()) {
            cov.resize(pos/cov_bucket_size);
        }
        cov[pos / cov_bucket_size] += hit;
        if (pos < pattern_len_) {
            patterns_[0] |= ((uint64_t) hit) << (63 - pos);
        } else {
            patterns_[pos - pattern_len_ + 1] = (patterns_[pos - pattern_len_] << 1) | ((uint64_t) hit << (64 - pattern_len_));
        }
    }

    inline double GetHitCov() {
        int hit = 0;

        for (int i = 0; i < 10; i++) {
            hit += cov[i];
        }

        return (double) hit / (read_length_ - snp_detector_.shape_size + 1);
//        return (double) hit / (read_length_ - pattern_db->shape_size + 1);
    }


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


//    int Evaluate3(classification::Hit& hit, size_t sequence_length, size_t shape_size) {
//        auto &&[offset, hits, forward] = hit.GetOffset().GetHits(shape_size, sequence_length);
////        OverlappingKmers()
//        return 0;
//    }

    int Evaluate2(int taxid, int geneid, int abs_pos, bool forward, size_t shape_length, size_t read_length, int start, int n) {
        std::string lookup_str = "";
        int rel_pos = forward ? abs_pos : abs_pos + read_length - shape_length;
        for (int pos = 0; pos < current_pos_; pos++) {
            // Iterate all k-mer lookups
            auto& lookup = hitmiss_[pos];

            // Does that work?
            if (lookup.IsMiss()) {
                lookup_str += '0';
                continue;
            }

            if (hitmiss_.size() > 200) {
                std::cout << "asdsad" << std::endl;
                exit(9);
            }

//            std::cout << " (lookup.genepos + ((forward * -1) + (!forward * 1)) * pos)): " << "rel: " << rel_pos << " " << (lookup.genepos + ((forward * -1) + (!forward * 1)) * pos) << std::endl;

            bool leaf_match = (lookup.IsGeneSet() && lookup.geneid == geneid && lookup.taxid == taxid && (rel_pos == (lookup.genepos + ((forward * -1) + (!forward * 1)) * pos)));

            if (hitmiss_.size() > 200) std::cout << "2. loop: " << pos << " size: " << hitmiss_.size() << std::endl;

            bool hit = lookup.taxid >= 0 && // First requirement taxid must not be -1
                       (leaf_match || // Second, geneid and genepos may not be 0
                        (!lookup.IsGeneSet() && (taxid == lookup.taxid || taxonomy_.IsNodeAncestor(lookup.taxid, taxid))));


            if (hit) {
                lookup_str += '1';
            } else {
                lookup_str += 'O';
            }

            if (leaf_match || (sensitive_mode_ && hit)) {
                if (first_hit_ == INT_MAX) first_hit_ = pos;
                last_hit_ = pos - pattern_len_ + 1;
            }

            Add(pos , hit);
        }

        std::cout << "full:  " << lookup_str << std::endl;
        std::cout << "short: " << lookup_str.substr(start, n) << std::endl;

        for (int i = first_hit_; i <= last_hit_; i++) {
            if (!patterns_[i]) continue;
            uint32_t pattern = patterns_[i] >> 32;

            auto& mutation = pattern_db->map_array[pattern >> (32llu - pattern_db->pattern_size)];

            for (int m = 0; m < mutation.Size(); m++) {
                int mutation_pos = mutation.mutations[m] + i;
                int gene_pos = forward ? rel_pos + mutation_pos : rel_pos - mutation_pos + shape_length - 1;

                auto mut = Mutation{ mutation_pos, gene_pos };
                mutations.insert(mut);
            }
        }

        return mutations.size();
    }

    int Evaluate(int taxid, int geneid, int rel_pos, bool forward, int shape_length) {
        std::string lookup_str = "";
        for (int pos = 0; pos < current_pos_; pos++) {
            // Iterate all k-mer lookups
            auto& lookup = hitmiss_[pos];

            // Does that work?
            if (lookup.IsMiss()) {
                lookup_str += '0';
                continue;
            }

            if (hitmiss_.size() > 200) {
                std::cout << "asdsad" << std::endl;
                exit(9);
            }

//            std::cout << " (lookup.genepos + ((forward * -1) + (!forward * 1)) * pos)): " << "rel: " << rel_pos << " " << (lookup.genepos + ((forward * -1) + (!forward * 1)) * pos) << std::endl;

            bool leaf_match = (lookup.IsGeneSet() && lookup.geneid == geneid && lookup.taxid == taxid && (rel_pos == (lookup.genepos + ((forward * -1) + (!forward * 1)) * pos)));

            if (hitmiss_.size() > 200) std::cout << "2. loop: " << pos << " size: " << hitmiss_.size() << std::endl;

            bool hit = lookup.taxid >= 0 && // First requirement taxid must not be -1
                    (leaf_match || // Second, geneid and genepos may not be 0
                     (!lookup.IsGeneSet() && (taxid == lookup.taxid || taxonomy_.IsNodeAncestor(lookup.taxid, taxid))));


            if (hit) {
                lookup_str += '1';
            } else {
                lookup_str += 'O';
            }

            if (leaf_match || (sensitive_mode_ && hit)) {
                if (first_hit_ == INT_MAX) first_hit_ = pos;
                last_hit_ = pos - pattern_len_ + 1;
            }

            Add(pos , hit);
        }

        std::cout << lookup_str << std::endl;

        for (int i = first_hit_; i <= last_hit_; i++) {
            if (!patterns_[i]) continue;
            uint32_t pattern = patterns_[i] >> 32;

            auto& mutation = pattern_db->map_array[pattern >> (32llu - pattern_db->pattern_size)];

            for (int m = 0; m < mutation.Size(); m++) {
                int mutation_pos = mutation.mutations[m] + i;
                int gene_pos = forward ? rel_pos + mutation_pos : rel_pos - mutation_pos + shape_length - 1;

                auto mut = Mutation{ mutation_pos, gene_pos };
                mutations.insert(mut);
            }
        }

        return mutations.size();
    }

    PatternDB* pattern_db = nullptr;
    ShapeBuilder snp_detector_;
};

struct TreeNode {
    int id = -1, parent_id = -1, count = 0, cumulated_count = 0;
    bool visited = false;

    void ResetCounts() {
        count = 0;
        cumulated_count = 0;
        visited = false;
    }

    bool IsReset() {
        return count == 0 && cumulated_count == 0 && !visited;
    }
    std::string ToString() {
        std::string result = "";
        result += "id: " + std::to_string(id) + '\n';
        result += "pid: " + std::to_string(parent_id) + '\n';
        result += "count: " + std::to_string(count) + '\n';
        result += "cumulated: " + std::to_string(cumulated_count) + '\n';
        result += "visited: " + std::to_string(visited) + '\n';
        return result;
    }

    ~TreeNode() {
//        std::cout << id << ": WTF_________________________________________________________________" << std::endl;
    }
};

struct LeafNode {
    int taxid = -1;
    int geneid = -1;
    int readpos_a = -1;
    int readpos_b = -1;
    int counts = 0;
    int counts_a = 0;
    int counts_b = 0;
    int leaf_hits = 0;
    int cumulated = 0;
    int pattern_idx = -1;

    std::vector<bool> pattern;
    size_t pattern_size = 0;

    std::string ToString() {
        std::string result = "";
        result += "taxid: " + to_string(taxid);
        result += "  geneid: " + to_string(geneid);
        result += "  rpos_a: " + to_string(readpos_a);
        result += "  rpos_b: " + to_string(readpos_b);
        result += "  counts: " + to_string(counts);
        result += "  rcnt_a: " + to_string(counts_a);
        result += "  rcnt_b: " + to_string(counts_b);
        result += "  cumulated: " + to_string(cumulated);
        return result;
    }

    LeafNode(){
        pattern = std::vector<bool> (200, false);
    };

    bool Equals(int _taxid, int _geneid, int _genepos, int _readpos) {
        return
                taxid == _taxid &&
                geneid == _geneid &&
                (readpos_a == _genepos - _readpos ||
                 readpos_b == _genepos + _readpos);
    }

    bool Supports(Taxonomy::IntTaxonomy& taxonomy, int _taxid, int _geneid, int _genepos, int _readpos) {
        return
                taxonomy.IsNodeAncestor(_taxid ,taxid) &&
                geneid == _geneid &&
                (readpos_a == _genepos - _readpos ||
                 readpos_b == _genepos + _readpos);
    }

    void Increment(int _genepos, int _readpos) {
        if (readpos_a == _genepos - _readpos) {
            counts_a++;
        } else if (readpos_b == _genepos + _readpos) {
            counts_b++;
        }
        counts = std::max(counts_a, counts_b);
        pattern[_readpos] = true;
    }

    void IncrementLeafs() {
        leaf_hits++;
    }

    int GetLeafHits() {
        return leaf_hits;
    }

    void Set(int _taxid, int _geneid, int _genepos, int _readpos) {
        this->taxid = _taxid;
        this->geneid = _geneid;
        this->readpos_a = _genepos - _readpos;
        this->readpos_b = _genepos + _readpos;

//        if (this->readpos_a == this->readpos_b) {
//            std::cout << "Detect_____" << std::endl;
//            std::cout << "rposa: " << this->readpos_a << std::endl;
//            std::cout << "rposb: " << this->readpos_b << std::endl;
//            std::cout << "_geneid: " << _geneid << std::endl;
//            std::cout << "_genepos: " << _genepos << std::endl;
//            std::cout << "_readpos: " << _readpos << std::endl;
//        }

        this->counts_a = 1;
        this->counts_b = 1;
        this->leaf_hits = 1;
        this->counts = 1;
        this->cumulated = 0;
        this->pattern_idx = -1;

        // Apparently with -O3 and gcc std::fill will use memset.
        // So it wont iterate over the elements to set them false
        pattern_size = _readpos;
        std::fill(pattern.begin(), pattern.end(), false);
    }
};

class ReadClassifierMM {
    std::vector<LeafNode> leaves_;
    size_t leaves_size_ = 0;

    size_t read_length_ = 0;
    size_t read_position_ = 0;
    ValueExtractorMG extractor_;

    int max_leaf_ = -1;
    int max_count_ = 0;


public:
    ReadClassifierMM(ValueExtractorMG& extractor, int shape) :
        extractor_(extractor) {
        leaves_ = std::vector<LeafNode>(30000, LeafNode());
    }

//    void InitProcessor(size_t pattern_len, std::string snp_finder, bool sensitive_mode) {
//        pattern_processor_ = make_unique<PatternProcessor2>(snp_finder, sensitive_mode);
//    }

    void Reset() {
        read_position_ = 0;
        leaves_size_ = 0;
        max_leaf_ = -1;
        max_count_ = 0;
    }

    void SetLength(size_t read_length) {
        read_length_ = read_length;
    }

    void AddHit(size_t value) {
        auto gene_pos = extractor_.GetGenePos(value);
        auto taxid = extractor_.GetTaxid(value);
        auto gene_id = extractor_.GetGeneId(value);

        auto i = 0;
        bool found = false;
        for (i = 0; i < leaves_size_; i++) {
            if (leaves_[i].Equals(taxid, gene_id, gene_pos, read_position_)) {
                found = true;
                leaves_[i].Increment(gene_pos, read_position_);
                leaves_[i].IncrementLeafs();
            }
        }
        if (!found) {
            leaves_[i].Set(taxid, gene_id, gene_pos, read_position_);
            leaves_size_++;
        }

        // Find the highest hit fast
        if (leaves_[i].counts > max_count_) {
            max_count_ = leaves_[i].counts;
            max_leaf_ = i;
        }

        if (leaves_size_ >= leaves_.size()) {
            std::cout << "resize: " << leaves_.size() * 2 << std::endl;
            leaves_.resize(leaves_.size() * 2);
        }
    }

    size_t GetLeavesSize() {
        return leaves_size_;
    }

    size_t GetMaxCount() {
        return max_count_;
    }

    void PrintLeafs() {
        if (leaves_size_ == 0)
            return;
        std::cout << ">" << leaves_size_ << "---------------------------------------------" << std::endl;
        for (auto i = 0; i < leaves_size_; i++) {
            auto& leaf = leaves_[i];

            std::cout << leaf.ToString() << std::endl;
        }

    }

    void Evaluate() {

    }

    void Move() {
        read_position_++;
    }
};



class ReadClassifier {
    // This part is for non leaf counts
    std::vector<TreeNode> nodes;
    size_t size;

    std::vector<int> hits;
    std::vector<int> clear;

//    int hits[1000];
//    int clear[1000];
    int hits_size_ = 0;
    int clear_size_ = 0;

    unique_ptr<PatternProcessor2> processor = nullptr;

private:
    int shape_size_;
    int read_length_ = 0;

    // This part is for leaf counts (include gene_id, and gene_pos)
    inline static const size_t leaves_capacity = 1000;
    size_t leaves_size_ = 0;
    LeafNode leaves[leaves_capacity];

    int max_index_ = -1;
    size_t max_count_ = 0;

    ValueExtractorMG extractor_;
    Taxonomy::IntTaxonomy taxonomy_;


    void Init() {
        int max = 0;
        for (auto iterator = taxonomy_.map.begin(); iterator != taxonomy_.map.end(); iterator++) {
            if (iterator->first > max) max = iterator->first;
        }

        size = max + 1;
        nodes.resize(size);
//        nodes = new TreeNode[size];

        for (auto iterator = taxonomy_.map.begin(); iterator != taxonomy_.map.end(); iterator++) {
            nodes[iterator->first].id = iterator->second.id;
            nodes[iterator->first].parent_id = iterator->second.parent_id;
        }
        assert(AreNodesReset());
    }


    int AddGenePosSupport(int taxid, int geneid, int genepos, int readpos) {
        // Check all previous counts if one matches (linear search)
        for (int i = 0; i < leaves_size_; i++) {
            if (leaves[i].Supports(taxonomy_, taxid, geneid, genepos, readpos)) {
                leaves[i].Increment(genepos, readpos);
                if (leaves[i].counts > max_count_) {
                    max_count_ = leaves[i].counts;
                    max_index_ = i;
                }
                return leaves[i].pattern_idx;
            }
        }
        return -1;
    }

    int AddLeaf(int taxid, int geneid, int genepos, int readpos) {

        // Check if it was the last hit
        if (leaves_size_ && leaves[max_index_].Equals(taxid, geneid, genepos, readpos)) {
            leaves[max_index_].Increment(genepos, readpos);
            max_count_++;
            return leaves[max_index_].pattern_idx;
        }

        // Check all previous counts if one matches (linear search)
        for (int i = 0; i < leaves_size_; i++) {
            if (leaves[i].Equals(taxid, geneid, genepos, readpos)) {
                leaves[i].Increment(genepos, readpos);
                if (leaves[i].counts > max_count_) {
                    max_count_ = leaves[i].counts;
                    max_index_ = i;
                }
                return leaves[i].pattern_idx;
            }
        }

        // new taxid + gene combination.
        auto& leaf = leaves[leaves_size_];

        leaf.Set(taxid, geneid, genepos, readpos);
        if (leaf.counts > max_count_) {
            max_count_ = leaf.counts;
            max_index_ = leaves_size_;
        }
        leaves_size_++;

        if (leaves_size_ >= leaves_capacity) {
            std::cout << "leaves_size exceeded " << leaves_size_ << " >= " << leaves_capacity << std::endl;
            exit(9);
        }

        return leaf.pattern_idx;
    }

    LeafNode& GetLastAddedLeaf() {
        return leaves[leaves_size_ - 1];
    }

    void AddNonLeaf(int taxid) {
        if (!nodes[taxid].count) {
            // New inner node
            if (hits_size_ >= hits.size())
                hits.resize((hits.size() + 1) * 2);
            hits[hits_size_++] = taxid;
        }
        nodes[taxid].count++;
    }

public:
    size_t pos_ = 0;


    ReadClassifier(ValueExtractorMG &extractor, Taxonomy::IntTaxonomy &taxonomy, int shape_size) :
            extractor_(extractor),
            taxonomy_(taxonomy),
            shape_size_(shape_size) {
        Init();
    }

    ~ReadClassifier() {
//        delete[] nodes;
    }

    void InitProcessor(std::string snp_finder, bool sensitive_mode) {
        processor = make_unique<PatternProcessor2>(taxonomy_, snp_finder, sensitive_mode);
    }

    void InitProcessor(PatternDB* pattern_db, bool sensitive_mode) {
        processor = make_unique<PatternProcessor2>(taxonomy_, pattern_db, sensitive_mode);
        assert(processor);
    }

    void PrintLeafsP(int index = -1) {
        if (index != -1) {
            auto leaf = leaves[index];
            std::cout << index << " taxon: " << taxonomy_.Get(leaf.taxid).scientific_name << " geneid: " << leaf.geneid << " readpos_a: " << leaf.readpos_a << " readpos_b: " << leaf.readpos_b << " COUNT: " << leaf.counts << " ACCUM: " << leaf.cumulated << " counts a/b " << leaf.counts_a << "/" << leaf.counts_b << std::endl;
            return;
        }

        for (int i = 0; i < leaves_size_; i++) {
            auto leaf = leaves[i];
            std::cout << i << " taxon: " << taxonomy_.Get(leaf.taxid).scientific_name << " geneid: " << leaf.geneid << " readpos_a: " << leaf.readpos_a << " readpos_b: " << leaf.readpos_b << " COUNT: " << leaf.counts << " ACCUM: " << leaf.cumulated << std::endl;
        }
    }

    void AddMiss() {
        assert(processor);
        processor->AddMiss();
        pos_++;
    }

    void AddHit(size_t value) {
        // Add part

        auto taxid = extractor_.GetTaxid(value);
        auto geneid = extractor_.GetGeneId(value);
        auto genepos = extractor_.GetGenePos(value);

//        if (taxid == 29202) {
//            printf("taxid: %lu  geneid: %lu    genepos: %lu \n", taxid, geneid, genepos);
//        }

        if (taxonomy_.IsLeaf(taxid) && (geneid || genepos)) {
            AddLeaf(taxid, geneid, genepos, pos_);
            processor->AddLeaf(taxid, geneid, genepos);
        } else if (geneid) {
//            std::cout << "Add gene pos support: " << taxid << " " << geneid << std::endl;
            if (AddGenePosSupport(taxid, geneid, genepos, pos_) == -1) {
                AddNonLeaf(taxid);
                processor->AddNode(taxid);
            }
        } else {
            AddNonLeaf(taxid);
            processor->AddNode(taxid);
        }
        pos_++;
    }

    bool AreNodesReset() {
        for (int i = 0; i < nodes.size(); i++) {
            if (!nodes[i].IsReset()) {
                std::cout << "culprit: " << nodes[i].ToString() << std::endl;
                return false;
            }
        }
        return true;
    }

    int EvaluatePatterns(int index) {
        // index == -1 means no hit on species level
        if (index < 0) return 0;

        auto& leaf = leaves[index];

//        std::cout << "print leaf: " << leaf.ToString() << std::endl;

        bool forward = leaf.counts_a > leaf.counts_b;
        auto rel_pos = 0;
        if (forward) {
            rel_pos = leaf.readpos_a;
        } else {
            rel_pos = leaf.readpos_b;
        }

        return processor->Evaluate(leaf.taxid, leaf.geneid, rel_pos, forward, shape_size_);
    }

    int AddUp(TreeNode &node) {
        static int addup_counter = 0;
        if (node.cumulated_count || node.visited) return node.cumulated_count;

        if (clear_size_ >= clear.size())
            clear.resize((clear.size() + 1) * 2);
        clear[clear_size_++] = node.id;
        if (node.id == node.parent_id) return node.count;

        node.cumulated_count = AddUp(nodes[node.parent_id]) + node.count;
        node.visited = true;
        return node.cumulated_count;
    }

    void EvaluateTreeReturnBest(size_t &taxid, int &gene_id, int &gene_pos, bool &forward, size_t &count, size_t &leaf_hits, size_t &supportive_counts, int &pattern_idx) {
        int max_counts = 0, max_id = -1, max_idx = -1;

        // First go over leafs.
        for (int i = 0; i < leaves_size_; i++) {
            // Get leaf hit (points to genome in index)
            auto& leaf = leaves[i];

            // Get parent node, and calculate the cumulated hits for that node
            // If parent has already been visited AddUp just returns the hits
            auto& parent = nodes[taxonomy_.Get(leaf.taxid).parent_id];
            auto node_counts = leaf.counts + AddUp(parent);


            leaf.cumulated = node_counts;

            if (node_counts > max_counts) {
                max_counts = node_counts;
                max_idx = i;
                max_id = leaf.taxid;
            }
        }

        if (max_counts) {
            if (clear_size_ >= clear.size())
                clear.resize((clear.size() + 1) * 2);

            clear[clear_size_++] = max_id;
            nodes[max_id].cumulated_count = leaves[max_idx].cumulated;
            nodes[max_id].count = leaves[max_idx].counts;
        }

        // Second go over internal nodes.
        for (int i = 0; i < hits_size_; i++) {
            auto& node = nodes[hits[i]];
            auto node_counts = AddUp(node);
            if (node_counts > max_counts) {
                max_counts = node_counts;
                max_idx = -1;
                max_id = node.id;
            }
        }

        if (max_counts) {
            taxid = max_id;

            count = nodes[max_id].count;
            supportive_counts = nodes[max_id].cumulated_count;

            if (max_idx > -1 && max_id > 0 && taxonomy_.IsLeaf(max_id)) {
                auto& leaf_max = leaves[max_idx];
                leaf_hits = leaves[max_idx].GetLeafHits();

                forward = leaf_max.counts_a > leaf_max.counts_b;
                gene_id = leaf_max.geneid;

                gene_pos = leaf_max.counts_a > leaf_max.counts_b ? leaf_max.readpos_a : leaf_max.readpos_b;

                if (!forward) gene_pos += shape_size_;

                pattern_idx = max_idx;

            }
        }
    }


    void Evaluate(ReadClassifierResult &result) {
        int max_counts = 0, max_id = -1, max_idx = -1;

        // First go over leafs.
        for (int i = 0; i < leaves_size_; i++) {
            // Get leaf hit (points to genome in index)
            auto& leaf = leaves[i];

            // Get parent node, and calculate the cumulated hits for that node
            // If parent has already been visited AddUp just returns the hits
            auto& parent = nodes[taxonomy_.Get(leaf.taxid).parent_id];
            auto node_counts = leaf.counts + AddUp(parent);


            leaf.cumulated = node_counts;

            if (node_counts > max_counts) {
                max_counts = node_counts;
                max_idx = i;
                max_id = leaf.taxid;
            }
        }

        if (max_counts) {
            if (clear_size_ >= clear.size())
                clear.resize((clear.size() + 1) * 2);

            clear[clear_size_++] = max_id;
            nodes[max_id].cumulated_count = leaves[max_idx].cumulated;
            nodes[max_id].count = leaves[max_idx].counts;
        }

        // Second go over internal nodes.
        for (int i = 0; i < hits_size_; i++) {
            auto& node = nodes[hits[i]];
            auto node_counts = AddUp(node);
            if (node_counts > max_counts) {
                max_counts = node_counts;
                max_idx = -1;
                max_id = node.id;
            }
        }

        if (max_counts) {
            result.taxid = max_id;
            result.best_count = nodes[max_id].count;
            result.best_total_counts = nodes[max_id].cumulated_count;

            if (max_idx > -1 && max_id > 0 && taxonomy_.IsLeaf(max_id)) {
                auto& leaf_max = leaves[max_idx];
                result.leaf_hits = leaves[max_idx].GetLeafHits();

                result.forward = leaf_max.counts_a > leaf_max.counts_b;
                result.gene_id = leaf_max.geneid;
                result.gene_pos = leaf_max.counts_a > leaf_max.counts_b ? leaf_max.readpos_a : leaf_max.readpos_b;
                if (!result.forward) result.gene_pos += shape_size_;
                result.pattern_idx = max_idx;
            }
        }
    }

    pair<int, float> GetBestIdForConfidence(int id) {
        auto node = nodes[id];

        auto total_count = node.cumulated_count;
        auto node_count = node.count;

        return { id, node_count / total_count };
    }


    void Reset() {

        hits_size_ = 0;
        for (int i = 0; i < clear_size_; i++) {
            nodes[clear[i]].ResetCounts();
        }

        clear_size_ = 0;

        pos_ = 0;
        max_index_ = -1;
        max_count_ = 0;
        leaves_size_ = 0;

//        assert(AreNodesReset());

        processor->Reset();
    }

    const unique_ptr<PatternProcessor2> &GetProcessor() const {
        return processor;
    }

    void PrintLeafs() {
        for (int i = 0; i < leaves_size_; i++) {
            auto leaf = leaves[i];
            std::cout << leaf.ToString() << std::endl;
//            std::cout << i << " " << leaf.taxid << " " << leaf.geneid << " " << leaf.readpos_a << " " << leaf.readpos_b << " countsa: " << leaf.counts_a << " countsb: " << leaf.counts_b << std::endl;
        }
    }

    void SetLength(unsigned long i) {
        read_length_ = i;
        processor->read_length_ = i;
        processor->cov_bucket_size = (i - processor->snp_detector_.shape_size + 1) / 10;
//        processor->cov_bucket_size = (i - processor->pattern_db->shape_size + 1) / 10;
    }
};
