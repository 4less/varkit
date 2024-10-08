//
// Created by fritsche on 21/07/22.
//

#ifndef VARKIT_ALIGNER_H
#define VARKIT_ALIGNER_H

#include <vector>
#include "strain_handler.h"

struct AlignmentResult {
    struct SNP {
        size_t pos;
        char ref_base;
        char obs_base;
        SNP (size_t pos, char ref_base, char obs_base) :
                pos(pos), ref_base(ref_base), obs_base(obs_base) {};
    };
    using SNPList = std::vector<SNP>;
    double similarity = 0;
    size_t overlap = 0;
    SNPList snps;

};

class Aligner {
//    using Seed = std::pair<size_t, size_t>;
    struct Seed {
        const int64_t reference = -1;
        const int64_t query = -1;

        Seed(int64_t reference, int64_t query) : reference(reference), query(query) {};
    };
    static void SimpleAlignment(AlignmentResult &result, Seed const& seed, FastxRecord query, std::string& reference) {
        auto min = std::min(seed.reference, seed.query);
        auto rpos = seed.reference - min;
        auto qpos = seed.query - min;

        // Most basic, no alignment, no range checks, single seed.
        size_t total_bases = std::min(query.sequence.size() - qpos, reference.size() - rpos);
        for (; rpos < reference.size() && qpos < query.sequence.size(); rpos++, qpos++) {
            if (query.sequence[qpos] == reference[rpos]) {
                result.snps.emplace_back(AlignmentResult::SNP(rpos, reference[rpos], query.sequence[qpos]));
            }
        }
        result.snps.clear();
        result.overlap = total_bases;
        result.similarity = (double)result.snps.size()/total_bases;
    }
};


#endif //VARKIT_ALIGNER_H
