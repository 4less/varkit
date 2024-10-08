//
// Created by fritsche on 04/04/2022.
//

#include "HitsClassifier.h"

template< typename tPair >
struct second_t {
    typename tPair::second_type operator()( const tPair& p ) const { return p.second; }
};

template< typename tMap >
second_t< typename tMap::value_type > second( const tMap& m ) { return second_t< typename tMap::value_type >(); }



void HitsClassifier::Evaluate() {
//    std::vector<Hit> hits;
//    std::transform(m_hits.begin(), m_hits.end(), std::back_inserter(hits), second(m_hits));
//    std::sort(hits.begin(), hits.end(), [](Hit const& a, Hit const& b) { return std::max(a.hits_fwd, a.hits_rev) > std::max(b.hits_fwd, b.hits_rev); });
//
//    for (auto& hit : hits) {
//        std::cout << hit.ToString() << " <<< " << m_taxonomy.Get(hit.m_node_id).parent_id << std::endl;
//    }
//    std::cout << std::string(20, '-') << std::endl;
    EvaluateTree();
//    std::cout << std::string(20, '_') << std::endl;
}
