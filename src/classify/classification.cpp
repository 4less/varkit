//
// Created by fritsche on 25/04/2022.
//

#include "classification.h"

namespace classification {

    Classifier::Classifier(Taxonomy &taxonomy) :
            m_taxonomy(taxonomy),
            m_tree(m_taxonomy.MaxLeafId()) {
        LoadLineage();
    }

    void Classifier::LoadLineage() {
        m_lineages = new Lineage[m_taxonomy.MaxTaxid() + 1];
        for (const auto& [id, node] : m_taxonomy.map) {
            auto cid = node.id;

            m_lineages[id].SetRankId(Taxonomy::Rank2Id(node.rank));
            if (Taxonomy::Rank2Id(node.rank) == -1 && !m_taxonomy.IsRoot(node.id)) {
                std::cout << node.ToString() << std::endl;
                exit(2);
            }
            while (!m_taxonomy.IsRoot(cid)) {
                auto n = m_taxonomy.Get(cid);
                m_lineages[id].Set(Taxonomy::Rank2Id(n.rank), n.id);
                cid = n.parent_id;
            }
        }
    }

    bool Classifier::HasNode(UNID node_id) {
        return m_tree.HasNode(node_id);
    }

    void Classifier::ListNodes() {
        std::cout << "Size: " << m_tree.Data().size() << std::endl;
//        std::cout << m_tree.ToString() << std::endl;

        for (auto& [id, node] : m_tree.Data()) {
            assert(id == node.GetNodeId());
            std::cout << NodeId::ToString(id) << ", Node: ";
            std::cout << NodeId::ToString(node.GetNodeId()) << "\tParent: ";
            std::cout << NodeId::ToString(node.GetParentNodeId()) << "    ";

            std::cout << node.ToString() << std::endl;
        }
    }

    void Classifier::PrintTree(C1Node& node, unsigned int level= 0) {
        assert(m_tree.HasRoot());

        if (level > 15) return;
        std::cout << std::string(level, '\t') << NodeId::ToString(node.GetNodeId()) << '\t' << node.ToString() << std::endl;
        for (auto& child_id : node.Children())
            PrintTree(m_tree.GetNode(child_id), level + 1);
    }

    void Classifier::PrintTree() {
        if (m_tree.HasRoot())
            PrintTree(m_tree.Root());
    }

    void Classifier::Clear() {
        m_tree.Clear();
        m_consistent_strains = true;
        if (processor) processor->Reset();
    }

    size_t Classifier::Dummy() {
        size_t ret = 0;
        for (auto& [key, node] : m_tree.Data()) {
            ret += node.Data().GetMaxHits();
        }
        return ret;
    }

    std::string Classifier::Newick() {
        assert(m_tree.HasRoot());
        return Newick(m_tree.Root(), 0);
    }

    std::string Classifier::Newick(C1Node& node, unsigned int level) {
        std::string result = "";
        // Recursion
        if (!node.Children().empty()) {
            result += '(';
            result += Newick(m_tree.GetNode(node.Children()[0]), level+1);
            for (auto i = 1; i < node.Children().size(); i++) {
                result += ',';
                result += Newick(m_tree.GetNode(node.Children()[i]), level+1);
            }
            result += ')';
        }
        // Append current
        auto taxonomic_id = NodeId::GetTaxonomicId(node.GetNodeId());
        auto gene_id = NodeId::GetGeneId(node.GetNodeId());
        std::string node_label = std::to_string(taxonomic_id);
        if (gene_id) node_label += "_" + std::to_string(gene_id);
        node_label += "_hits=";
        node_label += node.ToString();
//        node_label += ']';
        result += node_label;
        return result;
    }


    Classifier::~Classifier() {
        delete[] m_lineages;
    }

    C1Tree &Classifier::Tree() {
        return m_tree;
    }


    void Classifier::AddHit(TID taxonomic_id, GID gene_id, GPOS gene_position, GPOS read_position) {
        if (m_tree.IsLeaf(taxonomic_id, gene_id, gene_position)) {
            // Hit is a tree leaf
            auto leaf_id = NodeId::ToNodeId(taxonomic_id, gene_id);
            if (m_tree.HasLeaf(leaf_id)) {
                auto& leaf = m_tree.GetLeaf(leaf_id);
                leaf.Data().Add(gene_position, read_position);
            } else {
                AddLeaf(leaf_id, taxonomic_id, gene_id, gene_position, read_position);
            }
            processor->AddLeaf(taxonomic_id, gene_id, gene_position);
        } else {
            // Hit is an inner node
            auto node_id = NodeId::ToNodeId(taxonomic_id);
            if (m_tree.HasNode(node_id)) {
                auto& node = m_tree.GetNode(node_id);
                node.Data().Increment();
            } else {
                AddNode(node_id, taxonomic_id, 0);
            }
            processor->AddNode(taxonomic_id);
        }
    }


    void Classifier::AddMiss() {
        assert(processor);
        processor->AddMiss();
    }

    void Classifier::AddLeaf(UNID id, TID taxonomic_id, GID gene_id, GPOS gene_pos, GPOS read_pos) {
        auto& lineage = m_lineages[taxonomic_id];

        auto leaf_node_id = NodeId::ToNodeId(taxonomic_id, gene_id);
        auto parent_node_id = NodeId::ToNodeId(taxonomic_id);
        auto leaf = C1Node(leaf_node_id, parent_node_id);

        leaf.Data().Set(gene_pos, read_pos);
        m_tree.Data().insert( { leaf_node_id, leaf } );

        if (!HasNode(parent_node_id)) {
            // Elicit Cascade to insert all non-existing nodes
            AddNode(parent_node_id, taxonomic_id, 0, lineage);
        }

        // Add new leaf to direct parent
        m_tree.GetNode(parent_node_id).AddChild(leaf_node_id);
    }

    void Classifier::AddNode(UNID id, TID taxonomic_id, GID gene_id) {
        // Node is new. Insert and increment, then make sure
        // Every node in the lineage is present as well.
        auto& lineage = m_lineages[taxonomic_id];
//        std::cout << "Add: " << lineage.ToString() << std::endl;
        AddNode(id, taxonomic_id, gene_id, lineage, 1);
    }

    void Classifier::AddNode(UNID id, TID taxonomic_id, GID gene_id, Lineage& lineage, size_t init_hits) {

        // In a line of taxa: a b c d e with the current being e
        // child_iterator points to e
        // current iterator points to d
        // parent iterator points to c

//        std::cout << "Create node cascade for " << taxonomic_id << " (" << NodeId::ToString(id) << ")" << std::endl;

        // Current iterator
        auto current_iterator = lineage.rbegin() - !gene_id; // d
        while (*current_iterator == -1) current_iterator--;
        if (current_iterator < lineage.begin()) current_iterator = &m_tree.RootId();

        // Child iterator
        auto child_node_id = id;

        auto current_node_id = NodeId::ToNodeId(*current_iterator);

        auto node = C1Node(child_node_id, current_node_id);
        if(init_hits)
            node.Data().Increment();
        m_tree.Insert(child_node_id, node);

        auto parent_iterator = current_iterator;
        while (*(--parent_iterator) == -1);
        if (parent_iterator < lineage.begin()) parent_iterator = &m_tree.RootId();

        auto parent_node_id = NodeId::ToNodeId(*parent_iterator);

//        std::cout << "while (!HasNode(" << NodeId::ToString(current_node_id) << ") { ..." << std::endl;

        while (!HasNode(current_node_id)) {
            // Insert current node
            m_tree.Insert(current_node_id, C1Node(current_node_id, parent_node_id));
            m_tree.GetNode(current_node_id).AddChild(child_node_id);

//            std::cout << "\tInsert: " << NodeId::ToString(current_node_id) << std::endl;

            child_node_id = current_node_id;

            current_iterator = parent_iterator;
            current_node_id = parent_node_id;

            while (*(--parent_iterator) == -1);
            if (parent_iterator < lineage.begin()) parent_iterator = &m_tree.RootId();
            if (*current_iterator == m_tree.RootId()) parent_iterator = &m_tree.RootId();
            parent_node_id = NodeId::ToNodeId(*parent_iterator);
        }
//        std::cout << "}" << std::endl;

        if (current_node_id != child_node_id)
            m_tree.GetNode(current_node_id).AddChild(child_node_id);
    }

    bool Classifier::IsAmbiguous() {
        for (auto& [id, node] : m_tree.Data()) {
            if (node.Data().HasAlternativeOffsets()) return true;
        }
        return false;
    }

    void Classifier::PrintAmbiguous() {
        const auto& map = m_tree.Data();
        for (auto&& [ id, node] : m_tree.Data()) {
            if (node.Data().HasAlternativeOffsets()) {
                std::cout << "Ambiguous Node: (alternative offsets: " << node.Data().m_offsets.size() << ")" << std::endl;
                std::cout << node.ToString() << std::endl;
                std::cout << node.Data().m_offset.ToString() << std::endl;
                for (auto& offset : node.Data().m_offsets) {
                    std::cout << offset.ToString2() << std::endl;
                }
            }
        }
    }

    void Classifier::InitializePatternProcessor(ds::PatternMap* pattern_db, bool sensitive_mode) {
        processor = std::make_unique<PatternHandler>(m_taxonomy, pattern_db, sensitive_mode);
    }

    ReadOffset &Classifier::Evaluate() {

//        SumUpNodes();
        return m_tree.Root().Data().m_offset;
    }

    Lineage *Classifier::GetLineages() {
        return m_lineages;
    }


    TID *Lineage::begin() {
        return m_lineage;
    }

    TID *Lineage::end() {
        return m_lineage + m_rank_id + 1;
    }

    TID *Lineage::rbegin() {
        return m_lineage + m_rank_id;
    }

    TID *Lineage::rend() {
        return m_lineage - 1;
    }

    void Lineage::SetRankId(int rank_id) {
        m_rank_id = rank_id;
    }

    void Lineage::Set(int rank_id, TID rank_taxonomic_id) {
        //TODO
        assert(rank_id < 10);
        if (rank_id < 0) {
            return;
        }
        m_lineage[rank_id] = rank_taxonomic_id;
    }

    std::string Lineage::ToString() const {
        std::string result;
        result += "(" + std::to_string(m_rank_id) + ") ";
        for (int i = 0; i < LINEAGE_SIZE-1; i++) {
            result += std::to_string(m_lineage[i]);
            result += "\t";
        }
        result += std::to_string(m_lineage[LINEAGE_SIZE-1]);
        result += " taxid: " + std::to_string(m_lineage[m_rank_id]);

        return result;
    }

    bool Lineage::NotInit() {
        std::cout << ToString() << std::endl;
        for (int i = 0; i < LINEAGE_SIZE; i++) {
            if (m_lineage[i] != -1) return false;
        }
        return true;
    }

    std::pair<size_t, TID> Lineage::GetRankIdAndTaxId() const {
        return { m_rank_id, m_lineage[m_rank_id] };
    }

    auto Lineage::GetRankId() const {
        return m_rank_id;
    }

    std::pair<size_t, TID> Lineage::GetLCARankIdAndTaxId(const Lineage& other) const {
        int i = std::min(GetRankId(), other.GetRankId());
        for (; i >= 0; i--) {
            if (m_lineage[i] == other.m_lineage[i]) {
                if (m_lineage[i] == -1) continue;
                else {
                    return { i, m_lineage[i] };
                }
            }
        }
//        auto taxonomic_id = i == -1 || m_lineage[i] == -1 ? ROOT_NODE_ID : m_lineage[i];
        return { -1, ROOT_NODE_ID };
    }

    TID Lineage::GetTaxonomicId() const {
        return m_lineage[m_rank_id];
    }

    bool Lineage::IsRoot() const {
        return m_rank_id == -1;
    }

    const TID Lineage::GetNextSetTaxonomicId(size_t rank_id) const {
        while(rank_id < LINEAGE_SIZE && m_lineage[rank_id] == -1) rank_id++;
        return(m_lineage[rank_id]);
    }

    const TID Lineage::GetPrevSetTaxonomicId(size_t rank_id) const {
        rank_id -= rank_id != 0;
        while(rank_id > 0 && m_lineage[rank_id] == -1) rank_id--;
        if (m_lineage[rank_id] == -1) {
            return ROOT_NODE_ID;
        }
        return m_lineage[rank_id];
    }

    bool Lineage::IsIdAncestor(int taxid) {
        for (int i = m_rank_id; i >= 0; i--) {
            if (m_lineage[i] == taxid) return true;
        }
        return taxid == 1;
    }


    //########################################################################
    // Classifier 2
    //########################################################################

    Classifier2::Classifier2(Taxonomy &taxonomy) :
            m_taxonomy(taxonomy),
            m_tree(m_taxonomy.MaxLeafId()),
            m_max_leaf_id(m_taxonomy.MaxLeafId()) {
        LoadLineage();
    }

    void Classifier2::LoadLineage() {
        m_lineages_size = m_taxonomy.MaxTaxid() + 1;
        m_lineages = new Lineage[m_lineages_size];
//        std::cout << "LoadLineage with size:  " << m_lineages_size << std::endl;
//        std::cout << "Taxonomysize: " << m_taxonomy.map.size() << std::endl;
//        std::cout << "MaxLeafId:    " << m_taxonomy.MaxLeafId() << std::endl;
//        std::cout << "MaxTaxid:     " << m_taxonomy.MaxTaxid() << std::endl;
//        std::cout << "RootNode: " << m_taxonomy.Get(1).ToString() << std::endl;
        for (const auto& [id, node] : m_taxonomy.map) {
            auto cid = node.id;
            assert(id < m_lineages_size);

            m_lineages[id].SetRankId(Taxonomy::Rank2Id(node.rank));
            if (Taxonomy::Rank2Id(node.rank) == -1 && !m_taxonomy.IsRoot(node.id)) {
                std::cout << node.ToString() << std::endl;
                std::cout << "id:          " << node.id << std::endl;
                std::cout << "rank:        " << node.rank << std::endl;
                std::cout << "parent_id:   " << node.parent_id << std::endl;
                std::cout << "external_id: " << node.external_id << std::endl;
                std::cout << "level:       " << node.level << std::endl;
                std::cout << "Check taxonomy file. 'no rank' is reserved for  root node" << std::endl;
                exit(2);
            }
            while (!m_taxonomy.IsRoot(cid)) {
                if (cid == -1) std::cout << "HERE I AM CID" << std::endl;
                const auto& n = m_taxonomy.Get(cid);
                assert(id < m_lineages_size);
                m_lineages[id].Set(Taxonomy::Rank2Id(n.rank), n.id);
                cid = n.parent_id;
            }
        }
    }


    LineageTIDMap Classifier2::GetLineages() {
        return m_lineages;
    }


    void Classifier2::AddMiss(bool reverse) {
        m_processor[reverse]->AddMiss();
    }

    bool Classifier2::IsLeaf(TID taxonomic_id, GID gene_id, GPOS gene_position) {
        return taxonomic_id > m_tree.RootId() && taxonomic_id <= m_max_leaf_id &&
               gene_id > 0 && gene_position > 0;
    }

    void Classifier2::AddHit(TID taxonomic_id, GID gene_id, GPOS gene_position, GPOS read_position, bool second) {
        assert(taxonomic_id >= 1 && taxonomic_id < m_lineages_size);
        if (IsLeaf(taxonomic_id, gene_id, gene_position)) {
            AddLeaf(taxonomic_id, gene_id, gene_position, read_position, second);
            m_processor[second]->AddLeaf(taxonomic_id, gene_id, gene_position);

        } else {
            AddNode(taxonomic_id, 0, second);
//            AddNode(taxonomic_id, gene_id, second);
            m_processor[second]->AddNode(taxonomic_id);
        }
    }

    void Classifier2::AddLeaf(TID taxonomic_id, GID gene_id, GPOS gene_position, GPOS read_position, bool second) {
        auto node_id = NodeId::ToNodeId(taxonomic_id, gene_id);
        auto iterator = m_lookup_idx_map.find(node_id);

        size_t index;
        if (iterator == m_lookup_idx_map.end()) {
            index = m_lookup_idx_map.size();
            m_lookup_idx_map.insert( { node_id, index } );
            m_lookup_results.emplace_back(
                    // Add new LookupResult to vector.
                    LookupResult(m_lineages[taxonomic_id],
                                 taxonomic_id,
                                 gene_id) );

            assert(m_lookup_results.at(index).GetHit(second).GetMaxHits() == 0);
        } else {
            index = iterator->second;
        }

        auto& hit = m_lookup_results.at(index).GetHit(second);

        if (!hit.IsInit()) {
            hit.Set(gene_position, read_position);
        } else {
            hit.IncrementIfMatch(gene_position, read_position);
        }



//        if (iterator != m_lookup_idx_map.end()) {
//            // Node already exists
//            auto& node = m_lookup_results.at(iterator->second);
//
//            node.GetHit(second).IncrementIfMatch(gene_position, read_position);
//        } else {
//            m_lookup_idx_map.insert( { node_id, m_lookup_idx_map.size() } );
//            m_lookup_results.emplace_back(
//                    // Add new LookupResult to vector.
//                    LookupResult(m_lineages[taxonomic_id],
//                                 taxonomic_id,
//                                 gene_id,
//                                 gene_position,
//                                 read_position,
//                                 second) );
//        }
    }



    void Classifier2::AddNode(TID taxonomic_id, GID gene_id, bool second) {
        auto node_id = NodeId::ToNodeId(taxonomic_id);
        auto iterator = m_lookup_idx_map.find(node_id);


        if (iterator != m_lookup_idx_map.end()) {
            // Node already exists
            auto& node = m_lookup_results.at(iterator->second);
            node.GetHit(second).Increment();
        } else {
            auto index = m_lookup_idx_map.size();
            m_lookup_idx_map.insert( { node_id, index } );
            m_lookup_results.emplace_back(
                    // Add new LookupResult to vector.
                    LookupResult(m_lineages[taxonomic_id],
                                 taxonomic_id,
                                 gene_id) );
            m_lookup_results.at(index).GetHit(second).Increment();
        }
    }

    void Classifier2::ProcessWorker() {
        static const auto comparator = [](const LookupResult& a, const LookupResult& b) -> bool {
            for (uint32_t i = 0; i < Lineage::LINEAGE_SIZE; i++) {
                //If values are same continue
                if (a.GetLineage().GetRankId() > b.GetLineage().GetRankId()) return true;

                if (a.GetLineage().m_lineage[i] != b.GetLineage().m_lineage[i]) {
                    return a.GetLineage().m_lineage[i] > b.GetLineage().m_lineage[i];
                }
            }
            return a.GeneId() > b.GeneId();
        };

        static const auto comparator2 = [](const LookupResult& a, const LookupResult& b) -> bool {
            //std::cout << "Memcmp result: " << std::flush;
            static constexpr auto lin_bytes = Lineage::LINEAGE_SIZE * sizeof(TID);
//            if (a.GetLineage().GetRankId() != b.GetLineage().GetRankId()) return a.GetLineage().GetRankId() < b.GetLineage().GetRankId();

            auto cmp_lin = std::memcmp(a.GetLineage().m_lineage, b.GetLineage().m_lineage, lin_bytes);
            if (cmp_lin != 0) return cmp_lin > 0;

            return a.GeneId() > b.GeneId();
        };

        std::sort(m_lookup_results.begin(), m_lookup_results.end(), comparator2);

//        std::cout << "Sorted lookups" << std::endl;
//        for (auto& lr : GetLookupResults()) {
//            std::cout << lr.ToString() << std::endl;
//        }

    }

    void Classifier2::Clear() {
        m_lookup_idx_map.clear();
        m_lookup_results.clear();
        m_tree.Clear();
        m_processor[0]->Reset();
        m_processor[1]->Reset();
    }

    std::vector<LookupResult>& Classifier2::GetLookupResults() {
        return m_lookup_results;
    }

    void Classifier2::BuildTree() {
//        std::cout << std::string(79, '#') << std::endl;
        Lineage working_lineage;

        // Early return if no hits
        if (m_lookup_results.empty())
            return;

        // Insert root
        if (m_lookup_results.at(0).GetTaxonomicId() == m_tree.RootId()) {
            auto node_id = m_tree.RootNodeId();
            m_tree.Insert(node_id, CNode(node_id, node_id, m_lookup_results.at(0).GetHitPE()));
        } else {
            m_tree.InsertRoot(HitPE(1));
        }

        // p* parent, l* lca, t* current taxon
        // *r rank, *t taxid
        // p_ previous
        // e.g. pr = parent rank
        int32_t current_lca_rank_id = 0, current_rank_id = 0;
        int32_t parent_rank_id = 1, previous_lca_rank_id = 0, p_tr = 0;
        TID current_lca_taxon_id = 1, current_taxonomic_id = 1;
        TID parent_taxon_id = 1, previous_lca_taxon_id = 1, p_tt = 1;

        UNID previous_tax_node_id = 0;


//        m_tree.PrettyPrint(PrintCNode);
//        std::cout << "Enter loop. " << std::endl;

        bool is_root = m_lookup_results.at(0).GetLineage().IsRoot();
        auto* previous_lineage = is_root ?
                &m_lookup_results.at(0).GetLineage() :
                &working_lineage;


        // Iterate all lookup results

        for (auto i = is_root ? 1 : 0; i < m_lookup_results.size(); i++) {
//            std::cout << std::string(30, '^') << std::endl;
            auto& result = m_lookup_results.at(i);                      // result
            auto& lineage = result.GetLineage();                        // result->lineage

            std::tie(current_rank_id, current_taxonomic_id) = lineage.GetRankIdAndTaxId();
            std::tie(current_lca_rank_id, current_lca_taxon_id) = lineage.GetLCARankIdAndTaxId(*previous_lineage);


            auto previous_lca_node_id = NodeId::ToNodeId(previous_lca_taxon_id);
            auto current_lca_node_id = NodeId::ToNodeId(current_lca_taxon_id);

            auto current_tax_node_id = IsLeaf(current_taxonomic_id, result.GeneId(), 1) ?
                    NodeId::ToNodeId(current_taxonomic_id, result.GeneId()) :
                    NodeId::ToNodeId(current_taxonomic_id);

//            std::cout << "___________________________" << std::endl;
//            std::cout << "Working lineage:  " << working_lineage.ToString() << std::endl;
//            std::cout << "Previous lineage: " << previous_lineage->ToString() << std::endl;
//            std::cout << "Current lineage:  " << lineage.ToString() << std::endl;
//            std::cout << "LCA, node: " << current_lca_taxon_id << "(" << current_lca_node_id << ")\t\tTaxonomic Id: " << current_taxonomic_id << "(" << current_rank_id << ")" << std::endl;

//            std::cout << "1:" << working_lineage.ToString() << std::endl;
            // LCA is root.
            if (current_lca_taxon_id == 1) {
                parent_rank_id = -1;
            } else {
                if (working_lineage.m_lineage[current_lca_rank_id] == -1 || working_lineage.m_lineage[current_lca_rank_id] != current_taxonomic_id) {
//                    std::cout << "_________________________" << std::endl;
                    // Create LCA if not existent.
//                    std::cout << "Previous LCA: " << previous_lca_taxon_id << "(" << previous_lca_node_id << ")" << std::endl;
                    m_tree.Insert(current_lca_node_id, CNode(current_lca_node_id, previous_lca_node_id, HitPE(current_lca_taxon_id)));

                    //                    std::cout << "Insert: " << current_lca_taxon_id << "(" << current_lca_node_id << ")" << std::endl;
                    working_lineage.m_lineage[current_lca_rank_id] = current_lca_taxon_id;
                    // Reset any working_lineage above current lca rank id (current_lca_rank_id)
                    for (int r = current_lca_rank_id + 1; r < Lineage::LINEAGE_SIZE; r++) working_lineage.m_lineage[r] = -1;
//                    std::cout << "_________________________" << std::endl;
                }

//                std::cout << "working_lineage.m_lineage[current_lca_rank_id] == -1 " << (working_lineage.m_lineage[current_lca_rank_id] == -1) << std::endl;
//                std::cout << "working_lineage.m_lineage[current_lca_rank_id] != current_taxonomic_id " << (working_lineage.m_lineage[current_lca_rank_id] != current_taxonomic_id) << std::endl;
//                std::cout << "current_taxonomic_id: " <<  current_taxonomic_id << "  working_lineage.m_lineage[current_lca_rank_id]: "<< working_lineage.m_lineage[current_lca_rank_id] << std::endl;
//                std::cout << working_lineage.ToString() << std::endl;
                auto& current_lca_node = m_tree.GetNode(current_lca_node_id);
//                m_tree.PrettyPrint(PrintCNode);

//                std::cout << (current_lca_rank_id > previous_lca_rank_id ? "LCA moved up" : "LCA moved down") << std::endl;
//                std::cout << "Step in? " << (current_lca_node_id != previous_tax_node_id) << std::endl;
                if (current_lca_node_id != previous_tax_node_id && current_lca_rank_id > previous_lca_rank_id) { // current lca rank is higher than previous lca rank
                    // unlink previous lca and previous taxon _________________
                    auto& previous_lca_node = m_tree.GetNode(NodeId::ToNodeId(previous_lca_taxon_id));

                    auto __assert_size = previous_lca_node.Children().size();

//                    std::cout << "(" << previous_lca_taxon_id << ").RemoveChild(" << NodeId::ToString(previous_tax_node_id) << ")" << std::endl;
//                    std::cout << "Children previous_lca_taxon_id: " << NodeId::ToString(previous_lca_node.GetNodeId()) << "  Children(";
//                    for (auto& child : previous_lca_node.Children()) { std::cout << NodeId::ToString(child) << "|"; };
//                    std::cout << ")" << std::endl;


                    previous_lca_node.RemoveChild(previous_tax_node_id);

                    // after child removal Children list must be smaller by one
                    assert(__assert_size - 1 == previous_lca_node.Children().size());

                    // Fit previous LCA node is before current lca node
                    previous_lca_node.AddChild(current_lca_node_id);
                    // Relink old taxid with new intermediate
                    current_lca_node.AddChild(previous_tax_node_id);
//                    std::cout << "(" << previous_lca_taxon_id << ").AddChild(" << NodeId::ToString(current_lca_node_id) << "  Children(" << std::endl;
//                    std::cout << "Children previous_lca_taxon_id: " << NodeId::ToString(previous_lca_node.GetNodeId()) << "  Children(";
//                    for (auto& child : previous_lca_node.Children()) { std::cout << NodeId::ToString(child) << "|"; };
//                    std::cout << ")" << std::endl;
//                    std::cout << "Children current_lca_node: ";
//                    for (auto& child : current_lca_node.Children()) { std::cout << NodeId::ToString(child) << "|"; };
//                    std::cout << ")" << std::endl;
//
//                    std::cout << ">>>>>>>> ListNodes" << std::endl;
//                    m_tree.ListNodes();
//                    std::cout << "<<<<<<<< ListNodes" << std::endl;
//
//                    std::string stop;
//                    std::cin >> stop;
//                    m_tree.PrettyPrint(PrintCNode);

                    current_lca_node.SetParent(previous_lca_node_id);
//                    std::cout << "Add child and set parent" << std::endl;
                    assert(__assert_size == previous_lca_node.Children().size());
                    // __________________unlink previous lca and previous taxon
                } else {
                    parent_rank_id = current_lca_rank_id;
                    for (; parent_rank_id >= 0; parent_rank_id--) {
                        if (working_lineage.m_lineage[parent_rank_id] != -1) break;
                    }
                    // previous_lca_taxon_id last created node left to
                    parent_taxon_id = parent_rank_id == -1 ? 1 : working_lineage.m_lineage[parent_rank_id];

                    auto parent_node_id = NodeId::ToNodeId(parent_taxon_id);
                    auto parent_node = m_tree.GetNode(parent_node_id);
                    parent_node.AddChild(current_lca_node_id);
                    current_lca_node.SetParent(parent_node_id);
                }
            }

//            std::cout << std::string(20, '-') << std::endl;
//            std::cout << "Current LCA: " << NodeId::ToString(current_lca_node_id) << std::endl;
//
//            std::cout << ">>>>>>>> ListNodes" << std::endl;
//            m_tree.ListNodes();
//            std::cout << "<<<<<<<< ListNodes" << std::endl;

            // Insert current taxon finally
//            auto is_leaf = IsLeaf(current_taxonomic_id, result.GeneId(), 1);
            auto is_leaf = result.GetHitPE().GetHit(0).IsLeaf() || result.GetHitPE().GetHit(1).IsLeaf();
            auto node_id = is_leaf ?
                    NodeId::ToNodeId(current_taxonomic_id, result.GeneId()) :
                    NodeId::ToNodeId(current_taxonomic_id);


            if (node_id != current_lca_node_id) {
                // current node has not been inserted yet
                m_tree.Insert(node_id, CNode(node_id, current_lca_node_id, result.GetHitPE()));

                if (!is_leaf) working_lineage.m_lineage[current_rank_id] = current_taxonomic_id;
                m_tree.GetNode(current_lca_node_id).AddChild(node_id);
            } else {
                // if current node id == lca node id that means we have already inserted that node
//                std::cout << "(" << NodeId::ToString(current_lca_node_id) << ").AddChild(" << NodeId::ToString(node_id) << "  Children(" << std::endl;
//                m_tree.PrettyPrint(PrintCNode);
            }

            assert(m_tree.HasNode(node_id));

            previous_lca_rank_id = current_lca_rank_id;
            previous_lca_taxon_id = current_lca_taxon_id;
            previous_tax_node_id = node_id;
            previous_lineage = &result.GetLineage();

//            std::string stop;
//            std::cin >> stop;
        }
    }

    CTree &Classifier2::GetTree() {
        return m_tree;
    }

    bool Classifier2::HasNode(TID taxonomic_id) {
        return m_lookup_idx_map.contains(NodeId::ToNodeId(taxonomic_id));
    }

    void Classifier2::InitializePatternProcessor(ds::PatternMap *pattern_db, bool sensitive_mode) {
        if (!pattern_db) {
            std::cerr << "Pattern DB not initialized" << std::endl;
            exit(9);
        }
        m_processor[0] = std::make_unique<PatternHandler>(m_taxonomy, pattern_db, sensitive_mode);
        m_processor[1] = std::make_unique<PatternHandler>(m_taxonomy, pattern_db, sensitive_mode);
    }

    void Classifier2::AddUpTree(ClassificationResult &result) {
        auto& root = m_tree.Root();
        root.Data().GetHit(0).AddToTotal();
        root.Data().GetHit(1).AddToTotal();
        root.Data().m_total_count = root.Data().GetHit(0).Total() + root.Data().GetHit(1).Total();
        if (root.Children().empty()) result.m_hits.emplace_back(&root);

        for (auto& child : root.Children()) {
            auto& child_node = GetTree().GetNode(child);
            AddUpTree(result, root, child_node);
        }
//        AddUpTree(result, root, root);
    }

    void Classifier2::AddUpTree(ClassificationResult &result, CNode& parent, CNode& target) {
        target.Data().GetHit(0).AddToTotalFromOther(parent.Data().GetHit(0));
        target.Data().GetHit(1).AddToTotalFromOther(parent.Data().GetHit(1));
        target.Data().m_total_count = target.Data().GetHit(0).Total() + target.Data().GetHit(1).Total();

        target.Data().m_rank_confidence = (double) target.Data().NodeHitsOnly() / target.Data().m_total_count;

        if (target.Children().empty()) {
            result.m_hits.emplace_back(&target);
        }

        for (auto& child : target.Children()) {
            auto& child_node = GetTree().GetNode(child);
            AddUpTree(result, target, child_node);
        }
    }

    std::pair<size_t,size_t> Classifier2::GetSumOfHitsForTaxon(TID taxid, GID gene_id) {
        std::cout << "\n______GetSumOfHitsForTaxon(" << taxid << ", " << gene_id << ")" << std::endl;
        size_t sum_a = 0;
        size_t sum_b = 0;
        for (auto& lr : GetLookupResults()) {
            std::cout << lr.ToString();
//            auto gene_id = lr.GetHitPE().m_gene_id;
//            std::cout << "m_taxonomy.IsNodeAncestor(lr.GetTaxonomicId(), taxid): " << m_taxonomy.IsNodeAncestor(lr.GetTaxonomicId(), taxid) << std::endl;
//            std::cout << "(!m_taxonomy.IsLeaf(lr.GetTaxonomicId()) || gene_id < 1 || gene_id == lr.GeneId()): " << (!m_taxonomy.IsLeaf(lr.GetTaxonomicId()) || gene_id < 1 || gene_id == lr.GeneId()) << std::endl;
//            std::cout << gene_id << " == " << lr.GeneId() << std::endl;
            if (m_taxonomy.IsNodeAncestor(lr.GetTaxonomicId(), taxid) && (!m_taxonomy.IsLeaf(lr.GetTaxonomicId()) || gene_id < 1 || gene_id == lr.GeneId())) {
                sum_a += lr.GetHit(0).GetMaxHits();
                sum_b += lr.GetHit(1).GetMaxHits();
                std::cout << "            sumup: " << sum_a << ", " << sum_b;
            }
            std::cout << std::endl;
        }
        return { sum_a, sum_b };
    }

    bool Classifier2::ProcessHits3(ClassificationResult& result, ClassificationResult& result2, double min_confidence, size_t min_hits) {
        constexpr bool debug = false;

        // Sort lookup results (this is mandatory)
        ProcessWorker();

        if constexpr (debug) {
            std::cout << "\n______Process Hits:" << std::endl;
            PrintLookupResults();
            std::cout << "----------" << std::endl;
        }

        result.Clear();
        result2.Clear();

        size_t sum_hits = 0;
        size_t sum_hits_a = 0;
        size_t sum_hits_b = 0;

        int range_start = 0;
        int range_end = GetLookupResults().size();
        int rank_level = 0;

        struct Range {
            uint16_t start, end;
            uint32_t count;
            std::string ToString() const {
                return(std::to_string(start) + ',' + std::to_string(end) + ": " + std::to_string(count));
            }
            Range(uint16_t start, uint16_t end, uint32_t count) : start(start), end(end), count(count){};
        };
        std::vector<Range> sums;

        TID lca_id = 1;
        int lca_index = -1;

        auto& lrs = GetLookupResults();

        // If first lookup result is root change deal with it now and start processing at index 1
        if (lrs[0].GetLineage().IsRoot()) {
            sum_hits = GetLookupResults()[0].GetHitPE().NodeHitsOnly();
            sum_hits_a = GetLookupResults()[0].GetHitPE().GetHit(0).GetMaxHits();
            sum_hits_b = GetLookupResults()[0].GetHitPE().GetHit(1).GetMaxHits();
            range_start = 1;
            lca_index = 0;
        }

        int loop_num = 0;
        while (range_start < range_end - 1) { // Criterion for keeping on looking
            loop_num++;
            // GO UP RANK LEVELS OR INCREASE RANGE START UNTIL TAXID OF START AND END OF RANGE
            // AT RANK LEVEL ARE NOT THE SAME ANYMORE
            while (rank_level < Lineage::LINEAGE_SIZE &&
                    lrs[range_start].GetLineage().m_lineage[rank_level] == lrs[range_end-1].GetLineage().m_lineage[rank_level] &&
                    (range_end-range_start > 1)) {

                assert(range_end > range_start + 1);

                // If current the current range start's RankId matches the current rank level
                // That means it is an inner node of the next nodes.
                // In a leaf hit case, we might have two subsequent hits on the same taxonomic id but on different genes
                // E.g. Read1 maps to gene 1, Read2 maps to gene 2

                // if (CURRENT NODE IS NOT START OF A RANGE BUT ITS LCA) -> move range start by one
                // else (CURRENT NODE IS START OF A RANGE) -> move rank up by one
                if (lrs[range_start].GetLineage().GetRankId() == rank_level) {
                    if  (lrs[range_start].GetTaxonomicId() == lrs[range_start+1].GetTaxonomicId()) break;

                    size_t sum = 0;
                    bool pass = false;
                    for (auto i = range_start+1; i < range_end && sum < min_hits; i++) {
                        sum += lrs[i].GetHitPE().NodeHitsOnly();
                        assert(lrs[range_start].GetLineage().GetTaxonomicId() == lrs[range_start].GetLineage().m_lineage[rank_level]);
                    }

                    // Sum beneath node don't go down further (not enough evidence)
                    if (sum >= min_hits) {
                        if constexpr (debug) {
                            sum_hits += lrs[range_start].GetHitPE().NodeHitsOnly();
                            sum_hits_a += lrs[range_start].GetHitPE().GetHit(0).GetMaxHits();
                            sum_hits_b += lrs[range_start].GetHitPE().GetHit(1).GetMaxHits();
                            std::cout << "1 sum_hits_a: " << sum_hits_a << std::endl;
                            std::cout << "1 sum_hits_b: " << sum_hits_b << std::endl;
                        }
                        range_start++;
                    } else {
                        // sum_hits is above node (have there been enough ancestry hits supporting this node)
                        if (sum_hits >= min_hits) {
                            lca_id = lrs[range_start].GetTaxonomicId();
                            lca_index = range_start;
                            range_end = range_start + 1;

                            if constexpr (debug) {
                                std::cout
                                        << "Goto end____________________________________________________________________________"
                                        << std::endl;
                            }
                            goto end_loop;
                        } else {
                            return false;
                        }
                    }
                } else {
                    rank_level++;
                }
            }

            if (rank_level == Lineage::LINEAGE_SIZE) {
                std::cout << "(rank_level == Lineage::LINEAGE_SIZE)" << std::endl;
                std::cout << range_start << " " << range_end << std::endl;
                for (auto i = range_start; i < range_end; i++) {
                    std::cout << lrs[i].GetHitPE().ToVerboseString() << std::endl;
                }
                std::cout << "_______________________________" << std::endl;
                PrintLookupResults();

                std::cout << "## DEBUGDEBUGDEBUG" << std::endl;
                exit(12);
            }


            // LCA id is based off the last element in the range at the current rank level
            // rank level is lifted such that
            lca_id = rank_level > 0 ? lrs[range_end-1].GetLineage().GetPrevSetTaxonomicId(rank_level) : 1;
            if (lca_id == -1) {
                std::cerr << "rank level: " << rank_level << std::endl;
                std::cerr << "lrs[range_end-1].GetLineage().GetPrevSetTaxonomicId(rank_level):  " << lrs[range_end-1].GetLineage().GetPrevSetTaxonomicId(rank_level) << std::endl;
                std::cerr << "range_end:  " << range_end << std::endl;
            }

            lca_index = range_start;
            // If ID == LCA_ID that means the current node is the LCA node.
            // However, if there are two leafs with the same taxid, then

            if constexpr (debug) {
                std::cout << "\n____Check if node above range is LCA: " << std::endl;
                if (range_start > 0) {
                    std::cout << "ID above: " << lrs[range_start - 1].GetLineage().GetTaxonomicId() << std::endl;
                }
                std::cout << "Current lca: " << lca_index << std::endl;
            }
            // Set lca index to previous item, if previous item is the lca current range.
            lca_index -= range_start > 0 && lrs[range_start-1].GetLineage().GetTaxonomicId() == lca_id && (lrs[range_start].GetTaxonomicId() != lrs[range_start+1].GetTaxonomicId() || range_end-range_start == 1);


            TID previous_id = lrs[range_start].GetLineage().GetNextSetTaxonomicId(rank_level);

            uint16_t current_start = range_start;
            uint32_t count = lrs[range_start].GetHitPE().NodeHitsOnly();
            uint32_t total_count = count;

            sums.clear();

            // Get the count for different branching options to assess confidence
            for (uint16_t i = range_start+1; i < range_end; i++) {
                if (previous_id != lrs[i].GetLineage().GetNextSetTaxonomicId(rank_level) || rank_level == lrs[i].GetLineage().GetRankId()) {

                    previous_id = lrs[i].GetLineage().GetNextSetTaxonomicId(rank_level);
                    sums.emplace_back(Range{current_start, i, count});
                    count = 0;
                    current_start = i;
                }
                count += lrs[i].GetHitPE().NodeHitsOnly();
                total_count += lrs[i].GetHitPE().NodeHitsOnly();
            }
            sums.emplace_back(Range{current_start, (uint16_t) range_end, count});



            // Decide if tere is enough confidence for one of the branches to continue
            bool exit = true;
            if constexpr (debug) {
                std::cout << "\n______Check Range: " << std::endl;
                std::cout << "Fallback to : " << lca_id << " " << lca_index << std::endl;
                std::cout << "-----" << std::endl;
            }
            for (auto& range : sums) {
                double ratio = ((double) range.count)/total_count;

                if constexpr (debug) {
                    std::cout << "range: " << range.ToString() << " " << ratio << " "
                              << (ratio > min_confidence && range.count >= min_hits) << std::endl;
                }
                if (ratio > min_confidence && range.count >= min_hits) {
                    range_start = range.start;
                    range_end = range.end;
                    exit = false;
                    break;
                }
            }
            if constexpr (debug) {
                std::cout << "-----" << std::endl;
                std::cout << "exit: " << exit << " current lca = " << lca_id << " (" << (lca_index >= 0 ? lrs[lca_index].ToString() : std::to_string(lca_index)) << ")" << std::endl;
            }
            if (exit) break;
        }
        end_loop:

        if constexpr (debug) {
            std::cout << "\n_______End loop sum up hit counts next" << std::endl;
            std::cout << "current: " <<  sum_hits_a << ", " << sum_hits_b << std::endl;
            std::cout << "Range: " << range_start << " -> " << range_end << std::endl;
        }

        // Add up final counts for range if range_size = 1
        if (range_end - range_start == 1) {
            lca_index = range_start;
            lca_id = lrs[lca_index].GetLineage().GetTaxonomicId();

            sum_hits += lrs[range_start].GetHitPE().NodeHitsOnly();
            sum_hits_a += lrs[range_start].GetHitPE().GetHit(0).GetMaxHits();
            sum_hits_b += lrs[range_start].GetHitPE().GetHit(1).GetMaxHits();

            if constexpr (debug) {
                std::cout << "\n_______Range is one" << std::endl;
                std::cout << "taxid:  " << lca_id << std::endl;
                std::cout << "sum_hits_a: " << sum_hits_a;
                std::cout << ", sum_hits_b: " << sum_hits_b << std::endl;
            }
        }
        assert(sum_hits == sum_hits_a + sum_hits_b);

        result.m_taxonomic_id = lca_id;
        result.m_gene_id = 0;

        if (lca_id == -1) {
            std::cout << "LOCATION HERE" << std::endl; // REMOVE
            PrintLookupResults();
            std::cout << "LOCATION HERE" << std::endl; // REMOVE
            lca_id = 1;
        }
        bool leaf = m_taxonomy.IsLeaf(lca_id);
        bool classified = false;


        if constexpr(debug) {
            std::cout << "\t1. Read1 pass? " << sum_hits_a << " > " << sum_hits << " = " <<  (sum_hits_a > min_hits) << std::endl;
            std::cout << "\t1. Read2 pass? " << sum_hits_b << " > " << sum_hits << " = " <<  (sum_hits_b > min_hits) << std::endl;
        }

        // First read has passed the threshold, set stats
        if (sum_hits_a > min_hits) {
            //TODO computation of sum_hits_a seems to be flawed
            result.SetNodeHit(lca_id, sum_hits_a);

            auto& hit = m_lookup_results[lca_index].GetHit(false);

            if (leaf && hit.IsLeaf()) {
                if constexpr(debug) {
                    std::cout << "\t\t\tSet hit for read1:   " << hit.ToString() << std::endl;
                }
                result.SetLeafHit(&hit, m_lookup_results[lca_index].GeneId());
            }
            classified = true;
        }

        // Second read has passed the threshold, set stats
        if (sum_hits_b > min_hits) {
            result2.SetNodeHit(lca_id, sum_hits_b);

            auto& hit = m_lookup_results[lca_index].GetHit(true);

            if (leaf && hit.IsLeaf()) {
                if constexpr(debug) {
                    std::cout << "\t\t\tSet hit for read2:   " << hit.ToString() << std::endl;
                }
                result2.SetLeafHit(&hit, m_lookup_results[lca_index].GeneId());
            }
            classified = true;
        }

//        auto [sum_a, sum_b] = GetSumOfHitsForTaxon(result.m_taxonomic_id, result.m_gene_id);

//        if (!(sum_a == sum_hits_a && sum_b == sum_hits_b)) {
//            std::cout << "\nError in ProcessHits3" << std::endl;
//            PrintLookupResults();
//            std::cout << "\n" << result.m_taxonomic_id << std::endl;
//            std::cout << "Geneid:       " << result.m_gene_id << std::endl;
//            std::cout << "sum_hits_a:   " << sum_hits_a << " != " << sum_a << std::endl;
//            std::cout << "sum_hits_b:   " << sum_hits_b << " != " << sum_b << std::endl;
//
//            std::cout << std::endl;
//        }

        if constexpr(debug) {
            std::cout << "\t2. Read1: " << result.Success() << " and Leaf? " << result.IsLeafHit() << "    "
                      << (result.m_hit ? result.m_hit->ToString() : "") << std::endl;
            std::cout << "\t2. Read2: " << result2.Success() << " and Leaf? " << result2.IsLeafHit() << "    "
                      << (result2.m_hit ? result2.m_hit->ToString() : "") << std::endl;
        }
        if (!classified) return false;
        return true;

        if constexpr (debug) {
            auto [sum_a_tmp, sum_b_tmp] = GetSumOfHitsForTaxon(result.m_taxonomic_id, result.m_gene_id);

//        if (true) {
            if (sum_a_tmp != sum_hits_a || sum_b_tmp != sum_hits_b) {
                std::cout << result.m_taxonomic_id << " gene: " << result.m_gene_id << std::endl;
                std::cout << "taxid: " << lca_id << m_lookup_results[lca_index].ToString() << std::endl;
                std::cout << "sum_hits_a: " << sum_hits_a << ", " << sum_a_tmp << " : sum_a_tmp" << std::endl;
                std::cout << "sum_hits_b: " << sum_hits_b << ", " << sum_b_tmp << " : sum_b_tmp" << std::endl;
                GetTree().PrettyPrint(classification::PrintCNode);
                Utils::Input();
            }
        }
        assert(GetSumOfHitsForTaxon(result.m_taxonomic_id, result.m_gene_id).first == sum_hits_a &&
               GetSumOfHitsForTaxon(result.m_taxonomic_id, result.m_gene_id).second == sum_hits_b);



        return true;
    }

    bool Classifier2::ProcessHits2(ClassificationResult& result, double min_confidence, size_t min_hits) {
//        ProcessWorker();

//        PrintLookupResults();

        size_t sum_hits = 0;
        int range_start = 0;
        int range_end = GetLookupResults().size();
        int rank_level = 0;

        struct Range {
            uint16_t start, end;
            uint32_t count;
            std::string ToString() {
                return(std::to_string(start) + ',' + std::to_string(end) + ": " + std::to_string(count));
            }
            Range(uint16_t start, uint16_t end, uint32_t count) : start(start), end(end), count(count){};
        };
        std::vector<Range> sums;

        if (GetLookupResults()[0].GetLineage().IsRoot()) {
            sum_hits = GetLookupResults()[0].GetHitPE().NodeHitsOnly();
            range_start = 1;
        }

        TID lca_id = 1;
        int lca_index = -1;

        auto& lrs = GetLookupResults();

        while (range_start < range_end - 1) { // Criterion for keeping on looking
//            std::cout << std::string(79, '-') << " " << range_start << " - " << range_end << std::endl;
            while (rank_level < Lineage::LINEAGE_SIZE && lrs[range_start].GetLineage().m_lineage[rank_level] ==
                    lrs[range_end-1].GetLineage().m_lineage[rank_level] && (range_end-range_start > 1)) {
                // If current the current range start's RankId matches the current rank level
                // That means it is an inner node of the next nodes.
                if (lrs[range_start].GetLineage().GetRankId() == rank_level) {
                    size_t sum = 0;
                    bool pass = false;
                    for (auto i = range_start+1; i < range_end; i++) {
                        sum += lrs[i].GetHitPE().NodeHitsOnly();
                        if (lrs[range_start].GetLineage().GetTaxonomicId() != lrs[range_start].GetLineage().m_lineage[rank_level]) {
                            std::cout << range_start << " " << lrs[range_start].GetLineage().ToString() << " " << lrs[range_start].ToString() << std::endl;
                            std::cout << i << " " << lrs[i].GetLineage().ToString() << " " << lrs[i].ToString() << std::endl;
                            exit(9);
                        }
                    }
//                    std::cout << "X Sum: " << sum << std::endl;
                    if (sum >= min_hits) {
//                        std::cout << "sum >= min_hits  " << sum << " > " << min_hits << std::endl;
                        sum_hits += lrs[range_start].GetHitPE().NodeHitsOnly();
                        range_start++;
                    } else {
//                        std::cout << "STOP _______________________ " << std::endl;
//                        std::cout << range_start << " " << lrs[range_start].GetLineage().ToString() << " " << lrs[range_start].ToString() << std::endl;

                        if (sum_hits >= min_hits) {
                            result.m_taxonomic_id = lrs[range_start].GetTaxonomicId();
                            result.pattern_idx = range_start;
                            return true;
                        } else {
                            return false;
                        }
                    }
                } else {
                    rank_level++;
                }
            }
//            std::cout << "Final: " << range_start << " - " << range_end << std::endl;
            if (rank_level == Lineage::LINEAGE_SIZE) {
                std::cout << range_start << " " << range_end << std::endl;
                for (auto i = range_start; i < range_end; i++) {
                    std::cout << lrs[i].GetHitPE().ToString() << std::endl;
                }
                std::cout << "DEBUGDEBUGDEBUG" << std::endl;
                exit(9);
            }

//            std::cout << range_start << " " << lrs[range_start].GetLineage().m_lineage[rank_level] << " " << lrs[range_end -1].GetLineage().m_lineage[rank_level] << std::endl;

            lca_id = rank_level > 0 ? lrs[range_end-1].GetLineage().GetPrevSetTaxonomicId(rank_level) : 1;
            lca_index = range_start;
            lca_index -= range_start > 0 && lrs[range_start-1].GetLineage().GetTaxonomicId() == lca_id;
//            std::cout << "lca_id: "  << lca_id << " rank level: " << rank_level << std::endl;
//            std::cout << "lca_index: " << lca_index << std::endl;

//            // Is node lca
//            if (rank_level >= 0 && lrs[range_start].GetLineage().m_lineage[rank_level-1] == lrs[range_start].GetLineage().GetTaxonomicId()) {
//                std::cout << "Sum hits from " << sum_hits << " to ";
//                sum_hits += lrs[range_start].GetHitPE().NodeHitsOnly();
//                std::cout << sum_hits << " node: " << lrs[range_start].GetHitPE().ToVerboseString() << std::endl;
//                range_start++;
//                std::cout << "shift start to: " << range_start << std::endl;
//                continue;
//            }

            TID previous_id = lrs[range_start].GetLineage().GetNextSetTaxonomicId(rank_level);

            uint16_t current_start = range_start;
            uint32_t count = lrs[range_start].GetHitPE().NodeHitsOnly();
            uint32_t total_count = count;

            sums.clear();

//            std::cout << "Rank level:  " << rank_level << " -> id: " << lrs[range_start].GetLineage().m_lineage[rank_level] << std::endl;
//            std::cout << "Range start: " << range_start << " " << lrs[range_start].ToString() << std::endl;
//            std::cout << "Range end:   " << range_end << " " << lrs[range_end-1].ToString() << std::endl;

            // Get the count for different options
            for (uint16_t i = range_start+1; i < range_end; i++) {
                if (previous_id != lrs[i].GetLineage().GetNextSetTaxonomicId(rank_level)) {
                    previous_id = lrs[i].GetLineage().GetNextSetTaxonomicId(rank_level);
                    sums.emplace_back(Range{current_start, i, count});
                    count = 0;
                    current_start = i;
                }
                count += lrs[i].GetHitPE().NodeHitsOnly();
                total_count += lrs[i].GetHitPE().NodeHitsOnly();
            }
            sums.emplace_back(Range{current_start, (uint16_t) range_end, count});

//            std::cout << "total count : " << total_count << std::endl;
//            for (auto& range : sums) {
//                double ratio = ((double) range.count)/total_count;
//                std::cout << range.ToString() << " " << ratio << std::endl;
//            }
//            std::cout << "Decide?!" << std::endl;
            bool exit = true;
//            std::cout << "------------------------------------------" << std::endl;
            for (auto& range : sums) {
                double ratio = ((double) range.count)/total_count;
                if (ratio > min_confidence && range.count >= min_hits) {
//                    std::cout << "Yalla: " << range.count << " -> " << ratio << " -> " << range.ToString() << std::endl;
                    range_start = range.start;
                    range_end = range.end;
//                    std::cout << "range:start : " << range_start << std::endl;
                    exit = false;
                    break;
                }
            }
//            std::cout << "exit? " << exit << std::endl;
//            std::string stop;
//            std::cin >> stop;
            if (exit) break;
        }

        if (range_end - 1 == range_start) {
            lca_index = range_start;
            lca_id = lrs[lca_index].GetLineage().GetTaxonomicId();

            for (auto i = range_start; i < range_end; i++) {
                sum_hits += lrs[i].GetHitPE().NodeHitsOnly();
            }
        }


        result.m_taxonomic_id = lca_id;
        result.pattern_idx = lca_index;
//        std::cout << "End result: " << lca_index << "," << lca_id << std::endl;
//        if (lca_index > -1 && lca_index < GetLookupResults().size())
//            std::cout << lrs[lca_index].ToString() << std::endl;
//        std::cout << "Sum hits: " << sum_hits << std::endl;
//
//        std::string stop;
//        std::cin >> stop;

        return true;
    }

    void Classifier2::ProcessHits(ClassificationResult &result) {
//        std::cout << std::string(79, '-') << " sus double hit?" << std::endl;
//        for (auto& x : m_lookup_results) {
//            std::cout << x.ToString() << " " << std::endl;
//        }
//        for (auto& x : m_lookup_idx_map) {
//            std::cout << x.first << " " << x.second << " " << std::endl;
//        }
        ProcessWorker();
        BuildTree();
        result.Clear();
        AddUpTree(result);
        result.SortHits();

        assert(!result.m_hits.empty());

        // Stuff for candidate confidence ///////////
        int32_t min = INT32_MAX;
        auto max = 0;
        auto max_idx = 0;
        auto total = 0.0;

        for (auto i = 0; i < result.m_hits.size(); i++) {
            auto& node = result.m_hits.at(i);
            total += node->Data().TotalHits();
            if (node->Data().TotalHits() < min) min = node->Data().TotalHits();
        }
        min -= 1;
//        std::cout << "total: " << total << std::endl;
//        std::cout << "min: " << min << std::endl;
//        std::cout << "result.m_hits.size() * min: " << result.m_hits.size() * min << std::endl;
        total -= (double) result.m_hits.size() * min;
//        std::cout << total << std::endl;
//        std::cout << "List nodes.." << std::endl;

        for (auto& node : result.m_hits) {
            node->Data().m_candidate_confidence = (double) (node->Data().TotalHits() - min)/total;
//            std::cout << "+" << node->Data().TotalHits() << "+\t\t" << node->Data().m_taxonomic_id << " " << node->Data().m_gene_id << " " << node->Data().GetHit(0).ToString() <<  " " << node->Data().GetHit(1).ToString() << std::endl;
//            std::cout << node->Data().TotalHits() - min << " " << node->Data().m_candidate_confidence << std::endl;
        }

        // //////////////////////////////////////////


//        auto& best_node = result.m_hits.at(0);
        auto& best_node = result.GetBestNode();

        auto& best_hit_pe = best_node.Data();

//        m_tree.PrettyPrint(classification::PrintCNode);
//
//        std::cout << best_node->ToString() << std::endl;
//        std::cout << best_hit_pe.HasHit(0) << " " << best_hit_pe.HasHit(1) << std::endl;
//
//        std::string stop;
//        std::cin >> stop;

        bool has_hit_pe_first = best_hit_pe.HasHit(0);
        bool has_hit_pe_second = best_hit_pe.HasHit(1);

        if (has_hit_pe_first && has_hit_pe_second) {
            result.m_best = &best_hit_pe;
        } else {
            if (has_hit_pe_first) {
                result.m_best_first = &best_hit_pe;
            }
            if (has_hit_pe_second) {
                result.m_best_second = &best_hit_pe;
            }
        }
    }

    void Classifier2::SetReadLength(size_t l, bool second) {
        m_read_length = l;
        m_processor[second]->read_length_ = l;
//        m_processor[1]->read_length_ = l;
        m_processor[second]->cov_bucket_size = (l - m_processor[0]->snp_detector_.shape_size + 1) / 10;
//        m_processor[1]->cov_bucket_size = (l - m_processor[1]->snp_detector_.shape_size + 1) / 10;
    }

    const std::unique_ptr<PatternHandler> &Classifier2::GetProcessor(bool second) const {
        return m_processor[second];
    }

    void Classifier2::PrintLookupResults() const {
        for (auto& lookup_result : m_lookup_results) {
            std::cout << lookup_result.ToString() << " " << lookup_result.GetHitPE().m_taxonomic_id << ", " << lookup_result.GetHitPE().m_gene_id << std::endl;
        }
    }

    std::string Classifier2::Newick() {
        assert(m_tree.HasRoot());
        return Newick(m_tree.Root(), 0);
    }

    std::string Classifier2::Newick(CNode& node, unsigned int level) {
        std::string result = "";
        // Recursion
        if (!node.Children().empty()) {
            result += '(';
            result += Newick(m_tree.GetNode(node.Children()[0]), level+1);
            for (auto i = 1; i < node.Children().size(); i++) {
                result += ',';
                result += Newick(m_tree.GetNode(node.Children()[i]), level+1);
            }
            result += ')';
        }
        // Append current
        auto taxonomic_id = NodeId::GetTaxonomicId(node.GetNodeId());
        auto gene_id = NodeId::GetGeneId(node.GetNodeId());
        std::string node_label = std::to_string(taxonomic_id);
        if (gene_id) node_label += "_" + std::to_string(gene_id);
        node_label += "_hits=";
        node_label += node.ToString();
//        node_label += ']';
        result += node_label;
        return result;
    }



    //########################################################################
    // LookupResult
    //########################################################################

    LookupResult::LookupResult(const class Lineage &lineage, TID taxonomic_id, GID gene_id, GPOS gene_position, GPOS read_position, bool second) :
            m_lineage(lineage),
            m_taxonomic_id(taxonomic_id),
            m_gene_id(gene_id),
            m_hit{taxonomic_id, gene_id} {
        m_hit.GetHit(second).Set(gene_position, read_position);
    }

    LookupResult::LookupResult(const class Lineage &lineage, TID taxonomic_id, GID gene_id, bool second) :
            m_lineage(lineage),
            m_taxonomic_id(taxonomic_id),
            m_gene_id(gene_id),
            m_hit{taxonomic_id, gene_id} {
        m_hit.GetHit(second).Increment();
    }

    LookupResult::LookupResult(const class Lineage &lineage, TID taxonomic_id, GID gene_id) :
            m_lineage(lineage),
            m_taxonomic_id(taxonomic_id),
            m_gene_id(gene_id),
            m_hit({taxonomic_id, gene_id}) {
    }

    Hit &LookupResult::GetHit(bool second) {
        return m_hit.GetHit(second);
    }

    HitPE &LookupResult::GetHitPE() {
        return m_hit;
    }

    const HitPE &LookupResult::GetHitPE() const {
        return m_hit;
    }

    const Lineage &LookupResult::GetLineage() const {
        return m_lineage;
    }

    const GID &LookupResult::GeneId() const {
        return m_gene_id;
    }

    const TID &LookupResult::GetTaxonomicId() const {
        return m_taxonomic_id;
    }

    const std::string LookupResult::ToString() const {
        return m_lineage.ToString() + "\t\t" + std::to_string(m_taxonomic_id) + "\t" + std::to_string(m_gene_id) + "\t\t1: " + m_hit.GetHit(false).ToString() + "\t\t2: " + m_hit.GetHit(true).ToString();
    }

    PatternHandler::PatternHandler(Taxonomy &taxonomy, std::string shape_path, bool sensitive_mode) :
            taxonomy_(taxonomy),
            sensitive_mode_(sensitive_mode) {

        hitmiss_.resize(INIT_HITMISS_SIZE);
        patterns_.resize(INIT_PATTERNS_SIZE);
        cov.resize(10);

        if (!Utils::exists(shape_path)) {
            std::cerr << "File " << shape_path << " does not exist. Abort." << std::endl;
            exit(7);
        }
        snp_detector_.LoadArray(shape_path);

        pattern_len_ = snp_detector_.pattern_size;
    }

    PatternHandler::PatternHandler(Taxonomy &taxonomy, ds::PatternMap *db, bool sensitive_mode) :
            taxonomy_(taxonomy),
            sensitive_mode_(sensitive_mode) {

        hitmiss_.resize(200);
        patterns_.resize(200);
        cov.resize(10);

        pattern_db = db;
        pattern_len_ = db->pattern_size;
    }

    void PatternHandler::Reset() {
        patterns_[0] = 0;
        current_pos_ = 0;
        first_hit_ = INT_MAX;
        last_hit_ = 0;
        std::fill(cov.begin(), cov.end(), 0);
    }

    void PatternHandler::AddNode(int taxid) {
        last_size = hitmiss_.size();

        hitmiss_[current_pos_].Reset();
        hitmiss_[current_pos_++].taxid = taxid;


        if (current_pos_ == hitmiss_.size()) {
            hitmiss_.resize(current_pos_ * 2);
        }
    }

    void PatternHandler::AddLeaf(int taxid, int geneid, int genepos) {
        last_size = hitmiss_.size();

        if (current_pos_ >= hitmiss_.size()) {
            std::cout << "" << (size_t)this << std::endl;
            std::cerr << "current_pos_ >= hitmiss_.size()  " << current_pos_ << " >= " << hitmiss_.size() << std::endl;
            exit(0);
        }

        hitmiss_[current_pos_].taxid = taxid;
        hitmiss_[current_pos_].geneid = geneid;
        hitmiss_[current_pos_++].genepos = genepos;


        if (current_pos_ == hitmiss_.size()) {
            hitmiss_.resize(current_pos_ * 2);
        }
    }

    void PatternHandler::AddMiss() {
        hitmiss_[current_pos_].Reset();
        current_pos_++;


        if (current_pos_ == hitmiss_.size()) {
            hitmiss_.resize(current_pos_ * 2);
        }
    }

    void PatternHandler::Add(size_t pos, bool hit) {
        if (pos/cov_bucket_size > cov.size()) {
            cov.resize(pos/cov_bucket_size);
        }



        cov[pos / cov_bucket_size] += hit;
        if (pos < pattern_len_) {
            patterns_[0] |= ((uint64_t) hit) << (63 - pos);
        } else {
            if ( (pos - pattern_len_ + 1) >= patterns_.size()) {
                patterns_.resize(patterns_.size() * 2);
            }
            if (pos - pattern_len_ + 1 >= patterns_.size() || pos - pattern_len_ + 1 < 0) {
                std::cout << (pos - pattern_len_ + 1) << std::endl;
                std::cout << patterns_.size() << std::endl;
                exit(8);
            }
            patterns_[pos - pattern_len_ + 1] = (patterns_[pos - pattern_len_] << 1) | ((uint64_t) hit << (64 - pattern_len_));
        }
    }

    double PatternHandler::GetHitCov() {
        int hit = 0;

        for (int i = 0; i < 10; i++) {
            hit += cov[i];
        }

        return (double) hit / (read_length_ - snp_detector_.shape_size + 1);
//        return (double) hit / (read_length_ - pattern_db->shape_size + 1);
    }

    int PatternHandler::Evaluate2(int taxid, int geneid, int abs_pos, bool forward, size_t shape_length, size_t read_length,
                                  int start, int n) {

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
                auto mut = ds::SNP{mutation_pos, gene_pos };
                mutations_.insert(mut);
            }
        }

        return mutations_.size();
    }

    int PatternHandler::Evaluate3(ClassificationResult &result, size_t sequence_length, size_t shape_size, size_t gene_length, LineageTIDMap lineages) {
        constexpr bool debug = false;

//        std::string lookup_str = "";

        mutations_.clear();
        auto& hit = *result.m_hit;
        auto &&[offset, hits, forward] = hit.GetOffset().GetHits(shape_size, sequence_length);
        auto [start, n] = OverlappingKmers(forward, offset, sequence_length,
                                             shape_size, gene_length);


//        if (result.m_result_kmers > n) {
//            std::cerr << "_________________" << std::endl;
//            std::cerr << "result_kmers:  " << result.m_result_kmers << std::endl;
//            std::cerr << "overlap:       " << n << std::endl;
//            std::cerr << "genelength:    " << gene_length << std::endl;
//            std::cerr << "genepos:       " << offset << std::endl;
//            std::cerr << "forward:       " << forward << std::endl;
//            std::cerr << result.m_hit->ToString() << std::endl;
//            std::cerr << std::endl;
//        }

        result.m_gene_pos = offset;
        result.m_forward = forward;

        auto ov_start = std::max(0, result.m_gene_pos);
        auto ov_end = std::min(gene_length, result.m_gene_pos + sequence_length);
        result.m_overlap = ov_end - ov_start;

        //TODO Switch loop from 0 -> end to start -> start + n as this is the portion that is interesting to us
        if constexpr(debug) {
            std::cout << "\n______Evaluate3(" << sequence_length << ", " << shape_size << ", " << gene_length << ")" << std::endl;
            std::cout << "Start: " << start << "  N: " << n << std::endl;
            std::cout << "current_pos_: " << current_pos_ << std::endl;
            std::cout << "Range of hit: [ " << offset << ", " << (offset + sequence_length) << " ]" << (forward ? ">>" : "<<") << std::endl;
        }


        TID taxonomic_id = result.m_taxonomic_id;
        GID gene_id = result.m_gene_id;
        Lineage& lineage = lineages[taxonomic_id];


        int rel_pos = forward ? offset : offset + sequence_length - shape_size;
        auto total_hits = 0;
        for (int pos = 0; pos < current_pos_; pos++) {
            // Iterate all k-mer lookups
            auto& lookup = hitmiss_[pos];

            if constexpr(debug) {
                std::cout << lookup.ToString() << std::endl;
            }

            // Does that work?
            if (lookup.IsMiss()) {
//                if constexpr(debug) lookup_str += '0';
                Add(pos , false);
                continue;
            }

            if (hitmiss_.size() > 200) {
                std::cout << "asdsad" << std::endl;
                exit(9);
            }

            bool leaf_match = (lookup.IsGeneSet() && lookup.geneid == gene_id && lookup.taxid == taxonomic_id && (rel_pos == (lookup.genepos + ((forward * -1) + (!forward * 1)) * pos)));


//            bool hit = lookup.taxid >= 0 && // First requirement taxid must not be -1
//                       (leaf_match || // Second, geneid and genepos may not be 0
//                        (!lookup.IsGeneSet() && (taxonomic_id == lookup.taxid || taxonomy_.IsNodeAncestor(lookup.taxid, taxonomic_id))));

            bool hit = lookup.taxid >= 0 && // First requirement taxid must not be -1
                       (leaf_match || // Second, geneid and genepos may not be 0
                        (!lookup.IsGeneSet() && (taxonomic_id == lookup.taxid || lineage.IsIdAncestor(lookup.taxid))));

            assert(lineage.IsIdAncestor(lookup.taxid) == taxonomy_.IsNodeAncestor(lookup.taxid, taxonomic_id));

            if constexpr(debug) {
//                if (hit) {
//                    lookup_str += leaf_match ? '1' : '.';
//                } else {
//                    lookup_str += 'O';
//                }
            }

            if (leaf_match || (sensitive_mode_ && hit)) {
                if (first_hit_ == INT_MAX) first_hit_ = pos;
                last_hit_ = pos - pattern_len_ + 1;
            }

            if (pos >= start && pos < start+n) total_hits++;
            Add(pos , hit);
        }

        if constexpr(debug) {
//            std::cout << "full:  " << lookup_str << " " << total_hits << std::endl;
//            std::cout << "short: " << lookup_str.substr(start, n) << " " << total_hits << std::endl;
        }
//        for (int i = first_hit_; i <= last_hit_; i++) {

        result.SetTotalKmers(n);
        result.SetKmerHits(total_hits);

        if (total_hits > result.m_overlap) {
            std::cerr << "\n\nOverlap:         " << result.m_overlap << " > " << total_hits << std::endl;
            std::cerr << "Start:         " << start << std::endl;
            std::cerr << "n:             " << n << std::endl;
            std::cerr << "gene_length:   " << gene_length << std::endl;
            std::cerr << "_________________" << std::endl;
            std::cerr << "result_kmers:  " << result.m_result_kmers << std::endl;
            std::cerr << "overlap:       " << n << std::endl;
            std::cerr << "genelength:    " << gene_length << std::endl;
            std::cerr << "genepos:       " << offset << std::endl;
            std::cerr << "forward:       " << forward << std::endl;
            std::cerr << result.m_hit->ToString() << std::endl;
            std::cerr << std::endl;
        }

        if (n < pattern_db->pattern_size) {
//            std::cout << "pattern_db->pattern_size   " << n << ", " << pattern_db->pattern_size << std::endl;
//            std::cerr << "ERROR: " << n << "," << total_hits << "\tGL: " << gene_length  << " " << offset << "," << hits  << "," << forward << "\t\t" << (result.m_hit ? result.m_hit->ToString() : "") << std::endl;
//
//            result.SetTotalKmers(n);
//            result.SetKmerHits(total_hits);
            return 0;
        }

        auto total_patterns = n - pattern_db->pattern_size + 1;
        for (int i = start; i < start + total_patterns; i++) {
            uint32_t pattern = patterns_[i] >> 32;

            if constexpr(debug) {
                std::cout << "pattern(" << i << ")\t" << std::bitset<32>(pattern).to_string().substr(0, pattern_db->pattern_size);
            }

            if (!patterns_[i]) {
                if constexpr(debug) {
                    std::cout << std::endl;
                }
                continue;
            }

            auto& mutation = pattern_db->map_array[pattern >> (32llu - pattern_db->pattern_size)];

//            if (!mutation.Empty()) Utils::Input();
            if constexpr(debug) {
                std::cout << " " << mutation.ToString() << " (" << mutation.Size() << ")\t";
            }

            for (int m = 0; m < mutation.Size(); m++) {
                int mutation_pos = mutation.mutations[m] + i;
                if constexpr(debug) {
                    std::cout << mutation_pos << ",";
                }
                int gene_pos = forward ? rel_pos + mutation_pos : rel_pos - mutation_pos + shape_size - 1;
                auto mut = ds::SNP{ mutation_pos, gene_pos };
                mutations_.insert(mut);
            }

            if constexpr(debug) {
                std::cout << std::endl;
            }
        }
        if constexpr(debug) {
            std::cout << "Mutation count: " << mutations_.size() << std::endl;
        }
        return mutations_.size();
    }

    int PatternHandler::Evaluate(int taxid, int geneid, int rel_pos, bool forward, int shape_length) {
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

        for (int i = first_hit_; i <= last_hit_; i++) {
            if (!patterns_[i]) continue;
            uint32_t pattern = patterns_[i] >> 32;

            auto& mutation = pattern_db->map_array[pattern >> (32llu - pattern_db->pattern_size)];

            for (int m = 0; m < mutation.Size(); m++) {
                int mutation_pos = mutation.mutations[m] + i;
                int gene_pos = forward ? rel_pos + mutation_pos : rel_pos - mutation_pos + shape_length - 1;

                auto mut = ds::SNP{mutation_pos, gene_pos };
                this->mutations_.insert(mut);
            }
        }

        return mutations_.size();
    }

//    void LineFromHitPE(IO::ClassificationLine &line, const HitPE &hit_pe, size_t record_id, FastxRecord &record1,
//                       FastxRecord &record2, size_t shape_size, size_t total_hits) {
//
//    }
    void ClassificationResult::SetNodeHit(TID taxonomic_id, size_t kmer_hits_result) {
        m_taxonomic_id = taxonomic_id;
        m_result_kmers = kmer_hits_result;
        m_success = true;
    }

    void ClassificationResult::SetLeafHit(Hit* hit, GID gene_id) {
        assert(m_taxonomic_id != 0);
        m_hit = hit;
        m_gene_id = gene_id;
    }

    void ClassificationResult::SetTotalKmers(size_t total_kmers) {
        m_total_kmers = total_kmers;
    }

    void ClassificationResult::SetKmerHits(size_t kmer_hits) {
        m_result_kmers = kmer_hits;
    }
}