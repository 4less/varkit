//
// Created by joachim on 20/07/2020.
//

#pragma once

#include <vector>
#include <fstream>
#include <c++/7/unordered_map>
#include <c++/7/bits/locale_facets.tcc>
#include "Utils.h"
#include "csv.h"

class TaxonomyInterface {
public:
    virtual int lca(int a, int b) = 0;
};

class TaxonomyNode {
public:
    int id = -1;
    TaxonomyNode * parent;
    std::string rank;
    int level;
    std::vector<TaxonomyNode*> children;
    
    void print() {
        cout << "id:\t" << id << endl;
        if (parent)
            cout << "pid:\t" << parent->id << endl;
        else cout << "root" << endl;
        cout << "level:\t" << level << endl;
    }
    
    ~TaxonomyNode() {
    }
};

struct NCBITaxonomy : public TaxonomyInterface {
    std::unordered_map<int, TaxonomyNode*> map_;
    std::unordered_map<int, string> names_;
    std::unordered_map<int,int> to_ncbi_;
    std::unordered_map<int,int> from_ncbi_;
    
    bool custom = false;
    
    inline static const string SCIENTIFIC_NAME = "scientific name";
    
    NCBITaxonomy(std::string path) {
        loadNodes(path);
    }
    
    NCBITaxonomy(){};
    
    
    void loadNodes(std::string path) {
        int count = 0;
        
        auto str = csv::load_file(path.c_str());
        auto parser = csv::make_parser( str , '\t');
        
        for (auto&& row : parser ) {
            auto it = row.begin();
            int id = (*it).to_int();
            int parent_id = (*(++++it)).to_int();
            string rank = (*(++++it)).to_string();
    
            TaxonomyNode* node = getOrNew(id);
            node->id = id;
            if (id != 1) {
                node->rank = rank;
                TaxonomyNode* parent = getOrNew(parent_id);
                node->parent = parent;
                parent->id = parent_id;
                parent->children.push_back(node);
            } else {
                node->rank = "root";
                node->parent = nullptr;
                node->level = 0;
        
            }
        }
        relevel(map_.at(1));
    }
    
    void loadNames(std::string path) {
        auto str = csv::load_file(path.c_str());
        auto parser = csv::make_parser( str , '\t');
        for (auto&& row : parser ) {
            auto it = row.begin();
            
            int id = (*it).to_int();

            if (names_.find(id) != names_.end()) continue;
            
            string name = (*(++++it)).to_string();
            string type = (*(++++++++it)).to_string();
    
            if (type == SCIENTIFIC_NAME)
                names_.insert({id, name});
        }
    }

    void subsetByTaxa(string path) {
        vector<int> taxa;
        
        ifstream taxa_in;
        string line;
        taxa_in.open(path);
        if (taxa_in.is_open()) {
            while (getline(taxa_in, line)) {
                taxa.push_back(stoi(line));
            }
        }
        taxa_in.close();
        
        auto taxa_set = getCompleteSetOfTaxa(taxa);
    
    
        cout << "set has 1315283: " << (taxa_set.find(1315283) != taxa_set.end()) << endl;
        
        int root_id = -1;
        for (auto id : taxa_set) {
            if (root_id == -1) root_id = id;
            else root_id = lca(root_id, id);
        }
        cout << "root: " << root_id << endl;
        
        TaxonomyNode* root = map_.at(root_id);
        root->parent = nullptr;
    
        cout << "before worker size: " << map_.size() << endl;
        subsetWorker(root, taxa_set);
        cout << "after worker size: " << map_.size() << endl;
        deleteUnlinkedNodes(root);
        cout << "after delete unlinked size: " << map_.size() << endl;
        
        cout << "map has 1315283: " << (map_.find( 1315283) != map_.end()) << endl;
        
        int rel = 0;
        int ire = 0;
        countRelevantNodesInSubtree(root, rel, ire);
        cout << "relevant: " << rel << endl;
        cout << "irrelevant: " << ire << endl;
    
        auto children = root->children;
        for (auto n : children) {
            cleanUp(n);
        }
        relabel(root);
        relevel(root);
        custom = true;
    }
    
    void subsetWorker(TaxonomyNode* node, unordered_set<int> &taxa) {
        auto children = node->children;
        for (auto child : children) {
            if (taxa.find(child->id) == taxa.end()) {
                deleteFromVector(node->children, child);
            }
        }
        for (auto child : node->children) {
            subsetWorker(child, taxa);
        }
    }
    
    unordered_set<int> getCompleteSetOfTaxa(vector<int> &taxa) {
        unordered_set<int> t_set;
        int root_id = taxa[0];
        for (int i = 0; i < taxa.size(); i++) {
            root_id = lca(root_id, taxa[i]);
            t_set.insert(taxa[i]);
            for (int j = i+1; j < taxa.size(); j++) {
                t_set.insert(lca(taxa[i], taxa[j]));
            }
        }
        auto copy = t_set;
        for (auto taxid : copy) {
            if (taxid == root_id) continue;
            auto node = map_.at(taxid)->parent;
            while (node && node->id != root_id && t_set.find(node->id) == t_set.end()) {
                t_set.insert(node->id);
                auto cache = node;
                node = node->parent;
            }
        }
        return t_set;
    }
    
    void print(TaxonomyNode* node) {
        cout << tabs(node->level-1) << node->id;
        if (!names_.empty())
            cout << ": " << names_.at(getNCBI(node->id));
        if (node->parent)
            cout << " - pid: " << node->parent->id;
        cout << endl;
        
        for (auto child : node->children)
            print(child);
    }
    
    void deleteUnlinkedNodes(TaxonomyNode* node) {
        unordered_set<int> valid_taxa;
        unordered_set<int> invalid_taxa;
        extractTaxaInSubtree(node, valid_taxa);
        for (auto pair : map_) {
            if (valid_taxa.find(pair.first) == valid_taxa.end()) {
                invalid_taxa.insert(pair.first);
            }
        }
        
        for (auto id : invalid_taxa) {
            auto cache = map_.at(id);
            map_.erase(id);
            delete(cache);
        }
    }
    
    string tabs(int count) {
        string tabs = "";
        for (int i = 0; i < count; i++)
            tabs += '\t';
        return tabs;
    }
    
    void relevel(TaxonomyNode* node, int level = 0) {
        node->level = ++level;
        for (auto child : node->children)
            relevel(child, level);
    }
    
    void relabel(TaxonomyNode* root) {
        int custom_id = 0;
        int ncbi_id;
        queue<TaxonomyNode*> q;
        q.push(root);
        map_.clear();
        
        int count = 0;
        while (!q.empty()) {
            count++;
            auto node = q.front();
            q.pop();
            
            int ncbi_id = node->id;
            node->id = ++custom_id;
            
            map_.insert({custom_id, node});
            // fill conversion maps
            to_ncbi_.insert({custom_id,ncbi_id});
            from_ncbi_.insert({ncbi_id, custom_id});
    
            for (auto child : node->children) {
                q.push(child);
            }
        }
    }
    
    void countRelevantNodesInSubtree(TaxonomyNode* node, int &sum, int &ir) {
        if (node->children.size() == 1) {
            ir++;
            countRelevantNodesInSubtree(node->children[0], sum, ir);
        } else {
            sum++;
            for (auto n : node->children)
                countRelevantNodesInSubtree(n, sum, ir);
        }

    }
    
    void cleanUp(TaxonomyNode* node) {
        auto tmp = node;
    
        // if node has only one child it is an unnecessary intermediate node
        if (tmp->children.size() == 1) {
            while (tmp->children.size() == 1)
                tmp = tmp->children[0];
            
            // rewire new parent and new child
            deleteFromVector(node->parent->children, node);
            node->parent->children.push_back(tmp);
            tmp->parent = node->parent;
            
            //delete all intermediate nodes
            while (node != tmp) {
                auto cache = node;
                node = node->children[0];
                
                //delete from map and finally delete node
                map_.erase(cache->id);
                delete cache;
            }
        }
        
        auto children = tmp->children;
        for (auto n : children) {
            cleanUp(n);
        }
    }
    
    void checkValidity(TaxonomyNode* node) {
        if (map_.find(node->id) == map_.end()) {
            cout << "danger" << endl;
            node->print();
            exit(0);
        }
        for (auto n : node->children)
            checkValidity(n);
    }
    
    static bool compareNodesByLevel(TaxonomyNode* n1, TaxonomyNode* n2)
    {
        return (n1->level < n2->level);
    }

    int subset(vector<int> taxids) {
        int root_id;
        if (taxids.size() == 1) {
            root_id = taxids[0];
            auto root = map_.at(root_id);
            deleteTreeUpstream(root);
            return root_id;
        }
    
        root_id = taxids[0];
        for (int i = 1; i < taxids.size(); i++)
            root_id = lca(root_id, taxids[i]);
        
        auto root = map_.at(root_id);
        
        deleteTreeUpstream(root);
        
        vector<TaxonomyNode*> nodes;
        for (auto id : taxids) nodes.push_back(map_.at(id));
        std::sort(nodes.begin(), nodes.end(), compareNodesByLevel);
        
        for (auto node : nodes) {
            cout << node->id << " " << node->level << endl;
        }
        
        // remove taxid that are descendants of other taxa in the list
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = i+1; j < nodes.size(); j++) {
                if (nodes[i]->level == nodes[j]->level) continue;
                if (isAncestorOf(nodes[i], nodes[j]->id)) {
                    deleteFromVector(nodes, nodes[j]);
                    j--;
                }
            }
        }
        
        // only keep subtrees of taxons in the list
        deleteWorker(root, nodes);
        
        auto children = root->children;
        for (auto n : children) {
            cleanUp(n);
        }
        
        return root->id;
    }
    
    void extractTaxaInSubtree(TaxonomyNode* node, unordered_set<int> &taxa) {
        taxa.insert(node->id);
        for (auto child : node->children)
            extractTaxaInSubtree(child, taxa);
    }
    
    void deleteWorker(TaxonomyNode* node, vector<TaxonomyNode*> &nodes) {
        bool isRelevant = false;
        for (auto n : nodes) {
            if (node == n) return;
            if (isAncestorOf(node, n->id)) {
                isRelevant = true;
                break;
            }
        }
        
        if (isRelevant) {
            // iterate over copy of children object to ensure not to modify the children vector
            // while accessing it
            auto children = node->children;
            
            for (auto n : children)
                deleteWorker(n, nodes);
        } else {
            // remove from parent
            if (node->parent)
                deleteFromVector(node->parent->children, node);
            
            // remove subtree
            deleteSubtree(node);
        }
    }
    
    void deleteTreeUpstream(TaxonomyNode* newroot) {
        TaxonomyNode* tmp = newroot;
        
        while (tmp->parent) {
            //cout << "delete siblings: "  << tmp->id << endl;
            for (auto subt : tmp->parent->children) {
                if (subt == tmp) continue;
                deleteSubtree(subt);
            }
            auto cache = tmp;
            tmp = tmp->parent;
            if (cache != newroot) {
                map_.erase(cache->id);
                delete cache;
            }
            else cache->parent = nullptr;
        }
    }
    
    void deleteSubtree(TaxonomyNode* node) {
        if (!node->children.empty()) {
            for (auto child : node->children)
                deleteSubtree(child);
        }
        //remove from map
        map_.erase(node->id);
        delete node;
    }
    
    void deleteAsChild(TaxonomyNode* node) {
        TaxonomyNode* parent = node->parent;
        if (parent) {
            deleteFromVector(parent->children, node);
            if (parent->children.empty()) {
                deleteAsChild(parent);
                if (map_.find(parent->id) != map_.end()) {
                    map_.erase(parent->id);
                    delete parent;
                }
            }
        }
    }
    
    void clear() {
        map_.clear();
        from_ncbi_.clear();
        to_ncbi_.clear();
    }
    
    void deleteFromVector(vector<TaxonomyNode*> &v, TaxonomyNode* rm) {
        if (std::find(v.begin(), v.end(), rm) != v.end())
            v.erase(std::remove(v.begin(), v.end(), rm), v.end());
    }
    
    bool isAncestorOf(TaxonomyNode* node, int i) {
        auto tmp = map_.at(i);
        if (node->level > tmp->level) return false;
        while (tmp->level > node->level)
            tmp = tmp->parent;
        return tmp->id == node->id;
    }
    
    bool isDescendantOf(TaxonomyNode* node, int i) {
        auto tmp = map_.at(i);
        if (node->level < tmp->level) return false;
        while (tmp->level < node->level)
            node = node->parent;
        
        cout << endl << "tmp: " << tmp->level << endl;
        cout << "tmp: " << tmp->id << endl;
        cout << "nod: " << node->level << endl;
        cout << "nod: " << node->id << endl;
        if (node->parent) {
            cout << "nodparent: " << node->parent->id << endl;
            cout << "nodparent: " << node->parent->level << endl;
        }
    
        return tmp->id == node->id;
    }
    
    inline TaxonomyNode * getOrNew(int taxid) {
        if (map_.count(taxid) == 1) {
            return map_.at(taxid);
        } else {
            auto node = new TaxonomyNode();
            map_.insert({taxid, node});
            return node;
        }
    }


public:
    int lca(int t1, int t2) override {
        if (custom) {
            t1 = from_ncbi_.at(t1);
            t2 = from_ncbi_.at(t2);
        }
        auto node1 = map_.at(t1);
        auto node2 = map_.at(t2);
        
        int min = std::min(node1->level, node2->level);
        
        while (node1->level > min)
            node1 = node1->parent;
        while (node2->level > min)
            node2 = node2->parent;
        
        while (node1 != node2) {
            node1 = node1->parent;
            node2 = node2->parent;
        }
        return node1->id;
    }
    
    void shrinkToSubtreesOf(vector<int> taxids) {
        int root = subset(taxids);
        
        //clean up
        cout << "clean up" << endl;
        deleteUnlinkedNodes(map_.at(root));
        
        //relabel
        cout << "relabel" << endl;
        relabel(map_.at(root));
        
        //relevel
        cout << "relevel" << endl;
        relevel(map_.at(1));
    }
    
    void loadCustomNames(string path) {
    
        auto str = csv::load_file(path.c_str());
        auto parser = csv::make_parser(str);
    
        bool header = true;
    
        for (auto&& row : parser ) {
            if (header) {
                header = false;
                continue;
            }
        
            auto it = row.begin();
        
            int taxid = (*it).to_int();
            string scientific_name = (*(++it)).to_string();
            names_.insert({taxid, scientific_name});
        }
        custom = true;
    }
    
    void loadCustomNodes(string path) {
        clear();
        
        cout << path << endl;
        
        auto str = csv::load_file(path.c_str());
        auto parser = csv::make_parser(str);
    
        bool header = true;
        
        for (auto&& row : parser ) {
            if (header) {
                header = false;
                continue;
            }
            
            auto it = row.begin();
            
            int taxid = (*it).to_int();
            int ncbi_taxid = (*(++it)).to_int();
            
            to_ncbi_.insert({taxid, ncbi_taxid});
            from_ncbi_.insert({ncbi_taxid, taxid});
            
            int parent_id = (*(++it)).to_int();
            string rank = (*(++it)).to_string();
            int level = (*(++it)).to_int();
        
            TaxonomyNode* node = getOrNew(taxid);
            node->id = taxid;
            node->rank = rank;
            node->level = level;
            if (taxid != 1) {
                TaxonomyNode* parent = getOrNew(parent_id);
                node->parent = parent;
                parent->id = parent_id;
                parent->children.push_back(node);
            } else {
                node->parent = nullptr;
            }
        }
        custom = true;
    }
    
    int getNCBI(int custom) {
        return to_ncbi_.at(custom);
    }
    
    int getCustom(int ncbi) {
        return from_ncbi_.at(ncbi);
    }
    
    TaxonomyNode* getNode(int custom) {
        if (map_.find(custom) != map_.end())
            return map_.at(custom);
        return nullptr;
    }
    
    void saveCustomNodes(string path) {
        std::ofstream out(path);
        
        //auto err = map_.at(1734091392);
        
        
        if (!out)
            cerr << "unable to open outfile " << path << endl;
        
        auto writer = csv::make_writer(out);
        
        writer.write_row("taxid","ncbi_taxid","parent_taxid","rank","level");
        for (auto pair : map_) {
            auto node = pair.second;
            
            writer.write_row(
                    node->id,
                    to_ncbi_.at(node->id),
                    (node->parent ? to_string(node->parent->id) : ""),
                    node->rank,
                    node->level);
        }
        out.close();
    }
    
    void subsetNames() {
        auto copy = names_;
        for (auto row : copy) {
            if (from_ncbi_.find(row.first) == from_ncbi_.end())
                names_.erase(row.first);
        }
    }
    
    void saveCustomNames(string path) {
        std::ofstream out(path);
        if (!out)
            cerr << "unable to open outfile " << path << endl;
    
        auto writer = csv::make_writer(out);
        
        writer.write_row("taxid","scientific_name");
        for (auto pair : names_) {
            writer.write_row(
                    pair.first,
                   pair.second);
        }
        out.close();
    }
};

