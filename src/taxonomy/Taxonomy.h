//
// Created by joachim on 20/07/2020.
//

#pragma once

#include <vector>
#include <fstream>
//#include <c++/7/unordered_map>
//#include <c++/7/bits/locale_facets.tcc>
#include "Utils.h"
#include "csv.h"
#include <unordered_set>
#include <unordered_map>
#include <queue>

using namespace std;

class TaxonomyNode {
public:
    int id = -1;
    TaxonomyNode * parent;
    std::string rank;
    std::string name = "";
    int level;
    std::vector<TaxonomyNode*> children;
    
    void print() {
        std::cout << "id:\t" << id << endl;
        if (!name.empty())
            cout << "name:\t" << name << endl;
        if (parent) {
            cout << "pid:\t" << parent->id << endl;
            cout << "pname:\t" << parent->name << endl;
        } else cout << "root" << endl;
        cout << "level:\t" << level << endl;
        cout << "child_count:\t" << children.size() << endl;
    }
    
    ~TaxonomyNode() {
//        cout << "destroy TaxonomyNode" << endl;
    }
};

class TaxonomyInterface {
public:
    virtual int lca(int a, int b) = 0;
    virtual int getCustom(int t) = 0;
    virtual int getCustom(std::string t) = 0;
    virtual std::string getOriginal(int taxid) = 0;
    
    virtual void loadCustomNodes(std::string basicString) = 0;
    virtual void loadCustomNames(std::string path) = 0;
    
    virtual bool hasNode(int i) = 0;
    virtual TaxonomyNode* getNode(int i) = 0;
};

class GTDBTaxonomy : public TaxonomyInterface {
public:
    GTDBTaxonomy() {};
    ~GTDBTaxonomy() {
//        cout << "destroy GTDBTaxonomy" << endl;
    }
    
    
    std::string getOriginal(int taxid) override {
        return getNode(taxid)->name;
    }
    
    
    int lca(int t1, int t2) override {
        
        if (map_.find(t1) == map_.end() || map_.find(t2) == map_.end()) {
            cerr << t1 << " or " << t2 << " are unknown custom taxids in the map at lca(t1, t2)" << endl;
        }
        
        auto node1 = map_.at(t1);
        auto node2 = map_.at(t2);
        
        int min = std::min(node1->level, node2->level);
        
        while (node1->level > min) {
            node1 = node1->parent;
        }
        
        while (node2->level > min) {
            node2 = node2->parent;
        }
        
        while (node1 != node2) {
            node1 = node1->parent;
            node2 = node2->parent;
        }
        return node1->id;
    }
    
    std::string lca(std::string a, std::string b) {
        if (name_map_.find(a) == name_map_.end() || name_map_.find(b) == name_map_.end()) {
            cerr << a << " or " << b << " are unknown taxids in the map at lca(a, b)" << endl;
        }

        auto node1 = name_map_.at(a);
        auto node2 = name_map_.at(b);
    
        int min = std::min(node1->level, node2->level);
    
        while (node1->level > min)
            node1 = node1->parent;
        while (node2->level > min)
            node2 = node2->parent;
    
        while (node1 != node2) {
            node1 = node1->parent;
            node2 = node2->parent;
        }
        return node1->name;
    }

    std::string lca(TaxonomyNode* node1, TaxonomyNode* node2) {
        int min = std::min(node1->level, node2->level);

        while (node1->level > min)
            node1 = node1->parent;
        while (node2->level > min)
            node2 = node2->parent;

        while (node1 != node2) {
            node1 = node1->parent;
            node2 = node2->parent;
        }
        return node1->name;
    }
    
    int getCustom(int t) override {
        if (from_gtdb_.find(to_string(t)) != from_gtdb_.end())
            return from_gtdb_.at(to_string(t));
        else {
            cerr << "GTDB: " << t << " missing at getCustom(int gtdb)" << endl;
            return -1;
        }
    }
    
    int getCustom(std::string t) override {
        return from_gtdb_.at(t);
    }
    
    
    
    void loadNodes(std::string path) {
        int count = 0;
        
        if (!Utils::exists(path)) {
            cerr << "Node file " << path << " does not exist." << endl;
            exit(0);
        }
        
        auto str = csv::load_file(path.c_str());
        auto parser = csv::make_parser( str , '\t');
        
        std::string gtdb = "";
        int id = 0;
        int parent_id = 0;
        std::string name = "";
        std::string parent_name = "";
        
        std::string lineage = "";
        std::string rank = "";
        
        bool header = true;
    
        TaxonomyNode* root = getOrNew("root");
        root->id = 0;
        
        int row_count = 0;
        
        for (auto&& row : parser ) {
            if (header) {
                header = false;
                continue;
            }
            row_count++;
            
            auto it = row.begin();
            gtdb = (*it).to_string();


            if (gtdb.substr(0, 10).compare("GUT_GENOME") == 0) {
                id = stoi(gtdb.substr(10));
            }
            if (gtdb.substr(0, 10).compare("MGYG-HGUT-") == 0) {
                id = stoi(gtdb.substr(10));
            }
    
            TaxonomyNode* node;
            TaxonomyNode* parent;
            int i = 0;
            
            lineage = (*(++it)).to_string();
            std::vector<std::string> lineage_v;
            Utils::split(lineage_v, lineage, ";");


//            if (strcmp("MGYG-HGUT-04643", gtdb.c_str()) == 0) {
//                std::cout << "MGYG-HGUT-04643 yes" << std::endl;
//                std::cout << lineage << std::endl;
//                exit(0);
//            }

            for (i = 0; i < lineage_v.size(); i++) {
                string token = lineage_v[i];
                if (token.length() == 3) {
//                    if (id_string.empty()) {
//                        id_string = Utils::split(lineage_v[i-1], "__")[1];
//                    }
                    lineage_v[i] = lineage_v[i] + gtdb;
                }
                
                node = getOrNew(lineage_v[i]);
                parent = i == 0 ? root : getOrNew(lineage_v[i-1]);
                node->parent = parent;

                if (std::find(parent->children.begin(), parent->children.end(), node) == parent->children.end()) {
                    parent->children.push_back(node);
                }
            }
            //node = getOrNew(to_string(id));
            //parent = getOrNew(lineage_v[--i]);
            //node->parent = parent;
            //node->id = id;
            //parent->children.push_back(node);
            gtdb_map_.insert({gtdb, node});
            map_.insert({id, node});
        }
        
        relevel(name_map_.at("root"));

        std::cout << "size: " << map_.size() << std::endl;
        std::cout << "size: " << name_map_.size() << std::endl;
        std::cout << "size: " << gtdb_map_.size() << std::endl;
    }

    TaxonomyNode* get(std::string name) {
        if (name_map_.find(name) != name_map_.end()) {
            return(name_map_.at(name));
        } else if (gtdb_map_.find(name) != gtdb_map_.end()) {
            return(gtdb_map_.at(name));
        }
        return nullptr;
    }
    
    
    void subsetByTaxa(string path) {
        cout << "subsetByTaxa" << endl;
        ifstream taxa_in;
        string line;
        taxa_in.open(path);
        
        if (!taxa_in.is_open()) {
            cerr << "file " << path << " could not be loaded." << endl;
            exit(0);
        }
        string id;
        unordered_set<std::string> taxa_set;
        int lines = 0;
        while (getline(taxa_in, line)) {
            ++lines;

            if (line.substr(0, 10).compare("GUT_GENOME") == 0)
                id = to_string(stoi(line.substr(10,16)));
            if (line.substr(0, 10).compare("MGYG-HGUT-") == 0)
                id = to_string(stoi(line.substr(10,16)));

//            if (id[1] != '_') id = "s__" + line;

//            if (name_map_.find(id) != name_map_.end()) {
            if (gtdb_map_.find(line) != gtdb_map_.end()) {
                auto name = gtdb_map_.at(line)->name;
                if (strcmp("s__Collinsella sp002232035", name.c_str()) == 0) {
                    std::cout << "name: " << name << std::endl;
                    std::cout << "gtdb: " << line << std::endl;
                }
                taxa_set.insert(name);
//                std::cout << gtdb_map_.at(line)->name << " <- " << line << std::endl;
            } else {
                cerr << "discard " << line << endl;
            }
//            std::cout << "lines: " << lines << " taxa: " << taxa_set.size() << std::endl;
        }
        cout << lines << endl;
        taxa_in.close();
        
        
        auto orig_set = taxa_set;
        cout << "final set needs to contain at least " << taxa_set.size() << " ids." << endl;

        TaxonomyNode* node1;
        TaxonomyNode* node2;

        std::string root_name = "";

        std::cout << "taxasize: " << taxa_set.size() << std::endl;
        for (auto id : taxa_set) {
            std::cout << "id: " << id << std::endl;
            if (root_name.empty()) root_name = id;
            else {
                root_name = lca(name_map_.at(root_name),name_map_.at(id));
            }
            if (name_map_.find(id) == name_map_.end()) {
                std::cout << "sayonara" << std::endl;
                exit(0);
            }
        }
        cout << "root: " << root_name << endl;


        TaxonomyNode* root = get(root_name);
        root->parent = nullptr;

        cout << "make subset complete" << endl;
        makeSubsetComplete(taxa_set, root_name);
        cout << "is complete check" << endl;
        isSubsetComplete(taxa_set, root_name);



        subsetWorker(root, taxa_set);
        cout << "after worker size: " << name_map_.size() << endl;

        deleteUnlinkedNodes(root);
        cout << "after delete unlinked size: " << name_map_.size() << endl;


        int rel = 0;
        int ire = 0;
        countRelevantNodesInSubtree(root, rel, ire);

        cout << "relevant: " << rel << endl;
        cout << "irrelevant: " << ire << endl;


        relabel(root);
        relevel(root);
//
//
//
        //custom = true;
    }

    void relabel() {
        auto root = name_map_["root"];
        std::cout << root->name << std::endl;

        relabel(root);
        relevel(root);
    }
    
    void saveCustomNodes(string path) {
        std::ofstream out(path);
        
        //auto err = map_.at(1734091392);
        
        
        if (!out)
            cerr << "unable to open outfile " << path << endl;
        
        auto writer = csv::make_writer(out);
        
        writer.write_row("taxid","gtdb_taxid","parent_taxid","rank","level");
        for (auto pair : map_) {
            auto node = pair.second;
            
            writer.write_row(
                    node->id,
                    to_gtdb_.at(node->id),
                    (node->parent ? to_string(node->parent->id) : ""),
                    node->rank,
                    node->level);
        }
        out.close();
    }

    void saveCustomNodes2(string path) {
        std::ofstream out(path);

        //auto err = map_.at(1734091392);


        if (!out)
            cerr << "unable to open outfile " << path << endl;

        auto writer = csv::make_writer(out);

        writer.write_row("taxid","parent_taxid","rank","level");
        for (auto pair : map_) {
            auto node = pair.second;

            writer.write_row(
                    node->id,
                    (node->parent ? to_string(node->parent->id) : ""),
                    node->rank,
                    node->level);
        }
        out.close();
    }
    
    void saveCustomNames(string path) {
        std::ofstream out(path);
        if (!out)
            cerr << "unable to open outfile " << path << endl;
        
        auto writer = csv::make_writer(out);
        
        writer.write_row("taxid","scientific_name");
        for (auto pair : map_) {
            writer.write_row(
                    pair.first,
                    pair.second->name);
        }
        out.close();
    }


    void saveCustomNames2(string path) {
        std::ofstream out(path);
        if (!out)
            cerr << "unable to open outfile " << path << endl;

        auto writer = csv::make_writer(out);

        writer.write_row("taxid","scientific_name","type");
        for (auto pair : map_) {
            writer.write_row(
                    pair.first,
                    pair.second->name,
                    "name");
        }
        out.close();
    }
    
    
    void loadCustomNodes(string path) override {
//        clear();
        
        cout << "custom nodes: " << path << endl;
        
        if (!Utils::exists(path)) {
            cerr << "file " << path << " does not exist." << endl;
            exit(0);
        }
        
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
            std::string gtdb_id = (*(++it)).to_string();
            
            
            to_gtdb_.insert({taxid, gtdb_id});
            from_gtdb_.insert({gtdb_id, taxid});
            
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
        
//        custom = true;
    }
    
    std::string rank(char c) {
        switch(c) {
            case ('r'):
                return "root";
                break;
            case ('p'):
                return "phylum";
                break;
            case ('c'):
                return "class";
                break;
            case ('o'):
                return "order";
                break;
            case ('f'):
                return "family";
                break;
            case ('g'):
                return "genus";
                break;
            case ('s'):
                return "species";
                break;
            default:
                return "genome";
                break;
        }
    }
    
    void loadCustomNames(string path) override {
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
            string name = (*(++it)).to_string();
            auto node = map_.at(taxid);
            node->name = name;
            node->rank = rank(name[0]);
        }
    }
    
    bool hasNode(int t) override {
        return (map_.find(t) != map_.end());
    }
    
    TaxonomyNode* getNode(int custom) override {
        if (map_.find(custom) != map_.end())
            return map_.at(custom);
        else {
            cerr << "id unknown: " << custom << " for getNode(custom)" << endl;
            return nullptr;
        }
    }

private:
    std::unordered_map<std::string, TaxonomyNode*> name_map_;
    std::unordered_map<std::string, TaxonomyNode*> gtdb_map_;
    std::unordered_map<int, TaxonomyNode*> map_;
    std::vector<std::string> names_;
    std::unordered_map<int,std::string> to_gtdb_;
    std::unordered_map<std::string,int> from_gtdb_;
    
    inline TaxonomyNode * getOrNew(int taxid) {
        static int new_count = 0;
        if (map_.count(taxid) == 1) {
            return map_.at(taxid);
        } else {
            auto node = new TaxonomyNode();
            map_.insert({taxid, node});
            return node;
        }
    }
    
    void relevel(TaxonomyNode* node, int level = 0) {
        node->level = ++level;
        for (auto child : node->children)
            relevel(child, level);
    }
    
    inline TaxonomyNode * getOrNew(std::string taxid) {
        static int new_count = 0;
        if (name_map_.find(taxid) != name_map_.end()) {
            return name_map_.at(taxid);
        } else {
            //cout << "new: " << ++new_count << endl;
            auto node = new TaxonomyNode();
            node->name = taxid;
            names_.push_back(taxid);
            name_map_.insert({taxid, node});
            return node;
        }
    }
    
    void makeSubsetComplete(unordered_set<std::string> &taxa, std::string root_name) {
        cout << "after subset complete" << taxa.size() << endl;
        TaxonomyNode* node;
        TaxonomyNode* root = get(root_name);
        taxa.insert(root->name);
        auto copy = taxa;
        for (auto id : copy) {
            node = get(id);
            auto cache = node;
            while (node != root) {
                taxa.insert(node->name);
                node = node->parent;
            }
        }
        cout << "after subset complete" << taxa.size() << endl;
    }
    
    void isSubsetComplete(unordered_set<std::string> &taxa, std::string root_name) {
        TaxonomyNode* node;
        TaxonomyNode* root = get(root_name);
        if (taxa.find(root->name) == taxa.end()) {
            cerr << "violation: " << root->name << "(root) missing in set." << endl;
            exit(0);
        }
        int c = 0;
        for (auto id : taxa) {
            c++;
            node = get(id);
            node = name_map_.at(node->name);
            while (node != root) {
                if (taxa.find(node->name) == taxa.end()) {
                    cerr << "violation: " << node->name << " missing in set." << endl;
                    exit(0);
                }
                node = node->parent;
            }
        }
        std::cout << "c: " << c << std::endl;
    }
    
    void subsetWorker(TaxonomyNode* node, unordered_set<std::string> &taxa) {
        auto children = node->children;
        for (auto child : children) {
            if (taxa.find(child->name) == taxa.end()) {
                deleteFromVector(node->children, child);
            }
        }
        for (auto child : node->children) {
            subsetWorker(child, taxa);
        }
    }
    
    void deleteFromVector(vector<TaxonomyNode*> &v, TaxonomyNode* rm) {
        if (std::find(v.begin(), v.end(), rm) != v.end())
            v.erase(std::remove(v.begin(), v.end(), rm), v.end());
    }
    
    void deleteUnlinkedNodes(TaxonomyNode* node) {
        unordered_set<std::string> valid_taxa;
        unordered_set<std::string> invalid_taxa;
        extractTaxaInSubtree(node, valid_taxa);
        for (auto pair : name_map_) {
            if (valid_taxa.find(pair.first) == valid_taxa.end()) {
                invalid_taxa.insert(pair.first);
            }
        }
        
        for (auto id : invalid_taxa) {
            auto cache = name_map_.at(id);
            name_map_.erase(id);
            delete(cache);
        }
    }
    
    void extractTaxaInSubtree(TaxonomyNode* node, unordered_set<std::string> &taxa) {
        taxa.insert(node->name);
        for (auto child : node->children)
            extractTaxaInSubtree(child, taxa);
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
    
    void cleanUp(TaxonomyNode* node, unordered_set<std::string> &taxa) {
        auto tmp = node;
        
        // if node has only one child it is an unnecessary intermediate node
        if (tmp->children.size() == 1) {
            while (tmp->children.size() == 1 && taxa.find(tmp->name) == taxa.end())
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
                name_map_.erase(cache->name);
                delete cache;
            }
        }
        
        auto children = tmp->children;
        for (auto n : children) {
            cleanUp(n, taxa);
        }
    }
    
    void relabel(TaxonomyNode* root) {
        int custom_id = 0;
        int ncbi_id;
        std::queue<TaxonomyNode*> q;
        q.push(root);
        name_map_.clear();
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
            to_gtdb_.insert({custom_id, node->name});
            from_gtdb_.insert({node->name, custom_id});
        
            for (auto child : node->children) {
                q.push(child);
            }
        }
    }
};


class NCBITaxonomy : public TaxonomyInterface {
    std::unordered_map<int, TaxonomyNode*> map_;
    std::unordered_map<int, string> names_;
    std::unordered_map<int,int> to_ncbi_;
    std::unordered_map<int,int> from_ncbi_;
    
    bool custom = false;
    
    inline static const string SCIENTIFIC_NAME = "scientific name";
    
    NCBITaxonomy(std::string path) {
        loadNodes(path);
    };
    
    
    void printPath(TaxonomyNode* node, TaxonomyNode* root) {
        while (node != root) {
            cout << node->id << " >> ";
            node = node->parent;
        }
        cout << root->id << endl;
    }
    
    void makeSubsetComplete(unordered_set<int> &taxa, int root_id) {
        cout << taxa.size() << endl;
        TaxonomyNode* node;
        TaxonomyNode* root = map_.at(root_id);
        taxa.insert(root->id);
        auto copy = taxa;
        for (auto id : copy) {
            node = map_.at(id);
            auto cache = node;
            while (node != root) {
                taxa.insert(node->id);
                node = node->parent;
            }
        }
    }
    
    void isSubsetComplete(unordered_set<int> &taxa, int root_id) {
        
        TaxonomyNode* node;
        TaxonomyNode* root = map_.at(root_id);
        if (taxa.find(root->id) == taxa.end()) {
            cerr << "violation: " << root->id << "(root) missing in set." << endl;
            exit(0);
        }
        for (auto id : taxa) {
            node = map_.at(id);
            while (node != root) {
                if (taxa.find(node->id) == taxa.end()) {
                    cerr << "violation: " << node->id << " missing in set." << endl;
                    exit(0);
                }
                node = node->parent;
            }
        }
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
        cout << "get comP " << endl;
        unordered_set<int> t_set;
        int root_id = taxa[0];
        cout << "add to set ... ";
        
        for (int i = 0; i < taxa.size(); i++) {
            root_id = lca(root_id, taxa[i]);
            t_set.insert(taxa[i]);

            for (int j = i+1; j < taxa.size(); j++) {
                t_set.insert(lca(taxa[i], taxa[j]));
            }
        }
        cout << "ok" << endl;
        
        /*
        auto copy = t_set;
        for (auto taxid : copy) {
            if (taxid == root_id) continue;
            cout << taxid;
            auto node = map_.at(taxid)->parent;
            cout << endl;
            while (node && node->id != root_id && t_set.find(node->id) == t_set.end()) {
                t_set.insert(node->id);
                auto cache = node;
                node = node->parent;
            }
        }*/
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
    
    void testNames(std::unordered_set<int> taxa) {
        bool err = false;
        for (auto t : taxa) {
            if (names_.find(t) == names_.end()) {
                cerr << "taxon " << t << " is missing in names_." << endl;
                err = true;
            }
            if (from_ncbi_.find(t) == from_ncbi_.end()) {
                cerr << "taxon " << t << " is missing in from_ncbi_." << endl;
                err = true;
            }
            if (err) exit(1);
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
    
    void cleanUp(TaxonomyNode* node, unordered_set<int> &taxa) {
        auto tmp = node;
        
        // if node has only one child it is an unnecessary intermediate node
        if (tmp->children.size() == 1) {
            while (tmp->children.size() == 1 && taxa.find(tmp->id) == taxa.end())
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
            cleanUp(n, taxa);
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
        
        unordered_set<int> dummy;
        
        auto children = root->children;
        for (auto n : children) {
            cleanUp(n, dummy);
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
    
        if (map_.find(t1) == map_.end() || map_.find(t2) == map_.end()) {
            cerr << t1 << " or " << t2 << " are unknown taxids in the map at lca(t1, t2)" << endl;
        }
        
        auto node1 = map_.at(t1);
        auto node2 = map_.at(t2);
        
        int min = std::min(node1->level, node2->level);
        
        while (node1->level > min) {
            node1 = node1->parent;
        }
    
        while (node2->level > min) {
            node2 = node2->parent;
        }
        
        while (node1 != node2) {
            node1 = node1->parent;
            node2 = node2->parent;
        }
        return node1->id;
    }
    
    std::string getOriginal(int custom) override {
        if (to_ncbi_.find(custom) != to_ncbi_.end())
            return to_string(to_ncbi_.at(custom));
        else {
            cerr << "custom: " << custom << " missing at getNCBI(int custom)" << endl;
            return "";
        }
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
    
    bool hasNode(int ncbi) override {
        return (from_ncbi_.find(ncbi) != from_ncbi_.end());
    }
    
    void loadCustomNames(string path) {
        cout << "custom names: " << path << endl;
    
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
    
    void loadCustomNodes(string path) override {
        clear();
        
        cout << "custom nodes: " << path << endl;
        
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
        cout << "nodes: " << map_.size() << endl;
    }
    
    int getNCBI(int custom) {
        if (to_ncbi_.find(custom) != to_ncbi_.end())
            return to_ncbi_.at(custom);
        else {
            cerr << "custom: " << custom << " missing at getNCBI(int custom)" << endl;
            return -1;
        }
    }
    
    int getCustom(int ncbi) {
        if (from_ncbi_.find(ncbi) != from_ncbi_.end())
            return from_ncbi_.at(ncbi);
        else {
            cerr << "NCBI: " << ncbi << " missing at getCustom(int ncbi)" << endl;
            return -1;
        }
    }
    
    int getCustom(std::string ncbi) {
        int tid = stoi(ncbi);
        if (from_ncbi_.find(tid) != from_ncbi_.end())
            return from_ncbi_.at(tid);
        else {
            cerr << "NCBI: " << ncbi << " missing at getCustom(int ncbi)" << endl;
            return -1;
        }
        return -1;
    }
    
    TaxonomyNode* getNode(int custom) override {
        if (map_.find(custom) != map_.end())
            return map_.at(custom);
        else {
            cerr << "id unknown: " << custom << " for getNode(custom)" << endl;
            return nullptr;
        }
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
    
    NCBITaxonomy(){}
    
    void loadNodes(std::string path) {
        int count = 0;
        
        if (!Utils::exists(path)) {
            cerr << "Node file " << path << " does not exist." << endl;
            exit(0);
        }
        
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
    
    void subsetByTaxa(string path) {
        cout << "subsetByTaxa" << endl;
        vector<int> taxa;
        
        ifstream taxa_in;
        string line;
        taxa_in.open(path);
    
        if (!taxa_in.is_open()) {
            cerr << "file " << path << " could not be loaded." << endl;
            exit(0);
        }
        
        unordered_set<int> taxa_set;
        while (getline(taxa_in, line)) {
            std::cout << "line: " << line << std::endl;
            if (map_.find(stoi(line)) != map_.end() && std::find(taxa.begin(), taxa.end(), stoi(line)) == taxa.end()) {
                taxa_set.insert(stoi(line));
            } else {
                cout << "discard " << line << endl;
            }
        }
        taxa_in.close();
    
        
        //auto taxa_set = getCompleteSetOfTaxa(taxa);
        //cout << "ok" << endl;
        auto orig_set = taxa_set;
        
        cout << "final set needs to contain at least " << taxa_set.size() << " ids." << endl;

        int root_id = -1;
        for (auto id : taxa_set) {
            if (root_id == -1) root_id = id;
            else root_id = lca(root_id, id);
        }
        cout << "root: " << root_id << endl;
        
        
        TaxonomyNode* root = map_.at(root_id);
        root->parent = nullptr;
        
        cout << "makesubsetcomplete" << endl;
        makeSubsetComplete(taxa_set, root_id);
        cout << "iscompletecheck" << endl;
        isSubsetComplete(taxa_set, root_id);
    
    
        
        subsetWorker(root, taxa_set);
        cout << "after worker size: " << map_.size() << endl;
    
    
        
        deleteUnlinkedNodes(root);
        cout << "after delete unlinked size: " << map_.size() << endl;
        
        
        int rel = 0;
        int ire = 0;
        countRelevantNodesInSubtree(root, rel, ire);
    
        cout << "relevant: " << rel << endl;
        cout << "irrelevant: " << ire << endl;
    
        auto children = root->children;
        for (auto n : children) {
            cout << "cleanup: " << n->id << "    ->    ";
            cleanUp(n, taxa_set);
            cout << "cleanup. 82993: " << (map_.find(82993) != map_.end()) << endl;
        }
        
        cout << "7. 82993: " << (map_.find(82993) != map_.end()) << endl;
    
        cout << "final check" << endl;
        for (auto id : orig_set) {
            if (taxa_set.find(id) == taxa_set.end()) {
                cout << "invalid." << id << endl;
                exit(0);
            }
        }
    
        cout << "Validity check. .";
        for (auto t : taxa_set) {
            if (map_.find(t) == map_.end()) {
                cout << " . failed for node: " << t << endl;
                exit(1);
            }
        }
        cout << " . passed" << endl;
        
        relabel(root);
        relevel(root);
    
        
        
        custom = true;
    }
    
    void loadNames(std::string path) {
        if (!Utils::exists(path)) {
            cerr << "Names file " << path << " does not exist." << endl;
            exit(0);
        }
        
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
};

class TaxonomyNodeComparer {
public:
    TaxonomyInterface* taxonomy;
    explicit TaxonomyNodeComparer(TaxonomyInterface* taxonomy) : taxonomy(taxonomy) {};
//    ~TaxonomyNodeComparer() {
//        cout << omp_get_thread_num << ": destroy TaxonomyNodeComparer" << endl;
//    }
    bool operator() (const int &a, const int &b) const {
        return taxonomy->getNode(a)->level > taxonomy->getNode(b)->level;
    }
};

