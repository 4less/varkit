//
// Created by fritsche on 01/03/2021.
//

#pragma once

#include <iostream>
#include "robin_map.h"
#include "Taxonomy.h"
#include "Utils.h"

// Use robin map as hashmap (same interface as std::unordered_map)
typedef tsl::robin_map<uint32_t, TaxonomyNode*> TaxonomyMap;
typedef tsl::robin_map<uint32_t, tsl::robin_map<std::string, std::string>> AnnotationMap;
typedef tsl::robin_map<uint32_t, std::string> IntToRankMap;
typedef tsl::robin_map<std::string, uint32_t> RankToIntMap;


class InternalTaxonomy {
public:

    InternalTaxonomy(std::string nodes_path, std::string annotation_path) {
        if (!Utils::exists(nodes_path)) {
            std::cerr << "Taxonomy needs a node file." << "\nGiven file: " << nodes_path << " does not exist." << std::endl;
            exit(0);
        }
        if (!Utils::exists(annotation_path)) {
            std::cerr << "Taxonomy needs an annotation file." << "\nGiven file: " << annotation_path << " does not exist." << std::endl;
            exit(0);
        }

        initRanks();
        readNodes(nodes_path);
        readAnnotation(annotation_path);
    }

    std::string getRank(uint32_t id) {
        if (int_to_ranks_.find(id) != int_to_ranks_.end()) {
            return int_to_ranks_.at(id);
        } else {
            return "unknown";
        }
    }

    uint32_t getRankId(std::string rank) {
        if (ranks_to_int_.find(rank) != ranks_to_int_.end()) {
            return ranks_to_int_.at(rank);
        } else {
            return -1;
        }
    }

    std::string getAnnotation(uint32_t id) {
        if (annotation_.find(id) != annotation_.end()) {
            auto map = annotation_.at(id);
            if (map.find("name") != map.end()) {
                return map["name"];
            }
        }
        return "";
    }

    TaxonomyNode* getNode(uint32_t id) {
        return nodes_.find(id) != nodes_.end() ? nodes_[id] : nullptr;
    }

    static bool isNumber(const std::string& s)
    {
        std::string::const_iterator it = s.begin();
        while (it != s.end() && std::isdigit(*it)) ++it;
        return !s.empty() && it == s.end();
    }


    static void testTaxonomy() {
        std::string nodes = "/usr/users/QIB_fr017/fritsche/Projects/varkit/data/taxonomy/gtdb/nodes2.dmp";
        std::string annotation = "/usr/users/QIB_fr017/fritsche/Projects/varkit/data/taxonomy/gtdb/names2.dmp";

        InternalTaxonomy taxonomy(nodes, annotation);

        std::string input;
        std::string input2;
        std::cin >> input;

        while (strcmp(input.c_str(), "stop") != 0) {
            if (isNumber(input)) {
                std::cout << taxonomy.getAnnotation(stoi(input)) << std::endl;
                std::cout << "second number for lca" << std::endl;
                std::cin >> input2;
                if (isNumber(input2)) {
                    std::cout << taxonomy.getAnnotation(stoi(input2)) << std::endl;

                    int tid1 = stoi(input);
                    int tid2 = stoi(input2);

                    auto lca = taxonomy.lca(tid1, tid2);
                    std::cout << "lca: " << lca << std::endl;
                    std::cout << taxonomy.getAnnotation(lca) << std::endl;
                }

            }
            std::cin >> input;
        }

    }

    int lca(int t1, int t2) {

        if (nodes_.find(t1) == nodes_.end() || nodes_.find(t2) == nodes_.end()) {
            cerr << t1 << " or " << t2 << " are unknown taxids in the map at lca(t1, t2)" << endl;
        }

        auto node1 = nodes_.at(t1);
        auto node2 = nodes_.at(t2);

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

private:
    TaxonomyMap nodes_;
    AnnotationMap annotation_;

    IntToRankMap int_to_ranks_;
    RankToIntMap ranks_to_int_;

    void initRanks() {
        int_to_ranks_[0] = "root";
        int_to_ranks_[1] = "kingdom";
        int_to_ranks_[2] = "phylum";
        int_to_ranks_[3] = "class";
        int_to_ranks_[4] = "order";
        int_to_ranks_[5] = "family";
        int_to_ranks_[6] = "genus";
        int_to_ranks_[7] = "species";
        int_to_ranks_[8] = "strain";
        int_to_ranks_[9] = "genome";

        ranks_to_int_["root"] = 0;
        ranks_to_int_["kingdom"] = 1;
        ranks_to_int_["phylum"] = 2;
        ranks_to_int_["class"] = 3;
        ranks_to_int_["order"] = 4;
        ranks_to_int_["family"] = 5;
        ranks_to_int_["genus"] = 6;
        ranks_to_int_["species"] = 7;
        ranks_to_int_["strain"] = 8;
        ranks_to_int_["genome"] = 9;
    }

    void readNodes(std::string path) {
        auto str = csv::load_file(path.c_str());
        auto parser = csv::make_parser(str);

        bool header = true;

        for (auto&& row : parser ) {
            if (header) {
                header = false;
                continue;
            }

            auto it = row.begin();

            uint32_t taxid = (*it).to_int();
            uint32_t parent_taxid = (*(++it)).to_int();
            uint32_t rank_id = (*(++it)).to_int();
            uint32_t level = (*(++it)).to_int();

            TaxonomyNode* node = getOrNew(taxid);
            node->id = taxid;
            node->rank = getRank(rank_id);
            node->level = level;

            if (taxid != 1) {
                TaxonomyNode* parent = getOrNew(parent_taxid);
                node->parent = parent;
                parent->id = parent_taxid;
                parent->children.push_back(node);
            } else {
                node->parent = nullptr;
            }
        }
    }



    inline TaxonomyNode * getOrNew(uint32_t taxid) {
        if (nodes_.count(taxid) == 1) {
            return nodes_.at(taxid);
        } else {
            auto node = new TaxonomyNode();
            nodes_.insert({taxid, node});
            return node;
        }
    }

    void readAnnotation(std::string path) {
        auto str = csv::load_file(path.c_str());
        auto parser = csv::make_parser(str);

        bool header = true;

        for (auto&& row : parser ) {
            if (header) {
                header = false;
                continue;
            }

            auto it = row.begin();

            uint32_t taxid = (*it).to_int();
            std::string name = (*(++it)).to_string();
            std::string type = (*(++it)).to_string();

            annotation_[taxid][type] = name;
        }
    }


};