//
// Created by fritsche on 07/10/2021.
//

#pragma once

#include <TaxonomyNew.h>



class TestTaxonomy {
    Taxonomy::StdTaxonomy taxonomy;

    TestTaxonomy() {}

public:
    static void Test() {
        std::string nodes = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/humgut/pipeline/gtdb/nodes.dmp";
        std::string names = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/humgut/pipeline/gtdb/names.dmp";
        std::string subset = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/humgut/pipeline/genome2iid.tsv";
        std::string output = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/humgut/pipeline/gtdb/internal_taxonomy.dmp";

        Taxonomy::StdTaxonomy taxonomy;
        taxonomy.LoadNodes(nodes);
        taxonomy.LoadNames(names);
        taxonomy.SubsetAndRelabel(subset);
        taxonomy.Export(output);

        std::string nodes2 = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/humgut/pipeline/ncbi/nodes.dmp";
        std::string names2 = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/humgut/pipeline/ncbi/names.dmp";
        std::string output2 = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/humgut/pipeline/ncbi/internal_taxonomy_ncbi.dmp";

        Taxonomy::StdTaxonomy taxonomy_ncbi;
        taxonomy_ncbi.LoadNodes(nodes2);
        taxonomy_ncbi.LoadNames(names2);
        taxonomy_ncbi.SubsetAndRelabel(subset);
        taxonomy_ncbi.Export(output2);

    }

    static void LoadInteralTest() {
        std::string internal = "/usr/users/QIB_fr017/fritsche/Projects/data/taxonomy/uhgg/ncbi/internal_taxonomy.dmp";

        Taxonomy::IntTaxonomy taxonomy(internal);


    }
};