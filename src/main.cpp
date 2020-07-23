#include <iostream>
#include <cxxopts.hpp>
#include <assert.h>
#include "BHashMap.h"
#include "VarkitExecutor.h"
#include "OptionsContainer.h"
#include <filesystem>
#include "iotest.h"
#include "Taxonomy.h"

#define TEST 0

using namespace std;

template class BHashMap<44,20,CustomHash>;

// define program modes
std::vector<pair<string,string>> modes;

const static string build_mode = "build";
const static string db_stats_mode = "db-stats";
const static string classify_mode = "classify";
const static string new_shape_mode = "new-shape";

cxxopts::Options build_options("varkit build", "Varkit is a program for classifying metagenomic reads, identify strains by using unknown variants and link and track these unknown strains accross samples.");
cxxopts::Options db_stats_options("varkit db_stats", "Coming soon..");
cxxopts::Options classify_options("varkit classify", "Coming soon..");
cxxopts::Options new_shape_options("varkit new_shape", "Coming soon..");


void initModes() {
    modes.push_back({
        "build",
        "Build your custom database to then use for classifying reads. If you do not want to use one of the shapes provided, you might want to look at mode <new-shape> and <train-shape> before. Type varkit new-shape or varkit train-shape to display the help for these modes."
    });
    
    modes.push_back({
        "db-stats",
        "Get stats on a database that has either been built using the mode <build> or built using the supplied script for the varkit default database."
    });
    
    modes.push_back({
        "classify",
        "Classify unknown reads and extract the variants by using an existing database."
    });
    
    modes.push_back({
        "new-shape",
        "Chose this mode if you want to use a custom shape."
    });
    
    modes.push_back({
        "train-shape",
        "Chose this mode if you have created a custom shape with <new-shape> and now want to train the pattern to variant map."
    });
    
    modes.push_back({
        "help",
        "Display this help. For further help on the different modes, run 'varkit <mode>'."
    });
}

void initBuildOptions() {
    build_options.add_options()
            //("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("c,capacity", "Initial db capacity (to avoid many resizes)", cxxopts::value<int>())
            ("d,db", "Folder to database. Must contain shape information in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("h,help", "Display help.")
            ("fasta", "reference files", cxxopts::value<std::vector<std::string>>());
}

void initClassifyOptions() {
    
    
    classify_options.positional_help("[reads...]\n\n  [reads...] input one to several read files. Accepted formats are  ")
            //.show_positional_help()
            .add_options()
            ("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("d,db", "Specify path to database (folder)",cxxopts::value<std::string>())
            ("n,no_variants", "If specified, varkit does not identify variants.")
            ("h,help", "Display help.")
            ("reads", "", cxxopts::value<std::vector<std::string>>());
}

void initOptions() {
    initBuildOptions();
    initClassifyOptions();
}

static bool stringComparator(pair<string,string> &i, pair<string,string> &j) {
    return i.first.size() < j.first.size();
}

static string padStringTo(string s, int to) {
    assert(s.size() <= to);
    
    int spaces = to - s.size();
    return s + string(spaces, ' ');
}

static string breakString (string s, int line_length, int indent_by = 0, bool indent_first = false) {
    string result = "";
    
    int last_index = 0;
    int cur_line_len = 0;
    
    if (indent_first) {
        result += string(indent_by, ' ');
        cur_line_len = indent_by;
    }
    auto split = Utils::split(s, " ");
    for (auto word : split) {
        int cur_line_len_tmp = cur_line_len + word.length();
        if (cur_line_len_tmp < line_length - 1) {
        
        }
        if (cur_line_len_tmp == line_length || cur_line_len_tmp == line_length -1) {
            result += word + "\n";
        }
        if ((cur_line_len_tmp) > line_length) {
            if (cur_line_len == indent_by) {
                int wraps = word.length() / (line_length-indent_by);
                for (int i = 0; i < wraps; i++) {
                    result += word.substr(i * (line_length-indent_by), line_length-indent_by) + "\n";
                }
                result += word.substr(wraps * (line_length-indent_by), line_length-indent_by) + "\n";
                
            }
        }
    }
}

static void printGeneralHelp(int min_dist, int max_cols) {
    // Formatting
    const auto p = max_element(modes.begin(), modes.end(), stringComparator);
    
    int mode_col_size = p->first.size() + min_dist;
    int desc_col_size = max_cols - mode_col_size;
    
    cout << "Usage: " << endl;
    cout << "varkit <mode> <...> (e.g. 'varkit build' to display help for <build> mode" << endl << endl;
    
    for (auto p : modes) {
        cout << padStringTo(p.first, mode_col_size);
        cout << p.second.substr(0, desc_col_size) << endl;
        
        for (int i = desc_col_size; i < p.second.size(); i += desc_col_size) {
            while (p.second[i] == ' ') i++;
            int len = (i + desc_col_size) > p.second.size() ? p.second.size() - i : desc_col_size;
            cout << string(mode_col_size, ' ') << p.second.substr(i, len) << endl;
        }
        cout << endl;
    }
}

#if TEST == 1
#include "test/test.h"
#include "Benchmark.h"
int main() {
    /*cached_computation(4);
    cached_computation(7);
    cached_computation(7);
    cached_computation(0);*/
    
    //threadTest();
    
    /*Benchmark bm("fasta_test");
    cout << "start fastaTest" << endl;
    bm.start();
    fastaTest(1);
    bm.stop();
    bm.printResults();*/
    
    //iotest::compareFastaIO("/home/joachim/CLionProjects/varkit/data/test/ultra_test2.fa");
    Benchmark bm("taxonomy");
    bm.start();
    //TaxonomyInterface* taxonomy = new NCBITaxonomy("/home/joachim/CLionProjects/varkit/data/taxonomy/nodes.dmp");
    
    NCBITaxonomy* taxonomy = new NCBITaxonomy("/home/joachim/CLionProjects/varkit/data/taxonomy/nodes.dmp");
    taxonomy->subsetByTaxa("/home/joachim/CLionProjects/varkit/data/test/bacteriax.taxid");
    taxonomy->saveCustomNodes("/home/joachim/CLionProjects/varkit/data/taxonomy/custom.dmp");
    taxonomy->loadNames("/home/joachim/CLionProjects/varkit/data/taxonomy/names_wo.dmp");
    taxonomy->saveCustomNames("/home/joachim/CLionProjects/varkit/data/taxonomy/custom_names.dmp");
    
    return 0;
    //NCBITaxonomy* taxonomy = new NCBITaxonomy();
    taxonomy->loadCustomNodes("/home/joachim/CLionProjects/varkit/data/taxonomy/custom.dmp");
    cout << taxonomy->getNCBI(1) << endl;
    
    //taxonomy->loadNames("/home/joachim/CLionProjects/varkit/data/taxonomy/names_wo.dmp");
    //cout << "subset names" << endl;
    //taxonomy->subsetNames();
    cout << "savecustom names" << endl;
    //taxonomy->saveCustomNames("/home/joachim/CLionProjects/varkit/data/taxonomy/custom_names.dmp");
    
    taxonomy->loadCustomNames("/home/joachim/CLionProjects/varkit/data/taxonomy/custom_names.dmp");
    taxonomy->print(taxonomy->getNode(1));
    
    cout << taxonomy->getCustom(1315283) << endl;
    
    return 0;
    //neglect:
    //vector<int> taxids = {85991, 810, 1444190, 813, 562};
    //taxonomy->shrinkToSubtreesOf(taxids);
    //
    
    taxonomy->loadCustomNodes("/home/joachim/CLionProjects/varkit/data/taxonomy/custom.dmp");
    bm.stop();
    bm.printResults();
    int lca = taxonomy->lca(3514,3515);
    cout << lca << endl;
    cout << taxonomy->getNCBI(lca) << endl;
}
#else
int main(int argc, char *argv[]) {
    initModes();
    initOptions();
    
    int min_dist = 3;
    int max_cols = 79;
    
    //printGeneralHelp(min_dist, max_cols);
    
    if (argc < 2) {
        cout << "Please choose runmode." << endl << endl;
        printGeneralHelp(min_dist, max_cols);
    }
    
    if (argc >= 2) {
        int argc_new = argc-1;
        char** argv_copy = new char*[argc];
        
        for (int i = 0; i < argc_new; i++)
            argv_copy[i] = argv[i+1];
        
        if (string(argv[1]) == build_mode) {
            
            build_options.parse_positional({"fasta"});
            
            auto result = build_options.parse(argc_new, argv_copy);
    
            
            if (result.count("help")) {
                cout << "BUILD MODE" << endl;
                cout << build_options.help() << endl;
                exit(0);
            }
            if (!result.count("fasta")) {
                cout << "Please provide fasta files from which to build the db" << endl;
                exit(0);
            }
            
            BuildOptionsContainer build_options_container(
                    1, //result["threads"].as<int>(),
                    result["db"].as<string>(),
                    result["fasta"].as<vector<string>>(),
                    result["capacity"].as<int>()
                    );
            
            VarkitExecutor::runBuild(build_options_container);
            
        } else if (string(argv[1]) == db_stats_mode) {
        
        } else if (string(argv[1]) == classify_mode) {
            cout << "CLASSIFY MODE" << endl;
            classify_options.parse_positional({"reads"});
            auto result = classify_options.parse(argc_new, argv_copy);
    
            if (result.count("help")) {
                cout << classify_options.help() << endl;
                exit(0);
            }
    
            if (!result.count("db")) {
                cout << "You must specify a path to the database." << endl;
            }
            
            MetaDataDB meta_db = loadMetaDataDB(result["db"].as<string>());
            
            ClassifyOptionsContainer classify_options_container(
                    meta_db,
                    result["reads"].as<vector<string>>().at(0),
                    result["threads"].as<int>()
                    );
    
            VarkitExecutor::runClassify(classify_options_container);
    
    
        } else if (string(argv[1]) == new_shape_mode) {
        
        } else {
            cout << "'" << string(argv[1]) << "' is not a valid mode. Please provide a valid mode." << endl << endl;
            printGeneralHelp(min_dist, max_cols);
        }
    }
    
    return 0;
}
#endif

