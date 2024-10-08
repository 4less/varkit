//
// Created by fritsche on 05/04/2022.
//

#pragma once

#include <iostream>
#include <cxxopts.hpp>
#include <assert.h>
#include "BHashMap.h"
#include "VarkitExecutor.h"
#include "OptionsContainer.h"
#include <experimental/filesystem>
#include <kmer_iterator_test.h>
#include <ShapeBuilder.h>
#include "iotest.h"
#include "test2.h"
#include "Taxonomy.h"
#include "TaxonomyNew.h"
#include "CustomTaxonomy.h"
//#include "build_scripts.h"
#include "classify/classify_scripts.h"
#include <filesystem>
#include <thread_pool.h>
#include <modes/build.h>
#include <modes/classify.h>
#include "ShapeHistory.h"
#include "FFReader.h"
#include "TestTaxonomy.h"
#include "VCFTest.h"
#include "AbundanceEstimator.h"
#include "modes/profile.h"

#include <taxonomy/HitsClassifier.h>
#include <build/build_scripts.h>

#define TEST 0

using namespace std;

template class BHashMap<44,20,CustomHash>;

// define program modes
std::vector<pair<string,string>> modes;

const static string build_mode = "build";
const static string build_mg_mode = "build-mg";
const static string build_mmg_mode = "build-mmg";
const static string build_pre_mode = "build-pre";
const static string build_mg_buckets_mode = "build-mg-buckets";
const static string taxonomy_mode = "taxonomy";
const static string db_stats_mode = "db-stats";
const static string query_db_mode = "query-db";
const static string classify_mode = "classify";
const static string classify_mg_mode = "classify-mg";
const static string classify_mmg_mode = "classify-mmg";
const static string profile_mode = "profile";
const static string new_shape_mode = "new-shape";
const static string find_shape_mode = "find-shape";
//const static string train_shape_mode = "train-shape";
const static string test_mode = "test-mode";
const static string help_mode = "help";

cxxopts::Options build_options("varkit build", "Varkit is a program for classifying metagenomic reads, identify strains by using unknown variants and link and track these unknown strains accross samples.");
cxxopts::Options build_mg_options("varkit " + build_mg_mode, "Build database for marker genes. ");
cxxopts::Options build_mmg_options("varkit " + build_mmg_mode, "Build database for marker genes. ");
cxxopts::Options build_pre_options("varkit " + build_pre_mode, "Build database for pre_classification. ");
cxxopts::Options build_mg_buckets_options("varkit " + build_mg_buckets_mode, "Build database for marker genes. ");
cxxopts::Options taxonomy_options("varkit taxonomy", "Create a taxonomy subset for set of taxa (NCBI).");
cxxopts::Options query_db_options("varkit query-db", "Coming soon..");
cxxopts::Options db_stats_options("varkit db-stats", "Coming soon..");
cxxopts::Options classify_options("varkit classify", "Coming soon..");
cxxopts::Options profile_options("varkit profile", "Coming soon..");
cxxopts::Options classify_mg_options("varkit classify-mg", "Coming soon..");
cxxopts::Options classify_mmg_options("varkit classify-mg", "Coming soon..");
cxxopts::Options new_shape_options("varkit new-shape", "Coming soon..");
cxxopts::Options find_shape_options("varkit find-shape", "Coming soon..");
//cxxopts::Options train_shape_options("varkit train-shape", "Coming soon..");

void initModes() {
    modes.push_back({
                            build_mode,
                            "Build your custom database to then use for classifying reads. If you do not want to use one of the shapes provided, you might want to look at mode <new-shape> and <train-shape> before. Type varkit new-shape or varkit train-shape to display the help for these modes."
                    });

    modes.push_back({
                            build_mg_mode,
                            "Build your custom marker gene database to then use for classifying reads. If you do not want to use one of the shapes provided, you might want to look at mode <new-shape> and <train-shape> before. Type varkit new-shape or varkit train-shape to display the help for these modes."
                    });

    modes.push_back({
                            build_mmg_mode,
                            "Build your custom marker gene database to then use for classifying reads. If you do not want to use one of the shapes provided, you might want to look at mode <new-shape> and <train-shape> before. Type varkit new-shape or varkit train-shape to display the help for these modes."
                    });

    modes.push_back({
                            build_pre_mode,
                            "Build your custom mg pre database to then use for classifying reads. If you do not want to use one of the shapes provided, you might want to look at mode <new-shape> and <train-shape> before. Type varkit new-shape or varkit train-shape to display the help for these modes."
                    });

    modes.push_back({
                            taxonomy_mode,
                            "Prepare taxonomy for database building."
                    });

    modes.push_back({
                            db_stats_mode,
                            "Get stats on a database that has either been built using the mode <build> or built using the supplied script for the varkit default database."
                    });

    modes.push_back({
                            query_db_mode,
                            "Get stats on a database that has either been built using the mode <build> or built using the supplied script for the varkit default database."
                    });

    modes.push_back({
                            classify_mode,
                            "Classify unknown reads and extract the variants by using an existing database."
                    });

    modes.push_back({
                            classify_mg_mode,
                            "Classify unknown reads and extract the variants by using an existing database."
                    });

    modes.push_back({
                            classify_mmg_mode,
                            "Classify unknown reads and extract the variants by using an existing database."
                    });

    modes.push_back({
                            profile_mode,
                            "Profile classification output."
                    });

    modes.push_back({
                            new_shape_mode,
                            "Chose this mode if you want to use a custom shape."
                    });

    modes.push_back({
                            find_shape_mode,
                            "Chose this mode if you want to use a custom shape."
                    });

//    modes.push_back({
//        train_shape_mode,
//        "Chose this mode if you have created a custom shape with <new-shape> and now want to train the pattern to variant map."
//    });

    modes.push_back({
                            "help",
                            "Display this help. For further help on the different modes, run 'varkit <mode>'."
                    });
}

void initBuildOptions() {
    build_options.add_options()
            //("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("c,capacity", "Initial db capacity (to avoid many resizes)", cxxopts::value<std::string>())
            ("d,db", "Folder to database. Must contain shape information in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("t,threads", "How many threads to use.", cxxopts::value<size_t>())
            ("v,value_bits", "How bits does the value take up.", cxxopts::value<size_t>())
            ("s,shape", "Kmer shape", cxxopts::value<std::string>())
            ("h,help", "Display help.")
            ("fasta", "reference files", cxxopts::value<std::vector<std::string>>());
}

void initBuildMGOptions() {
    build_mg_options.add_options()
            //("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("d,db", "Folder to database. Must contain shape information and taxonomy in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("t,threads", "How many threads to use.", cxxopts::value<size_t>())
            ("v,value", "Value size.", cxxopts::value<int>())
            ("y,taxonomy", "Taxonomy for which to build database.", cxxopts::value<std::string>())
            ("g,max_ram", "Max ram available. This is only for estimating the bucket sizes prior to building the database. If the estimated database size is larger than the available memory, the database construction will fail.", cxxopts::value<size_t>())
            ("l,load_factor", "Initial load factor.", cxxopts::value<double>())
            ("r,target_rank", "If the fasta files have e.g. internal taxids on species level in their headers but you want genus level then provide --target_rank genus.", cxxopts::value<std::string>())
            ("b,bucket", "Bucket table size (bits off key)", cxxopts::value<size_t>())
            ("e,eval", "Evaluate DB after building.")
            ("h,help", "Display help.")
            ("f,force", "Force redo every step.")
            ("fasta", "reference files", cxxopts::value<std::vector<std::string>>());
}

void initBuildMMGOptions() {
    build_mmg_options.add_options()
            //("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("d,db", "Folder to database. Must contain shape information and taxonomy in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("t,threads", "How many threads to use.", cxxopts::value<size_t>())
            ("v,value", "Value size.", cxxopts::value<int>())
            ("l,load_factor", "Initial load factor.", cxxopts::value<double>())
            ("b,bucket", "Bucket table size (bits off key)", cxxopts::value<size_t>())
            ("h,help", "Display help.")
            ("fasta", "reference files", cxxopts::value<std::vector<std::string>>());
}

void initBuildPreOptions() {
    build_pre_options.add_options()
            //("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("d,db", "Folder to database. Must contain shape information and taxonomy in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("t,threads", "How many threads to use.", cxxopts::value<size_t>())
            ("v,value", "Value size.", cxxopts::value<int>())
            ("l,load_factor", "Initial load factor.", cxxopts::value<double>())
            ("b,bucket", "Bucket table size (bits off key)", cxxopts::value<size_t>())
            ("h,help", "Display help.")
            ("fasta", "reference files", cxxopts::value<std::vector<std::string>>());
}

void initBuildMGBucketsOptions() {
    build_mg_buckets_options.add_options()
            ("b,buckets", "Folder to database. Must contain shape information and taxonomy in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("f,change_factor", "Initial load factor.", cxxopts::value<double>())
            ("h,help", "Display help.")
            ("targets", "reference files", cxxopts::value<std::vector<size_t>>());
}

void initNewShapeOptions() {
    new_shape_options.add_options()
            //("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("d,db", "Folder to database. Must contain shape information and taxonomy in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("s,shape_string", "Shape string in Format X_X_XX_X_X (example). X denote positions that contribute to the k-mer, _ positions are ignored.", cxxopts::value<std::string>())
            ("t,threads", "How many threads to use.", cxxopts::value<size_t>())
            ("i,iterations", "How many threads to use.", cxxopts::value<size_t>())
            ("m,mutations", "How many threads to use.", cxxopts::value<size_t>())
            ("r,read_length", "How many threads to use.", cxxopts::value<size_t>())
            ("y,test", "How many threads to use.")
            ("p,pattern_size", "Pattern size in bits (default 22).", cxxopts::value<size_t>())
            ("z,train_stats", "Train pattern and output detailed statistics.")
            ("l, limit", "Train up to <limit> mutations per pattern (default=22). Larger values drastically increase computation time.", cxxopts::value<size_t>())
            ("h,help", "Display help.");
}

void initFindShapeOptions() {
    find_shape_options.add_options()
            //("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("o,output", "Output.", cxxopts::value<std::string>())
            ("t,threads", "How many threads to use.", cxxopts::value<size_t>())
            ("i,iterations", "How many threads to use.", cxxopts::value<size_t>())
            ("m,limit_inc", "How many threads to use.", cxxopts::value<size_t>())
            ("k,kmer", "How many threads to use.", cxxopts::value<size_t>())
            ("s,spaces", "How many spaces to use.", cxxopts::value<size_t>())
            ("r,random", "Probability for random read.", cxxopts::value<double>())
            ("p,pattern_size", "Pattern size in bits (default 16).", cxxopts::value<size_t>())
            ("l, limit", "Train up to <limit> mutations per pattern (default=16). Larger values drastically increase computation time.", cxxopts::value<size_t>())
            ("h,help", "Display help.");
}

void initTaxonomyOptions() {
    taxonomy_options.add_options()
            ("o,output", "output folder for internal_taxonomy.dmp'", cxxopts::value<std::string>())
            ("t,taxa", "Taxa mapping to ids.", cxxopts::value<std::string>())
            ("n,nodes", "Folder to ncbi nodes file", cxxopts::value<std::string>())
            ("m,names", "Folder to ncbi names file.", cxxopts::value<std::string>())
            ("h,help", "Display help.");
}

void initQueryDBOptions() {
    query_db_options.add_options()
            ("d,db", "Specify path to database (folder)",cxxopts::value<std::string>())
//            ("r,rank_composition", "Bool",cxxopts::value<bool>())
            ("h,help", "Display help.");
//            ("reads", "", cxxopts::value<std::vector<std::string>>());
}

void initClassifyOptions() {

    classify_options.positional_help("[reads...]\n\n  [reads...] input one to several read files. Accepted formats are  ")
                    //.show_positional_help()
            .add_options()
                    ("o,output", "output prefix.", cxxopts::value<std::string>())
                    ("d,db", "Specify path to database (folder)",cxxopts::value<std::string>())
                    ("t,threads", "Number of threads to use. (Multithreading based on OpenMP library)", cxxopts::value<int>())
                    ("f,force", "Force rerun classification.")
                    ("h,help", "Display help.")
                    ("reads", "", cxxopts::value<std::vector<std::string>>());
}

void initProfileOptions() {

    profile_options.positional_help("...")
            .add_options()
                    ("f,classification_file", "Path to classification.", cxxopts::value<std::string>())
                    ("s,sample_name", "Name of sample (will be stored in output header)", cxxopts::value<std::string>())
                    ("d,db", "Specify path to database (folder) used for classification",cxxopts::value<std::string>())
                    ("h,help", "Display help.")
                    ("r,reads_threshold", "Path to classification.", cxxopts::value<size_t>())
                    ("a,ani_threshold", "Path to classification.", cxxopts::value<double>())
                    ("e,expected_threshold", "Path to classification.", cxxopts::value<double>())
                    ("v,verbose", "Verbose output file.");
}

void initClassifyMGOptions() {
    classify_mg_options.positional_help("[reads...]\n\n  [reads...] input one to several read files. Accepted formats are  ")
                    //.show_positional_help()
            .add_options()
                    ("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
                    ("d,db", "Specify path to database (folder)",cxxopts::value<std::string>())
                    ("y,taxonomy", "Specify taxonomy (if database has different builds for different taxonomies).",cxxopts::value<std::string>())
                    ("t,threads", "Number of threads to use. (Multithreading based on OpenMP library)", cxxopts::value<int>())
                    ("n,no_variants", "If specified, varkit does not identify variants.")
                    ("h,help", "Display help.")
                    ("f,force", "Force rerun classification.")
                    ("l,list", "List available databases.")
                    ("s,sensitive_snps", "Call SNPs between hits, not only species level hits.")
                    ("reads", "", cxxopts::value<std::vector<std::string>>());
}

void initClassifyMMGOptions() {
    classify_mmg_options.positional_help("[reads...]\n\n  [reads...] input one to several read files. Accepted formats are  ")
                    //.show_positional_help()
            .add_options()
                    ("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
                    ("d,db", "Specify path to database (folder)",cxxopts::value<std::string>())
                    ("t,threads", "Number of threads to use. (Multithreading based on OpenMP library)", cxxopts::value<int>())
                    ("n,no_variants", "If specified, varkit does not identify variants.")
                    ("h,help", "Display help.")
                    ("reads", "", cxxopts::value<std::vector<std::string>>());
}


void initOptions() {
    initBuildOptions();
    initBuildMGOptions();
    initBuildMMGOptions();
    initBuildMGBucketsOptions();
    initClassifyOptions();
    initClassifyMGOptions();
    initClassifyMMGOptions();
    initProfileOptions();
    initBuildPreOptions();
    initTaxonomyOptions();
    initQueryDBOptions();
    initNewShapeOptions();
    initFindShapeOptions();
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
    return "";
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


static const string formatDNAInput(string in) {
    string out = "";

    for (int i = 0; i < in.length(); i++) {
        char c = in.c_str()[i];
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'X')
            out += c;
    }
    return out;
}

static int RunTests() {
    TestTaxonomy::LoadInteralTest();
//    TestTaxonomy::Test();

    // Subset according to list with species int pairs

//    auto str = csv::load_file(path.c_str());
//    auto parser = csv::make_parser( str , '\t');
//
//
//    bool header = false;
//
//    tsl::sparse_map<std::string, Taxonomy::TaxId> species2id;
//    std::unordered_set<int> blocked_ids;
//
//    for (auto&& row : parser ) {
//        if (header) {
//            header = false;
//            continue;
//        }
//
//        auto it = row.begin();
//
//        std::string species = (*(it)).to_string();
//        std::cout << "species : " << species << std::endl;
//        int id = (*(++it)).to_int();
//
//        species2id.insert({ species, id });
//        blocked_ids.insert(id);
//    }
//    exit(9);
}

static int run(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);


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

        // BUILD MODE
        if (string(argv[1]) == build_mode) {

            /* ############################################################################
             * Build database using marker genes and Normal hashmap
             * ############################################################################ */

            build_mg_options.parse_positional({ "fasta" });
            auto result = build_mg_options.parse(argc_new, argv_copy);

            // Print Help
            if (result.count("help")) {
                std::cout << "BUILD MG MODE" << endl;
                std::cout << build_mg_options.help() << endl;
                exit(0);
            }

            // Handle cases
            if (!result.count("fasta")) {
                cout << "Please provide fasta files from which to build the db" << endl;
                exit(0);
            }

            std::string db_dir = result["db"].as<std::string>();

            size_t bucket_bits = result["bucket"].as<size_t>();

            bool force = result.count("force");
            uint32_t threads = result.count("threads") ? result["threads"].as<size_t>() : 1;
            double load_factor = result.count("load_factor") ? result["load_factor"].as<double>() : 0.75f;
            auto references = result["fasta"].as<vector<std::string>>();

            std::string taxonomy = result.count("taxonomy") ? result["taxonomy"].as<std::string>() : "";

            VarkitOptionsContainer oc;
            oc.references_ = references;
            oc.db_dir_ = db_dir;
            oc.threads_ = threads;
            oc.load_factor_ = load_factor;
            oc.force_rebuild_ = force;
            oc.offset_key_bits_ = bucket_bits;

            std::cout << "taxonomy: " << taxonomy << std::endl;
            oc.Init(taxonomy);
            oc.LoadInternalTaxonomy();

            BuildMode::Run(oc);


            // TAXONOMY MODE
        } else if (string(argv[1]) == build_pre_mode) {
            build_pre_options.parse_positional({"fasta"});
            auto result = build_pre_options.parse(argc_new, argv_copy);


            if (result.count("help")) {
                std::cout << "BUILD MG MODE" << endl;
                std::cout << build_mg_options.help() << endl;
                exit(0);
            }
            if (!result.count("fasta")) {
                cout << "Please provide fasta files from which to build the db" << endl;
                exit(0);
            }

            bool eval = result.count("eval");
            bool update = result.count("update");


            std::string db_dir = result["db"].as<std::string>();

            if (!Utils::exists(db_dir + "/shape.txt")) {
                exit(10);
            }
            size_t kmer_size = ShapeUtils::GetK(ShapeUtils::LoadShape(db_dir + "/shape.txt"));

            size_t key_bits_total = kmer_size * 2;
            size_t bucket_bits = result["bucket"].as<size_t>();
            size_t key_bits_internal = key_bits_total - bucket_bits;
            size_t value_bits = 64 - key_bits_internal;

            uint32_t threads = result.count("threads") ? result["threads"].as<size_t>() : 1;

            double load_factor = result.count("load_factor") ? result["load_factor"].as<double>() : 0.75f;
            auto references = result["fasta"].as<vector<std::string>>();


            Build::BuildOptionsMG options(
                    key_bits_internal,
                    value_bits,
                    key_bits_total,
                    kmer_size,
                    threads,
                    load_factor,
                    0,
                    bucket_bits,
                    db_dir,
                    true,
                    "",
                    false,
                    references);
//            Build::BuildPreMap(options);

        } else if (string(argv[1]) == build_mg_mode) {
            /* ############################################################################
             * Build database using marker genes and Normal hashmap
             * ############################################################################ */

            build_mg_options.parse_positional({ "fasta" });
            auto result = build_mg_options.parse(argc_new, argv_copy);

            // Print Help
            if (result.count("help")) {
                std::cout << "BUILD MG MODE" << endl;
                std::cout << build_mg_options.help() << endl;
                exit(0);
            }

            // Handle cases
            if (!result.count("fasta")) {
                cout << "Please provide fasta files from which to build the db" << endl;
                exit(0);
            }

            std::string db_dir = result["db"].as<std::string>();

            // Shape file must exist
            if (!Utils::exists(db_dir + "/shape.txt")) {
                std::cerr << (db_dir + "/shape.txt") << " does not exist. Shape file must exist." << std::endl;
                exit(10);
            }


            size_t kmer_size = ShapeUtils::GetK(ShapeUtils::LoadShape(db_dir + "/shape.txt"));

            size_t bucket_bits = result["bucket"].as<size_t>();

            bool force = result.count("help");
            uint32_t threads = result.count("threads") ? result["threads"].as<size_t>() : 1;

            double load_factor = result.count("load_factor") ? result["load_factor"].as<double>() : 0.75f;
            double bloom_limit_gb = result.count("max_ram") ? result["max_ram"].as<size_t>() : 100'000'000;

            auto references = result["fasta"].as<vector<std::string>>();

            size_t key_bits_total = kmer_size * 2;
            size_t key_bits_internal = key_bits_total - bucket_bits;
            size_t value_bits = 64 - key_bits_internal;

            Build::BuildOptionsMG options(
                    key_bits_internal,
                    value_bits,
                    key_bits_total,
                    kmer_size,
                    threads,
                    load_factor,
                    bloom_limit_gb,
                    bucket_bits,
                    db_dir,
                    true,
                    "",
                    false,
                    references);

            VarkitOptionsContainer oc;
            oc.references_ = references;
            oc.db_dir_ = db_dir;
            oc.k_ = kmer_size;
            oc.threads_ = threads;
            oc.internal_key_bits_ = key_bits_internal;
            oc.key_bits_ = key_bits_total;
            oc.load_factor_ = load_factor;
            oc.force_rebuild_ = force;

            oc.Init();

            // Verify variables...
            options.Check();

            Build::BuildMarkerGenes(options);


        } else if (string(argv[1]) == build_mmg_mode) {
            /* ############################################################################
             * Build database using marker genes and multimap
             * ############################################################################ */

            build_mg_options.parse_positional({ "fasta" });
            auto result = build_mg_options.parse(argc_new, argv_copy);

            // Print Help
            if (result.count("help")) {
                std::cout << "BUILD MG MODE" << endl;
                std::cout << build_mg_options.help() << endl;
                exit(0);
            }

            // Handle cases
            if (!result.count("fasta")) {
                cout << "Please provide fasta files from which to build the db" << endl;
                exit(0);
            }

            std::string db_dir = result["db"].as<std::string>();

            // Shape file must exist
            if (!Utils::exists(db_dir + "/shape.txt")) {
                std::cerr << (db_dir + "/shape.txt") << " does not exist. Shape file must exist." << std::endl;
                exit(10);
            }


            size_t kmer_size = ShapeUtils::GetK(ShapeUtils::LoadShape(db_dir + "/shape.txt"));
            size_t key_bits_total = kmer_size * 2;
            size_t bucket_bits = result["bucket"].as<size_t>();
            size_t key_bits_internal = key_bits_total - bucket_bits;
            size_t value_bits = 62 - key_bits_internal; // MULTIMAP NEEDS 2 CONTROL BITS PER ENTRY

            uint32_t threads = result.count("threads") ? result["threads"].as<size_t>() : 1;

            double load_factor = result.count("load_factor") ? result["load_factor"].as<double>() : 0.75f;
            double bloom_limit_gb = result.count("max_ram") ? result["max_ram"].as<size_t>() : 100'000'000;

            auto references = result["fasta"].as<vector<std::string>>();

            Build::BuildOptionsMMG options(
                    key_bits_internal,
                    value_bits,
                    key_bits_total,
                    kmer_size,
                    threads,
                    load_factor,
                    bloom_limit_gb,
                    bucket_bits,
                    db_dir,
                    references);

            // Verify variables...
            options.Check();

            Build::BuildIndex(options);


        } else if (string(argv[1]) == build_mg_buckets_mode) {

            build_mg_buckets_options.parse_positional({"targets"});
            auto result = build_mg_buckets_options.parse(argc_new, argv_copy);


            if (result.count("help")) {
                std::cout << "BUILD MG MODE" << endl;
                std::cout << build_mg_buckets_options.help() << endl;
                exit(0);
            }

            std::string buckets_dir = result["buckets"].as<std::string>();
            double factor = result.count("f") ? result["f"].as<double>() : 2.0f;


            auto targets = result["targets"].as<vector<std::size_t>>();
//        Build::ChangeBucketSizes(buckets_dir, targets, factor);

//        auto buckets = IndexedMap::LoadBucketSizes(buckets_dir);
//
//        for (int i = 0; i < buckets.second; i++) {
//            std::cout << i << ": " << buckets.first[i] << std::endl;
//        }
//        delete[] buckets.first;

            // TAXONOMY MODE
        } else if (string(argv[1]) == taxonomy_mode) {
            auto result = taxonomy_options.parse(argc_new, argv_copy);

            if (result.count("help")) {
                cout << "TAXONOMY MODE" << endl;
                cout << taxonomy_options.help() << endl;
                exit(0);
            }

            string nodes = result["nodes"].as<string>();
            string names = result["names"].as<string>();
            string taxa = result["taxa"].as<string>();
            string output = result["output"].as<string>() + "/internal_taxonomy.dmp";

            Taxonomy::StdTaxonomy taxonomy;
            cout << "load nodes ..." << endl;
            taxonomy.LoadNodes(nodes);
            cout << "load names ..." << endl;
            taxonomy.LoadNames(names);
            cout << "create subset ..." << endl;
            taxonomy.SubsetAndRelabel(taxa);
            cout << "export ..." << endl;
            taxonomy.Export(output);

        } else if (string(argv[1]) == db_stats_mode) {

        } else if (string(argv[1]) == classify_mode) {
            cout << "CLASSIFY MODE" << endl;

            /* ############################################################################
             * Classify with marker gene and lca based database
             * ############################################################################ */

            classify_mg_options.parse_positional({ "reads" });
            auto result = classify_mg_options.parse(argc_new, argv_copy);

            // Print Help
            if (result.count("help")) {
                std::cout << classify_mg_options.help() << endl;
                exit(0);
            }

            if (result.count("l")) {
                std::string db_dir = result["db"].as<std::string>();
                std::cout << "List available databases:" << std::endl;
                VarkitOptionsContainer oc;
                oc.db_dir_ = db_dir;
                oc.LoadIndices();
                oc.ListIndices();
                exit(0);
            }


            // Handle cases
            if (!result.count("reads")) {
                cout << "Please provide fasta files from which to build the db" << endl;
                exit(0);
            }

            std::string db_dir = result["db"].as<std::string>();
            std::string output_prefix = result["output"].as<std::string>();

            std::string taxonomy = result.count("taxonomy") ? result["taxonomy"].as<std::string>() : "";

            bool force = result.count("force");

            int threads = result.count("threads") ? result["threads"].as<int>() : 1;
            auto reads = result["reads"].as<vector<std::string>>();

            varkit::Benchmark varkit_bm("Varkit");

            VarkitOptionsContainer oc;
            oc.LoadReadstrings(reads);
            oc.reads_ = reads;
            oc.db_dir_ = db_dir;
            oc.output_prefix_ = output_prefix;
            oc.output_dir_ = output_prefix;
            oc.threads_ = threads;
            oc.force_rebuild_ = force;
            oc.Init(taxonomy);

            std::vector<Sample> samples_classify;
            std::vector<Sample> samples_profile;

            std::cout << std::string(60, '-') << std::endl;
            for (auto& sample : oc.GetSamples()) {
                std::cout << sample.GetSampleString() << std::endl;
                std::cout << "Classification output: " << Utils::exists(oc.ClassificationOutputFile(sample)) << std::endl;
                std::cout << "Profile output:        " << Utils::exists(oc.ProfileOutputFile(sample)) << std::endl;
                std::cout << std::string(60, '-') << std::endl;
            }

            if (force) {
                samples_classify = oc.GetSamples();
                samples_profile = oc.GetSamples();
            } else {
                for (int i = 0; i < oc.GetSamples().size(); i++) {
                    auto &sample = oc.GetSamples()[i];

                    if (!Utils::exists(oc.ClassificationOutputFile(sample))) {
                        samples_classify.emplace_back(sample);
                    }
                    if (!Utils::exists(oc.ProfileOutputFile(sample))) {
                        samples_profile.emplace_back(sample);
                    }
                }
            }

            if (!samples_classify.empty()) {
                varkit::Benchmark load_bm("Loading resources");
                oc.LoadResources();
                load_bm.stop();
                load_bm.printResults();
                oc.IsValidClassifyOptions();
            } else {
                std::cout << "All samples have already been classified" << std::endl;
            }

            for (auto& sample : samples_classify) {
                ClassifyMode::Run(oc,sample);
            }

            // PROFILE
            if (!samples_profile.empty()) {
                if (oc.taxonomy == nullptr) {
                    oc.LoadInternalTaxonomy();
                }

                auto read_filter = ReadFilter(0.9, oc.shape_str.length(), oc.k_);

                size_t r = 100;
                double a = 0.9;
                double e = 0.1;

                auto taxon_filter = TaxonFilter(r, a, e);

                AbundanceEstimator ae(oc.k_, oc.ShapeLength(), read_filter, taxon_filter, *oc.InternalTaxonomy());
                ae.LoadMarkerGenome( oc.MarkerGeneLengthsFile().c_str() );

                for (auto& sample : samples_profile) {
                    auto classification_output = oc.ClassificationOutputFile(sample);
                    ProfileMode::Run(ae, oc, sample);
                }
            } else {
                std::cout << "All samples have already been profiled" << std::endl;
            }


        } else if (string(argv[1]) == classify_mg_mode) {
            cout << "CLASSIFY MG MODE" << endl;
            classify_mg_options.parse_positional({"reads"});
            auto result = classify_mg_options.parse(argc_new, argv_copy);

            if (result.count("help")) {
                cout << classify_mg_options.help() << endl;
                exit(0);
            }

            if (!result.count("db")) {
                cout << "You must specify a path to the database." << endl;
                exit(9);
            }

            int threads = result["threads"].as<int>();
            std::string db_path = result["db"].as<string>();
            std::string output_dir = result["output"].as<string>();
            auto reads = result["reads"].as<vector<string>>();

            bool output_unclassified = false;
            bool sensitive_snps = result.count("sensitive_snps");

            Classify::ClassifyOptions options(
                    threads,
                    db_path,
                    output_dir,
                    output_unclassified,
                    sensitive_snps,
                    reads
            );

            // Verify variables...
            // options.Check();

            Classify::ClassifyMarkerGenes(options);

        } else if (string(argv[1]) == profile_mode) {
            cout << "PROFILE MODE" << endl;
            auto result = profile_options.parse(argc_new, argv_copy);

            if (result.count("help")) {
                std::cout << profile_options.help() << std::endl;
                exit(0);
            }

            if (!result.count("classification_file")) {
                std::cout << "You must specify a path to a classification output." << std::endl;
            }
            if (!result.count("db")) {
                std::cout << "You must specify a path to the database." << std::endl;
            }


            bool verbose = result.count("verbose");

            std::string db_path = result["db"].as<string>();
            std::string sample_name = result["sample_name"].as<string>();
            std::string file_path = result["classification_file"].as<string>();

            size_t r = result.count("r") ? result["r"].as<size_t>() : 100;
            double a = result.count("a") ? result["a"].as<double>() : 0.9;
            double e = result.count("e") ? result["e"].as<double>() : 0.1;

            std::string shape_str = ShapeUtils::LoadShape(db_path + "/shape.txt");
            size_t k = ShapeUtils::GetK(shape_str);

            Taxonomy::IntTaxonomy taxonomy(db_path + "/ncbi/internal_taxonomy.dmp");

            auto read_filter = ReadFilter(0.9, shape_str.length(), k);
            auto taxon_filter = TaxonFilter(r, e, a);

            AbundanceEstimator ae(k, shape_str.length(), read_filter, taxon_filter, taxonomy);
            ae.LoadMarkerGenome( (db_path + "/mg_lengths.tsv").c_str() );
//            ae.LoadAbundanceCorrection( (db_path + "/abundance_correction.tsv").c_str() );

            std::cout << "Paths:" << std::endl;
            std::cout << db_path + "/shape.txt" << std::endl;
            std::cout << shape_str << std::endl;
            std::cout << shape_str << std::endl;
            std::cout << db_path + "/ncbi/internal_taxonomy.dmp" << std::endl;
            std::cout << file_path << std::endl;
            std::cout << db_path + "/mg_lengths.tsv" << std::endl;
            std::cout << db_path + "/abundance_correction.tsv" << std::endl;

            Benchmark profile_bm("Profiling");
            profile_bm.start();
            std::cout << "ReadClassificationFile" << std::endl;
            ae.ReadClassificationFile(file_path.c_str());
            std::cout << "Profile" << std::endl;
            ae.SetVerbose(verbose);
            ae.SetSampleName(sample_name);


            ae.Profile(file_path + ".profile");
            std::cout << std::flush;
            profile_bm.stop();
            profile_bm.printResults();

            std::cout << "r: " << r << " " << result.count("r") << std::endl;
            std::cout << "a: " << a << " " << result.count("a") << std::endl;
            std::cout << "e: " << e << " " << result.count("e") << std::endl;

//            std::cout << "pass:    " << ae.GetReadFilter().GetPassCounter() << std::endl;
//            std::cout << "decline: " << ae.GetReadFilter().GetDeclineCounter() << std::endl;
//            std::cout << "pass:    " << read_filter.GetPassCounter() << std::endl;
//            std::cout << "decline: " << read_filter.GetDeclineCounter() << std::endl;

        } else if (string(argv[1]) == classify_mmg_mode) {
            cout << "CLASSIFY MG MODE" << endl;
            classify_mmg_options.parse_positional({"reads"});
            auto result = classify_mmg_options.parse(argc_new, argv_copy);

            if (result.count("help")) {
                cout << classify_mmg_options.help() << endl;
                exit(0);
            }

            if (!result.count("db")) {
                cout << "You must specify a path to the database." << endl;
            }

            int threads = result["threads"].as<int>();
            std::string db_path = result["db"].as<string>();
            std::string output_dir = result["output"].as<string>();
            auto reads = result["reads"].as<vector<string>>();

            bool output_unclassified = false;

            Classify::ClassifyOptionsMM options(
                    threads,
                    db_path,
                    output_dir,
                    output_unclassified,
                    reads
            );

            Classify::ClassifyMarkerGenesMM(options);


        } else if (string(argv[1]) == new_shape_mode) {
            Benchmark bm("total");
            bm.start();
            // HERE
            auto result = new_shape_options.parse(argc_new, argv_copy);


            if (result.count("help")) {
                std::cout << "NEW SHAPE MODE" << endl;
                std::cout << new_shape_options.help() << endl;
                exit(0);
            }

            std::string db = result["db"].as<std::string>();
            std::string shape_db = db + "/shape.db";

            size_t iterations = result.count("iterations") ? result["iterations"].as<size_t>() : 10000;
            size_t read_length = result.count("read_length") ? result["read_length"].as<size_t>() : 100;
            size_t mutations = result.count("mutations") ? result["mutations"].as<size_t>() : 10;
            int threads = result.count("threads") ? result["threads"].as<size_t>() : 1;


            if (!result.count("test")) {

                std::string shape_str = result["shape_string"].as<std::string>();
                size_t pattern_size = result["p"].as<size_t>();
                size_t limit = result["limit"].as<size_t>();

//            std::cout << result["threads"].as<int>() << std::endl;

                std::cout << shape_str << std::endl;
                std::cout << "psize: " << pattern_size << std::endl;
                std::cout << "limit: " << limit << std::endl;
                std::cout << "threads: " << threads << std::endl;
                std::cout << "db: " << shape_db << std::endl;

                // Train but with stats.
                if (result.count("train_stats")) {
                    ShapeBuilder::DeepTrain(shape_db, db + "/stats.csv", db + "/distr.csv", shape_str, threads, pattern_size, limit);
                    return 0;
                }


                ShapeBuilder* builder;


                if (Utils::exists(shape_db)) {
                    std::cout << "Shape exists.. " << std::endl;
                    builder = new ShapeBuilder();
                    builder->Load(shape_db);


                    std::cout << "pattern_size " << builder->pattern_size << std::endl;
                    std::cout << "builder_limit " << builder->limit << std::endl;
                    std::cout << "builder_limit " << builder->lower_limit << std::endl;

                    if (builder->pattern_size != pattern_size) {
                        std::cerr << "pattern already exists but has been trained with a pattern size of " << builder->pattern_size << " instead of the requested " << pattern_size << std::endl;
                        exit(8);
                    }
                    if (builder->limit >= limit) {
                        std::cerr << "Pattern exists and has been trained to limit " << builder->limit << std::endl;
                        exit(8);
                    }
                    builder->limit = limit;
                } else {
                    builder = new ShapeBuilder(shape_str);
                    builder->pattern_size = pattern_size;
                    builder->limit = limit;
                }

                std::cout << "lower_limit " << builder->lower_limit << std::endl;
                std::cout << "upper_limit " << builder->limit << std::endl;


                ProgressSingleton *progress = ProgressSingleton::GetInstance();
                builder->progress = progress;

                builder->Train(threads);
//                builder->TrainByPatternTP(threads * 3);


                size_t max_pot_mut = 0;
                for (auto i = 0; i < 100; i++) {
                    if (builder->pot_mut_counter[i] != 0)
                        max_pot_mut = i;
                }
                for (auto i = 0; i <= max_pot_mut; i++) {
                    std::cout << i << ": " << builder->pot_mut_counter[i] << std::endl;
                }
                auto median = ShapeBuilder::median(builder->pot_mut_counter, max_pot_mut);
                std::cout << "median: " << median << std::endl;


                std::cout << "save.." << std::endl;
                builder->Save(shape_db);

                bm.stop();
                bm.printResults();

                size_t muts[10] = { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
                double sens[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

                std::cout << "Test..." << std::endl;
                builder->TestTP(read_length, iterations, muts, sens, 10, threads * 3);
                for (auto i = 0; i < 10; i++) {
                    std::cout << 90 + i << " ANI: " << sens[i] << std::endl;
                }
                double mean = ShapeBuilder::Sum(sens, 10) / 10;
                std::cout << "mean: " << mean << std::endl;

                delete builder;
//                builder.Test(shape_db, read_length, iterations, mutations);
            } else {
                std::cout << "Test pattern." << std::endl;
                std::cout << "read_length: " << read_length << std::endl;
                std::cout << "iterations: " << iterations << std::endl;
                std::cout << "mutations: " << mutations << std::endl;
                ShapeBuilder builder;
                builder.Load(shape_db);

//                size_t rl = 100;
//                size_t probe_length = 10000;
//                bool *read = new bool[rl];
//                size_t probes_sum = 0;
//
//                for (int i = 0; i < probe_length; i++) {
//                    builder.GenerateReadRandom(read, rl, 0.05, 0.8);
//
//                    std::cout << builder.ReadToString(read, rl) << std::endl;
//
//                    size_t mut_count = 0;
//                    for (int j = 0; j < rl; j++) {
//                        mut_count += read[j];
//                    }
//                    std::cout << "Mutations: " << mut_count << std::endl;
//                    probes_sum += mut_count;
//                }
//                std::cout << "mean: " << (double) probes_sum / probe_length << std::endl;
//
//
//                delete[] read;
//
//                exit(9);

                size_t muts[10] = { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
                double sens[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

                builder.TestTP(read_length, 10000, muts, sens, 10, threads * 3);

                for (auto i = 0; i < 10; i++) {
                    std::cout << 90 + i << " ANI: " << sens[i] << std::endl;
                }
                double mean = ShapeBuilder::Sum(sens, 10) / 10;
                std::cout << "mean: " << mean << std::endl;

//                builder.Test(shape_db, read_length, iterations, mutations);
            }


        }  else if (string(argv[1]) == find_shape_mode) {
            // HERE
            auto result = find_shape_options.parse(argc_new, argv_copy);

            if (result.count("help")) {
                std::cout << "FIND SHAPE MODE" << endl;
                std::cout << find_shape_options.help() << endl;
                exit(0);
            }

            std::string output_folder = result["o"].as<std::string>();

            size_t iterations = result.count("iterations") ? result["iterations"].as<size_t>() : 1000;

            int threads = result.count("threads") ? result["threads"].as<size_t>() : 1;
            size_t pattern_size = result.count("p") ? result["p"].as<size_t>() : 20;
            size_t limit = result.count("l") ? result["l"].as<size_t>() : 16;
            size_t limit_increment = result.count("m") ? result["m"].as<size_t>() : 10;

            double random = result["r"].as<double>();

            size_t k = result["k"].as<size_t>();
            size_t s = result["s"].as<size_t>();

            FinderOptions options;

//            options.threads = threads;
            options.threads = threads;
            options.training_max = limit;
            options.iterations = iterations;
            options.space = s;
            options.take = k;
            options.pattern_size = pattern_size;
            options.previous_findings = output_folder;


            ShapeFinder finder { options };
            finder.Find();


//            ShapeBuilder::ShapeSpaceExplorer(k, s, output_folder, random, threads, pattern_size, limit, limit_increment, iterations);


        } else if (string(argv[1]) == query_db_mode) {
//            query_db_options.parse_positional({"reads"});
            auto result = query_db_options.parse(argc_new, argv_copy);

            string db = result["db"].as<string>();
//            auto references = result["reads"].as<vector<std::string>>();
//            bool rank_composition = result.count("r");

//            Classify::QueryDB(db, references, rank_composition);
//            Build::RankLevelComposition(db);

            Build::DBCompositionMM(db);


        } else if (string(argv[1]) == test_mode) {

        } else {
            cout << "'" << string(argv[1]) << "' is not a valid mode. Please provide a valid mode." << endl << endl;
            printGeneralHelp(min_dist, max_cols);
        }

        delete[] argv_copy;
    }


    return 0;
}
//#endif

