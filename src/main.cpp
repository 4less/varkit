#include <iostream>
#include <cxxopts.hpp>
#include <assert.h>
#include "BHashMap.h"
#include "VarkitExecutor.h"
#include "OptionsContainer.h"
#include <filesystem>
#include "iotest.h"
#include "test2.h"
#include "Taxonomy.h"

#define TEST 0

using namespace std;

template class BHashMap<44,20,CustomHash>;

// define program modes
std::vector<pair<string,string>> modes;

const static string build_mode = "build";
const static string taxonomy_mode = "taxonomy";
const static string db_stats_mode = "db-stats";
const static string query_db_mode = "query-db";
const static string classify_mode = "classify";
const static string new_shape_mode = "new-shape";
const static string train_shape_mode = "train-shape";
const static string test_mode = "test-mode";
const static string help_mode = "help";

cxxopts::Options build_options("varkit build", "Varkit is a program for classifying metagenomic reads, identify strains by using unknown variants and link and track these unknown strains accross samples.");
cxxopts::Options taxonomy_options("varkit taxonomy", "Create a taxonomy subset for set of taxa (NCBI).");
cxxopts::Options query_db_options("varkit query-db", "Coming soon..");
cxxopts::Options db_stats_options("varkit db-stats", "Coming soon..");
cxxopts::Options classify_options("varkit classify", "Coming soon..");
cxxopts::Options new_shape_options("varkit new-shape", "Coming soon..");
cxxopts::Options train_shape_options("varkit train-shape", "Coming soon..");

void initModes() {
    modes.push_back({
        build_mode,
        "Build your custom database to then use for classifying reads. If you do not want to use one of the shapes provided, you might want to look at mode <new-shape> and <train-shape> before. Type varkit new-shape or varkit train-shape to display the help for these modes."
    });
    
    modes.push_back({
        taxonomy_mode,
        "Build your custom database to then use for classifying reads. If you do not want to use one of the shapes provided, you might want to look at mode <new-shape> and <train-shape> before. Type varkit new-shape or varkit train-shape to display the help for these modes."
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
        new_shape_mode,
        "Chose this mode if you want to use a custom shape."
    });
    
    modes.push_back({
        train_shape_mode,
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
            ("c,capacity", "Initial db capacity (to avoid many resizes)", cxxopts::value<std::string>())
            ("t,taxonomy", "Set the taxonomy. Options are NCBI and GTDB. The taxonomy must be available in the database, a.k.a it must have been build with the taxonomy mode and the specific options.", cxxopts::value<std::string>())
            ("d,db", "Folder to database. Must contain shape information in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>())
            ("s,threads", "How many threads to use.", cxxopts::value<int>())
            ("v,validate", "Validate DB afterwards by extracting the k-mers again and checking if they are in the build index.")
            ("h,help", "Display help.")
            ("fasta", "reference files", cxxopts::value<std::vector<std::string>>());
}

void initTaxonomyOptions() {
    taxonomy_options.add_options()
            ("o,output", "output folder for subsetted files 'nodes.dmp' and 'names.dmp'", cxxopts::value<std::string>())
            ("t,taxa", "Specify taxa to build the tree.", cxxopts::value<std::string>())
            ("n,nodes", "Folder to ncbi nodes file", cxxopts::value<std::string>())
            ("m,names", "Folder to ncbi names file.", cxxopts::value<std::string>())
            ("h,help", "Display help.");
}

void initQueryDBOptions() {
    query_db_options.add_options()
            ("d,db", "Folder to database. Must contain shape information in folder shape (can be created with mode <new_shape> or downloaded)", cxxopts::value<std::string>());
}

void initClassifyOptions() {
    
    
    classify_options.positional_help("[reads...]\n\n  [reads...] input one to several read files. Accepted formats are  ")
            //.show_positional_help()
            .add_options()
            ("o,output", "output folder (must be empty). If it does not exist, it will be created.", cxxopts::value<std::string>())
            ("d,db", "Specify path to database (folder)",cxxopts::value<std::string>())
            ("t,threads", "Number of threads to use. (Multithreading based on OpenMP library)", cxxopts::value<int>())
            ("n,no_variants", "If specified, varkit does not identify variants.")
            ("l,vmer_length", "Length of variant mers with variant position in the middle.", cxxopts::value<size_t>())
            ("h,help", "Display help.")
            ("reads", "", cxxopts::value<std::vector<std::string>>());
}

void initOptions() {
    initBuildOptions();
    initClassifyOptions();
    initTaxonomyOptions();
    initQueryDBOptions();
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

#if TEST == 1
int main(int argc, char *argv[]) {
    cout << "test start" << endl;
    
    test_kmer_buffer();
    
    cout << "TEst end" << endl;
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
        
        // BUILD MODE
        if (string(argv[1]) == build_mode) {
            
            build_options.parse_positional({"fasta"});
            // DEFAULT TAXONOMY CHANGE THIS MECHANISM IN THE FUTURE
            std::string taxonomy = "NCBI";
            
            auto result = build_options.parse(argc_new, argv_copy);
    
            int threads = result.count("threads") ? result["threads"].as<int>() : 1;
            
            if (result.count("help")) {
                cout << "BUILD MODE" << endl;
                cout << build_options.help() << endl;
                exit(0);
            }
            if (!result.count("fasta")) {
                cout << "Please provide fasta files from which to build the db" << endl;
                exit(0);
            }
            if (result.count("taxonomy")) {
                taxonomy = result["taxonomy"].as<std::string>();
            }
            
            BuildOptionsContainer build_options_container(
                    threads,
                    result["db"].as<string>(),
                    result["fasta"].as<vector<std::string>>(),
                    stoll(result["capacity"].as<std::string>()),
                    taxonomy,
                    result.count("validate")
                    );
            
            VarkitExecutor::runBuild(build_options_container);
            
            cout << "done" << endl;
            
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
            string output = result["output"].as<string>();
            
            NCBITaxonomy taxonomy;
            cout << "load nodes ..." << endl;
            taxonomy.loadNodes(nodes);
            
            cout << "create subset ..." << endl;
            taxonomy.subsetByTaxa(taxa);

            
            cout << "save custom nodes ..." << endl;
            taxonomy.saveCustomNodes(output + "/nodes.dmp");
            
            cout << "load names ..." << endl;
            taxonomy.loadNames(names);
            
            cout << "subset names ..." << endl;
            taxonomy.subsetNames();
        
            cout << "save custom names ..." << endl;
            taxonomy.saveCustomNames(output + "/names.dmp");
            
        
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
    
            size_t vmer_length = 15;
            if (!result.count("vmer_length")) {
                vmer_length = result["vmer_length"].as<size_t>();
            }
            
            MetaDataDB meta_db = loadMetaDataDB(result["db"].as<string>());
            
            ClassifyOptionsContainer classify_options_container(
                    meta_db,
                    result["reads"].as<vector<string>>(),
                    0,//result["threads"].as<int>(),
                    result["output"].as<string>(),
                    vmer_length
                    );
    
            
            VarkitExecutor::runClassifyOMP(classify_options_container);
    
    
        } else if (string(argv[1]) == new_shape_mode) {
        
        } else if (string(argv[1]) == query_db_mode) {
            cout << "QUERY DB MODE" << endl;
            auto result = query_db_options.parse(argc_new, argv_copy);
            
            string db = result["db"].as<string>();
            
            cout << "metapath: " << db << endl;
            MetaDataDB meta_db = loadMetaDataDB(db);
    
            // Shape information
            cout << "shape_path: " << db + MetaDataDB::SHAPE_FILE << endl;
            string shape_str = MetaDataDB::loadShape(db + MetaDataDB::SHAPE_FILE);
            const size_t shape_length = shape_str.length();
            bool * shape = MetaDataDB::getShape(shape_str);
            cout << "shape: " << shape_str << endl;
            cout << "shape_length: " << shape_str.length() << endl;
            
    
            BHashMap<44, 20, CustomHash> map (meta_db.capacity, meta_db.load_factor, meta_db.growth_factor);
            map.load(meta_db.path + "/index.db", meta_db.path + "/index.meta");
            cout << "done loading." << endl;
            
            int k = 22;
            const int key_bytes = (2*k + 8 - 1) /  8;
            

            SpacedKmerIterator iterator(k, shape, shape_length);
    
            static uint8_t *key = new uint8_t[key_bytes];
            static uint8_t *value_cache = new uint8_t[8];
            
            string input;
            cin >> input;
            input = formatDNAInput(input);
            FastxRecord record;
            record.header = ">dummy";
            
            bool testfile = false;
            cout << "input: " << input << endl;
            while (input.compare("stop") != 0) {
                if (input.compare("file") == 0) {
                    testfile = true;
                    break;
                }
                cout << "format: " << input << endl;
                input = formatDNAInput(input);
                cout << "query: " << input << endl;
                if (input.length() >= k + meta_db.shape.length()) {
                    record.sequence = input;
                    
                    iterator.setRecord(record);
                    
                    cout << input << endl;
                    cout << shape_str << endl;
                    
                    while (iterator.hasNext()) {
                        iterator.operator() (key);
                        cout << "key: " << KmerUtils::bytesToDNA(key, k) << ": ";
                        
                        uint64_t tid = map.search(key);
                        *((uint64_t*)value_cache) = tid;
                        Utils::swapEndianess(value_cache, 3);
                        cout << *((uint64_t*)value_cache) << endl;
                        //uint64_t hash = map.Hash(key);
                        //map.printMapPSL(hash, 5);
                    }
                }
                cin >> input;
                
            }
            
            if (testfile) {
                cout << "input filepath" << endl;
                cin >> input;
    
//                // Load Taxonomy
//                NCBITaxonomy taxonomy = NCBITaxonomy();
//                taxonomy.loadCustomNodes(db + MetaDataDB::TAX_NODES_FILE);
//                taxonomy.loadCustomNames(db + MetaDataDB::TAX_NAMES_FILE);
                
                BufferedFastxReader reader = BufferedFastxReader();
                std::istream* is = new std::ifstream(input);
                uint64_t value = 1;
                uint64_t tax_id = 0;
                
                while (true) {
                    bool ok = false;
                    ok = reader.LoadBlock(*is, 200000);
                    if (!ok) break;
                    
                    while (true) {
                        auto valid_fragment = reader.NextSequence(record);
                        if (!valid_fragment) break;
                        iterator.setRecord(record);
                        cout << record.to_string() << endl;
    
//                        // Extract taxonomic identifier from sequence header (has to be a ncbi identifier e.g.: >813)
//                        value = stoll(record.header.substr(1));
//                        tax_id = taxonomy.getCustom(value);
//                        if (tax_id == -1) {
//                            cerr << "unknown taxid: " << tax_id << " for value " << value << endl;
//                            continue;
//                        }
    
                        while (iterator.hasNext()) {
                            // Extract key from k-mer
                            iterator.operator()(key);
                            cout << map.Hash(key) << endl;
                            cout << KmerUtils::bytesToDNA(key, 22) << ": ";
                            uint64_t tid = map.search(key);
                            *((uint64_t*)value_cache) = tid;
                            Utils::swapEndianess(value_cache, 3);
                            cout << *((uint64_t*)value_cache) << endl;
                            value = map.search(key);
                            if (*((uint64_t*)value_cache) == 0) {
                                cerr << "missing k-mer at pos in sequence" << iterator.getPos() << endl;
                                cerr << iterator.getSubstring(iterator.getPos()-1, 100);
                                cout << record.to_string();
                                exit (50);
                            }
                        }
                    }
                }
            }
            
            delete[] key;
            delete[] value_cache;
            
        } else if (string(argv[1]) == test_mode) {
//            std::string folder = "/home/joachim/CLionProjects/varkit/data/test/";
//            std::string file = folder + "gzip_test.txt";
//            std::string out = folder + "gzip_test.txt.gz";
//            std::string unzip_out = folder + "gzip_test2_double.txt";
//            std::string test = folder + "test.txt";
//
//            std::string dna_in = folder + "ultra_test.fa";
//            std::string dna_out = folder + "test.fa.gz";
            
//            test_key();
//            test_kmer_buffer();
//            cout << "zip" << endl;
//            gzip_test(file.c_str(), out.c_str());
//            cout << "done zippin" << endl;
//            gunzip_test(out.c_str(), unzip_out.c_str());
            
//            cout << file << endl;
//            read_auto(file.c_str());
//            cout << out << endl;
//            read_auto(out.c_str());
//            cout << unzip_out << endl;
//            read_auto(unzip_out.c_str());
//            cout << "auto detect: " << test << endl;
//            read_auto(test.c_str());
//            cout << "txt only: " << test << endl;
//            read_unzipped(test.c_str());
//            cout << "getline: " << test << endl;
//            read_getline(test.c_str());
            
            //gzip_test(dna_in.c_str(), dna_out.c_str());
            
//            cout << "read" << endl;
//            read_auto(dna_out.c_str());

            string nodes =  "/media/joachim/TOSHIBA EXT/Bioinformatics/gtdb_taxonomy.dmp";

            string subset_path = "/home/joachim/CLionProjects/varkit/data/taxonomy/gtdb/taxa.txt";
            string names_out = "/home/joachim/CLionProjects/varkit/data/taxonomy/gtdb/names.dmp";
            string nodes_out = "/home/joachim/CLionProjects/varkit/data/taxonomy/gtdb/nodes.dmp";
            
            GTDBTaxonomy taxonomy;
            taxonomy.loadNodes(nodes);
            taxonomy.subsetByTaxa(subset_path);
            taxonomy.saveCustomNodes(nodes_out);
            taxonomy.saveCustomNames(names_out);
            
            
        } else {
            cout << "'" << string(argv[1]) << "' is not a valid mode. Please provide a valid mode." << endl << endl;
            printGeneralHelp(min_dist, max_cols);
        }
        
        delete[] argv_copy;
    }

    
    return 0;
}
#endif

