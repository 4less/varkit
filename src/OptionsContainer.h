//
// Created by joachim on 06/07/2020.
//
#pragma once

#include <compact_map/compact_map.h>
#include <taxonomy/TaxonomyNew.h>
#include <taxonomy/ReadClassifier.h>
#include "io/MetaDataDB.h"
#include "filesystem"
#include <regex>
#include "FastxReader.h"

class VarkitOptionsContainer;

struct Sample {
private:
    std::string se_read_path = "";
    std::string pe1_read_path = "";
    std::string pe2_read_path = "";

    std::string base_name = "";
    size_t index = 0;

    void CheckFile(std::string file_path, std::string read_type) {
        errx(EX_IOERR,
             "%s file %s does not exist.", read_type.c_str(), file_path.c_str());
    }

public:
    Sample(std::string se_read_path, bool check=false) :
                se_read_path(se_read_path) {

        base_name = Utils::GetBasename(se_read_path);

        std::cout << "basename: " << base_name << " (" << se_read_path << ")" << std::endl;

        if (!Utils::HasSuffix(base_name, ".fastq")  && !Utils::HasSuffix(base_name, ".fq")) {
            errx(EX_IOERR, "Only fastq files are accepted. %s does not end on .fastq or .fq", se_read_path.c_str());
        }
        base_name = Utils::RemoveSuffix(base_name, ".fastq");
        base_name = Utils::RemoveSuffix(base_name, ".fq");

        if (!check) return;

        CheckFile(se_read_path, "Single-end");
    }

    Sample() {};

    Sample(std::string pe1_read_path, std::string pe2_read_path, bool check=false, bool use_first_as_base=false) :
                pe1_read_path(pe1_read_path),
                pe2_read_path(pe2_read_path) {

        auto filename1 = Utils::GetBasename(pe1_read_path);
        auto filename2 = Utils::GetBasename(pe2_read_path);

        if (!Utils::HasSuffix(filename1, ".fastq")  && !Utils::HasSuffix(filename1, ".fq")) {
            errx(EX_IOERR, "Only fastq files are accepted. %s does not end on .fastq or .fq", pe1_read_path.c_str());
        }
        if (!Utils::HasSuffix(filename2, ".fastq")  && !Utils::HasSuffix(filename2, ".fq")) {
            errx(EX_IOERR, "Only fastq files are accepted. %s does not end on .fastq or .fq", pe2_read_path.c_str());
        }

        filename1 = Utils::RemoveSuffix(filename1, ".fastq");
        filename1 = Utils::RemoveSuffix(filename1, ".fq");
        filename2 = Utils::RemoveSuffix(filename2, ".fastq");
        filename2 = Utils::RemoveSuffix(filename2, ".fq");

        if (use_first_as_base) {
            base_name = filename1;
        } else {
            base_name = Utils::GetSharedPrefix(filename1, filename2);
            if (base_name[base_name.size()-1] == '_')
                base_name = Utils::RemoveSuffix(base_name, "_");
        }

        if (base_name.empty()) {
            std::cerr << "Warning: fastq files do not share common prefix." << std::endl;
            std::cerr << pe1_read_path << " " << pe2_read_path << std::endl;
            base_name = filename1;
        }
//        std::cout << "basename1: " << filename1 << " (" << pe1_read_path << ")" << std::endl;
//        std::cout << "basename2: " << filename2 << " (" << pe2_read_path << ")" << std::endl;
//        std::cout << "basename: " << base_name << std::endl;

        if (!check) return;

        CheckFile(pe1_read_path, "Paired-end (1)");
        CheckFile(pe2_read_path, "Paired-end (2)");
    }

    Sample(std::string se_read_path, std::string pe1_read_path, std::string pe2_read_path, bool check=false) :
                se_read_path(se_read_path),
                pe1_read_path(pe1_read_path),
                pe2_read_path(pe2_read_path) {
        if (!check) return;

        CheckFile(se_read_path, "Single-end");
        CheckFile(pe1_read_path, "Paired-end (1)");
        CheckFile(pe2_read_path, "Paired-end (2)");
    }

    std::pair<std::string, std::string> GetPairedEndPaths() {
        return { pe1_read_path, pe2_read_path };
    }

    std::string GetSingleEndPath() {
        return se_read_path;
    }

//    std::string GetDirectory() {
//        return GetBasenameNoExtension()
//    }

    bool HasSingleEnd() {
        return !se_read_path.empty();
    }

    void SetBasename(std::string base_name) {
        this->base_name = base_name;
    }

    std::string GetBasenameNoExtension() const {
        return base_name;
    }

    void SetSingleEnd(std::string se_read_path) {
        this->se_read_path = se_read_path;
    }

    bool HasPairedEnd() {
        return !pe1_read_path.empty() && !pe2_read_path.empty();
    }

    void SetIndex(size_t index) {
        this->index = index;
    }

    size_t GetIndex() const {
        return index;
    }

    std::string GetSampleName() const {
        return GetBasenameNoExtension();
    }

    std::string GetSampleString() {
        std::string result = "Sample: \n";

        size_t col_width = 20;

        if (HasSingleEnd()) {
            Utils::exists(se_read_path);
            std::string exists = "";
            result += Utils::PadString("Single-end:", col_width, ' ');
            result += (Utils::exists(se_read_path) ? "X" : " ");
            result += "    " + se_read_path + '\n';
        }
        if (HasPairedEnd()) {
            Utils::exists(pe1_read_path);
            result += Utils::PadString("Paired-end (1):", col_width, ' ');
            result += (Utils::exists(pe1_read_path) ? "X" : " ");
            result += "   " + pe1_read_path + '\n';
            result += Utils::PadString("Paired-end (2):", col_width, ' ');
            result += (Utils::exists(pe1_read_path) ? "X" : " ");
            result += "   " + pe2_read_path + '\n';
        }
        result += "Basename: " + base_name;

        return result;
    }
};

class VarkitOptionsContainer {
public:
    const std::string SHAPE_FILE = "/shape.txt";
    const std::string SHAPE_DB_FILE = "/shape.db";
    const std::string META_FILE = "/entry.meta";
//    const std::string INDEX_META_FILE = "/db.meta";
    const std::string DATABASE_FILE = "/database.bin";
    const std::string TAXONOMY_FILE = "/internal_taxonomy.dmp";
    const std::string LIBRARY_FOLDER = "/library/";
    const std::string TAXONOMY_FOLDER = "/taxonomy/";
    static inline const std::string BUCKETS_FILE = "/buckets.bin";
    static inline const std::string MARKER_GENE_LENGTHS_FILE = "/mg_lengths.tsv";
    static inline const std::string ABUNDANCE_CORRECTION_FILE = "/abundance_correction.tsv";

    static inline const std::string CLASS_EXTENSION = ".class.bin";
    static inline const std::string PROFILE_EXTENSION = ".profile.tsv";
    static inline const std::string SNPS_EXTENSION = ".snps.bin";
    static inline const std::string STRAIN_OUTPUT = "/strain_distance/";

    static inline const std::string REFERENCE_FILE = "/reference.fna";
    static inline const std::string REFERENCE_MAP_FILE = "/reference.map";


//    std::string index_folder = "";
//    std::string index_name = "";
//    std::string default_name = "";

//    std::unordered_map<std::string, std::string> indices;

    size_t threads_ = 0;
    size_t memory_upper_gb_ = 8;
    std::string output_dir_ = "";
    std::string output_prefix_ = "";
    std::string db_dir_ = "";


    std::vector<std::string> reads_;
    std::vector<std::string> reads1_;
    std::vector<std::string> reads2_;

    std::vector<Sample> samples_;

    std::vector<std::string> output_;
    std::vector<std::string> references_;


    double load_factor_ = -1.;

    bool force_rebuild_ = false;
    bool raw_taxa = false;

    size_t k_ = 0;
    size_t key_bits_ = 0;
    size_t value_bits_ = 0;

    size_t offset_key_bits_ = 0;
    size_t bucket_count_ = 0;
    size_t internal_key_bits_ = 0;

    size_t internal_taxid_bits_ = 0;
    size_t geneid_bits_ = 0;
    size_t genepos_bits_ = 0;

    size_t max_marker_gene_pos_ = 0;

    std::string shape_str = "";
    bool* shape;

    IndexedMap *map = nullptr;
    Taxonomy::IntTaxonomy *taxonomy = nullptr;

    void LoadResources(bool force=false) {
        std::cout << "Load resources" << std::endl;
        if (InternalTaxonomyFileExists()) {
            LoadInternalTaxonomy();
        }

        if (HasDatabase() && (!IsDatabaseLoaded() || force)) {
            std::cout << "Load " << DatabaseFile() << std::endl;
            map = IndexedMap::Load(DatabaseFile());
            if (offset_key_bits_ != map->offset_bits_) {
                std::cerr << "Offset bits are different" << std::endl;
                exit(9);
            }
        }
        LoadPatternDB();
    }

    static bool IsInterlaced(std::string& read, size_t limit=0) {
        FastxRecord record1;
        FastxRecord record2;
        BufferedFastxReader reader;

        if (!Utils::exists(read)) {
            std::cerr << "______________________________" << std::endl;
            std::cerr << "Read: " << read << std::endl;
            std::cerr << "Read file does not exist" << std::endl;
            exit(12);
        }

        ifstream is(read, ios::in);
        const size_t batch_size = 1024;

        if (limit == 0) limit = UINT64_MAX;
        int i = 0;
        while (true) {
            bool ok = false;
            ok = reader.LoadBatch(is, batch_size);
            if (!ok) break;

            for (; i < limit; i++) {
                auto valid_fragment = reader.NextSequence(record1);
                if (!valid_fragment) break;
                auto valid_fragment2 = reader.NextSequence(record2);
                if (!valid_fragment2) return false;

                if (!IsInterlaced(record1.id, record2.id)) {
                    return false;
                }
            }
            if (i == limit) break;
        }

        return true;
    }

    static bool IsInterlaced(std::string& header1, std::string& header2) {
        if (header1.length() != header2.length()) return false;

        std::string_view prefix1(header1.c_str(), header1.length() - 2);
        std::string_view prefix2(header2.c_str(), header1.length() - 2);

        if (prefix1 != prefix2) return false;

        std::string_view suffix1(header1.c_str() + prefix1.length(), 2);
        std::string_view suffix2(header2.c_str() + prefix2.length(), 2);

        if (suffix1 == "/1" && suffix2 == "/2") {
            return true;
        }

        return false;
    }

    static bool IsPair(std::string read1, std::string read2) {
        std::smatch m1;
        std::smatch m2;
        std::regex e1 ("[_|\\.]1\\.");
        std::regex e2 ("[_|\\.]2\\.");

        std::string name1 = Utils::strip_path(read1);
        std::string name2 = Utils::strip_path(read2);

        std::string prefix1 = "";
        string::const_iterator searchStart1( read1.cbegin() );

        while ( std::regex_search( searchStart1, read1.cend(), m1, e1 ) ) {
            prefix1 += m1.prefix();
            std::string suffix1 = m1.suffix();

            std::string prefix2 = "";
            string::const_iterator searchStart2( read2.cbegin() );
            while ( std::regex_search( searchStart2, read2.cend(), m2, e2 ) ) {
                prefix2 += m2.prefix();
                std::string suffix2 = m2.suffix();

                if (suffix1 == suffix2 && prefix1 == prefix2)
                    return true;

                prefix2 += m2[0];
                searchStart2 = m2.suffix().first;
            }
            prefix1 += m1[0];
            searchStart1 = m1.suffix().first;
        }

        return false;
    }

    static std::vector<std::string> GetHeaders(std::string read, size_t limit) {
        std::vector<std::string> headers;

        FastxRecord record;
        BufferedFastxReader reader;
        ifstream is(read, ios::in);
        const size_t batch_size = 1024;

        if (limit == 0) limit = UINT64_MAX;
        int i = 0;
        while (true) {
            bool ok = false;
            ok = reader.LoadBatch(is, batch_size);
            if (!ok) break;

            for (; i < limit; i++) {
                auto valid_fragment = reader.NextSequence(record);
                if (!valid_fragment) break;

                headers.push_back(record.id);
            }
            if (i == limit) break;
        }

        return headers;
    }

    void LoadReadstrings(std::vector<std::string> &reads, bool auto_sort=false) {
        if (IsInterlaced(reads[0])) {
            std::cerr << "Interlaced reads are currently not supported. Please deinterlace your reads." << std::endl;
            std::cerr << "Interlaced: " << reads[0] << std::endl;

            exit(2);
        }

        if (reads.size() == 1) {
            std::cout << "Singled end files" << std::endl;
            samples_.push_back(Sample(reads[0]));
            return;
        }

        if ((reads.size() % 2) != 0) {
            exit(3);
        }

        if (auto_sort) {
            std::sort(reads.begin(), reads.end());
        }


        reads_.clear();
        reads2_.clear();
        for (int i = 0; i < reads.size(); ++i) {
            //TODO: fix
            if (IsPair(reads[i], reads[i+1])) {
//                std::cout << "paired? " << std::endl;
//                std::cout << "One: " << reads[i] << std::endl;
//                std::cout << "Two: " << reads[i+1] << std::endl;
                Sample sample(Sample(reads[i], reads[i+1]));
                samples_.push_back(sample);
                i++;
            } else {
                Sample sample(reads[i]);
                samples_.push_back(sample);
            }
        }

    }

    void SetCustomBasenames(std::vector<std::string> &output_names) {
        if (output_names.size() != samples_.size()) {
            errx(EX_IOERR, "If custom basenames are provided, m_sample count and basenames must be the same.");
        }
        for (auto i = 0; i < samples_.size(); i++) {
            samples_[i].SetBasename(output_names[i]);
        }
    }
    void Init() {
//        LoadIndices();
        LoadShape();
        LoadMeta();

        // Set valuebits if not set.
        if (offset_key_bits_ == 0) {
            // offset_key_bits not set..
            std::cerr << "Warning: offset key bits not set. Can lead to undefined behaviour." << std::endl;
            int missing_bits = key_bits_ + value_bits_ - 64;
            offset_key_bits_ = missing_bits > 0 ? missing_bits : 0;
        } else {
            value_bits_ = 64 - internal_key_bits_;
        }

        bucket_count_ = 1llu << offset_key_bits_;
        max_marker_gene_pos_ = (1 << genepos_bits_) - 1;
    }

//    void OldInit(std::string taxonomy_name="") {
////        LoadIndices();
//        LoadShape();
//        LoadMeta();
//
//        if (taxonomy_name.empty()) {
//            ChangeDatabase(default_name);
//        } else {
//            std::cout << "Use database " << taxonomy_name << std::endl;
//            ChangeDatabase(taxonomy_name);
//        }
//
//        // Set valuebits if not set.
//        if (offset_key_bits_ == 0) {
//            // offset_key_bits not set..
//            std::cerr << "Warning: offset key bits not set. Can lead to undefined behaviour." << std::endl;
//            int missing_bits = key_bits_ + value_bits_ - 64;
//            offset_key_bits_ = missing_bits > 0 ? missing_bits : 0;
//        } else {
//            value_bits_ = 64 - internal_key_bits_;
//        }
//
//        bucket_count_ = 1llu << offset_key_bits_;
//        max_marker_gene_pos_ = (1 << genepos_bits_) - 1;
//    }

//    void LoadIndices() {
//        if (DBFolderExists() && IndexFileExists()) {
//            std::string index_file = IndexFile();
//
//            std::ifstream is(index_file.c_str());
//
//            std::string line;
//            std::vector<std::string> tokens;
//
//            bool is_default = false;
//
//            while (getline(is, line)) {
//                if (line.starts_with("#")) {
//                    if (line.erase(0, 1) == "DEFAULT_INDEX") {
//                        is_default = true;
//                    }
//                    continue;
//                }
//
//                Utils::split(tokens, line, "=");
//                if (tokens.size() != 2) {
//                    std::cout << "Error in reading meta file" << std::endl;
//                    exit(8);
//                }
//
//                auto key = tokens[0];
//                auto value = tokens[1];
//
//                if (is_default) {
//                    default_name = key;
//                    is_default = false;
//                }
//
//                indices.insert ( { key, value } );
//            }
//
//            is.close();
//        } else {
//            std::cerr << "Index file does not exist." << std::endl;
//            exit(9);
//        }
//    }
//
//    void ListIndices() {
//        if (!indices.empty()) {
//            for (auto& pair : indices) {
//                std::cout << pair.first;
//                if (pair.first == default_name) {
//                    std::cout << "*";
//                }
//                std::cout << "\t" << pair.second << std::endl;
//            }
//        } else {
//            std::cout << "No databases detected. " << std::endl;
//        }
//    }

    bool HasPatternDB() {
        return pattern_db != nullptr;
    }

    void LoadPatternDB() {
//        pattern_db = new PatternDB(ShapeDatabaseFile());
        pattern_db = new ds::PatternMap(ShapeDatabaseFile());
    }

    void LoadMeta() {
        if (Utils::exists(MetaFile())) {
            std::ifstream is(MetaFile());
            std::string line;

            if (getline(is, line)) {
                internal_taxid_bits_ = stoull(line);
            }
            if (getline(is, line)) {
                geneid_bits_ = stoull(line);
            }
            if (getline(is, line)) {
                genepos_bits_ = stoull(line);
            }
            if (getline(is, line)) {
                offset_key_bits_ = stoull(line);
            }
            is.close();

//            value_bits_ = internal_taxid_bits_ + geneid_bits_ + genepos_bits_;
//            value_bits_ = 64 - internal_taxid_bits_;

        } else {
            std::cerr << "Warning: Meta file does not exist." << std::endl;
        }
    }

    void LoadShape() {
        if (Utils::exists(ShapeFile())) {
            this->shape_str = ShapeUtils::LoadShape(ShapeFile());
            this->shape = ShapeUtils::GetShape(this->shape_str);
            this->k_ = ShapeUtils::GetK(this->shape_str);
            this->key_bits_ = this->k_ * 2;

            internal_key_bits_ = key_bits_ - offset_key_bits_;


            value_bits_ = 64 - internal_key_bits_;


        } else {
            std::cerr << "Warning" << std::endl;
        }
    }
//
//    std::string OldInternalTaxonomyFile() {
//        return db_dir_ + '/' + index_folder + '/' + TAXONOMY_FILE;
//    }

    std::string InternalTaxonomyFile() {
        return db_dir_ + '/' + TAXONOMY_FOLDER + '/' + TAXONOMY_FILE;
    }

    std::vector<Sample>& GetSamples() {
        return samples_;
    }

    void LoadInternalTaxonomy() {
        std::cout << "Load internal taxonomy file:  " << InternalTaxonomyFile() << std::endl;
        taxonomy = new Taxonomy::IntTaxonomy(InternalTaxonomyFile());
    }


    Taxonomy::IntTaxonomy* InternalTaxonomy() {
        return taxonomy;
    }

//    void SaveIndices() {
//        std::ofstream os(MetaFile());
//
//        for (auto kv : indices) {
//            auto name = kv.first;
//            auto path = kv.second;
//
//            if (name == default_name) {
//                os << "#DEFAULT_INDEX" << std::endl;
//            }
//            os << name << "=" << path << std::endl;
//        }
//
//        os.close();
//    }

    bool IsValidCellBits() {
        return false;
    }

    bool DBFolderExists() {
        if (!filesystem::is_directory(db_dir_)) {
            // Exception
            std::cerr << "DB folder does not exist." << std::endl;
            exit(9);
        }
        return true;
    }

    bool IsMapInitialized() {
        return map != nullptr;
    }

    void SetMap(IndexedMap *map) {
        this->map = map;
    }

//    void ChangeDatabase(std::string taxonomy) {
//        if (!indices.empty()) {
//            if (indices.contains(taxonomy)) {
//                index_name = taxonomy;
//                index_folder = indices[index_name];
//            } else {
//                std::cerr << "Taxonomy " << taxonomy << " not found." << std::endl;
//            }
//        }
//    }

    std::string MetaFile() const {
        return db_dir_ + META_FILE;
    }

//    std::string OldDatabaseFile() const {
//        return db_dir_ + '/' + index_folder + '/' + DATABASE_FILE;
//    }

    std::string DatabaseFile() const {
        return db_dir_ + '/' + DATABASE_FILE;
    }

    std::string ReferenceFile() const {
        return db_dir_ + '/' + LIBRARY_FOLDER + '/' + REFERENCE_FILE;
    }
    std::string ReferenceMapFile() const {
        return db_dir_ + '/' + LIBRARY_FOLDER + '/' + REFERENCE_MAP_FILE;
    }
//
//    std::string IndexFile() const {
//        return db_dir_ + INDEX_META_FILE;
//    }

    std::string ShapeFile() const {
        return db_dir_ + SHAPE_FILE;
    }

    std::string ShapeDatabaseFile() const {
        return db_dir_ + SHAPE_DB_FILE;
    }

    std::string ClassificationOutputFile() const {
        return output_prefix_ + ".csv";
    }

    std::string MarkerGeneLengthsFile() const {
        return db_dir_ + '/' + MARKER_GENE_LENGTHS_FILE;
    }

    std::string OutputPrefix() const {
        return output_prefix_;
    }

    std::string ClassificationOutputFile(Sample &sample) const {
        return output_dir_ + '/' + sample.GetBasenameNoExtension() + CLASS_EXTENSION;
    }

    std::string ProfileOutputFile(Sample& sample) const {
        return output_dir_ + '/' + sample.GetBasenameNoExtension() + PROFILE_EXTENSION;
    }

    std::string RawTaxaOutputFile(Sample& sample) const {
        return output_dir_ + '/' + sample.GetBasenameNoExtension() + ".taxa.tsv";
    }

    std::string OutputDir() const {
        return output_dir_ + '/';
    }

    std::string StrainOutputDir() const {
        return output_dir_ + '/' + STRAIN_OUTPUT;
    }

    string ProfileOutputFileFromClass(string classification_file) const {
        auto basename = Utils::RemoveSuffix(classification_file, CLASS_EXTENSION);
        return basename + PROFILE_EXTENSION;
    }

    std::string SNPOutputFile(Sample& sample) const {
        return output_dir_ + '/' + sample.GetBasenameNoExtension() + SNPS_EXTENSION;
    }
//
//    std::string AbundanceOutputPrefix() {
//        return output_prefix_ + ".abundance";
//    }

    std::string AbundanceCorrectionFile() const {
        return db_dir_ + "/abundance_correction.tsv";
    }

    std::string SNPOutputFile() {
        return output_prefix_ + SNPS_EXTENSION;
    }

//    std::string OldBucketsFile() const {
//        return db_dir_ + '/' + index_folder + '/' + BUCKETS_FILE;
//    }

    std::string BucketsFile() const {
        return db_dir_ + '/' + BUCKETS_FILE;
    }


    bool* Shape() {
        return shape;
    }

    size_t ReferencesTotalSize() const {
        // First estimate bucket size
        size_t total_size = 0;

        for (auto file : references_) {
            total_size += Utils::GetFileSize(file);
        }
        return total_size;
    }

    size_t ReferenceSize() const {
        return Utils::GetFileSize(ReferenceFile());
    }

    size_t ReadsTotalSize() {
        // First estimate bucket size
        size_t total_size = 0;

        for (auto file : reads_) {
            total_size += Utils::GetFileSize(file);
        }
        return total_size;
    }

    size_t ShapeLength() {
        return shape_str.size();
    }

    bool MetaFileExists() {
        return Utils::exists(MetaFile());
    }

//    bool IndexFileExists() {
//        return Utils::exists(IndexFile());
//    }

    bool InternalTaxonomyFileExists() {
        return Utils::exists(InternalTaxonomyFile());
    }

    bool IsShapeValid() {
        return false;
    }

    bool IsValidClassifyOptions() {
        if (!std::filesystem::is_directory(db_dir_)) {
            errx(EX_IOERR,
                 "Database folder %s is no directory or does not exist.", db_dir_.c_str());
        }
        if (!Utils::exists(db_dir_)) {
            errx(EX_IOERR,
                 "Database file %s does not exist.", db_dir_.c_str());
        }
        if (!MetaFileExists()){
            errx(EX_IOERR,
                 "Meta file %s does not exist. Abort.", MetaFile().c_str());
        }
        if (!ReadsExist()) {
            errx(EX_IOERR,
                 "Reads do not exist. Abort.");
        }
        if (!InternalTaxonomyFileExists()) {
            errx(EX_IOERR,
                 "Taxonomy file %s do not exist.", InternalTaxonomyFile().c_str());
        }
        if (!Utils::exists(MarkerGeneLengthsFile())) {
            errx(EX_IOERR,
                 "Marker gene lengths file %s do not exist.", MarkerGeneLengthsFile().c_str());
        }
        if (!HasShapeDatabase()) {
            errx(EX_IOERR,
                 "Shape Database file %s do not exist.", ShapeDatabaseFile().c_str());
        }
        if ((64 - internal_key_bits_) < value_bits_) {
            errx(EX_IOERR,
                 "Value bits %lu and internal key bits %lu are bigger than cell size 64.", value_bits_, internal_key_bits_);
        }
        return true;
    }

    bool IsValidReads() {
        return false;
    }

    bool ReferencesExist() {
        if (references_.empty()) {
            std::cerr << "References are empty." << std::endl;
        }
        for (auto& ref : references_) {
            if (!Utils::exists(ref)) {
                std::cerr << "Reference " << ref << " does not exist." << std::endl;
                return false;
            }
            if (Utils::GetFileSize(ref) == 0) {
                std::cerr << "File is empty" << std::endl;
                return false;
            }
        }
        return true;
    }

    bool ReferenceExist() {
        return Utils::exists(ReferenceFile());
    }

    bool ReadsExist() {
        for (auto& read : reads_) {
            if (!Utils::exists(read)) {
                std::cerr << "Read " << read << " does not exist." << std::endl;
                return false;
            }
            if (Utils::GetFileSize(read) == 0) {
                std::cerr << "File is empty" << std::endl;
                return false;
            }
        }
        return true;
    }

    bool IsValidReferences() {
        if (references_.empty()) {
            std::cerr << "References are empty. Abort" << std::endl;
            exit(12);
        }

        for (auto& ref : references_) {
            if (!Utils::exists(ref)) {
                return false;
            }
            if (Utils::GetFileSize(ref) == 0) {
                return false;
            }
        }
        return true;
    }

    bool IsValidDBPreBuild() {
        if (!std::filesystem::is_directory(db_dir_)) {
            errx(EX_IOERR,
                 "Database folder %s is no directory or does not exist.", db_dir_.c_str());
        }
        if (!Utils::exists(db_dir_)) {
            errx(EX_IOERR,
                 "Database file %s does not exist.", db_dir_.c_str());
        }
        if (!MetaFileExists()){
            errx(EX_IOERR,
                 "Meta file %s does not exist. Abort.", MetaFile().c_str());
        }
//        if (!ReferencesExist()) {
//            errx(EX_IOERR,
//                 "References do not exist. Abort.");
//        }
        if (!ReferenceExist()) {
            errx(EX_IOERR,
                 "Reference file does not exist. Abort.");
        }
        if (!InternalTaxonomyFileExists()) {
            errx(EX_IOERR,
                 "Taxonomy file %s do not exist.", InternalTaxonomyFile().c_str());
        }
        if ((64 - internal_key_bits_) < value_bits_) {
            errx(EX_IOERR,
                 "Value bits %lu and internal key bits %lu are bigger than cell size 64.", value_bits_, internal_key_bits_);
        }
        return true;
    }

    bool IsValidDBPostBuild() {
        if (!std::filesystem::is_directory(db_dir_)) {
            errx(EX_IOERR,
                 "Database folder %s is no directory or does not exist.", db_dir_.c_str());
        }
        if (!Utils::exists(db_dir_)) {
            errx(EX_IOERR,
                 "Database file %s does not exist.", db_dir_.c_str());
        }
        if (!MetaFileExists()){
            errx(EX_IOERR,
                 "Meta file %s does not exist. Abort.", MetaFile().c_str());
        }
        if (!ReadsExist()) {
            errx(EX_IOERR,
                 "Reads do not exist. Abort.");
        }
        if (!InternalTaxonomyFileExists()) {
            errx(EX_IOERR,
                 "Taxonomy file %s do not exist.", InternalTaxonomyFile().c_str());
        }
        if ((64 - internal_key_bits_) < value_bits_) {
            errx(EX_IOERR,
                 "Value bits %lu and internal key bits %lu are bigger than cell size 64.", value_bits_, internal_key_bits_);
        }
        return true;
    }

    bool HasBuckets() {
        return Utils::exists(BucketsFile());
    }

    bool HasDatabase() {
        return Utils::exists(DatabaseFile());
    }

    bool IsDatabaseLoaded() {
        return map != nullptr;
    }

    bool HasShapeDatabase() {
        return Utils::exists(ShapeDatabaseFile());
    }

    bool IsPairedEnd() {
        return !reads2_.empty();
    }

//    PatternDB *pattern_db = nullptr;
    ds::PatternMap *pattern_db = nullptr;

};


class ClassifyOptionsContainer {
public:
    const MetaDataDB meta_db;
    const int threads;
//    const string read;
    const std::vector<std::string> reads;
    const string output;
    const size_t vmer_length;
    
    ClassifyOptionsContainer(MetaDataDB &meta_db, std::vector<std::string> reads, int threads, string output, size_t vmer_length) :
            meta_db(meta_db),
            reads(reads),
            threads(threads),
            output(output),
            vmer_length(vmer_length) {};
};

class BuildOptionsContainer {
public:
    const int threads;
    const string db;
    const int64_t initial_capacity;
    const string taxonomy;
    const bool validate;
    const vector<string> reference;
    
    BuildOptionsContainer(int threads, string db, vector<string> reference, int64_t initial_capacity, string taxonomy, bool validate) :
            reference(reference),
            db(db),
            threads(threads),
            initial_capacity(initial_capacity),
            taxonomy(taxonomy),
            validate(validate) {};
};
