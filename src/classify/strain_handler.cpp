//
// Created by fritsche on 03/07/2022.
//

#include "strain_handler.h"


bool varkit::SequenceRange::Overlap(size_t start, size_t end) const {
    // MS###S#####ME Start between m_start and m_end or
    // MS##E######ME End between m_start and m_end
    // Start is inclusive, end is exclusive
    return start >= m_start && start <= m_end || end >= m_start && end <= m_end || m_start >= start && m_start <= end || m_end >= start && m_end <= end;
}

void varkit::SampleHandler::AddSample(std::string sample_name) {
//    m_samples.insert( { sample_name, SampleStrains{ } } );
}

void varkit::SampleHandler::AddSample(Sample &sample) {
    m_samples.insert( {sample.GetBasenameNoExtension(), SampleStrains(sample) } );
}

varkit::SampleStrains &varkit::SampleHandler::GetSample(std::string sample_name) {
    auto it = m_samples.find(sample_name);
    if (it == m_samples.end()) {
        errx(EX_IOERR, "GetSample() in SampleHandler in strain_handler.h - sample_name not present");
    }
    return it.value();
}

varkit::StrainSet varkit::SampleHandler::GetStrainList() const {
    StrainSet strain_keys;
    for (const auto &[sample_key, sample] : m_samples) {
        for (const auto& [strain_key, strain] : sample) {
            strain_keys.insert(strain_key);
        }
    }
    return strain_keys;
}

varkit::StrainSet varkit::SampleHandler::GetPredictedStrainList() const {
    StrainSet strain_keys;
    for (const auto &[sample_key, sample] : m_samples) {
        for (const auto& strain_key : sample.GetPredictedStrains()) {
            strain_keys.insert(strain_key);
        }
    }
    return strain_keys;
}



void varkit::SampleHandler::SetGenomeLoader(varkit::SampleHandler::GenomeLoaderPtr genome_loader) {
    m_loader = genome_loader;
}

void varkit::SampleHandler::SetInternalTaxonomy(InternalTaxonomyPtr taxonomy) {
    m_taxonomy = taxonomy;
}


void varkit::SampleHandler::CompareSNPsToReference() {
    if (!m_loader) return;
    for (auto &[sample_key, sample] : m_samples) {
        for (auto &[strain_key, strain] : sample) {
            m_loader->GetGenome(strain_key).LoadGenome();
            auto& genome = m_loader->GetGenome(strain_key);

            for (auto& [gene_key, gene] : strain.GetGenes().GetGeneMap()) {
                auto& lgene = genome.GetGene(gene_key);
                gene.CompareSNPsToReference(genome.GetGene(gene_key).Sequence());
            }
        }
    }
}

varkit::SampleHandler::GenomeLoaderPtr varkit::SampleHandler::GetGenomeLoader() {
    return m_loader;
}

//varkit::SampleStrains &varkit::SampleHandler::GetOrAddSample(std::string sample_name) {
//    // [] operator for map creates entry if it does not exist as [] operator
//    // returns a reference.
//    return m_samples[sample_name];
//}


/* ####################################################################################
 * SAMPLE STRAINS
 * ####################################################################################
 */

//void varkit::SampleStrains::LoadSNPs(const std::string path, tsl::sparse_set<size_t>& target_species, tsl::sparse_set<size_t>& valid_read_ids) {
void varkit::SampleStrains::LoadSNPs(const std::string path, tsl::sparse_set<size_t>& target_species, tsl::sparse_set<std::pair<size_t, uint8_t>, Utils::pair_hash>& valid_read_ids) {
    std::ifstream is(path, std::ios::binary);
    IO::IOSNP snp;
    while (IO::SNPReader::Read(is, snp)) {
        if (snp.snp_position == (uint32_t)-1) {
            std::cerr << "Malformed " << snp.ToString() << std::endl;
            continue;
        }
        auto taxonomic_id = snp.tax_id;

        assert(snp.gene_id < 121);
        // TODO Remove in final version
        if (snp.gene_id > 120) exit(12);

        if (!target_species.contains(taxonomic_id)) continue;
        if (!valid_read_ids.contains( { snp.read_id, snp.read_num } )) continue;
//        if (!valid_read_ids.contains( snp.read_id )) continue;

//
//        if (snp.tax_id == 130 && snp.gene_id == 38 && snp.snp_position == 377) {
//            std::cout << snp.ToString() << std::endl;
//        }

        auto& strain = GetStrain(taxonomic_id);
        strain.GetGenes().GetGene(snp.gene_id).AddSNP(snp);

    }

//            for (const auto& [strain_key, strain] : m_strains) {
    for (const auto& strain_key : target_species) {
        auto& strain = GetStrain(strain_key);
        auto& genes = strain.GetGenes();
        for (auto& [gene_key, gene] : genes.GetGeneMap()) {
//                    std::cout << strain_key << " " << gene_key << std::endl;
            gene.SNPsToVector();
//            std::cout << "Genekey: " << gene_key << std::endl;
            gene.AddReferenceCountToSNPs();
        }
    }
}



const tsl::sparse_map<size_t, varkit::Strain>::const_iterator varkit::SampleStrains::begin() const noexcept {
    return m_strains.begin();
}

const tsl::sparse_map<size_t, varkit::Strain>::const_iterator varkit::SampleStrains::end() const noexcept {
    return m_strains.end();
}

varkit::StrainMap &varkit::SampleStrains::GetStrains() {
    return m_strains;
}

const varkit::StrainMap &varkit::SampleStrains::GetStrains() const {
    return m_strains;
}

varkit::Strain &varkit::SampleStrains::GetStrain(size_t internal_id) {
    auto it = m_strains.find(internal_id);
    if (it == m_strains.end()) {
        errx(EX_IOERR, "GetStrain() in SampleStrains in strain_handler.h - Gene id not present");
    }
    return it.value();
}

varkit::Strain &varkit::SampleStrains::GetOrAddStrain(size_t internal_id) noexcept {
    return m_strains[internal_id];
}

const varkit::StrainSet &varkit::SampleStrains::GetPredictedStrains() const {
    return m_predicted_strains;
}

void varkit::SequenceRangeHandler::Add(size_t start, size_t end) {
    SequenceRange query_range(start, end);
    auto range_it = GetSequenceRangeIterator(start, end);

    if (range_it == m_ranges.end() || *range_it != query_range) {
        m_ranges.insert(range_it, query_range);

    } else if (*range_it == query_range) {
        range_it->Union(query_range);
        if (range_it+1 != m_ranges.end() && *(range_it+1) == *range_it) {
            auto del_it = range_it + 1;
            range_it->Union(*del_it);
            m_ranges.erase(del_it);
        }
    }
}

void varkit::SequenceRangeHandler::Add(const IO::ClassificationLine &line) {
    size_t start = std::max(line.genepos, 0);
    size_t end = start + line.overlap_with_gene;

    SequenceRange query_range(start, end);
    ReadInfo rinfo;
    rinfo.read_id = line.record_id;
    rinfo.length = line.overlap_with_gene;
    rinfo.start = start;
    rinfo.forward = line.forward;
    query_range.AddReadInfo(rinfo);
    auto range_it = GetSequenceRangeIterator(start, end);

    if (range_it == m_ranges.end() || *range_it != query_range) {
        m_ranges.insert(range_it, query_range);

    } else if (*range_it == query_range) {
        range_it->Union(query_range);
        if (range_it+1 != m_ranges.end() && *(range_it+1) == *range_it) {
            auto del_it = range_it + 1;
            range_it->Union(*del_it);
            m_ranges.erase(del_it);
        }
    }
}

varkit::SequenceRangeHandler varkit::SequenceRangeHandler::Unify(const varkit::SequenceRangeHandler &handler) const {
    SequenceRangeHandler union_ranges;

    int this_i = 0;
    int other_i = 0;

    while (this_i < m_ranges.size() && other_i < handler.m_ranges.size()) {
        if (m_ranges[this_i] == handler.m_ranges[other_i]) {
            // create new copy
            auto union_range = m_ranges[this_i];
            // Intersect with other
            union_range.Union(handler.m_ranges[other_i]);
            // Push to new intersection handler
            union_ranges.m_ranges.emplace_back(union_range);
        }


        bool increment_this = m_ranges[this_i].m_end < handler.m_ranges[other_i].m_end;
        this_i += increment_this;
        other_i += !increment_this;

    }
    return union_ranges;
}

varkit::SequenceRangeList &varkit::SequenceRangeHandler::GetRanges() {
    return m_ranges;
}

const varkit::SequenceRangeList &varkit::SequenceRangeHandler::GetRanges() const {
    return m_ranges;
}


varkit::SequenceRangeHandler varkit::SequenceRangeHandler::Intersect(const varkit::SequenceRangeHandler &handler) const {
    SequenceRangeHandler intersection;

    int this_i = 0;
    int other_i = 0;

    while (this_i < m_ranges.size() && other_i < handler.m_ranges.size()) {
        if (m_ranges[this_i] == handler.m_ranges[other_i]) {
            // create new copy
            auto new_range = m_ranges[this_i];
            // Intersect with other
            new_range.Intersect(handler.m_ranges[other_i]);
            // Push to new intersection handler
            intersection.m_ranges.emplace_back(new_range);
        }


        bool increment_this = m_ranges[this_i].m_end < handler.m_ranges[other_i].m_end;
        this_i += increment_this;
        other_i += !increment_this;

    }
    return intersection;
}

varkit::SequenceRangeHandler::SRIterator varkit::SequenceRangeHandler::FindSequenceRange(size_t start, size_t end) {
    auto range_it = m_ranges.begin();
    SequenceRange query_range(start, end);

    // ranges are sorted in ascending order by their start.
    // increment iterator until range is smaller
    while (range_it != m_ranges.end() && *range_it < query_range) {
        range_it++;
    }
    if (*range_it == query_range) return range_it;
    return m_ranges.end();
}

std::string varkit::SequenceRangeHandler::ToString() const {
    if (m_ranges.empty()) return "[]";
    std::string str = "[";
    str += m_ranges[0].ToString();
    for (int i = 1; i < m_ranges.size(); i++) {
        str += ", " + m_ranges[i].ToString();
    }
    str += "]";
    return str;
}

std::string varkit::SequenceRangeHandler::ToVerboseString() const {
    if (m_ranges.empty()) return "[]";
    std::string str = "";
    for (int i = 0; i < m_ranges.size(); i++) {
        str += m_ranges[i].ToVerboseString();
    }
    return str;
}


size_t varkit::SequenceRangeHandler::SequenceLength() const {
    return std::accumulate(m_ranges.begin(), m_ranges.end(), 0, [](int a, SequenceRange const& range){ return range.Length(); });
}

size_t varkit::SequenceRangeHandler::Size() const {
    return m_ranges.size();
}

varkit::SequenceRangeHandler::SRIterator varkit::SequenceRangeHandler::begin() noexcept {
    return m_ranges.begin();
}

varkit::SequenceRangeHandler::SRIterator varkit::SequenceRangeHandler::end() noexcept {
    return m_ranges.end();
}

varkit::SequenceRangeHandler::SRIterator
varkit::SequenceRangeHandler::GetOrAddSequenceRangeIterator(const size_t start, const size_t end) {
    auto range_it = m_ranges.begin();
    SequenceRange query_range(start, end);

    // ranges are sorted in ascending order by their start.
    // increment iterator until range is smaller
    while (range_it != m_ranges.end() && *range_it < query_range) {
        std::advance(range_it, 1);
    }
    if (range_it == m_ranges.end() || *range_it != query_range) {
        // Insert new range at end
        m_ranges.insert(range_it, query_range);
    } else if (*range_it == query_range) {
        range_it->Union(query_range);
        if (range_it+1 != m_ranges.end() && *(range_it+1) == *range_it) {
            auto del_it = range_it + 1;
            range_it->Union(*del_it);
            m_ranges.erase(del_it);
        }
    }
    return range_it;
}

varkit::SequenceRangeHandler::SRIterator
varkit::SequenceRangeHandler::GetSequenceRangeIterator(const size_t start, const size_t end) {
    auto range_it = m_ranges.begin();
    SequenceRange query_range(start, end);

    // ranges are sorted in ascending order by their start.
    // increment iterator until range is smaller
    while (range_it != m_ranges.end() && *range_it < query_range) {
//                std::cout << range_it->ToString() << " < " << query_range.ToString() << " " << (query_range > *range_it) << std::endl;
        range_it++;
    }
    return range_it;
}

const varkit::SequenceRange &varkit::SequenceRangeHandler::GetRange(size_t pos) const {
    auto find = std::find_if(m_ranges.begin(), m_ranges.end(), [&pos](SequenceRange const& range){
        return range.Contains(pos);
    });
    //TODO implement exception
    return *find;
}

varkit::SequenceRange &varkit::SequenceRangeHandler::GetRange(size_t pos) {
    auto find = std::find_if(m_ranges.begin(), m_ranges.end(), [&pos](SequenceRange const& range){
        return range.Contains(pos);
    });
    //TODO implement exception
    return const_cast<varkit::SequenceRange&>(*find);
}

bool varkit::SequenceRangeHandler::HasRange(size_t pos) const {
    auto find = std::find_if(m_ranges.begin(), m_ranges.end(), [&pos](SequenceRange const& range) {
        return range.Contains(pos);
    });
    //TODO implement exception
    return find != m_ranges.end();
}

void varkit::SequenceRangeHandler::Add(varkit::SequenceRange &range) {
    m_ranges.emplace_back(range);
}

void varkit::SequenceRangeHandler::Add(varkit::SequenceRange &&range) {
    m_ranges.emplace_back(range);
}

/*
 * GENE Class
 */

varkit::SNP &varkit::Gene::GetSNP(varkit::SNPKey key) {
    return m_snps_map[key];
}

void varkit::Gene::AddSNP(const IO::IOSNP &iosnp) {
    constexpr bool debug = false;
//    SNPKey key = { iosnp.snp_position, iosnp.snp_base };
    SNPKey key = iosnp.snp_position;
//    if constexpr(debug) {
//        std::cout << iosnp.ToString() << " -> ";
//    }

//    auto& range = m_ranges.GetRange(iosnp.snp_position);
//    bool contains = std::any_of(range.GetReadInfo().begin(), range.GetReadInfo().end(), [&iosnp](ReadInfo const& ri) { return ri.read_id == iosnp.read_id; });
//    if (!contains) {
//        std::cerr << iosnp.ToString() << std::endl;
//        return;
//    } else {
//        std::cerr << "BONBONBON" << std::endl;
//    }

    auto& snp = GetSNP(key);
    snp.Add(iosnp);

//    if (iosnp.gene_id == 119 && iosnp.snp_position == 285) {
//        std::cout << snp.ToString() << "  ->  " << iosnp.ToString() << std::endl;
//    }
    if constexpr(debug) {
//        m_ranges.GetRange(iosnp.snp_position);
        std::cout << iosnp.ToString() << " -> ";
        std::cout << GetSNP(key).ToString() << std::endl;
    }

}

void varkit::Gene::SNPsToVector() const {
    constexpr bool debug = false;

    if constexpr(debug) {
        std::cout << "Print all snps" << std::endl;
    }

    for (const auto& [key, snp] : m_snps_map) {

        if constexpr(debug) {
            std::cout << snp.ToString() << std::endl;
        }
        m_snps.emplace_back(snp);
    }
    std::sort(m_snps.begin(), m_snps.end());

    if constexpr(debug) {
        std::cout << "Sorted SNPs" << std::endl;
        for (auto& snp : m_snps) {
            std::cout << snp.ToString() << std::endl;
        }
    }

    m_snps_map.clear();
}

void varkit::Gene::AddRange(const IO::ClassificationLine &line) {
    size_t start = std::max(line.genepos, 0);
    size_t end = start + line.overlap_with_gene;
//            std::cout << "Add " << start << ", " << end << std::endl;
//    m_ranges.Add(start, end);
    m_ranges.Add(line);
}

const varkit::SequenceRangeHandler &varkit::Gene::GetRangeHandler() const {
    return m_ranges;
}

varkit::SequenceRangeHandler varkit::Gene::GetRangeHandlerForCoverage(size_t minimum_coverage, size_t min_length) const {
    SequenceRangeHandler new_range_handler;

    size_t minimum_stretch_size = min_length;

    for (auto& range : m_ranges.GetRanges()) {
        range.SortReadInfo();
        auto& read_info = range.GetReadInfo();
        auto cov = range.CoverageVector();

        auto range_start = range.Start();
        auto sub_range_start = -1;
        auto sub_range_length = -1;
        size_t added = 0;
        for (auto i = 0; i < cov.size(); i++) {
            if (cov[i] >= minimum_coverage) {
                if (sub_range_start == -1) sub_range_start = i;
            } else {
                sub_range_length = i - sub_range_start;
                if (sub_range_start >= 0 && sub_range_length >= minimum_stretch_size) {
                    added++;
                    new_range_handler.Add(SequenceRange(sub_range_start + range_start, sub_range_start + range_start + sub_range_length));
                }
                sub_range_start = -1;
            }
        }
//        std::cout << "sub_range_start: " << sub_range_start << std::endl;
        sub_range_length = cov.size() - sub_range_start;
        if (sub_range_start >= 0 && sub_range_length > minimum_stretch_size) {
            added++;
            new_range_handler.Add(SequenceRange(sub_range_start + range_start, sub_range_start + range_start + sub_range_length));
        }
//
//            std::cout << "Original range: " << range.ToString() << std::endl;
//        std::cout << "______\nNew ranges: " << std::endl;
//        for (auto& range : new_range_handler.GetRanges()) {
//            std::cout << range.ToString() << std::endl;
//        }
//        if (added > 1) {
//            std::cout << m_ranges.ToString() << std::endl;
//            for (auto c : cov) std::cout << c << " ";
//            std::cout << std::endl;
//            Utils::Input();
//        }

    }
//    std::cout << m_ranges.ToString() << std::endl;
//    std::cout << new_range_handler.ToString() << std::endl;
//    Utils::Input();

    return new_range_handler;
}

const varkit::SNPList &varkit::Gene::GetSNPList() const {
    return m_snps;
}

varkit::SNPList &varkit::Gene::GetSNPList() {
    return m_snps;
}

varkit::SequenceRangeHandler &varkit::Gene::GetRangeHandler() {
    return m_ranges;
}

void varkit::Gene::SetGeneLength(size_t length) {
    m_gene_length = length;
}

size_t varkit::Gene::GetGeneLength() const {
    return m_gene_length;
}

const void varkit::Gene::CompareSNPsToReference(std::string const& reference) const {
    auto it = m_snps.begin();

    while (it != m_snps.end()) {
        SNP& snp = *it;

        size_t pos = snp.Position();
        char ref_base = reference.at(pos);
        snp.SetReference(ref_base);

        bool true_positive = true;
        for (auto& allele : snp.GetAlleles()) {
            if (allele.Base() == 'R') continue;

            true_positive &= allele.Base() != ref_base;
        }

//        std::cout << snp.ToString();
        if(!true_positive) {
//            std::cout << " remove";
            it = m_snps.erase(it);
        }
        else ++it;
//        std::cout << std::endl;
    }
}


/*
 * Strain Class
 */

varkit::Genes &varkit::Strain::GetGenes() {
    return m_genes;
}

const varkit::Genes &varkit::Strain::GetGenes() const {
    return m_genes;
}



/*
 * GENES Class
 */

bool varkit::Genes::HasGene(varkit::GeneID id) const {
    return m_genes.find(id) != m_genes.end();
}

varkit::GeneMap &varkit::Genes::GetGeneMap() {
    return m_genes;
}

const varkit::GeneMap &varkit::Genes::GetGeneMap() const {
    return m_genes;
}

varkit::Gene &varkit::Genes::GetOrAddGene(varkit::GeneID id) noexcept {
    return m_genes[id];
}

void varkit::Genes::AddGene(varkit::GeneID id, size_t gene_length) {
    m_genes.insert( { id, Gene(gene_length) } );
}

varkit::Gene &varkit::Genes::GetGene(varkit::GeneID id) {
    auto it = m_genes.find(id);
    if (it == m_genes.end()) {
        errx(EX_IOERR, "GetGene() in strain_handler.h - Gene id (%lu) not present", id);
    }
//            return it->second;
    return m_genes[id];
}

const varkit::Gene &varkit::Genes::GetGene(varkit::GeneID id) const {
    auto it = m_genes.find(id);
    if (it == m_genes.end()) {
        errx(EX_IOERR, "GetGene() in strain_handler.h - Gene id (%lu) not present", id);
    }
//            return it->second;
    return it->second;
}
