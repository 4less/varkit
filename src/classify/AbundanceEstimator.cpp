//
// Created by fritsche on 07/12/2021.
//

#include <iostream>
#include "AbundanceEstimator.h"
#include <csv.h>
#include <cmath>
#include <progress_bar.h>

void AbundanceEstimator::Clear() {
    this->taxon2depth.clear();
    this->taxa.clear();
}

void AbundanceEstimator::ReadClassificationFile(const char* file_path) {
    std::ifstream is(file_path, std::ios::binary);

    // Debug and remove later
//    m_valid_readids.clear();
    m_valid_readuids.clear();

    auto filesize = Utils::GetFileSize(file_path);

    ProgressBar bar;
    bar.reset(filesize);

    size_t progress = 0;
    size_t progress_chunk_size = filesize / 100 + 1;
    size_t progress_chunk_count = 0;

    bar.Update(progress);

    bool has_truth = false;

    // For
    size_t special_counter = 0;

    size_t counter = 0llu;
    while (!is.eof()) {
        IO::ClassificationLine line;
        is.read((char*) &line, sizeof(IO::ClassificationLine));
        if (is.eof()) break;

//        //4100	71	1691
//        if (line.taxid == 4100 && line.geneid == 71 && line.genepos >= (1691 - 150) && line.genepos < (1691 + 150)) {
//            special_counter++;
////            std::cout << line.ToString() << std::endl;
//        }

        if (line.geneid > 120) {
            std::cout << "o.g. whoops: " << line.geneid << std::endl;
            std::cout << line.ToString() << std::endl;
            exit(9);
        }
        counter++;

//        std::cout << line.ToString() << std::endl;

        if (has_truth) {
            std::string header{line.buffer};
            auto true_id = stol(Utils::split(header, '_'));
        }

        AddLine(line);

        progress += sizeof(IO::ClassificationLine);

        if (progress > (progress_chunk_count + 1) * progress_chunk_size) {
            while (progress > (progress_chunk_count + 1) * progress_chunk_size)
                progress_chunk_count++;
            bar.Update(progress);
        }
    }

//    std::cout << "special_counter: " << special_counter << " " << m_active_sample->GetStrain(4100).GetGenes().GetGene(71).GetRangeHandler().GetRange(1691).Coverage(1691) << std::endl;

    if (has_truth)
        bce.WriteStats();
    is.close();
}

AbundanceEstimator::AbundanceEstimator(size_t k, size_t shape_length, ReadFilter &read_filter, TaxonFilter &taxon_filter, Taxonomy::IntTaxonomy &taxonomy) :
        k(k),
        shape_length(shape_length),
        read_filter(read_filter),
        taxon_filter(taxon_filter),
        taxonomy(taxonomy) {

}

bool AbundanceEstimator::AddLine(IO::ClassificationLine &line, bool force_pass, int true_internal_id) {
    // Calculate overlap of read with gene
    line.overlap_with_gene = line.read_length;

    int correct_classification = 0;
    int incorrect_classification = 0;

    if (taxonomy.HasNode(true_internal_id)) {
        correct_classification = true_internal_id == line.taxid;
        incorrect_classification = !taxonomy.IsNodeAncestor(line.taxid, true_internal_id);
    }

    if (line.geneid > 0) {
        if (line.genepos < 0) {
            // G: Gene, R: Read  genepos < 0
            //     RRRRRRRRRRRR
            // -------GGGGGGGGGGGGG------------------ Genome
            //        rrrrrrrrr <-- overlap
            line.overlap_with_gene = line.read_length + line.genepos;
        } else if (line.genepos + line.read_length > m_marker_genomes[line.taxid].Get(line.geneid)) {
            // G: Gene, R: Read  genepos + readlen > gene_len
            //               RRRRRRRRRRRR
            // -------GGGGGGGGGGGGG------------------ Genome
            //               rrrrrr <-- overlap
            line.overlap_with_gene = line.read_length - (line.genepos + line.read_length - m_marker_genomes[line.taxid].Get(line.geneid));

        }
        // Debug
        if (line.overlap_with_gene < shape_length) return false;


        // Predict ani not by snps but by hits
        if ((line.overlap_with_gene - shape_length + 1) == 0) {
            std::cerr << "Blurp: " << line.ToString() << std::endl;
        }
        if (m_estimate_ani)
            line.predicted_ani = AbundanceEstimator::EstimateAni(line.total_hits_best, line.overlap_with_gene, shape_length, k);


        if (line.predicted_ani < 0 || line.predicted_ani > 1) {
            std::cerr << "ANI:" << line.predicted_ani << '\t' << line.ToString() << std::endl;
            std::cerr << "MGLEN: " << m_marker_genomes[line.taxid].Get(line.geneid) << std::endl;
            Utils::Input();
        }
    }



    // Filter reads if hit is not good enough
    auto pass = read_filter.Pass(line, force_pass);
//
//    if (line.record_id == 617014 || line.record_id == 866054 ||
//            line.record_id == 505041 || line.record_id == 655207 || line.record_id == 526329 || line.record_id == 487458 ||
//            line.record_id == 636123 || line.record_id == 142820 || line.record_id == 733553 || line.record_id == 891181) {
//        std::cout << pass << " " << line.ToString() << std::endl;
//        Utils::Input();
//    }


    if (pass) {
        // PROCESS PASSED READ FOR DISTANCE MATRIX TO SAMPLE HANDLER
//        std::cout << "line.HasGene() " << line.HasGene() << " -> " << line.geneid << std::endl;
//        std::cout << "Go in if: " << (taxonomy.IsLeaf(line.taxid) && line.HasGene()) << std::endl;

//        std::cout << "Load into sample handler" << std::endl;
        if (taxonomy.IsLeaf(line.taxid) && line.HasGene()) {
            auto &genes = m_active_sample->GetOrAddStrain(line.taxid).GetGenes();
            if (!genes.HasGene(line.geneid)) {
                auto marker_genome_it = m_marker_genomes.find(line.taxid);
                auto length = marker_genome_it->second.Get(line.geneid);
                genes.AddGene(line.geneid, marker_genome_it->second.Get(line.geneid));
            }
            if (line.genepos < 0) line.genepos = 0;

            size_t old = genes.GetGene(line.geneid).GetRangeHandler().HasRange(line.genepos) ?
                    genes.GetGene(line.geneid).GetRangeHandler().GetRange(line.genepos).Coverage(line.genepos) : 0;
            genes.GetGene(line.geneid).AddRange(line);

            // Debug eventuelly remove
//            m_valid_readids.insert(line.record_id);
            m_valid_readuids.insert( { line.record_id, line.read_num } );



            if (!genes.GetGene(line.geneid).GetRangeHandler().HasRange(line.genepos) || old+1 != genes.GetGene(line.geneid).GetRangeHandler().GetRange(line.genepos).Coverage(line.genepos)) {
                std::cout << "Pos: " << line.genepos << std::endl;
                std::cout << "Old: " << old << " -> " << genes.GetGene(line.geneid).GetRangeHandler().GetRange(line.genepos).Coverage(line.genepos) + 1 << std::endl;
                std::cout << genes.GetGene(line.geneid).GetRangeHandler().GetRange(line.genepos).ToVerboseString() << std::endl;
                auto& range = genes.GetGene(line.geneid).GetRangeHandler().GetRange(line.genepos);
                std::cout << range.Coverage(line.genepos) << std::endl;
                for (auto& rinfo : range.GetReadInfo()) {
                    std::cout << rinfo.ToString() << " " << rinfo.Contains(line.genepos) << std::endl;
                }
                exit(81);
            }
            //           genes.GetOrAddGene(line.geneid).AddRange(line);
        } else {
//            std::cout << "Not registered" << std::endl;
        }

        if (!taxa.contains(line.taxid)) {
            // If hit not in taxa create new taxon to record hits
            // If hit is on leaf level (to genome) load Information on marker gene
            taxa.insert(
                    {
                        line.taxid,
                        taxonomy.IsLeaf(line.taxid) ? Taxon(line.taxid, m_marker_genomes[line.taxid]) : Taxon(line.taxid)
                    });
        }
        // Add line to taxon
        taxa[line.taxid].Add(line);

        if (true_internal_id > 0) {
            bce.tp += correct_classification;
            bce.fp += incorrect_classification;
        }

    } else {
        // Record hits that are not true
        // Install mechanism to differentiate true positive from true negative if reference is available
//        std::cout << "Not pass" << std::endl;

        if (true_internal_id > 0) {
            bce.tn += incorrect_classification;
            bce.fn += correct_classification;
        }
    }
    return pass;
}

void AbundanceEstimator::LoadMarkerGenome(const char *file_path) {

    auto str = csv::load_file(file_path);
    auto parser = csv::make_parser( str , '\t');

    for (auto&& row : parser ) {
        auto it = row.begin();

        size_t taxid = (*it).to_long();
        size_t geneid = (*(++it)).to_long();
        int mg_length = (*(++it)).to_int();

        if (!m_marker_genomes.contains(taxid)) {
//            std::cout << "Insert: " << taxid << ", " << geneid << ", " << mg_length << std::endl;
            m_marker_genomes.insert({taxid, MarkerGenome(120)});
        }
        m_marker_genomes[taxid].Set(geneid, mg_length);
    }

    for (auto& mg : m_marker_genomes) {
        mg.second.CalculateMGSize();
    }
}

void AbundanceEstimator::LoadAbundanceCorrection(const char *file_path) {

    auto str = csv::load_file(file_path);
    auto parser = csv::make_parser( str , '\t');

    for (auto&& row : parser ) {
        auto it = row.begin();

        int taxid = (*it).to_int();
        int geneid = (*(++it)).to_int();
        double factor = (*(++it)).to_double();

        if (!taxon2abundance_factor.contains(taxid)) {
            std::unordered_map<size_t, double> inner;
            taxon2abundance_factor.insert( { (size_t)taxid, inner });
        }
        taxon2abundance_factor[taxid].insert( std::make_pair( geneid, factor) );
    }
    correct_abundances = true;
}

ReadFilter &AbundanceEstimator::GetReadFilter() {
    return read_filter;
}

std::string AbundanceEstimator::GetBioboxHeader(std::string sample_name) const {
    return "@SampleID: " + sample_name +
            "\n@Version: 0.9.1" +
            "\n@Ranks:superkingdom|phylum|class|order|family|genus|species|strain" +
            "\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n";
}

std::string AbundanceEstimator::BioboxString(size_t id, double abundance, bool as_leaf, std::vector<std::string> ranks) const {
    //id, rank, idstring, namestring, abundance
    auto& node = taxonomy.Get(id);

    std::string rank = node.rank;
    std::string rep_genome = "";
    auto external_str = std::to_string(node.external_id);

    if (node.IsLeaf()) {
        if (as_leaf) {
            rank = "strain";
            external_str += ".1";
        }
        if (!node.rep_genome.empty())
            rep_genome = node.rep_genome;
        else {
            rep_genome = node.scientific_name;
        }
    }

    auto lin = taxonomy.LineageStr(id, ranks, "|");
    auto external_lin = taxonomy.LineageExternalIds(id, ranks, "|");
    auto name = node.scientific_name;

    // HACKY HACKY HACKY

    // TODO make unhacky
    lin = "Bacteria" + Utils::trim_end(lin, '|');
    external_lin = "2" + Utils::trim_end(external_lin, '|');


    return external_str + '\t' + rank + '\t' + external_lin + '\t' + lin + '\t' + std::to_string((abundance*100)) + '\t' + rep_genome;
}


void AbundanceEstimator::SetVerbose(bool verbose) {
    this->verbose = verbose;
}

void AbundanceEstimator::Profile(std::string output_file) {
//    double expected_abundance = 216.0;
    double expected_abundance = 0.862;
    double se = 0;
    size_t counter = 0;

    bool tobiom = true;

    if (verbose) {
        verbose_output = new std::ofstream(output_file + ".log", std::ios::out);
    }

    for (auto& pair : taxa) {
        auto &id = pair.first;

        auto &taxon = pair.second;
        if (taxon.leaf) {
            if (raw_taxa) {
                *raw_taxa << taxon.VerboseString() << std::endl;
            }
//            std::cout << pair.first << " " << pair.second.VerboseString() << std::endl;


            auto mean_ani = Statistics::Mean(taxon.anis);
            auto median_ani = Statistics::Median(taxon.anis);

            if (mean_ani < 0 || median_ani < 0 || mean_ani > 1 || median_ani > 1) {
                std::cerr << mean_ani << "," << median_ani << "\t" << taxon.VerboseString() << std::endl;
                std::cerr << "ANIS: ";
                for (auto& ani : taxon.anis) {
                    std::cout << ani << ",";
                }
                std::cout << "\n" << std::endl;
            }

            double presence = ((double) taxon.PresentGenes() / taxon.marker_genome->MarkerGeneCount());
            double expected_presence = ((double) taxon_filter.ExpectedGenes(taxon.MarkerGeneCount(), taxon.mapped_reads)) /
                                       taxon.MarkerGeneCount();

            auto gene_abundances = correct_abundances ?
                    taxon.GetBucketAbundances3(taxon2abundance_factor[id], 0) :
                    taxon.GetBucketAbundances2();

            auto median = Statistics::Median(gene_abundances);

            if (gene_abundances.empty()) {
                if (taxon.mapped_reads > 5)
                    std::cerr << taxon.VerboseString() << "   abundance: " << median << std::endl;
            }
//            auto abundance = median;
            auto abundance = median;

            se += pow(expected_abundance - abundance, 2);
            counter++;


//            if (abundance == 0)
//                std::cout << "abundance: " << abundance << " " << taxonomy.LineageStr(id) << std::endl;
            taxon2depth.insert( { id, abundance } );
        }
    }



    total_depth = 0;
    size_t leaf_depth = 0;
    for (auto t : taxon2depth) {
        auto taxon = taxonomy.Get(t.first);
        if (taxon.IsLeaf())
            leaf_depth += t.second;
        else {
            std::cout << taxon.rank << " " << t.second << std::endl;
        }

        total_depth += t.second;
    }
//    std::cout << "total_depth: " << total_depth << std::endl;
//    std::cout << "leaf_depth: " << leaf_depth << std::endl;
    total_depth = leaf_depth;


    // Filter taxa and output into requested formats
    size_t passed_taxa = 0;
    size_t total_taxa = 0;

//    std::cout << "Output profile... " << output_file << std::endl;
//    std::ofstream of(output_file, std::ios::out);


//    std::unordered_map<size_t, double> abundances;
//    std::unordered_map<std::string, double> total_depths;
    m_abundances.clear();
    m_total_depths.clear();
    m_predicted_strains.clear();

    // Calculate the total sequencing depth for each rank by summing up
    // depths for all taxa that passed the filter.
    for (auto& entry : taxa) {
        auto internal_id = entry.first;
        auto &taxon = entry.second;
        auto node = taxonomy.Get(internal_id);
        if (taxon_filter.Pass(taxon)) {
            auto sequencing_depth = taxon2depth[internal_id];

            m_predicted_strains.insert(taxon.id);

            // Add to predicted strains in sample handler
            m_active_sample->AddPredictedStrain(taxon.id);

            if (node.IsLeaf()) {
                auto rank = "strain";
                if (!m_total_depths.contains(rank)) {
                    m_total_depths.insert( { rank, 0.0 } );
                }
                m_total_depths[rank] += sequencing_depth;
            }
            while (internal_id != 1) {
                auto node = taxonomy.Get(internal_id);
                if (!m_total_depths.contains(node.rank)) {
                    m_total_depths.insert( { node.rank, 0.0 } );
                }
                m_total_depths[node.rank] += sequencing_depth;
                internal_id = taxonomy.Get(internal_id).parent_id;
            }
        }
    }

//    std::cout << "Print total per-rank sequencing depths" << std::endl;
//    for (auto &entry : total_depths) {
//        std::cout << entry.first << " " << entry.second << std::endl;
//    }

//    std::cout << "ITERATE ALL TAXA AND FILTER " << taxa.size() << std::endl;
    std::vector<std::string> lineage_ranks={ "superkingdom", "phylum", "class", "order", "family", "genus", "species" };
    for (auto& entry : taxa) {
        auto internal_id = entry.first;
        auto& taxon = entry.second;

        if (!taxon.leaf) {
            continue;
        }
        if (!taxon2depth.contains((uint64_t)internal_id)) {
            std::cout << internal_id << " not in taxon2depth" << std::endl;
        }

        auto sequencing_depth = taxon2depth[internal_id];

        if (verbose && taxon.leaf) {
            auto external_lin = taxonomy.LineageExternalIds(taxon.id, lineage_ranks, "|");
            auto expected_genes = TaxonFilter::ExpectedGenes(taxon);
//            std::cout << "XXX" << taxon.id << " " << taxon.mapped_reads << " " << expected_genes << std::endl;
            *verbose_output << taxonomy.Get(taxon.id).external_id << '\t' << expected_genes << '\t' << external_lin << '\t' << taxon.VerboseString() << '\n';
        }

        total_taxa++;
        if (taxon_filter.Pass(taxon)) {
            passed_taxa++;

            while (internal_id != 1) {
                auto node = taxonomy.Get(internal_id);
                auto rank_depth = m_total_depths[node.rank];

                auto abundance = (double) sequencing_depth/rank_depth;

                if (!m_abundances.contains(internal_id))
                    m_abundances.insert( { internal_id, 0.0 } );

                m_abundances[internal_id] += abundance;

                internal_id = taxonomy.Get(internal_id).parent_id;
            }

        } else {

        }
    }


    if (verbose) {
        verbose_output->close();
        delete verbose_output;
    }

    // OUTPUT BIOM FORMAT

//    bool output_strain = true;
//    of << GetBioboxHeader(sample_name) << std::endl;
//    if (output_strain) {
//        for (auto &entry : m_abundances) {
//            auto internal_id = entry.first;
//            auto abundance = entry.second;
//            auto taxon = taxonomy.Get(internal_id);
//
//            auto rank = taxon.rank;
//
//            if (taxon.IsLeaf()) {
//                auto taxon2 = taxa[internal_id];
////                std::cout << taxon.scientific_name << " ani: " << Statistics::Mean(taxon2.anis) << " " << Statistics::Median(taxon2.anis);
////                std::cout << " " << BioboxString(internal_id, abundance, true) << std::endl;
//                of << BioboxString(internal_id, abundance, true) << std::endl;
//            }
//        }
//    }
//
//    auto ranks = { "species", "genus", "family", "order", "class", "phylum", "root" };
////    auto ranks = { "phylum" };
//    for (auto rank : ranks) {
////        std::cout << "RANK:  " << rank << std::endl;
//
//        for (auto& entry : m_abundances) {
//
//            auto internal_id = entry.first;
//            auto abundance = entry.second;
//            auto taxon = taxonomy.Get(internal_id);
//
//
//            if (taxon.rank == rank) {
////                std::cout << abundance << std::endl;
////                std::cout << BioboxString(internal_id, abundance) << std::endl;
//                of << BioboxString(internal_id, abundance) << std::endl;
//            }
//        }
//    }
//    of.close();

    WriteBiobox(output_file);



    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Report:" << std::endl;
    std::cout << "passed rate: " << double(passed_taxa)/double(total_taxa) << std::endl;
    std::cout << "total passed: " << passed_taxa << std::endl;




//    std::cout << "------" << std::endl;
//    auto hit = taxon2depth[10233];
//    std::cout << BioboxString(10233, (double)hit/total_depth) << std::endl;


//    std::cout << "total: " << total_depth << std::endl;
    std::cout << std::flush;
//    std::cerr << "Mean squared error is: " << (se/counter) << std::endl;
}

void AbundanceEstimator::SetSample(std::string basicString) {
    this->sample_name = basicString;
    this->sample_idx = m_current_sample_idx++;
    m_strain_handler.AddSample(basicString);
    m_active_sample = &m_strain_handler.GetSample(basicString);
}

void AbundanceEstimator::SetSample(Sample& sample) {
    this->sample_name = sample.GetBasenameNoExtension();
    this->sample_idx = sample.GetIndex();
    m_strain_handler.AddSample(sample);
    m_active_sample = &m_strain_handler.GetSample(sample.GetBasenameNoExtension());
}

AbundanceEstimator::~AbundanceEstimator() {
    delete raw_taxa;
    delete verbose_output;
}

void AbundanceEstimator::SetRawOutput(std::string path) {
    std::cout << "RAW OUTPUT: " << path << std::endl;
    raw_taxa = new std::ofstream(path, std::ios::out);
}


// DEBUG FUNCTION
void AbundanceEstimator::DebugFunction() {
    auto& sample = *m_active_sample;


//    for (auto pair1 : m_sample.GetStrains()) {
//        for (auto pair2 : m_sample.GetStrains()) {
//
//        }
//    }

    for (auto pair : sample.GetStrains()) {

        auto& strain = pair.second;
        auto& genes = strain.GetGenes();

        for (auto& gene_pair : genes.GetGeneMap()) {
            std::cout << "GeneId: " << gene_pair.first;
            std::cout << "  length: " << gene_pair.second.GetGeneLength();
            std::cout << "    Range: " << gene_pair.second.GetRangeHandler().ToString() << std::endl;
        }
        std::cout << "Strain: " << pair.first << std::endl;
        std::cout << std::string(79, '-') << std::endl;


//
//        for (auto& gene_pair1 : genes.GetGeneMap()) {
//            auto& handler1 = gene_pair1.second.GetRangeHandler();
//            for (auto& gene_pair2 : genes.GetGeneMap()) {
//                auto& handler2 = gene_pair2.second.GetRangeHandler();
//
//                if (&handler1 == &handler2) continue;
////                std::cout << "GeneId: " << gene_pair1.first;
////                std::cout << "  length: " << gene_pair1.second.GetGeneLength();
////                std::cout << "    Range: " << gene_pair1.second.GetRangeHandler().ToString() << std::endl;
//                std::cout << "Compare " << gene_pair1.first << " with " << gene_pair2.first << std::endl;
//
//                std::cout << "Intersect: " << gene_pair1.first << " with " << gene_pair2.first << std::endl;
//                std::cout << "Range1: " << handler1.ToString() << std::endl;
//                std::cout << "Range2: " << handler2.ToString() << std::endl;
//
//                auto intersection = handler1.Intersect(handler2);
//                std::cout << "RangeI: " << intersection.ToString() << std::endl;
//                if (intersection.Size() > 1)
//                    Utils::Input();
//            }
//        }



    }
}

void AbundanceEstimator::WriteBiobox(const std::string path) const {
    // OUTPUT BIOM FORMAT

    ofstream of(path, std::ios::out);
    bool output_strain = true;
    of << GetBioboxHeader(sample_name) << std::endl;
    if (output_strain) {
        for (auto &entry : m_abundances) {
            auto internal_id = entry.first;
            auto abundance = entry.second;
            auto taxon = taxonomy.Get(internal_id);

            auto rank = taxon.rank;

            if (taxon.IsLeaf()) {

                auto taxon2 = taxa.at(internal_id);
//                std::cout << taxon.scientific_name << " ani: " << Statistics::Mean(taxon2.anis) << " " << Statistics::Median(taxon2.anis);
//                std::cout << " " << BioboxString(internal_id, abundance, true) << std::endl;
                of << BioboxString(internal_id, abundance, true) << std::endl;
            }
        }
    }

    auto ranks = { "species", "genus", "family", "order", "class", "phylum", "root" };
//    auto ranks = { "phylum" };
    for (auto rank : ranks) {
//        std::cout << "RANK:  " << rank << std::endl;

        for (auto& entry : m_abundances) {

            auto internal_id = entry.first;
            auto abundance = entry.second;
            auto taxon = taxonomy.Get(internal_id);


            if (taxon.rank == rank) {
//                std::cout << abundance << std::endl;
//                std::cout << BioboxString(internal_id, abundance) << std::endl;
                of << BioboxString(internal_id, abundance) << std::endl;
            }
        }
    }
    of.close();
}


void Taxon::Add(IO::ClassificationLine &line) {
    // Add taxon to line
    mapped_reads += 1;
    mapped_length += line.overlap_with_gene;

    hits.push_back(line.total_hits_best);
    overlap.push_back(line.overlap_with_gene);

    if (line.predicted_ani != -1)
        anis.push_back(line.predicted_ani);

    if (line.predicted_ani != -1 && (line.predicted_ani < 0 && line.predicted_ani > 1)) {
        std::cerr << "\nANI:" << line.predicted_ani << '\t' << line.ToString() << std::endl;
        Utils::Input();
    }

    if (line.geneid > 0) {
        if (line.geneid >= marker_gene_buckets.size()) {
            std::cout << "whoops geneid: " << line.geneid << " >= " << marker_gene_buckets.size() << " taxid: " << line.taxid << std::endl;
            exit(123);
        }
        marker_gene_buckets[line.geneid].Add(line);
    }
}

size_t Taxon::PresentGenes() {
//    std::cout << "present genes: " << std::count_if(marker_gene_buckets.begin(), marker_gene_buckets.end(), [](int i){return i > 0;}) << std::endl;
    return std::count_if(marker_gene_buckets.begin(), marker_gene_buckets.end(), [](Gene gene){return gene.mapped_bases > 0;});
}

size_t Taxon::MarkerGeneCount() {
    if (!marker_genome) {
        std::cerr << this->id << " has no marker genome" << std::endl;
        exit(8);
    }
    return marker_genome->MarkerGeneCount();
}

double Taxon::NaiveVerticalCoverage() {
//    std::cout << "\nMapped SequenceLength: " << mapped_length << " MGLEN: " << marker_genome.MarkerGenomeLength() << " mapped reads: " << mapped_reads << std::endl;
    return (double)mapped_length/marker_genome->MarkerGenomeLength();
}

std::vector<double> Taxon::GetBucketAbundances() {
    std::vector<double> result;
    for (auto& gene : marker_gene_buckets) {
        if (!gene.HasGene()) continue;
        if (!gene.IsPresent()) continue;

        //debug
        auto vc = gene.VerticalCoverage();
        if (isnanl(vc)) {
            std::cout << "isnan: " << gene.TotalBuckets() << " " << gene.mapped_bases << " " << gene.total_reads << std::endl;
        }

        result.push_back(gene.VerticalCoverage());
    }
    return result;
}

std::vector<double> Taxon::GetBucketAbundances2() {
    std::vector<double> result;
    for (auto& gene : marker_gene_buckets) {
        if (!gene.HasGene()) continue;
        if (!gene.IsPresent()) continue;
//        result.push_back(gene.VerticalCoverage2());
        result.push_back(gene.VerticalCoverage());
    }
    return result;
}

double Taxon::VerticalCoverage() {
    if (m_vertical_coverage == -1) {
        m_vertical_coverage = Statistics::Median(GetBucketAbundances2());
    }
    return m_vertical_coverage;
//    return 0.0;
}

std::vector<double> Taxon::GetBucketAbundances3(std::unordered_map<size_t, double> &map, double threshold) {
    std::vector<double> result;
    for (auto& gene : marker_gene_buckets) {
        if (!gene.HasGene()) continue;
        if (!gene.IsPresent()) continue;
        if (!map[gene.geneid]) continue;
        if ((map[gene.geneid]) == 0) {
            std::cerr << "map[gene.geneid] is 0" << std::endl;
            exit(51);
        }
        if (map[gene.geneid] < threshold) continue;
        auto corrected_abundance = gene.VerticalCoverage() * (1/map[gene.geneid]);
        result.push_back(corrected_abundance);
    }
    if (result.empty()) result.push_back(0);
    return result;
}

std::string Taxon::VerboseString() {
    std::string return_str = "";


    auto mean_ani = Statistics::Mean(anis);
    auto median_ani = Statistics::Median(anis);

    return_str += std::to_string(this->id);
    return_str += '\t';
    return_str += std::to_string(this->VerticalCoverage());
    return_str += '\t';
    return_str += std::to_string(mean_ani);
    return_str += '\t';
    return_str += std::to_string(median_ani);
    return_str += '\t';
    return_str += std::to_string(this->mapped_length);
    return_str += '\t';
    return_str += std::to_string(this->mapped_reads);
    return_str += '\t';
    return_str += std::to_string(this->PresentGenes());
    return_str += '\t';
    return_str += std::to_string(this->MarkerGeneCount());
    return_str += '\t';
    return_str += std::to_string((double) this->PresentGenes()/this->MarkerGeneCount());
    return_str += '\t';

    std::for_each(this->marker_gene_buckets.begin(), this->marker_gene_buckets.end(), [&return_str](Gene &gene){return_str += std::to_string(gene.total_reads) + ',';});

    return return_str;
}


Gene::Gene(MarkerGenome &marker_genome, size_t geneid) : geneid(geneid) {
    gene_length = marker_genome.Get(geneid);
    size_t min_size = 66;
    size_t target_size = 100;
    size_t max_size = 133;

    size_t bucket_count_upper = marker_genome.Get(geneid)/target_size + 1;
    size_t bucket_count_lower = marker_genome.Get(geneid)/target_size;

    if (bucket_count_lower == 0) {
        bucket_count = bucket_count_upper;
    } else if (target_size - (gene_length/bucket_count_upper) < (gene_length/bucket_count_lower) - target_size) {
        bucket_count = bucket_count_upper;
    } else
        bucket_count = bucket_count_lower;
    if ((bucket_count) == 0) {
        std::cerr << "bucket_count is 0" << std::endl;
        exit(51);
    }
//    std::cout << "gene_length: " << gene_length << std::endl;
//    std::cout << "bucket_count: " << bucket_count << std::endl;
    bucket_length = gene_length / bucket_count;
    bucket_count++; // check if correct
    if (bucket_length > 10000) {
        std::cout << bucket_length << std::endl;
        exit(8);
    } else if (bucket_length == 0) {
        std::cout << "Gene length: " << gene_length << std::endl;
        std::cout << "bucket_count: " << bucket_count << std::endl;
        std::cout << "geneid: " << geneid << std::endl;
        exit(2);
    }

    if (bucket_length == 0) {
        std::cerr << "gene_length: " << gene_length << std::endl;
        std::cerr << "bucket_count: " << bucket_count << std::endl;
        exit(51);
    }

    if ((bucket_count * bucket_length) < gene_length) {
        std::cout << "genelength: " << gene_length << std::endl;
        std::cout << "bucket_count * bucket_length: " << (bucket_length * bucket_count) << std::endl;
    }
//    std::cout << gene_length/bucket_count_upper << " " << (gene_length/bucket_count_lower) << std::endl;
//    std::cout << "bucket_count_upper " << bucket_count_upper << " bucket_count_lower: " << bucket_count_lower << std::endl;
//    std::cout << "bucket_count " << bucket_count << " bucket_length: " << bucket_length << std::endl;
    buckets.resize(bucket_count);
}

double Gene::ExpectedActiveBuckets() {
    return (1 - pow((bucket_count - 1)/bucket_count, total_reads)) * bucket_count * correction_factor;
}

double Gene::VerticalCoverage() {
    if (!IsPresent()) return 0;
    if (ExpectedActiveBuckets() * bucket_length == 0) {
        std::cerr << "(ExpectedActiveBuckets() * bucket_length is 0" << std::endl;
        exit(51);
    }
    return (double) mapped_bases / (ExpectedActiveBuckets() * bucket_length);
}



double Gene::VerticalCoverage2() {
    if (!IsPresent()) return 0;
    if (ActiveBuckets() * bucket_length == 0) {
        std::cerr << "(ActiveBuckets() * bucket_length is 0" << std::endl;
        exit(51);
    }
    return (double) mapped_bases / (ActiveBuckets() * bucket_length);
}


size_t Gene::ActiveBuckets() {
    if (buckets.empty()) {
        exit(123);
    }
    size_t active_buckets = 0;
    for (size_t i = 0; i < bucket_count; i++) {
        if (!buckets.empty() && buckets[i] > 0) {
            active_buckets++;
        }
    }
    return active_buckets;
}

size_t Gene::TotalBuckets() {
    return bucket_count;
}

double Gene::ActiveRate() {
    return (double) ActiveBuckets()/TotalBuckets();
}

void Gene::Add(IO::ClassificationLine &line) {
    auto pos = line.genepos < 0 ? 0 : line.genepos;
//    if (pos >= bucket_count) {
//        std::cout << line.genepos << std::endl;
//        std::cout << bucket_length << std::endl;
//        std::cout << pos << " >= " << bucket_count << std::endl;
//        exit(23);
//    }
    total_reads++;
    if (bucket_length == 0) {
        std::cerr << line.ToString() << std::endl;
        std::cerr << "bucket length is 0" << std::endl;
        exit(51);
    }

    if ((pos / bucket_length) >= bucket_count || bucket_length > 10000) {
        std::cout << "abort" << std::endl;
        std::cout << "pos: " << pos << std::endl;
        std::cout << "bucket_count: " << bucket_count << std::endl;
        std::cout << "bucket_length: " << bucket_length << std::endl;
        std::cout << "gene_length: " << gene_length << std::endl;
        std::cout << "(pos / bucket_length): " << (pos / bucket_length) << std::endl;
        std::cout << bucket_length << std::endl;
        exit(34);
    }

    buckets[pos / bucket_length]++;

    mapped_bases += line.overlap_with_gene;
}

bool Gene::HasGene() {
    return !buckets.empty();
}

bool Gene::IsPresent() {
    return total_reads > 0;
}

std::vector<int>& MarkerGenome::MarkerGeneLengths() {
    return marker_gene_lengths;
}
