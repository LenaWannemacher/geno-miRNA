#include "json.hpp"

#include "Candidate.h"
#include "miRNA.h"
#include "algorithm.h"
#include "output.h"
#include "intervallTree.h"
#include "GenomeProcessor.h"
#include "GFFReader.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <direct.h>
#include <map>
#include <vector>
#include <algorithm>


using namespace std;
using json = nlohmann::json;


json readParameters(const char* path)
{
    std::ifstream istrm(path);

    if(!istrm) {
        std::cerr << "Could not open file " << path << " for reading.\n";

        // Call exit here to avoid complicated error handling. Only use this
        // method in the main file as it makes program abortion unrecoverable
        // and hard to predict.
        std::exit(-2);
    }

    try {
        json parameters;
        istrm >> parameters;

        return parameters;
    } catch(const std::invalid_argument& e) {
        std::cerr << "Error while parsing config file: " << e.what() << '\n';
    }

    // See above
    std::exit(-3);
}


int main(int argc, char* argv[]) {

    string option = argv[1];

    if(option == "-h" || option == "--h"){

        cout << "Usage:\n\t"
             << "--align   miRNA.fa   genome.fa   genome_annotation.gff3   config.json   tmp_dir   output_dir\nor\n\t"
             << "further instructions can be found in the README file\n";
        return 0;

    }

    if (argc < 2){
        cerr << "Not enough arguments!\nUsage:\n\t"
             << "--align   miRNA.fa   genome.fa   genome_annotation.gff3   config.json   tmp_dir   output_dir\nor\n\t"
             << "type in   --help or -h   for detailed help\n";
        return -1;
    }

    /**
     * align miRNA with genome
     */
    if(option == "--align"){

        /// 1. Processing the miRNA file
        double duration_overall;
        std::clock_t start_clock_overall = std::clock();

        auto parameters = readParameters(argv[5]);

        size_t match = parameters["match"];
        size_t gu = parameters["g_u_wobble"];
        int penalty = parameters["penalty"];
        int gap_open = parameters["gap_open_pen"];
        int gap_ext = parameters["gap_ext_pen"];
        float scal_f = parameters["scaling_factor"];
        size_t threshold = parameters["threshold"];
        size_t token_size = parameters["token_size"];
        size_t overhang = parameters["overhang"];
        size_t wobble_region = parameters["wobble_region"];
        size_t output_alignment = parameters["output_alignment"];
        size_t output_stats = parameters["output_statistics"];

        const char* tmp_dir = argv[6];
        const char* output = argv[7];

        std::string tmp_dir_genome_path = std::string(tmp_dir).append("\\temp_genome.fa");
        std::string tmp_dir_candidate_path = std::string(tmp_dir).append("\\all_candidates.fa");

        std::string output_table_path = std::string(output);
        std::string output_alignment_path = std::string(output);

        fstream mirna_file;
        mirna_file.open(argv[2]);

        std::string stream_line;
        std::vector<miRNA> miRNA_list;

        std::string mirna_seq, mirna_id;

        if(!mirna_file.is_open()){
            cerr << "Could not open the given miRNA file!" << endl;
            exit(EXIT_FAILURE);
        }

        cout << "Parsing the given miRNA-file..." << endl;

        while(getline(mirna_file, stream_line)){

            if(stream_line.substr(0,1) == ">"){
                mirna_id = stream_line.substr(1,stream_line.find(' '));
            }
            else{
                mirna_seq = stream_line;

                miRNA mirna(mirna_seq, mirna_id);
                mirna.generate_complementary_RNA(mirna.get_reversed_miRNA());
                mirna.convert_to_DNA(mirna.get_changed_miRNASequence());
                mirna.set_miRNAToken(token_size);
                mirna.set_miRNAOverhang(overhang);

                miRNA_list.push_back(mirna);
            }
        }

        cout << "--- Parsing complete! ---\n-------------------------\n" << endl;

        /// 2. Processing the genome file
        GenomeProcessor genpro(argv[3]);

        cout << "\nNow processing genome...\n";

        genpro.setTempGenFilePath(tmp_dir_genome_path.c_str());
        genpro.setCandidateFilePath(tmp_dir_candidate_path.c_str());

        genpro.process_organism();

        /// 3. Processing the GFF3 file
        cout << "\nNow generating IT's for every Chromosome...\n";
        std::string organism;
        std::map<std::string, ITNode*> node_map;
        organism = genpro.getOrganism();
        GFFReader gff_reader = GFFReader(argv[4], organism);

        gff_reader.processGFFFile(node_map);
        cout << "--- Generating complete! ---\n\n";

        Statistics stats = {0,0,0,0};

        std::string mirna_ident;

        std::vector<SeqAlignment> best_alignment_list;

        int mirna_counter = 1;

        for(auto it: miRNA_list){

            /// 4. Search for candidates
            genpro.search_for_candidates(genpro.getTempGenFilePath(), it);

            /// 5. Apply algorithm
            fstream candidate_file;
            candidate_file.open(genpro.getCandidateFile());

            std::string line, trash;
            std::string candidate_sequence;
            std::string overhang_str;
            std::string chr_id;
            size_t start, end;
            size_t count = 0;
            size_t key = 1;

            std::map<size_t, SeqAlignment> all_alignments;

            if(!candidate_file.is_open()){
                cerr << "Could not open the temporary candidate file!" << endl;
                exit(EXIT_FAILURE);
            }

            while(getline(candidate_file, line)){

                istringstream iss(line);

                iss >> chr_id >> start >> end >> candidate_sequence >> overhang_str;

                if(evaluate_candidates(candidate_sequence, overhang_str, it, match, gu, penalty, gap_open, gap_ext, threshold, scal_f, all_alignments, key, wobble_region)){
                    count++;
                    key++;

                    all_alignments.at(key-1).setStart(start);
                    all_alignments.at(key-1).setEnd(end);
                    all_alignments.at(key-1).setChrID(chr_id);

                }
            }

            candidate_file.close();

            if(remove(genpro.getCandidateFile()) != 0) {
                perror("Error deleting temporary candidate file");
            }

            /// 6. Search IT's for precise location
            std::string loc;
            int i_start, i_end;
            Interval res;
            Interval i;
            std::vector<Interval> possible_hits;
            std::vector<ITNode*> possible_node_hits;

            cout << all_alignments.size() << " made it!\n";

            for(auto it_align: all_alignments){

                i_start = (int)it_align.second.getStart();
                i_end = (int)it_align.second.getEnd();

                i = {i_start, i_end};

                std::string align_key = it_align.second.getChrID();

                process_for_overlapSearch_with_nodes(node_map[align_key], i, possible_node_hits);

                for(auto node: possible_node_hits){
                    process_for_overlapSearch(node->gene_root, i, possible_hits);
                }
                res = best_Overlap(possible_hits, i);

                if (res.loc == "no_hit") {
                    loc = "intergenic";
                    all_alignments.at(it_align.first).setLocOnGenome(loc);
                } else {
                    loc = res.loc;
                    all_alignments.at(it_align.first).setLocOnGenome(loc);
                }

                possible_hits.clear();
                possible_hits.shrink_to_fit();

                possible_node_hits.clear();
                possible_node_hits.shrink_to_fit();

            }


            /// 7. Output
            output_table_path.append("\\output_tables\\table_" + it.get_miRNAId() + ".txt");
            output_alignment_path.append("\\output_alignments\\alignment_" + it.get_miRNAId() + ".txt");

            mirna_ident = it.get_miRNAId();

            if(output_stats == 1) {
                print_table(all_alignments, output_table_path.c_str(), mirna_ident, all_alignments.size(), stats);
            }
            if(output_alignment == 1) {
                output_results(all_alignments, output_alignment_path.c_str(), mirna_ident, best_alignment_list);
            }

            output_table_path.clear(), output_alignment_path.clear();
            cout << mirna_counter << " done!\n";

            mirna_counter++;
        }

        if(remove(genpro.getTempGenFilePath()) != 0 ) {
            perror("Error deleting temporary genome file file");
        }

        std::string output_best_alignments = std::string(output);
        std::string output_best_alignments_path = output_best_alignments.append("\\best_alignments.txt");

        if(output_alignment == 1) {
            print_best_alignments(best_alignment_list, output_best_alignments_path.c_str());
        }

        std::string output_overall = std::string(output);
        std::string output_overall_path = output_overall.append("\\overall_table.txt");

        overall_stats_output(stats, miRNA_list.size(), output_overall_path.c_str());

        duration_overall = (std::clock() - start_clock_overall) / (double) CLOCKS_PER_SEC;

        size_t minutes_overall = floor(duration_overall/60);
        double seconds_overall = (duration_overall/60 - minutes_overall) * 60;

        cout << "\nTime needed: " << duration_overall << " seconds, which equals " << minutes_overall << ":" << seconds_overall << " minutes" << endl;

        return 0;

    }

    else {
        cerr << "wrong/unknown command!\nUsage:\n\t"
                << "--align   miRNA.fa   genome.fa   genome_annotation.gff3   config.json   tmp_dir   output_dir\nor\n\t"
                << "type in   --help or -h   for detailed help\n";
        return -1;
    }

 }

