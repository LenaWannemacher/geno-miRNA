#include "GenomeProcessor.h"
#include "miRNA.h"
#include "Candidate.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <vector>
#include <algorithm>
#include <regex>

using namespace std;

GenomeProcessor::GenomeProcessor(const char *genome_file_path) {
    this->genome_file_path = genome_file_path;
}

/**
 * Processes the given genome .fa-file and creates a temporary genome .fa-file with no '\n' to make it easier to search
 * for our Token-Sequence
 */
void GenomeProcessor::process_organism() {

    fstream genome_file;
    ofstream temp_genome;
    temp_genome.open(temp_genome_file_path);

    string line;
    string org;
    string chr_id;
    string sequence;

    size_t count = 1;
    size_t pos;
    size_t org_begin;

    bool first_found = true;
    bool chr_to_create = false;

    genome_file.open(genome_file_path);

    if (!genome_file.is_open()){
        cerr << "Could not open genome file!" << endl;
        exit(EXIT_FAILURE);
    }

    while(getline(genome_file, line)){
        if(line.substr(0,1) == ">"){
            pos = line.find("chromosome");
            if(pos != string::npos){
                if(first_found){
                    org_begin = line.find_first_of(':', 15);
                    organism = line.substr(org_begin+1, (line.find_first_of(':',org_begin+1)) - org_begin-1);
                    first_found = false;
                }
                if (chr_to_create){
                    temp_genome << ">" << chr_id << "\n" << sequence << endl;
                    count++;
                }
                chr_id = line.substr(1, line.find(' ',1)); // hier evtl die 1 bei find wegmachen
                chr_to_create = true;
                sequence.clear();
            }
            else {
                break;
            }
        }
        else {
            sequence.append(line);
        }
    }

    temp_genome << ">" << chr_id << "\n" << sequence;
    temp_genome.close();

    cout << "--- Processing successful! ---" << endl;
}

/**
 * Searches for the Token-Sequence in the temp_genome file and creates a temporary file containing all the Candidates
 *
 * @param temp_gen_file_path    Path to the temporary genome file
 * @param mirna The given miRNA
 */
void GenomeProcessor::search_for_candidates(const char *temp_gen_file_path, miRNA &mirna) {

    string chr_seq;
    string chr_id;

    string dna_seq;
    string token = mirna.getToken();

    size_t start, end, pos;
    size_t mirna_length = mirna.get_miRNASequence().size();
    size_t overlap = mirna.getOverhang();
    string overlap_str;

    ofstream all_candidates_file;
    all_candidates_file.open(candidates_file_path);

    fstream temp_genome;
    temp_genome.open(temp_gen_file_path);

    if(!temp_genome.is_open()){
        cerr << "Could not open the temporary genome storage file!" << endl;
        exit(EXIT_FAILURE);
    }

    while (getline(temp_genome, chr_seq)){

        if(chr_seq.substr(0,1) == ">"){
            chr_id = chr_seq.substr(1, chr_seq.find(' ',1));
        }
        else {
            /// experiment
            pos = chr_seq.find(token);
            while (pos != std::string::npos){

                if(pos < (mirna_length + overlap - token.size()-1)){
                    pos = chr_seq.find(token, pos + token.size());
                    continue;
                }

                start = pos - mirna_length - overlap + token.size() + 1;
                end = pos + token.size();

                dna_seq = chr_seq.substr(start+overlap, mirna_length);

                overlap_str = chr_seq.substr(start, overlap);

                all_candidates_file << chr_id << " " << start << " " << end << " " << dna_seq << " " << overlap_str << "\n";

                pos = chr_seq.find(token, pos + token.size());
                dna_seq.clear();
            }
        }
    }

    temp_genome.close();
    all_candidates_file.close();
}
