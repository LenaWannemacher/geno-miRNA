#ifndef MIRNA_GENOMEPROCESSOR_H
#define MIRNA_GENOMEPROCESSOR_H

#include "Candidate.h"
#include "miRNA.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <regex>
#include <map>

class GenomeProcessor{

private:

    std::string organism;
    const char* genome_file_path;
    const char* temp_genome_file_path;
    std::map<size_t, Candidate*> all_candidates;
    const char* candidates_file_path;

public:

    explicit GenomeProcessor(const char* genome_file_path);

    /**
     * Getter and Setter
     */
    std::string getOrganism(){return this->organism;}
    const char* getFilePath(){return this->genome_file_path;}
    const char* getTempGenFilePath(){return this->temp_genome_file_path;}
    const char* getCandidateFile(){return this->candidates_file_path;}
    std::map<size_t, Candidate*> getAllCandidates(){return this->all_candidates;}

    void setOrganism(std::string& org){this->organism = org;}
    void setFilePath(const char* file_path){this->genome_file_path = file_path;}
    void setTempGenFilePath(const char* temp_file_path){this->temp_genome_file_path = temp_file_path;}
    void setCandidateFilePath(const char* can_file_path){this->candidates_file_path = can_file_path;}

    /**
     * Processes genome .fa-file and creates a temporary genome file with no '\n'
     */
    void process_organism();

    /**
     * Searches for the Token-Sequence in the temp_genome file and creates a temporary file containing all the Candidates
     *
     * @param temp_gen_file_path    Path to the temporary genome file
     * @param mirna The given miRNA
     */
    void search_for_candidates(const char* temp_gen_file_path, miRNA& mirna);

};

#endif //MIRNA_GENOMEPROCESSOR_H
