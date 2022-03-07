#ifndef MIRNA_SEQALIGNMENT_H
#define MIRNA_SEQALIGNMENT_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <vector>
#include <algorithm>

/// think about making it a struct??

class SeqAlignment{

private:

    size_t seq_start{}, seq_end{};
    std::string sequence_1, sequence_2, between;
    std::string chr_id;
    std::string bs_motif;                            //6mer, 7mer-A1, 7mer-m8, 8mer
    std::string loc_on_genome;
    std::string complementary_site;         // seed, 3' or seed + 3'  (seed will be always true)
    size_t mismatches{};
    float score{};

    std::string mirna_id{};         // usually only important for statistics -> best alignments (either all together or only of "chunk")

public:

    SeqAlignment(std::string& seq1, std::string& seq2, std::string& betweeen, std::string& bs_motif, std::string& compl_site, size_t mismatch, float align_score){
        this->sequence_1 = seq1;
        this->sequence_2 = seq2;
        this->between = betweeen;
        this->bs_motif = bs_motif;
        this->complementary_site = compl_site;
        this->mismatches = mismatch;
        this->score = align_score;
    }

    /**
     * Generic Getter and Setter
     */
    void setStart(size_t start){this->seq_start = start;}
    void setEnd(size_t end){this->seq_end = end;}
    void setSeq1(std::string& seq){this->sequence_1 = seq;}
    void setSeq2(std::string& seq){this->sequence_2 = seq;}
    void setBetweenSeq(std::string& between_seq){this->between = between_seq;}
    void setChrID(std::string& chrID){this->chr_id = chrID;}
    void setBSMotif(std::string& binding_site_motif){this->bs_motif = binding_site_motif;}
    void setLocOnGenome(std::string& loc_on_gen){this->loc_on_genome = loc_on_gen;}
    void setComplementarySite(std::string& compl_site){this->complementary_site = compl_site;}
    void setMismatches(size_t mismatch){this->mismatches = mismatch;}
    void setScore(float align_score){this->score = align_score;}
    void setmiRNAID(std::string& m_id){this->mirna_id = m_id;}

    size_t getStart(){return this->seq_start;}
    size_t getEnd(){return this->seq_end;}
    size_t getMismatches(){return this->mismatches;}
    std::string getSeq1(){return this->sequence_1;}
    std::string getSeq2(){return this->sequence_2;}
    std::string getBetweenSeq(){return this->between;}
    std::string getChrID(){return this->chr_id;}
    std::string getBSMotif(){return this->bs_motif;}
    std::string getLocOnGenome(){return this->loc_on_genome;}
    std::string getComplementarySite(){return this->complementary_site;}
    float getScore(){return this->score;}
    std::string getmiRNAID(){return this->mirna_id;}

};


#endif //MIRNA_SEQALIGNMENT_H
