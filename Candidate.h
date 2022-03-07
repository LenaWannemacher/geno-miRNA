#ifndef MIRNA_CANDIDATE_H
#define MIRNA_CANDIDATE_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <map>

class Candidate{

private:

    /**
     * Contains the start and end of the sequence as well as the sequence itself and the chromosome ID
     */

    size_t can_start;
    size_t can_end;
    std::string can_seq;
    std::string can_chr_id;

public:

    Candidate(size_t start, size_t end, std::string& seq, std::string& chr_id){
        can_start = start;
        can_end = end;
        can_seq = seq;
        can_chr_id = chr_id;
    }

    /**
     * Getter and Setter functions
     */

    size_t getStart(){return this->can_start;}
    size_t getEnd(){return this->can_end;}
    std::string getSeq(){return this->can_seq;}
    std::string getChrID(){return this->can_chr_id;}

    void setStart(size_t var){can_start = var;}
    void setEnd(size_t var){can_end = var;}
    void setSeq(std::string& var_seq){can_seq = var_seq;}
    void setChrID(std::string& chr_id){can_chr_id = chr_id;}

    ~Candidate() = default;

};




#endif //MIRNA_CANDIDATE_H