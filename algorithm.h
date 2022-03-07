#ifndef MIRNA_ALIGNMENT_H
#define MIRNA_ALIGNMENT_H

#include "InvalidCharacter.h"
#include "Matrix.h"
#include "miRNA.h"
#include "Candidate.h"
#include "SeqAlignment.h"
#include "StringMatrix.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <vector>
#include <direct.h>


/**
 * generates the alignments of two given strings
 *
 * @param mirna The given miRNA string
 * @param dna   The given DNA string
 * @return The string between the alignments, showing matches, GU wobbles and mismatches
 */
std::string getAlignment(string& mirna, string& dna){

    std::string between;

    for(size_t x = 0; x < mirna.size(); x++){

        if(mirna.at(x) == 'A' && dna.at(x) == 'T'){
            between += '|';
        }
        else if (mirna.at(x) == 'C' && dna.at(x) == 'G'){
            between += '|';
        }
        else if (mirna.at(x) == 'G' && dna.at(x) == 'C'){
            between += '|';
        }
        else if (mirna.at(x) == 'U' && dna.at(x) == 'A'){
            between += '|';
        }
        else if (mirna.at(x) == 'U' && dna.at(x) == 'G'){
            between += ':';
        }
        else if(mirna.at(x) == 'G' && dna.at(x) == 'T'){
            between += ':';
        }
        else {
            between += ' ';
        }
    }

    return between;

}

/**
 * counts mismatches of an alignment
 *
 * @param between   The string containing information about matches, mismatches and GU wobbles
 * @return The amount of mismatches
 */
size_t count_mismatches(std::string& between){

    size_t mismatches = 0;

    for(char x : between){
        if(x == ' '){
            mismatches++;
        }
    }

    return mismatches;
}

/**
 * checks for site complementary
 *
 * @param between   The string containing information about matches, mismatches and GU wobbles
 * @return The complementary sites (seed, 3' or seed + 3') || seed match is always given
 */
std::string complementary_sites(std::string& between){

    bool seed_match = true;
    bool threeprime_match = true;

    std::string rev_between = string(between.rbegin(),between.rend());

    for(size_t x = 1; x <= 6; x++){
        if(rev_between[x] == ' '){
            seed_match = false;
        }
    }

    for(size_t x = 11; x <= 16; x++){
        if(rev_between[x] == ' '){
            threeprime_match = false;
        }
    }

    if(seed_match && threeprime_match){return "seed + 3'";}
    else if(seed_match){return "seed";}
    else if(threeprime_match){return "3'";}
    else {return "none";}

}

/**
 * computes the gap penalty
 *
 * @param mirna The given miRNA string
 * @param dna   The given DNA string
 * @param gap_open  Penalty for opening a gap
 * @param gap_ext   Penalty for extending a gap
 * @param scal_f    Scaling factor for the first n positions
 * @return The gap penalty score of the given alignment
 */
float computeGapPenalty(string& mirna, string& dna, int gap_open, int gap_ext, int penalty, float scal_f){

    bool gap_opened = true;
    bool gap_extended = false;
    float score = 0;

    std::string rev_mirna = string(mirna.rbegin(),mirna.rend());
    std::string rev_dna = string(dna.rbegin(), dna.rend());

    for(size_t x = 0; x < rev_mirna.size(); x++){

        if(x<11) {
            if (rev_mirna.at(x) == '-') {
                if (gap_opened) {
                    score += (((float)gap_open * scal_f) + (float)penalty);
                    gap_opened = false;
                    gap_extended = true;
                } else if (gap_extended) {
                    score += (((float)gap_ext * scal_f) + (float)penalty);
                }
            } else {
                gap_opened = true;
                gap_extended = false;
            }
        }
        else {
            if (rev_mirna.at(x) == '-') {
                if (gap_opened) {
                    score += (float) gap_open + (float)penalty;
                    gap_opened = false;
                    gap_extended = true;
                } else if (gap_extended) {
                    score += (float) gap_ext;
                }
            } else {
                gap_opened = true;
                gap_extended = false;
            }
        }

    }

    gap_opened = true;
    gap_extended = false;

    for(size_t x = 0; x < rev_dna.size(); x++){

        if(x<11) {
            if (rev_dna.at(x) == '-') {
                if (gap_opened) {
                    score += (((float)gap_open * scal_f) + (float)penalty);
                    gap_opened = false;
                    gap_extended = true;
                } else if (gap_extended) {
                    score += (((float)gap_ext * scal_f) + (float)penalty);
                }
            } else {
                gap_opened = true;
                gap_extended = false;
            }
        }
        else {
            if (rev_dna.at(x) == '-') {
                if (gap_opened) {
                    score += (float) gap_open + (float)penalty;
                    gap_opened = false;
                    gap_extended = true;
                } else if (gap_extended) {
                    score += (float) gap_ext;
                }
            } else {
                gap_opened = true;
                gap_extended = false;
            }
        }

    }

    return score;

}

/**
 * compares two chars (miRNA and DNA) to determine if it is a match, mismatch or GU wobble
 *
 * @param rna   Char of the given miRNA sequence
 * @param dna   Char of the given DNA sequence
 * r and c used to check where we are in the alignment and if it is necessary to use our scaling factor
 * @param r
 * @param c
 * @param match The match score
 * @param gu    The GU wobble score
 * @param penalty   The penalty score
 * @param scal_f    Scaling factor for the first n positions
 * @return The score of the two given chars
 */
float compareChars(char rna, char dna, size_t r, size_t c, size_t match, size_t gu, int penalty, float scal_f){

    /**
     * last change: added G-T wobble (RNA-DNA)
     *
     * r,c >= 11 für scaling factor Abfrage
     */

    float score = 0;

    if(rna == 'A'){
        if(dna == 'T'){
            score = (float)match;
        } else { score = (float)penalty; }
    }
    else if(rna == 'C'){
        if(dna == 'G'){
            score = (float)match;
        } else { score = (float)penalty; }
    }
    else if(rna == 'G'){
        if (dna == 'C'){
            score = (float)match;
        }
        else if (dna == 'T'){
            score = (float)gu;
        }
        else { score = (float)penalty; }
    }
    else if(rna == 'U'){
        if(dna == 'A'){
            score = (float)match;
        }
        else if (dna == 'G'){
            score = (float)gu;
        } else { score = (float)penalty; }
    }
    else if (rna == '-' || dna == '-'){
        score = (float)penalty;
    }
    if (r>=11 && c>=11){
        score *= scal_f;
    }
    return score;
}

/**
 * generates our scoring matrix
 *
 * @param rows
 * @param cols
 * @param penalty   Used to initialize the matrix
 * @return Empty, initialized scoring matrix
 */
Matrix generateScoringMat(size_t rows, size_t cols, int penalty){

    // set default score to 0.0
    Matrix mat = Matrix(rows,cols,0.0);

    // Initializing the scoring matrix
    for(int c = 0; c < cols; c++){
        mat(0,c) = (float)(penalty*c);
    }
    for(int r = 0; r < rows; r++){
        mat(r,0) = (float)(penalty*r);
    }

    return mat;
}

/**
 * generates the memory matrix, used to keep track of where we came from in the scoring matrix.
 * (needed to perform the traceback)
 *
 * @param rows
 * @param cols
 * @return Empty, initialized memory matrix
 */
Matrix generateMemoryMat(size_t rows, size_t cols){

    Matrix mat = Matrix(rows, cols, 0);

    // 1 = left, 2 = up, 3 = up-left, 0 = left and up equal
    for(size_t c = 1; c < cols; c++){
        mat(0,c) = 1;
    }
    for(size_t r = 1; r < rows; r++){
        mat(r,0) = 2;
    }

    return mat;
}

/**
 * checks if an alignment passes certain rules
 *
 * @param between   The string containing information about matches, mismatches and GU wobbles
 * @return True if the alignment passes the rules, otherwise False
 */
bool rules_check(string& between, size_t wobble_region){

    std::string rev_between = string(between.rbegin(), between.rend());

    // rule 1: fewer than 5 mismatches at positions 3-12
    size_t count_mismatches = 0;
    for(size_t x = 2; x < 12; x++){
        if(rev_between[x] == ' '){
            count_mismatches++;
        }
    }
    if(count_mismatches >= 5){
        return false;
    }

    count_mismatches = 0;

    // rule 2: at least one mismatch 9 - (length-5)
    for(size_t x = 8; x < rev_between.size()-5; x++){
        if(rev_between[x] == ' '){
            count_mismatches++;
        }
    }
    if(count_mismatches == 0){
        return false;
    }

    count_mismatches = 0;

    // rule 3: fewer than two mismatches in the last five positions
    for(size_t x = rev_between.size()-5; x < rev_between.size(); x++){
        if(rev_between[x] == ' '){
            count_mismatches++;
        }
    }
    if(count_mismatches >= 2){
        return false;
    }

    count_mismatches = 0;

    // rule 4: wobble region (optional, default off)
    if(wobble_region == 1) {
        for (size_t x = 7; x < 14; x++) {
            if (rev_between[x] == ' ') {
                count_mismatches++;
            }
        }
        if(count_mismatches != 0){
            return false;
        }
    }

    return true;
}

/**
 * determines the binding site motif of an alignment
 *
 * @param dna   The given DNA string
 * @param between   The string containing information about matches, mismatches and GU wobbles
 * @return The binding site motif
 */
std::string binding_site_motif(string& dna, string& between){

    std::string binding_site;

    size_t match_count = 0;

    std::string rev_dna = string(dna.rbegin(), dna.rend());
    std::string rev_between = string(between.rbegin(), between.rend());

    for(size_t x = 1; x < 7; x++){
        if(rev_between[x] == '|'){
            match_count++;
        }
    }

    if(match_count == 6){                               // seed positions match (2-7)

        if(rev_between[7] == '|'){                          // position 8 matches

            if(rev_dna[0] == 'A'){                          // seed positions + 8 + A at position 1
                return binding_site = "8mer";
            }
            else{
                return binding_site = "7mer-m8";        // seed position + 8 matches
            }
        }
        else{
            if(rev_dna[0] == 'A'){                          // seed positions match + A at position 1, but no match at 8
                return binding_site = "7mer-A1";
            }
        }

        return binding_site = "6mer";                   // only seed positions match

    }
    else {
        return binding_site = "none";                   // no canonical binding site
    }

}
/**
 * computes the final alignment score
 *
 * @param mirna The given miRNA sequence
 * @param dna   The given DNA sequence
 * @param match The match score
 * @param gu    The GU wobble score
 * @param penalty   The penalty score
 * @param gap_open  Penalty for opening a gap
 * @param gap_ext   Penalty for extending a gap
 * @param scal_f    Scaling factor for the first n positions
 * @return  The final alignment score
 */
float computeAlignmentScore(std::string& mirna, std::string& dna, size_t match, size_t gu, int penalty, int gap_open, int gap_ext, float scal_f){

    std::string rev_mirna = string(mirna.rbegin(),mirna.rend());
    std::string rev_dna = string(dna.rbegin(),dna.rend());

    bool gap_opened = true;
    bool gap_extended = false;

    float score = 0;
    float scaling_score = 0;

    for(size_t x = 0; x < rev_mirna.length(); x++){

        if(x == 11){
            scaling_score = score;
            score = 0;
        }

        if(rev_mirna[x] == '-' || rev_dna[x] == '-'){
            if(gap_opened){
                score += (float)gap_open;
                gap_opened = false;
                gap_extended = true;
            }
            else if (gap_extended){
                score += (float)gap_ext;
            }
        }
        else if(rev_mirna[x] == 'A') {
            if(rev_dna[x] == 'T') {
                score += (float)match;
                gap_opened = true;
                gap_extended = false;
            } else {
                score += (float)penalty;
                gap_opened = true;
                gap_extended = false;
            }
        }
        else if(rev_mirna[x] == 'C'){
            if(rev_dna[x] == 'G'){
                score += (float)match;
                gap_opened = true;
                gap_extended = false;
            } else {
                score += (float)penalty;
                gap_opened = true;
                gap_extended = false;
            }
        }
        else if(rev_mirna[x] == 'G'){
            if(rev_dna[x] == 'C'){
                score += (float)match;
                gap_opened = true;
                gap_extended = false;
            } else if (rev_dna[x] == 'T') {
                score += (float)gu;
                gap_opened = true;
                gap_extended = false;
            } else {
                score += (float)penalty;
                gap_opened = true;
                gap_extended = false;
            }
        }
        else if(rev_mirna[x] == 'U'){
            if(rev_dna[x] == 'A'){
                score += (float)match;
                gap_opened = true;
                gap_extended = false;
            } else if (rev_dna[x] == 'G'){
                score += (float)gu;
                gap_opened = true;
                gap_extended = false;
            } else {
                score += (float)penalty;
                gap_opened = true;
                gap_extended = false;
            }
        }
        else {
            score += (float)penalty;
            gap_opened = true;
            gap_extended = false;
        }
    }

    score += (scal_f * scaling_score);
    return score;

}

/**
 * The algorithm. Same as above, but instead of the object "Candidate" we only get the sequence needed
 *
 * @param candidate_seq The DNA sequence of our Candidate
 * @param mirna The given miRNA sequence
 * @param match The match score
 * @param gu    The GU wobble score
 * @param penalty   The penalty score
 * @param gap_open  Penalty for opening a gap
 * @param gap_ext   Penalty for extending a gap
 * @param threshold Threshold to determine whether or not an alignment is good enough
 * @param scal_f    Scaling factor for the first n positions
 * @param seqalignments The map of all final sequence alignments
 * @param key   Key for the map item
 * @return True, if the given Candidate is good enough, else False
 */
bool evaluate_candidates(string& candidate_seq, string& overhang_str, miRNA& mirna, size_t match, size_t gu, int penalty, int gap_open, int gap_ext, size_t threshold, float scal_f, std::map<size_t, SeqAlignment>& seqalignments, size_t key, size_t wobble_region){

    string seq1 = mirna.get_reversed_miRNA();
    string seq2 = candidate_seq;

    size_t cols = seq1.size()+1;        // miRNA (top)
    size_t rows = seq2.size()+1;        // DNA (side)

    Matrix scoring_mat = generateScoringMat(rows, cols, penalty);
    Matrix memory_mat = generateMemoryMat(rows, cols);

    float score, max, left, up, upleft;

    // filling scoring matrix and memory matrix
    for(size_t r = 1; r < rows; r++){
        for(size_t c = 1; c < cols; c++){

            char rna = seq1.at(c-1);
            char dna = seq2.at(r-1);

            upleft = scoring_mat(r-1,c-1) + (compareChars(rna, dna, r, c, match, gu, penalty, scal_f));
            up = scoring_mat(r-1,c) + (compareChars('-',dna, r, c, match, gu, penalty, scal_f));
            left = scoring_mat(r,c-1) + (compareChars(rna , '-', r, c, match, gu, penalty, scal_f));

            if (upleft >= up && upleft >= left) {
                score = upleft;
                memory_mat(r,c) = 3;
                max = upleft;
            }
            else if (up > upleft && up > left){
                score = up;
                memory_mat(r,c) = 2;
                max = up;
            }
            else if (left > upleft && left > up){
                score = left;
                memory_mat(r,c) = 1;
                max = left;
            } else if (left == up) {
                score = left;
                memory_mat(r,c) = 1;
                max = left;
            }

            scoring_mat(r,c) = score;

        }
    }

    std::string aligned_miRNA, aligned_DNA;

    size_t r = rows-1; size_t c = cols-1;

    bool fillup = false;

    // traceback
    for(r; r >= 0; --r){

        if(r == 0 && c > 0){
            fillup = true;
            break;
        }

        for(c; c > 0; --c){

            if(memory_mat(r,c) == 3){
                aligned_miRNA += seq1.at(c-1);
                aligned_DNA += seq2.at(r-1);
                c--;

                break;
            }
            else if(memory_mat(r,c) == 2){
                aligned_miRNA += '-';
                aligned_DNA += seq2.at(r-1);

                break;
            }
            else if (memory_mat(r,c) == 1 || memory_mat(r,c) == 0){
                aligned_miRNA += seq1.at(c-1);
                aligned_DNA += '-';

            }
        }
        if(c == 0){
            break;
        }
    }

    if(fillup){
        for(c; c > 0; --c){
            aligned_miRNA += seq1.at(c-1);

            if(c <= overhang_str.length()){
                aligned_DNA += overhang_str.back();
                overhang_str.pop_back();
            } else {
                aligned_DNA += '-';
            }
        }
    }

    reverse(aligned_DNA.begin(), aligned_DNA.end());
    reverse(aligned_miRNA.begin(), aligned_miRNA.end());

    std::string between = getAlignment(aligned_miRNA, aligned_DNA);

    std::string bs_motif = binding_site_motif(aligned_DNA, between);

    std::string compl_site = complementary_sites(between);
    size_t mismatches = count_mismatches(between);

    /// changed to computeAlignmentScore !!!
    float align_score = computeAlignmentScore(aligned_miRNA, aligned_DNA, match, gu, penalty, gap_open, gap_ext, scal_f);

    /// rules check new
    if(rules_check(between, wobble_region)){
        if(align_score >= (float)threshold) {

            if(bs_motif == "none") {
                return false;
            }

            SeqAlignment my_seqalignment = SeqAlignment(aligned_miRNA, aligned_DNA, between, bs_motif, compl_site, mismatches, align_score);
            seqalignments.insert(pair<size_t, SeqAlignment>(key, my_seqalignment));

            return true;

        }
        else{
            return false;
        }
    }
    else {
        return false;
    }
}


#endif //MIRNA_ALIGNMENT_H
