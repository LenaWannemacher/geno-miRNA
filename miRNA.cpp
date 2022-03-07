#include "InvalidCharacter.h"
#include "miRNA.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>

using namespace std;

/**
 * Constructor getting all the information from the miRNA.fa file
 *
 * @param miRNA_file_path
 */
miRNA::miRNA(const char *miRNA_file_path) {

    fstream miRNA_file;
    string line;

    miRNA_file.open(miRNA_file_path);

    if(!miRNA_file.is_open()){
        cerr << "Could not open miRNA file!" << endl;
    }

    while (getline(miRNA_file, line)){

        if (line.substr(0,1) == ">"){
            miRNA_ID = line.substr(1, line.find(' '));
        }
        else {
            miRNA_sequence.append(line);
        }
    }

    changed_miRNA_seq = miRNA_sequence;

    reversed_miRNA = miRNA_sequence;
    reverse(reversed_miRNA.begin(), reversed_miRNA.end());

}

miRNA::miRNA(std::string& seq, std::string& id){

    miRNA_ID = id;
    miRNA_sequence = seq;
    changed_miRNA_seq = seq;
    reversed_miRNA = seq;
    reverse(reversed_miRNA.begin(), reversed_miRNA.end());

}


/**
 * Generates complementary strand-sequence
 */
void miRNA::generate_complementary_RNA(std::string& seq){

    string origin_seq = seq;
    string complementary_seq(origin_seq.length(), ' ');

    for (int i = 0; i < origin_seq.length(); ++i){
        if (origin_seq[i] == 'A'){
            complementary_seq[i] = 'U';
        }
        else if (origin_seq[i] == 'C'){
            complementary_seq[i] = 'G';
        }
        else if (origin_seq[i] == 'G'){
            complementary_seq[i] = 'C';
        }
        else if (origin_seq[i] == 'U'){
            complementary_seq[i] = 'A';
        }
        else if (origin_seq[i] == 'N'){
            complementary_seq[i] = 'N';
        }
        else {
            throw InvalidCharacter(origin_seq[i]);
        }
    }
    changed_miRNA_seq = complementary_seq;
}

/**
 * Generates complementary strand-sequence
 */
void miRNA::generate_complementary_DNA(std::string& seq){

    string origin_seq = seq;
    string complementary_seq(origin_seq.length(), ' ');

    for (int i = 0; i < origin_seq.length(); ++i){
        if (origin_seq[i] == 'A'){
            complementary_seq[i] = 'T';
        }
        else if (origin_seq[i] == 'C'){
            complementary_seq[i] = 'G';
        }
        else if (origin_seq[i] == 'G'){
            complementary_seq[i] = 'C';
        }
        else if (origin_seq[i] == 'T'){
            complementary_seq[i] = 'A';
        }
        else if (origin_seq[i] == 'N'){
            complementary_seq[i] = 'N';
        }
        else {
            throw InvalidCharacter(origin_seq[i]);
        }
    }
    changed_miRNA_seq = complementary_seq;
}

/**
 * converts a sequence into RNA
 */
void miRNA::convert_to_RNA(std::string& seq){

    string origin_seq = seq;
    string converted_seq(origin_seq.length(), ' ');

    for (int i = 0; i < origin_seq.length(); ++i){
        if (origin_seq[i] == 'A'){
            converted_seq[i] = 'A';
        }
        else if (origin_seq[i] == 'C'){
            converted_seq[i] = 'C';
        }
        else if (origin_seq[i] == 'G'){
            converted_seq[i] = 'G';
        }
        else if (origin_seq[i] == 'T'){
            converted_seq[i] = 'U';
        }
        else if (origin_seq[i] == 'N'){
            converted_seq[i] = 'N';
        }
        else {
            throw InvalidCharacter(origin_seq[i]);
        }
    }
    changed_miRNA_seq = converted_seq;
}

/**
 * converts a sequence into DNA
 */
void miRNA::convert_to_DNA(std::string& seq){

    string origin_seq = seq;
    string converted_seq(origin_seq.length(), ' ');

    for (int i = 0; i < origin_seq.length(); ++i){
        if (origin_seq[i] == 'A'){
            converted_seq[i] = 'A';
        }
        else if (origin_seq[i] == 'C'){
            converted_seq[i] = 'C';
        }
        else if (origin_seq[i] == 'G'){
            converted_seq[i] = 'G';
        }
        else if (origin_seq[i] == 'U'){
            converted_seq[i] = 'T';
        }
        else if (origin_seq[i] == 'N'){
            converted_seq[i] = 'N';
        }
        else {
            throw InvalidCharacter(origin_seq[i]);
        }
    }
    changed_miRNA_seq = converted_seq;
}


/**
 * generic Getter and Setter
 */

string &miRNA::get_miRNASequence() {
    return miRNA_sequence;
}

const string &miRNA::get_miRNAId() const {
    return miRNA_ID;
}

string& miRNA::get_changed_miRNASequence() {
    return changed_miRNA_seq;
}

string& miRNA::get_reversed_miRNA() {
    return reversed_miRNA;
}

string& miRNA::getToken() {
    return token;
}

size_t miRNA::getOverhang() {
    return overhang;
}

/**
 * Catch input errors
 * @param seq
 */
void miRNA::set_miRNASequence(string& seq){

    string new_seq;

    for (char i : seq){
        if (i == 'A'){
            new_seq.append("A");
        }
        else if (i == 'C'){
            new_seq.append("C");
        }
        else if (i == 'G'){
            new_seq.append("G");
        }
        else if (i == 'U'){
            new_seq.append("U");
        }
        else if (i == 'N'){
            new_seq.append("N");
        }
        else {
            throw InvalidCharacter(i);
        }
    }

    this->miRNA_sequence = new_seq;

}

void miRNA::set_miRNAId(string& id){
    miRNA_ID = id;
}

void miRNA::set_changed_miRNASequence(string &seq) {
    changed_miRNA_seq = seq;
}

void miRNA::set_miRNAToken(size_t token_size) {

    if (token_size > changed_miRNA_seq.length() || token_size < 1){
        cerr << "token_size must be in range of miRNA_sequence size!" << endl;
    }
    else{
        token = changed_miRNA_seq.substr(changed_miRNA_seq.size()-token_size-1, token_size);
    }

}

void miRNA::set_miRNAOverhang(size_t oh) {
    overhang = oh;
}

miRNA::miRNA() = default;
