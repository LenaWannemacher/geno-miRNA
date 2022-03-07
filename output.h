#ifndef MIRNA_OUTPUT_H
#define MIRNA_OUTPUT_H

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


struct Statistics{

    size_t total_hits{}, total_8mer{}, total_threeprime{}, total_8mer_threeprime{};
    size_t total_6mer{}, total_7merA1{}, total_7merm8{};

    size_t inter_six{}, inter_seven_a{}, inter_seven_m{}, inter_eight{};
    size_t lnc_six{}, lnc_seven_a{}, lnc_seven_m{}, lnc_eight{};
    size_t mrna_six{}, mrna_seven_a{}, mrna_seven_m{}, mrna_eight{};
    size_t three_six{}, three_seven_a{}, three_seven_m{}, three_eight{};
    size_t five_six{}, five_seven_a{}, five_seven_m{}, five_eight{};
    size_t cds_six{}, cds_seven_a{}, cds_seven_m{}, cds_eight{};
    size_t exon_six{}, exon_seven_a{}, exon_seven_m{}, exon_eight{};
    size_t v_gene_six{}, v_gene_seven_a{}, v_gene_seven_m{}, v_gene_eight{};
    size_t pseudo_six{}, pseudo_seven_a{}, pseudo_seven_m{}, pseudo_eight{};

    std::vector<SeqAlignment> best_alignment_list;

};

/**
 * Prints out the final results into given file (currently "output.txt")
 *
 * @param alignments    The map with all the final alignments
 */
void print_output(std::map<size_t, SeqAlignment>& alignments){

    /// setting counters for every case
    size_t inter_six = 0; size_t inter_seven_a = 0; size_t inter_seven_m = 0; size_t inter_eight = 0;
    size_t lnc_six = 0; size_t lnc_seven_a = 0; size_t lnc_seven_m = 0; size_t lnc_eight = 0;
    size_t mrna_six = 0; size_t mrna_seven_a = 0; size_t mrna_seven_m = 0; size_t mrna_eight = 0;
    size_t three_six = 0; size_t three_seven_a = 0; size_t three_seven_m = 0; size_t three_eight = 0;
    size_t five_six = 0; size_t five_seven_a = 0; size_t five_seven_m = 0; size_t five_eight = 0;
    size_t cds_six = 0; size_t cds_seven_a = 0; size_t cds_seven_m = 0; size_t cds_eight = 0;
    size_t exon_six = 0; size_t exon_seven_a = 0; size_t exon_seven_m = 0;  size_t exon_eight = 0;

    std::string loc_on_genome, bsm;

    for(auto it : alignments){

        loc_on_genome = it.second.getLocOnGenome();
        bsm = it.second.getBSMotif();

        if(loc_on_genome == "intergenic" && bsm == "6mer") { inter_six++; }
        if(loc_on_genome == "intergenic" && bsm == "7mer-A1") { inter_seven_a++; }
        if(loc_on_genome == "intergenic" && bsm == "7mer-m8") { inter_seven_m++; }
        if(loc_on_genome == "intergenic" && bsm == "8mer") { inter_eight++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "6mer") { lnc_six++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "7mer-A1") { lnc_seven_a++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "7mer-m8") { lnc_seven_m++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "8mer") { lnc_eight++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "6mer") { three_six++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "7mer-A1") { three_seven_a++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "7mer-m8") { three_seven_m++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "8mer") { three_eight++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "6mer") { five_six++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "7mer-A1") { five_seven_a++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "7mer-m8") { five_seven_m++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "8mer") { five_eight++; }
        if(loc_on_genome == "mRNA" && bsm == "6mer") { mrna_six++; }
        if(loc_on_genome == "mRNA" && bsm == "7mer-A1") { mrna_seven_a++; }
        if(loc_on_genome == "mRNA" && bsm == "7mer-m8") { mrna_seven_m++; }
        if(loc_on_genome == "mRNA" && bsm == "8mer") { mrna_eight++; }
        if(loc_on_genome == "CDS" && bsm == "6mer") { cds_six++; }
        if(loc_on_genome == "CDS" && bsm == "7mer-A1") { cds_seven_a++; }
        if(loc_on_genome == "CDS" && bsm == "7mer-m8") { cds_seven_m++; }
        if(loc_on_genome == "CDS" && bsm == "8mer") { cds_eight++; }
        if(loc_on_genome == "exon" && bsm == "6mer") { exon_six++; }
        if(loc_on_genome == "exon" && bsm == "7mer-A1") { exon_seven_a++; }
        if(loc_on_genome == "exon" && bsm == "7mer-m8") { exon_seven_m++; }
        if(loc_on_genome == "exon" && bsm == "8mer") { exon_eight++; }

    }

    /// initializing the "output_matrix"
    StringMatrix strmat = StringMatrix(11,17," ");      //13,19
    strmat(1,3) = "intergenic"; strmat(1,5) = "lnc_RNA"; strmat(1,7) = "mRNA"; strmat(1,9) = "3'UTR"; strmat(1,11) = "5'UTR"; strmat(1,13) = "CDS";
    strmat(1,15) = "exon";
    strmat(3,1) = "6mer"; strmat(5,1) = "7mer-A1"; strmat(7,1) = "7mer-m8"; strmat(9,1) = "8mer";

    /// filling the "output_matrix"
    strmat(3,3) = to_string(inter_six);      strmat(3,5) = to_string(lnc_six) ;    strmat(3,7) = to_string(mrna_six);     strmat(3,9) = to_string(three_six);     strmat(3,11) = to_string(five_six);     strmat(3,13) = to_string(cds_six);     strmat(3,15) = to_string(exon_six);
    strmat(5,3) = to_string(inter_seven_a) ; strmat(5,5) = to_string(lnc_seven_a); strmat(5,7) = to_string(mrna_seven_a); strmat(5,9) = to_string(three_seven_a); strmat(5,11) = to_string(five_seven_a); strmat(5,13) = to_string(cds_seven_a); strmat(5,15) = to_string(exon_seven_a);
    strmat(7,3) = to_string(inter_seven_m);  strmat(7,5) = to_string(lnc_seven_m); strmat(7,7) = to_string(mrna_seven_m); strmat(7,9) = to_string(three_seven_m); strmat(7,11) = to_string(five_seven_m); strmat(7,13) = to_string(cds_seven_m); strmat(7,15) = to_string(exon_seven_m);
    strmat(9,3) = to_string(inter_eight);    strmat(9,5) = to_string(lnc_eight);   strmat(9,7) = to_string(mrna_eight);   strmat(9,9) = to_string(three_eight);   strmat(9,11) = to_string(five_eight);   strmat(9,13) = to_string(cds_eight);   strmat(9,15) = to_string(exon_eight);

    /// calling the printToFile function
    strmat.printMatrixToFile(R"(C:\\Users\\Snappy\\Documents\\Bachelor\\output_table.txt)", true);

}
/**
 * Prints the statistics table
 *
 * @param alignments    The map containing all the alignments
 * @param output_path   The path to the output file
 * @param mirna_id  The miRNA ID
 * @param all_hits  The total number of hits for this miRNA
 * @param stats The statistics object, keeping track for the overall statistics
 */
void print_table(std::map<size_t, SeqAlignment>& alignments, const char* output_path, std::string& mirna_id, size_t all_hits, Statistics& stats){

    /// setting counters for every case
    size_t inter_six = 0; size_t inter_seven_a = 0; size_t inter_seven_m = 0; size_t inter_eight = 0;
    size_t lnc_six = 0; size_t lnc_seven_a = 0; size_t lnc_seven_m = 0; size_t lnc_eight = 0;
    size_t mrna_six = 0; size_t mrna_seven_a = 0; size_t mrna_seven_m = 0; size_t mrna_eight = 0;
    size_t three_six = 0; size_t three_seven_a = 0; size_t three_seven_m = 0; size_t three_eight = 0;
    size_t five_six = 0; size_t five_seven_a = 0; size_t five_seven_m = 0; size_t five_eight = 0;
    size_t cds_six = 0; size_t cds_seven_a = 0; size_t cds_seven_m = 0; size_t cds_eight = 0;
    size_t exon_six = 0; size_t exon_seven_a = 0; size_t exon_seven_m = 0;  size_t exon_eight = 0;
    size_t v_gene_six = 0; size_t v_gene_seven_a = 0; size_t v_gene_seven_m = 0; size_t v_gene_eight = 0;
    size_t pseudo_six = 0; size_t pseudo_seven_a = 0; size_t pseudo_seven_m = 0; size_t pseudo_eight = 0;

    std::string loc_on_genome, bsm;

    for(auto it : alignments){

        loc_on_genome = it.second.getLocOnGenome();
        bsm = it.second.getBSMotif();

        if(loc_on_genome == "intergenic" && bsm == "6mer") { inter_six++; }
        if(loc_on_genome == "intergenic" && bsm == "7mer-A1") { inter_seven_a++; }
        if(loc_on_genome == "intergenic" && bsm == "7mer-m8") { inter_seven_m++; }
        if(loc_on_genome == "intergenic" && bsm == "8mer") { inter_eight++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "6mer") { lnc_six++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "7mer-A1") { lnc_seven_a++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "7mer-m8") { lnc_seven_m++; }
        if(loc_on_genome == "lnc_RNA" && bsm == "8mer") { lnc_eight++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "6mer") { three_six++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "7mer-A1") { three_seven_a++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "7mer-m8") { three_seven_m++; }
        if(loc_on_genome == "three_prime_UTR" && bsm == "8mer") { three_eight++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "6mer") { five_six++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "7mer-A1") { five_seven_a++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "7mer-m8") { five_seven_m++; }
        if(loc_on_genome == "five_prime_UTR" && bsm == "8mer") { five_eight++; }
        if(loc_on_genome == "mRNA" && bsm == "6mer") { mrna_six++; }
        if(loc_on_genome == "mRNA" && bsm == "7mer-A1") { mrna_seven_a++; }
        if(loc_on_genome == "mRNA" && bsm == "7mer-m8") { mrna_seven_m++; }
        if(loc_on_genome == "mRNA" && bsm == "8mer") { mrna_eight++; }
        if(loc_on_genome == "CDS" && bsm == "6mer") { cds_six++; }
        if(loc_on_genome == "CDS" && bsm == "7mer-A1") { cds_seven_a++; }
        if(loc_on_genome == "CDS" && bsm == "7mer-m8") { cds_seven_m++; }
        if(loc_on_genome == "CDS" && bsm == "8mer") { cds_eight++; }
        if(loc_on_genome == "exon" && bsm == "6mer") { exon_six++; }
        if(loc_on_genome == "exon" && bsm == "7mer-A1") { exon_seven_a++; }
        if(loc_on_genome == "exon" && bsm == "7mer-m8") { exon_seven_m++; }
        if(loc_on_genome == "exon" && bsm == "8mer") { exon_eight++; }
        if(loc_on_genome == "V_gene_segment" && bsm == "6mer") { v_gene_six++; }
        if(loc_on_genome == "V_gene_segment" && bsm == "7mer-A1") { v_gene_seven_a++; }
        if(loc_on_genome == "V_gene_segment" && bsm == "7mer-m8") { v_gene_seven_m++; }
        if(loc_on_genome == "V_gene_segment" && bsm == "8mer") { v_gene_eight++; }
        if(loc_on_genome == "pseudogenic_transcript" && bsm == "6mer") { pseudo_six++; }
        if(loc_on_genome == "pseudogenic_transcript" && bsm == "7mer-A1") { pseudo_seven_a++; }
        if(loc_on_genome == "pseudogenic_transcript" && bsm == "7mer-m8") { pseudo_seven_m++; }
        if(loc_on_genome == "pseudogenic_transcript" && bsm == "8mer") { pseudo_eight++; }

    }

    /// initializing the "output_matrix"
    StringMatrix strmat = StringMatrix(11,21," ");      //13,19
    strmat(1,3) = "intergenic"; strmat(1,5) = "lnc_RNA"; strmat(1,7) = "mRNA"; strmat(1,9) = "3'UTR"; strmat(1,11) = "5'UTR"; strmat(1,13) = "CDS";
    strmat(1,15) = "exon"; strmat(1,17) = "V_gene"; strmat(1,19) = "pseudo";
    strmat(3,1) = "6mer"; strmat(5,1) = "7mer-A1"; strmat(7,1) = "7mer-m8"; strmat(9,1) = "8mer"; //strmat(11,1) = "none";

    /// filling the "output_matrix"
    strmat(3,3) = to_string(inter_six);      strmat(3,5) = to_string(lnc_six) ;    strmat(3,7) = to_string(mrna_six);     strmat(3,9) = to_string(three_six);     strmat(3,11) = to_string(five_six);     strmat(3,13) = to_string(cds_six);     strmat(3,15) = to_string(exon_six);     strmat(3,17) = to_string(v_gene_six);      strmat(3,19) = to_string(pseudo_six);
    strmat(5,3) = to_string(inter_seven_a) ; strmat(5,5) = to_string(lnc_seven_a); strmat(5,7) = to_string(mrna_seven_a); strmat(5,9) = to_string(three_seven_a); strmat(5,11) = to_string(five_seven_a); strmat(5,13) = to_string(cds_seven_a); strmat(5,15) = to_string(exon_seven_a); strmat(5,17) = to_string(v_gene_seven_a);  strmat(5,19) = to_string(pseudo_seven_a);
    strmat(7,3) = to_string(inter_seven_m);  strmat(7,5) = to_string(lnc_seven_m); strmat(7,7) = to_string(mrna_seven_m); strmat(7,9) = to_string(three_seven_m); strmat(7,11) = to_string(five_seven_m); strmat(7,13) = to_string(cds_seven_m); strmat(7,15) = to_string(exon_seven_m); strmat(7,17) = to_string(v_gene_seven_m);  strmat(7,19) = to_string(pseudo_seven_m);
    strmat(9,3) = to_string(inter_eight);    strmat(9,5) = to_string(lnc_eight);   strmat(9,7) = to_string(mrna_eight);   strmat(9,9) = to_string(three_eight);   strmat(9,11) = to_string(five_eight);   strmat(9,13) = to_string(cds_eight);   strmat(9,15) = to_string(exon_eight);   strmat(9,17) = to_string(v_gene_eight);    strmat(9,19) = to_string(pseudo_eight);


    strmat.setMatrix_mirna_id(mirna_id);
    strmat.setMatrix_allhits(all_hits);

    /// calling the printToFile function
    strmat.printMatrixToFile(output_path, true);

    /// For the TOTAL statistics part
    stats.total_hits += all_hits;
    stats.total_8mer += (inter_eight + lnc_eight + mrna_eight + three_eight + five_eight + cds_eight + exon_eight + v_gene_eight + pseudo_eight);
    stats.total_threeprime += (three_six + three_seven_a + three_seven_m + three_eight);
    stats.total_8mer_threeprime += (three_eight);
    stats.total_6mer += (inter_six + lnc_six + mrna_six + three_six + five_six + cds_six + exon_six + v_gene_six + pseudo_six);
    stats.total_7merA1 += (inter_seven_a + lnc_seven_a + mrna_seven_a + three_seven_a + five_seven_a + cds_seven_a + exon_seven_a + v_gene_seven_a + pseudo_seven_a);
    stats.total_7merm8 += (inter_seven_m + lnc_seven_m + mrna_seven_m + three_seven_m + five_seven_m + cds_seven_m + exon_seven_m + v_gene_seven_m + pseudo_seven_m);
    stats.inter_six += inter_six; stats.inter_seven_a += inter_seven_a; stats.inter_seven_m += inter_seven_m; stats.inter_eight += inter_eight;
    stats.lnc_six += lnc_six; stats.lnc_seven_a += lnc_seven_a; stats.lnc_seven_m += lnc_seven_m; stats.lnc_eight += lnc_eight;
    stats.mrna_six += mrna_six; stats.mrna_seven_a += mrna_seven_a; stats.mrna_seven_m += mrna_seven_m; stats.mrna_eight += mrna_eight;
    stats.three_six += three_six; stats.three_seven_a += three_seven_a; stats.three_seven_m += three_seven_m; stats.three_eight += three_eight;
    stats.five_six += five_six; stats.five_seven_a += five_seven_a; stats.five_seven_m += five_seven_m; stats.five_eight += five_eight;
    stats.cds_six += cds_six; stats.cds_seven_a += cds_seven_a; stats.cds_seven_m += cds_seven_m; stats.cds_eight += cds_eight;
    stats.exon_six += exon_six; stats.exon_seven_a += exon_seven_a; stats.exon_seven_m += exon_seven_m; stats.exon_eight += exon_eight;
    stats.v_gene_six += v_gene_six; stats.v_gene_seven_a += v_gene_seven_a; stats.v_gene_seven_m += v_gene_seven_m; stats.v_gene_eight += v_gene_eight;
    stats.pseudo_six += pseudo_six; stats.pseudo_seven_a += pseudo_seven_a; stats.pseudo_seven_m += pseudo_seven_m; stats.pseudo_eight += pseudo_eight;

    /// For the LOCAL statistics part (only this miRNA)
    size_t all_8mer = inter_eight + lnc_eight + mrna_eight + three_eight + five_eight + cds_eight + exon_eight + v_gene_eight + pseudo_eight;
    size_t all_threeprime = three_six + three_seven_a + three_seven_m + three_eight;

    double percent_eight = ceil(((double)(all_8mer*100) / all_hits) * 100) / 100;                           // Out of x hits, y are 8mer
    double percent_threeprime = ceil(((double)(all_threeprime*100) / all_hits) * 100) / 100;                // Out of x hits, y are 3'UTR
    double percent_8mer_of_threeprime = ceil(((double)(three_eight*100) / all_threeprime) * 100) / 100;     // Out of x 3'UTR hits, y are 8mer
    double percent_threeprime_of_8mer = ceil(((double)(three_eight*100) / all_8mer) * 100) / 100;           // Out of x 8mer hits, y are 3'UTR

    ofstream output;
    output.open(output_path, std::ios_base::app);

    output << endl << "Statistics\n" << "----------\n";
    output << "Out of " << all_hits << " hits, there are:\n" << "- " << all_8mer << " 8mer hits (" << percent_eight << " %)\n";
    output << "- " << all_threeprime << " 3'-UTR hits (" << percent_threeprime << " %)\n";
    output << "Out of the " << all_threeprime << " 3'-UTR hits, there are " << three_eight << " 8mer hits (" << percent_8mer_of_threeprime << " %)\n";
    output << "Out of the " << all_8mer << " 8mer hits, there are " << three_eight << " hits in 3'-UTR (" << percent_threeprime_of_8mer << " %)\n";

    output.close();

}

/**
 * Prints the alignment in the command line
 *
 * @param mirna The given miRNA string
 * @param dna   The given DNA string
 */
void printAlignment(string& mirna, string& dna){

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
        else {
            between += ' ';
        }
    }

    cout << "Alignment:\n" << mirna << "\n" << between << "\n" << dna << endl;

}



/**
 * Writes output_results into file.
 * Results will be sorted as follows:
 *      region (i.e. intergenic, 5'/3' UTR, exon etc)
 *      binding site motif (6mer, 7mer-m8, 7mer-A1, 8mer)
 * mismatches and the position of the alignment will also be shown
 *
 * @param alignments    The map containing all the alignments and information needed
 * @param output_path   The path to the output file
 * @param mirna_id  The mIRNA ID
 * @param best_alignments_list  A sequence alignment vector keeping track of the best alignments
 */
void output_results(std::map<size_t, SeqAlignment>& alignments, const char* output_path, std::string& mirna_id, std::vector<SeqAlignment>& best_alignments_list){

    std::string loc_on_genome, bsm;

    std::string what_happened;

    std::string int_six, int_seven1, int_seven8, int_eight, lnc_six, lnc_seven1, lnc_seven8, lnc_eight, three_six, three_seven1, three_seven8, three_eight;
    std::string five_six, five_seven1, five_seven8, five_eight, mRNA_six, mRNA_seven1, mRNA_seven8, mRNA_eight, cds_six, cds_seven1, cds_seven8, cds_eight, exon_six, exon_seven1, exon_seven8, exon_eight;
    std::string v_gene_six, v_gene_seven1, v_gene_seven8, v_gene_eight, pseudo_six, pseudo_seven1, pseudo_seven8, pseudo_eight;

    float best_score = 0;
    vector<SeqAlignment> best_alignments;

    for(auto it: alignments) {

        loc_on_genome = it.second.getLocOnGenome();
        bsm = it.second.getBSMotif();

        if (loc_on_genome == "intergenic" && bsm == "6mer") {

            int_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                           to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                           to_string(it.second.getScore()) + "\n");
            int_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            int_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            int_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "intergenic" && bsm == "7mer-A1") {
            int_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            int_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            int_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            int_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "intergenic" && bsm == "7mer-m8") {
            int_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            int_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            int_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            int_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "intergenic" && bsm == "8mer") {
            int_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                             to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                             to_string(it.second.getScore()) + "\n");
            int_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            int_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            int_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "lnc_RNA" && bsm == "6mer") {
            lnc_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                           to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                           to_string(it.second.getScore()) + "\n");
            lnc_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            lnc_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            lnc_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "lnc_RNA" && bsm == "7mer-A1") {
            lnc_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            lnc_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            lnc_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            lnc_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "lnc_RNA" && bsm == "7mer-m8") {
            lnc_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            lnc_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            lnc_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            lnc_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "lnc_RNA" && bsm == "8mer") {
            lnc_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                             to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                             to_string(it.second.getScore()) + "\n");
            lnc_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            lnc_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            lnc_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "three_prime_UTR" && bsm == "6mer") {
            three_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                             to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                             to_string(it.second.getScore()) + "\n");
            three_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            three_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            three_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "three_prime_UTR" && bsm == "7mer-A1") {
            three_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                                to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                                to_string(it.second.getScore()) + "\n");
            three_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            three_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            three_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "three_prime_UTR" && bsm == "7mer-m8") {
            three_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                                to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                                to_string(it.second.getScore()) + "\n");
            three_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            three_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            three_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "three_prime_UTR" && bsm == "8mer") {
            three_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                               to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                               to_string(it.second.getScore()) + "\n");
            three_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            three_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            three_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "five_prime_UTR" && bsm == "6mer") {
            five_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                            to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                            to_string(it.second.getScore()) + "\n");
            five_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            five_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            five_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "five_prime_UTR" && bsm == "7mer-A1") {
            five_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                               to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                               to_string(it.second.getScore()) + "\n");
            five_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            five_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            five_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "five_prime_UTR" && bsm == "7mer-m8") {
            five_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                               to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                               to_string(it.second.getScore()) + "\n");
            five_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            five_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            five_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "five_prime_UTR" && bsm == "8mer") {
            five_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            five_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            five_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            five_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "mRNA" && bsm == "6mer") {
            mRNA_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                            to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                            to_string(it.second.getScore()) + "\n");
            mRNA_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            mRNA_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            mRNA_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "mRNA" && bsm == "7mer-A1") {
            mRNA_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                               to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                               to_string(it.second.getScore()) + "\n");
            mRNA_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            mRNA_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            mRNA_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "mRNA" && bsm == "7mer-m8") {
            mRNA_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                               to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                               to_string(it.second.getScore()) + "\n");
            mRNA_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            mRNA_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            mRNA_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "mRNA" && bsm == "8mer") {
            mRNA_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            mRNA_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            mRNA_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            mRNA_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "CDS" && bsm == "6mer") {
            cds_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                           to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                           to_string(it.second.getScore()) + "\n");
            cds_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            cds_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            cds_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "CDS" && bsm == "7mer-A1") {
            cds_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            cds_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            cds_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            cds_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "CDS" && bsm == "7mer-m8") {
            cds_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            cds_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            cds_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            cds_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "CDS" && bsm == "8mer") {
            cds_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                             to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                             to_string(it.second.getScore()) + "\n");
            cds_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            cds_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            cds_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "exon" && bsm == "6mer") {
            exon_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                            to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                            to_string(it.second.getScore()) + "\n");
            exon_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            exon_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            exon_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "exon" && bsm == "7mer-A1") {
            exon_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                               to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                               to_string(it.second.getScore()) + "\n");
            exon_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            exon_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            exon_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "exon" && bsm == "7mer-m8") {
            exon_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                               to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                               to_string(it.second.getScore()) + "\n");
            exon_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            exon_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            exon_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "exon" && bsm == "8mer") {
            exon_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            exon_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            exon_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            exon_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "V_gene_segment" && bsm == "6mer") {
            v_gene_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                           to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                           to_string(it.second.getScore()) + "\n");
            v_gene_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            v_gene_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            v_gene_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "V_gene_segment" && bsm == "7mer-A1") {
            v_gene_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            v_gene_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            v_gene_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            v_gene_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "V_gene_segment" && bsm == "7mer-m8") {
            v_gene_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            v_gene_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            v_gene_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            v_gene_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "V_gene_segment" && bsm == "8mer") {
            v_gene_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                             to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                             to_string(it.second.getScore()) + "\n");
            v_gene_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            v_gene_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            v_gene_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "pseudogenic_transcript" && bsm == "6mer") {
            pseudo_six.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                              to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                              to_string(it.second.getScore()) + "\n");
            pseudo_six.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            pseudo_six.append("\t    " + it.second.getBetweenSeq() + "\n");
            pseudo_six.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "pseudogenic_transcript" && bsm == "7mer-A1") {
            pseudo_seven1.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                                 to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                                 to_string(it.second.getScore()) + "\n");
            pseudo_seven1.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            pseudo_seven1.append("\t    " + it.second.getBetweenSeq() + "\n");
            pseudo_seven1.append("\t5'  " + it.second.getSeq2() + "  3'\n");if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "pseudogenic_transcript" && bsm == "7mer-m8") {
            pseudo_seven8.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                                 to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                                 to_string(it.second.getScore()) + "\n");
            pseudo_seven8.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            pseudo_seven8.append("\t    " + it.second.getBetweenSeq() + "\n");
            pseudo_seven8.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
        else if (loc_on_genome == "pseudogenic_transcript" && bsm == "8mer") {
            pseudo_eight.append("Chr.: " + it.second.getChrID() + ", Pos.: " + to_string(it.second.getStart()) + " - " + to_string(it.second.getEnd()) + ", Mismatches: " +
                                to_string(it.second.getMismatches()) + ", " + it.second.getComplementarySite() + " match, Score: " +
                                to_string(it.second.getScore()) + "\n");
            pseudo_eight.append("\t3'  " + it.second.getSeq1() + "  5'\n");
            pseudo_eight.append("\t    " + it.second.getBetweenSeq() + "\n");
            pseudo_eight.append("\t5'  " + it.second.getSeq2() + "  3'\n");
            if(it.second.getScore() > best_score){
                best_score = it.second.getScore();
                best_alignments.clear();
                best_alignments.push_back(it.second);
            }
            else if(it.second.getScore() == best_score){
                best_alignments.push_back(it.second);
            }
        }
    }

    ofstream output;
    output.open(output_path);
    output << mirna_id << endl << endl;
    /// Best alignments/scores
    output << "# Best hits:\n";
    for(auto best: best_alignments){
        output << best.getBSMotif() << " / " << best.getLocOnGenome() << endl << "Chr.: " << best.getChrID() << ", Pos.: " << best.getStart() << " - " << best.getEnd() << ", Mismatches: " << best.getMismatches()
               << ", " << best.getComplementarySite() << " match, Score: " << best.getScore() << endl;
        output << "\t3'  " << best.getSeq1() << "  5'\n" << "\t    " << best.getBetweenSeq() << endl << "\t5'  " << best.getSeq2() << "  3'\n";

        best.setmiRNAID(mirna_id);
        best_alignments_list.push_back(best);
    }
    ///
    output << "#\n\n";
    output << "6mer:\n-----\n\n" << "Intergenic Region:\n" << int_six << "\nlnc_RNA Region:\n" << lnc_six << "\n3'-UTR Region:\n" << three_six << "\n5'-UTR Region:\n" << five_six <<
           "\nmRNA Region:\n" << mRNA_six << "\nCDS Region:\n" << cds_six << "\nExon Region:\n" << exon_six << "\nV_gene Region:\n" << v_gene_six << "\nPseudo Region:\n" << pseudo_six << endl << endl;
    output << "7mer-A1:\n-------\n\n" << "Intergenic Region:\n" << int_seven1 << "\nlnc_RNA Region:\n" << lnc_seven1 << "\n3'-UTR Region:\n" << three_seven1 << "\n5'-UTR Region:\n" << five_seven1 <<
           "\nmRNA Region:\n" << mRNA_seven1 << "\nCDS Region:\n" << cds_seven1 << "\nExon Region:\n" << exon_seven1 << "\nV_gene Region:\n" << v_gene_seven1 << "\nPseudo Rgeion:\n" << pseudo_seven1 << endl << endl;
    output << "7mer-m8:\n-------\n\n" << "Intergenic Region:\n" << int_seven8 << "\nlnc_RNA Region:\n" << lnc_seven8 << "\n3'-UTR Region:\n" << three_seven8 << "\n5'-UTR Region:\n" << five_seven8 <<
           "\nmRNA Region:\n" << mRNA_seven8 << "\nCDS Region:\n" << cds_seven8 << "\nExon Region:\n" << exon_seven8 << "\nV_gene Region:\n" << v_gene_seven8 << "\nPseudo Region:\n" << pseudo_seven8 << endl << endl;
    output << "8mer:\n-----\n\n" << "Intergenic Region:\n" << int_eight << "\nlnc_RNA Region:\n" << lnc_eight << "\n3'-UTR Region:\n" << three_eight << "\n5'-UTR Region:\n" << five_eight <<
           "\nmRNA Region:\n" << mRNA_eight << "\nCDS Region:\n" << cds_eight << "\nExon Region:\n" << exon_eight << "\nV_gene Region:\n" << v_gene_eight << "\nPseudo Region:\n" << endl << endl;

    output.close();

}

/**
 * Prints the overall statistics table
 *
 * @param stats The statistics object containig all the information needed
 * @param miRNA_size    The number of miRNAs this query
 * @param output_path   The path to the output file
 */
void overall_stats_output(Statistics stats, size_t miRNA_size, const char* output_path){

    ofstream stats_output;
    stats_output.open(output_path);

    stats_output << "Overall Statistics of last query\n--------------------------------\nWe had " << miRNA_size << " miRNA(s)\n--------------------------------\n\n";
    stats_output.close();

    StringMatrix strmat = StringMatrix(11,21," ");      //13,19
    strmat(1,3) = "intergenic"; strmat(1,5) = "lnc_RNA"; strmat(1,7) = "mRNA"; strmat(1,9) = "3'UTR"; strmat(1,11) = "5'UTR"; strmat(1,13) = "CDS";
    strmat(1,15) = "exon"; strmat(1,17) = "V_gene"; strmat(1,19) = "pseudo";
    strmat(3,1) = "6mer"; strmat(5,1) = "7mer-A1"; strmat(7,1) = "7mer-m8"; strmat(9,1) = "8mer";

    /// filling the "output_matrix"
    strmat(3,3) = to_string(stats.inter_six);      strmat(3,5) = to_string(stats.lnc_six) ;    strmat(3,7) = to_string(stats.mrna_six);     strmat(3,9) = to_string(stats.three_six);     strmat(3,11) = to_string(stats.five_six);     strmat(3,13) = to_string(stats.cds_six);     strmat(3,15) = to_string(stats.exon_six);     strmat(3,17) = to_string(stats.v_gene_six);      strmat(3,19) = to_string(stats.pseudo_six);
    strmat(5,3) = to_string(stats.inter_seven_a) ; strmat(5,5) = to_string(stats.lnc_seven_a); strmat(5,7) = to_string(stats.mrna_seven_a); strmat(5,9) = to_string(stats.three_seven_a); strmat(5,11) = to_string(stats.five_seven_a); strmat(5,13) = to_string(stats.cds_seven_a); strmat(5,15) = to_string(stats.exon_seven_a); strmat(5,17) = to_string(stats.v_gene_seven_a);  strmat(5,19) = to_string(stats.pseudo_seven_a);
    strmat(7,3) = to_string(stats.inter_seven_m);  strmat(7,5) = to_string(stats.lnc_seven_m); strmat(7,7) = to_string(stats.mrna_seven_m); strmat(7,9) = to_string(stats.three_seven_m); strmat(7,11) = to_string(stats.five_seven_m); strmat(7,13) = to_string(stats.cds_seven_m); strmat(7,15) = to_string(stats.exon_seven_m); strmat(7,17) = to_string(stats.v_gene_seven_m);  strmat(7,19) = to_string(stats.pseudo_seven_m);
    strmat(9,3) = to_string(stats.inter_eight);    strmat(9,5) = to_string(stats.lnc_eight);   strmat(9,7) = to_string(stats.mrna_eight);   strmat(9,9) = to_string(stats.three_eight);   strmat(9,11) = to_string(stats.five_eight);   strmat(9,13) = to_string(stats.cds_eight);   strmat(9,15) = to_string(stats.exon_eight);   strmat(9,17) = to_string(stats.v_gene_eight);    strmat(9,19) = to_string(stats.pseudo_eight);

    strmat.setMatrix_allhits(stats.total_hits);

    /// printing the overall normal string matrix
    strmat.printMatrixToFile(output_path, false);

    /// now calculating the %
    double percent_six = ceil(((double)(stats.total_6mer*100) / stats.total_hits) * 100) / 100;
    double percent_seven_a = ceil(((double)(stats.total_7merA1*100) / stats.total_hits) * 100) / 100;
    double percent_seven_m = ceil(((double)(stats.total_7merm8*100) / stats.total_hits) * 100) / 100;
    double percent_eight = ceil(((double)(stats.total_8mer*100) / stats.total_hits) * 100) / 100;

    /// set precision
    std::string per_six = to_string(percent_six);
    per_six.erase(per_six.find_last_not_of('0') + 1, std::string::npos); per_six.append(" %");
    std::string per_seven_a = to_string(percent_seven_a);
    per_seven_a.erase(per_seven_a.find_last_not_of('0') + 1, std::string::npos); per_seven_a.append(" %");
    std::string per_seven_m = to_string(percent_seven_m);
    per_seven_m.erase(per_seven_m.find_last_not_of('0') + 1, std::string::npos); per_seven_m.append(" %");
    std::string per_eight = to_string(percent_eight);
    per_eight.erase(per_eight.find_last_not_of('0') + 1, std::string::npos); per_eight.append(" %");

    /// now generating and filling the stat_table
    StringMatrix stat_mat = StringMatrix(9, 7, " ");
    stat_mat(1,1) = "6mer"; stat_mat(1, 3) = to_string(stats.total_6mer); stat_mat(1,5) = per_six;
    stat_mat(3,1) = "7mer-A1"; stat_mat(3,3) = to_string(stats.total_7merA1); stat_mat(3,5) = per_seven_a;
    stat_mat(5,1) = "7mer-m8"; stat_mat(5,3) = to_string(stats.total_7merm8); stat_mat(5,5) = per_seven_m;
    stat_mat(7,1) = "8mer"; stat_mat(7,3) = to_string(stats.total_8mer); stat_mat(7,5) = per_eight;

    /// print that matrix
    stat_mat.print_statmat_toFile(output_path);

    /// now adding the 3 specific cases

    double percent_threeprime = ceil(((double)(stats.total_threeprime*100) / stats.total_hits) * 100) / 100;
    double percent_8mer_of_threeprime = ceil(((double)(stats.total_8mer_threeprime*100) / stats.total_threeprime) * 100) / 100;
    double percent_threeprime_of_8mer = ceil(((double)(stats.total_8mer_threeprime*100) / stats.total_8mer) * 100) / 100;

    stats_output.open(output_path, std::ios_base::app);

    stats_output << "\n> 3'-UTR: " << stats.total_threeprime << " (" << percent_threeprime << " %)" << endl;
    stats_output << "\t> Out of " << stats.total_threeprime << " hits, " << stats.total_8mer_threeprime << " of these were 8mer (" << percent_8mer_of_threeprime << " %)\n";
    stats_output << "\t> " << percent_threeprime_of_8mer << " % of 8mer hits were in 3'-UTR\n";

    stats_output.close();
}

/**
 * Processes an alignment file
 *
 * @param input_path    Path to the alignment file
 * @param alignment_map A sequence alignment map to save all the alignments of the file in
 * @param best_alignments   A vector containing all the best alignments
 */
void process_alignment_input(const string& input_path, std::map<size_t, SeqAlignment>& alignment_map, vector<SeqAlignment>& best_alignments){

    fstream input;
    input.open(input_path);

    if(!input.is_open()){
        cerr << "Could not open the given file!" << endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    size_t hash_count = 0;
    bool first_hash = false;
    bool sec_hash = false;
    bool first = true;
    size_t key = 0;

    std::string mirna_id;

    // for SeqAlignment
    std::string seq1, seq2, between, chr_id, bsm, loc, compl_site;
    size_t seq_start, seq_end, mismatches;
    float score;
    //
    std::string one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirt, fourt;

    while(getline(input, line)){

        if(first){
            istringstream iss(line);
            iss >> mirna_id;
            first = false;
        }
        else {

            if (line.substr(0, 1) == "#") {
                hash_count++;
                if (hash_count == 1) {
                    first_hash = true;
                    continue;
                }
                if(hash_count == 2){
                    first_hash = false;
                    sec_hash = true;
                }
            }

            if (first_hash) {
                istringstream iss(line);
                iss >> one;
                if (one == "6mer" || one == "7mer-A1" || one == "7mer-m8" || one == "8mer"){
                    iss >> two >> three;
                    bsm = one;
                    loc = three;
                } else if (one == "Chr.:") {
                    iss >> two >> three >> four >> five >> six >> seven >> eight >> nine >> ten >> eleven >> twelve;
                    two.pop_back(); six.pop_back(); eight.pop_back();
                    chr_id = two;
                    seq_start = stoi(four);
                    seq_end = stoi(six);
                    mismatches = stoi(eight);
                    if (ten == "match,") {
                        compl_site = "seed";
                        score = stof(twelve);
                    } else {
                        iss >> thirt >> fourt;
                        compl_site = "seed + 3'";
                        score = stof(fourt);
                    }
                } else if(one == "3'"){
                    iss >> two;
                    seq1 = two;
                } else if(one == "5'"){
                    iss >> two;
                    seq2 = two;
                    SeqAlignment seq_ali = SeqAlignment(seq1, seq2, between, bsm, compl_site, mismatches, score);
                    seq_ali.setStart(seq_start); seq_ali.setEnd(seq_end);
                    seq_ali.setLocOnGenome(loc); seq_ali.setChrID(chr_id);
                    seq_ali.setmiRNAID(mirna_id);
                    best_alignments.push_back(seq_ali);
                } else {
                    between = line;
                }
            }
            if(sec_hash) {
                istringstream iss(line);
                iss >> one >> two;
                if (one == "6mer:" || one == "7mer-A1:" || one == "7mer-m8:" || one == "8mer:") {
                    one.pop_back();
                    bsm = one;
                }
                else if (two == "Region:") {
                    if (one == "Intergenic") { loc = "intergenic"; }
                    else if (one == "lnc_RNA") { loc = "lnc_RNA"; }
                    else if (one == "3'-UTR") { loc = "three_prime_UTR"; }
                    else if (one == "5'-UTR") { loc = "five_prime_UTR"; }
                    else if (one == "mRNA") { loc = "mRNA"; }
                    else if (one == "CDS") { loc = "CDS"; }
                    else if (one == "Exon") { loc = "exon"; }
                    else if (one == "V_gene") { loc = "V_gene_segment"; }
                    else if (one == "Pseudo") { loc = "pseudogenic_transcript"; }
                } else {
                    if (one == "Chr.:") {
                        iss >> three >> four >> five >> six >> seven >> eight >> nine >> ten >> eleven >> twelve;
                        two.pop_back();
                        six.pop_back();
                        eight.pop_back();
                        chr_id = two;
                        seq_start = stoi(four);
                        seq_end = stoi(six);
                        mismatches = stoi(eight);
                        if (ten == "match,") {
                            compl_site = "seed";
                            score = stof(twelve);
                        } else {
                            iss >> thirt >> fourt;
                            compl_site = "seed + 3'";
                            score = stof(fourt);
                        }
                    } else if (one == "3'") {
                        seq1 = two;
                    } else if (one == "5'") {
                        seq2 = two;
                        SeqAlignment seq_ali = SeqAlignment(seq1, seq2, between, bsm, compl_site, mismatches, score);
                        seq_ali.setStart(seq_start);
                        seq_ali.setEnd(seq_end);
                        seq_ali.setLocOnGenome(loc);
                        seq_ali.setChrID(chr_id);
                        seq_ali.setmiRNAID(mirna_id);
                        alignment_map.insert(pair<size_t, SeqAlignment>(key, seq_ali));
                        key++;
                        one.clear();
                    } else {
                        between = line;
                    }
                }
            }
        }
    }
}

/**
 * Prints the best alignments into a given output file
 *
 * @param best_alignments   The vector containing all the best alignments
 * @param output_path   The path to the output file
 */
void print_best_alignments(std::vector<SeqAlignment>& best_alignments, const char* output_path){

    float best = 0;

    vector<SeqAlignment> output_alignments;

    for(auto it: best_alignments){
        if(best < it.getScore()){
            best = it.getScore();
            output_alignments.clear();
            output_alignments.push_back(it);
        }
        else if(best == it.getScore()){
            output_alignments.push_back(it);
        }
    }

    ofstream best_alignments_output;
    best_alignments_output.open(output_path);

    best_alignments_output << "Best Alignments:\n\n";

    for(auto it: output_alignments){
        best_alignments_output << it.getmiRNAID() << " | " << it.getBSMotif() << " | " << it.getLocOnGenome() << endl;
        best_alignments_output << it.getChrID() << " | " << it.getStart() << " - " << it.getEnd() << " | " << it.getMismatches() << " MM | " << it.getComplementarySite() << " match | Score: " << it.getScore() << endl;
        best_alignments_output << "\t3'  " << it.getSeq1() << "  5'\n\t    " << it.getBetweenSeq() << "\n\t5'  " << it.getSeq2() << "  3'\n\n";
    }

}

/**
 * Extracts the best alignments out of an alignment file
 *
 * @param input_path    Path to the alignment file
 * @param all_alignments    A vector ontaining all the alignments
 */
void extract_best_alignments(const string& input_path, std::vector<SeqAlignment>& all_alignments){

    fstream input;
    input.open(input_path);

    if(!input.is_open()){
        cerr << "Could not open the given file!" << endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    size_t hash_count = 0;
    bool first = true;
    size_t key = 0;

    std::string mirna_id;

    // for SeqAlignment
    std::string seq1, seq2, between, chr_id, bsm, loc, compl_site;
    size_t seq_start, seq_end, mismatches;
    float score;
    //
    std::string one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirt, fourt;

    while(getline(input, line)){

        if(first){
            istringstream iss(line);
            iss >> mirna_id;
            first = false;
        }
        else {
            istringstream iss(line);
            iss >> one >> two;
            if (one == "6mer:" || one == "7mer-A1:" || one == "7mer-m8:" || one == "8mer:") {
                one.pop_back();
                bsm = one;
            }
            else if (two == "Region:") {
                if (one == "Intergenic") { loc = "intergenic"; }
                else if (one == "lnc_RNA") { loc = "lnc_RNA"; }
                else if (one == "3'-UTR") { loc = "three_prime_UTR"; }
                else if (one == "5'-UTR") { loc = "five_prime_UTR"; }
                else if (one == "mRNA") { loc = "mRNA"; }
                else if (one == "CDS") { loc = "CDS"; }
                else if (one == "Exon") { loc = "exon"; }
                else if (one == "V_gene") { loc = "V_gene_segment"; }
                else if (one == "Pseudo") { loc = "pseudogenic_transcript"; }
            } else {
                if (one == "Chr.:") {
                    iss >> three >> four >> five >> six >> seven >> eight >> nine >> ten >> eleven >> twelve;
                    two.pop_back();
                    six.pop_back();
                    eight.pop_back();
                    chr_id = two;
                    seq_start = stoi(four);
                    seq_end = stoi(six);
                    mismatches = stoi(eight);
                    if (ten == "match,") {
                        compl_site = "seed";
                        score = stof(twelve);
                    } else {
                        iss >> thirt >> fourt;
                        compl_site = "seed + 3'";
                        score = stof(fourt);
                    }
                } else if (one == "3'") {
                    seq1 = two;
                } else if (one == "5'") {
                    seq2 = two;
                    SeqAlignment seq_ali = SeqAlignment(seq1, seq2, between, bsm, compl_site, mismatches, score);
                    seq_ali.setStart(seq_start);
                    seq_ali.setEnd(seq_end);
                    seq_ali.setLocOnGenome(loc);
                    seq_ali.setChrID(chr_id);
                    seq_ali.setmiRNAID(mirna_id);
                    all_alignments.push_back(seq_ali);
                    key++;
                    one.clear();
                } else {
                    between = line;
                }
            }
        }
    }
}

/**
 * FInds the best alignments based on scores
 *
 * @param all_alignments    A vector containing all the alignments
 * @param best_alignments   A vector containing the best alignments
 */
void get_best_alignment(std::vector<SeqAlignment>& all_alignments, std::vector<SeqAlignment>& best_alignments){

    float best = 0;
    vector<SeqAlignment> temp_best_align;

    for(auto it: all_alignments){
        if(best < it.getScore() && it.getScore()){
            best = it.getScore();
            temp_best_align.clear();
            temp_best_align.push_back(it);
        }
        else if(best == it.getScore()){
            temp_best_align.push_back(it);
        }
    }

    for(const auto& it: temp_best_align){
        best_alignments.push_back(it);
    }

}

#endif //MIRNA_OUTPUT_H
