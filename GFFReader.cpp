//#include "json.hpp"
#include "GFFReader.h"

#include <iostream>
#include <random>
#include <string>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <direct.h>
#include <vector>

using namespace std;

/**
 * processes the .gff3 - file and generates the IntervalTree
 *
 * @param node_map string = key = chromosome-id and
 *                 ITNode* = root of that chromosome (every node represents a gene, every gene has its own gene_root_node with every location)
 */
void GFFReader::processGFFFile(std::map<std::string, ITNode*> &node_map) {

    std::string line, chr_id;

    std::string trash1, trash2, trash3, location;
    int start, end;

    size_t pos;

    bool chr_found = false;
    bool gene_found = false;
    bool find_next = false;
    bool first = true;

    std::vector<Interval> ints;

    ITNode *root = nullptr;

    ITNode *gene_root;
    Interval gene_i;
    std::vector<Gene> gene_nodes;
    Gene gene;

    Interval i;

    ofstream temp_gff;

    fstream gff_file;
    gff_file.open(this->gff_file_path);

    if(!gff_file.is_open()){
        cerr << "Could not open the given .gff3 file!" << endl;
        exit(EXIT_FAILURE);
    }

    while(getline(gff_file, line)){

        if(first) {

            pos = line.find("ID=chromosome");
            if (pos != std::string::npos) {

                chr_id = line.substr(0, line.find('\t'));
                chr_found = true;
                first = false;

            }
        }
        else{

            if(gene_found){

                if(line.substr(0,3) == "###"){
                    gene_found = false;

                    shuffle(ints.begin(), ints.end(), std::mt19937(std::random_device()()));
                    gene = {gene_i, NULL};

                    unsigned int n = ints.size();
                    for(int x = 0; x < n; x++){
                        gene.gene_list_node = insert(gene.gene_list_node, ints[x]);
                    }

                    gene_nodes.push_back(gene);
                    ints.clear();

                    continue;
                }
                else if(line.substr(0, line.find('\t')) != chr_id){
                    gene_found = false;
                }
                else {
                    stringstream iss(line);
                    iss >> chr_id >> trash1 >> location >> start >> end;

                    if(i.low == start && i.high == end){

                        if(!ints.empty()) {
                            ints.back().loc = location;
                            continue;
                        }

                    }

                    i = {start, end, location};

                    ints.push_back(i);

                    continue;
                }
            }

            if(chr_found){

                // we need to check if we're still on the same chromosome or if we "enter" another one
                if (line.substr(0, line.find('\t')) == chr_id) {
                    // we are still on the same chromosome

                    // since we are still on the same chr, we can look for genes
                    pos = line.find("ID=gene");
                    if (pos != std::string::npos) {
                        //we found the beginning of a gene, we can go to the next line
                        //and search for the next separator

                        gene_found = true;

                        std::stringstream iss2(line);
                        iss2 >> chr_id >> trash1 >> trash2 >> start >> end;

                        gene_i = {start, end};

                    }

                } else if (line.substr(0, 3) == "###") {
                    // separator found -> check next line
                } else {
                    // no separator, chr_id changed
                    // on the next chromosome
                    chr_found = false;
                    find_next = true;

                }
            }

            if(find_next){

                pos = line.find("ID=chromosome");
                if(pos != std::string::npos){

                    shuffle(gene_nodes.begin(), gene_nodes.end(), std::mt19937(std::random_device()()));

                    unsigned int n = gene_nodes.size();
                    for(int x = 0; x < n; x++){
                        root = insert_with_node(root, gene_nodes[x].gene_interval, gene_nodes[x].gene_list_node);

                    }

                    node_map.insert(std::pair<std::string, ITNode*>(chr_id, root));

                    chr_id = line.substr(0, line.find('\t'));

                    find_next = false;
                    chr_found = true;

                    gene_nodes.clear();
                    root = NULL;
                }
            }
        }
    }

    shuffle(gene_nodes.begin(), gene_nodes.end(), std::mt19937(std::random_device()()));

    unsigned int n = gene_nodes.size();
    for(int x = 0; x < n; x++){
        root = insert_with_node(root, gene_nodes[x].gene_interval, gene_nodes[x].gene_list_node);

    }

    node_map.insert(std::pair<std::string, ITNode*>(chr_id, root));

}