#ifndef MIRNA_GFFREADER_H
#define MIRNA_GFFREADER_H

#include "json.hpp"
#include "SeqAlignment.h"
#include "intervallTree.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <direct.h>
#include <random>

class GFFReader
{

private:

    const char* gff_file_path;
    std::string organism;

public:

    /**
     * constructor
     */
    GFFReader(const char* path, std::string& organism) {
        this->gff_file_path = path;
        this->organism = organism;
    }

    /**
     * processes the .gff3 - file and generates the IntervalTree
     *
     * @param node_map  String = key = chromosome-id and
     *                  ITNode* = root of that chromosome (every node represents a gene, every gene has its own gene_root_node with every location)
     */
    void processGFFFile(std::map<std::string, ITNode*> &node_map);

};

#endif //MIRNA_GFFREADER_H
