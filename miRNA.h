#ifndef MIRNA_MIRNA_H
#define MIRNA_MIRNA_H

#include "InvalidCharacter.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <vector>
#include <algorithm>


class miRNA{

private:

    std::string miRNA_sequence;      // original sequence
    std::string changed_miRNA_seq;   // sequence subject to change
    std::string miRNA_ID;
    std::string reversed_miRNA{};

    std::string token;
    size_t overhang{};

public:

  explicit miRNA(const char* miRNA_file_path);
  miRNA(std::string& seq, std::string& id);
  miRNA();

  /**
   * Getter and Setter
   */
  std::string &get_miRNASequence();
  const std::string &get_miRNAId() const;
  std::string& get_changed_miRNASequence();
  std::string& get_reversed_miRNA();
  std::string& getToken();
  size_t getOverhang();

  void set_miRNASequence(std::string& seq);
  void set_miRNAId(std::string& id);
  void set_changed_miRNASequence(std::string& seq);
  void set_miRNAToken(size_t token_size);
  void set_miRNAOverhang(size_t oh);

  /**
   * Generates complementary strand-sequences
   */
  void generate_complementary_RNA(std::string& seq);
  void generate_complementary_DNA(std::string& seq);

  /**
   * converts a sequence to either DNA or RNA
   */
  void convert_to_RNA(std::string& seq);
  void convert_to_DNA(std::string& seq);


};

#endif //MIRNA_MIRNA_H