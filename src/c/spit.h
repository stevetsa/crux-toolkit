#ifndef SPIT_H
#define SPIT_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <math.h>
#include "match_collection.h"
#include "utils.h"
#include "match.h"
#include "carp.h"
#include "objects.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "protein.h"
#include "database.h"

// A map from peptides to their scores
//typedef std::map<std::string, FLOAT_T> PeptideScore;


int main(int argc, char** argv);
//int spit_main(int argc, char** argv);

 /*
bool get_txt_matches(
	PeptideScore&,
	ProteinPeptides&,
	char* psm_folder
	);
// 
bool get_database_matches(
	PeptideScore&,
	ProteinPeptides&,
	char* psm_folder,
	char* database
 	);

//
void get_spit_scores(
	ProteinScores&,
	const PeptideScore&,
	const ProteinPeptides&
	);

//
bool print_spit_scores (
	const ProteinScores&,
	char* file 
	);
*/
// old
/*
class SpitMatch {
  public:
   SpitMatch(const std::string& peptide_sequence,
	     const std::string& parent_proteins,
	     FLOAT_T score)
   :peptide_sequence_(peptide_sequence),
   parent_proteins_(parent_proteins),
   score_(score)
   {}
   
   SpitMatch(const SpitMatch& that);
   
  const std::string& get_peptide_sequence() const {
    return peptide_sequence_;
  }
  const std::string& get_parent_proteins() const {
    return parent_proteins_;
  }
  FLOAT_T get_score() const {
    return score;
  }

  private:
  std::string peptide_sequence_;
  std::string parent_proteins_; // string of all parent proteins sep. by commas
  FLOAT_T score_; // best pvalue score for peptide
};
*/
#endif
