/**
 * \file q-ranker.cpp
 */
/*
 * AUTHOR: Barbara Frewen and Marina Spivak
 * CREATE DATE: November 25, 2008
 * DESCRIPTION: Copied from match_analysis.c with only the percolator
 *         functionality kept.
 *         Given as input a directory containing binary psm files and
 *         a protein database, run q-ranker and return a txt file
 *         with results.
 *
 *         Handles at most 4 files (target and decoy).  Looks for .csm
 *         files in the input directory and for corresponding
 *         -decoy[123].csm files.  Multiple target files in the given
 *         directory are concatinated together and presumed to be
 *         non-overlaping parts of the same ms2 file. 
 * 
 * $Revision: 1.1.2.3 $
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "q-ranker.h"
#include "protein.h"
#include "peptide.h"
#include "spectrum.h"
#include "parse_arguments.h" 
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "scorer.h"
#include "match.h"
#include "match_collection.h"
#include "QRankerCInterface.h"
#include "output-files.h"

#include "DelimitedFile.h"

#include "ChargeIndex.h"

using namespace std;

/* 
 * Private function declarations.  Details below
 */


void qcInitiate2(NSet sets, unsigned int numFeatures, 
  std::vector<unsigned int>& numSpectra, char ** featureNames, double pi0);


void run_mpsm_q(
  char* psm_result_folder, 
  char* fasta_file, 
  char* feature_file,
  DelimitedFile& qranker_results); 


void collapseScans(DelimitedFile& matches_in, DelimitedFile& matches_out, const string& column_name);

/**
 * \brief crux-analyze-matches: takes in a directory containing binary
 * psm files and a protein index and analyzes the psms.
 */
int mpsm_qranker_main(int argc, char** argv){


  /* Define command line arguments */
  const char* option_list[] = {
    "version",
    "verbosity",
    "parameter-file",
    "fileroot",
    "feature-file",
    "output-dir",
    "overwrite",
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  const char* argument_list[] = {
    "protein input",
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);
  carp(CARP_INFO,"Initialize run");
  initialize_run(MPSM_QRANKER_COMMAND, argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  /* Get arguments */
  char* psm_dir = get_string_parameter("output-dir");
  char* protein_input_name = get_string_parameter("protein input");
  // TODO (BF oct-22-09): consider adding feature file to OutputFiles
  char* feature_file = get_string_parameter("feature-file");


  if (feature_file != NULL) {
    prefix_fileroot_to_name(&feature_file);
  }


  /* Perform the analysis */
  //MATCH_COLLECTION_T* match_collection = NULL;

  DelimitedFile qranker_results;

  carp(CARP_INFO,"run_mpsm_q");
  run_mpsm_q(psm_dir,
             protein_input_name,
             feature_file,
             qranker_results);

  string output_path = string(psm_dir) + string("/q-ranker.target.txt");

  qranker_results.saveData(output_path);

  //collapse scans to lowest q-value.
  DelimitedFile collapsed_results;
  collapseScans(qranker_results, collapsed_results, "q-ranker q-value");

  output_path = string(psm_dir) + string("/q-ranker.target.collapsed.txt");
  collapsed_results.saveData(output_path);

  // clean up
  free(psm_dir);
  free(protein_input_name);
  free(feature_file);


  carp(CARP_INFO, "crux q-ranker finished.");
  exit(0);

}

//const int number_features = 11;
const char* feature_names[number_features] = {
  "xcorr",
  "max_diff",
  "abs_max_diff",
  "neutral_mass",
  "lnSM",
  "max_pep_len",
  "ncharge1",
  "ncharge2",
  "ncharge3",
  "num_peptides",
  "relative_rtime"
};

const int number_features = sizeof(feature_names) / sizeof(char*);


double getMaxDiff(double a, vector<double>& bs) {

  double ans_max = 0;
  double abs_max = 0;
  for (unsigned int i=0;i<bs.size();i++) {
    double diff = a - bs[i];
    double abs_diff = fabs(diff);
    if (abs_diff > abs_max) {
      abs_max = abs_diff;
      ans_max = diff;
    }
  }
  return ans_max;
}

unsigned int getMaxLength(vector<string> sequences) {
  
  unsigned int ans = 0;
  for (unsigned int i=0;i<sequences.size();i++) {
    if (sequences[i].length() > ans) {
      ans = sequences[i].length();
    }
  }
  return ans;
}


double* get_mpsm_features(DelimitedFile& matches) {


  double* features = (double*)mycalloc(number_features, sizeof(double));
  
  double xcorr_score = matches.getDouble("xcorr score");
  double neutral_mass = matches.getDouble("spectrum neutral mass");
  vector<double> peptide_masses;
  matches.getDoubleVectorFromCell("peptide mass", peptide_masses);
  double max_diff = getMaxDiff(neutral_mass, peptide_masses);
  double abs_max_diff = fabs(max_diff);
  

  int matches_spectrum = matches.getInteger("matches/spectrum");

  
  double lnSM = 0;

  if (matches_spectrum > 0) {
    lnSM = log(matches_spectrum);
  }
  //cout <<"Matches:"<<matches_spectrum<<" lnSM:"<<lnSM<<endl;

  vector<string> peptide_sequences;
  matches.getStringVectorFromCell("sequence", peptide_sequences);
  int max_pep_len = getMaxLength(peptide_sequences); 
    
  int num_peptides = peptide_sequences.size();

  
  string string_charge = matches.getString("charge");
  ChargeIndex charge(string_charge);

  int ncharge1 = charge.numCharge(1);
  int ncharge2 = charge.numCharge(2);
  int ncharge3 = charge.numCharge(3);


  

  double relative_rtime = 0;//matches.getDouble("corr");

  features[0] = xcorr_score;
  features[1] = max_diff;//max_diff;
  features[2] = abs_max_diff;//abs_max_diff;
  features[3] = neutral_mass;//neutral_mass;
  features[4] = lnSM;//lnSM;
  features[5] = max_pep_len;//max_pep_len;
  features[6] = ncharge1;
  features[7] = ncharge2;
  features[8] = ncharge3;
  features[9] = num_peptides;
  features[10] = relative_rtime;


  return features;

}

void registerMatches(DelimitedFile& matches, SetType set_idx, FILE* feature_file) {
  double* features = NULL;
  matches.reset();
  while (matches.hasNext()) {
    //carp(CARP_INFO,"Getting features:%d",set_idx);
    features = get_mpsm_features(matches);
    //carp(CARP_INFO,"Registering features");
    /*
    for (int i=0;i<11;i++) {
      carp(CARP_INFO,"feature[%d]=%lf",i,features[i]);
    }
  */
    if ((int)set_idx == 0) {
       string id = matches.getString("scan")+"-"+matches.getString("charge"); //unique;
     
       qcRegisterPSM(set_idx, (char*)id.c_str(), features);
    } else {
       qcRegisterPSM(set_idx, NULL, features);
    }
    
    if (feature_file != NULL) {
      fprintf(feature_file,"%d\t",matches.getInteger("scan"));
      fprintf(feature_file,"%d\t",(int)set_idx);

      fprintf(feature_file,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",features[0],features[1],features[2],features[3],features[4],features[5],features[6],features[7],features[8],features[9],features[10]);
    }
    //carp(CARP_INFO,"freeing features");
    free(features);
    matches.next();
  }
}


/*  ****************** Subroutines ****************/

/**
 * \brief Analyze matches using the q-ranker algorithm
 * 
 * Runs Marina Spivak's q-ranker algorithm on the PSMs in the psm_result_folder
 * for a search against the sequence database fasta_file. Optionally 
 * puts the q-ranker PSM feature vectors into feature_file, if it is 
 * not NULL.
 * \returns a pointer to a MATCH_COLLECTION_T object
 * \callgraph
 */
void run_mpsm_q(
  char* psm_result_folder, 
  char* fasta_file, 
  char* feature_file,
  DelimitedFile& qranker_results){ 

  DelimitedFile target_search_results;
  vector<DelimitedFile> decoys_search_results;   
  vector<unsigned int> counts;

  carp(CARP_INFO,"Fasta file:%s",fasta_file);

  double* results_q = NULL;
  double* results_score = NULL;
  double pi0 = get_double_parameter("pi0");
  FILE* feature_fh = NULL;
//  int set_idx = 0;
  
  // optional feature_file
  if(feature_file != NULL){  
    if((feature_fh = fopen(feature_file, "w")) == NULL){
      carp(CARP_FATAL, "Problem opening output file %s", feature_file);
      return;
    }
  }

  carp(CARP_DETAILED_DEBUG, "Created feature file");

  //load the results file.
  carp(CARP_INFO,"Loading target");
  string target_path = psm_result_folder+string("/search.target.txt");
  qranker_results.loadData(target_path, true);

  counts.push_back(qranker_results.numRows());

  for (int i=1;i<=2;i++) {
    DelimitedFile decoy_search_results;
    string decoy_path = 
      psm_result_folder + 
      string("/search.decoy-") + 
      DelimitedFile::to_string<int>(i) + 
      string(".txt");
    carp(CARP_INFO,"registering %s",decoy_path.c_str());

    decoy_search_results.loadData(decoy_path);
    counts.push_back(decoy_search_results.numRows());

    decoys_search_results.push_back(decoy_search_results);  
  }

  int num_target_matches = qranker_results.numRows();

    


  carp(CARP_INFO,"Num target matches:%d", num_target_matches);

  carp(CARP_INFO,"Allocating result memory");
  results_q = (double*)mycalloc(num_target_matches, sizeof(double));
  results_score = (double*)mycalloc(num_target_matches, sizeof(double));




  carp(CARP_INFO,"Calling qcInitiate");
  // Call that initiates q-ranker
  qcInitiate2((NSet)3, number_features, counts, (char**)feature_names, pi0);
  //carp(CARP_INFO,"%d",pi0);
  // Call that sets verbosity level
  // 0 is quiet, 2 is default, 5 is more than you want
  if(verbosity < CARP_ERROR){
    qcSetVerbosity(0);
  }    
  else if(verbosity < CARP_INFO){
    qcSetVerbosity(1);
  }
  else{
    qcSetVerbosity(5);
  }
    
  // create iterator, to register each PSM feature to q-ranker.
  carp(CARP_INFO,"Registering targets");
  registerMatches(qranker_results, SetType(0), feature_fh);

  for (int i=1;i<=2;i++) {
    DelimitedFile decoy_search_results;
    string decoy_path = 
      psm_result_folder + 
      string("/search.decoy-") + 
      DelimitedFile::to_string<int>(i) + 
      string(".txt");
    carp(CARP_INFO,"registering %s",decoy_path.c_str());

    decoy_search_results.loadData(decoy_path);
    carp(CARP_INFO,"Number of decoy features:%d", decoy_search_results.numRows());
    registerMatches(decoy_search_results, SetType(i), feature_fh);  
  }

  /***** Q-RANKER run *********/

    carp(CARP_INFO, "got to here");
    
    // Start processing
  qcExecute(); 
  carp(CARP_INFO," Done executing q-ranker");  
  /* Retrieving target scores and qvalues after 
   * processing, the array should be numSpectra long and will be filled in 
   * the same order as the features were inserted */
  carp(CARP_INFO," Getting results");
  qcGetScores(results_score, results_q); 



  carp(CARP_INFO,"Filling in results");
  carp(CARP_INFO,"Number of columns:%d",qranker_results.numCols());
  carp(CARP_INFO,"Number of rows:%d",qranker_results.numRows());
  // fill results for QRANKER_Q_VALUE and QRANKER_SCORE
  int qranker_q_col = qranker_results.findColumn("q-ranker q-value");
  carp(CARP_INFO,"q_col:%d",qranker_q_col);

  int qranker_s_col = qranker_results.findColumn("q-ranker score");
  carp(CARP_INFO,"s_col:%d",qranker_s_col);
 
  carp(CARP_INFO,"Filling in results2");
  for (unsigned int row_idx=0;row_idx<qranker_results.numRows();row_idx++) {
    qranker_results.setValue(qranker_q_col, row_idx, results_q[row_idx]);
    qranker_results.setValue(qranker_s_col, row_idx, results_score[row_idx]);
  }

  carp(CARP_INFO,"Calling cleanup");
  // Function that should be called after processing finished
  //qcCleanUp();

  carp(CARP_INFO,"freeing results arrays");
  free(results_q);
  free(results_score);


  // TODO put free back in. took out because glibc claimed it was corrupted
  // double linked list
  // free_parameters();
}


void getBestScore(DelimitedFile& matches_in, int start_row, int last_row, 
                  const string& column_name, int& best_row, FLOAT_T& best_score) {
  
  best_row = start_row;
  best_score = matches_in.getFloat(column_name.c_str(), best_row);

  for (int idx=start_row+1;idx<=last_row;idx++) {
    FLOAT_T current_score = matches_in.getFloat(column_name.c_str(), idx);
    if (current_score < best_score) {
      best_row = idx;
      best_score = current_score;
    }
  }
}

void collapseScans(DelimitedFile& matches_in, DelimitedFile& matches_out, const string& column_name) {
  //TODO: figure out what to do with duplicates. (Take the least complex?)
  matches_out.clear();
  //matches_out.addColumns(matches_in.getColumnNames());
  vector<string>& column_names = matches_in.getColumnNames();
  for (unsigned int idx=0;idx<column_names.size();idx++) {
    matches_out.addColumn(column_names[idx]);
  }
  
  //make sure scans are together.
  matches_in.sortByIntegerColumn("scan");

  int last_scan = matches_in.getInteger("scan", 0);
  int first_row = 0;
  int best_row = 0;
  FLOAT_T best_score = 0;
  

  for (unsigned int match_idx = 0;match_idx < matches_in.numRows(); match_idx++) {
    int current_scan = matches_in.getInteger("scan", match_idx);
    if (last_scan != current_scan) {
      //process the scans between the first and match_idx-1.
      //find the best row and calculate the bonferroni corrected p-value.
      carp(CARP_DEBUG,"Collaping %d %d %d",last_scan, first_row, match_idx-1);
      getBestScore(matches_in, first_row, match_idx-1, column_name, best_row, best_score);  
      carp(CARP_DEBUG,"Copying row best: %d %lf", best_row, best_score);
      //update matches out
      int new_row = matches_out.addRow();
      matches_in.copyToRow(matches_out, best_row, new_row);

      //update first_row and last scan.
      first_row = match_idx;
      last_scan = current_scan;
    }
  }

  carp(CARP_DEBUG,"Doing last entry");
  //finish the last entry
  getBestScore(matches_in, first_row, matches_in.numRows() - 1, column_name, best_row, best_score);
  carp(CARP_DEBUG,"best: %d %lf", best_row, best_score);
  int new_row = matches_out.addRow();
  matches_in.copyToRow(matches_out, best_row, new_row);
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
