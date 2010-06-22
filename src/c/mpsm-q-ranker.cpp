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
#include <fstream>

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
#include "DelimitedFileReader.h"

#include "KrokhinRetentionPredictor.h"

using namespace std;
#include "./q-ranker/PSMDescription.h"
#include "./q-ranker/DataSet.h"
#include "./q-ranker/Scores.h"
#include "./q-ranker/SetHandler.h"
//#include "./q-ranker/Net.h"
#include "./q-ranker/Caller.h"
//#include "./q-ranker/Globals.h"

#include "ChargeIndex.h"

using namespace std;

/* 
 * Private function declarations.  Details below
 */

void writeResults(DelimitedFileReader& input, const char* output_path);
void writeSubsamples(DelimitedFileReader& matches, vector<unsigned int>& subsample_indices, const char* output_path);
void collapseResults(const char* input_path, const char* output_path, const char* column_name);

int argmax(vector<double>& a) {
  if (a.size() == 0) return -1;
  int ans = 0;

  double max = a[0];
  for (unsigned int idx=1;idx < a.size();idx++) {
    if (a[idx] > max) {
      ans = idx;
      max = a[idx];
    }
  }

  return ans;
}


void qcInitiate2(NSet sets, unsigned int numFeatures, 
  std::vector<unsigned int>& numSpectra, char ** featureNames, double pi0);

Caller* getCallerQR();


void run_mpsm_q(
  char* psm_result_folder, 
  char* fasta_file, 
  char* feature_file); 


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
    "subsample-percent",
    "random-sample"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  const char* argument_list[] = {
    "protein input",
    "search results directory"
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);
  carp(CARP_INFO,"Initialize run");
  initialize_run(MPSM_QRANKER_COMMAND, argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  /* Get arguments */
  char* psm_dir = get_string_parameter("search results directory");
  char* protein_input_name = get_string_parameter("protein input");

  // TODO (BF oct-22-09): consider adding feature file to OutputFiles
  char* feature_file = "features.txt";


  if (feature_file != NULL) {
    prefix_fileroot_to_name(&feature_file);
  }


  /* Perform the analysis */
  //MATCH_COLLECTION_T* match_collection = NULL;

  DelimitedFile qranker_results;

  carp(CARP_INFO,"run_mpsm_q");
  run_mpsm_q(psm_dir,
             protein_input_name,
             feature_file);

  carp(CARP_INFO,"free psm_dir");
  // clean up
  free(psm_dir);
  carp(CARP_INFO,"free protein input name");
  free(protein_input_name);
  carp(CARP_INFO,"free feature_file");
  //free(feature_file);


  carp(CARP_INFO, "crux q-ranker finished.");
  exit(0);

}

//const int number_features = 11;
const char* feature_names[] = {
  "XCorr",
  "avedM",
  "maxdM",
  "theoryFrac",
  "obsFrac",
  "numTests",
  "aveMissed",
  "maxMissed",
  "aveLength",
  "maxLength",
  "charge1",
  "charge2",
  "charge3",
  "numPep",
  "aveRT",
  "maxRT"
};

const int number_features = sizeof(feature_names) / sizeof(char*);


double getAveDiff(vector<double>& as, vector<double>& bs) {
  double ans = 0;

  for (unsigned int i=0;i<bs.size();i++) {
    double abs_diff = fabs(as[i] - bs[i]);
    ans += abs_diff;
  }

  return ans / (double)bs.size();

}

double getMaxDiff(vector<double>& as, vector<double>& bs) {

  double ans = 0;

  for (unsigned int i=0;i<bs.size();i++) {
    double abs_diff = fabs(as[i] - bs[i]);

    if (abs_diff > ans) {
      ans = abs_diff;
    }
  }
  return ans;
}

double getTheoryFrac() {
  return 0.0;
}

double getObsFrac() {
  return 0.0;
}

int get_missed_cleavage_sites(string& sequence) {
  unsigned int aa_idx = 0;
  int missed_count = 0;
  for (;aa_idx < sequence.length()-1; ++aa_idx) {
    if (sequence[aa_idx] == 'K' ||
        sequence[aa_idx] == 'R'){

      if (sequence[aa_idx+1] == 'P') {
        continue;
      } else {
        ++missed_count;
      }
    }
  }
  return missed_count;
}

double getAveMissed(vector<string>& sequences) {
  
  double ans = 0;

  for (unsigned int i=0;i<sequences.size();i++) {
    ans += get_missed_cleavage_sites(sequences[i]);
  }
  
  return ans / (double)sequences.size();
}

int getMaxMissed(vector<string>& sequences) {

  int ans = 0;

  for (unsigned int i=0;i<sequences.size();i++) {
    int current = get_missed_cleavage_sites(sequences[i]);

    if (current > ans) {
      ans = current;
    }
  }
  return ans;
}


double getAveLength(vector<string>& sequences) {

  double ans = 0.0;
  for (unsigned int idx=0;idx<sequences.size();idx++) {
    ans += sequences[idx].length();
  }

  return ans / (double)sequences.size();

}

unsigned int getMaxLength(vector<string>& sequences) {
  
  unsigned int ans = 0;
  for (unsigned int i=0;i<sequences.size();i++) {
    if (sequences[i].length() > ans) {
      ans = sequences[i].length();
    }
  }
  return ans;
}

void getRTimes(vector<string>& sequences, vector<double>& rtimes) {

  KrokhinRetentionPredictor rt_predictor;

  rtimes.clear();

  for (unsigned int idx=0;idx < sequences.size();idx++) {
    double rt = 
      rt_predictor.predictRTimeS(sequences[idx].c_str());
    rtimes.push_back(rt);
  }
}

double getAveRT(vector<double>& rtimes) {

  if (rtimes.size() <= 1) return 0.0;
  double ans = 0.0;
  int count = 0;
  for (unsigned int idx1 = 0;idx1 < rtimes.size() - 1; idx1++) {
    for (unsigned int idx2 = idx1 + 1;idx2 < rtimes.size() ;idx2++) {
      count++;
      ans += fabs(rtimes[idx1] - rtimes[idx2]);
    }
  }
  

  return ans / (double)count;
}

double getMaxRT(vector<double>& rtimes) {
  if (rtimes.size() <= 1) return 0.0;
  double ans = -1;
  for (unsigned int idx1 = 0;idx1 < rtimes.size() - 1; idx1++) {
    for (unsigned int idx2 = idx1 + 1;idx2 < rtimes.size() ;idx2++) {
      double current = fabs(rtimes[idx1] - rtimes[idx2]);
      if (current > ans) {
        ans = current;
      }
    }
  }
  return ans;
}




double* get_mpsm_features(DelimitedFileReader& matches) {

  double* features = (double*)mycalloc(number_features, sizeof(double));

  vector<double> neutral_masses;
  matches.getDoubleVectorFromCell("spectrum neutral mass", neutral_masses);
  vector<double> peptide_masses;
  matches.getDoubleVectorFromCell("peptide mass", peptide_masses);

  int matches_spectrum = matches.getInteger("matches/spectrum");
  double lnSM = 0;
  if (matches_spectrum > 0) {
    lnSM = log(matches_spectrum);
  }

  vector<string> peptide_sequences;
  matches.getStringVectorFromCell("sequence", peptide_sequences);
  
  vector<double> rtimes;
  getRTimes(peptide_sequences, rtimes);

  string string_charge = matches.getString("charge");
  ChargeIndex charge(string_charge);

  features[0] = matches.getDouble("xcorr score");         //xcorr_score;
  features[1] = getAveDiff(neutral_masses, peptide_masses); //abs_ave_diff;
  features[2] = getMaxDiff(neutral_masses, peptide_masses); //abs_max_diff;
  features[3] = getTheoryFrac();
  features[4] = getObsFrac();
  features[5] = lnSM;//lnSM;
  features[6] = getAveMissed(peptide_sequences);
  features[7] = getMaxMissed(peptide_sequences);
  features[8] = getAveLength(peptide_sequences);//ave_pep_len
  features[9] = getMaxLength(peptide_sequences);//max_pep_len;
  features[10] = charge.numCharge(1); //ncharge1;
  features[11] = charge.numCharge(2);//ncharge2;
  features[12] = charge.numCharge(3);//ncharge3;
  features[13] = peptide_sequences.size();//num_peptides;
  features[14] = getAveRT(rtimes);
  features[15] = getMaxRT(rtimes);

  return features;

}

/**
 * get a sorted list of subsamples from an 0..total_samples-1 indices.
 */
void get_random_indices(unsigned int total_samples, 
  unsigned int num_subsamples, 
  vector<unsigned int>& indices) {

  indices.clear();

  cout <<"Generating list:"<<total_samples<<endl;
  //generate a list of 0..total_samples-1 indices.
  for (unsigned int idx = 0;idx < total_samples; idx++) {
    indices.push_back(idx);
  }

  cout <<"Shuffling"<<endl;
  //shuffle the indices.
  random_shuffle(indices.begin(), indices.end());

  cout <<"Resizing to "<<num_subsamples<<endl;
  //truncate to the num_subsamples
  indices.resize(num_subsamples);

  cout <<"Sorting"<<endl;
  //sort the subsamples in ascending order.
  sort(indices.begin(), indices.end());

}

void getTopIndices(DelimitedFileReader& matches, 
  vector<unsigned int>& indices) {
  
  matches.reset();
  indices.clear();

  string last_charge = matches.getString("charge");
  int current_index = 0;
  indices.push_back(0);
  while (matches.hasNext()) {

    string current_charge = matches.getString("charge");
    if (current_charge != last_charge) {
      indices.push_back(current_index);
      last_charge = current_charge;
    }
    current_index++;
    matches.next();
  }

}

void getSubsampleIndices(DelimitedFileReader &matches,
  vector<unsigned int>& subsample_indices,
  double subsample_percent) {

  subsample_indices.clear();
  if (subsample_percent >= 100.0) {
    for (unsigned int idx=0;idx<matches.numRows();idx++) {
      subsample_indices.push_back(idx);
    }
    return;
  }

  if (get_boolean_parameter("random-sample")) {
    cout <<"Generating random indices:"<<endl;
    get_random_indices(matches.numRows(), 
      (unsigned int)(matches.numRows()*subsample_percent/100.0), 
      subsample_indices);
  } else {

    vector<unsigned int> top_indices;
    getTopIndices(matches, top_indices);
    cout <<"There are "<<top_indices.size()<<" of "<<matches.numRows()<<endl;
    vector<unsigned int> temp_indices;
    get_random_indices(top_indices.size(),
      (unsigned int)(top_indices.size()*subsample_percent/100.0),
      temp_indices);
    for (unsigned int idx=0;idx<temp_indices.size();idx++) {
      subsample_indices.push_back(top_indices[temp_indices[idx]]);
    }
  }  
}


void registerMatches(DelimitedFileReader& matches, 
  SetType set_idx, 
  vector<unsigned int>& subsample_indices,
  FILE* feature_file) {

  cout <<"Inside registerMatches"<<endl;
  double* features = NULL;
  unsigned int sample_count = 0;

  matches.reset();
  cout<<"Building feature set"<<endl;
  for (unsigned int sub_iidx = 0; 
    sub_iidx < subsample_indices.size(); 
    sub_iidx++) {

    unsigned int sub_idx = subsample_indices[sub_iidx];
    //cout<<"going to :"<<sub_idx<<" "<<sample_count<<" "<<matches.numRows()<<endl;
    for (;sample_count < sub_idx;sample_count++) {
      matches.next();
    }
    features = get_mpsm_features(matches);

    if ((int)set_idx == 0) {
       string id = matches.getString("scan")+"-"+matches.getString("charge")+"-"+matches.getString("xcorr rank"); //unique;
     
       qcRegisterPSM(set_idx, (char*)id.c_str(), features);
    } else {
       qcRegisterPSM(set_idx, NULL, features);
    }
    
    if (feature_file != NULL) {
      fprintf(feature_file,"%d\t",matches.getInteger("scan"));
      fprintf(feature_file,"%s\t",matches.getString("charge").c_str());
      fprintf(feature_file,"%d",(int)set_idx);
      for (int idx=0;idx < number_features;idx++) {
        fprintf(feature_file,"\t%f",features[idx]);
      }
      fprintf(feature_file,"\n");
    }
    //carp(CARP_INFO,"freeing features");
    free(features);
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
  char* feature_file){ 

  DelimitedFileReader target_search_results;
  vector<DelimitedFileReader> decoy_search_results;   
  vector<unsigned int> num_matches;
  vector<unsigned int> num_subsamples;
  vector<vector<unsigned int> > subsample_indices_indices;
    
  char* output_dir = get_string_parameter("output-dir");

  double subsample_percent = get_double_parameter("subsample-percent");

  carp(CARP_INFO,"Subsample percent:%g",get_double_parameter("subsample-percent"));

  carp(CARP_INFO,"Fasta file:%s",fasta_file);

  double pi0 = get_double_parameter("pi0");
  FILE* feature_fh = NULL;
//  int set_idx = 0;
  
  // optional feature_file
  if(feature_file != NULL){  
    if((feature_fh = fopen(feature_file, "w")) == NULL){
      carp(CARP_FATAL, "Problem opening output file %s", feature_file);
      return;
    }
    fprintf(feature_fh, "scan\tcharge\tset");
    for (int idx=0;idx<number_features;idx++) {
      fprintf(feature_fh,"\t%s", feature_names[idx]);
    }
    fprintf(feature_fh, "\n");


  }

  carp(CARP_DETAILED_DEBUG, "Created feature file");

  //load the results file.
  carp(CARP_INFO,"counting target");



  string target_path = psm_result_folder+string("/search.target.txt");

  vector<DelimitedFileReader> search_results;
  search_results.push_back(target_search_results);
  search_results.push_back(target_search_results);
  search_results.push_back(target_search_results);
  search_results[0].loadData(target_path.c_str());

  //num_matches.push_back(search_results[0].numRows());
  //num_subsamples.push_back((unsigned int)(search_results[0].numRows() * subsample_percent));
  //carp(CARP_INFO,"Target count:%u",search_results[0].numRows());

  for (int idx=1;idx<=2;idx++) {
    carp(CARP_INFO,"idx:%d",idx);
    carp(CARP_INFO,"Building");
    string decoy_path = 
      psm_result_folder + 
      string("/search.decoy-") + 
      DelimitedFile::to_string<int>(idx) + 
      string(".txt");
    //carp(CARP_INFO,"reading %s",decoy_path.c_str());
    search_results[idx].loadData(decoy_path);
    //num_matches.push_back(search_results[idx].numRows());
    //num_subsamples.push_back((unsigned int)(search_results[idx].numRows() * subsample_percent));
    //carp(CARP_INFO,"count:%u", search_results[idx].numRows());
  }

  for (int idx=0;idx<3;idx++) {
    num_matches.push_back(search_results[idx].numRows());
    vector<unsigned int> subsample_indices;
    subsample_indices_indices.push_back(subsample_indices);
    getSubsampleIndices(search_results[idx], subsample_indices_indices[idx], subsample_percent);
    num_subsamples.push_back(subsample_indices_indices[idx].size());
  }

  carp(CARP_INFO,"Calling qcInitiate");
  // Call that initiates q-ranker
  qcInitiate2((NSet)3, number_features, num_subsamples, (char**)feature_names, pi0);
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

  for (int i=0;i<3;i++) {
    carp(CARP_INFO,"Registering :%i %d of %d", i, num_subsamples[i], num_matches[i]);
    string output_path = string(output_dir) + string("/subsample.") + DelimitedFile::to_string<int>(i) + string(".txt");
    writeSubsamples(search_results[i], subsample_indices_indices[i], output_path.c_str());

    registerMatches(search_results[i], SetType(i), subsample_indices_indices[i], feature_fh);  
  }

  /***** Q-RANKER run *********/

  if (feature_fh != NULL) {
    fclose(feature_fh);
  }
    carp(CARP_INFO, "got to here");
    
    // Start processing
  qcExecute();  
  carp(CARP_INFO," Done executing q-ranker");  
  /* Retrieving target scores and qvalues after 
   * processing, the array should be numSpectra long and will be filled in 
   * the same order as the features were inserted */
  carp(CARP_INFO," Getting results");

  string output_path = string(output_dir) + string("/qranker.target.txt");
  writeResults(search_results[0], output_path.c_str()); 

  output_path = string(output_dir) + string("/qranker.decoy-1.txt");
  writeResults(search_results[1], output_path.c_str());
  
  output_path = string(output_dir) + string("/qranker.decoy-2.txt");
  writeResults(search_results[2], output_path.c_str());

  string in_path = string(output_dir) + string("/qranker.target.txt");
  output_path = string(output_dir) + string("/qranker.target.collapsed.txt");
  collapseResults(in_path.c_str(), output_path.c_str(), "q-ranker score");

  in_path = string(output_dir) + string("/qranker.decoy-1.txt");
  output_path = string(output_dir) + string("/qranker.decoy-1.collapsed.txt");
  collapseResults(in_path.c_str(), output_path.c_str(), "q-ranker score");

  in_path = string(output_dir) + string("/qranker.decoy-2.txt");
  output_path = string(output_dir) + string("/qranker.decoy-2.collapsed.txt");
  collapseResults(in_path.c_str(), output_path.c_str(), "q-ranker score");

  carp(CARP_INFO, "mpsm-q-ranker: done.");
}

void writeSubsamples(DelimitedFileReader& matches, vector<unsigned int>& subsample_indices, const char* output_path) {
  matches.reset();

  ofstream fout(output_path);

  for (unsigned int idx = 0; idx < matches.numCols() ; idx++) {
    if (idx != 0) {
      fout << "\t";
    }
    fout << matches.getColumnName(idx);
  }
  fout << endl;

  unsigned int sample_count = 0;
  for (unsigned int sub_iidx = 0; 
    sub_iidx < subsample_indices.size(); 
    sub_iidx++) {

    unsigned int sub_idx = subsample_indices[sub_iidx];
    //cout<<"going to :"<<sub_idx<<" "<<sample_count<<" "<<matches.numRows()<<endl;
    for (;sample_count < sub_idx;sample_count++) {
      matches.next();
    }

    fout << matches.getDataString() << endl;
  }  

  fout.close();

}

void writeResults(DelimitedFileReader& input, const char* output_path) {
  input.reset();
  ofstream fout(output_path);
  
  for (unsigned int idx = 0; idx < input.numCols() ; idx++) {
    if (idx != 0) {
      fout << "\t";
    }
    fout << input.getColumnName(idx);
  }
  fout << endl;
  
  int qranker_score_col = input.findColumn("q-ranker score");
  if (qranker_score_col == -1) {
    carp(CARP_FATAL,"Couldn't find q-ranker score in input");
  }

  Caller* pCaller = getCallerQR();


  while (input.hasNext()) {
    //create feature vector.
    //calculate q-ranker score.
    //carp(CARP_INFO,"get_mpsm_features");
    double* features = get_mpsm_features(input);
    //carp(CARP_INFO,"getQRankerScore");
    double qranker_score = pCaller -> getQRankerScore(features);

    //carp(CARP_INFO,"free features");
    free(features);
    //write everything out.
    //carp(CARP_INFO, "write out");
    for (unsigned int idx=0;idx<input.numCols();idx++) {
      //carp(CARP_INFO,"col %u",idx);
      if (idx != 0) {
        fout << "\t";
      }
      if (idx == (unsigned int)qranker_score_col) {
        //carp(CARP_INFO,"outputting qranker score:%lf",qranker_score);
        fout << qranker_score;
      } else {
        fout << input.getString(idx);
      }
    }
    fout << endl;
    //carp(CARP_INFO,"getting next record");
    input.next();
  }
  fout.close();
}

void collapseResults(const char* input_path, const char* output_path, const char* column_name) {
  DelimitedFileReader input(input_path, true);

  ofstream fout(output_path);

  //print the header.
  
  for (unsigned int idx = 0; idx < input.numCols() ; idx++) {
    if (idx != 0) {
      fout << "\t";
    }
    fout << input.getColumnName(idx);
  }
  fout << endl;

  vector<int> best_indices;

  double best_score = input.getFloat(column_name);
  string best_row = input.getDataString();
  
  int last_scan = input.getInteger("scan");

  while(input.hasNext()) {
    int current_scan = input.getInteger("scan");
    double current_score = input.getFloat(column_name);
    if (current_scan != last_scan) {
      //output best_row.
      fout << best_row << endl;
      best_row = input.getDataString();
      best_score = current_score;
      last_scan = current_scan;
    } else {
      if (current_score > best_score) {
        best_row = input.getDataString();
        best_score = current_score;
      }
    }
    input.next();
  }
  //output the last best row.
  fout << best_row << endl;
  fout.close();
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
