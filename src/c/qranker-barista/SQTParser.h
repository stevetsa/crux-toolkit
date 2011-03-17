#ifndef SQTPARSER_H
#define SQTPARSER_H

#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <string>
#include <math.h>
#include <map>
#include <set>
#include <cstring>
#include "SpecFeatures.h"
#include "BipartiteGraph.h"
using namespace std;

typedef enum {TRYPSIN_ENZ,CHYMOTRYPSIN_ENZ,ELASTASE_ENZ} enzyme;

class SQTParser{
 public:
  struct sqt_match{
    int scan;
    int charge;
    double precursor_mass;
    int num_sequence_comparisons;
    vector <int> xcorr_rank;
    vector<int> sp_rank;
    vector<double> calc_mass;
    vector<double> delta_cn;
    vector<double> xcorr_score;
    vector<double> sp_score;
    vector <double> num_ions_matched;
    vector<double> num_total_ions;
    vector<string> peptides;
    vector<int> num_proteins_in_match;
    vector<string> proteins;
  };
  /*
 SQTParser() : num_spectra(0),num_psm(0),num_pos_psm(0),num_neg_psm(0),num_features(0),
    num_spec_features(0), num_pep(0),num_pos_pep(0), num_neg_pep(0),num_prot(0),num_pos_prot(0),num_neg_prot(0),
   num_mixed_labels(0),psmind(0){}
  */
  SQTParser(int capacity);
  ~SQTParser();
  int run();
  void erase_matches();
  inline void set_num_features(int nf) {num_features = nf;}
  inline void set_num_spec_features(int nsf) {num_spec_features = nsf;}
  inline int get_num_features() const {return num_features;}
  inline void set_output_dir(string output_dir){out_dir = output_dir;}
  inline string& get_output_dir(){return out_dir;}
  inline void set_input_dir(string input_dir){in_dir = input_dir;}
  inline void set_input_dir_ms2(string input_dir){in_dir_ms2 = input_dir;}
  inline string& get_input_dir(){return in_dir;}
  inline void set_db_name(string database){db_name = database;}
  inline void set_decoy_prefix(string prefix){decoy_prefix = prefix;}
  inline void set_enzyme(enzyme enz){e = enz;}
  int set_command_line_options(int argc, char *argv[]);
  
  void read_sqt_file(ifstream &is, string &decoy_prefix, int final_hits_per_spectrum, enzyme enz, int pass);
  int parse_sqt_spectrum_matches(ifstream &is, sqt_match &m);
  void read_S_line(ifstream &is, sqt_match &m);
  void read_M_line(ifstream &is, sqt_match &m);
  void digest_database(ifstream &f_db, enzyme e);
  int cntEnzConstraints(string& seq,enzyme enz);
  
  static int cntEnz(const string& peptide, enzyme enz);
  static double isTryptic(const char n,const char c);
  static double isChymoTryptic(const char n,const char c);
  static double isElastasic(const char n,const char c);
  static double isEnz(const char n,const char c, enzyme enz);
  void extract_psm_features(sqt_match &m, enzyme enz, double *x, int i);
  void extract_features(sqt_match &m, string &decoy_prefix, int hits_read, int final_hits,enzyme enz);
  void add_matches_to_tables(sqt_match &m, string &decoy_prefix, int hits_read, int final_hits);
  void allocate_feature_space();
  void fill_graphs();
  
  void save_data_in_binary(string out_dir);
  void load_data_in_binary(string dir_in);
  void clean_up(string dir);
  
  //spec features generator
  SpecFeaturesGenerator sfg;
  

 protected:
  //auxiliary variables
  sqt_match m;
  int num_mixed_labels;
  map<int,set<int> > pepind_to_protinds_map;
  map<int,set<int> > protind_to_pepinds_map;
  map<int,set<int> > pepind_to_psminds_map;
  map<int,int> psmind_to_label_map;
  map<int,int> pepind_to_label_map;
  map<int,int> protind_to_label_map;

  //summary of the dataset
  int num_features;
  int num_spec_features;
  int num_spectra;
  int num_psm;
  int num_pos_psm;
  int num_neg_psm;
  int num_pep;
  int num_pos_pep;
  int num_neg_pep;
  int num_prot;
  int num_pos_prot;
  int num_neg_prot;
  
  int psmind;

  //psm feature vector
  double *x;
  //spec feature vector
  double *xs;
  //psm info
  int* psmind_to_pepind;
  int* psmind_to_scan;
  int* psmind_to_charge;
  int* psmind_to_label;

  //peptide info
  map<string,int> pep_to_ind;
  map<int,string> ind_to_pep;
  int* pepind_to_label;
  BipartiteGraph pepind_to_protinds;
  BipartiteGraph pepind_to_psminds;
  
  //protein info
  map<string,int> prot_to_ind;
  map<int,string> ind_to_prot;
  BipartiteGraph protind_to_pepinds;
  int* protind_to_label;

  //digested database info
  map<string,int> protein_to_num_all_pep_map;  
  map<int,int> protind_to_num_all_pep_map;  
  int *protind_to_num_all_pep;

  //writing out data
  string in_dir;
  string in_dir_ms2;
  string out_dir;
  string db_name;
  vector<string> sqt_file_names;
  vector<string> ms2_file_names;

  ofstream f_psm;
  ofstream f_spec_feat;
  ofstream f_psm_spec_feat;

  
  //final hits per spectrum
  int fhps;
  //decoy prefix
  string decoy_prefix;
  //enzyme
  enzyme e;
  //max peptide length to be considered
  int max_len;
  //min peptide length to be considered
  int min_len;
  
};

#endif
