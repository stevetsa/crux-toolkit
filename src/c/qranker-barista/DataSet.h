#ifndef DATASET_H_
#define DATASET_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <cmath>
#include <map>
#include "BipartiteGraph.h"
using namespace std;



class Dataset
{
 public:

  Dataset();
  ~Dataset();
  void load_data();
  void load_data(string &summary_fn, string &psm_fn);
  void load_psm_data_for_training(string &summary_fn, string &psm_fn);
  void load_psm_data_for_reporting_results();
  void load_prot_data();
  
  inline void set_input_dir(string input_dir){in_dir = input_dir;}
  void normalize_psms();
  
  inline double* psmind2features(int psmind){return (psmind_to_features+num_features*psmind);}
  inline int psmind2label(int psmind){return psmind_to_label[psmind];}
  inline int psmind2scan(int psmind){return psmind_to_scan[psmind];}
  inline int psmind2charge(int psmind){return psmind_to_charge[psmind];}
  inline int* psmind2charges(int psmind){return (psmind_to_charge+psmind_to_ofst[psmind]);}
  inline int psmind2pepind(int psmind){return psmind_to_pepind[psmind];}
  inline int* psmind2pepinds(int psmind){return (psmind_to_pepind+psmind_to_ofst[psmind]);}
  inline int psmind2num_pep(int psmind){return psmind_to_num_pep[psmind];}
  inline double* psmind2neutral_mass(int psmind){return (psmind_to_neutral_mass+psmind_to_ofst[psmind]);}
  inline double* psmind2peptide_mass(int psmind){return (psmind_to_peptide_mass+psmind_to_ofst[psmind]);}

  inline string ind2pep(int ind){return ind_to_pep[ind];}
  inline int get_num_psms(){return num_psms;}
  inline int get_num_features(){return num_features;}

  inline int get_num_peptides(){return num_pep;}
  inline int pepind2num_psm(int pepind){return pepind_to_psminds.get_range_length(pepind);}
  inline int* pepind2psminds(int pepind){return pepind_to_psminds.get_range_indices(pepind);}

  inline int pepind2num_prot(int pepind){return pepind_to_protinds.get_range_length(pepind);}
  inline int* pepind2protinds(int pepind){return pepind_to_protinds.get_range_indices(pepind);}

  inline int get_num_proteins(){return num_prot;}
  inline string ind2prot(int ind){return ind_to_prot[ind];}
  inline int protind2label(int protind){return protind_to_label[protind];}
  inline int protind2num_pep(int protind){return protind_to_pepinds.get_range_length(protind);}
  inline int* protind2pepinds(int protind){return protind_to_pepinds.get_range_indices(protind);}
  inline int protind2num_all_pep(int protind){return protind_to_num_all_pep[protind];}
  //returns false if not subset, true if yes, subset
  inline bool is_prot_subset(int protind1, int protind2){return protind_to_pepinds.is_subset(protind1, protind2);}

 protected:
  int num_psms;
  int num_pos_psms;
  int num_neg_psms;
  int num_features;
  int num_all_pep_in_psms;
  double* psmind_to_features;
  int* psmind_to_label;
  int *psmind_to_pepind;
  int *psmind_to_num_pep;
  int *psmind_to_ofst;
  int *psmind_to_scan;
  int *psmind_to_charge;
  double *psmind_to_neutral_mass;
  double *psmind_to_peptide_mass;
  map <int, string> ind_to_pep;


  int num_pep;
  int num_pos_pep;
  int num_neg_pep;
  BipartiteGraph pepind_to_psminds;
  BipartiteGraph pepind_to_protinds;


  int num_prot;
  int num_pos_prot;
  int num_neg_prot;
  int *protind_to_label;
  BipartiteGraph protind_to_pepinds;
  int *protind_to_num_all_pep;
  map <int, string> ind_to_prot;

  string in_dir;
};



#endif /*DATASET_H */
