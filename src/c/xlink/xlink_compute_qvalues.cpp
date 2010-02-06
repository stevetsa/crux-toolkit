#include "xhhc.h"
#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"

extern "C" {
#include "objects.h"
#include "scorer.h"
#include "spectrum_collection.h"
}

#include "DelimitedFile.h"

#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
//#define bin_width_mono 1.0005079


#define NUM_ARGUMENTS 8
#define NUM_OPTIONS 5

void getBestBonf(DelimitedFile& matches, int start, int stop, 
  int& best_index, double& best_bonf) {

  int numScans = stop - start + 1;
  //if there is only 1 scan, then return it.
  if (numScans == 1) {
    best_index = start;
    double pvalue = matches.getDouble("p-value");
    int ntests = matches.getInteger("matches/spectrum");
    best_bonf = bonferroni_correction(pvalue, ntests);
  } else {
    map<int, pair<int, double> > charge_best_score; 
    map<int, int> charge_ntests;

    for (int match_idx = start;match_idx <= stop;match_idx++) {
      int charge = matches.getInteger("charge", match_idx);
      double pvalue = matches.getDouble("p-value", match_idx);

      if (charge_best_scan_idx.find(charge) == charge_best_scan_idx.end()) {
        charge_best_scan_idx[charge] = make_pair(match_idx, pvalue);
      } else if (pvalue < charge_best_scan_idx[charge].second) {
        charge_best_scan_idx[charge] = make_pair(match_idx, pvalue);
      }
      if (charge_ntests.find(charge) == charge_ntests.end()) {
        int ntests = matches.getDouble("matches/spectrum", match_idx);
        charge_ntests[charge] = ntests;
      }
    }
  }
  
  int best_charge = charge_best_score.begin() -> first;
  double best_charge_pvalue = charge_best_score[best_charge] -> second.second;
  double best_charge_pvalue_bonf = 
    bonferroni_correct(best_charge_pvalue, charge_ntests[best_charge]);

  int best_charge_idx = charge_best_score[best_charge] -> second.first;

  for (map<int, pair<int, double> >::iterator iter = 
    charge_best_score.begin();
    iter != charge_best_score.end();
    ++iter) {

    double current_charge_pvalue_bonf =
      bonferroni_correct(iter -> second.second, charge_ntests[iter -> first]);
    if (current_charge_pvalue_bonf < best_charge_pvalue_bonf) {
      best_charge = iter -> first;
      best_charge_pvalue = iter -> second.second;
      best_charge_pvalue_bonf = current_charge_pvalue_bonf;
      best_charge_idx = charge_best_score[best_charge] -> second.first;  
    }
    
  }

  int ntests_total = 0;
  for (map<int, int>::iterator iter =
    charge_ntests.begin();
    iter != charge_ntests.end();
    ++iter) {

    ntests_total += iter -> second;
  }

  best_bonf = bonferonni_correct(best_charge_pvalue, ntests_total);
  best_index = best_charge_idx;
}

void collapseScans(DelimitedFile& matches_in, DelimitedFile& matches_out) {

  matches_out.clear();
  matches_out.addColumns(matches_in.getColumnNames());

  matches_in.sortByColumn("scan", INTEGER_TYPE, ASCENDING);

  matches_out.addColumn("p-value bonf.");

  int last_scan = matches_in.getInteger("scan", 0);
  int first_row = 0;

  

  for (int match_idx = 0;match_idx < matches_in.numRows(); match_idx++) {
    int current_scan = matches_in.getInteger("scan");
    if (last_scan != current_scan) {
      //process the scans between the first and match_idx-1.
      int best_row;
      double best_bonf;
      //find the best row and calculate the bonferroni corrected p-value
      getBestBonf(matches_in, first_row, match_idx-1, best_index, best_bonf);  

      //update matches out
      int new_row = matches_out.addRow();
      matches_in.copyRow(best_index, matches_out, new_row);
      matches_out.setValue("p-value bonf.", best_bonf);


      //update first_row and last scan.
      first_row = match_idx;
      last_scan = current_scan;
    }
  }

  //finish the last entry
  getBestBonf(matches_in, first_row, matches_in.numRows() - 1, best_row, best_bonf);
  new_row = matches_out.addRow();
  matches_in.copyRow(best_index, matches_out, new_row);
  matches_out.setValue("p-value bonf.", new_row, best_bonf);




}


int main(int argc, char** argv){

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  int num_options = NUM_OPTIONS;
  const char* option_list[NUM_OPTIONS] = {
    "verbosity",
    "version",
    "use-mgf",
    "ion-tolerance",
    "precision"
  };

  

  /* Define required command line arguments */
  int num_arguments = NUM_ARGUMENTS;
  const char* argument_list[NUM_ARGUMENTS] = {"peptide A",
                                              "peptide B",
					      "pos A",
					      "pos B",
					      "link mass",
					      "charge state",
					      "scan number",
					      "ms2 file"};

  
  /* for debugging of parameter processing */
  //set_verbosity_level( CARP_DETAILED_DEBUG );
  set_verbosity_level( CARP_ERROR );
  
  /* Set default values for parameters in parameter.c */
  initialize_parameters();

  /* Define optional and required command line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  parse_cmd_line_into_params_hash(argc, argv, "xlink-assign-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
 
  // create new ion series


  //Read in targets.
  string target_file = output_dir + "/search.target.txt";
  DelimitedFile target_matches(target_file);

  //Read in decoys.
  string decoy_file = output_dir + "/search.decoy.txt";
  DelimitedFile decoy_matches(decoy_file);

  //Collapse to the best target/decoy and calculate bonferonni corrected p-value.
  DelimitedFile target_matches_bonf;
  collapseScans(target_matches, target_matches_bonf);
  
  DelimitedFile decoy_matches_bonf;
  collapseScans(decoy_matches, decoy_matches_bonf);



  //Sort by increasing p-value.
  target_matches_bonf.sortByColumn("p-value bonf.", DOUBLE_TYPE, ASCENDING);
  decoy_matches_bonf.sortByColumn("p-value bonf.", DOUBLE_TYPE, ASCENDING); 
  
  target_matches.addColumn("q-value b-h");
  target_matches.addColumn("q-value decoy");

  int decoy_idx = 0;

  for (int target_idx = 0; target_idx < target_matches_bonf.numRows(); target_idx++) {
    //calculate q-value by b-h
    double current_pvalue = target_matches_bonf.getDouble("p-value bonf.", target_idx);
    double q_value_bh = (double)(target_idx + 1) / 
      (double) target_matches_bonf.size() * 
      current_pvalue;

    while ((decoy_idx < decoy_matches_bonf.size()) && 
      (decoy_matches_bonf.getDouble("p-value bonf.") <= current_pvalue)) {
      decoy_idx++;
    }

    double q_value_decoy = 0;
    if (decoy_idx != 0) {
      q_value_decoy = decoy_idx / (target_idx + 1);
    }

    target_matches_bonf.setDouble("q-value b-h", target_idx, q_value_bh);
    target_matches_bonf.setDouble("q-value decoy", target_idx, q_value_decoy);
  }
  
  string result_file = outputdir + "/qvalues.target.txt";
  target_matches_bonf.saveData(outputdir);
}
