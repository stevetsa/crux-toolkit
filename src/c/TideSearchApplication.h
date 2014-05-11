#ifndef TIDESEARCHAPPLICATION_H
#define TIDESEARCHAPPLICATION_H

#include "CruxApplication.h"
#include "TideMatchSet.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <gflags/gflags.h>
#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "tide/theoretical_peak_set.h"
#include "tide/max_mz.h"

using namespace std;

class TideSearchApplication : public CruxApplication {

protected:

  static bool HAS_DECOYS;

  /**
   * Free all existing mods
   */
  void cleanMods();
  
  // with thread.
  void search(
    const vector<SpectrumCollection::SpecCharge>* spec_charges,
    ActivePeptideQueue* active_peptide_queue[],
    ProteinVec& proteins,
    vector<const pb::AuxLocation*>& locations,
    double precursor_window,
    WINDOW_TYPE_T window_type[],
    double spectrum_min_mz,
    double spectrum_max_mz,
    int min_scan,
    int max_scan,
    int min_peaks,
    int search_charge,
    int top_matches,
    double highest_mz,
    OutputFiles* output_files,
    ofstream* target_file[],
    ofstream* decoy_file[],
    bool compute_sp
  );

  void collectScoresCompiled(
    ActivePeptideQueue* active_peptide_queue,
    const Spectrum* spectrum,
    const ObservedPeakSet& observed,
    TideMatchSet::Arr* match_arr,
    int queue_size,
    int charge
  );

  void computeWindow(
    const SpectrumCollection::SpecCharge& sc,
    WINDOW_TYPE_T window_type,
    double precursor_window,
    double* out_min,
    double* out_max
  );

  struct ScSortByMz {
    ScSortByMz(double precursor_window) { precursor_window_ = precursor_window; }
    bool operator() (const SpectrumCollection::SpecCharge x,
                     const SpectrumCollection::SpecCharge y) {
      return (x.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_) * x.charge <
             (y.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_) * y.charge;
    }
    double precursor_window_;
  };

public:

  /**
   * Constructor
   */
  TideSearchApplication();

  /**
   * Destructor
   */
  ~TideSearchApplication();

  /**
   * Main method
   */

  int NUM_THREADS;

  virtual int main(int argc, char** argv);

  static bool hasDecoys();

  /**
   * Returns the command name
   */
  virtual string getName();

  /**
   * Returns the command description
   */
  virtual string getDescription();

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory();

  virtual COMMAND_T getCommand();

  // for threading spec charges
  void * searchthread(void *threadarg);

  // struct holding necessary information for each thread to run.
  struct thread_data {

    const vector<SpectrumCollection::SpecCharge>* spec_charges;
    ActivePeptideQueue* active_peptide_queue;
    ProteinVec proteins;
    vector<const pb::AuxLocation*> locations;
    double precursor_window;
    WINDOW_TYPE_T window_type;
    double spectrum_min_mz;
    double spectrum_max_mz;
    int min_scan;
    int max_scan;
    int min_peaks;
    int search_charge;
    int top_matches;
    double highest_mz;
    OutputFiles* output_files;
    ofstream* target_file;
    ofstream* decoy_file;
    bool compute_sp;
    long lo;
    long hi;

    void setParams(const vector<SpectrumCollection::SpecCharge>* spec_charges_,
    	ActivePeptideQueue* active_peptide_queue_, ProteinVec proteins_,
    	vector<const pb::AuxLocation*> locations_, double precursor_window_,
    	WINDOW_TYPE_T window_type_, double spectrum_min_mz_, double spectrum_max_mz_,
    	int min_scan_, int max_scan_, int min_peaks_, int search_charge_, int top_matches_,
    	double highest_mz_, OutputFiles* output_files_, ofstream* target_file_,
    	ofstream* decoy_file_, bool compute_sp_, long lo_, long hi_) {
      spec_charges = spec_charges_;
      active_peptide_queue = active_peptide_queue_;
      proteins = proteins_;
      locations = locations_;
      precursor_window = precursor_window_;
      window_type = window_type_;
      spectrum_min_mz = spectrum_min_mz_;
      spectrum_max_mz = spectrum_max_mz_;
      min_scan = min_scan_;
      max_scan = max_scan_;
      min_peaks = min_peaks_;
      search_charge = search_charge_;
      top_matches = top_matches_;
      highest_mz = highest_mz_;
      output_files = output_files_;
      target_file = target_file_;
      decoy_file = decoy_file_;
      compute_sp = compute_sp_;
      lo = lo_;
      hi = hi_;
    }
  };
  
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
