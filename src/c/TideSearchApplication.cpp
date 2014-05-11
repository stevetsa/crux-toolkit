#include <cstdio>
#include "tide/abspath.h"
#include "tide/records_to_vector-inl.h"

#include "carp.h"
#include "parameter.h"
#include "SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "TideSearchApplication.h"
#include <pthread.h>
#include <sstream>

extern HASH_T* parameters;
extern AA_MOD_T* list_of_mods[MAX_AA_MODS];
extern int num_mods;

bool TideSearchApplication::HAS_DECOYS = false;

TideSearchApplication::TideSearchApplication() {
}

TideSearchApplication::~TideSearchApplication() {
}

int TideSearchApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "precursor-window",
    "precursor-window-type",
    "spectrum-min-mz",
    "spectrum-max-mz",
    "min-peaks",
    "spectrum-charge",
    "scan-number",
    "top-match",
    "store-spectra",
    "concat",
    "compute-sp",
    "remove-precursor-peak",
    "remove-precursor-tolerance",
    "txt-output",
    "sqt-output",
    "pepxml-output",
    "mzid-output",
    "pinxml-output",
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "tide spectra file",
    "tide database index"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-search...");

  string cmd_line = "crux tide-search";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  NUM_THREADS = 2; // MINIMUM # = 1. (Meaning just main thread) Do not make
                   // this value below 1.
  carp(CARP_INFO, "Number of Threads: %d", NUM_THREADS); // is this legal ?

  string index_dir = get_string_parameter_pointer("tide database index");
  string peptides_file = index_dir + "/pepix";
  string proteins_file = index_dir + "/protix";
  string auxlocs_file = index_dir + "/auxlocs";
  string spectra_file = get_string_parameter_pointer("tide spectra file");

  double window = get_double_parameter("precursor-window");
  WINDOW_TYPE_T window_type[NUM_THREADS];
  for (int i = 0; i < NUM_THREADS; i++) {
      window_type[i] = get_window_type_parameter("precursor-window-type");
  }

  // Check spectrum-charge parameter
  string charge_string = get_string_parameter_pointer("spectrum-charge");
  int charge_to_search;
  if (charge_string == "all") {
    carp(CARP_DEBUG, "Searching all charge states");
    charge_to_search = 0;
  } else {
    charge_to_search = atoi(charge_string.c_str());
    if (charge_to_search < 1 || charge_to_search > 6) {
      carp(CARP_FATAL, "Invalid spectrum-charge value %s",
           charge_string.c_str());
    }
    carp(CARP_DEBUG, "Searching charge state %d", charge_to_search);
  }

  // Check scan-number parameter
  string scan_range = get_string_parameter_pointer("scan-number");
  int min_scan, max_scan;
  if (scan_range == "__NULL_STR") {
    min_scan = 0;
    max_scan = BILLION;
    carp(CARP_DEBUG, "Searching all scans");
  }
  else if (scan_range.find('-') == string::npos) {
    // Single scan
    min_scan = max_scan = atoi(scan_range.c_str());
    carp(CARP_DEBUG, "Searching single scan %d", min_scan);
  } else {
    if (!get_range_from_string(scan_range.c_str(), min_scan, max_scan)) {
      carp(CARP_FATAL, "The scan number range '%s' is invalid. "
           "Must be of the form <first>-<last>", scan_range.c_str());
    } else {
      if (min_scan > max_scan) {
        int tmp_scan = min_scan;
        min_scan = max_scan;
        max_scan = tmp_scan;
        carp(CARP_DEBUG, "Switched scan range min and max");
      }
      carp(CARP_DEBUG, "Searching scan range %d-%d", min_scan, max_scan);
    }
  }

  // Check concat parameter
  bool concat = get_boolean_parameter("concat");
  if (concat) {
    OutputFiles::setConcat();
  }

  // Check compute-sp parameter
  bool compute_sp = get_boolean_parameter("compute-sp");
  if (get_boolean_parameter("sqt-output") && !compute_sp){
    compute_sp = true;
    carp(CARP_INFO, "Enabling parameter compute-sp since SQT output is enabled "
                    " (this will increase runtime).");
  }

  const double max_mz = 0.0;

  carp(CARP_INFO, "Reading index %s", index_dir.c_str());
  // Read proteins index file
  ProteinVec proteins;
  pb::Header protein_header;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins,
      proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
  }
  carp(CARP_DEBUG, "Read %d proteins", proteins.size());

  // Read auxlocs index file
  vector<const pb::AuxLocation*> locations;
  if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file)) {
    carp(CARP_FATAL, "Error reading index (%s)", auxlocs_file.c_str());
  }
  carp(CARP_DEBUG, "Read %d auxlocs", locations.size());

  // Read peptides index file
  pb::Header peptides_header;
  HeadedRecordReader* peptide_reader[NUM_THREADS];
  for (int i = 0; i < NUM_THREADS; i++) {
    peptide_reader[i] = new HeadedRecordReader(peptides_file, &peptides_header);
  }
  if (!peptides_header.file_type() == pb::Header::PEPTIDES ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }
  
  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader.decoys();
  if (headerDecoyType != NO_DECOYS) {
    HAS_DECOYS = true;
    if (headerDecoyType == PROTEIN_REVERSE_DECOYS) {
      OutputFiles::setProteinLevelDecoys();
    }
  }
  MassConstants::Init(&pepHeader.mods());
  TideMatchSet::initModMap(pepHeader.mods());

  ActivePeptideQueue * active_peptide_queue[NUM_THREADS];
  for (int i = 0; i < NUM_THREADS; i++) {
    active_peptide_queue[i] = new ActivePeptideQueue(peptide_reader[i]->Reader(), proteins);
  }

  carp(CARP_INFO, "Reading spectra file %s", spectra_file.c_str());
  // Try to read file as spectrumrecords file
  SpectrumCollection spectra;
  pb::Header spectrum_header;
  string delete_spectra_file = "";
  if (!spectra.ReadSpectrumRecords(spectra_file, &spectrum_header)) {
    // Failed, try converting to spectrumrecords file
    carp(CARP_INFO, "Converting %s to spectrumrecords format",
                    spectra_file.c_str());
    string converted_spectra_file = get_string_parameter_pointer("store-spectra");
    if (converted_spectra_file.empty()) {
      char tmpnam_buffer[L_tmpnam];
      tmpnam(tmpnam_buffer);
      delete_spectra_file = converted_spectra_file = tmpnam_buffer;
    }
    carp(CARP_DEBUG, "New spectrumrecords filename: %s",
                     converted_spectra_file.c_str());
    if (!SpectrumRecordWriter::convert(spectra_file, converted_spectra_file)) {
      carp(CARP_FATAL, "Error converting %s to spectrumrecords format",
                       spectra_file.c_str());
    }
    carp(CARP_DEBUG, "Reading converted spectra file %s",
                     spectra_file.c_str());
    // Re-read converted file as spectrumrecords file
    if (!spectra.ReadSpectrumRecords(converted_spectra_file, &spectrum_header)) {
      carp(CARP_DEBUG, "Deleting %s", converted_spectra_file.c_str());
      remove(converted_spectra_file.c_str());
      carp(CARP_FATAL, "Error reading spectra file %s",
                       converted_spectra_file.c_str());
    }
  } else {
    // Successfully read file as spectrumrecords format
    if (get_int_parameter("remove_precursor_peak") != 0) {
      carp(CARP_FATAL, "remove_precursor_peak can only be used during conversion "
                       "to spectrumrecords format.");
    }
  }

  carp(CARP_INFO, "Sorting spectra");
  if (window_type[0] != WINDOW_MZ) {
    spectra.Sort();
  } else {
    spectra.Sort<ScSortByMz>(ScSortByMz(window));
  }
  if (max_mz == 0) {
    double highest_mz = spectra.FindHighestMZ();
    carp(CARP_DEBUG, "Max m/z %f", highest_mz);
    MaxMZ::SetGlobalMax(highest_mz);
  } else {
    MaxMZ::SetGlobalMaxFromFlag();
  }

  char* digestString =
    digest_type_to_string(pepHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);

  bool txt_only = !get_boolean_parameter("sqt-output") &&
                  !get_boolean_parameter("pepxml-output") &&
                  !get_boolean_parameter("mzid-output") &&
                  !get_boolean_parameter("pinxml-output");
  OutputFiles* output_files = NULL;
  ofstream* target_file[NUM_THREADS];
  ofstream* decoy_file[NUM_THREADS];

  if (!txt_only) {
    carp(CARP_DEBUG, "Using OutputFiles to write matches");
    // Overwrite enzyme/digestion parameters in the hash
    // TODO Find a better way to do this?
    add_or_update_hash(parameters, "enzyme", pepHeader.enzyme().c_str());
    add_or_update_hash(parameters, "digestion", digestString);
    free(digestString);
    output_files = new OutputFiles(this);
  } else {
    carp(CARP_DEBUG, "Using TideMatchSet to write matches");
    bool overwrite = get_boolean_parameter("overwrite");
    stringstream ss;
    ss << pepHeader.enzyme() << '-' << digestString;
    free(digestString);
    TideMatchSet::setCleavageType(ss.str());
    if (!concat) {
      string target_file_name = make_file_path("tide-search.target.txt");
      target_file[0] = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
      for (int i = 1; i < NUM_THREADS; i++) {
        ostringstream target_file_name_string;
        target_file_name_string << "tide-search.target" << i << ".txt";
        target_file_name = make_file_path(target_file_name_string.str());
        target_file[i] = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
      }

      if (HAS_DECOYS) {
        string decoy_file_name = make_file_path("tide-search.decoy.txt");
        decoy_file[0] = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
        for (int i = 1; i < NUM_THREADS; i++) {
          ostringstream decoy_file_name_string;
          decoy_file_name_string << "tide-search.decoy" << i << ".txt";
          decoy_file_name = make_file_path(decoy_file_name_string.str());
          decoy_file[i] = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
        }
      }
    } else {
      string concat_file_name = make_file_path("tide-search.txt");
      target_file[0] = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
    }
  }

  // Do the search
  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6); // may need utils.cpp ?
  carp(CARP_INFO, "Running search");
  cleanMods();
  search(spectra.SpecCharges(), active_peptide_queue, proteins, locations,
         window, window_type, get_double_parameter("spectrum-min-mz"),
         get_double_parameter("spectrum-max-mz"), min_scan, max_scan,
         get_int_parameter("min-peaks"), charge_to_search,
         get_int_parameter("top-match"), spectra.FindHighestMZ(),
         output_files, target_file, decoy_file, compute_sp);

  // Delete temporary spectrumrecords file
  if (!delete_spectra_file.empty()) {
    carp(CARP_DEBUG, "Deleting %s", delete_spectra_file.c_str());
    remove(delete_spectra_file.c_str());
  }

  // Clean up
  for (int i = 0; i < NUM_THREADS; i++) {
    delete peptide_reader[i];
    delete active_peptide_queue[i];
  }

  for (ProteinVec::const_iterator i = proteins.begin(); i != proteins.end(); ++i) {
    delete const_cast<pb::Protein*>(*i);
  }
  if (output_files) {
    delete output_files;
  }
  if (target_file[0]) {
    for (int i = 0; i < NUM_THREADS; i++) {
      target_file[i]->close();
      delete target_file[i];
    }

    if (decoy_file[0]) {
      for (int i = 0; i < NUM_THREADS; i++) {
        decoy_file[i]->close();
        delete decoy_file[i];
      }
    }
  }

  if (!concat) {
    string target_file_name = make_file_path("tide-search.target.txt");
    ofstream target_file(target_file_name.c_str(), ios::out | ios::app);

    for (int i = 1; i < NUM_THREADS; i++) {
      ostringstream target_file_name_string;
      target_file_name_string << "tide-search.target" << i << ".txt";
      target_file_name = make_file_path(target_file_name_string.str());
      ifstream target_file2(target_file_name.c_str(), ios::in);

      if (target_file.is_open()) {
        target_file << target_file2.rdbuf();
      }
      int remtarget = remove(target_file_name.c_str());
    }
    target_file.close();

    if (HAS_DECOYS) {
      string decoy_file_name = make_file_path("tide-search.decoy.txt");
      ofstream decoy_file(decoy_file_name.c_str(), ios::out | ios::app);

      for (int i = 1; i < NUM_THREADS; i++) {
        ostringstream decoy_file_name_string;
        decoy_file_name_string << "tide-search.decoy" << i << ".txt";
        decoy_file_name = make_file_path(decoy_file_name_string.str());
        ifstream decoy_file2(decoy_file_name.c_str(), ios::in);
        if (decoy_file.is_open()) {
          decoy_file << decoy_file2.rdbuf();
        }
        int remdecoy = remove(decoy_file_name.c_str());
      }
      decoy_file.close();
    }
  } else {
        // Concat option ?
//      string concat_file_name = make_file_path("tide-search.txt");
//      target_file = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
  }

  return 0;
}

/**
 * Free all existing mods
 */
void TideSearchApplication::cleanMods() {
  for (int i = 0; i < MAX_AA_MODS; ++i) {
    free_aa_mod(list_of_mods[i]);
    list_of_mods[i] = NULL;
  }
  num_mods = 0;
}

// method for each spec charge vector segment
void * TideSearchApplication::searchthread(void* threadarg)
{
  struct thread_data *my_data = (struct thread_data *) threadarg;

  const vector<SpectrumCollection::SpecCharge>* spec_charges = my_data->spec_charges;
  ActivePeptideQueue* active_peptide_queue = my_data->active_peptide_queue;
  ProteinVec& proteins = my_data->proteins;
  vector<const pb::AuxLocation*>& locations = my_data->locations;
  double precursor_window = my_data->precursor_window;
  WINDOW_TYPE_T window_type = my_data->window_type;
  double spectrum_min_mz = my_data->spectrum_min_mz;
  double spectrum_max_mz = my_data->spectrum_max_mz;
  int min_scan = my_data->min_scan;
  int max_scan = my_data->max_scan;
  int min_peaks = my_data->min_peaks;
  int search_charge = my_data->search_charge;
  int top_matches = my_data->top_matches;
  double highest_mz = my_data->highest_mz;
  OutputFiles* output_files = my_data->output_files;
  ofstream* target_file = my_data -> target_file;
  ofstream* decoy_file = my_data->decoy_file;
  bool compute_sp = my_data->compute_sp;
  long lo = my_data->lo;
  long hi = my_data->hi;

  // This is the main search loop.
  ObservedPeakSet observed;
  // cycle through spectrum-charge pairs, sorted by neutral mass
  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin()+lo;
       sc != spec_charges->begin()+hi;
       ++sc) {
//  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();
//       sc != spec_charges->end();
//       ++sc) {
    const Spectrum* spectrum = sc->spectrum;

    double precursor_mz = spectrum->PrecursorMZ();
    int charge = sc->charge;
    int scan_num = spectrum->SpectrumNumber();

    if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz ||
        scan_num < min_scan || scan_num > max_scan ||
        spectrum->Size() < min_peaks ||
        (search_charge != 0 && charge != search_charge)) {
      continue;
    }

    // Normalize the observed spectrum and compute the cache of
    // frequently-needed values for taking dot products with theoretical
    // spectra.
    observed.PreprocessSpectrum(*spectrum, charge);

    // The active peptide queue holds the candidate peptides for spectrum.
    // Calculate and set the window, depending on the window type.
    double min_mass, max_mass;
    computeWindow(*sc, window_type, precursor_window, &min_mass, &max_mass);

    int size = active_peptide_queue->SetActiveRange(min_mass, max_mass);
    TideMatchSet::Arr match_arr(size); // Scored peptides will go here.

    // Programs for taking the dot-product with the observed spectrum are laid
    // out in memory managed by the active_peptide_queue, one program for each
    // candidate peptide. The programs will store the results directly into
    // match_arr. We now pass control to those programs.
    collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr,
                          size, charge);

    // matches will arrange the results in a heap by score, return the top
    // few, and recover the association between counter and peptide. We output
    // the top matches.
    TideMatchSet matches(&match_arr, highest_mz);

    if (output_files) {
      matches.report(output_files, top_matches, spectrum, charge,
                     active_peptide_queue, proteins, locations, compute_sp);
    } else {
      matches.report(target_file, decoy_file, top_matches, spectrum, charge,
                     active_peptide_queue, proteins, locations, compute_sp);
    }
  }
  pthread_exit(NULL);
}

void * searchthread_helper(void * threadarg)
    {
        return ((TideSearchApplication *)threadarg)->searchthread(threadarg);
    }

void TideSearchApplication::search(
  const vector<SpectrumCollection::SpecCharge>* spec_charges,
  ActivePeptideQueue * active_peptide_queue[],
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
) {

  if (output_files) {
    output_files->writeHeaders();
  } else if (target_file[0]) {
    TideMatchSet::writeHeaders(target_file[0], false, compute_sp);
    TideMatchSet::writeHeaders(decoy_file[0], true, compute_sp);
  }

  // Attempting to parallelize main search loop...
  
  pthread_t threads[NUM_THREADS - 1];
  struct thread_data thread_data_array[NUM_THREADS - 1];
  for (int i= 0; i < NUM_THREADS - 1; i++) {
    thread_data_array[i].setParams(spec_charges, active_peptide_queue[i+1],
      proteins, locations, precursor_window, window_type[i+1], spectrum_min_mz,
      spectrum_max_mz, min_scan, max_scan, min_peaks, search_charge, top_matches,
      highest_mz, output_files, target_file[i+1], decoy_file[i+1], compute_sp,
      ((i+1) * spec_charges->size()) / NUM_THREADS, ((i+2) * spec_charges->size()) / NUM_THREADS);
  }

  int rc;
  for (long t = 0; t < NUM_THREADS - 1; t++) {
    rc = pthread_create(&threads[t], NULL, &searchthread_helper, 
        (void *) &(thread_data_array[t]));
  }

  ObservedPeakSet observed;
  // cycle through spectrum-charge pairs, sorted by neutral mass
  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();
       sc != spec_charges->begin() + (spec_charges->size()/NUM_THREADS);
       ++sc) {
    const Spectrum* spectrum = sc->spectrum;

    double precursor_mz = spectrum->PrecursorMZ();
    int charge = sc->charge;
    int scan_num = spectrum->SpectrumNumber();

    if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz ||
        scan_num < min_scan || scan_num > max_scan ||
        spectrum->Size() < min_peaks ||
        (search_charge != 0 && charge != search_charge)) {
      continue;
    }

    // Normalize the observed spectrum and compute the cache of
    // frequently-needed values for taking dot products with theoretical
    // spectra.
    observed.PreprocessSpectrum(*spectrum, charge);

    // The active peptide queue holds the candidate peptides for spectrum.
    // Calculate and set the window, depending on the window type.
    double min_mass, max_mass;
    computeWindow(*sc, window_type[0], precursor_window, &min_mass, &max_mass);

    int size = active_peptide_queue[0]->SetActiveRange(min_mass, max_mass);
    TideMatchSet::Arr match_arr(size); // Scored peptides will go here.

    // Programs for taking the dot-product with the observed spectrum are laid
    // out in memory managed by the active_peptide_queue, one program for each
    // candidate peptide. The programs will store the results directly into
    // match_arr. We now pass control to those programs.
    collectScoresCompiled(active_peptide_queue[0], spectrum, observed, &match_arr,
                          size, charge);

    // matches will arrange the results in a heap by score, return the top
    // few, and recover the association between counter and peptide. We output
    // the top matches.
    TideMatchSet matches(&match_arr, highest_mz);
    if (output_files) {
      matches.report(output_files, top_matches, spectrum, charge,
                     active_peptide_queue[0], proteins, locations, compute_sp);
    } else {
      matches.report(target_file[0], decoy_file[0], top_matches, spectrum, charge,
                     active_peptide_queue[0], proteins, locations, compute_sp);
    }
  }
  for (int i = 0; i < NUM_THREADS - 1; i++) {
    pthread_join(threads[i], NULL);
  }

/*
  // This is the main search loop.
  ObservedPeakSet observed;
  // cycle through spectrum-charge pairs, sorted by neutral mass
  for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();
       sc != spec_charges->end();
       ++sc) {
    const Spectrum* spectrum = sc->spectrum;

    double precursor_mz = spectrum->PrecursorMZ();
    int charge = sc->charge;
    int scan_num = spectrum->SpectrumNumber();

    if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz ||
        scan_num < min_scan || scan_num > max_scan ||
        spectrum->Size() < min_peaks ||
        (search_charge != 0 && charge != search_charge)) {
      continue;
    }

    // Normalize the observed spectrum and compute the cache of
    // frequently-needed values for taking dot products with theoretical
    // spectra.
    observed.PreprocessSpectrum(*spectrum, charge);

    // The active peptide queue holds the candidate peptides for spectrum.
    // Calculate and set the window, depending on the window type.
    double min_mass, max_mass;
    computeWindow(*sc, window_type, precursor_window, &min_mass, &max_mass);

    int size = active_peptide_queue->SetActiveRange(min_mass, max_mass);
    TideMatchSet::Arr match_arr(size); // Scored peptides will go here.

    // Programs for taking the dot-product with the observed spectrum are laid
    // out in memory managed by the active_peptide_queue, one program for each
    // candidate peptide. The programs will store the results directly into
    // match_arr. We now pass control to those programs.
    collectScoresCompiled(active_peptide_queue, spectrum, observed, &match_arr,
                          size, charge);

    // matches will arrange the results in a heap by score, return the top
    // few, and recover the association between counter and peptide. We output
    // the top matches.
    TideMatchSet matches(&match_arr, highest_mz);
    if (output_files) {
      matches.report(output_files, top_matches, spectrum, charge,
                     active_peptide_queue, proteins, locations, compute_sp);
    } else {
      matches.report(target_file, decoy_file, top_matches, spectrum, charge,
                     active_peptide_queue, proteins, locations, compute_sp);
    }
  }
*/
  if (output_files) {
    output_files->writeFooters();
  }
}

void TideSearchApplication::collectScoresCompiled(
  ActivePeptideQueue* active_peptide_queue,
  const Spectrum* spectrum,
  const ObservedPeakSet& observed,
  TideMatchSet::Arr* match_arr,
  int queue_size,
  int charge
) {
  if (!active_peptide_queue->HasNext()) {
    return;
  }
  // prog gets the address of the dot-product program for the first peptide
  // in the active queue.
  const void* prog = active_peptide_queue->NextPeptide()->Prog(charge);
  const int* cache = observed.GetCache();
  // results will get (score, counter) pairs, where score is the dot product
  // of the observed peak set with a candidate peptide. The candidate
  // peptide is given by counter, which refers to the index within the
  // ActivePeptideQueue, counting from the back. This complication
  // simplifies the generated programs, which now simply dump the counter.
  pair<int, int>* results = match_arr->data();

  // See compiler.h for a description of the programs beginning at prog and
  // how they are generated. Here we initialize certain registers to the
  // values expected by the programs and call the first one (*prog).
  //
  // See gnu assembler format for more on this format. We tell the compiler
  // to set these registers:
  // edx/rdx points to the cache.
  // eax/rax points to the first program.
  // ecx/rcx is the counter and gets the size of the active queue.
  // edi/rdi points to the results buffer.
  //
  // The push and pop operations are a workaround for a compiler that
  // doesn't understand that %ecx and %edi (or %rcx and %rdi) get
  // clobbered. Since they're already input registers, they can't be
  // included in the clobber list.
  __asm__ __volatile__("cld\n" // stos operations increment edi
#ifdef __x86_64__
                       "push %%rcx\n"
                       "push %%rdi\n"
                       "call *%%rax\n"
                       "pop %%rdi\n"
                       "pop %%rcx\n"
#else
                       "push %%ecx\n"
                       "push %%edi\n"
                       "call *%%eax\n"
                       "pop %%edi\n"
                       "pop %%ecx\n"
#endif
                       : // no outputs
                       : "d" (cache),
                         "a" (prog),
                         "c" (queue_size),
                         "D" (results)
  );

  // match_arr is filled by the compiled programs, not by calls to
  // push_back(). We have to set the final size explicitly.
  match_arr->set_size(queue_size);
}

void TideSearchApplication::computeWindow(
  const SpectrumCollection::SpecCharge& sc,
  WINDOW_TYPE_T window_type,
  double precursor_window,
  double* out_min,
  double* out_max
) {
  switch (window_type) {
  case WINDOW_MASS:
    *out_min = sc.neutral_mass - precursor_window;
    *out_max = sc.neutral_mass + precursor_window;
    break;
  case WINDOW_MZ: {
    double mz_minus_proton = sc.spectrum->PrecursorMZ() - MASS_PROTON;
    *out_min = (mz_minus_proton - precursor_window) * sc.charge;
    *out_max = (mz_minus_proton + precursor_window) * sc.charge;
    break;
    }
  case WINDOW_PPM: {
    double tiny_precursor = precursor_window * 1e-6;
    *out_min = sc.neutral_mass * (1.0 - tiny_precursor);
    *out_max = sc.neutral_mass * (1.0 + tiny_precursor);
    break;
    }
  default:
    carp(CARP_FATAL, "Invalid window type");
  }
  carp(CARP_DETAILED_DEBUG, "Scan %d.%d mass window is [%f, %f]",
       sc.spectrum->SpectrumNumber(), sc.charge, *out_min, *out_max);
}

bool TideSearchApplication::hasDecoys() {
  return HAS_DECOYS;
}

string TideSearchApplication::getName() {
  return "tide-search";
}

string TideSearchApplication::getDescription() {
  return
  "Search a collection of spectra against a sequence "
  "database, returning a collection of peptide-spectrum "
  "matches (PSMs) scored by XCorr.";
}

bool TideSearchApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideSearchApplication::getCommand() {
  return TIDE_SEARCH_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
