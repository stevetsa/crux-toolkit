#include "spit.h"

#define NUM_SPIT_OPTIONS 8
#define NUM_SPIT_ARGUMENTS 1

// crux spit (simple protein identification tool)
// given a colection of scored PSMs, produces a
// ranked list of proteins


using namespace std;

// A map from peptides to their scores
typedef std::map<std::string, FLOAT_T> PeptideScore;

// A map of protein names to all peptides it contains
typedef std::map<std::string, std::set<std::string> > ProteinPeptides;

// list of all pairs of proteins and their scores
typedef std::vector< std::pair<FLOAT_T, std::string> > ProteinScore;

/* private method headers */

// prints tab-delimited file of proteins and their
// scores in sorted order
bool print_spit_scores (
	const ProteinScore&,
	char* file 		// output file
	);

// assigns a score to each protein
void get_spit_scores(
	ProteinScore&,
	const PeptideScore&,
	ProteinPeptides&,
	map<string, int>&,
        bool
	);

// post process to create unique protein-peptide links
// by removing peptides from lower scored proteins for filter option.
// recalculates spit scores afterwards.
void filter_proteins(ProteinScore&,
	ProteinPeptides&);

// gets top scoring psms from a txt file created by
// crux search-for-matches
bool get_txt_matches(
	PeptideScore&,
	ProteinPeptides&,
	char* psm_folder
	);

// gets top scoring psms from a database
bool get_database_matches(
	PeptideScore&,
	ProteinPeptides&,
	map<string, int>&,
	char* psm_folder,
	char* database
 	);

// splits a string, returning a vector of tokens
void string_split(const std::string&, std::vector<std::string>&, const std::string& delimiters = std::string(" "));

int spit_main(int argc, char** argv) {
//int main(int argc, char** argv) {

  // Define optional and required command line arguments
  int num_options = NUM_SPIT_OPTIONS;
  char* option_list[NUM_SPIT_OPTIONS] = {
    "database",
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity",
    "filter",
    "version"};

  int num_arguments = NUM_SPIT_ARGUMENTS;
  char* argument_list[NUM_SPIT_ARGUMENTS] = {
    "psm-folder"
  };

  initialize_parameters();

  // Define optional command line arguments in parameter.c 
  select_cmd_line_arguments(argument_list, num_arguments);
  select_cmd_line_options(option_list, num_options );

  // Parse the command line and optional parameter file
  // does sytnax, type, and bounds checking and dies on error
  parse_cmd_line_into_params_hash(argc, argv, "spit");

  char* psm_dir = get_string_parameter("psm-folder");
  char* database = get_string_parameter("database");
  char* output_dir = get_string_parameter("output-dir");
  char* output_file = get_string_parameter("spit-output-file");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  BOOLEAN_T filter = get_boolean_parameter("filter");

  prefix_fileroot_to_name(&output_file);

  char* full_output_file = get_full_filename(output_dir, output_file); 
  create_output_directory(output_dir, TRUE, FALSE);
  create_file_in_path(output_file, output_dir, overwrite);
  
  carp(CARP_DETAILED_DEBUG, "psm-folder: %s", psm_dir);
  carp(CARP_DETAILED_DEBUG, "database: %s", database);
  
  // initialize peptide-score and protein-peptides maps
  PeptideScore peptideScore;
  ProteinPeptides proteinPeptides;
  ProteinScore proteinScores;
  map<string, int> numPeptides; // keeps counts of peptides for database option

  // get matches from .txt file or database
  if (database) {
    if (!get_database_matches(peptideScore, proteinPeptides, numPeptides, psm_dir, database)) {
      carp(CARP_ERROR, "error getting matches from database"); 
      exit(1);
    }
  } else {
    if (!get_txt_matches(peptideScore, proteinPeptides, psm_dir)) {
      carp(CARP_ERROR, "error getting matches from .txt file");
      exit(1);
    }
  }

  // calculate final protein scores
  get_spit_scores(proteinScores, peptideScore, proteinPeptides, numPeptides, filter);

  // print out scores in crux output directory
  if (!print_spit_scores(proteinScores, full_output_file)) {
    carp(CARP_ERROR, "error printing protein scores");
    exit(1);
  }
  carp(CARP_INFO, "completed succesfully");
  return 0;
}

// does post-processing on the proteins. 
void filter_proteins(ProteinScore& proteinScores, ProteinPeptides& proteinPeptides) {
  carp(CARP_INFO, "filtering protein-peptide links");
  set<string> top_peptides;
  for (ProteinScore::iterator score_pair = proteinScores.begin(); score_pair != proteinScores.end(); ++score_pair) {
    for (set<string>::iterator peptide = proteinPeptides[score_pair->second].begin(); 
	                       peptide != proteinPeptides[score_pair->second].end(); ++peptide) {
      // if a unique peptide so far
      if (top_peptides.find(*peptide) == top_peptides.end()) {
	top_peptides.insert(*peptide);
      } else {
        proteinPeptides[score_pair->second].erase(peptide);
      }
    }
  } 
}


// takes a map of peptide scores and map from proteins to peptides and
// stores a list of pairs of protein IDs and scores
void get_spit_scores (ProteinScore& proteinScores, const PeptideScore& peptideScore, ProteinPeptides& proteinPeptides, map<string, int>& numPeptides, bool filter) {
  carp(CARP_INFO, "computing protein scores");
  
  char* database = get_string_parameter("database");
 
  for (ProteinPeptides::const_iterator protein = proteinPeptides.begin(); protein != proteinPeptides.end(); ++protein) { 
    // calculate score
    FLOAT_T score = 1.0;
    set<string> peptides = protein->second;
    // for each peptide belonging to the protein
    for (set<string>::const_iterator peptide = peptides.begin(); peptide != peptides.end(); ++peptide) {
      PeptideScore::const_iterator it = peptideScore.find(*peptide);
      // multiply score by peptides score
      score *= it->second;
    }
    // take the nth root. if using a database, n is total number of peptides found in the database,
    // otherwise n is number of peptides found in the text file.
    if (database) {
      string id = protein->first;
      score = pow(score, (1.0 / numPeptides[protein->first.c_str()]));
    } else {
      score = pow(score, (1.0 / peptides.size()));
    }
    pair<FLOAT_T, string> score_pair (score, protein->first);
    proteinScores.push_back(score_pair);
  }
  // sort the pairs by their scores
  sort(proteinScores.begin(), proteinScores.end()); 
  reverse(proteinScores.begin(), proteinScores.end());
  if (filter) {
    filter_proteins(proteinScores, proteinPeptides);
    ProteinScore new_proteinScores;
    // recalculate spit scores, skipping filter step
    get_spit_scores(new_proteinScores, peptideScore, proteinPeptides, numPeptides, false); 
  }
}

// prints the tab-delimited file of protein ids and their
// scores to the target file
bool print_spit_scores(const ProteinScore& proteinScores, char* file) {
  ofstream output_stream(file);
  if (!output_stream.is_open()) {
    carp(CARP_ERROR, "error opening output file");
    return false;
  }
  // for each pair in proteinScores
  for (ProteinScore::const_iterator score_pair = proteinScores.begin(); score_pair != proteinScores.end(); ++score_pair) {
    output_stream << score_pair->second << "\t" << score_pair->first << endl;
  } 
  return true;
}

// gets all the spit matches from a protein database
// and counts number of peptides for each parent protein for scoring
// returns true on success 
bool get_database_matches(
	PeptideScore& peptideScore,
	ProteinPeptides& proteinPeptides,
	map<string, int>& numPeptides,
	char* psm_dir,
	char* database
 	) {

  // declarations
  PEPTIDE_SRC_ITERATOR_T* src_iterator = NULL;
  PROTEIN_T* protein = NULL;  
  PEPTIDE_SRC_T* peptide_src; 
  PEPTIDE_CONSTRAINT_T* peptide_constraint = NULL;
  PROTEIN_PEPTIDE_ITERATOR_T* peptide_iterator = NULL;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_T* match = NULL;
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator = 
	new_match_collection_iterator(psm_dir, database);
  PEPTIDE_T* peptide;
  
  // string representations
  char* peptide_sequence;
  char* protein_id;

  carp(CARP_INFO, "reading matches from %s", psm_dir);
  while(match_collection_iterator_has_next(match_collection_iterator)) {
    peptide_constraint = new_peptide_constraint_from_parameters();
    // get the next match_collection
    match_collection =
      match_collection_iterator_next(match_collection_iterator);
   
    // create iterator
    match_iterator = new_match_iterator(match_collection, LOGP_BONF_WEIBULL_XCORR, FALSE);

    // for each match    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      // skip if not the best match 
      if (get_match_rank(match, LOGP_BONF_WEIBULL_XCORR) != 1) {
	continue;
      }	
      // get p-value
      FLOAT_T score = get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
      // ignore unscored psms
      if (score == P_VALUE_NA) {
        carp(CARP_INFO, "skipping unscored psm");
        continue;
      }

      // get peptide and parent protein information from match
      peptide = get_match_peptide(match); 
      peptide_sequence = get_match_sequence_sqt(match);

      src_iterator = new_peptide_src_iterator(peptide);
      // iterate over all parent proteins
      while (peptide_src_iterator_has_next(src_iterator)) {
        
	// get the id for the protein
        peptide_src = peptide_src_iterator_next(src_iterator);
        protein = get_peptide_src_parent_protein(peptide_src);
	protein_id = get_protein_id_pointer(protein);

        // if peptides haven't been counted yet for the protein
	if (!numPeptides[protein_id]) {
	  peptide_iterator = new_protein_peptide_iterator(protein, peptide_constraint);
          int peptide_count = 0;
          // count the peptides
          while (protein_peptide_iterator_has_next(peptide_iterator)) {
	    protein_peptide_iterator_next(peptide_iterator);
	    peptide_count++;
          }
 	  numPeptides[protein_id] = peptide_count;
          free_protein_peptide_iterator(peptide_iterator);
        }
	PeptideScore::iterator it = peptideScore.find(peptide_sequence);
      	// if a new peptide, or a larger score, update peptide score
      	if (it == peptideScore.end() || (it->second < score)) {
	  peptideScore[peptide_sequence] = score;
      	}	
	// add the peptide to the current protein
        proteinPeptides[protein_id].insert(peptide_sequence);
      } // get next peptide
      free_peptide_src_iterator(src_iterator);
      carp(CARP_DEBUG, "num proteins: %d", get_match_collection_num_proteins(match_collection));
    } // get the next match
    free_match_iterator(match_iterator); 
  } // get the next match collection
  free_match_collection_iterator(match_collection_iterator);
  return true;
}


// gets all spit matches from a txt file created by search-for-matches
// returns true on success
bool get_txt_matches(PeptideScore& peptideScore, ProteinPeptides& proteinPeptides, char* psm_dir) {
 
  // file handling declarations
  DIR *dir = opendir(psm_dir);
  struct dirent* directory_entry = NULL;
  char* file;
  ifstream txt_file;
  string line;
  vector<string> tokens;

  int scan_number;
  int next_scan;

  FLOAT_T pvalue;

  // indeces for each psm
  int scan_index, pvalue_index, sequence_index, protein_index;

  while ((directory_entry = readdir(dir)) != NULL) {
    if (!suffix_compare(directory_entry->d_name, "search.target.txt")) {
      continue;
    }
    
    file = get_full_filename(psm_dir, directory_entry->d_name);

    carp(CARP_INFO, "reading matches from %s", file); 
    // open the file stream 
    txt_file.open(file, ios::in);

    scan_index = -1;
    pvalue_index = -1;
    sequence_index = -1;
    protein_index = -1;

    // get header and find indeces 
    tokens.clear();
    getline(txt_file, line);
    string_split(line, tokens, "\t");
    for (int i = 0; i < tokens.size(); ++i) {
      if (tokens[i] == "-log(p-value)") pvalue_index = i;
      else if (tokens[i] == "sequence") sequence_index = i;
      else if (tokens[i] == "protein id") protein_index = i;
      else if (tokens[i] == "scan") scan_index = i;
    }

    carp(CARP_DETAILED_DEBUG, "header indeces: pvalue %d sequence %d protein %d",
			          pvalue_index, sequence_index, protein_index);
    if (pvalue_index   < 0
    ||  sequence_index < 0
    ||  protein_index  < 0
    ||  scan_index     < 0) {
      carp(CARP_INFO, "found invalid .txt peptide-spectrum match file, skipping");
      continue;
    }
  
    scan_number = 0;
    // for each psm in file
    while (getline(txt_file, line)) {
      tokens.clear();
      string_split(line, tokens, "\t");
      next_scan = atoi(tokens[scan_index].c_str());

      // skip if not the best match
      if (scan_number == next_scan) continue;
      else scan_number = next_scan;
   
      carp(CARP_DETAILED_DEBUG, "scan number %d", scan_number);

      pvalue = atof(tokens[pvalue_index].c_str());
      // skip if nan or unscored
      if (pvalue == 0 || pvalue != pvalue) {
        carp(CARP_INFO, "No p-value found for scan %d, skipping", scan_number);
        continue;
      }
   
      PeptideScore::iterator it = peptideScore.find(tokens[sequence_index]);
      // if a new peptide, or a larger score, add the peptide and its score to the map
      if (it == peptideScore.end() || (it->second < pvalue)) {
	peptideScore[tokens[sequence_index]] = pvalue;
      }	
      // split the protein string by commas
      vector<string> parent_proteins;
      string_split(tokens[protein_index], parent_proteins, ",");
      // for each parent protein, add the peptide 
      for (vector<string>::const_iterator protein = parent_proteins.begin(); protein != parent_proteins.end(); ++protein) {
	proteinPeptides[*protein].insert(tokens[sequence_index]);	
      }      
    } // get next psm
    txt_file.close();
  } // get next txt file
  closedir(dir);
  // die if no p-values found
  if (peptideScore.size() == 0) {
    carp(CARP_ERROR, "No p-values for peptide-spectrum matches found.");
    carp(CARP_ERROR, "spit requires search-for-matches to be run with --compute-p-values option");
    return false;
  }
  return true;
}

// splits a string into a vector of tokens
void string_split(const string& str, vector<string>& tokens, const string& delimiters) {
    // skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos)
    {
        // found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
 	if (pos != string::npos) lastPos = pos+1;
        else lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
