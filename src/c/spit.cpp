#include "spit.h"

#define NUM_SPIT_OPTIONS 7
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
typedef std::vector< std::pair<FLOAT_T, std::string> > ProteinScores;

bool print_spit_scores (
	const ProteinScores&,
	char* file 
	);

void get_spit_scores(
	ProteinScores&,
	const PeptideScore&,
	const ProteinPeptides&
	);

bool get_txt_matches(
	PeptideScore&,
	ProteinPeptides&,
	char* psm_folder
	);

bool get_database_matches(
	PeptideScore&,
	ProteinPeptides&,
	char* psm_folder,
	char* database
 	);

// splits a string, returning a vector of tokens
void string_split(const std::string&, std::vector<std::string>&, const std::string& delimiters = std::string(" "));

//int spit_main(int argc, char** argv) {
int main(int argc, char** argv) {

  // Define optional and required command line arguments
  int num_options = NUM_SPIT_OPTIONS;
  char* option_list[NUM_SPIT_OPTIONS] = {
    "database",
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity",
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
  char* output_dir = get_string_parameter("output-dir");
  char* database = get_string_parameter("database");
  
  // get the output file
  char* output_file = "proteins.target.txt";
  prefix_fileroot_to_name(&output_file);
  
  //set_verbosity_level(CARP_DETAILED_DEBUG);
  //set_verbosity_level(CARP_ERROR);
  //set_verbosity_level(CARP_INFO);

  carp(CARP_DETAILED_DEBUG, "psm-folder: %s", psm_dir);
  carp(CARP_DETAILED_DEBUG, "output-dir: %s", output_dir);
  carp(CARP_DETAILED_DEBUG, "database: %s", database);
  
  // initialized peptide-score and protein-peptides maps
  PeptideScore peptideScore;
  ProteinPeptides proteinPeptides;
  
  ProteinScores proteinScores;
  
  if (database) { 
    carp(CARP_ERROR, "database option not working (yet)");
    exit(1);
  }

  // get matches from .txt file or database
  if (database) {
    if (!get_database_matches(peptideScore, proteinPeptides, psm_dir, database)) {
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
  get_spit_scores(proteinScores, peptideScore, proteinPeptides);

  // print out scores to proteins.txt in output directory
  if (!print_spit_scores(proteinScores, output_file)) {
    carp(CARP_ERROR, "error printing protein scores");
    exit(1);
  }
  return 0;
}


// gets all the spit matches from a database
// returns TRUE on success 
bool get_database_matches(
	PeptideScore& peptideScore,
	ProteinPeptides& proteinPeptides,
	char* psm_dir,
	char* database_file
 	) {

  // declarations
  //PROTEIN_PEPTIDE_ITERATOR_T* peptide_iterator = NULL;
  PEPTIDE_SRC_ITERATOR_T* src_iterator = NULL;
  PROTEIN_T* protein = NULL;  
  PEPTIDE_SRC_T* peptide_src;
  PEPTIDE_CONSTRAINT_T* peptide_constraint = NULL;
  DATABASE_PEPTIDE_ITERATOR_T* peptide_iterator = NULL;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_T* match = NULL;
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator = 
	new_match_collection_iterator(psm_dir, database_file);
  PEPTIDE_T* peptide;
  PEPTIDE_T* match_peptide;

  //TODO: how to check if memory mapped database?
  DATABASE_T* database = new_database(database_file, FALSE);
  
  char* parent_proteins;
  char* peptide_sequence;

  while(match_collection_iterator_has_next(match_collection_iterator)) {
    
    // get the next match_collection
    match_collection =
      match_collection_iterator_next(match_collection_iterator);
   
    // create iterator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);
	
    // for each match    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      
      // get p-value
      FLOAT_T score = get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
      // ignore unscored psms
      if (score == P_VALUE_NA) {
        carp(CARP_INFO, "skipping unscored psm");
        continue;
      }
      // get sequence and parent information from match
      match_peptide = get_match_peptide(match); 
      //parent_proteins = get_protein_ids(match_peptide);
      peptide_sequence = get_match_sequence_sqt(match);
      
      peptide_constraint = new_peptide_constraint_from_parameters();

      src_iterator = new_peptide_src_iterator(match_peptide);
    
      // iterate over all parent proteins
      while (peptide_src_iterator_has_next(src_iterator)) {
        peptide_src = peptide_src_iterator_next(src_iterator);
        protein = get_peptide_src_parent_protein(peptide_src);
        peptide_iterator = new_database_peptide_iterator(database, peptide_constraint);
      // iterate over all peptides
      while (database_peptide_iterator_has_next(peptide_iterator)) {
	// if peptide not found yet, add it to the map
	peptide = database_peptide_iterator_next(peptide_iterator);
	//cout << get_peptide_sequence(match_peptide) << endl;
	//cout << get_peptide_sequence(peptide) << endl;
	if (get_peptide_sequence(match_peptide) == get_peptide_sequence(peptide))
        cout << "matched" << endl;
	//peptideScore[

      }
      free_database_peptide_iterator(peptide_iterator);
    }
    free_peptide_src_iterator(src_iterator);
     // peptide_iterator = new_protein_peptide_iterator(PROTEIN_T* protein, 
/*    
    src_iterator = new_peptide_src_iterator(match->peptide);
    // iterate overall parent proteins
    // find largest numProt and pepSite among the parent proteins
    while(peptide_src_iterator_has_next(src_iterator)){
      peptide_src = peptide_src_iterator_next(src_iterator);
      protein = get_peptide_src_parent_protein(peptide_src);
      protein_idx = get_protein_protein_idx(protein);
    } 
*/
  
     /* 
      PeptideScore::iterator it = peptideScore.find(fields[sequence_index]);
      // if a new peptide, or a larger score, add the peptide and its score to the map
      if (it == peptideScore.end() || (it->second < atof(fields[pvalue_index]))) {
	peptideScore[fields[sequence_index]] = atof(fields[pvalue_index]);
      }	
      // split the protein string by commas
      vector<string> parent_proteins;
      string_split(string(fields[protein_index]), parent_proteins, ",");
      // for each protein
      for (vector<string>::const_iterator protein = parent_proteins.begin(); protein != parent_proteins.end(); ++protein) {
	proteinPeptides[*protein].insert(fields[sequence_index]);	
      }
      */
      carp(CARP_DEBUG, "num proteins: %d", get_match_collection_num_proteins(match_collection));
      //free peptide_str 
      // free parent_proteins	
    } // get the next match
    free_match_iterator(match_iterator); 
  } // get the next match_collection
  free_match_collection_iterator(match_collection_iterator);
  return true;
}

bool print_spit_scores(const ProteinScores& proteinScores, char* file) {
  ofstream output_file(file);
  if (!output_file.is_open()) {
    carp(CARP_ERROR, "error opening output file");
    return false;
  }
  // for each pair in proteinScores
  for (ProteinScores::const_iterator score_pair = proteinScores.begin(); score_pair != proteinScores.end(); ++score_pair) {
    output_file << score_pair->second << "\t" << score_pair->first << endl;
  } 
  return true;
}

// takes a map of peptide scores and map from proteins to peptides and
// stores a list of pairs of protein names and scores
void get_spit_scores (ProteinScores& proteinScores, const PeptideScore& peptideScore, const ProteinPeptides& proteinPeptides) {
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
    // take the nth root
    score = pow(score, (1.0 / peptides.size()));
    pair<FLOAT_T, string> score_pair (score, protein->first);
    proteinScores.push_back(score_pair);
  }
  // sort the pairs by their scores
  sort(proteinScores.begin(), proteinScores.end()); 
}


// gets all spit matches from a txt file created by search-for-matches
// returns TRUE on success

// TODO: make it c++
bool get_txt_matches(PeptideScore& peptideScore, ProteinPeptides& proteinPeptides, char* psm_dir) {
    DIR *dir = opendir(psm_dir);
    struct dirent* directory_entry = NULL;
    FILE* txt_file;
    
    // which index the values are in the tab delimited txt file
    int scan_index, pvalue_index, sequence_index, protein_index, field_index;

    int scan_number; 

    char line[500]; 
    char field[50]; 
    int field_length;
    
    char *ptr; // pointer to the current position in the line 
    while ((directory_entry = readdir(dir)) != NULL) {
      carp(CARP_DEBUG, "reading file: %s", directory_entry->d_name);
      if (!suffix_compare(directory_entry->d_name, "txt")) {
	continue;
      }
      txt_file = fopen(directory_entry->d_name, "r");
      if (!txt_file) {
	carp(CARP_ERROR, "Error reading %s", directory_entry->d_name);
    	return false;
      }
      
      //initialize indeces
      scan_index = -1;
      pvalue_index = -1;
      sequence_index = -1;
      protein_index = -1;
      field_index = 0;

      // store the header
      fgets(line, 500, txt_file);

      // extract the positions of the useful information from header
      ptr = line;
      while (sscanf(ptr, "%[^\t]%n", field, &field_length) == 1) {
        //carp(CARP_INFO, "%s", field);
        if (strstr("-log(p-value)", field)) {
	  pvalue_index = field_index;
        } else if (strstr("sequence", field)) {
	  sequence_index = field_index;
        } else if (strstr("protein id", field)) {
          protein_index = field_index;
        } else if (strstr("scan", field)) {
	  scan_index = field_index;
	}
        ptr += field_length; // move to next field
        if (*ptr != '\t') break; // if done
        while (*ptr == '\t') ++ptr; // adjust for multiple delimiters
        ++field_index;
      } // get next field in header 

	// die if the fields not found in the header
      if (pvalue_index   < 0
      ||  sequence_index < 0
      ||  protein_index  < 0
      ||  scan_index     < 0) {
        carp(CARP_INFO, "found invalid .txt file, skipping");
        continue;
      }
 
      carp(CARP_DETAILED_DEBUG, "header indeces: pvalue %d sequence %d protein %d",
				          pvalue_index, sequence_index, protein_index);
      scan_number = 0;    
      // get every line from txt file, adding top matches
      // to spit_matches array
      while (fgets(line, 500, txt_file)) {
        char* fields[field_index + 1]; // store fields of every line
        field_index = 0;
	ptr = line; 
        while (sscanf(ptr, "%[^\t]%n", field, &field_length)) {
	  fields[field_index++] = ptr;
	  ptr += field_length;
	  if (*ptr != '\t') break; // no more fields
	    // split the string, move to next field
          *ptr = '\0';
 	  ++ptr; 
	  // adjust for missing fields
 	  while (*ptr == '\t') {
	    ++ptr;
	    ++field_index; 
	  }
        }

	// skip if not best match for the scan
 	if (scan_number == atoi(fields[scan_index]) ) {
	  continue;
	}else {
	  scan_number = atoi(fields[scan_index]);
	}
	
	carp(CARP_DETAILED_DEBUG, "scan %d", scan_number);

	if (strcmp(fields[pvalue_index], "nan") == 0 ||
            strcmp(fields[pvalue_index], "NaN") == 0 ) {
	  carp(CARP_INFO, "skipping unscored PSM");
	  continue;
        }
	
      PeptideScore::iterator it = peptideScore.find(fields[sequence_index]);
      // if a new peptide, or a larger score, add the peptide and its score to the map
      if (it == peptideScore.end() || (it->second < atof(fields[pvalue_index]))) {
	peptideScore[fields[sequence_index]] = atof(fields[pvalue_index]);
      }	
      // split the protein string by commas
      vector<string> parent_proteins;
      string_split(string(fields[protein_index]), parent_proteins, ",");
      // for each protein
      for (vector<string>::const_iterator protein = parent_proteins.begin(); protein != parent_proteins.end(); ++protein) {
	proteinPeptides[*protein].insert(fields[sequence_index]);	
      }      
    } // get next scan
       
    fclose(txt_file);
     
  } // get another file from psm_dir
    closedir(dir);
    return true;
}

// splits a string into a vector of tokens
void string_split(const string& str,
                  vector<string>& tokens,
                  const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
