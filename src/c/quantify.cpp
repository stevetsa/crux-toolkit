/**
 * @file quantify.cpp
 *
 * This is a crux command to use after finding peptide spectrum matches with 
 * p values and returns a file with a list of proteins with scores 
 *
 *
 **/

#include "quantify.h"

using namespace std;


typedef map<string, int> MappingToSpc;                /* maps protein or peptide to 
							 thier spectrum count */
typedef map<string, int> ProteinToPepc;               /* maps protein to thier
							 peptide count */
typedef map<string, int> ProteinToLength;             /* maps protein to thier 
							 length */
typedef map<string, set<string> > ProteinToPeptides;  /* maps protein to the set
							 of peptides matched */
typedef map<string, FLOAT_T> MappingToScore;          /* maps protein or peptide
							 to thier score 
							 (nsaf or sin)*/
typedef map<string, FLOAT_T> PeptideToIntensity;      /* maps peptide thier ion 
							 intensity */
typedef vector<pair<double, string> > ScoreToMapping; /* vector of score to
							 protein or peptide*/
typedef map<pair<int, int>, string> ScanChargeToPeptides; /* maps scan/charge to
							     peptide */

// for nsaf, prints tab-delimited file of proteins and thier nsaf
// scores in sorted order from biggest to smallest
bool print_nsaf_scores(
		       ScoreToMapping& scoreToMapping,
		       string full_output_file,
		       string measure,
		       BOOLEAN_T use_prot_level
		       );

// gets the indices of all the values based off of header 
void getIndices(
		int *score_index, 
		int *proteins_index, 
		int *peptide_index, 
		int *rank_index, 
		int *scan_index, 
		int *charge_index,
		bool isPerc,
		string header
		);


// for sin, prints tab-delimited file of proteins and thier nsaf
// scores in sorted order from biggest to smallest
bool print_scores(
		  ScoreToMapping& scoreToMapping,
		  string full_output_file,
		  string measure,
		  MappingToScore mappingToScore,
		  ProteinToLength proteinToLength,
		  double score_summation,
		  map<string, int>& proteinToNumSpectra,
		  map<string, int>& proteinToPn,
		  BOOLEAN_T use_prot_level
		 );



// gets the Spc / L for each protein
bool get_nsaf_scores(
		     MappingToScore& mappingToScore,
		     MappingToSpc& mappingToSpc,
		     ProteinToLength& proteinToLength,
		     double * score_summation,
		     BOOLEAN_T use_prot_level
		     );

// gets the total ion intensity for the protein
bool get_sin_scores(
		     MappingToScore& mappingToScore,
		     ProteinToPeptides& proteinToPeptides,
		     PeptideToIntensity& peptideToIntensity,
		     double * score_summation,
		     map<string, int>& proteinToNumSpectra,
		     map<string, int>& pepToNumSpectra,
		     map<string, int>& proteinToPn
		     ); 

// get the ion intensity score by running scripts and reading
// the output
bool get_ion_intensity_scores(
			      PeptideToIntensity& peptideToIntensity, 
			      float threshold,
			      char * txt_file,
			      char * output_dir,
			      string input_ms2,
			      map<string, int>& pepToNumSpectra
			      );

// find the normalized score by dividing all the scores by the 
// summation (and by length of the protein if its for sin score)
bool get_normalized_score(
			  MappingToScore& mappingToScore,
			  ScoreToMapping& normScoreToMapping,
			  ProteinToLength& proteinToLength,
			  string measure,
			  double summation,
			  BOOLEAN_T use_prot_level
			  );

// creates and returns a match_collection from an sqt file
bool get_matches_from_txt(
			  char* txt_file,
			  char* database,
			  ProteinToLength& proteinToData,
			  MappingToSpc& mappingToSpc,
			  ProteinToPepc& proteinToPepc,
			  ProteinToPeptides& proteinToPeptides,
			  ScanChargeToPeptides& scanChargeToPeptides,
			  float threshold,
			  BOOLEAN_T only_unique,
			  BOOLEAN_T use_prot_level
			  );

/* 
 * Takes PeptideToIntensity mapping and fills them with 
 * mapping of peptides to areas of bullseye results
 * NOTE: Name of mappping is inaccurate
 */
bool get_bullseye_areas(
			PeptideToIntensity& peptideToIntensity,
			char * bullseye_file,
			ScanChargeToPeptides& scanChargeToPeptides
			);


// retreives lengths of protiens from the database
bool get_data_from_db(
		      ProteinToLength& proteinToData,
		      DATABASE_T *database
		      );

bool file_exists(string file_name);


// take in the threshold value and returns the output file name
string create_file_name(float number);

// takes in a string and returns a vector with the string split by the delimiter
void tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters);

// prints a tab delimited protein to peptide mapping with the comma delimited
// peptide list
bool printProteinToPeptides(map<string, set<string> > proteinToPeptides);

// Takes a string a returns a list of the string split by the delimiters
void split(const string& str, 
	   vector<string>& tokens, 
	   const string& delimiters);

int quantify_main(int argc, char** argv){
  // Define optional and required command line arguments
  const char* option_list[] = {
    "threshold",
    "input-ms2",
    "fileroot",
    "output-dir",
    "overwrite",
    "unique-mapping",
    "input-bullseye",
    "measure",
    "version"
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* argument_list[] = {
    "input-PSM",
    "database"
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);
  
  initialize_run(QUANTIFY_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);

  string  measure = get_string_parameter("measure");
  char* psm_dir = get_string_parameter("input-PSM");
  char* database = get_string_parameter("database");
  double threshold = get_double_parameter("threshold");
  char* output_dir = get_string_parameter("output-dir");
  char *  ms2 = get_string_parameter("input-ms2");
  char * bullseye = get_string_parameter("input-bullseye");
  string quant_level = get_string_parameter("quant-level");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  BOOLEAN_T only_unique = get_boolean_parameter("unique-mapping");

  
  // check validity of parameters
  transform(measure.begin(), measure.end(), measure.begin(), ::toupper); 
  transform(quant_level.begin(), quant_level.end(), quant_level.begin(), ::toupper); 

  if (!(measure.compare("NSAF") || measure.compare("SIN"))){
    carp(CARP_FATAL, "Must pick NSAF or SIN for measure option");
  } 
  if (!measure.compare("SIN") && ms2==NULL){
    carp(CARP_FATAL, "Must use input-ms2 option if you want to find SIN");
  }
  if (!(quant_level=="PROTEIN" || quant_level=="PEPTIDE")){
    carp(CARP_FATAL, "Must pick PROTEIN or PEPTIDE for quant-level option");
  }



  BOOLEAN_T use_prot_level = (quant_level=="PROTEIN");
  string input_ms2 = ms2;
  char* output_file = (char*)malloc(sizeof(char)*20);
  strcpy(output_file, "quantify.target.txt");
  output_file[19] = '\0';
  prefix_fileroot_to_name(&output_file);
  

  carp(CARP_INFO, "debugging");



  // creates the full path name
  string temp1 = output_dir;
  string temp2 = output_file;
  string full_output_file = temp1 +'/'+ temp2;
  
  /* Open the log file to record carp messages */
  char * log_file_name = (char*)malloc(sizeof(char)*17);
  strcpy(log_file_name, "quantify.log.txt");
  log_file_name[16] = '\0';

  open_log_file(&log_file_name);
  log_command_line(argc, argv);

  carp(CARP_INFO, "Running quantify with threshold of %f", threshold);

  create_output_directory(output_dir, overwrite);
  create_file_in_path(output_file, output_dir, overwrite);

  carp(CARP_DETAILED_DEBUG, "psm-folder: %s", psm_dir);
  carp(CARP_DETAILED_DEBUG, "database: %s", database);

  // initialize variables
  MappingToScore mappingToScore;         /* maps protein to thier score 
					   (nsaf or sin) */
  MappingToSpc mappingToSpc;             /* maps proteins to spectrum count */
  ProteinToPepc proteinToPepc;           /* maps proteins to peptide count */
  ProteinToLength proteinToLength;       /* maps proteins to length (number of 
					    amino acids). Contains all proteins
					    in database */
  ProteinToPeptides proteinToPeptides;   /* maps proteins to set of matched 
					    peptides */
  PeptideToIntensity peptideToIntensity; /* maps peptides to sum of ion 
					    intensities for all matched 
					    spectra */
  ScoreToMapping normScoreToMapping;     /* maps normalized score to protein */
  ScanChargeToPeptides scanChargeToPeptides; /* maps scan/charge to peptide */
  double score_summation;                /* sums up the scores for all proteins 
					    to use for normalization */
  map<string, int> pepToNumSpectra;      /* number of spectra with ion 
					    intensity for each peptide */
  map<string, int> proteinToNumSpectra;  /* number of spectra with ion 
					    intensity for each protein */
  map<string, int> proteinToPn;           /* number of peptides used with 
					     ion intensity for the protein */

  /* read from the database and target.txt file */
  if (!get_matches_from_txt(psm_dir, 
			    database, 
			    proteinToLength, 
			    mappingToSpc, 
			    proteinToPepc, 
			    proteinToPeptides, 
			    scanChargeToPeptides,
			    threshold, 
			    only_unique, 
			    use_prot_level)){
    carp(CARP_ERROR, "error retrieving matches");
    exit(1);
  }

  /* This part gets the (un)normalized scores depending on whether the user 
     asked for nsaf or sin score */
  if (!measure.compare("NSAF")){
    if (!get_nsaf_scores(mappingToScore, 
			 mappingToSpc, 
			 proteinToLength, 
			 &score_summation, 
			 use_prot_level
			 )){
      carp(CARP_ERROR, "error finding nsaf scores"); 
      exit(1);
    }
    carp(CARP_INFO, "%i", mappingToScore.size());
    carp(CARP_INFO, "%i", mappingToSpc.size());
  } else {

    /* if bullseye parameter was passed, then use bullseye */
    if (bullseye != NULL){
      
      /* calculate areas under peak in bullseye file */
      if (!get_bullseye_areas(
			peptideToIntensity,
			bullseye,
			scanChargeToPeptides)){
	    carp(CARP_ERROR, "error finding area under peaks in bullseye");
	    exit(1);
	  }

    } else {
      /* calculate the ion intensity scores for each protein */
      if (!get_ion_intensity_scores(peptideToIntensity, 
				    threshold ,
				    psm_dir, 
				    output_dir, 
				    input_ms2, 
				    pepToNumSpectra)){
	carp(CARP_ERROR, "error finding ion intensities");
	exit(1);
      }
    }
    
    /* using the ion intensities, find the sin scores */
    if (!use_prot_level) {
      mappingToScore = peptideToIntensity;
    } else {
      if (!get_sin_scores(mappingToScore, 
			  proteinToPeptides, 
			  peptideToIntensity, 
			  &score_summation, 
			  proteinToNumSpectra, 
			  pepToNumSpectra, 
			  proteinToPn)){
	carp(CARP_ERROR, "error finding sin scores");
	exit(1);
      }  
    }
  }
  
  /* find the normalized scores */
  if (!get_normalized_score(mappingToScore ,
			    normScoreToMapping, 
			    proteinToLength , 
			    measure ,
			    score_summation,
			    use_prot_level
			    )){
    carp(CARP_ERROR, "error finding the normalized scores");
    exit(1);
  }


  
  /* print the results */
  if (!measure.compare("NSAF")){
    if (!print_nsaf_scores(normScoreToMapping, 
			   full_output_file, 
			   measure,
			   use_prot_level
			   )){
      carp(CARP_ERROR, "error printing scores");
      exit(1);
    }
  } else{
    /* SIN has more data to print */
    if (!print_scores(normScoreToMapping, 
		      full_output_file, 
		      measure, 
		      mappingToScore, 
		      proteinToLength, 
		      score_summation, 
		      proteinToNumSpectra, 
		      proteinToPn,
		      use_prot_level
		      )){
      carp(CARP_ERROR, "error printing scores");
      exit(1);
    }
  }

  carp(CARP_INFO, "Crux quantify finished");

  return 1;  
}



// (just for nsaf) Iterate through all the proteins and output its id and nsaf scores 
// into a file. Sorted by scores and delimited by tab
bool print_nsaf_scores(
		       ScoreToMapping& scoreToMapping,
		       string full_output_file,
		       string measure,
		       BOOLEAN_T use_prot_level
		       ){
  carp(CARP_INFO, "Outputting matches");
  ofstream output_stream(full_output_file.c_str());
  string prot_level;
  if (use_prot_level){
    prot_level = "Protein id";
  } else {
    prot_level = "Peptide seq";
  }
  
  if (!output_stream.is_open()){
    carp(CARP_ERROR, "error opening output file");
    return false;
  }
  
  
  /* print header */
  output_stream<< prot_level+"\t "+measure+" score" <<endl;
  // for each protein, print out its id and score
  for (ScoreToMapping::const_iterator score_pair = scoreToMapping.begin();
	 score_pair != scoreToMapping.end(); ++score_pair){
    output_stream << score_pair->second << "\t" << score_pair->first << endl;
  }
  
  
  

  return true;
}

/* Iterate through all the proteins and output its id and nsaf scores 
   into a file. Sorted by scores and delimited by tab */
bool print_scores(
		  ScoreToMapping& scoreToMapping,
		  string full_output_file,
		  string measure,
		  MappingToScore mappingToScore,
		  ProteinToLength proteinToLength,
		  double score_summation,
		  map<string, int>& proteinToNumSpectra,
		  map<string, int>& proteinToPn,
		  BOOLEAN_T use_prot_level
		  ){
  /* create a stream to the output file */
  carp(CARP_INFO, "Outputting matches");
  ofstream output_stream(full_output_file.c_str());
  if (!output_stream.is_open()){
    carp(CARP_ERROR, "error opening output file");
    return false;
  }
  /* print header */
  if (use_prot_level){
    /* For Protein level quantification */
    output_stream<< "Protein id \t "+measure+" score" 
		 << "\ttotal_intensity" 
		 << "\tlength" 
		 << "\tnum_spectra"
		 << "\tnum_peptides"
		 << "\tsum_intensities=" 
		 << score_summation 
		 << endl;
    /* for each protein, print out its id and score */
    for (ScoreToMapping::const_iterator score_pair = scoreToMapping.begin();
	 score_pair != scoreToMapping.end(); ++score_pair){
      output_stream << score_pair->second << "\t" 
		    << score_pair->first  << "\t" 
		    << mappingToScore[score_pair->second] << "\t" 
		    << proteinToLength[score_pair->second] <<"\t" 
		    << proteinToNumSpectra[score_pair->second] << "\t" 
		    << proteinToPn[score_pair->second]
		    << endl;
    }
  } else {
    /* For Peptide level quantification */
    output_stream<< "Peptide seq \t "+measure+" score" 
		 << "\ttotal_intensity" 
		 << "\tnum_spectra"
		 << "\tsum_intensities=" 
		 << score_summation 
		 << endl;
    /* for each protein, print out its id and score */
    for (ScoreToMapping::const_iterator score_pair = scoreToMapping.begin();
	 score_pair != scoreToMapping.end(); ++score_pair){
      output_stream << score_pair->second << "\t" 
		    << score_pair->first  << "\t" 
		    << mappingToScore[score_pair->second] << "\t" 
		      << proteinToNumSpectra[score_pair->second] << "\t" 
		    << endl;
    }
  }
  
  
  return true;
}

/* Takes tne summation of spc/length and gets the nsaf scoresv for each 
 *  protein 
 */
bool get_normalized_score(
		     MappingToScore& mappingToScore,
		     ScoreToMapping& normScoreToMapping,
		     ProteinToLength& proteinToLength,
		     string measure,
		     double summation,
		     BOOLEAN_T use_prot_level
		     ){
  /* initialize variables */
  string protein_id;
  double normalized_score;
  int length;

  carp(CARP_INFO, "Finding normalized scores");
  /* for each protein */
  for (MappingToScore::const_iterator score_pair = mappingToScore.begin();
       score_pair != mappingToScore.end(); ++score_pair){
    /* get nessasary variables from maps */
    protein_id = (*score_pair).first;
    normalized_score = (*score_pair).second / summation;
    if (!measure.compare("SIN") && use_prot_level){
      length = proteinToLength[protein_id];
      normalized_score = normalized_score / (FLOAT_T)length;
    }
    /* record score into vector */
    normScoreToMapping.push_back(make_pair(normalized_score, protein_id));
  }
  /* sort the vector and reverse to order them from biggest to largest */
  sort(normScoreToMapping.begin(), normScoreToMapping.end());
  reverse(normScoreToMapping.begin(), normScoreToMapping.end());
  return true;
}



/* Gets the sin scores and enters them into mappingToScore map */
bool get_nsaf_scores(
		     MappingToScore& mappingToScore,
		     MappingToSpc& mappingToSpc,
		     ProteinToLength& proteinToLength,
		     double * score_summation,
		     BOOLEAN_T use_prot_level
		     ){
  /* initialize variables */
  string protein_id;
  int count;
  int length;
  
  carp(CARP_INFO, "Finding summation");
  /* start summation at zero */
  (*score_summation) = 0;
  /* for each protein */
  for (map<string,int>::const_iterator score_pair = mappingToSpc.begin();
       score_pair != mappingToSpc.end(); ++score_pair){
    /* get nessasary varaibles for maps */
    protein_id = (*score_pair).first;
    count = (*score_pair).second;
    length = proteinToLength[protein_id];

    /* add score to mapping and incremement normalization constant */
    if (use_prot_level){
      mappingToScore.insert(make_pair(protein_id , (float)count / (float)length));
      (*score_summation) += (float)count / (float)length;
    } else {
      mappingToScore.insert(make_pair(protein_id , (float)count));
      (*score_summation) += (float)count;
    }

  }
  return true;
}

/* Gets the sin scores and enters them into mappingToScore map. Assumes 
   that peptideToIntensity is filled out */
bool get_sin_scores(
		     MappingToScore& mappingToScore,
		     ProteinToPeptides& proteinToPeptides,
		     PeptideToIntensity& peptideToIntensity,
		     double * score_summation,
		     map<string, int>& proteinToNumSpectra,
		     map<string, int>& pepToNumSpectra,
		     map<string, int>& proteinToPn
		    ){
  string protein_id;
  string peptide_id;
  FLOAT_T total_intensity;
  int total_spectra;
  int total_peptide;
  (*score_summation) = 0;

  /* iterate through all peptides of all proteins while adding up the intensity 
     of each spectra and storing them into peptideToIntensity */
  for (ProteinToPeptides::const_iterator proteinIterator = proteinToPeptides.begin();
       proteinIterator != proteinToPeptides.end(); ++proteinIterator){
    protein_id = (*proteinIterator).first;
    //carp(CARP_INFO, "For Portein %s", protein_id.c_str());
    total_intensity = 0.0;
    total_spectra = 0;
    total_peptide = 0;
    /* for each peptide */
    for (set<string>::const_iterator peptideIterator = 
	   proteinToPeptides[protein_id].begin(); peptideIterator 
	   != proteinToPeptides[protein_id].end(); ++peptideIterator){
      peptide_id = (*peptideIterator);
      /* add to sum if the intensity exists */
      if (peptideToIntensity.find(peptide_id) != peptideToIntensity.end()){
	total_intensity += peptideToIntensity[peptide_id];
	total_spectra += pepToNumSpectra[peptide_id];
	total_peptide += 1;
	//carp (CARP_INFO, "successfully added ion score");
      } else {
	//carp (CARP_INFO, "\tMissing ion intensity score for peptide %s ", 
	//peptide_id.c_str());
      }
    }
    /* extra test since bullseye scores do not match up perfectly */
    if (total_intensity != 0){
      mappingToScore.insert(make_pair(protein_id, total_intensity));
      proteinToNumSpectra.insert(make_pair(protein_id, total_spectra));
      proteinToPn.insert(make_pair(protein_id, total_peptide));
      (*score_summation) += total_intensity;
    }
  }
  return true;
}



// converts a float to a string and returns it
string create_file_name(float number){
   stringstream ss;//create a stringstream
   ss << "nsaf_"; 
   ss << number ;//add number to the stream
   ss << ".target.txt";
   return ss.str();//return a string with the contents of the stream
}



/* Iterates through database and the text file to retrieve
 * protein lengths and spc scores 
 */
bool  get_matches_from_txt(
			   char* txt_file,
			   char* database_file,
			   ProteinToLength& proteinToData,
			   MappingToSpc& mappingToSpc,
			   ProteinToPepc& proteinToPepc,
			   ProteinToPeptides& proteinToPeptides,
			   ScanChargeToPeptides& scanChargeToPeptides,
			   float threshold,
			   BOOLEAN_T only_unique,
			   BOOLEAN_T use_prot_level
			   ){
  DATABASE_T* database = NULL;
  BOOLEAN_T use_index = is_directory(database_file);
  char* binary_fasta = NULL;
  /* TODO check the validity of files */
 
  /* create new database object */
  if (use_index){ 
    binary_fasta = get_index_binary_fasta_name(database_file);
  } else {
    binary_fasta = get_binary_fasta_name(database_file);
    carp(CARP_DEBUG, "Looking for binary fasta %s", binary_fasta);
    if (access(binary_fasta, F_OK)){
      carp(CARP_DEBUG, "Could not find binary fasta %s", binary_fasta);
      if (!create_binary_fasta_here(binary_fasta, binary_fasta)){
	carp(CARP_FATAL, "Could not create binary fasta file %s", binary_fasta);
      };
    }
  }
  database = new_database(binary_fasta, TRUE);
  /* check if already parsed */
  if(!get_database_is_parsed(database)){
    carp(CARP_DETAILED_DEBUG,"Parsing database");
    if(!parse_database(database)){
      carp(CARP_FATAL, "Failed to parse database, cannot create new index");
    }
  }
  free(binary_fasta);


  /* find all lengths from database */
  get_data_from_db(proteinToData, database);
  free_database(database);
  


  /* find indices of which element are the scores for a particular file (temporary solution) */
  int score_index, proteins_index, 
    peptide_index, rank_index,
    scan_index, charge_index;
  string file_name = txt_file;

  carp(CARP_INFO, "Getting match scores from text file");
  ifstream myfile (txt_file);
  string line;
  getline(myfile, line);
  /* get indices from header file */
  getIndices(&score_index, 
	     &proteins_index, 
	     &peptide_index, 
	     &rank_index, 
	     &scan_index, 
	     &charge_index,
	     file_name.find("percolator")!=string::npos,
	     line);
  string peptide_sequence;
  string protein_id;
  float score = 0.0;
  int rank, scan, charge;
  vector<string> match_fields;
  vector<string> matched_proteins;
  string PSMId;

  while (getline(myfile, line) && score < threshold){
    split(line, match_fields, "\t");
    rank = atoi(match_fields.at(rank_index).c_str());
    score = atof(match_fields.at(score_index).c_str());
    carp(CARP_INFO, "rank %i score %f ", rank, score);
    /* only record if score is below threshold and rank is one */
    if (score < threshold && rank == 1){
      scan = atoi(match_fields.at(scan_index).c_str());
      charge = atoi(match_fields.at(charge_index).c_str());
      peptide_sequence = match_fields.at(peptide_index);
      scanChargeToPeptides.insert(make_pair(make_pair(scan ,charge ),
					    peptide_sequence));
      tokenize(match_fields.at(proteins_index), matched_proteins, ",");
      carp(CARP_INFO, "scan %i charge %i %s", scan, charge, peptide_sequence.c_str());
      
      /* for each match: */
      if (use_prot_level){
	/* increment spectra/peptide counts for protein level */
	if (matched_proteins.size() == 1 || !only_unique){
	  for (vector<string>::iterator it = matched_proteins.begin(); 
	       it != matched_proteins.end(); ++it){
	    protein_id = *it;
	    carp(CARP_INFO, "protein %s", protein_id.c_str());
	    /* check if peptide has already been counted */
	    if (proteinToPeptides[protein_id].find(peptide_sequence) == 
		proteinToPeptides[protein_id].end()){
	      proteinToPeptides[protein_id].insert(peptide_sequence);
	      if (proteinToPepc.find(protein_id) == proteinToPepc.end()){
		proteinToPepc.insert(make_pair(protein_id, 0));
		mappingToSpc.insert(make_pair(protein_id, 0));
	      }
	      proteinToPepc[protein_id]++;
	    }
	    mappingToSpc[protein_id]++;
	  }
	}
      }else {
	/*incrememt spectra counts for peptide level*/
	if (mappingToSpc.find(peptide_sequence) == mappingToSpc.end()){
	  mappingToSpc.insert(make_pair(peptide_sequence, 0));
	}
	mappingToSpc[peptide_sequence]++;
      }
      
    }
    matched_proteins.clear();
    match_fields.clear();
  }
  myfile.close();
  /*
  for (MappingToSpc::const_iterator it = mappingToSpc.begin();
       it != mappingToSpc.end(); ++it){
    protein_id = (*it).first;
    carp(CARP_INFO, "%s  %d", protein_id.c_str(), (*it).second);
  }
  */

  return true;
}



/* Gets the ion_intensity_scores for each peptide by first running
 * couple scripts and then reading from the last output file
 */
bool get_ion_intensity_scores(
			      PeptideToIntensity& peptideToIntensity, 
			      float threshold,
			      char * txt_file,
			      char * output_path,
			      string input_ms2,
			      map<string, int>& pepToNumSpectra
			      ){
  carp(CARP_INFO, "Getting ion intensity scores");
  /* TODO: Somehow not rely on knowing the paths...*/
  string output_dir = output_path;
  output_dir += "/";
  string cluster_python = "/net/gs/vol3/software/Python-2.5.2/bin/python";
  string output_mgf = output_dir+"output.mgf";
  string output_seq = output_dir+"seq_output.mgf";
  string output_annotate = output_dir+"annotated-spectrum.mgf2";
  string path_to_mgf_script = "/net/noble/vol2/home/mkm24/proj/rubel/results_Michael"
    "/2010-01-04/crux-quantify/src/perl/ms22mgf.pl";
  string path_to_seq_script = "/net/noble/vol2/home/mkm24/proj/rubel/results_Michael/"
    "2010-01-04/crux-quantify/src/python/seq_into_mgf_exclu.py";
  string path_to_annotate_script = "/net/noble/vol2/home/mkm24/proj/rubel/results_Michael/"
    "2010-01-04/crux-quantify/src/python/psm-annotator.py";
  
  stringstream ss (stringstream::in | stringstream::out);
  ss << threshold;
  string threshold_s = ss.str();
  cout << threshold_s <<endl;

  /* make sure python files exists */
  if (!(file_exists(path_to_mgf_script) 
	&& file_exists(path_to_seq_script)
	&& file_exists(path_to_annotate_script))){
    carp(CARP_FATAL, "Missing script(s) to find ion intensity scores.");
  }
  
  string  get_mgf_cmd =  path_to_mgf_script+" "+input_ms2+" "+output_mgf;
  string get_seq_cmd = cluster_python+" "+path_to_seq_script+" "
    +output_mgf+" "+txt_file+ " "+output_dir+ " "+threshold_s;
  string annotate_cmd = cluster_python+" "+path_to_annotate_script+" "
    +output_seq+" "+output_dir;
  

  /* run all nessasary scripts */
  carp(CARP_INFO, "Running script: %s" , get_mgf_cmd.c_str());
  system(get_mgf_cmd.c_str());
  carp(CARP_INFO, "Running script: %s" , get_seq_cmd.c_str());
  system(get_seq_cmd.c_str());
  carp(CARP_INFO, "Running script: %s", annotate_cmd.c_str());
  system(annotate_cmd.c_str());

  /* iterate through the file */
  ifstream mgf2_file ((char *)output_annotate.c_str());
  FLOAT_T totalIntensity;
  string line;
  string ion;
  string sequence_line;
  string sequence;
  vector<string> all_sequences;
  vector<string> match_fields;
  while (1){
    /* keep reading lines while its not begin ions */
    while(getline(mgf2_file, line) && line.compare("BEGIN IONS")){}
    if (line == ""){ 
      /* stop if there are no more lines */
      break;
    }
    totalIntensity = 0.0;
    getline(mgf2_file, line);
    getline(mgf2_file, line);
    sequence_line = line.substr(4,string::npos);
    /* keep reading lines while its not end ions */
    while(line.find("END IONS") == string::npos){
      getline(mgf2_file, line);
      tokenize(line, match_fields, " "); 
      if ( (int)match_fields.size() > 2){  // must have Y or B
	ion = match_fields[2];
	/* must match [b|y]  */
	if ((ion.find("b")!=string::npos || ion.find("y")!=string::npos) 
	    && ion.find("H2O") == string::npos){
	  //carp(CARP_INFO, "\tadding ion intensity: %s", (char *)match_fields[1].c_str());
	  totalIntensity += (FLOAT_T) atof((char *)match_fields[1].c_str());
	}
      }
      /* reset the vector */
      match_fields.clear();
    }

    tokenize(sequence_line, all_sequences, ",");
    /* record the total ion intensity for all peptides */
    for (vector<string>::iterator it = all_sequences.begin(); it != all_sequences.end(); ++it){
      sequence = *it;
      if (peptideToIntensity.find(sequence) == peptideToIntensity.end()){
	/* if doesn't exist, insert a new entry */
	peptideToIntensity.insert(make_pair(sequence, totalIntensity));
	pepToNumSpectra.insert(make_pair(sequence, 1));
	carp(CARP_INFO, "peptide: %s\t addingScore: %f", sequence.c_str(), totalIntensity);
      }else{
	/* if peptide exists, then add to the totalIntensity */
	peptideToIntensity[sequence]+= totalIntensity;
	pepToNumSpectra[sequence]+=1;
      }
    }
    /* reset vector */
    all_sequences.clear();
  }
  
  return true;
}


bool printProteinToPeptides(DATABASE_T* database){
    ofstream output_stream("pep_to_pro.txt");
    PEPTIDE_T* peptide;
    PEPTIDE_CONSTRAINT_T* peptide_constraint = new_peptide_constraint_from_parameters();
    string proteins_list;
    string peptide_str;
    DATABASE_PEPTIDE_ITERATOR_T* it = new_database_peptide_iterator(database, peptide_constraint);
    while (database_peptide_iterator_has_next(it)){
      peptide = database_peptide_iterator_next(it);
      proteins_list = get_protein_ids(peptide);
      peptide_str = get_peptide_sequence(peptide);
      output_stream << peptide_str << "\t" << proteins_list << endl;
    }
    output_stream.close();
    return true;
}



/* Takes in the database and a map as input and 
 * iterates through the database and records all lengths
 * into the map 
 */
bool get_data_from_db(
		      ProteinToLength& proteinToData,
		      DATABASE_T*  database
		      ){
  DATABASE_PROTEIN_ITERATOR_T* it =  new_database_protein_iterator(database);
  int find_protein_mapping = false;

  PROTEIN_T* protein;
  carp(CARP_INFO, "Getting protein lengths from database");
  while (database_protein_iterator_has_next(it)){
    protein = database_protein_iterator_next(it);

    proteinToData.insert(make_pair(string(get_protein_id_pointer(protein)), 
				   get_protein_length(protein)));
  }

  /* printing out protein to peptide mapping */
  if (find_protein_mapping){
    printProteinToPeptides(database);
  }
  return true;
}



/* 
 * Takes PeptideToIntensity mapping and fills them with 
 * mapping of peptides to areas of bullseye results
 * NOTE: Name of mappping is inaccurate
 */
bool get_bullseye_areas(
			PeptideToIntensity& peptideToIntensity,
			char * bullseye_file,
			ScanChargeToPeptides& scanChargeToPeptides
			){
  carp(CARP_INFO, "Finding area under peaks from bullseye");
  
  if (!file_exists(bullseye_file)){
    carp(CARP_FATAL, "File doesn't exist, exiting quantify");
    exit(1);
  }
  
  ifstream myfile (bullseye_file);
  string line, prev;
  vector<string> split_line;

  /* get rid of Header lines */
  getline(myfile, line);



  while (line != "" && line.at(0) == 'H'){
    prev = line;
    getline(myfile, line);
  }

  /* Bullseye accidently outputs first scan on the last Header line */
  float best_area, pre_mz, best_mz, mz;
  int scan, charge, temp_charge;
  ScanChargeToPeptides::iterator it;
  tokenize(prev, split_line, "\t");
  pre_mz = atof(split_line.at(5).c_str());
  scan = atoi(split_line.at(3).c_str());
  best_area = 0.0;
  best_mz = 0.0;
  charge = 0;
  while(1){
    /* finding the area for the mz closest to pre_mz */
    while (line != "" && line.at(0) != 'S'){
      if (line.find("EZ") != string::npos){
	split_line.clear();
	tokenize(line, split_line ,"\t");
	temp_charge = atoi(split_line.at(2).c_str());
	mz = atof(split_line.at(3).c_str()) / (float) temp_charge;
	if (fabs(pre_mz-mz) < fabs(best_mz-pre_mz)){
	  best_mz = mz;
	  best_area = atoi(split_line.at(5).c_str());
	  charge = temp_charge;
	}
      }
      getline(myfile, line);
    }
    /* record best results into the mapping */
    it = scanChargeToPeptides.find(make_pair(scan, charge));
    if (it != scanChargeToPeptides.end()){
      //carp(CARP_INFO, "%f", best_area);
      peptideToIntensity.insert(make_pair((*it).second, best_area));
    }
    /* end of the file */
    if (line == ""){
      break;
    }
    /* next scan*/
    split_line.clear();
    tokenize(line, split_line, "\t");
    scan = atoi(split_line.at(1).c_str());
    pre_mz = atof(split_line.at(3).c_str());
    getline(myfile, line);
  }
    
  myfile.close();
  return true;
}


/* Takes a string, vector and delimiters and splits the string by delimiters and
 * puts them into the vector 
 */
void tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters)
{
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos){
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(" ", 
				    str.find_first_not_of(delimiters, pos));
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void getIndices(int *score_index, 
		int *proteins_index, 
		int *peptide_index, 
		int *rank_index, 
		int *scan_index, 
		int *charge_index,
		bool isPerc,
		string header){
  int i=0;
  vector<string> elements;
  split(header, elements, "\t");
  
  for (vector<string>::iterator iter = elements.begin(); 
       iter != elements.end(); ++iter){
    if ("scan"==*iter){
      (*scan_index)=i;
    } else if ("charge"==*iter){
      (*charge_index)=i;
    } else if ("percolator q-value"==*iter &&
	       isPerc){
      (*score_index)=i;
    } else if ("q-ranker q-value"==*iter &&
	       !isPerc){
      (*score_index)=i;
    } else if ("xcorr rank"==*iter){
      (*rank_index)= i;
    } else if ("protein id"==*iter){
      (*proteins_index) = i;
    } else if ("sequence"==*iter){
      (*peptide_index) = i;
    }
    i++;
  }
}




void split(const string& str, vector<string>& tokens, const string& delimiters)
{
  string::size_type lastPos = (size_t) 0;
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  while (string::npos != pos && string::npos != lastPos){
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = ++pos;
    pos = str.find_first_of(delimiters, lastPos);
  }
}

/* returns true if file exists and false otherwise */
bool file_exists(string filePath){
  ifstream file(filePath.c_str());
  bool result = file.is_open();
  file.close();
  return result;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
