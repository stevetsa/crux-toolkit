#include "spectral-counts.h"

using namespace std;

/**
 * \typedef peptideToScore
 * \brief Mapping of peptide object to scores
 */
typedef map<PEPTIDE_T*, FLOAT_T, bool(*)(PEPTIDE_T*, PEPTIDE_T*)> PeptideToScore;
/**
 * \typedef ProteinToScore
 * \brief Mapping of protein object to scores
 */
typedef map<PROTEIN_T*, FLOAT_T, bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToScore;
/**
 * \typedef MetaProtein
 * \brief Collection of protein objects
 */
typedef set<PROTEIN_T*, bool(*)(PROTEIN_T*, PROTEIN_T*)> MetaProtein;
/**
 * \typedef PeptideSet
 * \brief Collection of peptide objects (not a meta-peptide)
 */
typedef set<PEPTIDE_T*, bool(*)(PEPTIDE_T*, PEPTIDE_T*)> PeptideSet;
/**
 * \typedef MetaMapping
 * \brief Mapping of peptideSet to MetaProtein
 * Each entry is a set of peptides mapped to a set of proteins of which
 * all contain the set of peptides
 */
typedef map<PeptideSet, MetaProtein, bool(*)(PeptideSet, PeptideSet) > MetaMapping;
/**
 * \typedef ProteinToPeptides
 * Mapping of Protein objects to a set of peptides that are part
 * of the protein sequence
 */
typedef map<PROTEIN_T*, PeptideSet , bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToPeptides;
/**
 * \typedef MetaToRan
 * \brief Mapping of MetaProtein to ranks to the rank asigned to it
 */
typedef map<MetaProtein, int, bool(*)(MetaProtein, MetaProtein)> MetaToRank;
/**
 * \typedef MetaToScore
 * \brief Mapping of MetaProtein to the score assigned to it
 */
typedef map<MetaProtein, FLOAT_T, bool(*)(MetaProtein, MetaProtein)> MetaToScore;
/**
 * \typedef ProteinToMeta
 * \brief Mapping of Protein to MetaProtein to which it belongs
 */
typedef map<PROTEIN_T*, MetaProtein, bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToMetaProtein;

/* private function declarations */
string pepsToString(PeptideSet s);
void getDirPath(char* path, char** dir);
set<MATCH_T*> filterMatches(
  MATCH_COLLECTION_ITERATOR_T* match_collection_it
);
void get_peptide_scores(
  set<MATCH_T*>*  matches, 
  PeptideToScore* peptideToScore
);
void get_protein_scores(
  PeptideToScore* peptideToScore, 
  ProteinToScore* proteinToScore
);
void get_protein_to_peptides(
  PeptideToScore* peptideToScore, 
  ProteinToPeptides* proteinToPeptides
);
void get_protein_to_meta_protein(
  MetaMapping* metaMapping, 
  ProteinToMetaProtein* proteinToMetaProtein
);
void get_meta_mapping(
  ProteinToPeptides* proteinToPeptides, 
  MetaMapping* metaMapping
);
void get_meta_ranks(
  MetaToScore* metaToScore, 
  MetaToRank* metaToRank
);
void get_meta_scores(
  MetaMapping* metaMapping, 
  ProteinToScore* proteinToScore, 
  MetaToScore* metaToScore
);
void perform_parsimony_analysis(
  MetaMapping* metaMapping
);
void normalize_peptide_scores(
  PeptideToScore* peptideToScore
);
void normalize_protein_scores(
  ProteinToScore* proteinToScore
);
void print_protein_results(
  ProteinToScore* proteinToScore,
  MetaToRank* metaToRank,
  ProteinToMetaProtein* proteinToMeta,
  char * pathToOutput
);
void print_peptide_results(
  PeptideToScore* peptideToScore,
  char * output_path
);
void make_unique_mapping(
  PeptideToScore* peptideToScore
);
string meta_protein_to_string(MetaProtein s);
string peps_to_string(PeptideSet s);
/* comparison function declarations */
bool compare_pep(PEPTIDE_T*, PEPTIDE_T*);
bool compare_prot(PROTEIN_T*, PROTEIN_T*);
bool compare_peptide_sets(PeptideSet, PeptideSet);
bool compare_meta_proteins(MetaProtein, MetaProtein);
bool compare_sets(
  pair<PeptideSet, MetaProtein>, 
  pair<PeptideSet, MetaProtein> 
);




int spectral_counts_main(int argc, char** argv){
  const char* option_list[] = {
    "parsimony",
    "threshold",
    "input-ms2",
    "fileroot",
    "output-dir",
    "overwrite",
    "unique-mapping",
    "input-bullseye",
    "quant-level",
    "measure",
    "version",
    "verbosity"
  };
  const char* argument_list[] = {
    "input PSM",
    "protein database"
  };

  int num_options = sizeof(option_list) / sizeof(char*);
  int num_arguments = sizeof(argument_list) / sizeof(char*);


  initialize_run(SPECTRAL_COUNTS_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);

  bool unique_mapping = get_boolean_parameter("unique-mapping");
  char * psm_file = get_string_parameter("input-PSM");
  char * database = get_string_parameter("database");
  char * parsimony = get_string_parameter("parsimony");
  char * quant_level = get_string_parameter("quant-level");
  char * output_dir = get_string_parameter("output-dir");
  char * measure = get_string_parameter("measure"); 
  char * ms2 = get_string_parameter("input-ms2");
  char * dir_path = NULL;

  /* Check input validity */
  if (strcmp(measure, "SIN") != 0 &&
      strcmp(measure, "NSAF") != 0 ){
    carp(CARP_FATAL, "SIN or NSAF are the only "
	 "valid input for --measure");
  }
  if (strcmp(measure, "SIN") == 0 &&
      ms2 == NULL){
    carp(CARP_FATAL, "Must specify ms2 file "
	 "with --input-ms2 option if SIN is "
	 "to be measured.");
  }
  if ( strcmp(quant_level,"PROTEIN") !=0 && 
       strcmp(quant_level,"PEPTIDE") != 0 ){
    carp(CARP_FATAL, "Must pick PROTEIN or PEPTIDE for "
	 "quant-level option");
  } 
  if (strcmp(parsimony, "NONE") != 0 && 
      strcmp(parsimony, "SIMPLE") != 0 &&
      strcmp(parsimony, "GREEDY") != 0 ){
    carp(CARP_FATAL, "Parsimony option, "
	 "must be set to NONE|SIMPLE|GREEDY");
  }
  if (!is_directory(output_dir)){
    carp(CARP_FATAL, 
	 "Value passed for output-dir is not a directory");
  }

  char * output_path = 
    get_full_filename(output_dir, "spectral-count.target.txt");
  int decoy_count = 0;
  getDirPath(psm_file, &dir_path);
  // create match collection
  MATCH_COLLECTION_ITERATOR_T* match_collection_it 
    = new_match_collection_iterator(dir_path, database, &decoy_count);
   
  // get a set of matches that pass threshold
  set<MATCH_T*> matches = filterMatches(match_collection_it);
  carp(CARP_INFO, "Number of matches passed the threshold %i", 
       matches.size());

  
  PeptideToScore peptideToScore(compare_pep);
  get_peptide_scores(&matches, &peptideToScore);
  if (unique_mapping){
    make_unique_mapping(&peptideToScore);
  }
  carp(CARP_INFO, "Number of Unique Peptides %i", peptideToScore.size());

  if (strcmp(quant_level, "PEPTIDE") == 0){ // peptide level
    normalize_peptide_scores(&peptideToScore);
    print_peptide_results(&peptideToScore, output_path);
  } else { // protein level
    ProteinToScore proteinToScore(compare_prot);
    ProteinToPeptides proteinToPeptides(compare_prot);
    ProteinToMetaProtein proteinToMeta(compare_prot);
    MetaMapping metaMapping(compare_peptide_sets);
    MetaToScore metaToScore(compare_meta_proteins);
    MetaToRank metaToRank(compare_meta_proteins);
    
    get_protein_scores(&peptideToScore, &proteinToScore);
    normalize_protein_scores(&proteinToScore);
    carp(CARP_INFO, "Number of Proteins %i", proteinToScore.size());
    
    
    if (strcmp(parsimony, "NONE") != 0){ //if parsimony is not none
      get_protein_to_peptides(&peptideToScore, &proteinToPeptides);
      get_meta_mapping(&proteinToPeptides, &metaMapping);
      get_protein_to_meta_protein(&metaMapping, &proteinToMeta);
      carp(CARP_INFO, "Number of meta proteins %i", metaMapping.size());
      
      if (strcmp(parsimony, "GREEDY") == 0){ //if parsimony is greedy
	perform_parsimony_analysis(&metaMapping);
      }
      get_meta_scores(&metaMapping, &proteinToScore, &metaToScore);
      get_meta_ranks(&metaToScore, &metaToRank);
    }
    
    print_protein_results(&proteinToScore,
			&metaToRank,
			&proteinToMeta,
			output_path);
  }
  
  free(dir_path);
  return 1;
}

/**
 * For every protein that can be mapped from the set of 
 * peptides in PeptideToScore map, enter the protein and
 * the set of identified peptides it maps to, into 
 * ProteinToPeptide
 */
void get_protein_to_peptides(
  PeptideToScore* peptideToScore, 
  ProteinToPeptides* proteinToPeptides
){
  for (PeptideToScore::iterator pep_it = peptideToScore->begin();
       pep_it != peptideToScore->end(); ++pep_it){
    PEPTIDE_T* peptide = pep_it->first;
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator =
      new_peptide_src_iterator(peptide);
    while( peptide_src_iterator_has_next(peptide_src_iterator)) {
      PEPTIDE_SRC_T* peptide_src = peptide_src_iterator_next(peptide_src_iterator);
      PROTEIN_T* protein = get_peptide_src_parent_protein(peptide_src);
      if (proteinToPeptides->find(protein) == proteinToPeptides->end()){
	PeptideSet newset(compare_pep);
	proteinToPeptides->insert(make_pair(protein, newset));
      }
      (*proteinToPeptides)[protein].insert(peptide);
    }
    free(peptide_src_iterator);   
  }
}


/**
 * Enters the mapping of protein to its metaProtein
 * into ProteinToMetaProtein. MetaProteins are retreieved
 * from MetaMapping
 *
 */
void get_protein_to_meta_protein(
			     MetaMapping* metaMap, 
			     ProteinToMetaProtein* proteinToMeta
			     ){
  // for every meta protein
  for (MetaMapping::iterator meta_protein_it = metaMap->begin();
       meta_protein_it != metaMap->end(); ++meta_protein_it){
    MetaProtein proteins = meta_protein_it->second;
    // for every protein in the meta protein
    for (MetaProtein::iterator proteins_it = proteins.begin();
	 proteins_it != proteins.end(); ++proteins_it){
      // create a mapping of protein to meta protein
      proteinToMeta->insert(make_pair((*proteins_it), proteins));
    }
  }
}


/**
 * A score for each protein is calculated by summing 
 * the scores of each peptide that belongs to a protein
 *
 */
void get_protein_scores(
		      PeptideToScore* peptideToScore,
		      ProteinToScore* proteinToScore
		      ){
  
  // iterate through each peptide
  for (PeptideToScore::iterator pep_it = peptideToScore->begin();
       pep_it != peptideToScore->end(); ++pep_it){
    PEPTIDE_T* peptide = pep_it->first;
    FLOAT_T pep_score = pep_it->second;
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator =
      new_peptide_src_iterator(peptide);
    while( peptide_src_iterator_has_next(peptide_src_iterator)) {
      PEPTIDE_SRC_T* peptide_src = peptide_src_iterator_next(peptide_src_iterator);
      PROTEIN_T* protein = get_peptide_src_parent_protein(peptide_src);
      if (proteinToScore->find(protein) == proteinToScore->end()){
	proteinToScore->insert(make_pair(protein, 0.0));
      }
      (*proteinToScore)[protein] += pep_score;
    }
    free(peptide_src_iterator);
  }
}


/**
 * Takes the PeptideToScores map and updates all
 * values with normalized values. Normalized by sum
 * of all scores and then by the peptide length
 *
 */
void normalize_peptide_scores(PeptideToScore* peptideToScore){
  carp(CARP_INFO, "Normalizing peptide scores");
  FLOAT_T total = 0.0;

  // calculate sum of all scores
  for (PeptideToScore::iterator it = peptideToScore->begin();
       it != peptideToScore->end(); ++it){
    FLOAT_T score = it->second;
    total += score;
  }
  
  // normalize by sum of scores and length
  for (PeptideToScore::iterator it = peptideToScore->begin();
       it != peptideToScore->end(); ++it){
    FLOAT_T score = it->second;
    PEPTIDE_T* peptide = it->first;
    it->second = score / total / strlen(get_peptide_sequence(peptide));
  }
  
}


/**
 * Takes ProteinToScore mapping and updates all the scores
 * with normalized values. Normalized by sum of all scores,
 * then by the protein length
 *
 */

void normalize_protein_scores(ProteinToScore* proteinToScore){
  carp(CARP_INFO, "Normalizing protein scores");
  FLOAT_T total = 0.0;

  // calculate sum of all scores
  for (ProteinToScore::iterator it = proteinToScore->begin();
       it != proteinToScore->end(); ++it){
    FLOAT_T score = it->second;
    total += score;
  }

  // normalize by sum of all scores and length
  for (ProteinToScore::iterator it = proteinToScore->begin();
       it != proteinToScore->end(); ++it){
    FLOAT_T score = it->second;
    PROTEIN_T* protein = it->first;
    it->second = score / total / strlen(get_protein_sequence(protein));
  }

}


/**
 * For SIN, it will parse the spectras from ms2 to
 * find the total ion intensity from b and y ions without
 * h20 modifications. These ion intensities are summed
 * for each peptide 
 *
 * For NSAF, it will sum up the counts of matches for each
 * peptide
 *
 */
void get_peptide_scores(
		      set<MATCH_T*> * matches,
		      PeptideToScore* peptideToScore
		      ){
  
  char * measure = get_string_parameter("measure"); //TODO check this value
  char * ms2 = get_string_parameter("input-ms2");
  FLOAT_T bin_width = get_double_parameter("mz-bin-width");
  bool isSin = (measure == NULL || !strcmp(measure, "SIN"));
  map<pair<int,int>, Spectrum*> spectras;
  
  SCORER_TYPE_T score_type = XCORR;
    
  // for SIN, parse out spectrum collection from ms2 fiel
  if (isSin){ 
    IonSeries ion_series;
    SpectrumCollection* spectrumCollection = new SpectrumCollection(ms2);
    if (!spectrumCollection->parse()){
      carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2);
    } 
    // get spectras from ms2 file
    FilteredSpectrumChargeIterator* spectrum_iterator = 
      new FilteredSpectrumChargeIterator(spectrumCollection);
    while (spectrum_iterator->hasNext()){
      int charge = 0;
      Spectrum* spectrum = spectrum_iterator->next(&charge);
      spectras.insert(make_pair(make_pair(spectrum->getFirstScan(), 
					  charge), spectrum));
    }
    carp(CARP_INFO, "Number of Spectras %i", spectras.size());
  }
  

  for(set<MATCH_T*>::iterator match_it = matches->begin();
      match_it != matches->end(); ++match_it){
    int match_intensity = 0;
    MATCH_T* match = (*match_it);
    // for sin, calculate total ion intensity for match by
    // summing up peak intensities
    if (isSin){ 
      char* peptide_seq = get_match_sequence(match);
      MODIFIED_AA_T* modified_sequence = get_match_mod_sequence(match);
      int charge = get_match_charge(match);
      Spectrum* temp = get_match_spectrum(match);
      int scan = temp -> getFirstScan();
      Spectrum* spectrum = spectras[make_pair(scan, charge)];
      Ion* ion;
      IonConstraint* ion_constraint = 
	IonConstraint::newIonConstraintSmart(score_type, charge);
      IonSeries* ion_series = new IonSeries(ion_constraint, charge);
      ion_series->update(peptide_seq, modified_sequence);
      ion_series->predictIons();
      for (IonIterator ion_it = ion_series->begin(); 
	   ion_it != ion_series->end(); ++ion_it){
	ion = (*ion_it);
	if (ion -> getType() == B_ION || ion -> getType() == Y_ION){
	  if (!ion->isModified()){
	    PEAK_T* peak = spectrum->getNearestPeak(ion->getMassZ(),
						    bin_width);
	    if (peak != NULL){
	      match_intensity += get_peak_intensity(peak);
	    }
	  }
	}
      }
    } else { // for nsaf every match is counted as 1 
      match_intensity = 1;
    }
    // add ion_intensity to peptide scores
    PEPTIDE_T* peptide = get_match_peptide(match);
    if (peptideToScore->find(peptide) ==  peptideToScore->end()){
      peptideToScore->insert(make_pair(peptide, 0.0));
    } 
    (*peptideToScore)[peptide]+=match_intensity;
  }
}



/**
 * For every match in the match collection, return a set of 
 * matches that have a qvalue score lower than user-specified
 * threshold and has a XCORR ranking of 1
 *
 */
set<MATCH_T*> filterMatches(MATCH_COLLECTION_ITERATOR_T* match_collection_it){
  set<MATCH_T*> match_set;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  FLOAT_T threshold = get_double_parameter("threshold");
  bool qualify = FALSE;
  while (match_collection_iterator_has_next(match_collection_it)){

    
    match_collection = match_collection_iterator_next(match_collection_it);
    match_iterator = new_match_iterator(match_collection, XCORR, TRUE);

    while(match_iterator_has_next(match_iterator)){
      MATCH_T* match = match_iterator_next(match_iterator);
      qualify = FALSE;
      if (get_match_rank(match, XCORR) != 1){
	continue;
      }
      // find a qvalue score lower than threshold
      if (get_match_score(match, PERCOLATOR_QVALUE) != FLT_MIN &&
	  get_match_score(match, PERCOLATOR_QVALUE) <= threshold)  {
	qualify = TRUE;
      } else if (get_match_score(match, QRANKER_QVALUE) != FLT_MIN &&
	       get_match_score(match, QRANKER_QVALUE) <= threshold)  {
	qualify = TRUE;
      } else if (get_match_score(match, DECOY_XCORR_QVALUE) != FLT_MIN &&
		 get_match_score(match, DECOY_XCORR_QVALUE) <= threshold)  {
	qualify = TRUE;
      } 
      
      if (qualify == TRUE){
	  match_set.insert(match);
	}
      }
    }
  return match_set;
}




/**
 * Fills in the MetaMapping with entries of set of 
 * peptides that can be found in every protein in
 * the meta protein
 *
 */
void get_meta_mapping(
		    ProteinToPeptides* proteinToPeptides,
		    MetaMapping* metaMapping
		    ){
  carp(CARP_INFO, "Creating a mapping of meta protein to peptides");
  int count = 0;
  for (ProteinToPeptides::iterator prot_it = proteinToPeptides->begin();
       prot_it != proteinToPeptides->end(); ++prot_it){
    PROTEIN_T* protein = prot_it->first;
    PeptideSet pep_set = prot_it->second;

    if (metaMapping->find(pep_set) == metaMapping->end()){
      MetaProtein meta_protein(compare_prot);
      count++;
      metaMapping->insert(make_pair(pep_set, meta_protein));
    }
    (*metaMapping)[pep_set].insert(protein);
  }

}


/**
 * Takes a mapping of set of peptides to meta proteins and 
 * a mapping of protein to scores, and finds the largest score
 * of the meta proteins for each protein. The object returned
 * is a mapping of MetaProteins to the highest score
 *
 */
void get_meta_scores(
		   MetaMapping* metaMapping, 
		   ProteinToScore* proteinToScore,
		   MetaToScore* metaToScore
		   ){
  carp(CARP_INFO, "Finding scores of meta proteins");
  for (MetaMapping::iterator meta_it = metaMapping->begin(); 
       meta_it != metaMapping->end(); ++meta_it ){
    MetaProtein proteins = (*meta_it).second;
    FLOAT_T top_score = -1.0;
    for (MetaProtein::iterator protein_it = proteins.begin();
	 protein_it != proteins.end(); ++protein_it){
      PROTEIN_T* protein = (*protein_it);
      FLOAT_T score = (*proteinToScore)[protein];
      top_score = max(score, top_score);
    }
    (*metaToScore).insert(make_pair(proteins, top_score));
  }

}


/**
 * Takes a mapping of MetaProteins to scores and returns 
 * a mapping of set of peptides to rank
 *
 */
void get_meta_ranks(
		  MetaToScore* metaToScore,
		  MetaToRank* metaToRank
		  ){
  carp(CARP_INFO, "Finding ranks of meta proteins");
  vector< pair<FLOAT_T, MetaProtein> > metaVector;
  for (MetaToScore::iterator meta_it = metaToScore->begin();
       meta_it != metaToScore->end(); ++meta_it){
    MetaProtein proteins = (*meta_it).first;
    FLOAT_T score = (*meta_it).second;
    metaVector.push_back(make_pair(score, proteins));
  }
  sort(metaVector.begin(), metaVector.end());
  reverse(metaVector.begin(), metaVector.end());
 
  int cur_rank = 1;
  for (vector< pair<FLOAT_T, MetaProtein> >::iterator 
	 vector_it = metaVector.begin();
       vector_it != metaVector.end(); ++vector_it){
    MetaProtein proteins = (*vector_it).second;
    (*metaToRank).insert(make_pair(proteins, cur_rank));
    cur_rank++;
  }

}



/**
 * Greedily finds a peptide to protein mapping where each
 * peptide is only mapped to a single meta-protein. 
 *
 * Would of been better to implement with priority queue w/
 * adjancency lists: O(n*log(n)) but input size should be
 * small enough where performance would not be an issue
 *
 */
void perform_parsimony_analysis(MetaMapping* metaMapping){
  carp(CARP_INFO, "Performing Greedy Parsimony analysis");
  MetaMapping result(compare_peptide_sets);
  vector< pair<PeptideSet, MetaProtein > > peps_vector;

  // get all meta mappings into a vector 
  for (MetaMapping::iterator meta_iter = (*metaMapping).begin();
       meta_iter != (*metaMapping).end(); ++meta_iter){
    peps_vector.push_back((*meta_iter));
  }
  
  
  // greedy algorithm to pick off the meta proteins with
  // most peptide mappings
  while (!peps_vector.empty()){
    sort(peps_vector.begin(), peps_vector.end(), compare_sets);
    pair<PeptideSet, MetaProtein> node = peps_vector.back();
    peps_vector.pop_back();
    if (node.first.size() == 0) break; // do not enter anything without peptide sizes
    result.insert(node);
    PeptideSet cur_peptides = node.first;
    // update the peptide sets for the rest of meta proteins
    for (vector< pair<PeptideSet, MetaProtein > >::iterator  
	   iter= peps_vector.begin();
	 iter != peps_vector.end(); ++iter){
      PeptideSet peptides = (*iter).first;
      PeptideSet difference(compare_pep);
      set_difference(peptides.begin(), peptides.end(), 
		     cur_peptides.begin(), cur_peptides.end(), 
		     inserter(difference, difference.end()), compare_pep);
      (*iter).first = difference;
    }
  }
  (*metaMapping) = result;
}



/** 
 * Outputs the protein id and their corresponding quantifcation
 * score in a tab delimited format. Outputs the meta protein
 * rank if parsimony was called.
 */
void print_protein_results(
			 ProteinToScore* proteinToScore,
			 MetaToRank* metaToRank,
			 ProteinToMetaProtein* proteinToMeta,
			 char * output_path
			 ){
  carp(CARP_INFO, "Outputting results");
  ofstream targetFile;
  targetFile.open(output_path);
  bool isParsimony = (proteinToMeta->size() != 0);


  vector<pair<FLOAT_T, PROTEIN_T*> > scoreToProtein;
  for (ProteinToScore::iterator it = proteinToScore->begin(); 
       it != proteinToScore->end(); ++it){
    PROTEIN_T* protein = it->first;
    FLOAT_T score = it->second;
    scoreToProtein.push_back(make_pair(score, protein));
  }
  
  /* sort and reverse the list */
  sort(scoreToProtein.begin(), scoreToProtein.end());
  reverse(scoreToProtein.begin(), scoreToProtein.end());

  /* write header */
  if (isParsimony){
    targetFile << "ProteinId\tScore\tParsimonyRank" << endl;
  } else {
    targetFile << "ProteinId\tScore" << endl;
  }

  for (vector<pair<FLOAT_T, PROTEIN_T*> >::iterator
	 it = scoreToProtein.begin(); 
       it != scoreToProtein.end(); ++it){
    FLOAT_T score = it->first;
    PROTEIN_T* protein = it->second;
    if (isParsimony){
      MetaProtein metaProtein = (*proteinToMeta)[protein];
      int rank = -1;
      if (metaToRank->find(metaProtein) != metaToRank->end()){
	rank = (*metaToRank)[metaProtein];
      } 
      targetFile << get_protein_id(protein) << "\t" 
		 << score << "\t" <<rank << endl;
    } else {
      targetFile << get_protein_id(protein) << "\t" 
		 << score << endl;
    }
  }
  targetFile.close();
}


/**
 * Outputs the peptide sequence and its corresponding
 * score in a tab delimited format
 */
void print_peptide_results(
			 PeptideToScore* peptideToScore,
			 char * output_path
			 ){
  carp(CARP_INFO, "Outputting results");
  ofstream targetFile;
  targetFile.open(output_path);
  
  vector<pair<FLOAT_T, PEPTIDE_T*> > scoreToPeptide;
  for (PeptideToScore::iterator it = peptideToScore->begin();
       it != peptideToScore->end(); ++it){
    PEPTIDE_T* peptide = it->first;
    FLOAT_T score = it->second;
    scoreToPeptide.push_back(make_pair(score, peptide));
  }
  
  targetFile << "PeptideSequence\tScore" << endl;
  
  sort(scoreToPeptide.begin(), scoreToPeptide.end());
  reverse(scoreToPeptide.begin(), scoreToPeptide.end());
  for (vector<pair<FLOAT_T, PEPTIDE_T*> >::iterator 
	 it = scoreToPeptide.begin();
       it != scoreToPeptide.end(); ++it){
    PEPTIDE_T* peptide = it->second;
    FLOAT_T score = it->first;
    targetFile << get_peptide_sequence(peptide) 
	       << "\t"  << score << endl;
  }
  targetFile.close();
}


/**
 * Removes peptides from the map if the peptide
 * sequence belongs in more than one protein
 */
void make_unique_mapping(
		       PeptideToScore* peptideToScore
		       ){
  carp(CARP_INFO, "Filtering peptides that have more"
       "than one protein source");
  for (PeptideToScore::iterator it = peptideToScore->begin();
       it != peptideToScore->end(); ++it){
    PEPTIDE_T* peptide = it->first;
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator =
      new_peptide_src_iterator(peptide);
    int num_proteins = 0;
    while(peptide_src_iterator_has_next(peptide_src_iterator)){
      peptide_src_iterator_next(peptide_src_iterator);
      if (++num_proteins > 1) break;
    }
    if (num_proteins > 1){
      peptideToScore->erase(it);
    }
  }
  
}




/* comparison and helper functions */

/** 
 *Comparison function for comparing set of peptides.
 */
bool compare_peptide_sets(PeptideSet set_one, PeptideSet set_two){
  PeptideSet pep_union(compare_pep);
  for (PeptideSet::iterator it = set_one.begin(); it != set_one.end();
       ++it){
    PEPTIDE_T* peptide = (*it);
    if (pep_union.find(peptide) == pep_union.end()){
      pep_union.insert((*it));
    }
  }
  for (PeptideSet::iterator it = set_two.begin(); it != set_two.end();
       ++it){
    if (pep_union.find((*it)) == pep_union.end()){
      pep_union.insert((*it));
    }
  }
  
  if (pep_union.size() == set_one.size() && set_one.size() == set_two.size()){
    return false;
  } else {
    string string_one = peps_to_string(set_one);
    string string_two = peps_to_string(set_two);
    return string_one.compare(string_two) > 0;
  }
}
/** 
 * Comparison function for comparing the equality of two proteins
 * based on thier protein ids'
 */
bool compare_prot(PROTEIN_T* protein_one, PROTEIN_T* protein_two){
  int compare = strcmp(get_protein_id(protein_one), get_protein_id(protein_two));
  if (compare == 0){
    return false;
  } else {
    return (compare > 0);
  }
}

/**
 * Comparison function for comparing the equality of two peptides
 * based on their sequence
 */
bool compare_pep(PEPTIDE_T* peptide_one, PEPTIDE_T* peptide_two){
  int compare = strcmp(get_peptide_sequence(peptide_one), get_peptide_sequence(peptide_two));
  if (compare == 0){
    return false;
  } else {
    return (compare > 0);
  }
}

/**
 * Comparison function for comparing the equality of MetaProteins
 * based off if they have the same set of proteins
 */
bool compare_meta_proteins(MetaProtein set_one, MetaProtein set_two){
  MetaProtein prot_union(compare_prot);
  for (MetaProtein::iterator it = set_one.begin(); it != set_one.end();
       ++it){
    if (prot_union.find((*it)) == prot_union.end()){
      prot_union.insert((*it));
    }
  }
  for (MetaProtein::iterator it = set_two.begin(); it != set_two.end();
       ++it){
    if (prot_union.find((*it)) == prot_union.end()){
      prot_union.insert((*it));
    }
  }
  
  if (prot_union.size() == set_one.size() && set_one.size() == set_two.size()){
    return false;
  } else {
    string string_one = meta_protein_to_string(set_one);
    string string_two = meta_protein_to_string(set_two);
    return string_one.compare(string_two) > 0;
  }
}

/**
 * returns a string of all peptide sequences separated
 * by space.
 */
string peps_to_string(PeptideSet s){
  stringstream ss (stringstream::in | stringstream::out);
  for (PeptideSet::iterator p_it = s.begin();
       p_it != s.end(); ++p_it){
    ss << get_peptide_sequence((*p_it)) << " ";
  }
  return ss.str();
}

/**
 * returns a string of all protein ids separated by
 * space
 */
string meta_protein_to_string(MetaProtein s){
  stringstream ss (stringstream::in | stringstream::out);
  for (MetaProtein::iterator p_it = s.begin();
       p_it != s.end(); ++p_it){
    ss << get_protein_id((*p_it)) << " ";
  }
  return ss.str();
}

/**
 *Comparison function that just compares the size
 * of PeptideSet
 */
bool compare_sets(
  pair<PeptideSet, MetaProtein > peps_one , 
  pair<PeptideSet, MetaProtein > peps_two){
  return ((peps_one).first.size() < (peps_two).first.size());
}

/**
 * Helper function to find the path of directory
 * for a given path to a file
 */
void getDirPath(char* path, char** dir){
  int pos = strlen(path)-1;
  while (pos >0 && *(path+pos) != '/'){pos--;}
  if (pos == 0){ // no backslashes exist
    *dir = (char*) malloc(sizeof(char)*2);
    strcpy(*dir, ".");
  } else {
    *dir = (char*) malloc((sizeof(char)*pos)+1);
    strncpy(*dir, path, pos);
    (*dir)[pos] = '\0';
  }
  carp(CARP_INFO, "directory: %s", *dir);
}
