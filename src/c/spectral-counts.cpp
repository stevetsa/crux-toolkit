#include "spectral-counts.h"

using namespace std;

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
void filter_matches(MATCH_COLLECTION_ITERATOR_T* match_collection_it,
                    set<MATCH_T*>& match_set);
void get_peptide_scores(
  set<MATCH_T*>&  matches, 
  PeptideToScore& peptideToScore
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
  MetaMapping& metaMapping
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
void getSpectra(map<pair<int,int>, Spectrum*>& spectra);
/* comparison function declarations */
bool compare_peptide_sets(PeptideSet, PeptideSet);
bool compare_meta_proteins(MetaProtein, MetaProtein);
bool sets_are_equal_size(
  pair<PeptideSet, MetaProtein>, 
  pair<PeptideSet, MetaProtein> 
);

/**
 * Given a collection of scored PSMs, print a list of proteins
 * ranked by their a specified score. Spectral-counts supports two
 * types of quantification: Normalized Spectral Abundance Factor (NSAF)
 * and Normalized Spectral Index (SIN). 
 * \returns 0 on successful completion.
 */
int spectral_counts_main(int argc, char** argv){

  const char* option_list[] = {
    "verbosity",
    "parameter-file",
    "parsimony",
    "threshold",
    "input-ms2",
    "fileroot",
    "output-dir",
    "overwrite",
    "unique-mapping",
    "input-bullseye",
    "quant-level",
    "measure"
  };
  const char* argument_list[] = {
    "input PSM",
    "protein database"
  };

  int num_options = sizeof(option_list) / sizeof(char*);
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize_run(SPECTRAL_COUNTS_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);

  // open output files
  OutputFiles output(SPECTRAL_COUNTS_COMMAND);
  output.writeHeaders();

  // get input file directory
  char* psm_file = get_string_parameter("input PSM");
  char** path_info = parse_filename_path(psm_file);
  if( path_info[1] == NULL ){
    path_info[1] = my_copy_string(".");
  }

  // create match collection
  int decoy_count = 0;
  char* database = get_string_parameter("protein database");
  MATCH_COLLECTION_ITERATOR_T* match_collection_it 
    = new_match_collection_iterator(path_info[1], database, &decoy_count);
   
  // get a set of matches that pass threshold
  set<MATCH_T*> matches;
  filter_matches(match_collection_it, matches);
  carp(CARP_INFO, "Number of matches passed the threshold %i", 
       matches.size());

  // get a set of peptides
  PeptideToScore peptideToScore(peptide_less_than);
  get_peptide_scores(matches, peptideToScore);
  if( get_boolean_parameter("unique-mapping") ){
    make_unique_mapping(&peptideToScore);
  }
  carp(CARP_INFO, "Number of Unique Peptides %i", peptideToScore.size());

  // quantify at either the peptide or protein level
  QUANT_LEVEL_TYPE_T quant_level =get_quant_level_type_parameter("quant-level");
  if( quant_level == PEPTIDE_QUANT_LEVEL ){ // peptide level
    normalize_peptide_scores(&peptideToScore);
    output.writeRankedPeptides(peptideToScore);

  } else if( quant_level == PROTEIN_QUANT_LEVEL ){ // protein level
    ProteinToScore proteinToScore(protein_id_less_than);
    ProteinToPeptides proteinToPeptides(protein_id_less_than);
    ProteinToMetaProtein proteinToMeta(protein_id_less_than);
    MetaMapping metaMapping(compare_peptide_sets);
    MetaToScore metaToScore(compare_meta_proteins);
    MetaToRank metaToRank(compare_meta_proteins);
    
    get_protein_scores(&peptideToScore, &proteinToScore);
    normalize_protein_scores(&proteinToScore);
    carp(CARP_INFO, "Number of Proteins %i", proteinToScore.size());
        
    PARSIMONY_TYPE_T parsimony = get_parsimony_type_parameter("parsimony");
    if( parsimony != PARSIMONY_NONE ){ //if parsimony is not none
      get_protein_to_peptides(&peptideToScore, &proteinToPeptides);
      get_meta_mapping(&proteinToPeptides, metaMapping);
      get_protein_to_meta_protein(&metaMapping, &proteinToMeta);
      carp(CARP_INFO, "Number of meta proteins %i", metaMapping.size());
      
      if( parsimony == PARSIMONY_GREEDY ){ //if parsimony is greedy
	perform_parsimony_analysis(&metaMapping);
      }
      get_meta_scores(&metaMapping, &proteinToScore, &metaToScore);
      get_meta_ranks(&metaToScore, &metaToRank);
    }
    
    output.writeRankedProteins(proteinToScore, metaToRank, proteinToMeta);

  } else {
    carp(CARP_FATAL, "Invalid quantification level.");
  }
  
  free(psm_file);
  free(path_info[1]);
  free(path_info[0]);
  free(path_info);
  free(database);

  return 0;
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
	PeptideSet newset(peptide_less_than);
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
    it->second = score / total / get_peptide_length(peptide);

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
 * For the spectrum associated with the match, sum the intensities of
 * all b and y ions that are not modified.
 * \return The sum of unmodified b and y ions.
 */
int sum_match_intensity(MATCH_T* match, 
                        map<pair<int,int>, Spectrum*>& spectra,
                        FLOAT_T bin_width)
{
  int match_intensity = 0;
  char* peptide_seq = get_match_sequence(match);
  MODIFIED_AA_T* modified_sequence = get_match_mod_sequence(match);
  int charge = get_match_charge(match);
  Spectrum* temp = get_match_spectrum(match);
  int scan = temp -> getFirstScan();
  Spectrum* spectrum = spectra[make_pair(scan, charge)];
  Ion* ion;
  SCORER_TYPE_T score_type = XCORR;
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
  delete ion_series;
  free(peptide_seq);
  
  return match_intensity;
}


/**
 * Generate a score for each peptide in the set of matches.  Store the
 * scores in the peptideToScore object.
 *
 * For SIN the score is the sum of intensites of b and y ions (without
 * H2O modifications).  Intensites are taken from the .ms2 file.
 *
 * For NSAF, the score is the number of matches for each peptide.
 */
void get_peptide_scores(
  set<MATCH_T*>& matches, ///< get peptides from these matches
  PeptideToScore& peptideToScore ///< store peptides and their scores here
){
  
  MEASURE_TYPE_T measure = get_measure_type_parameter("measure");
  FLOAT_T bin_width = get_double_parameter("mz-bin-width");
  map<pair<int,int>, Spectrum*> spectra;
  
    
  // for SIN, parse out spectrum collection from ms2 fiel
  if( measure == MEASURE_SIN ){ 
    getSpectra(spectra);
  }

  for(set<MATCH_T*>::iterator match_it = matches.begin();
      match_it != matches.end(); ++match_it){

    FLOAT_T match_intensity = 1; // for NSAF every match counted as 1

    MATCH_T* match = (*match_it);
    // for sin, calculate total ion intensity for match by
    // summing up peak intensities
    if( measure == MEASURE_SIN ){ 
      match_intensity = sum_match_intensity(match, spectra, bin_width);
    }

    // add ion_intensity to peptide scores
    PEPTIDE_T* peptide = get_match_peptide(match);
    if (peptideToScore.find(peptide) ==  peptideToScore.end()){
      peptideToScore.insert(make_pair(peptide, 0.0));
    } 
    peptideToScore[peptide] += match_intensity;
  }

}


/**
 * Parse spectra from the file given in the ms2-file parameter.  Store
 * spectra indexed by a scan-number,charge pair.
 */
void getSpectra(map<pair<int,int>, Spectrum*>& spectra){
  char* ms2 = get_string_parameter("input-ms2");
  SpectrumCollection* spectrum_collection = new SpectrumCollection(ms2);
  if (!spectrum_collection->parse()){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2);
  } 
  // get spectra from ms2 file
  FilteredSpectrumChargeIterator* spectrum_iterator = 
    new FilteredSpectrumChargeIterator(spectrum_collection);
  while (spectrum_iterator->hasNext()){
    int charge = 0;
    Spectrum* spectrum = spectrum_iterator->next(&charge);
    spectra.insert(make_pair(make_pair(spectrum->getFirstScan(), 
                                       charge), spectrum));
  }
  carp(CARP_INFO, "Number of Spectra %i", spectra.size());
  delete spectrum_iterator;
  free(ms2);
}

/**
 * Create a set of matches, all with an XCORR rank == 1 and all of which
 * have a qvalue score lower than user-specified threshold.
 */
void filter_matches(MATCH_COLLECTION_ITERATOR_T* match_collection_it, 
               set<MATCH_T*>& match_set )
{
  match_set.clear();
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  FLOAT_T threshold = get_double_parameter("threshold");
  bool qualify = false;
  while (match_collection_iterator_has_next(match_collection_it)){

    
    match_collection = match_collection_iterator_next(match_collection_it);
    match_iterator = new_match_iterator(match_collection, XCORR, TRUE);
    
    while(match_iterator_has_next(match_iterator)){
      MATCH_T* match = match_iterator_next(match_iterator);
      qualify = false;
      if (get_match_rank(match, XCORR) != 1){
	continue;
      }
      // find a qvalue score lower than threshold
      if (get_match_score(match, PERCOLATOR_QVALUE) != FLT_MIN &&
	  get_match_score(match, PERCOLATOR_QVALUE) <= threshold)  {
	qualify = true;
      } else if (get_match_score(match, QRANKER_QVALUE) != FLT_MIN &&
                 get_match_score(match, QRANKER_QVALUE) <= threshold)  {
	qualify = true;
      } else if (get_match_score(match, DECOY_XCORR_QVALUE) != FLT_MIN &&
		 get_match_score(match, DECOY_XCORR_QVALUE) <= threshold)  {
	qualify = true;
      } 

      if (qualify == true){
        match_set.insert(match);
      }
    } // next match
  } // next file
}




/**
 * Fills in the MetaMapping with entries of set of 
 * peptides that can be found in every protein in
 * the meta protein
 *
 */
void get_meta_mapping(
		    ProteinToPeptides* proteinToPeptides,
		    MetaMapping& metaMapping
		    ){
  carp(CARP_INFO, "Creating a mapping of meta protein to peptides");
  int count = 0;
  for (ProteinToPeptides::iterator prot_it = proteinToPeptides->begin();
       prot_it != proteinToPeptides->end(); ++prot_it){
    PROTEIN_T* protein = prot_it->first;
    PeptideSet pep_set = prot_it->second;

    if (metaMapping.find(pep_set) == metaMapping.end()){
      MetaProtein meta_protein(protein_id_less_than);
      count++;
      metaMapping.insert(make_pair(pep_set, meta_protein));
    }
    metaMapping[pep_set].insert(protein);
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
 * Greedily finds a peptide-to-protein mapping where each
 * peptide is only mapped to a single meta-protein. 
 *
 * Would of been better to implement with priority queue w/
 * adjancency lists: O(n*log(n)) but input size should be
 * small enough that performance should not be an issue.
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
    sort(peps_vector.begin(), peps_vector.end(), sets_are_equal_size);
    pair<PeptideSet, MetaProtein> node = peps_vector.back();
    peps_vector.pop_back();
    if (node.first.size() == 0){ break; }// do not enter anything without peptide sizes
    result.insert(node);
    PeptideSet cur_peptides = node.first;
    // update the peptide sets for the rest of meta proteins
    for (vector< pair<PeptideSet, MetaProtein > >::iterator  
	   iter= peps_vector.begin();
	 iter != peps_vector.end(); ++iter){
      PeptideSet peptides = (*iter).first;
      PeptideSet difference(peptide_less_than);
      set_difference(peptides.begin(), peptides.end(), 
		     cur_peptides.begin(), cur_peptides.end(), 
		     inserter(difference, difference.end()), 
                     peptide_less_than);
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
                           char * output_path)
{
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
			 char* output_path
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
    char* seq = get_peptide_sequence(peptide);
    targetFile << seq << "\t"  << score << endl;
    free(seq);
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
    int num_proteins = get_peptide_num_peptide_src(peptide);
    if (num_proteins > 1){
      peptideToScore->erase(it);
    }
  }
  
}




/* comparison and helper functions */

/** 
 * Compare two sets of peptides and return true if the first unshared
 * peptide sequence in set one is lexically less than that in set
 * two.
 */
bool compare_peptide_sets(PeptideSet set_one, PeptideSet set_two){

  // compare each peptides in the two (sorted) sets
  PeptideSet::iterator iter1 = set_one.begin();
  PeptideSet::iterator iter2 = set_two.begin();

  while( iter1 != set_one.end() && iter2 != set_two.end() ){
    int diff = tri_compare_peptide_sequence(*iter1, *iter2);
    if( diff < 0 ){
      return true;
    } else if(diff > 0 ){
      return false;
    } // else they are equal, compare the next

    ++iter1;
    ++iter2;
  }

  // all peptides were the same; are the sets the same size?
  if( set_one.size() == set_two.size() ){
    return false;
  } else if( set_one.size() > set_two.size() ){
    return false;
  } else { // one < two
    return true;
  }
}

/**
 * Comparison function for MetaProteins.  MetaProtein one is less than
 * MetaProtein two if the first non-matching protein id of one is less than
 * that of two.  
 * \returns True if one < two, false if one == two or one > two.
 */
bool compare_meta_proteins(MetaProtein set_one, MetaProtein set_two){

  // compare each protein in the two (sorted) sets
  MetaProtein::iterator iter1 = set_one.begin();
  MetaProtein::iterator iter2 = set_two.begin();

  while( iter1 != set_one.end() && iter2 != set_two.end() ){
    bool one_less_than_two = protein_id_less_than(*iter1, *iter2);
    bool two_less_than_one = protein_id_less_than(*iter1, *iter2);
    // different proteins one is less than the other
    if( one_less_than_two || two_less_than_one ){
      return one_less_than_two;
    }
    // else, they are the same, keep comparing
    ++iter1;
    ++iter2;
  } 

  // all proteins were the same, are the sets the same size?
  if( set_one.size() == set_two.size() ){
    return false;
  } else if( set_one.size() > set_two.size() ){
    return false;
  } else { // one < two
    return true;
  }
}

/**
 * Compare the size of the PeptideSets in the two given pairs.
 * \returns True if the PeptideSets are the same size, else false.
 */
bool sets_are_equal_size(
  pair<PeptideSet, MetaProtein > peps_one , 
  pair<PeptideSet, MetaProtein > peps_two){
  return ((peps_one).first.size() < (peps_two).first.size());
}

