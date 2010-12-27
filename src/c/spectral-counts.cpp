#include "spectral-counts.h"

using namespace std;

typedef map<PEPTIDE_T*, FLOAT_T, bool(*)(PEPTIDE_T*, PEPTIDE_T*)> PeptideToScore;
typedef map<PROTEIN_T*, FLOAT_T, bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToScore;
typedef set<PROTEIN_T*, bool(*)(PROTEIN_T*, PROTEIN_T*)> MetaProtein;
typedef set<PEPTIDE_T*, bool(*)(PEPTIDE_T*, PEPTIDE_T*)> PeptideSet;
typedef map<PeptideSet, MetaProtein, bool(*)(PeptideSet, PeptideSet) > MetaMapping;
typedef map<PROTEIN_T*, PeptideSet , bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToPeptides;
typedef map<MetaProtein, int, bool(*)(MetaProtein, MetaProtein)> MetaToRank;
typedef map<MetaProtein, FLOAT_T, bool(*)(MetaProtein, MetaProtein)> MetaToScore;
typedef map<PROTEIN_T*, MetaProtein, bool(*)(PROTEIN_T*, PROTEIN_T*)> ProteinToMetaProtein;

string pepsToString(PeptideSet s);
void getDirPath(char* path, char** dir);
set<MATCH_T*> filterMatches(MATCH_COLLECTION_ITERATOR_T* match_collection_it);
void getPeptideScores(set<MATCH_T*>  matches, PeptideToScore* peptideToScore);
void getProteinScores(PeptideToScore peptideToScore, ProteinToScore* proteinToScore);
void getProteinToPeptides(PeptideToScore peptideToScore, ProteinToPeptides* proteinToPeptides);
void getProteinToMetaProtein(MetaMapping, ProteinToMetaProtein* proteinToMetaProtein);
void getMetaMapping(ProteinToPeptides proteinToPeptides, MetaMapping* metaMapping);
void getMetaRanks(MetaToScore metaToScore, MetaToRank* metaToRank);
void getMetaScores(MetaMapping metaMapping, ProteinToScore proteinToScore, MetaToScore* metaToScore);
void performParsimonyAnalysis(MetaMapping* metaMapping);
void normalizePeptideScores(PeptideToScore* peptideToScore);
void normalizeProteinScores(ProteinToScore* proteinToScore);

bool compare_pep(PEPTIDE_T*, PEPTIDE_T*);
bool compare_prot(PROTEIN_T*, PROTEIN_T*);
bool compare_peptide_sets(PeptideSet, PeptideSet);
bool compare_meta_proteins(MetaProtein, MetaProtein);
string metaProteinToString(MetaProtein s);
string pepsToString(PeptideSet s);

bool compare_sets(pair<PeptideSet, MetaProtein > , pair<PeptideSet, MetaProtein >);
bool compare_sets(pair<PeptideSet, MetaProtein > peps_one , pair<PeptideSet, MetaProtein > peps_two){
  return ((peps_one).first.size() < (peps_two).first.size());
}


int spectral_counts_main(int argc, char** argv){
  const char* option_list[] = {
    "threshold",
    "input-ms2",
    "fileroot",
    "output-dir",
    "overwrite",
    "unique-mapping",
    "input-bullseye",
    "quant-level",
    "measure",
    "average",
    "version",
    "verbosity"
  };
  const char* argument_list[] = {
    "input-PSM",
    "database"
  };

  int num_options = sizeof(option_list) / sizeof(char*);
  int num_arguments = sizeof(argument_list) / sizeof(char*);


  initialize_run(SPECTRAL_COUNTS_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);

  /* Create Match Collection from input-PSM */
  char * psm_file = get_string_parameter("input-PSM");
  char * database = get_string_parameter("database");
  char * dir_path = NULL;
  char * quant_level = get_string_parameter("quant-level");

  if ( strcmp(quant_level,"PROTEIN") !=0 && strcmp(quant_level,"PEPTIDE") != 0 ){
    carp(CARP_FATAL, "Must pick PROTEIN or PEPTIDE for quant-level option");
  }

  getDirPath(psm_file, &dir_path);
  
  int decoy_count = 0;

  MATCH_COLLECTION_ITERATOR_T* match_collection_it 
    = new_match_collection_iterator(dir_path, database, &decoy_count);
 

  set<MATCH_T*> matches = filterMatches(match_collection_it);
  carp(CARP_INFO, "Number of matches passed the threshold %i", 
       matches.size());
  
  

  
  PeptideToScore peptideToScore(compare_pep);
  ProteinToScore proteinToScore(compare_prot);
  getPeptideScores(matches, &peptideToScore);
  carp(CARP_INFO, "Number of Unique Peptides %i", peptideToScore.size());

  // if peptide_score
  //normalizePeptideScores(&peptideToScore);

  
  //if protein score
  if (strcmp(quant_level,"PROTEIN") == 0){
    ProteinToScore proteinToScore(compare_prot);
    ProteinToPeptides proteinToPeptides(compare_prot);
    ProteinToMetaProtein proteinToMeta(compare_prot);
    MetaMapping metaMapping(compare_peptide_sets);
    MetaToScore metaToScore(compare_meta_proteins);
    MetaToRank metaToRank(compare_meta_proteins);
    
    
    getProteinScores(peptideToScore, &proteinToScore);
    normalizeProteinScores(&proteinToScore);
    carp(CARP_INFO, "Number of Proteins %i", proteinToScore.size());
    
    
    // if parsimony is not none{
    getProteinToPeptides(peptideToScore, &proteinToPeptides);
    getMetaMapping(proteinToPeptides, &metaMapping);
    getProteinToMetaProtein(metaMapping, &proteinToMeta);
    carp(CARP_INFO, "Number of meta proteins %i", metaMapping.size());
    
        
    // if parsimony is greedy{
    performParsimonyAnalysis(&metaMapping);
    //}
    
    getMetaScores(metaMapping, proteinToScore, &metaToScore);
    getMetaRanks(metaToScore, &metaToRank);
    // }
    
    
    for (ProteinToScore::iterator it = proteinToScore.begin(); 
	 it != proteinToScore.end(); it++){
      PROTEIN_T* protein = (*it).first;
      FLOAT_T score = (*it).second;
      MetaProtein metaProtein = proteinToMeta[protein];
      int rank = -1;
      FLOAT_T meta_score = -100.0;
      if (metaToRank.find(metaProtein) != metaToRank.end()){
	rank = metaToRank[metaProtein];
      }
      if (metaToScore.find(metaProtein) != metaToScore.end()){
	meta_score = metaToScore[metaProtein];
      }
      cout << get_protein_id(protein) << "\t" << score << "\t"  << meta_score << "\t" << rank << endl;
    }
    
  }     

  

  free(dir_path);
  return 1;
}


/*
 * Creates a mapping from protein to a set of peptides that is part
 * of the protein. The peptides only include the peptides with scores
 *
 */
void getProteinToPeptides(PeptideToScore peptideToScore, ProteinToPeptides* proteinToPeptides){
  for (PeptideToScore::iterator pep_it = peptideToScore.begin();
       pep_it != peptideToScore.end(); pep_it++){
    PEPTIDE_T* peptide = (*pep_it).first;
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator =
      new_peptide_src_iterator(peptide);
    while( peptide_src_iterator_has_next(peptide_src_iterator)) {
      PEPTIDE_SRC_T* peptide_src = peptide_src_iterator_next(peptide_src_iterator);
      PROTEIN_T* protein = get_peptide_src_parent_protein(peptide_src);
      if ((*proteinToPeptides).find(protein) == (*proteinToPeptides).end()){
	PeptideSet newset(compare_pep);
	(*proteinToPeptides).insert(make_pair(protein, newset));
      }
      (*proteinToPeptides)[protein].insert(peptide);
    }
    free(peptide_src_iterator);   
  }
}


void getProteinToMetaProtein(
			     MetaMapping metaMap, 
			     ProteinToMetaProtein* proteinToMeta
			     ){
  for (MetaMapping::iterator meta_protein_it = metaMap.begin();
       meta_protein_it != metaMap.end(); meta_protein_it++){
    MetaProtein proteins = (*meta_protein_it).second;
    for (MetaProtein::iterator proteins_it = proteins.begin();
	 proteins_it != proteins.end(); proteins_it++){
      (*proteinToMeta).insert(make_pair((*proteins_it), proteins));
    }
  }
  /*
  for (ProteinToMetaProtein::iterator it = (*proteinToMeta).begin();
       it != (*proteinToMeta).end(); it++){
    PROTEIN_T* protein = (*it).first;
    MetaProtein proteins = (*it).second;
    cout << get_protein_id(protein) << "||"<<metaProteinToString(proteins) << endl;
  }
  */

}


/*
 * For every protein that can be mapped from the collection of
 * peptides in the peptideToScore mapping, the sum of all scores
 * from mapped peptides are calculated and returned as a map
 *
 */
void getProteinScores(
		      PeptideToScore peptideToScore,
		      ProteinToScore* proteinToScore
		      ){
  
  // iterate through each peptide
  for (PeptideToScore::iterator pep_it = peptideToScore.begin();
       pep_it != peptideToScore.end(); pep_it++){
    PEPTIDE_T* peptide = (*pep_it).first;
    FLOAT_T pep_score = (*pep_it).second;
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator =
      new_peptide_src_iterator(peptide);
    while( peptide_src_iterator_has_next(peptide_src_iterator)) {
      PEPTIDE_SRC_T* peptide_src = peptide_src_iterator_next(peptide_src_iterator);
      PROTEIN_T* protein = get_peptide_src_parent_protein(peptide_src);
      if ((*proteinToScore).find(protein) == (*proteinToScore).end()){
	(*proteinToScore).insert(make_pair(protein, 0.0));
      }
      (*proteinToScore)[protein] += pep_score;
    }
    free(peptide_src_iterator);
  }

}


/*
 * Takes the PeptideToScores map and updates all
 * values with normalized values. 
 *
 */
void normalizePeptideScores(PeptideToScore* peptideToScore){
  carp(CARP_INFO, "Normalizing peptide scores");
  FLOAT_T total = 0.0;

  // normalize by peptide length
  for (PeptideToScore::iterator it = (*peptideToScore).begin();
       it != (*peptideToScore).end(); it++){
    PEPTIDE_T* peptide = (*it).first;
    FLOAT_T score = (*it).second;
    FLOAT_T normalized_score = score / strlen(get_peptide_sequence(peptide));
    (*it).second = normalized_score;
    total += normalized_score;
  }
  
  // normalize by sum of scores
  for (PeptideToScore::iterator it = (*peptideToScore).begin();
       it != (*peptideToScore).end(); it++){
    FLOAT_T score = (*it).second;
    (*it).second = score / total;
  }
  
}


/*
 * Takes ProteinToScore mapping and updates all the scores
 * with normalized values
 *
 */

void normalizeProteinScores(ProteinToScore* proteinToScore){
  carp(CARP_INFO, "Normalizing protein scores");
  FLOAT_T total = 0.0;
  for (ProteinToScore::iterator it = (*proteinToScore).begin();
       it != (*proteinToScore).end(); it++){
    PROTEIN_T* protein = (*it).first;
    FLOAT_T score = (*it).second;
    FLOAT_T normalized = score / strlen(get_protein_sequence(protein));
    (*it).second = normalized;
    total += normalized;
  }

  for (ProteinToScore::iterator it = (*proteinToScore).begin();
       it != (*proteinToScore).end(); it++){
    FLOAT_T score = (*it).second;
    (*it).second = score / total;
  }

}


/*
 * For SIN, it will parse the spectras from ms2 to
 * find the total ion intensity from b and y ions without
 * h20 modifications. These ion intensities are summed
 * for each peptide 
 *
 * For NSAF, it will sum up the counts of matches for each
 * peptide
 *
 */
void getPeptideScores(
		      set<MATCH_T*> matches,
		      PeptideToScore* peptideToScore
		      ){
  
  char * measure = get_string_parameter("measure"); //TODO check this value
  char * ms2 = get_string_parameter("input-ms2");
  FLOAT_T bin_width = get_double_parameter("mz-bin-width");
  BOOLEAN_T isSin = (measure == NULL || !strcmp(measure, "SIN"));
  map<pair<int,int>, Spectrum*> spectras;
  
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

  for(set<MATCH_T*>::iterator match_it = matches.begin();
      match_it != matches.end(); match_it++){
    int match_intensity = 0;
    MATCH_T* match = (*match_it);
    if (isSin){ // for sin, calculate total ion intensity for match
      char* peptide_seq = get_match_sequence(match);
      int charge = get_match_charge(match);
      Spectrum* temp = get_match_spectrum(match);
      int scan = temp -> getFirstScan();
      // FREE temp
      Spectrum* spectrum = spectras[make_pair(scan, charge)];
      Ion* ion;
      IonSeries* ion_series
	= new IonSeries(peptide_seq, charge, new IonConstraint()); //What to do for ion constraint?
      for (IonIterator ion_it = ion_series->begin(); 
	   ion_it != ion_series->end(); ion_it++){
	ion = (*ion_it);
	if (ion -> getType() == B_ION || ion -> getType() == Y_ION){
	  // how to check if ion has H20 modification?
	  if (!ion->isModified()){
	    PEAK_T* peak = spectrum->getNearestPeak(ion->getMassZ(),
						    bin_width);
	    if (peak != NULL){
	      match_intensity += get_peak_intensity(peak);
	      carp(CARP_INFO, "peak intensity %i", match_intensity);
	    }
	  }
	}
      }
    } else { // for nsaf every match is counted as 1 
      match_intensity = 1;
    }
    // add ion_intensity to peptide scores
    PEPTIDE_T* peptide = get_match_peptide(match);
    if ((*peptideToScore).find(peptide) ==  (*peptideToScore).end()){
      (*peptideToScore).insert(make_pair(peptide, 0.0));
    } 
    (*peptideToScore)[peptide]+=match_intensity;
  }
}




/*
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
}


set<MATCH_T*> filterMatches(MATCH_COLLECTION_ITERATOR_T* match_collection_it){
  set<MATCH_T*> match_set;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  FLOAT_T threshold = get_double_parameter("threshold");
  BOOLEAN_T qualify = FALSE;
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
      } else if (get_match_score(match, QRANKER_QVALUE) != FLT_MIN ||
	       get_match_score(match, QRANKER_QVALUE) <= threshold)  {
	qualify = TRUE;
      } else if (get_match_score(match, DECOY_XCORR_QVALUE) != FLT_MIN ||
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



/** Methods to help with Parsimony Analysis **/


/*
 * Creates a mapping of set of peptides to a set of proteins that 
 * map to that exact set of peptides. The set of proteins can be
 * considered to be meta-proteins
 *
 */
void getMetaMapping(
		    ProteinToPeptides proteinToPeptides,
		    MetaMapping* metaMapping
		    ){
  carp(CARP_INFO, "Creating a mapping of meta protein to peptides");
  int count = 0;
  for (ProteinToPeptides::iterator prot_it = proteinToPeptides.begin();
       prot_it != proteinToPeptides.end(); prot_it++){
    PROTEIN_T* protein = (*prot_it).first;
    PeptideSet pep_set = (*prot_it).second;

    if ((*metaMapping).find(pep_set) == (*metaMapping).end()){
      MetaProtein meta_protein(compare_prot);
      count++;
      (*metaMapping).insert(make_pair(pep_set, meta_protein));
    }
    (*metaMapping)[pep_set].insert(protein);
  }

}


/*
 * Takes a mapping of set of peptides to meta proteins
 * a mapping of protein to scores, and finds the largest score
 * of the meta proteins for each protein. The object returned
 * is a mapping of set of peptides to the average score
 *
 */
void getMetaScores(
		   MetaMapping metaMapping, 
		   ProteinToScore proteinToScore,
		   MetaToScore* metaToScore
		   ){
  carp(CARP_INFO, "Finding scores of meta proteins");
  for (MetaMapping::iterator meta_it = metaMapping.begin(); 
       meta_it != metaMapping.end(); meta_it++ ){
    MetaProtein proteins = (*meta_it).second;
    FLOAT_T top_score = -1.0;
    for (MetaProtein::iterator protein_it = proteins.begin();
	 protein_it != proteins.end(); protein_it++){
      PROTEIN_T* protein = (*protein_it);
      FLOAT_T score = proteinToScore[protein];
      top_score = max(score, top_score);
    }
    (*metaToScore).insert(make_pair(proteins, top_score));
  }
  /*
  for (MetaToScore::iterator it = (*metaToScore).begin();
       it != (*metaToScore).end(); it++){
    FLOAT_T score = (*it).second;
    MetaProtein proteins = (*it).first;
    cout << metaProteinToString(proteins) << "\t" << score << endl;
  }
  */

}


/*
 * Takes a mapping of set of peptides to scores and returns 
 * a mapping of set of peptides to rank
 *
 */
void getMetaRanks(
		  MetaToScore metaToScore,
		  MetaToRank* metaToRank
		  ){
  carp(CARP_INFO, "Finding ranks of meta proteins");
  vector< pair<FLOAT_T, MetaProtein> > metaVector;
  for (MetaToScore::iterator meta_it = metaToScore.begin();
       meta_it != metaToScore.end(); meta_it++ ){
    MetaProtein proteins = (*meta_it).first;
    FLOAT_T score = (*meta_it).second;
    metaVector.push_back(make_pair(score, proteins));
  }
  sort(metaVector.begin(), metaVector.end());
  reverse(metaVector.begin(), metaVector.end());
 
  int cur_rank = 1;
  for (vector< pair<FLOAT_T, MetaProtein> >::iterator 
	 vector_it = metaVector.begin();
       vector_it != metaVector.end(); vector_it++){
    MetaProtein proteins = (*vector_it).second;
    (*metaToRank).insert(make_pair(proteins, cur_rank));
    cur_rank++;
  }
  
  /*
  for (MetaToRank::iterator it = (*metaToRank).begin();
       it != (*metaToRank).end(); it++){
    int rank= (*it).second;
    MetaProtein proteins = (*it).first;
    cout << metaProteinToString(proteins) << "\t" << rank << endl;
  }
  */

}



/*
 * Greedily finds a peptide to protein ranking where each
 * peptide is only mapped to a single meta-protein. 
 *
 * Would of been better to implement with priority queue w/
 * adjancency lists: O(n*log(n)) but input size should be
 * small enough where performance would not be an issue
 *
 */
void performParsimonyAnalysis(MetaMapping* metaMapping){
  carp(CARP_INFO, "performing Parsimony analysis");
  MetaMapping result(compare_peptide_sets);
  vector< pair<PeptideSet, MetaProtein > > peps_vector;

  // get all meta mappings into a vector 
  for (MetaMapping::iterator meta_iter = (*metaMapping).begin();
       meta_iter != (*metaMapping).end(); meta_iter++){
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
	 iter != peps_vector.end(); iter++){
      PeptideSet peptides = (*iter).first;
      PeptideSet difference(compare_pep);
      set_difference(peptides.begin(), peptides.end(), 
		     cur_peptides.begin(), cur_peptides.end(), 
		     inserter(difference, difference.end()), compare_pep);
      //if (difference.size() != peptides.size()) carp(CARP_INFO, "Difference size: %i Orignal %i", difference.size(), peptides.size());
      (*iter).first = difference;
    }
  }
  (*metaMapping) = result;
}




/* comparison and helper functions */



bool compare_peptide_sets(PeptideSet set_one, PeptideSet set_two){
  PeptideSet pep_union(compare_pep);
  for (PeptideSet::iterator it = set_one.begin(); it != set_one.end();
       it++){
    PEPTIDE_T* peptide = (*it);
    if (pep_union.find(peptide) == pep_union.end()){
      pep_union.insert((*it));
    }
  }
  for (PeptideSet::iterator it = set_two.begin(); it != set_two.end();
       it++){
    if (pep_union.find((*it)) == pep_union.end()){
      pep_union.insert((*it));
    }
  }
  
  if (pep_union.size() == set_one.size() && set_one.size() == set_two.size()){
    return false;
  } else {
    string string_one = pepsToString(set_one);
    string string_two = pepsToString(set_two);
    return string_one.compare(string_two) > 0;
  }
}

bool compare_prot(PROTEIN_T* protein_one, PROTEIN_T* protein_two){
  int compare = strcmp(get_protein_id(protein_one), get_protein_id(protein_two));
  if (compare == 0){
    return false;
  } else {
    return (compare > 0);
  }
}

bool compare_pep(PEPTIDE_T* peptide_one, PEPTIDE_T* peptide_two){
  int compare = strcmp(get_peptide_sequence(peptide_one), get_peptide_sequence(peptide_two));
  if (compare == 0){
    return false;
  } else {
    return (compare > 0);
  }
}

bool compare_meta_proteins(MetaProtein set_one, MetaProtein set_two){
  MetaProtein prot_union(compare_prot);
  for (MetaProtein::iterator it = set_one.begin(); it != set_one.end();
       it++){
    if (prot_union.find((*it)) == prot_union.end()){
      prot_union.insert((*it));
    }
  }
  for (MetaProtein::iterator it = set_two.begin(); it != set_two.end();
       it++){
    if (prot_union.find((*it)) == prot_union.end()){
      prot_union.insert((*it));
    }
  }
  
  if (prot_union.size() == set_one.size() && set_one.size() == set_two.size()){
    return false;
  } else {
    string string_one = metaProteinToString(set_one);
    string string_two = metaProteinToString(set_two);
    return string_one.compare(string_two) > 0;
  }
}


string pepsToString(PeptideSet s){
  stringstream ss (stringstream::in | stringstream::out);
  for (PeptideSet::iterator p_it = s.begin();
       p_it != s.end(); p_it++){
    ss << get_peptide_sequence((*p_it)) << " ";
  }

  return ss.str();
}


string metaProteinToString(MetaProtein s){
  stringstream ss (stringstream::in | stringstream::out);
  for (MetaProtein::iterator p_it = s.begin();
       p_it != s.end(); p_it++){
    ss << get_protein_id((*p_it)) << " ";
  }
  return ss.str();
}
