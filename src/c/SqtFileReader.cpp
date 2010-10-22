#include "SqtFileReader.h"

#include <fstream>
#include <sstream>

#include "objects.h"
#include "match.h"
#include "spectrum.h"
#include "peptide.h"
#include "carp.h"

using namespace std;


/**
 * \returns a blank SqtFileReader object
 */
SqtFileReader::SqtFileReader() : has_next_(false){
  file_ptr_ = NULL;
}


SqtFileReader::SqtFileReader(
  const char* file_name,
  DATABASE_T* database
) : has_next_(false){ 
  file_ptr_ = NULL;
  loadData(file_name);
  database_ = database;
}

SqtFileReader::SqtFileReader(
  const string& file_name,
  DATABASE_T* database
) : has_next_(false){
  file_ptr_ = NULL;
  loadData(file_name.c_str());
  database_ = database;
}

/**
 * Destructor
 */
SqtFileReader::~SqtFileReader(){
  
  if (file_ptr_ != NULL){
    file_ptr_ -> close();
    delete file_ptr_;
  }
  /*
  if (spectrum_ != NULL){
    free_spectrum(spectrum_);
  }
  */
}


/**
 * Loads the sqt file specified by the file_name and prepares the
 * iterator to iterate through the matches
 */
void SqtFileReader::loadData(
  const char* file_name		    
  ){
  file_name_ = string(file_name);
  
  // delete file pointer if one already exists
  if (file_ptr_ != NULL){
    file_ptr_ -> close();
    delete file_ptr_;
  }

  file_ptr_ = new fstream(file_name, ios::in);
  
  if (!file_ptr_ -> is_open()){
    carp(CARP_FATAL, "Opening %s or reading failed", file_name);
    return;
  }
  
  // Read all header lines
  has_next_ = getline(*file_ptr_, next_data_string_) != NULL;
  while (has_next_ && next_data_string_[0] != 'S'){
    if (next_data_string_[0] != 'H'){
      carp(CARP_FATAL, "Reading %s failed", file_name);
      return;
    }
    has_next_ = getline(*file_ptr_, next_data_string_) != NULL;
  }
  next();
}


/**
 * Updates the variables with the next match on iteration
 * and does nothing if iterated through entire file. Updates 
 * spectrum specific variables as it iterates.
 */
BOOLEAN_T  SqtFileReader::next(){
  BOOLEAN_T okay = false;
  if (has_next_){
    while (has_next_ && next_data_string_[0] == 'S'){ // new spectrum, load the S line values
      getSLine();
      has_next_ = (getline(*file_ptr_, next_data_string_) != NULL);
    } 
    if (has_next_ && next_data_string_[0] == 'M'){ // new match, load the M line values
      getMLine();
      has_next_ = (getline(*file_ptr_, next_data_string_) != NULL);
      okay = true;
    } else { 
      // In this case if last spectrum of the file did not have 
      // any matches
      okay = false;
    }
    
    protein_ids_.clear();
    while (has_next_  && next_data_string_[0] == 'L'){
      getLLine();
      has_next_ = (getline(*file_ptr_, next_data_string_) != NULL);
    }
  }
  return okay;
}



/**
 * Updates the protein variables with the next line string
 */
void SqtFileReader::getLLine(){
  istringstream linestr;
  string protein_id;
  linestr.str(next_data_string_);
  if (!(linestr >> letter_ >> protein_id)){
    carp(CARP_FATAL, "Could not parse the L line: %s", 
	 next_data_string_.c_str());
  }
  protein_ids_.push_back(protein_id);
}



/**
 * Updates the Spectrum variables with the next line string
 */
void SqtFileReader::getSLine(){
  istringstream linestr;
  linestr.str(next_data_string_);
  if (!(linestr >> letter_ >> low_scan_ >> high_scan_ >> charge_ >> process_time_ 
	>> server_ >> experimental_mass_ >> total_ion_intensity_
	>> lowest_sp_ >> num_seq_matched_)){
    carp(CARP_FATAL, "Could not parse the S line: %s", next_data_string_.c_str());
  }
  linestr.clear();
}

/**
 * Updates the Match variables with the next line string
 */
void SqtFileReader::getMLine(){
  istringstream linestr;
  linestr.str(next_data_string_);
  if (!(linestr >> letter_ >> x_corr_rank_ >> sp_rank_ >> calculated_mass_ >> delta_cn_
	>> xcorr_ >> sp_ >> matched_ions_ >> expected_ions_ 
	>> sequence_matched_ >> validation_status_)){
    carp(CARP_FATAL, "Could not parse the M line: %s", next_data_string_.c_str());
  }
  linestr.clear();
}


/**
 * resets the file pointer to the beginning of the file
 */
void SqtFileReader::reset(){
  if (file_ptr_ != NULL){
    file_ptr_ -> close();
    delete file_ptr_;
  }
  loadData(file_name_.c_str());
}


/**
 * \returns whether there are more rows to 
 * iterate through
 */
BOOLEAN_T SqtFileReader::hasNext(){
  return has_next_;
}



void SqtFileReader::fillMatch(
 MATCH_T* match
){
  SPECTRUM_T* spectrum = create_spectrum_sqt(low_scan_, high_scan_, experimental_mass_);
  PEPTIDE_T* peptide = parse_peptide_sqt(database_, sequence_matched_, calculated_mass_, protein_ids_, FULL_DIGEST, true);  
  match = create_match_sqt(spectrum, peptide, x_corr_rank_, sp_rank_, xcorr_, sp_, matched_ions_, expected_ions_);
}


/*
 *Getters
 */

int SqtFileReader::getLowScan(){
  return low_scan_;
}

int SqtFileReader::getHighScan(){
  return high_scan_;
}

int SqtFileReader::getCharge(){
  return charge_;
}

double SqtFileReader::getProcessTime(){
  return process_time_;
}

string SqtFileReader::getServer(){
  return server_;
}

double SqtFileReader::getExperimentalMass(){
  return experimental_mass_;
}

double SqtFileReader::getTotalIonIntensity(){
  return total_ion_intensity_;
}

double SqtFileReader::getLowestSp(){
  return lowest_sp_;
}

int SqtFileReader::getNumSeqMatched(){
  return num_seq_matched_;
}

int SqtFileReader::getXCorrRank(){
  return x_corr_rank_;
}

int SqtFileReader::getSpRank(){
  return sp_rank_;
}

double SqtFileReader::getCalculatedMass(){
  return calculated_mass_;
}

double SqtFileReader::getDeltaCn(){
  return delta_cn_;
}

double SqtFileReader::getXCorr(){
  return xcorr_;
}

double SqtFileReader::getSp(){
  return sp_;
}

int SqtFileReader::getMatchedIons(){
  return matched_ions_;
}

int SqtFileReader::getExpectedIons(){
  return expected_ions_;
}

string SqtFileReader::getSequenceMatched(){
  return sequence_matched_;
}

string SqtFileReader::getValidationStatus(){
  return validation_status_;
}

vector<string> SqtFileReader::getProteinIds(){
  return protein_ids_;
}
