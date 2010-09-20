#ifndef SQTFILEREADER_H
#define SQTFILEREADER_H

#include <set>
#include <sstream>
#include <vector>

#include "parameter.h"

class SqtFileReader{
  
 protected:
  std::string letter_;
  std::string next_data_string_;
  std::vector<std::string> data_;
  std::fstream* file_ptr_;
  std::string file_name_;
  BOOLEAN_T okay_;
  BOOLEAN_T has_next_;
  // S line
  int low_scan_;
  int high_scan_;
  int charge_;
  double process_time_;
  std::string server_;
  double experimental_mass_;
  double total_ion_intensity_;
  double lowest_sp_;
  int num_seq_matched_;
  // M line
  int x_corr_rank_;
  int sp_rank_;
  double calculated_mass_;
  double delta_cn_;
  double xcorr_;
  double sp_;
  int matched_ions_;
  int expected_ions_;
  std::string sequence_matched_;
  std::string validation_status_;
  // L line
  std::set<std::string> protein_ids_;
  
 public:
  /**
   * \returns a blank SqtFileReader object
   */
  SqtFileReader();
  
  /**
   * \returns a SqtFileReader object and loads the sqt file
   * specified by the file_name
   */
  SqtFileReader(
    const char* file_name
  );


  SqtFileReader(
    const std::string& file_name
  );

  /**
   * Desctructor
   */
  
  ~SqtFileReader();

  /**
   * Loads the sqt file specified by the file_name
   */
  void loadData(
    const char* file_name
  );
  
  /**
   * Updates the variables with the next match on iteration
   * and does nothing if iterated through entire file
   */
  void next();

  /**
   * Resets the file pointer to the beginning of the file
   */
  void reset();

  /**
   * \Returns TRUE if there are more rows to iterate and
   * false otherwise
   */
  BOOLEAN_T hasNext();


  /**
   * Getters
   */
  int getLowScan();
  int getHighScan();
  int getCharge();
  double getProcessTime();
  std::string getServer();
  double getExperimentalMass();
  double getTotalIonIntensity();
  double getLowestSp();
  int getNumSeqMatched();
  int getXCorrRank();
  int getSpRank();
  double getCalculatedMass();
  double getDeltaCn();
  double getXCorr();
  double getSp();
  int getMatchedIons();
  int getExpectedIons();
  std::string getSequenceMatched();
  std::string getValidationStatus();
  std::set<std::string>getProteinIds();
  BOOLEAN_T okay();
  
 private:
  void getLLine();
  void getSLine();
  void getMLine();
};



#endif //SQTFILEREADER_H
