#ifndef SEQUESTSEARCH_H
#define SEQUESTSEARCH_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"
#include "OutputFiles.h"


#include <string>
#include <vector>

class SequestSearch: public CruxApplication {

 protected:
  // Private functions, commented below at definition
  void print_matches(
    OutputFiles& output_files,       
    MATCH_COLLECTION_T* target_psms, 
    std::vector<MATCH_COLLECTION_T*>& decoy_psms,
    Spectrum* spectrum,             
    BOOLEAN_T combine_target_decoy,
    int num_decoy_files
    );

 public:

  SequestSearch();
  ~SequestSearch();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  virtual std::string getFileString();
};


#endif
