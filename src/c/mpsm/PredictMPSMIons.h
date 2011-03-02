/**
 * \file PredictMPSMIons.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running search-for-xlinks
 *****************************************************************************/

#ifndef PREDICTMPSMIONS_H
#define PREDICTMPSMIONS_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class PredictMPSMIons: public CruxApplication {

 public:

  PredictMPSMIons();
  ~PredictMPSMIons();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
