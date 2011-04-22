/**
 * \file PredictMPSMIons.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running search-for-xlinks
 *****************************************************************************/

#ifndef ASSIGNIONSMPSMSPECTRUM_H
#define ASSIGNIONSMPSMSPECTRUM_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class AssignIonsMPSMSpectrum: public CruxApplication {

 public:

  AssignIonsMPSMSpectrum();
  ~AssignIonsMPSMSpectrum();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
