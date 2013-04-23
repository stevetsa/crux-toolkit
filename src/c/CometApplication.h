/**
 * \file CruxHardklorApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef COMETAPPLICATION_H
#define COMETAPPLICATION_H

#include "CruxApplication.h"

#include <string>
#include <fstream>
#include <iostream>

class CometApplication: public CruxApplication {

 protected:

 public:

  /**
   * \returns a blank CometApplication object
   */
  CometApplication();

  /**
   * Destructor
   */
  ~CometApplication();

  /**
   * main method for CometApplication
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for CometApplication
   */
  virtual std::string getName();

  /**
   * \returns the description for CometApplication
   */
  virtual std::string getDescription();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();

  /**
   * \write parameters 
   */
  void writeParams(std::ofstream &fout, 
    std::string protein_database
  );
  
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
