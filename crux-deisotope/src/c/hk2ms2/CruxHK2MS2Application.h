/**
 * \file CruxBullseyeApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef CRUXHK2MS2APPLICATION_H
#define CRUXHK2MS2APPLICATION_H

#include "CruxApplication.h"

#include <string>
#include <fstream>

class CruxHK2MS2Application: public CruxApplication {

 protected:

  //Calls the main method in bullseye
  int hk2ms2Main(int argc, char* argv[]);

 public:

  /**
   * \returns a blank CruxHK2MS2Application object
   */
  CruxHK2MS2Application();

  /**
   * Destructor
   */
  ~CruxHK2MS2Application();

  /**
   * main method for CruxHK2MS2Application
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for CruxHK2MS2Application
   */
  virtual std::string getName();

  /**
   * \returns the description for CruxHK2MS2Application
   */
  virtual std::string getDescription();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
