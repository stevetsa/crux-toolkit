/**
 * \file StatColumn.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Given a delimited file and a column-name, print out statistics
 * for that column (n, min, max, sum, average, median).
 *****************************************************************************/
#ifndef CRUXHARDKLORAPPLICATION_H
#define CRUXHARDKLORAPPLICATION_H

#include "CruxApplication.h"

#include <string>

class CruxHardklorApplication: public CruxApplication {

 protected:
  std::string delimited_filename_;
  std::string column_name_string_;
  char delimiter_;
  bool header_;


 public:

  /**
   * \returns a blank CruxHardklorApplication object
   */
  CruxHardklorApplication();

  /**
   * Destructor
   */
  ~CruxHardklorApplication();

  /**
   * main method for CruxHardklorApplication
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for CruxHardklorApplication
   */
  virtual std::string getName();

  /**
   * \returns the description for CruxHardklorApplication
   */
  virtual std::string getDescription();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
