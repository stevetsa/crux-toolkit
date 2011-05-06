/**
 * \file StatColumn.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Give a tab delimited file and a comma-separated list of column names
 * print out a tab delimied file with only those columns
 *****************************************************************************/
#ifndef STATCOLUMN_H
#define STATCOLUMN_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class StatColumn: public CruxApplication {

 public:

  /**
   * \returns a blank ExtractRows object
   */
  StatColumn();

  /**
   * Destructor
   */
  ~StatColumn();

  /**
   * main method for StatColumn
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for StatColumn
   */
  virtual std::string getName();

  /**
   * \returns the description for StatColumn
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
