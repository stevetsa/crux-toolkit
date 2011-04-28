/**
 * \file SearchForXLinks.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running search-for-xlinks
 *****************************************************************************/

#ifndef SEARCHFORXLINKSMODS_H
#define SEARCHFORXLINKSMODS_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class SearchForXLinksMods: public CruxApplication {

 public:

  SearchForXLinksMods();
  ~SearchForXLinksMods();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  virtual bool needsOutputDirectory();

};


#endif
