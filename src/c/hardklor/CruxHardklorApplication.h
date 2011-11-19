/**
 * \file CruxHardklorApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef CRUXHARDKLORAPPLICATION_H
#define CRUXHARDKLORAPPLICATION_H

#include "CruxApplication.h"

#include <string>
#include <fstream>

class CruxHardklorApplication: public CruxApplication {

 protected:
  std::string delimited_filename_;
  std::string column_name_string_;
  char delimiter_;
  bool header_;

  //Calls the main method in HardklorApp
  int hardklorMain(int argc, char* argv[]);

  /**
   * writes the ISOTOPE.DAT file for hardklor
   */
  void writeIsotopeDat(
    std::string& filename ///<path for dat file
  );

  /**
   * write the ISOTOPE.DAT to an output stream
   */
  void writeIsotopeDat(
    std::ostream& os ///< the output stream to use
  );

  /**
   * writes the Hardklor.dat to a path
   */
  void writeHardklorDat(
    std::string& filename ///<path to write the Hardklor.dat to
    );  

  /**
   * writes the Hardklor.dat to a stream
   */
  void writeHardklorDat(
    std::ostream& os ///< stream to write to.
    );  

  /**
   * \returns the filename 'output_directory'/'fileroot'.'name' 
   */
  std::string makeFileName(
    const std::string& output_directory, ///< the directory
    const std::string& fileroot, ///< the fileroot
    const std::string& name ///< the name of the file
    );

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
