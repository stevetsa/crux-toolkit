/**
 * \file GetMs2Spectrum.h
 *
 * AUTHOR: Manijeh Naseri
 * CREATE DATE: February 2, 2012
 * DESCRIPTION: Main method for the generate-peptides command.
 *              Output all peptide sequences in the given fasta file
 *              that fall within all peptide constraints.
 */
#ifndef GETMS2SPECTRUM_H
#define GETMS2SPECTRUM_H

#include "CruxApplication.h"
#include "util/crux-utils.h"
#include "io/carp.h"
#include "parameter.h"

class GetMs2Spectrum: public CruxApplication {

 public:
  /**
   * \returns A blank GetMs2Spectrum object.
   */
  GetMs2Spectrum();
  
  /**
   * Destructor
   */
  ~GetMs2Spectrum();

  /**
   * Main method for GetMs2Spectrum.
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns The command name for GetMs2Spectrum.
   */
  virtual std::string getName() const;

  /**
   * \returns The description for GetMs2Spectrum.
   */
  virtual std::string getDescription() const;

  /**
   * \returns The command arguments
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * \returns The command options
   */
  virtual std::vector<std::string> getOptions() const;

  /**
   * \returns The command outputs
   */
  virtual std::map<std::string, std::string> getOutputs() const;

  /**
   * \returns The file stem of the application, default getName.
   */
 
  virtual COMMAND_T getCommand() const;

 
 protected:
  /**
   * Print header lines to stdout.
   */
  void printHeader();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

