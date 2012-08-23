/**
 * \file crux-main.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 * \brief The starting point for the main crux program.
 *
 * Usage is "crux [command] [options] [arguments]" where command
 * is one of the primary crux commands.
 **/

#include "crux-main.h"
#include "crux-utils.h" // Need to get definition of NUM_FEATURES.

#include "CruxApplicationList.h"
#include "CreateIndex.h"
#include "MatchSearch.h"
#include "SequestSearch.h"
#include "ComputeQValues.h"
#include "ComputeQValuesLegacy.h"
#include "Percolator.h"
#include "QRanker.h"
#include "Barista.h"
#include "PrintProcessedSpectra.h"
#include "GeneratePeptides.h"
#include "GetMs2Spectrum.h"
#include "PredictPeptideIons.h"
#include "SearchForXLinks.h"
//#include "ExtractScanNeutralMass.h"
#include "ExtractColumns.h"
#include "FilterSpectraByFragments.h"
#include "SpectralCounts.h"
#include "ExtractRows.h"
#include "PrintVersion.h"
#include "StatColumn.h"
#include "SortColumn.h"
#include "CruxHardklorApplication.h"
#include "CruxBullseyeApplication.h"

/**
 * The starting point for crux.  Prints a general usage statement when
 * given no arguments.  Runs one of the crux commands, including
 * printing the current version number.
 */
int main(int argc, char** argv){

#ifdef _MSC_VER
  // Turn off auto-tranlation of line-feed to 
  // carriage-return/line-feed
  _set_fmode(_O_BINARY);
#endif 

  CruxApplicationList applications("crux");

  applications.add(new CreateIndex());

  // search
  applications.add(new MatchSearch());
  applications.add(new SequestSearch());
  applications.add(new SearchForXLinks());

  // post-search
  applications.add(new ComputeQValues());
  applications.add(new ComputeQValuesLegacy()); // depricated name
  applications.add(new Percolator());
  applications.add(new QRanker());
  applications.add(new Barista());
  applications.add(new SpectralCounts());
  // fasta/ms2 utilities
  applications.add(new PrintProcessedSpectra());
  applications.add(new GeneratePeptides());
  applications.add(new PredictPeptideIons());
  applications.add(new FilterSpectraByFragments());
  //applications.add(new ExtractScanNeutralMass());
  applications.add(new GetMs2Spectrum());

  // delimited file utilities
  applications.add(new ExtractColumns());
  applications.add(new ExtractRows());
  applications.add(new StatColumn());
  applications.add(new SortColumn());

  applications.add(new CruxHardklorApplication());
  applications.add(new CruxBullseyeApplication());
  applications.add(new PrintVersion());
  



  int ret = applications.main(argc, argv);
  close_log_file();
  free_parameters();
  return ret;

}// end main

