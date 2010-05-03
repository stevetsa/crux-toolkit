/**
 * \file MPSM_OutputFiles.h
 */
/*
 * FILE: output-files.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: Aug 24, 2009
 * PROJECT: crux
 * DESCRIPTION: A class description for handling all the various
 * output files, excluding parameter and log files.  The filenames,
 * locations and overwrite status would be taken from parameter.c.
 * adds additional handling for mpsm matches
 */

#ifndef MPSM_OUTPUTFILES_H
#define MPSM_OUTPUTFILES_H

#include "output-files.h"

#include "MPSM_ChargeMap.h"
#include "MPSM_MatchCollection.h"
#include "MPSM_Match.h"

class MPSM_OutputFiles: public OutputFiles {

  public:
    
    MPSM_OutputFiles(COMMAND_T program_name): OutputFiles(program_name){};
    


    void writeMatches(MPSM_ChargeMap& charge_map);
  protected:
    void writeMatches(std::vector<MPSM_MatchCollection>& mpsm_match_collections);
    void writeMatches(FILE* file_ptr, MPSM_MatchCollection& mpsm_match_collection);
    void writeMatch(FILE* file_ptr, MPSM_Match& mpsm_match);
};

#endif //MPSM_OUTPUTFILES_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
