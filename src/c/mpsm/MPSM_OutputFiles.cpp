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

#include "MPSM_OutputFiles.h"

#include <iomanip>
#include <ios>
#include <limits>


using namespace std;


void MPSM_OutputFiles::writeMatches(MPSM_ZStateMap& charge_map) {

  for (MPSM_ZStateMap::iterator iter = charge_map.begin();
    iter != charge_map.end();
    ++iter) {
    ZStateIndex charge = iter -> first;
    cout <<"Printing matches for "<<charge<<endl;

    vector<MPSM_MatchCollection>& mpsm_match_collections = iter -> second;

    //if (charge.size() == 1) {
      for (int collection_idx = 0; collection_idx < mpsm_match_collections.size(); collection_idx++) {
        MatchFileWriter* file_ptr = getFilePtr(collection_idx);
        writeMatches(file_ptr, mpsm_match_collections[collection_idx]);
      }
    //} else {
    //  writeMatches(mpsm_match_collections);
    //}
  }
}


void MPSM_OutputFiles::writeMatches(vector<MPSM_MatchCollection>& mpsm_match_collections) {
  int collection_idx;

  cout <<"Size of collection:"<<mpsm_match_collections.size()<<endl;

  for (collection_idx = 0; collection_idx < mpsm_match_collections.size(); collection_idx++) {
    int file_idx;

    if (collection_idx == 0) {
      file_idx = 0;
    } else if (collection_idx <= 2) {
      file_idx = 1;
    } else {
      file_idx = 2;
    }

    //cout <<"collection_idx:"<<collection_idx<<" file_idx:"<<file_idx<<endl;
    
    MatchFileWriter* file_ptr = getFilePtr(file_idx);
    writeMatches(file_ptr, mpsm_match_collections[collection_idx]);
    
  }
}

void MPSM_OutputFiles::writeMatches(
  MatchFileWriter* file_ptr, 
  MPSM_MatchCollection& mpsm_match_collection) {

  bool do_sort = get_boolean_parameter("mpsm-do-sort");

  int top_match = get_int_parameter("mpsm-top-match");

  for (int match_idx = 0;match_idx < mpsm_match_collection.numMatches(); match_idx++) {
    //cout <<"Match "<<match_idx<<" out of "<<n<<endl;

    MPSM_Match& current_match = mpsm_match_collection.getMatch(match_idx);
    current_match.setParent(&mpsm_match_collection);
    int rank = current_match.getXCorrRank();
    
    

    if (do_sort) {
      if (rank <= top_match) {
        writeMatch(file_ptr, current_match);
      } else {
        break;
      }
    } else {
      if (match_idx < top_match) {
        writeMatch(file_ptr, current_match);
      } else {
        break;
      }
    }
  }
}

void MPSM_OutputFiles::writeMatch(
  MatchFileWriter* file_ptr,
  MPSM_Match& mpsm_match) {


  file_ptr->setColumnCurrentRow(SCAN_COL, mpsm_match.getFirstScan());
  file_ptr->setColumnCurrentRow(CHARGE_COL,  mpsm_match.getChargeString());    
  file_ptr->setColumnCurrentRow(SPECTRUM_PRECURSOR_MZ_COL, mpsm_match.getSpectrumPrecursorMZ());
  file_ptr->setColumnCurrentRow(SPECTRUM_NEUTRAL_MASS_COL, mpsm_match.getNeutralMassString());
  file_ptr->setColumnCurrentRow(PEPTIDE_MASS_COL, mpsm_match.getPeptideMassString());
  file_ptr->setColumnCurrentRow(XCORR_SCORE_COL, mpsm_match.getScore(XCORR));
  file_ptr->setColumnCurrentRow(XCORR_RANK_COL, mpsm_match.getXCorrRank());
  file_ptr->setColumnCurrentRow(MATCHES_SPECTRUM_COL, mpsm_match.getMatchesPerSpectrum());
  file_ptr->setColumnCurrentRow(SEQUENCE_COL, mpsm_match.getSequenceString());
  file_ptr->setColumnCurrentRow(NZSTATE_COL, mpsm_match.getSpectrum()->getNumZStates());
  file_ptr->setColumnCurrentRow(RTIME_MAX_DIFF_COL, mpsm_match.getRTimeMaxDiff());
  file_ptr->writeRow();

}
