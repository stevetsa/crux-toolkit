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


void MPSM_OutputFiles::writeMatches(MPSM_ChargeMap& charge_map) {

  for (MPSM_ChargeMap::iterator iter = charge_map.begin();
    iter != charge_map.end();
    ++iter) {
    ChargeIndex charge = iter -> first;
    cout <<"Printing matches for "<<charge<<endl;

    vector<MPSM_MatchCollection>& mpsm_match_collections = iter -> second;

    //if (charge.size() == 1) {
      for (int collection_idx = 0; collection_idx < mpsm_match_collections.size(); collection_idx++) {
        FILE* file_ptr = getFilePtr(collection_idx);
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

    cout <<"collection_idx:"<<collection_idx<<" file_idx:"<<file_idx<<endl;

    FILE* file_ptr = getFilePtr(file_idx);
    writeMatches(file_ptr, mpsm_match_collections[collection_idx]);
    
  }
}

void MPSM_OutputFiles::writeMatches(
  FILE* file_ptr, 
  MPSM_MatchCollection& mpsm_match_collection) {

  int match_idx;
  int n = min(mpsm_match_collection.numMatches(), getMatchesPerSpec());
  carp(CARP_INFO,"Printing %d matches", n);
  //TODO - fix this to rank after we have implemented the sort by XCORR.
  for (match_idx = 0;match_idx < n; match_idx++) {
    cout <<"Match "<<match_idx<<" out of "<<n<<endl;
    mpsm_match_collection.getMatch(match_idx).setParent(&mpsm_match_collection);
    writeMatch(file_ptr, mpsm_match_collection.getMatch(match_idx));
  }
}

void MPSM_OutputFiles::writeMatch(
  FILE* file_ptr,
  MPSM_Match& mpsm_match) {

  cout <<"writeMatch:start()"<<endl;
  ostringstream oss;

  oss << mpsm_match;

  std::string string_value = oss.str();

  fprintf(file_ptr, "%s\n", string_value.c_str());
  cout <<"writeMatch:end()"<<endl;
}
