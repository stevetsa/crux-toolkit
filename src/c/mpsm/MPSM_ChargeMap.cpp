#include "MPSM_ChargeMap.h"

#include <iostream>

using namespace std;

//assume that this is one charge.
void MPSM_ChargeMap::insert(vector<MPSM_MatchCollection>& match_collection) {


  MPSM_ChargeMap::iterator cur_iter;

  ChargeIndex& charge_index = match_collection[0][0].getChargeIndex();
  cur_iter = find(charge_index);

  if (cur_iter == end()) {

    (*this)[charge_index] = match_collection;

  } else {

    for (int collection_idx = 0;
      collection_idx < match_collection.size();
      collection_idx++) {
      MPSM_MatchCollection& this_collection = (*this)[charge_index][collection_idx];
      MPSM_MatchCollection& new_collection = match_collection[collection_idx];
      for (int match_idx=0;
        match_idx < new_collection.numMatches();
        match_idx++) {
          this_collection.addMatch(new_collection[match_idx]);
      }
        
    }
  }
}


void MPSM_ChargeMap::insert(MPSM_ChargeMap& mpsm_chargemap) {

  MPSM_ChargeMap::iterator new_iter;

  for (new_iter = mpsm_chargemap.begin();
    new_iter != mpsm_chargemap.end();
    ++new_iter) {

    insert(new_iter -> second);
  }
}


void MPSM_ChargeMap::sortMatches(SCORER_TYPE_T sort_type) {
  for (MPSM_ChargeMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;
    ChargeIndex current_charge = iter -> first;
    for (int collection_idx = 0; 
      collection_idx < match_collections.size();
      collection_idx++) {
      //cout <<"Sorting "<<current_charge<<" "<<collection_idx<<endl;
      match_collections[collection_idx].sortByScore(sort_type);
    }
     

  }

}

void MPSM_ChargeMap::clearMap() {
    for (MPSM_ChargeMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;
    ChargeIndex current_charge = iter -> first;
    /*
    for (int collection_idx = 0; 
      collection_idx < match_collections.size();
      collection_idx++) {
      match_collections[collection_idx].clear();
    } 
    */
    match_collections.clear();
  }
  this -> clear();

}


void MPSM_ChargeMap::calcDeltaCN() {
  for (MPSM_ChargeMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;

    for (int collection_idx = 0;
      collection_idx < match_collections.size();
      collection_idx++) {

      MPSM_MatchCollection& match_collection = match_collections[collection_idx];
      match_collection.calcDeltaCN();

    }


  }
}

void MPSM_ChargeMap::calcZScores() {

  for (MPSM_ChargeMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;

    for (int collection_idx = 0;
      collection_idx < match_collections.size();
      collection_idx++) {

      MPSM_MatchCollection& match_collection = match_collections[collection_idx];
      match_collection.calcZScores();

    }
  }
}


MPSM_ChargeMap::~MPSM_ChargeMap() {
  clearMap();
}
