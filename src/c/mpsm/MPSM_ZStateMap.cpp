#include "MPSM_ZStateMap.h"

#include <iostream>

using namespace std;


void MPSM_ZStateMap::insert(ZStateIndex& zstate_index) {
  MPSM_ZStateMap:: iterator cur_iter;
  
  cur_iter = find(zstate_index);

  if (cur_iter == end()) {
    vector<MPSM_MatchCollection> match_collections;
    (*this)[zstate_index] = match_collections;
  }

  vector<MPSM_MatchCollection>& match_collections = at(zstate_index);

  while (match_collections.size() < 3) {
    MPSM_MatchCollection new_collection;
    match_collections.push_back(new_collection);
  }
  

}

//assume that this is one charge.
void MPSM_ZStateMap::insert(vector<MPSM_MatchCollection>& match_collection) {


  MPSM_ZStateMap::iterator cur_iter;

  ZStateIndex& zstate_index = match_collection[0][0].getZStateIndex();
  cur_iter = find(zstate_index);

  if (cur_iter == end()) {

    (*this)[zstate_index] = match_collection;

  } else {

    for (int collection_idx = 0;
      collection_idx < match_collection.size();
      collection_idx++) {
      MPSM_MatchCollection& this_collection = (*this)[zstate_index][collection_idx];
      MPSM_MatchCollection& new_collection = match_collection[collection_idx];
      for (int match_idx=0;
        match_idx < new_collection.numMatches();
        match_idx++) {
          this_collection.addMatch(new_collection[match_idx]);
      }
        
    }
  }
}


void MPSM_ZStateMap::insert(MPSM_MatchCollection& match_collection, int match_collection_idx) {
  MPSM_ZStateMap::iterator cur_iter;
  
  ZStateIndex& zstate_index = match_collection[0].getZStateIndex();
  cur_iter = find(zstate_index);

  if (cur_iter == end()) {
    vector<MPSM_MatchCollection> new_match_collections;
    (*this)[zstate_index] = new_match_collections;
    cur_iter = find(zstate_index);
  }

  if (cur_iter -> second.size() <= match_collection_idx) {
    MPSM_MatchCollection new_collection;
    cur_iter -> second.push_back(new_collection);
  }

  MPSM_MatchCollection& new_collection = cur_iter -> second[match_collection_idx];

  for (int match_idx;match_idx < match_collection.numMatches();match_idx++) {
    new_collection.addMatch(match_collection[match_idx]);
  }

}


void MPSM_ZStateMap::insert(MPSM_Match& new_match, int match_collection_idx) {
  ZStateIndex& zstate_index = new_match.getZStateIndex();
  
  MPSM_ZStateMap::iterator cur_iter = find(zstate_index);

  if (cur_iter == end()) {
    vector<MPSM_MatchCollection> new_match_collections;
    (*this)[zstate_index] = new_match_collections;
    cur_iter = find(zstate_index);
  }

  if (cur_iter -> second.size() <= match_collection_idx) {
    MPSM_MatchCollection new_collection;
    cur_iter -> second.push_back(new_collection);
  }

  cur_iter -> second[match_collection_idx].addMatch(new_match);

}

void MPSM_ZStateMap::insert(
  MPSM_ZStateMap& mpsm_chargemap
  ) {

  MPSM_ZStateMap::iterator new_iter;

  for (new_iter = mpsm_chargemap.begin();
    new_iter != mpsm_chargemap.end();
    ++new_iter) {

    insert(new_iter -> second);
  }
}


void MPSM_ZStateMap::sortMatches(SCORER_TYPE_T sort_type) {
  for (MPSM_ZStateMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;
    ZStateIndex current_charge = iter -> first;
    for (int collection_idx = 0; 
      collection_idx < match_collections.size();
      collection_idx++) {
      //cout <<"Sorting "<<current_charge<<" "<<collection_idx<<endl;
      match_collections[collection_idx].sortByScore(sort_type);
    }
     

  }

}

void MPSM_ZStateMap::clearMap() {
    for (MPSM_ZStateMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;
    ZStateIndex current_zstate = iter -> first;
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


void MPSM_ZStateMap::calcDeltaCN() {
  for (MPSM_ZStateMap::iterator iter = begin();
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

void MPSM_ZStateMap::calcZScores() {

  for (MPSM_ZStateMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;

    double mean,std;
    match_collections.back().calcZParameters(mean, std);

    for (int collection_idx = 0;
      collection_idx < match_collections.size();
      collection_idx++) {

      match_collections[collection_idx].setZParameters(mean, std);
    }
  }
}

void MPSM_ZStateMap::calcXCorrRanks() {

  for (MPSM_ZStateMap::iterator iter = begin();
    iter != end();
    ++iter) {

    vector<MPSM_MatchCollection>& match_collections = iter -> second;
    for (int collection_idx = 0;
      collection_idx < match_collections.size();
      collection_idx++) {

      match_collections[collection_idx].calcXCorrRanks();
    }
  }
}


MPSM_ZStateMap::~MPSM_ZStateMap() {
  clearMap();
}
