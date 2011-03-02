#include <algorithm>
#include <iostream>

#include "ZStateIndex.h"

using namespace std;

static bool use_unique_windows = true;


ZStateIndex::ZStateIndex() {
  ;
}

ZStateIndex::ZStateIndex(SpectrumZState& zstate) {
  add(zstate);
}

bool ZStateIndex::add(SpectrumZState& zstate) {
/*
  if (use_unique_windows) {
    for (int idx=0;idx < size();idx++) {
      if (at(idx) == zstate) {
        return false;
      }
    }
  }
*/
  push_back(zstate);

  //sort(begin(), end());

  return true;

  

}

bool ZStateIndex::operator == (ZStateIndex& z) {
  return false;
}

bool ZStateIndex::operator != (ZStateIndex& z) {

  return !(*this == z);
}

bool ZStateIndex::operator < (ZStateIndex& z) {
  return false;
}

bool ZStateIndex::operator > (ZStateIndex& z) {
  return false;
}


std::ostream& operator<<(std::ostream &os, ZStateIndex &obj) {
  return os;
}

std::istream& operator>>(std::istream &is, ZStateIndex &obj) {
  return is;
}
