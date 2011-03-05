#include <algorithm>
#include <iostream>

#include <iomanip>
#include <ios>
#include <limits>
#include <sstream>

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

  if (use_unique_windows) {
    for (int idx=0;idx < size();idx++) {
      if (at(idx) == zstate) {
        return false;
      }
    }
  }

  push_back(zstate);

  sort(begin(), end());

  return true;

  

}

bool ZStateIndex::operator == (ZStateIndex& z) {

  //cerr <<"ZStateIndex::operator=="<<endl;
  if (size() != z.size()) {
    return false;
  }

  for (int idx=0;idx < size();idx++) {
    if (at(idx) != z.at(idx)) {
      return false;
    }
  }

  return true;
}

bool ZStateIndex::operator != (ZStateIndex& z) {
  //cerr <<"ZStateIndex::operator !="<<endl;
  return !(*this == z);
}

bool ZStateIndex::operator < (const ZStateIndex& z) const {
  //cerr <<"ZStateIndex::operator <"<<endl;
  if (size() != z.size()) {
    return size() < z.size();
  }

  for (int idx=0;idx<size();idx++) {

    if (at(idx) != z.at(idx)) {
      return at(idx) < z.at(idx);
    }
  }
  //so everything is equal, return false.
  return false;
  


}

bool ZStateIndex::operator > (ZStateIndex& z) {
  //cerr <<"ZStateIndex::operator >"<<endl;
  return false;
}

string ZStateIndex::getChargeString() {
  ostringstream oss;

  oss << at(0).getCharge();

  for (int idx=1;idx<size();idx++) {
    oss << "," << at(idx).getCharge();
  }

  return oss.str();
  
}


std::ostream& operator<<(std::ostream &os, ZStateIndex &obj) {

  for (int idx=0;idx<obj.size();idx++) {
    os << "(" << obj[idx].getCharge() << ","<< obj[idx].getNeutralMass() << ")"; 
  }

  return os;
}

std::istream& operator>>(std::istream &is, ZStateIndex &obj) {
  return is;
}

