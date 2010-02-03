#include <algorithm>
#include <iostream>

#include "ChargeIndex.h"

using namespace std;

ChargeIndex::ChargeIndex() {
  ;
}

ChargeIndex::ChargeIndex(int charge) {
  push_back(charge);
}

void ChargeIndex::add(int charge) {
  push_back(charge);
  sort(begin(), end());
}

bool ChargeIndex::operator==(ChargeIndex& c) {
  if (c.size() != size())
    return false;
  for (int i=0;i<size();i++) {
    if (c[i] != (*this)[i])
      return false;
  }
  return true;
}

bool ChargeIndex::operator < (ChargeIndex& c) {
  if (size() < c.size()) return true;
  if (size() > c.size()) return false;
  for (int i=0;i<size();i++) {
    if ((*this)[i] < c[i]) return true;
    if ((*this)[i] > c[i]) return false;
  }
  return false;
}

bool ChargeIndex::operator > (ChargeIndex& c) {
    if (size() > c.size()) return true;
    if (size() < c.size()) return false;
    for (int i=0;i<size();i++) {
      if ((*this)[i] > c[i]) return true;
      if ((*this)[i] < c[i]) return false; 
    }
    return false;
}

std::ostream& operator<<(std::ostream &os, ChargeIndex &obj) {
  if (obj.size() == 0) {os <<"Empty!";return os;}
  os << obj[0];
  for (int i=1;i<obj.size();i++) {
    os <<"," << obj[i];
  }
  return os;
}

std::istream& operator>>(std::istream &is, ChargeIndex &obj){
  // TODO.
}
