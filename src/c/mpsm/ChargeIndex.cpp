#include <algorithm>
#include <iostream>

#include "ChargeIndex.h"
#include "DelimitedFile.h"

using namespace std;

ChargeIndex::ChargeIndex() {
  ;
}

ChargeIndex::ChargeIndex(int charge) {
  push_back(charge);
}

ChargeIndex::ChargeIndex(string& charges, char delimiter) {
  vector<string> vcharges;

  DelimitedFile::tokenize(charges, vcharges, delimiter);

  for (int i=0;i<vcharges.size();i++) {
    int charge;
    DelimitedFile::from_string<int>(charge, vcharges[i]);
    add(charge);
  }


}


void ChargeIndex::add(int charge) {
  push_back(charge);
  sort(begin(), end());
}

int ChargeIndex::max() {
  return back();
}

int ChargeIndex::numCharge(int charge) {
  int ans = 0;

  for (int idx = 0; idx < size(); idx++) {
    if (at(idx) == charge) {
      ans++;
    }
  }
  return ans;
}

bool ChargeIndex::isHomogeneous() {

  if (size() == 1) return true;

  //since the charge mixture is sorted, we can just check the first and last.
  //elements if they are equal, then it is homogeneous, otherwise, it is not.
  return at(0) == at(size()-1);

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
