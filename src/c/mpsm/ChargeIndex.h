#ifndef _CHARGE_INDEX_H_
#define _CHARGE_INDEX_H_

#include <string>
#include <vector>

class ChargeIndex: public std::vector<int> {
  
public:

  ChargeIndex();
  ChargeIndex(int charge);
  ChargeIndex(std::string& charges, char delimiter=',');

  void add(int charge);

  int max();

  int numCharge(int charge);


  bool operator ==(ChargeIndex& c); 
  
  bool operator <(ChargeIndex& c);

  bool operator > (ChargeIndex& c);

  friend std::ostream& operator<<(std::ostream &os, ChargeIndex &obj);
  friend std::istream& operator>>(std::istream &is, ChargeIndex &obj);


};

#endif
