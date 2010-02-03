#ifndef _CHARGE_INDEX_H_
#define _CHARGE_INDEX_H_

#include <vector>

class ChargeIndex: public std::vector<int> {
  
public:

  ChargeIndex();
  ChargeIndex(int charge);

  void add(int charge);

  bool operator ==(ChargeIndex& c); 
  
  bool operator <(ChargeIndex& c);

  bool operator > (ChargeIndex& c);

  friend std::ostream& operator<<(std::ostream &os, ChargeIndex &obj);
  friend std::istream& operator>>(std::istream &is, ChargeIndex &obj);


};

#endif
