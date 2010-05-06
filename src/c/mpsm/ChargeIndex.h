#ifndef _CHARGE_INDEX_H_
#define _CHARGE_INDEX_H_

#include <string>
#include <vector>

class ChargeIndex: public std::vector<unsigned char> {
  
public:

  ChargeIndex();
  ChargeIndex(unsigned char charge);
  ChargeIndex(std::string& charges, char delimiter=',');

  void add(unsigned char charge);

  unsigned char max();

  int numCharge(unsigned char charge);

  bool isHomogeneous();



  bool operator ==(ChargeIndex& c); 
  
  bool operator <(ChargeIndex& c);

  bool operator > (ChargeIndex& c);

  friend std::ostream& operator<<(std::ostream &os, ChargeIndex &obj);
  friend std::istream& operator>>(std::istream &is, ChargeIndex &obj);


};

#endif
