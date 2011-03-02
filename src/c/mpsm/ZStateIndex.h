#ifndef ZSTATEINDEX_H
#define ZSTATEINDEX_H


#include <string>
#include <vector>

#include "SpectrumZState.h"


class ZStateIndex: public std::vector<SpectrumZState>{

 public:
  ZStateIndex();
  ZStateIndex(SpectrumZState& zstate);
  
  bool add(SpectrumZState& zstate);
  bool add(ZStateIndex zstates);


  int numZStates();
  
  bool operator == (ZStateIndex& z);
  bool operator != (ZStateIndex& z);
  bool operator < (ZStateIndex& z);
  bool operator > (ZStateIndex& z);

  friend std::ostream& operator<<(std::ostream &os, ZStateIndex &obj);
  friend std::istream& operator>>(std::istream &is, ZStateIndex &obj);

};

#endif
