#ifndef XLINK_H_
#define XLINK_H_
#include "objects.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class MatchCandidate;
class MatchCandidateVector;


void get_min_max_mass(
  FLOAT_T precursor_mz, 
  int charge, 
  BOOLEAN_T use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass);


  
#endif