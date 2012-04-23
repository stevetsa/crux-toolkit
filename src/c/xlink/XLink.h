#ifndef XLINK_H_
#define XLINK_H_
#include "objects.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>
#include <string>


class XLinkMatch;
class XLinkMatchCollection;


void get_min_max_mass(
  FLOAT_T precursor_mz, 
  int charge, 
  bool use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass);



namespace XLink {
std::string get_protein_ids_locations(Peptide* peptide);
void addAllocatedPeptide(Peptide* peptide);
void deleteAllocatedPeptides();
};
  
#endif
