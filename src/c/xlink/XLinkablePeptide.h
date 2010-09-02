#ifndef XLINKABLEPEPTIDE_H_
#define XLINKABLEPEPTIDE_H_

#include "objects.h"
#include "peptide.h"

#include <vector>
#include <string>

#include "XLinkBondMap.h"


class XLinkablePeptide {
 protected:
  PEPTIDE_T* peptide_;
  char* sequence_;
  std::vector<int> link_sites_;

 public:
  XLinkablePeptide();
  XLinkablePeptide(char* xlinkable);
  XLinkablePeptide(PEPTIDE_T* peptide, std::vector<int>& link_sites);
  XLinkablePeptide(PEPTIDE_T* peptide, XLinkBondMap& bondmap);
  
  virtual ~XLinkablePeptide();

  static void findLinkSites(
    PEPTIDE_T* peptide, 
    XLinkBondMap& bondmap, 
    std::vector<int>& link_sites
  );

  XLinkablePeptide shuffle();

  size_t numLinkSites();
  BOOLEAN_T isLinkable();
  int addLinkSite(int link_site);
  int getLinkSite(int link_site_idx);

  PEPTIDE_T* getPeptide();
  FLOAT_T getMass() const;
  FLOAT_T getMass(MASS_TYPE_T mass_type);
  char* getSequence();
  MODIFIED_AA_T* getModifiedSequence();
  //std::string getAASequence();
  std::string getModifiedSequenceString();
  
};

bool compareXLinkablePeptideMass(const XLinkablePeptide& xpep1, const XLinkablePeptide& xpep2);



#endif
