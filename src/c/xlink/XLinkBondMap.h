#ifndef XLINKBONDMAP_H_
#define XLINKBONDMAP_H_

#include <map>
#include <set>
#include <string>

#include "XLinkablePeptide.h"

class XLinkBondMap: public std::map<char, std::set<char> > {
 public:
  XLinkBondMap();
  virtual ~XLinkBondMap();
  XLinkBondMap(std::string& links_string);

  void setLinkString(std::string& links_string);

  BOOLEAN_T canLink(char aa);

  BOOLEAN_T canLink(XLinkablePeptide& pep1,
		    XLinkablePeptide& pep2,
		    int link1_site,
		    int link2_site);

  BOOLEAN_T canLink(XLinkablePeptide& pep,
		    int link1_site,
		    int link2_site);

  BOOLEAN_T canLinkIdx(XLinkablePeptide& pep1,
		       XLinkablePeptide& pep2,
		       int link1_site_idx,
		       int link2_site_idx);

};

#endif

