#ifndef CROSSLINKEDPEPTIDE_H_
#define CROSSLINKEDPEPTIDE_H_

class XLinkablePeptide {
 protected:
  PEPTIDE_T* peptide;
  vector<int> link_sites;
};


class XLinkPeptide {
 protected:
  static FLOAT_T linker_mass_;
  std::vector<XLinkablePeptide> linked_peptides_;
  std::vector<int> link_idx;

  bool mass_calculated_[NUMBER_MASS_TYPES];
  FLOAT_T mass_[NUMBER_MASS_TYPES];

  BOOLEAN_T is_decoy_;

 public:
};

class XLinkMatch {
 protected:
  CrossLinkedPeptide peptide_;
  FLOAT_T xcorr_score_;
  FLOAT_T ppm_error_;
  FLOAT_T pvalue_;
  
  

};




#endif



