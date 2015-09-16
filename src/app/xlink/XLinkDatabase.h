#include "objects.h"
#include "model/Peptide.h"
#include "model/Database.h"
#include "XLinkablePeptide.h"
#include "SelfLoopPeptide.h"
#include "LinearPeptide.h"

#include <vector>

class XLinkDatabase {

 protected:
  static Database* protein_database_;
  static XLinkBondMap bondmap_;

  static std::vector<std::vector<Crux::Peptide*> > target_peptides_; ///< all peptides generated, indexed by missed cleavages

  static std::vector<Crux::Peptide*> decoy_peptides0_;
  static std::vector<Crux::Peptide*> decoy_peptides1_;
  static std::vector<Crux::Peptide*> decoy_peptides2_;

  static std::vector<LinearPeptide> target_linear_peptides_;
  static std::vector<LinearPeptide> decoy_linear_peptides_;

  static std::vector<SelfLoopPeptide> target_selfloop_peptides_;
  static std::vector<SelfLoopPeptide> decoy_selfloop_peptides_;

  static std::vector<XLinkablePeptide> target_xlinkable_peptides_;
  static std::vector<XLinkablePeptide> target_xlinkable_peptides_flatten_;
  static std::vector<XLinkablePeptide> decoy_xlinkable_peptides_;
  static std::vector<XLinkablePeptide> target_xlinkable_peptides2_; //Peptides that could be selfloops.


  static void findLinearPeptides(vector<Crux::Peptide*>& peptides, vector<LinearPeptide>& linears);
  static void generateAllLinears(bool decoy);
  static void generateAllLinkablePeptides(std::vector<Crux::Peptide*>& peptides, 
    std::vector<XLinkablePeptide>& xpeptides);
  static void findSelfLoops(
			    std::vector<XLinkablePeptide>& linkable_peptides,
			    std::vector<SelfLoopPeptide>& ans);

  static void generateAllSelfLoops(bool decoy);
  static void flattenLinkablePeptides(
   std::vector<XLinkablePeptide>& xpeptides,
   std::vector<XLinkablePeptide>& flattened
   );
  

 public:
  XLinkDatabase() {;}
  virtual ~XLinkDatabase() {;}

  static void initialize();
  static void finalize();
  
  static XLinkBondMap& getXLinkBondMap();

  static std::vector<XLinkablePeptide>::iterator getXLinkableBegin();
  static std::vector<XLinkablePeptide>::iterator getXLinkableBegin(FLOAT_T min_mass);
  static std::vector<XLinkablePeptide>::iterator getXLinkableEnd();

  static std::vector<XLinkablePeptide>::iterator getXLinkableFlattenBegin();
  static std::vector<XLinkablePeptide>::iterator getXLinkableFlattenBegin(FLOAT_T min_mass);
  static std::vector<XLinkablePeptide>::iterator getXLinkableFlattenEnd();

  static std::vector<SelfLoopPeptide>::iterator getSelfLoopBegin();
  static std::vector<SelfLoopPeptide>::iterator getSelfLoopEnd();
  static std::vector<SelfLoopPeptide>::iterator getSelfLoopBegin(FLOAT_T min_mass);

  static std::vector<LinearPeptide>::iterator getLinearBegin();
  static std::vector<LinearPeptide>::iterator getLinearBegin(FLOAT_T min_mass);
  static std::vector<LinearPeptide>::iterator getLinearEnd();



  static void print();

};
