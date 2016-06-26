#include "XLinkDatabase.h"
#include "util/modifications.h"
#include "model/ModifiedPeptidesIterator.h"
#include "util/GlobalParams.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

using namespace std;

Database* XLinkDatabase::protein_database_;
XLinkBondMap XLinkDatabase::bondmap_;

vector<vector<Crux::Peptide*> > XLinkDatabase::target_peptides_; ///< all peptides generated with no additional missed cleavage;

vector<vector<Crux::Peptide*> > XLinkDatabase::decoy_peptides_;

std::vector<LinearPeptide> XLinkDatabase::target_linear_peptides_;
std::vector<LinearPeptide> XLinkDatabase::decoy_linear_peptides_;

std::vector<SelfLoopPeptide> XLinkDatabase::target_selfloop_peptides_;
std::vector<SelfLoopPeptide> XLinkDatabase::decoy_selfloop_peptides_;

std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides_;
std::vector<XLinkablePeptide> XLinkDatabase::decoy_xlinkable_peptides_;
std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides2_; //Peptides that could be selfloops.
std::vector<XLinkablePeptide> XLinkDatabase::decoy_xlinkable_peptides2_;

std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides_flatten_;
std::vector<XLinkablePeptide> XLinkDatabase::decoy_xlinkable_peptides_flatten_;
//vector<pair<int, vector<XLinkablePeptide> > > XLinkDatabase::protein_idx_to_xpeptides_;


void XLinkDatabase::initialize() {
  carp(CARP_INFO, "Initializing database");
  //Step one, load the database
  string input_file = Params::GetString("protein fasta file");
  string link_string = Params::GetString("link sites");
  bondmap_ = XLinkBondMap(link_string);
  protein_database_ = NULL;
  int num_protein = prepare_protein_input(input_file, &protein_database_);
  
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  //Step two create all peptides.
  int additional_cleavages = 0;

  bool generate_xlinkable = Params::GetBool("xlink-include-inter-intra") ||
           Params::GetBool("xlink-include-inter") ||
           Params::GetBool("xlink-include-intra") ||
           Params::GetBool("xlink-include-selfloops") ||
           Params::GetBool("xlink-include-deadends");
  
  if (generate_xlinkable) {
    additional_cleavages = 1;
  }
  
  if (Params::GetBool("xlink-include-selfloops")) {
    additional_cleavages = 2;
  }

  int total_missed_cleavages = GlobalParams::getMissedCleavages()+additional_cleavages;

  for (size_t idx=0;idx<=total_missed_cleavages;idx++) {
    target_peptides_.push_back(vector<Crux::Peptide*>());
    decoy_peptides_.push_back(vector<Crux::Peptide*>());
  }
  //Generate all possible target peptides
  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    double delta_mass = peptide_mod_get_mass_change(peptide_mod);
    carp(CARP_INFO, "delta mass:%g", delta_mass);
    //
    ModifiedPeptidesIterator* peptide_iterator =
      new ModifiedPeptidesIterator(
        GlobalParams::getMinMass(), 
        GlobalParams::getMaxMass(), 
        peptide_mod, 
        false, 
        protein_database_,
        additional_cleavages);

    //add the targets
    while (peptide_iterator->hasNext()) {
      
      Crux::Peptide* peptide = peptide_iterator->next();
//      carp(CARP_INFO, "Current peptide:%s", peptide->getUnshuffledSequence().c_str());
      int missed_cleavages = peptide->getMissedCleavageSites();
      target_peptides_[missed_cleavages].push_back(peptide);
      //Crux::Peptide* decoy_peptide = new Crux::Peptide(peptide);
      //decoy_peptide->transformToDecoy();
      //peptide->addDecoyPeptide(decoy_peptide);
      //decoy_peptides_[missed_cleavages].push_back(decoy_peptide);
    }

      
    delete peptide_iterator;
  }
  //Sort by mass
  //sort(target_peptides0_.begin(), target_peptides0_.end(), comparePeptideMass); 
  //sort(target_peptides1_.begin(), target_peptides1_.end(), comparePeptideMass);
  //sort(target_peptides2_.begin(), target_peptides2_.end(), comparePeptideMass);

  // If linear peptides were asked for, then generate them from the target_peptides0 list.
  if (Params::GetBool("xlink-include-linears")) {
    generateAllLinears(false);
    sort(target_linear_peptides_.begin(), 
      target_linear_peptides_.end(), 
      compareLinearPeptideMass);
    generateAllLinears(true);
    sort(decoy_linear_peptides_.begin(), 
      decoy_linear_peptides_.end(), 
      compareLinearPeptideMass);
    carp(CARP_INFO, "There are %d target linear peptides", target_linear_peptides_.size());
    carp(CARP_INFO, "There are %d decoy linear peptides", decoy_linear_peptides_.size());
  }

  if (generate_xlinkable) {
  
    //Generate all linkable peptides, allowing for suppressed cleavages
    carp(CARP_DEBUG, "Generating all linkable peptides");
    int max_cleavages = GlobalParams::getMissedCleavages()+1;
    for (size_t cleavage_idx=0;cleavage_idx<=max_cleavages;cleavage_idx++) {
      carp(CARP_DEBUG, "Generating target linkable peptides with cleavage:%d", cleavage_idx);
      generateAllLinkablePeptides(target_peptides_[cleavage_idx], target_xlinkable_peptides_);
      //carp(CARP_DEBUG, "Generating decoy linkable peptides with cleavage:%d", cleavage_idx);
      //generateAllLinkablePeptides(decoy_peptides_[cleavage_idx], decoy_xlinkable_peptides_);
    }

    sort(
      target_xlinkable_peptides_.begin(), 
      target_xlinkable_peptides_.end(), 
      compareXLinkablePeptideMass);
    for (size_t idx = 0;idx < target_xlinkable_peptides_.size();idx++) {
      target_xlinkable_peptides_[idx].setIndex(idx);
    }
    sort(
      decoy_xlinkable_peptides_.begin(), 
      decoy_xlinkable_peptides_.end(), 
      compareXLinkablePeptideMass);
    for (size_t idx = 0;idx < decoy_xlinkable_peptides_.size();idx++) {
      decoy_xlinkable_peptides_[idx].setIndex(idx);
      decoy_xlinkable_peptides_[idx].setDecoy(true);
    }
    carp(CARP_INFO, "There are %d xlinkable target peptides", target_xlinkable_peptides_.size());
    //carp(CARP_INFO, "There are %d xlinkable decoy peptides", decoy_xlinkable_peptides_.size());



    flattenLinkablePeptides(target_xlinkable_peptides_, target_xlinkable_peptides_flatten_);
    //flattenLinkablePeptides(decoy_xlinkable_peptides_, decoy_xlinkable_peptides_flatten_);
  }
  //TODO, generate all mono/dead link peptides

  if (Params::GetBool("xlink-include-selfloops")) {
    generateAllSelfLoops(true);
    generateAllSelfLoops(false);
    carp(CARP_INFO, "There are %d target self loop peptides", target_selfloop_peptides_.size());
    carp(CARP_INFO, "There are %d decoy self loop peptides", decoy_selfloop_peptides_.size());
  }

  //filter linkable peptides based upon user parameters
  /*
  if (!Params::GetBool("xlink-include-inter-intra") ||
      (!Params::GetBool("xlink-include-inter") && !Params::GetBool("xlink-include-intra"))) {
    carp(CARP_INFO, "Filtering out linkable peptides");
    vector<XLinkablePeptide> filtered;
    filterLinkablePeptides(target_xlinkable_peptides_, filtered);
    target_xlinkable_peptides_=filtered;
  }
  */

  //Build the protein_idx_to_xpeptides_;
  /*
  vector<vector<XLinkablePeptide> > protein_idx_to_xpeptides;

  for (size_t idx=0;idx < target_xlinkable_peptides_.size();idx++) {
    XLinkablePeptide& pep1 = target_xlinkable_peptides_[idx];
    pep1.setIndex(idx);
    carp(CARP_INFO, "%d %d", idx, pep1.getIndex());
    for (PeptideSrcIterator src_iterator1 = pep1.getPeptide()->getPeptideSrcBegin();
      src_iterator1 != pep1.getPeptide()->getPeptideSrcEnd();
      ++src_iterator1) {
      PeptideSrc* src1 = *src_iterator1;
      size_t id1 = src1->getParentProtein()->getProteinIdx();
      while(protein_idx_to_xpeptides.size() <= id1) {
        protein_idx_to_xpeptides.push_back(vector<XLinkablePeptide>());
      }
      protein_idx_to_xpeptides.at(id1).push_back(pep1);
      carp(CARP_INFO, "%d %d %d", idx, id1, protein_idx_to_xpeptides[id1].back().getIndex());
    }
   
  }
  //carp(CARP_FATAL, "stop!");
  for (size_t idx = 0;idx < protein_idx_to_xpeptides.size();idx++) {
    if (protein_idx_to_xpeptides[idx].size() > 0) {
      protein_idx_to_xpeptides_.push_back(make_pair(idx, protein_idx_to_xpeptides[idx]));
    }
  }
  */
  

  carp(CARP_INFO, "Done initializing database");
}

void XLinkDatabase::finalize() {

  target_linear_peptides_.clear();
  decoy_linear_peptides_.clear();
  target_selfloop_peptides_.clear();
  decoy_selfloop_peptides_.clear();
  target_xlinkable_peptides_.clear();
  decoy_xlinkable_peptides_.clear();

  for (size_t idx1=0;idx1<target_peptides_.size();idx1++) {
    for (size_t idx2=0;idx2<target_peptides_[idx1].size();idx2++) {
      delete target_peptides_[idx1][idx2];
    }
  }

  for (size_t idx1=0;idx1<decoy_peptides_.size();idx1++) {
    for (size_t idx2=0;idx2<decoy_peptides_[idx1].size();idx2++) {
      delete decoy_peptides_[idx1][idx2];
    }
  }

  Database::freeDatabase(protein_database_);
}

void XLinkDatabase::findSelfLoops(
  vector<XLinkablePeptide>& linkable_peptides, 
  vector<SelfLoopPeptide>& ans) {

  //Loop through all peptides, trying to find peptides with at least two link sites.
  for (vector<XLinkablePeptide>::iterator iter =
	 linkable_peptides.begin();
       iter != linkable_peptides.end();
       ++iter) {

    if (iter->numLinkSites() > 1) {
      for (size_t link1_idx = 0; link1_idx<iter->numLinkSites()-1;link1_idx++) {
	for (size_t link2_idx = link1_idx+1; link2_idx < iter->numLinkSites();link2_idx++) {
	  if (bondmap_.canLink(*iter, link1_idx, link2_idx)) {
            //create the candidate.
	    SelfLoopPeptide self_loop(*iter, iter->getLinkSite(link1_idx), iter->getLinkSite(link2_idx));
	    if (self_loop.getNumMissedCleavages() <= GlobalParams::getMissedCleavages()) {
	      self_loop.getMass(GlobalParams::getIsotopicMass());
              ans.push_back(self_loop);
	    }
	  }
	}
      }
    }
  }
}

void XLinkDatabase::generateAllSelfLoops(bool decoy) {

  int mc = GlobalParams::getMissedCleavages()+2;
  if (decoy) {
    findSelfLoops(decoy_xlinkable_peptides_, decoy_selfloop_peptides_);
    int mc = GlobalParams::getMissedCleavages()+2;
    generateAllLinkablePeptides(decoy_peptides_[mc], decoy_xlinkable_peptides2_);
    findSelfLoops(decoy_xlinkable_peptides2_, decoy_selfloop_peptides_);
    sort(decoy_selfloop_peptides_.begin(), 
      decoy_selfloop_peptides_.end(), 
      compareSelfLoopPeptideMass);
  } else {
    findSelfLoops(target_xlinkable_peptides_, target_selfloop_peptides_);
    generateAllLinkablePeptides(target_peptides_[mc], target_xlinkable_peptides2_);
    findSelfLoops(target_xlinkable_peptides2_, target_selfloop_peptides_);
    sort(
      target_selfloop_peptides_.begin(), 
      target_selfloop_peptides_.end(), 
      compareSelfLoopPeptideMass);
  }
}

void XLinkDatabase::findLinearPeptides(vector<Crux::Peptide*> &peptides, vector<LinearPeptide>& linear_peptides) {

  for (vector<Crux::Peptide*>::iterator iter = peptides.begin();
       iter != peptides.end();
       ++iter) {
    LinearPeptide lpeptide(*iter);
    lpeptide.getMass(GlobalParams::getIsotopicMass());
    linear_peptides.push_back(lpeptide);
    
  }
}

void XLinkDatabase::generateAllLinears(bool decoy) {
  //NOTE - Right now, we need to call this with decoy=false first.
  
  for (size_t cleavage_idx = 0; 
         cleavage_idx <= GlobalParams::getMissedCleavages(); 
         cleavage_idx++) {
      findLinearPeptides(target_peptides_[cleavage_idx], target_linear_peptides_);
    }
}


void XLinkDatabase::flattenLinkablePeptides(vector<XLinkablePeptide>& xpeptides,
					    vector<XLinkablePeptide>& flattened) {

  for (size_t idx=0;idx < xpeptides.size();idx++) {
    XLinkablePeptide& current = xpeptides[idx];
    for (size_t link1_idx=0;link1_idx < current.numLinkSites(); link1_idx++) {
      XLinkablePeptide onelink(current);
      onelink.getMass(MONO);
      onelink.clearSites();
      onelink.addLinkSite(current.getLinkSite(link1_idx));
      flattened.push_back(onelink);
    }
  }
}

void XLinkDatabase::filterLinkablePeptides(
  vector<XLinkablePeptide>& xpeptides,
  vector<XLinkablePeptide>& filtered_xpeptides
  ) {
  bool filter1 = !Params::GetBool("xlink-include-inter-intra");
  bool filter2 = !Params::GetBool("xlink-include-inter") && !Params::GetBool("xlink-include-intra");
  
  for (size_t idx = 0 ;idx < xpeptides.size();idx++) {
  //quick check to see if peptide can come from multiple protein sources.
  //if true, then this should return XLINK_INTER_INTRA_CANDIDATE
    XLINKMATCH_TYPE_T ctype = 
      XLink::getCrossLinkCandidateType(xpeptides[idx].getPeptide(), 
				xpeptides[idx].getPeptide());
    
    if ((filter1 && ctype == XLINK_INTER_INTRA_CANDIDATE) ||
      (filter2 && ctype != XLINK_INTER_INTRA_CANDIDATE)) {
    } else {
      filtered_xpeptides.push_back(xpeptides[idx]);
    }
  }

  carp(CARP_INFO, "kept %d out of %d", filtered_xpeptides.size(), xpeptides.size());

}

void XLinkDatabase::generateAllLinkablePeptides(
  vector<Crux::Peptide*>& peptides, 
  vector<XLinkablePeptide>& xpeptides) {

  carp(CARP_DEBUG, "XLinkDatabase::generateAllLinkablePeptides: start()");
  carp(CARP_DEBUG, "Number of peptides:%d", peptides.size());
  //Loop through peptides
  vector<int> link_sites;
  for (vector<Crux::Peptide*>::iterator iter = peptides.begin(); 
    iter != peptides.end(); 
       ++iter) {
    Crux::Peptide* peptide = *iter;

    if (peptide->countModifiedAAs() <= GlobalParams::getMaxXLinkMods()) {
      XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites); 
      if (!link_sites.empty()) {
	XLinkablePeptide xlp(peptide, link_sites);
	xlp.getMass(MONO);
        xpeptides.push_back(xlp);
      }
    }
  }
  carp(CARP_DEBUG, "XLinkDatabase::generateAllLinkablePeptides: done()");
}


void XLinkDatabase::print() {
   
    string output_directory = Params::GetString("output-dir");
    if (Params::GetBool("xlink-include-linears")) {
      ostringstream oss;
      oss << output_directory << "/" << "xlink_peptides.linear.txt";
      string temp = oss.str();
      ofstream peptides_file(temp.c_str());

      peptides_file << setprecision(8);

      peptides_file << "mass\tsequence\tprotein id\tunshuffled"<<endl;

      for (int idx=0;idx < target_linear_peptides_.size();idx++) {
     
        peptides_file << target_linear_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << target_linear_peptides_[idx].getSequenceString() << "\t";
        peptides_file << target_linear_peptides_[idx].getProteinIdString() << "\t";
        peptides_file << "" << endl;
      }
      
      for (int idx=0;idx < decoy_linear_peptides_.size();idx++) {
        peptides_file << decoy_linear_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << decoy_linear_peptides_[idx].getSequenceString() << "\t";
        peptides_file << decoy_linear_peptides_[idx].getProteinIdString() << "\t";
        peptides_file << decoy_linear_peptides_[idx].getPeptide(0)->getUnshuffledSequence() << endl;
      }
      
      peptides_file.flush();
    }

    if (Params::GetBool("xlink-include-selfloops")) {
      ostringstream oss;
      oss << output_directory << "/" << "xlink_peptides.selfloops.txt";
      string temp = oss.str();
      ofstream peptides_file(temp.c_str());

      peptides_file << setprecision(8);

      peptides_file << "mass\tsequence\tprotein id"<<endl;

      for (int idx=0;idx < target_selfloop_peptides_.size();idx++) {
     
        peptides_file << target_selfloop_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << target_selfloop_peptides_[idx].getSequenceString() << "\t";
        peptides_file << target_selfloop_peptides_[idx].getProteinIdString() << endl;
      }
      peptides_file.flush();
    }      
    
    ostringstream oss;
    oss << output_directory << "/" << "xlink_peptides.linkable.txt";
     
    string temp = oss.str();
    ofstream peptides_file(temp.c_str());

    peptides_file << setprecision(8);

    peptides_file << "mass\tsequence"<<endl;

    for (int idx=0;idx < target_xlinkable_peptides_.size();idx++) {
     
      peptides_file << target_xlinkable_peptides_[idx].getMass(MONO) << "\t";
      peptides_file << target_xlinkable_peptides_[idx].getModifiedSequenceString() << "\t";
      //peptides_file << target_xlinkable_peptides_[idx].getProteinIdString() << endl;
      peptides_file << endl;
    }
    peptides_file.flush();
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearBegin(
  bool decoy
  ) {

  if (decoy) {
    return (decoy_linear_peptides_.begin());
  } else {
    return (target_linear_peptides_.begin());
  }
}
vector<LinearPeptide>::iterator XLinkDatabase::getLinearBegin(
  bool decoy,
  FLOAT_T min_mass
  ) {

  if (decoy) {
    return (lower_bound(decoy_linear_peptides_.begin(),
		        decoy_linear_peptides_.end(),
                        min_mass,
                        compareLinearPeptideMassToFLOAT));
  } else {
    return (lower_bound(target_linear_peptides_.begin(), 
                        target_linear_peptides_.end(), 
                        min_mass, 
                        compareLinearPeptideMassToFLOAT)); 
  }
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearEnd(
  bool decoy
) {

  if (decoy) {
    return (decoy_linear_peptides_.end());
  } else {
    return (target_linear_peptides_.end());
  }
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearEnd(
  bool decoy,
  FLOAT_T max_mass
  ) {
  return(std::upper_bound(getLinearBegin(decoy),
                          getLinearEnd(decoy),
                          max_mass,
                            compareLinearPeptideMassToFLOAT2));  

 
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearEnd(
  bool decoy,
  std::vector<LinearPeptide>::iterator& siter,
  FLOAT_T max_mass
) {
   return(std::upper_bound(siter,
                            getLinearEnd(decoy),
                            max_mass,
                            compareLinearPeptideMassToFLOAT2));  
  
}

vector<SelfLoopPeptide>::iterator XLinkDatabase::getSelfLoopBegin(
  bool decoy,
  FLOAT_T min_mass
  ) {

  if (decoy) {
    return (lower_bound(decoy_selfloop_peptides_.begin(), 
                        decoy_selfloop_peptides_.end(), 
                        min_mass, 
                        compareSelfLoopPeptideMassToFLOAT));
  } else {

    return (lower_bound(target_selfloop_peptides_.begin(), 
                        target_selfloop_peptides_.end(), 
                        min_mass, 
                        compareSelfLoopPeptideMassToFLOAT));
  }
}

vector<SelfLoopPeptide>::iterator XLinkDatabase::getSelfLoopEnd(
  bool decoy
) {

  if (decoy) {
    return(decoy_selfloop_peptides_.end());
  } else {
    return(target_selfloop_peptides_.end());
  }
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableBegin() {
  return(target_xlinkable_peptides_.begin());

}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableBegin(FLOAT_T min_mass) {
  return (lower_bound(target_xlinkable_peptides_.begin(), target_xlinkable_peptides_.end(),
		      min_mass, compareXLinkablePeptideMassToFLOAT));
}
vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableEnd() {
    return(target_xlinkable_peptides_.end());
}

vector<XLinkablePeptide>& XLinkDatabase::getXLinkablePeptides(
  bool decoy
) {
  if (decoy) {
    return (decoy_xlinkable_peptides_);
  } else {
    return(target_xlinkable_peptides_);
  }
}



XLinkBondMap& XLinkDatabase::getXLinkBondMap() {
  return bondmap_;
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenBegin() {
  return(target_xlinkable_peptides_flatten_.begin());
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenBegin(
  bool decoy,
  FLOAT_T min_mass
  ) {
  if (decoy) {
    return(lower_bound(decoy_xlinkable_peptides_flatten_.begin(), 
                       decoy_xlinkable_peptides_flatten_.end(),
	  	     min_mass, compareXLinkablePeptideMassToFLOAT));
  } else {
    return(lower_bound(target_xlinkable_peptides_flatten_.begin(), 
                       target_xlinkable_peptides_flatten_.end(),
	  	     min_mass, compareXLinkablePeptideMassToFLOAT));
  }
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenEnd() {
  return(target_xlinkable_peptides_flatten_.end());
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenEnd(
									 bool decoy,
  FLOAT_T max_mass
  ) {

  if (decoy) {
    return(std::upper_bound(decoy_xlinkable_peptides_flatten_.begin(),
	  	     decoy_xlinkable_peptides_flatten_.end(),
		     max_mass,
		     compareXLinkablePeptideMassToFLOAT2));
  } else {

    return(std::upper_bound(target_xlinkable_peptides_flatten_.begin(),
	  	     target_xlinkable_peptides_flatten_.end(),
		     max_mass,
		     compareXLinkablePeptideMassToFLOAT2));
  }
}

/*
vector<pair<int, vector<XLinkablePeptide> > >& XLinkDatabase::getTargetProteinIdxToXPeptides() {
  return protein_idx_to_xpeptides_;

}
*/

int XLinkDatabase::getNLinkable() {
  return(target_xlinkable_peptides_.size());
}
