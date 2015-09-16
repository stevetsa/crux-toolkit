#include "XLinkDatabase.h"
#include "util/modifications.h"
#include "model/ModifiedPeptidesIterator.h"
#include "util/GlobalParams.h"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

Database* XLinkDatabase::protein_database_;
XLinkBondMap XLinkDatabase::bondmap_;

vector<vector<Crux::Peptide*> > XLinkDatabase::target_peptides_; ///< all peptides generated with no additional missed cleavage;

std::vector<Crux::Peptide*> XLinkDatabase::decoy_peptides0_;
std::vector<Crux::Peptide*> XLinkDatabase::decoy_peptides1_;
std::vector<Crux::Peptide*> XLinkDatabase::decoy_peptides2_;

std::vector<LinearPeptide> XLinkDatabase::target_linear_peptides_;
std::vector<LinearPeptide> XLinkDatabase::decoy_linear_peptides_;

std::vector<SelfLoopPeptide> XLinkDatabase::target_selfloop_peptides_;
std::vector<SelfLoopPeptide> XLinkDatabase::decoy_selfloop_peptides_;

std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides_;
std::vector<XLinkablePeptide> XLinkDatabase::decoy_xlinkable_peptides_;
std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides2_; //Peptides that could be selfloops.

std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides_flatten_;



void XLinkDatabase::initialize() {
  //Step one, load the database
  string input_file = get_string_parameter("protein fasta file");
  string link_string = get_string_parameter("link sites");
  bondmap_ = XLinkBondMap(link_string);
  protein_database_ = NULL;
  int num_protein = prepare_protein_input(input_file, &protein_database_);
  
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  //Step two create all peptides.
  int additional_cleavages = 1;

  if (get_boolean_parameter("xlink-include-selfloops")) {
    additional_cleavages = 2;
  }

  //Generate all possible target peptides
  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    double delta_mass = peptide_mod_get_mass_change(peptide_mod);
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
      int missed_cleavages = peptide->getMissedCleavageSites();
      while(target_peptides_.size() <= missed_cleavages) {
	target_peptides_.push_back(vector<Crux::Peptide*>());
      }
      target_peptides_[missed_cleavages].push_back(peptide);
    }
    delete peptide_iterator;
  }

  

  //Sort by mass
  //sort(target_peptides0_.begin(), target_peptides0_.end(), comparePeptideMass); 
  //sort(target_peptides1_.begin(), target_peptides1_.end(), comparePeptideMass);
  //sort(target_peptides2_.begin(), target_peptides2_.end(), comparePeptideMass);

  // If linear peptides were asked for, then generate them from the target_peptides0 list.
  if (get_boolean_parameter("xlink-include-linears")) {
    generateAllLinears(false);
    sort(target_linear_peptides_.begin(), target_linear_peptides_.end(), compareLinearPeptideMass);
    carp(CARP_INFO, "There are %d linear peptides", target_linear_peptides_.size());
  }

  //Generate all linkable peptides, allowing for suppressed cleavages
  int max_cleavages = GlobalParams::getMissedCleavages()+1;
  for (size_t cleavage_idx=0;cleavage_idx<=max_cleavages;cleavage_idx++) {
    generateAllLinkablePeptides(target_peptides_[cleavage_idx], target_xlinkable_peptides_);
  }
  sort(target_xlinkable_peptides_.begin(), target_xlinkable_peptides_.end(), compareXLinkablePeptideMass);
  carp(CARP_INFO, "There are %d xlinkable peptides", target_xlinkable_peptides_.size());

  flattenLinkablePeptides(target_xlinkable_peptides_, target_xlinkable_peptides_flatten_);

  //TODO, generate all mono/dead link peptides

  if (get_boolean_parameter("xlink-include-selfloops")) {
    generateAllSelfLoops(true);
    carp(CARP_INFO, "There are %d self loop peptides", target_selfloop_peptides_.size());
  }
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

  findSelfLoops(target_xlinkable_peptides_, target_selfloop_peptides_);
  
  
  int mc = GlobalParams::getMissedCleavages()+2;
  generateAllLinkablePeptides(target_peptides_[mc], target_xlinkable_peptides2_);
  findSelfLoops(target_xlinkable_peptides2_, target_selfloop_peptides_);

  sort(target_selfloop_peptides_.begin(), target_selfloop_peptides_.end(), compareSelfLoopPeptideMass);

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

  if (!decoy) {

    for (size_t cleavage_idx = 0; cleavage_idx <= GlobalParams::getMissedCleavages(); cleavage_idx++) {
      findLinearPeptides(target_peptides_[cleavage_idx], target_linear_peptides_);
    }
  } else {
    carp(CARP_FATAL, "XLinkDatabase::generateAllLinears(true) not implemented yet!");
  }
}


void XLinkDatabase::flattenLinkablePeptides(vector<XLinkablePeptide>& xpeptides,
					    vector<XLinkablePeptide>& flattened) {

  for (size_t idx=0;idx < xpeptides.size();idx++) {
    XLinkablePeptide& current = xpeptides[idx];
    for (size_t link1_idx=0;link1_idx < current.numLinkSites(); link1_idx++) {
      XLinkablePeptide onelink(current);
      onelink.clearSites();
      onelink.addLinkSite(current.getLinkSite(link1_idx));
      flattened.push_back(onelink);
    }
  }

}

void XLinkDatabase::generateAllLinkablePeptides(vector<Crux::Peptide*>& peptides, 
  vector<XLinkablePeptide>& xpeptides) {
  //Loop through peptides
  vector<int> link_sites;
  for (vector<Crux::Peptide*>::iterator iter = peptides.begin(); 
    iter != peptides.end(); 
       ++iter) {
    Crux::Peptide* peptide = *iter;
    if (peptide->countModifiedAAs() <= GlobalParams::getMaxXLinkMods()) {
      XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites); 
      if (!link_sites.empty()) {
        xpeptides.push_back(XLinkablePeptide(peptide, link_sites));
      }
    }
  }
}


void XLinkDatabase::print() {
   
    string output_directory = get_string_parameter("output-dir");
    if (get_boolean_parameter("xlink-include-linears")) {
      ostringstream oss;
      oss << output_directory << "/" << "xlink_peptides.linear.txt";
      string temp = oss.str();
      ofstream peptides_file(temp.c_str());

      peptides_file << setprecision(8);

      peptides_file << "mass\tsequence\tprotein id"<<endl;

      for (int idx=0;idx < target_linear_peptides_.size();idx++) {
     
        peptides_file << target_linear_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << target_linear_peptides_[idx].getSequenceString() << "\t";
        peptides_file << target_linear_peptides_[idx].getProteinIdString() << endl;
      }
      peptides_file.flush();
    }

    if (get_boolean_parameter("xlink-include-selfloops")) {
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

vector<LinearPeptide>::iterator XLinkDatabase::getLinearBegin() {

  return (target_linear_peptides_.begin());

}
vector<LinearPeptide>::iterator XLinkDatabase::getLinearBegin(FLOAT_T min_mass) {
  return (lower_bound(target_linear_peptides_.begin(), target_linear_peptides_.end(), min_mass, compareLinearPeptideMassToFLOAT)); 
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearEnd() {
  return (target_linear_peptides_.end());
}

vector<SelfLoopPeptide>::iterator XLinkDatabase::getSelfLoopBegin(FLOAT_T min_mass) {
  return (lower_bound(target_selfloop_peptides_.begin(), 
                      target_selfloop_peptides_.end(), 
                      min_mass, 
                      compareSelfLoopPeptideMassToFLOAT));
}

vector<SelfLoopPeptide>::iterator XLinkDatabase::getSelfLoopEnd() {
  return(target_selfloop_peptides_.end());
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

XLinkBondMap& XLinkDatabase::getXLinkBondMap() {
  return bondmap_;
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenBegin() {
  return(target_xlinkable_peptides_flatten_.begin());
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenBegin(FLOAT_T min_mass) {
  return(lower_bound(target_xlinkable_peptides_flatten_.begin(), target_xlinkable_peptides_flatten_.end(),
		     min_mass, compareXLinkablePeptideMassToFLOAT));
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenEnd() {
  return(target_xlinkable_peptides_flatten_.end());
}
