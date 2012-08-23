#include "XLinkPeptide.h"
#include "ModifiedPeptidesIterator.h"
#include "IonSeries.h"
#include "Ion.h"

#include "XLinkablePeptideIterator.h"


#include <iostream>
#include <sstream>

using namespace std;

FLOAT_T XLinkPeptide::linker_mass_;
set<Peptide*> XLinkPeptide::allocated_peptides_;
FLOAT_T XLinkPeptide::pmin_;


XLinkPeptide::XLinkPeptide() : XLinkMatch() {

}

XLinkPeptide::XLinkPeptide(XLinkablePeptide& peptideA,
			   XLinkablePeptide& peptideB,
			   int posA, int posB) : XLinkMatch() {
  //cout<<"XLinkPeptide"<<endl;
  linked_peptides_.push_back(peptideA);
  linked_peptides_.push_back(peptideB);
  link_pos_idx_.push_back(posA);
  link_pos_idx_.push_back(posB);

  doSort();

}

XLinkPeptide::XLinkPeptide(char* peptideA,
			   char* peptideB,
			   int posA, int posB) : XLinkMatch() {

  //cout <<"Creating peptideA"<<endl;
  XLinkablePeptide A(peptideA);
  linked_peptides_.push_back(A);
  //cout <<"Creating peptideB"<<endl;
  XLinkablePeptide B(peptideB);
  linked_peptides_.push_back(B);
  //cout <<"Adding links"<<endl;
  A.addLinkSite(posA);
  link_pos_idx_.push_back(0);
  B.addLinkSite(posB);
  link_pos_idx_.push_back(0);

  doSort();
}


/**
 * makes sure that sequence1 is smaller in alphanumeric value than
 * sequence 2
 */
void XLinkPeptide::doSort() {

  string seq1 = linked_peptides_[0].getModifiedSequenceString();
  
  string seq2 = linked_peptides_[1].getModifiedSequenceString();


  if (seq1 > seq2) {

    //swap peptides
    swap(linked_peptides_[0], linked_peptides_[1]);
    //swap links
    swap(link_pos_idx_[0], link_pos_idx_[1]);
  }

  seq1 = linked_peptides_[0].getModifiedSequenceString();
  
  seq2 = linked_peptides_[1].getModifiedSequenceString();

  assert(seq1 <= seq2);



}

XLinkPeptide::~XLinkPeptide() {
}

void XLinkPeptide::setLinkerMass(FLOAT_T linker_mass) {
  linker_mass_=linker_mass;
}

FLOAT_T XLinkPeptide::getLinkerMass() {
  return linker_mass_;
}

int XLinkPeptide::getLinkPos(int peptide_idx) {
  return linked_peptides_[peptide_idx].getLinkSite(link_pos_idx_[peptide_idx]);
}


bool XLinkPeptide::isInter() {

  Peptide* peptide_a = linked_peptides_[0].getPeptide();
  Peptide* peptide_b = linked_peptides_[1].getPeptide();

  vector<Protein*> proteins_a;
  vector<Protein*> proteins_b;


  for (PeptideSrcIterator src_iterator = peptide_a->getPeptideSrcBegin();
    src_iterator != peptide_a->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_a.push_back(src->getParentProtein());
  }

  for (PeptideSrcIterator src_iterator = peptide_b->getPeptideSrcBegin();
    src_iterator != peptide_b->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_b.push_back(src->getParentProtein());
  }
  
  sort(proteins_a.begin(), proteins_a.end());
  sort(proteins_b.begin(), proteins_b.end());


  for (unsigned int idx_a=0;idx_a<proteins_a.size();idx_a++) {

    for (unsigned int idx_b=0;idx_b<proteins_b.size();idx_b++) {

      if (proteins_a[idx_a] != proteins_b[idx_b]) {
        //Found an instance where the proteins are not equal, therefore it is
        //inter. but still could be intra....
        return true;
      }
    }
  }

  //all proteins equal, so it is definately intra
  return false;

}

bool XLinkPeptide::isIntra() {

  Peptide* peptide_a = linked_peptides_[0].getPeptide();
  Peptide* peptide_b = linked_peptides_[1].getPeptide();

  vector<Protein*> proteins_a;
  vector<Protein*> proteins_b;
  

  for (PeptideSrcIterator src_iterator = peptide_a->getPeptideSrcBegin();
    src_iterator != peptide_a->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_a.push_back(src->getParentProtein());
  }

  for (PeptideSrcIterator src_iterator = peptide_b->getPeptideSrcBegin();
    src_iterator != peptide_b->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_b.push_back(src->getParentProtein());
  }

  sort(proteins_a.begin(), proteins_a.end());
  sort(proteins_b.begin(), proteins_b.end());


  for (unsigned int idx_a=0;idx_a<proteins_a.size();idx_a++) {

    for (unsigned int idx_b=0;idx_b<proteins_b.size();idx_b++) {

      if (proteins_a[idx_a] == proteins_b[idx_b]) {
        //Found an instance where the protein are equal, therefore it is
        //intra. but still also be inter as well....
        return true;
      }
    }
  }

  //all proteins inequal, so it is definately inter
  return false;




}

void XLinkPeptide::addLinkablePeptides(double min_mass, double max_mass,
			 Index* index, Database* database,
			 PEPTIDE_MOD_T* peptide_mod, bool is_decoy, 
			 XLinkBondMap& bondmap, 
			 vector<XLinkablePeptide>& linkable_peptides) {

  ModifiedPeptidesIterator* peptide_iterator =
    new ModifiedPeptidesIterator(
      min_mass, 
      max_mass,
      peptide_mod, 
      is_decoy,
      index, 
      database);

  int max_mod_xlink = get_int_parameter("max-xlink-mods");

  while (peptide_iterator->hasNext()) {
    Peptide* peptide = peptide_iterator->next();
    vector<int> link_sites;
    
    char* seq= peptide->getModifiedSequenceWithMasses(MOD_MASSES_SEPARATE);
    string proteins = peptide->getProteinIdsLocations();
    char* useq = peptide->getUnshuffledModifiedSequence(); 
    //cerr <<"Finding sites for :"<<seq<<":";

    //cout << seq << "\t" << proteins << "\t" << useq << "\t" << is_decoy << endl;

    free(seq);
    free(useq);
    if (peptide->countModifiedAAs() > max_mod_xlink) {
      delete peptide;
    } else {
      XLinkablePeptide::findLinkSites(peptide, bondmap, link_sites);

      if (link_sites.size() > 0) {
        XLinkablePeptide xlinkable_peptide(peptide, link_sites);
        xlinkable_peptide.setDecoy(is_decoy);
        linkable_peptides.push_back(xlinkable_peptide);
        XLink::addAllocatedPeptide(peptide);
      } else {
        delete peptide;
      }
    }
  }
  
  delete peptide_iterator;
}

void XLinkPeptide::addCandidates(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  XLinkBondMap& bondmap,
  Index* index, 
  Database* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  XLinkMatchCollection& candidates
  ) {

  //bool mixed_target_decoys = get_boolean_parameter("mixed-target-decoys");

  FLOAT_T min_length_mass = get_mass_amino_acid('G', MONO) * 
    (FLOAT_T)get_int_parameter("min-length") + 
    MASS_H2O_MONO;

  pmin_ = max((FLOAT_T)get_double_parameter("min-mass"), min_length_mass);

  FLOAT_T peptide1_min_mass = pmin_;
  FLOAT_T peptide1_max_mass = max_mass-pmin_-linker_mass_;

  for (int mod_idx1 = 0; mod_idx1 < num_peptide_mods; mod_idx1++) {
    PEPTIDE_MOD_T* peptide_mod1 = peptide_mods[mod_idx1];
    XLinkablePeptideIterator iter1_target(peptide1_min_mass, peptide1_max_mass, index, database, peptide_mod1, false, bondmap);
    XLinkablePeptideIterator iter1_decoy(peptide1_min_mass, peptide1_max_mass, index, database, peptide_mod1, true, bondmap);
    for (int mod_idx2 = 0; mod_idx2 < num_peptide_mods; mod_idx2++) {
      PEPTIDE_MOD_T* peptide_mod2 = peptide_mods[mod_idx2];
      addCandidates(min_mass, max_mass, bondmap, index, database, peptide_mod2, false, iter1_target, candidates);
      addCandidates(min_mass, max_mass, bondmap, index, database, peptide_mod2, true, iter1_decoy, candidates);
    }
  } 

  if (get_boolean_parameter("mixed-target-decoys")) {

    for (int mod_idx1 = 0; mod_idx1 < num_peptide_mods; mod_idx1++) {
      PEPTIDE_MOD_T* peptide_mod1 = peptide_mods[mod_idx1];
      XLinkablePeptideIterator iter1_target(peptide1_min_mass, peptide1_max_mass, index, database, peptide_mod1, false, bondmap);
      for (int mod_idx2 = 0; mod_idx2 < num_peptide_mods; mod_idx2++) {
        PEPTIDE_MOD_T* peptide_mod2 = peptide_mods[mod_idx2];
        addCandidates(min_mass, max_mass, bondmap, index, database, peptide_mod2, true, iter1_target, candidates);
      }
    } 
  }


}


void XLinkPeptide::addCandidates(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  XLinkBondMap& bondmap,
  Index* index,
  Database* database,
  PEPTIDE_MOD_T* peptide_mod2,
  bool decoy2,
  XLinkablePeptideIterator& iter1,
  XLinkMatchCollection& candidates) {

  int max_mod_xlink = get_int_parameter("max-xlink-mods");

  set<XLinkablePeptide> visited;

  while (iter1.hasNext()) {
    XLinkablePeptide pep1 = iter1.next();
    //cerr <<"pep1:"<<pep1.getSequence()<<endl;


    FLOAT_T peptide2_min_mass = min_mass - pep1.getMass() - linker_mass_;
    FLOAT_T peptide2_max_mass = max_mass - pep1.getMass() - linker_mass_;
    //cerr << "peptide2 min mass:"<<peptide2_min_mass << endl;
    //cerr << "peptide2 max mass:"<<peptide2_max_mass << endl;

    assert (peptide2_min_mass <= peptide2_max_mass);
  

    if (peptide2_max_mass >= pmin_) {
      XLinkablePeptideIterator iter2(peptide2_min_mass, peptide2_max_mass, index, database, peptide_mod2, decoy2, bondmap);
      while (iter2.hasNext()) {
        XLinkablePeptide pep2 = iter2.next();

        //cerr << "  pep2 is:"<<pep2.getSequence()<<endl;
        if (visited.find(pep2) == visited.end()) {
          //cerr << "pep1:"<<pep1.getSequence()<<" pep2:"<<pep2.getSequence()<<endl;

          FLOAT_T mass = pep1.getMass() + pep2.getMass() + linker_mass_;
  
          if ((mass >= min_mass) && (mass <= max_mass)) {
  
            int mods = pep1.getPeptide()->countModifiedAAs() + pep2.getPeptide()->countModifiedAAs();
            if (mods <= max_mod_xlink) {
              //carp(CARP_INFO,"Generating xlink candidates");
  
              //for every linkable site, generate the candidate if it is legal.
              for (unsigned int link1_idx=0;link1_idx < pep1.numLinkSites(); link1_idx++) {
                for (unsigned int link2_idx=0;link2_idx < pep2.numLinkSites();link2_idx++) {
                //cerr<<"link1_idx:"<<link1_idx<<endl;
                //cerr<<"link2_idx:"<<link2_idx<<endl;
                //cerr<<"Testing link:"<<endl;
                  if (bondmap.canLink(pep1, pep2, link1_idx, link2_idx)) {
                  //create the candidate
                    XLinkMatch* newCandidate = 
                      new XLinkPeptide(pep1, pep2, link1_idx, link2_idx);
                    candidates.add(newCandidate);
                  }
                }
              } // for link1_idx 
            } // if (mods <= max_mod_xlink .. 
          } // if ((mass >= min_mass) && (mass <= max_mass))
 
        } // if (visited.find(pep2) == visited.end()

      } // while (iter2.hasNext()) 

    } // if (peptide2_max_mass >= pmin_)

    visited.insert(pep1);
  } // while (iter1.hasNext()) 
}  



void XLinkPeptide::addCandidatesOld(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  XLinkBondMap& bondmap,
  Index* index, 
  Database* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  XLinkMatchCollection& candidates
  ) {

  bool mixed_target_decoys = get_boolean_parameter("mixed-target-decoys");
  carp(CARP_INFO, "mixed_target_decoys:%i", mixed_target_decoys);
  //get all linkable peptides up to mass-linkermass.

  //cerr <<"XLinkPeptide::addCandidates() : start."<<endl;
  //int max_missed_cleavages = get_int_parameter("missed-cleavages");
  vector<XLinkablePeptide> linkable_peptides;
  
  //iterate through each modification, 
  //get all peptides that are linkable up to the max_mass-linker_mass.
  // assess scores after all pmods with x amods have been searched
  int cur_aa_mods = 0;

  int max_mod_xlink = get_int_parameter("max-xlink-mods");

  // for each peptide mod
  for(int mod_idx=0; mod_idx<num_peptide_mods; mod_idx++){
    // get peptide mod
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];

    // is it time to assess matches?
    int this_aa_mods = peptide_mod_get_num_aa_mods(peptide_mod);
    
    if( this_aa_mods > cur_aa_mods ){
      carp(CARP_DEBUG, "Finished searching %i mods", cur_aa_mods);
      /*
	TODO - do we need this?
      BOOLEAN_T passes = is_search_complete(match_collection, cur_aa_mods);
      if( passes ){
        carp(CARP_DETAILED_DEBUG, 
             "Ending search with %i modifications per peptide", cur_aa_mods);
        break;
      }// else, search with more mods
      */
      cur_aa_mods = this_aa_mods;
    }
    //carp(CARP_INFO,"Calling addLinkablePeptides:%d",mod_idx);

    //carp(CARP_INFO,"max mass-linker_mass:%f", (max_mass-linker_mass_));

    


    addLinkablePeptides(0, max_mass-linker_mass_, index, database,
			peptide_mod, false, bondmap, linkable_peptides);

    addLinkablePeptides(0, max_mass-linker_mass_, index, database,
                        peptide_mod, true, bondmap, linkable_peptides);

     

    carp(CARP_INFO,"Done calling addLinkablePeptides:%d",mod_idx);
    
  }//next peptide mod
/*
  cout << "sequence\tprotein id\tunshuffed"<<endl; 
  for (size_t idx = 0; idx < linkable_peptides.size() ;idx++) {

    Peptide* peptide = linkable_peptides[idx].getPeptide();
    cout << peptide->getSequence()<< "\t" << peptide->getProteinIdsLocations() << "\t" << peptide->getUnshuffledSequence() << endl; 

  }  
*/  
  //exit(-1);

  if (linkable_peptides.size() == 0) {
    carp(CARP_WARNING, "No linkable peptides found!");
    return;
  }

  carp(CARP_INFO,"Sorting by mass: %d", linkable_peptides.size());
  //sort by increasing mass.
  sort(linkable_peptides.begin(),
       linkable_peptides.end(), 
       compareXLinkablePeptideMass);
  /*
  for (unsigned int idx=0;idx<linkable_peptides.size();idx++) {
    cerr <<linkable_peptides[idx].getMass()<<" "
         <<linkable_peptides[idx].getModifiedSequenceString()<<endl;
  }
  cerr << "======================================="<<endl;
  */
  unsigned int first_idx = 0;
  unsigned int next_idx = 0;
  
  while(first_idx < linkable_peptides.size()) {

    FLOAT_T first_mass = linkable_peptides.at(first_idx).getMass() + linker_mass_;

    if ((first_mass + linkable_peptides.at(first_idx).getMass()) > max_mass) {
      break; //we are done.
    }

    //cerr<<" mass:"<<first_mass<<endl;
    next_idx = first_idx;
    //cerr<<"last:"<<last_idx<<endl;
    FLOAT_T current_mass = first_mass + 
      linkable_peptides.at(next_idx).getMass();
    //cerr<<"current_mass:"<<current_mass<<endl;
    while ((next_idx < (linkable_peptides.size()-1)) && (current_mass < min_mass)) {
      next_idx++;
      //cerr << "new last:" << last_idx << endl;
      current_mass = first_mass + linkable_peptides.at(next_idx).getMass();
      //cerr << "new current:" << current_mass << endl;
    }

    while (next_idx < linkable_peptides.size() && current_mass <= max_mass) {
      //cerr<<"Adding links for peptides:"<<first_idx<<":"<<next_idx<<endl;
      XLinkablePeptide& pep1 = linkable_peptides.at(first_idx);
      XLinkablePeptide& pep2 = linkable_peptides.at(next_idx);
      //cerr << pep1.getModifiedSequenceString() <<" "<<pep1.isDecoy()<<endl;
      //cerr << pep2.getModifiedSequenceString() <<" "<<pep2.isDecoy()<<endl;
      //cerr << "======================="<<endl;
      if (mixed_target_decoys || (pep1.isDecoy() == pep2.isDecoy())) {
        FLOAT_T candidate_mass = pep1.getMass() + pep2.getMass() + linker_mass_;
        if (candidate_mass >= min_mass && candidate_mass <= max_mass) {
          int mods = pep1.getPeptide()->countModifiedAAs() + pep2.getPeptide()->countModifiedAAs();
          if (mods <= max_mod_xlink) {
            //carp(CARP_INFO,"Generating xlink candidates");

            //for every linkable site, generate the candidate if it is legal.
            for (unsigned int link1_idx=0;link1_idx < pep1.numLinkSites(); link1_idx++) {
              for (unsigned int link2_idx=0;link2_idx < pep2.numLinkSites();link2_idx++) {
              //cerr<<"link1_idx:"<<link1_idx<<endl;
              //cerr<<"link2_idx:"<<link2_idx<<endl;
              //cerr<<"Testing link:"<<endl;
                if (bondmap.canLink(pep1, pep2, link1_idx, link2_idx)) {
                  //create the candidate
                  XLinkMatch* newCandidate = 
                    new XLinkPeptide(pep1, pep2, link1_idx, link2_idx);
                  candidates.add(newCandidate);
                }
              }
            } /* for link1_idx */
          } /* if (mods <= max_mod_xlink .. */
        } /* if (candidate_mass >= min_mass && candidate_mass <= max_mass) */
      } /* if (!mixed_target_decoys_ ... */
      next_idx++;
    
      if (next_idx < linkable_peptides.size()) {
        current_mass = first_mass + linkable_peptides[next_idx].getMass();
      }
    }

    //start with the next candidate.
    first_idx++;
    //carp(CARP_INFO,"Incremented first to:%d",first_idx);
  }

  //cerr <<"XLinkPeptide::addCandidates: done"<<endl;
}

XLINKMATCH_TYPE_T XLinkPeptide::getCandidateType() {
  return XLINK_CANDIDATE;
}

string XLinkPeptide::getSequenceString() {

  doSort();

  string seq1 = linked_peptides_[0].getModifiedSequenceString();
  
  string seq2 = linked_peptides_[1].getModifiedSequenceString();




  //assert(seq1 <= seq2);

  ostringstream oss;
  oss << seq1 << ", " << 
    seq2 << " (" <<
    (getLinkPos(0)+1) << "," <<
    (getLinkPos(1)+1) << ")";

  string svalue = oss.str();

  return svalue;
}

FLOAT_T XLinkPeptide::calcMass(MASS_TYPE_T mass_type) {
  return linked_peptides_[0].getMass(mass_type) + 
    linked_peptides_[1].getMass(mass_type) + 
    linker_mass_;
}

XLinkMatch* XLinkPeptide::shuffle() {

  XLinkPeptide* decoy = new XLinkPeptide();
  decoy->setZState(getZState());
  decoy->linked_peptides_.push_back(linked_peptides_[0].shuffle());
  decoy->linked_peptides_.push_back(linked_peptides_[1].shuffle());
  decoy->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy->link_pos_idx_.push_back(link_pos_idx_[1]);

  return (XLinkMatch*)decoy;


}

void XLinkPeptide::predictIons(IonSeries* ion_series, int charge) {
  //cerr << "Inside predictIons"<<endl;
  MASS_TYPE_T fragment_mass_type = get_mass_type_parameter("fragment-mass"); 
  //cerr << "Predicting "<<getSequenceString()<<" +"<<charge<<endl;
  //predict the ion_series of the first peptide.
  char* seq1 = linked_peptides_[0].getSequence();
  MODIFIED_AA_T* mod_seq1 = linked_peptides_[0].getModifiedSequence();

  //carp(CARP_INFO,"predictIONS %s",seq1);

  ion_series->setCharge(charge);
  ion_series->update(seq1, mod_seq1);
  ion_series->predictIons();

  //iterate through all of the ions, if the ion contains a link, then
  //add the mass of peptide2 + linker_mass.

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;

    unsigned int cleavage_idx = ion->getCleavageIdx();

    if (ion->isForwardType()) {
      if (cleavage_idx > (unsigned int)getLinkPos(0)) {
	FLOAT_T mass = ion->getMassFromMassZ();
	mass += linked_peptides_[1].getMass(fragment_mass_type) + linker_mass_;
	ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN1");
        }
      }
    } else {
      if (cleavage_idx >= (strlen(seq1) - (unsigned int)getLinkPos(0))) {
	FLOAT_T mass = ion->getMassFromMassZ();
	mass += linked_peptides_[1].getMass(fragment_mass_type) + linker_mass_;
	ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN2");
        }
      }
    }
  }

  //predict the ion_series of the second peptide.
  IonConstraint* ion_constraint = ion_series->getIonConstraint();

  IonSeries* ion_series2 = 
      new IonSeries(ion_constraint, charge);
  
  char* seq2 = linked_peptides_[1].getSequence();

  //carp(CARP_INFO,"seq2:%s",seq2);

  MODIFIED_AA_T* mod_seq2 = 
    linked_peptides_[1].getModifiedSequence();
  ion_series2->setCharge(charge);
  ion_series2->update(seq2, mod_seq2);
  ion_series2->predictIons();

  //modify the necessary ions and add to the ion_series  
  for (IonIterator ion_iter = ion_series2->begin();
    ion_iter != ion_series2->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;

    unsigned int cleavage_idx = ion->getCleavageIdx();
    if (ion->isForwardType()) {
      if (cleavage_idx > (unsigned int)getLinkPos(1)) {
	FLOAT_T mass = ion->getMassFromMassZ();
	mass += linked_peptides_[0].getMass(fragment_mass_type) + linker_mass_;
	ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN3");
        }
      }
    } else {
      if (cleavage_idx >= (strlen(seq2)-(unsigned int)getLinkPos(1))) {
	FLOAT_T mass = ion->getMassFromMassZ();
	mass += linked_peptides_[0].getMass(fragment_mass_type) + linker_mass_;
	ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN4");
        }
      }
    }
    ion_series->addIon(ion);
  }
  //carp(CARP_INFO,"free(seq1)");
  free(seq1);
  //carp(CARP_INFO,"free(seq2)");
  free(seq2);
  //carp(CARP_INFO,"free(mod_seq1)");
  free(mod_seq1);
  //carp(CARP_INFO,"free(mod_seq2)");
  free(mod_seq2);
  
  IonSeries::freeIonSeries(ion_series2, false);

  //carp(CARP_INFO,"Number of ions:%d",get_ion_series_num_ions(ion_series));
  
}

string XLinkPeptide::getIonSequence(Ion* ion) {

  int peptide_idx = 0;

  string ion_sequence = ion->getPeptideSequence();

  if (ion_sequence == linked_peptides_[0].getSequence()) {
    peptide_idx = 0;
  } else {
    peptide_idx = 1;
  }

  unsigned int cleavage_idx = ion->getCleavageIdx();

  bool is_linked = false;

  if (ion->isForwardType()) {
    is_linked = (cleavage_idx > (unsigned int)getLinkPos(peptide_idx)); 
  } else {
    is_linked = (cleavage_idx >= (ion_sequence.length() - getLinkPos(peptide_idx)));
  }

  string subseq;
  if (ion->isForwardType()) {
    subseq = ion_sequence.substr(0, cleavage_idx);
  } else {
    subseq = ion_sequence.substr(ion_sequence.length() - cleavage_idx, ion_sequence.length());
  }

  if (!is_linked) {
    return subseq;
  } else {
    string ans;
    if (peptide_idx == 0) {
      char* seq2 = linked_peptides_[1].getSequence();
      ans = subseq + string(",") + string(seq2);
      free(seq2);
    } else {
      char* seq1 = linked_peptides_[0].getSequence();
      ans = string(seq1) + string(",") + subseq;
      free(seq1);
    }
    return ans;
  }
}

Peptide* XLinkPeptide::getPeptide(int peptide_idx) {
  return linked_peptides_[peptide_idx].getPeptide();
}

int XLinkPeptide::getNumMissedCleavages() {

  char missed_cleavage_link_site = 'K';
  set<int> skip;

  int link1_site = getLinkPos(0);
  int link2_site = getLinkPos(1);
  
  Peptide* pep1 = linked_peptides_[0].getPeptide();
  Peptide* pep2 = linked_peptides_[1].getPeptide();
  
  char *seq1 = pep1->getSequencePointer();
  char *seq2 = pep2->getSequencePointer();

  if (seq1[link1_site] == missed_cleavage_link_site) {
    skip.insert(link1_site);
  }

  int missed1 = pep1->getMissedCleavageSites(skip);
  
  skip.clear();

  if (seq2[link2_site] == missed_cleavage_link_site) {
    skip.insert(link2_site);
  }
  
  int missed2 = pep2->getMissedCleavageSites(skip);

  //cerr<<getSequenceString()<<" "<<missed1<<" "<<missed2<<endl;

  return max(missed1, missed2);



}


bool XLinkPeptide::isModified() {

  return linked_peptides_[0].isModified() || linked_peptides_[1].isModified();
}

string XLinkPeptide::getProteinIdString() {

  doSort();

  ostringstream oss;

  Peptide* peptide = this -> getPeptide(0);

  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide : Null first peptide!");
  } else {
    oss << peptide->getProteinIdsLocations();
  }
  oss << ";";

  peptide = this -> getPeptide(1);

  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide : Null second peptide!");
  } else {
    oss << peptide->getProteinIdsLocations();
  }

  return oss.str();
}

string XLinkPeptide::getFlankingAAString() {

  doSort();

  ostringstream oss;

  Peptide* peptide = this -> getPeptide(0);
  
  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide::getFlankingAAString() : Null first peptide!");
  } else {

    char* flanking_aas = peptide->getFlankingAAs();
    oss << flanking_aas;

  }

  oss << ";";

  peptide = this->getPeptide(1);

  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide::getFlankingAAString() : Null second peptide!");
  } else {

    char* flanking_aas = peptide->getFlankingAAs();
    oss << flanking_aas;
  }
  
  return oss.str();
}
