#ifndef HHC_H
#define HHC_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <sstream>

#include "utils.h"
#include "carp.h"
#include "peptide_constraint.h"
#include "peptide.h"
#include "peptide_src.h"
#include "protein.h"
#include "database.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "objects.h"
#include "crux-utils.h"
#include "ion.h"

// private stuff, needs to be cleaned up and simplified
class Peptide;
//class LinkedIonSeries;
using namespace std;
class LinkedPeptide {
  public:
    // constructor for a linked peptide. If A or B null, then
    // a self link will be created. If an index is -1, a link to nothing
    // will be created.
    LinkedPeptide(char* peptide_A, char* peptide_B, int posA, int posB, FLOAT_T linker_mass, int charge);
    LinkedPeptide(int charge, FLOAT_T linkermass) : charge_(charge), linker_mass_(linkermass), decoy_(false) {}
    //LinkedPeptide(const LinkedPeptide& other); 
    std::vector<Peptide>& peptides() 	{ return peptides_;}
    int charge() 			{ return charge_;}
    void set_charge(int charge)		{ charge_ = charge; }
    void set_type(ION_TYPE_T type)   	{ type_ = type; }
    ION_TYPE_T type()			{ return type_;}
    FLOAT_T linker_mass()		{ return linker_mass_;}
    FLOAT_T get_mz();
    int size()				{ return peptides_.size(); }
    void calculate_mass();
    FLOAT_T mass() 			{ return mass_; }
    void set_decoy()			{ decoy_ = true; }
    bool is_decoy()			{ return decoy_;}
    void add_peptide(Peptide& peptide) { peptides_.push_back(peptide); } 

  
    //void split_many(std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs);
    void split(std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& ion_pairs);
    friend std::ostream &operator<< (std::ostream& os, const LinkedPeptide& lp); 
    // for sorting LinkedPeptides by mz
    friend bool operator< (const LinkedPeptide &lp1, const LinkedPeptide &lp2) {
      return lp1.mass_ < lp2.mass_;
      //return lp1.mz < lp2.mz;
    }
   private:
      bool decoy_;
      ION_TYPE_T type_;
      FLOAT_T mass_;
      FLOAT_T mz;
      FLOAT_T linker_mass_;
      int charge_;
      std::vector<Peptide> peptides_;
};

class Peptide {
  public:
    Peptide() {}
    Peptide(std::string sequence) : sequence_(sequence), length_(sequence.length())  {}
    //void move_link(Peptide old_pep, Peptide& new_pep);
    bool has_link_at(int index);
    Peptide& link_at(int index);
    int length();
    //Peptide& shuffle();
    bool has_link() { return !links_.empty(); }
    std::string sequence() {return sequence_;}	
    void add_link(int index, Peptide& peptide);
    void split_at(int index, std::vector<std::pair<LinkedPeptide, LinkedPeptide> >& pairs, int charge, LinkedPeptide& parent);
    int link_site(); 
    FLOAT_T mass();
private:
    std::map<int, Peptide> links_;
    int length_;
    std::string sequence_;
};
/*
void Peptide::shuffle() {
  string copy = string(sequence_);
  int start_idx = 0;
  int end_idx = length_-1;
  int switch_idx = 0;
  char temp_char = 0;

  ++start_idx;
  --end_idx;

  while(start_idx < end_idx) {
    switch_idx = get_random_number_interval(start_idx, end_idx);
    //temp_char = sequence_[start_idx];
    copy[start_idx] = sequence_[switch_idx];
    copy[switch_idx] = sequence_[start_idx];
  }
}
*/

 
LinkedPeptide::LinkedPeptide(char* peptide_A, char* peptide_B, int posA, int posB, FLOAT_T linkermass, int charge) {
  charge_ = charge;
  linker_mass_ = linkermass;
  decoy_ = false;
  Peptide pepA = Peptide(peptide_A);
  // if a self link
  if (peptide_B == NULL) {
     if (posB != -1) {
       pepA.add_link(posB, pepA);
       pepA.add_link(posA, pepA);
     } else {
       Peptide dead_end = Peptide();
       pepA.add_link(posA, dead_end);
     }
    peptides_.push_back(pepA);
  } else {
    Peptide pepB = Peptide(peptide_B);
    //cout << "adding links at " << posA << " and " << posB << endl;
    pepA.add_link(posA, pepB);
    pepB.add_link(posB, pepA);
    peptides_.push_back(pepA);
    peptides_.push_back(pepB);
  }
}

// calculates mass of linked peptide,
// remove H2O from mass if it's a b-ion
// only works for one or two peptides for now
bool Peptide::has_link_at(int index) { return (links_.find(index) != links_.end()); }

int Peptide::link_site() {
  for (int i = 0; i < length_; ++i) {
    if (has_link_at(i)) 
      return i;
  }
  return -1;
}

FLOAT_T Peptide::mass() {
  return calc_sequence_mass((char*)sequence_.c_str(), MONO);
}

void LinkedPeptide::calculate_mass() {
  mass_ = calc_sequence_mass((char*)peptides_[0].sequence().c_str(), MONO);   
  /*for (vector<Peptide>::iterator peptide = peptides_.begin(); peptide != peptides_.end(); ++peptide) {
    mass_ += calc_sequence_mass((char*)peptide->sequence().c_str(), MONO);
   // change this later!!
    //if (peptide->has_link()) cout << "has link" << endl;
    if (peptide != peptides_.begin() || peptide->has_link()) {
     // cout << " adding linker mass " << endl;
      mass_ += linker_mass_;
      cout << "adding linker mass ";
    }
  } */
  if (size() == 2  || peptides_[0].has_link()) {
    mass_ += linker_mass_;
    if (size() == 2) mass_ += calc_sequence_mass((char*)peptides_[1].sequence().c_str(), MONO);
  }
  if (type_ == B_ION) {
    mass_ = mass_ - MASS_H2O_MONO;
  } 
}

// 
FLOAT_T LinkedPeptide::get_mz() {
  if (mz < 1)
    mz = ((mass_ + MASS_H_MONO*charge_) / charge_);
  return mz;
}

int main(int argc, char** argv);


Peptide& Peptide::link_at(int index)     { return links_[index]; }

int Peptide::length()           { return length_;}

void Peptide::add_link(int index, Peptide& other) {
  links_[index] = other;
}
/*
void Peptide::move_link(Peptide pep_old, Peptide& pep_new) {
  for(map<int, Peptide>::iterator i = links_.begin(); i != links_.end(); ++i)
  {
    if (i->second.sequence() == pep_old.sequence()) {
      i->second = pep_new;
      break;
    }
  }
}
*/

void LinkedPeptide::split(vector<pair<LinkedPeptide, LinkedPeptide> >& ion_pairs) {
 //vector<pair<LinkedPeptide, LinkedPeptide> > pairs; 
 for (vector<Peptide>::iterator it = peptides_.begin(); it != peptides_.end(); ++it) {
   Peptide peptide = *it;
   for (int i = 1; i < peptide.length(); i++) {
      peptide.split_at(i, ion_pairs, charge_, *this);
    }
  }
      //ion_pairs.insert(ion_pairs.end(), pairs.begin(), pairs.end()); // concat
}

// temporary
std::ostream &operator<< (std::ostream& os, LinkedPeptide& lp) {
  vector<Peptide> peptides = lp.peptides();
  ostringstream link_positions;
  link_positions << "(";
  for (int i = 0; i < peptides[0].length(); ++i) {
	if (peptides[0].has_link_at(i)) link_positions << i << "," ;
  }
  if (peptides.size() == 2) {
    for (int i = 0; i < peptides[1].length(); ++i) {
	if (peptides[1].has_link_at(i)) link_positions << i << ")";
    }
    return os << peptides[0].sequence() << ", " << peptides[1].sequence() << " " << link_positions.str() << " +" << lp.charge();
  }
  string final = link_positions.str();
  if (final.length() > 1) final.erase(final.end()-1);
  return os << peptides[0].sequence() << " " << final << ") +" << lp.charge();
}

// this needs to change
void Peptide::split_at(int index, vector<pair<LinkedPeptide, LinkedPeptide> >& pairs, int charge, LinkedPeptide& parent) {
  bool self_flag = false;
  for (int c = 0; c <= charge; ++c) {
    Peptide pepA = Peptide(sequence_.substr(0, index));
    Peptide pepB = Peptide(sequence_.substr(index, length_ - index));
    LinkedPeptide linkedA = LinkedPeptide(c, parent.linker_mass());
    LinkedPeptide linkedB = LinkedPeptide(charge - c, parent.linker_mass());
    self_flag = true;
    for (int i = 0; i < index; i++) {
      if (has_link_at(i)) {
        Peptide& other = link_at(i);
        pepA.add_link(i, other); 
        //other.move_link(*this, pepA);
        //linkedA.add_peptide(pepA);
        if (parent.size() > 1) linkedA.add_peptide(other);
        else if (!other.sequence().empty()) (self_flag = !self_flag);
      }
    }
    if (!self_flag) continue;
    for (int i = index; i < length_; i++) {
      if (has_link_at(i)) {
        Peptide other = link_at(i);
        pepB.add_link(i - index, other);
        //other.move_link(*this, pepB);
        //linkedB.add_peptide(pepB);
	if (parent.size() > 1) linkedB.add_peptide(other);
	//else (self_flag = !self_flag);
      }
    } 
    linkedA.add_peptide(pepA);
    linkedB.add_peptide(pepB);
    pairs.push_back(pair<LinkedPeptide, LinkedPeptide> (linkedA, linkedB));
  }
}



#endif
