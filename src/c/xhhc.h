#ifndef XHHC_H
#define XHHC_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include "utils.h"
#include "carp.h"
#include "objects.h"
#include "peptide_constraint.h"
#include "peptide.h"
#include "protein.h"
#include "parameter.h"
#include "database.h"
#include "objects.h"
#include "crux-utils.h"
#include "utils.h"
#include "parse_arguments.h"

class Link;

class Peptide;

class LinkedPeptide {
  public:
    LinkedPeptide(char* peptide_A, char* peptide_B, int posA, int posB, FLOAT_T linker_mass, int charge);
    LinkedPeptide(Peptide& peptide, int charge) : charge_(charge) {
	peptides_.push_back(&peptide);
    }
    void add(Peptide& peptide) {
      peptides_.push_back(&peptide);
    }
    std::vector<Peptide*> peptides() { return peptides_;}
    int charge() {return charge_;}
    FLOAT_T get_mz(); // for each peptide, get mass
    LinkedPeptide(LinkedPeptide& other) {
      charge_ = other.charge();
      peptides_.assign(other.peptides().begin(), other.peptides().end());
    }  
    std::vector<std::pair<LinkedPeptide, LinkedPeptide> > split();
   private:
      int charge_;
      std::vector<Peptide*> peptides_;
}

class Peptide {
  public:
    Peptide(char* sequence);
    bool has_link_at(int index); 
    Link* link_at(int index);
    int length();	
    void add_link(int index, Link& link);
    std::vector<std::pair<LinkedPeptide, LinkedPeptide> > split_at(int index, int charge, LinkedPeptide* parent);
    
private:
    Link* links_;
    int length_;
    char* sequence_;
}
/*  
class Fragment {
public:
  // constructors
  Fragment(std::string& peptide, int charge) : peptide_(peptide), charge_(charge) {
    for (int i=0; i<peptide.length(); i++) links_.push_back(NULL);
    //for (int i=0; i<peptide.length(); i++) add_link(i, NULL);
  }
  
  Fragment(std::string peptide, std::vector<Link*> links, int charge) : 
	peptide_(peptide), charge_(charge), links_(links) {
  	for (int i=0; i < peptide.length(); i++) {
          if (has_link_at(i)) std::cout << "X"; else std::cout << " ";
        }
        std::cout << std::endl;     
  }
  //Fragment(Fragment& other);
  std::pair<Fragment, Fragment> split_at(int index);
  int length()                  	{ return peptide_.length(); }
  std::string& peptide()         	{ return peptide_; }
  bool has_link_at(int index)   	{ return links_[index] != NULL; }
  Link* link_at(int index)	     	{ return links_[index]; }
  void add_link(int index, Link& link);
  void add_link(Link& link) {links_.push_back(&link);}
  std::string print();
  std::string print_with_charge();
private:
  int charge_;
  std::vector<Link*> links_;
  std::string peptide_;
};
*/

class Link {
public:
  Link(Peptide& pepA, Peptide& pepB, FLOAT_T mass) :
	pepA_(pepA), pepB_(pepB), mass_(mass), traveled(false) {}
  Peptide* other(Peptide& p) { // 
    if (&p == &pepA_) return &pepB_;
    return &pepA_;
  }
  void move(Peptide& old_pep, Peptide& new_pep) {
    if (&old_pep == &pepA_) {
      pepA_ = new_pep;
    } else {
      pepB_ = new_pep;
   }
  } 
  FLOAT_T mass() { return mass_; }
  bool traveled;
private:
  FLOAT_T mass_;
  Peptide& pepA_;
  Peptide& pepB_;
};

//void split_all(std::vector<std::pair<Fragment, Fragment> >& all_pairs, Fragment& f);

int main(int argc, char** argv);

#endif
