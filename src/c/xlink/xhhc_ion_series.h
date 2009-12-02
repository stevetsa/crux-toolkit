#ifndef XHHC_ION_SERIES_H
#define XHHC_ION_SERIES_H

#include "xhhc.h"

//CRUX includes
extern "C" {
#include "ion_series.h"
#include "scorer.h"
#include "spectrum.h"
#include "spectrum_collection.h"
}

// should this be somewhere else?
#define bin_width_mono 1.0005079

using namespace std;

class LinkedIonSeries {
  public:
    // constructors
    LinkedIonSeries() : charge_(0) {}
    LinkedIonSeries(char* links, int charge);
    //LinkedIonSeries(char* sequenceA, char* sequenceB, int posA, int posB, int charge);
    //LinkedIonSeries(char* sequenceA, char* sequenceB, char* links, int charge);

    // getters
    int charge()                  { return charge_; }
    vector<LinkedPeptide>& ions() { return all_ions; }
    int size()                    { return all_ions.size(); }
    
    // other
    void set_charge(int charge) { charge_ = charge; }
    void clear()                { all_ions.clear(); }
    // splits
    void add_linked_ions(LinkedPeptide& linked_peptide);
    // print tab-delimited list of all ions
    void print();

  private:
    
    std::map<char, std::set<char> > bond_map;
    int charge_; 
    // a list of all the ions 
    std::vector<LinkedPeptide> all_ions;

};


void hhc_predict_ions(
  ION_SERIES_T* ion_series, ///< the ion series to predict ions for -in
  FLOAT_T linker_mass,
  int linker_site);

#endif
