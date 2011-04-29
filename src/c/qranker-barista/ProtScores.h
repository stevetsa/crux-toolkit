#ifndef PROTSCORES_H_
#define PROTSCORES_H_
#include <vector>
#include <algorithm>
#include <string>
using namespace std;
#include "DataSet.h"

class ProtScoreHolder{
public:
  double score; 
  int protind;
  vector<int> subset_protinds;
  double q;
  int label;
  int in_meta_prot;
  ProtScoreHolder():score(0.0),protind(0),q(0.0),label(0), in_meta_prot(0){;}
  ~ProtScoreHolder() {;}
};

class ProtScores
{
public:
	ProtScores();
	~ProtScores();
	void clear();
    vector<ProtScoreHolder>::iterator begin() {return scores.begin();}
    vector<ProtScoreHolder>::iterator end() {return scores.end();}    
    static double pi0;
    double factor;

    int calcOverFDR(double fdr);
    void calcMultiOverFDR(vector<double> &fdr, vector<int> &overFDR);
    inline ProtScoreHolder& operator[](int ix){return scores[ix];}
    void static fillProteinsFull(ProtScores& fullset, Dataset &d);
    int static traverse_peptide(Dataset &d, int pepind, int trn_tst, vector<int> &assignment_array);
    int static traverse_protein(Dataset &d, int protind, int trn_tst, vector<int> &assignment_array);
    void static fillProteinsSplit(ProtScores& train,ProtScores& test,Dataset &d, double ratio);
    void make_meta_set(Dataset &d);
    bool is_subset(Dataset &d, int protind1, int protind2);
    inline int size(){return scores.size();}
protected:
    int neg,pos,posNow;
    vector<ProtScoreHolder> scores;
};



#endif /*PROTSCORES_H_*/
