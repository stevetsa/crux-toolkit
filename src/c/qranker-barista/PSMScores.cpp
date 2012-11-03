#include <assert.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <set>
#include <vector>
#include <string>
#include <math.h>
using namespace std;
#include "PSMScores.h"

inline bool operator>(const PSMScoreHolder &one, const PSMScoreHolder &other) 
    {return (one.score>other.score);}

inline bool operator<(const PSMScoreHolder &one, const PSMScoreHolder &other) 
    {return (one.score<other.score);}



PSMScores::PSMScores()
{
    factor=1;
    neg=0;
    pos=0;
    posNow=0;
}

PSMScores::~PSMScores()
{

}

void PSMScores::clear()
{
  scores.clear();
  pos=neg=posNow = 0;
}

void PSMScores::add_psms(PSMScores& psms) {
  for (vector<PSMScoreHolder>::iterator iter = psms.begin();
    iter != psms.end();
    ++iter) {
    add_psm(*iter);
  }

}


/**
 * Percentage of target scores that are drawn according to the null.
 */
double PSMScores::pi0 = 1.0;

void PSMScores :: calc_factor()
{
  //cout << scores.size() << endl;
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      if(scores[i].label == 1)
	pos++;
      else
	neg++;
    }
  factor = (double)pos/(double)neg;
  //cout << pos << " " << neg << " " << factor << endl;
}

bool comparePValue(const PSMScoreHolder& s1, const PSMScoreHolder& s2) {

  if (s1.p == s2.p) {
    return s1.score > s2.score;
  }

  return s1.p < s2.p;
}


int PSMScores::calcOverFDRBH(double fdr) {

  sort(scores.begin(), scores.end(), comparePValue);

  int ans = 0;
  int posNow = 0;
  for (size_t idx = 0 ;idx < scores.size();idx++) {
    if (scores[idx].label == 1) {
      posNow++;
      scores[idx].q = (double)pos / (double)(posNow) * pi0 * scores[idx].p;
      //cout << "pos:" << pos<<" posNow:"<<posNow<<" p:"<<scores[idx].p<<" q:"<<scores[idx].q<<endl;
      if (scores[idx].q <= fdr) {
        ans = posNow;
      }
    }
  }

  double lastFDR = pi0;
  


  for (int ix=scores.size();ix>=0;--ix) {
    if (scores[ix].label == 1) {
      if (scores[ix].q < lastFDR) {
        lastFDR = scores[ix].q;
      } else {
        scores[ix].q = lastFDR;
      }
    }
  }

  return ans;

}


/**
 * Calculate the number of targets that score above a specified FDR.
 */
int PSMScores::calcOverFDR(double fdr) {
  
  vector<PSMScoreHolder>::iterator it = scores.begin();

  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
  
  int positives=0,nulls=0;
  double efp=0.0,q;
  posNow = 0;
  register int ix=0;

  map<double,double> score_to_fdr;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      positives++;
    if (it->label==-1) {
      nulls++;
      efp=pi0*nulls*factor;
    }
    if (positives)
      q=efp/(double)positives;
    else
      q=pi0;
    if (q>pi0)
      q=pi0;
    score_to_fdr[it->score] = q;
/*
    if (fdr>=q)
      posNow = positives;
*/
  }


  double min_fdr = pi0;
  //cerr << "assigning q-values"<<endl;
  for (ix=scores.size()-1;ix >=0;ix--) {
    //cerr << "ix:"<<ix << endl;
    double current_score = scores[ix].score;
    
    //cerr << "score:"<< current_score << endl;
    double current_fdr = score_to_fdr[current_score];
    //cerr << "fdr:"<<current_fdr<<endl;
    if (current_fdr < min_fdr) {
      min_fdr = current_fdr;
    }
    scores[ix].q = min_fdr;
    if (min_fdr <= fdr) {
      posNow++;
    }
  }
  

  return posNow;
}

void PSMScores::calcPValues() {

  vector<PSMScoreHolder>::iterator it = scores.begin();

  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());

  int nulls = 0;
  


  for ( vector<PSMScoreHolder>::iterator it = scores.begin();
    it != scores.end(); it++) {
      if (it->label == -1) {
        nulls++;
      }
      it->p = nulls;
   }

  for ( vector<PSMScoreHolder>::iterator it = scores.begin();
    it != scores.end();it++) {
    it-> p = it->p / (double)nulls;
  }



}


void PSMScores::getMaxPerScan(Dataset& d, PSMScores& max) {
  vector<PSMScoreHolder>::iterator it = scores.begin();

  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());

  max.clear();
  set<int>scans;

  for ( vector<PSMScoreHolder>::iterator it = scores.begin();
    it != scores.end(); it++) {
      int scan = d.psmind2scan(it->psmind);
      if (scans.find(scan) == scans.end()) {
        max.add_psm(*it);
        scans.insert(scan);
      }
  }

  max.calc_factor();

}


void PSMScores::calcMultiOverFDR(vector<double> &fdr, vector<int> &overFDR) {
  
  vector<PSMScoreHolder>::iterator it = scores.begin();

  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
  
  int positives=0,nulls=0;
  double efp=0.0,q;
  register unsigned int ix=0;
  
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label==1)
      positives++;
    if (it->label==-1) {
      nulls++;
      efp=pi0*nulls*factor;
    }
    if (positives)
      q=efp/(double)positives;
    else
      q=pi0;
    it->q=q;
    if (q>pi0)
      q=pi0;
       
    for(unsigned int ct = 0; ct < fdr.size(); ct++)
      if (fdr[ct]>=q)
	overFDR[ct] = positives;
  }
  for (ix=scores.size();--ix;) {
    if (scores[ix-1].q > scores[ix].q)
      scores[ix-1].q = scores[ix].q;  
  }
}

void PSMScores::fillFeaturesSplitScan(
  PSMScores& in, 
  Dataset& d, 
  PSMScores& train, 
  PSMScores& test) {

  map<int, pair<int,int> > scan_to_pos_neg;

  for (int idx = 0; idx < in.size(); idx++) {
    int scan = d.psmind2scan(in[idx].psmind);
    int label = d.psmind2label(in[idx].psmind);
    if (scan_to_pos_neg.find(scan) == scan_to_pos_neg.end()) {
      pair<int,int> pos_neg;
      scan_to_pos_neg[scan] = pos_neg;
    }

    if (label == 1) {
      scan_to_pos_neg[scan].first++;
    } else {
      scan_to_pos_neg[scan].second++;
    }
  }

  split(train, test, d, scan_to_pos_neg);

}

void PSMScores::fillFeaturesSplitScan(PSMScores& train, PSMScores& test, Dataset& d) {

  map<int, pair<int,int> > scan_to_pos_neg;

  for (int idx = 0; idx < d.get_num_psms();idx++) {
    int scan = d.psmind2scan(idx);
    int label = d.psmind2label(idx);
    if (scan_to_pos_neg.find(scan) == scan_to_pos_neg.end()) {
      pair<int,int> pos_neg;
      scan_to_pos_neg[scan] = pos_neg;
    }

    if (label == 1) {
      scan_to_pos_neg[scan].first++;
    } else {
      scan_to_pos_neg[scan].second++;
    }
  }

  split(train, test, d, scan_to_pos_neg);
}


void PSMScores::split(PSMScores& train, PSMScores& test, Dataset& d, map<int, pair<int,int> >& scan_to_pos_neg) {

  cerr << "number of scans:"<<scan_to_pos_neg.size()<<endl;

  map<pair<int,int>, vector<int> > pos_neg_to_scans;

  for (map<int, pair<int,int> >::iterator iter = scan_to_pos_neg.begin();
    iter != scan_to_pos_neg.end();
    ++iter) {

    pair<int,int> pos_neg = iter->second;
    int scan = iter->first;

    if (pos_neg_to_scans.find(pos_neg) == pos_neg_to_scans.end()) {
      vector<int> scans;
      pos_neg_to_scans[pos_neg] = scans;
    }
    pos_neg_to_scans[pos_neg].push_back(scan);
  }

  set<int> train_scans;
  set<int> test_scans;
  for (map<pair<int,int>, vector<int> >::iterator iter = pos_neg_to_scans.begin();
    iter != pos_neg_to_scans.end();
    ++iter) {

    pair<int,int> pos_neg = iter->first;

    vector<int> scans = iter->second;

    cerr << "pos: "<<pos_neg.first<<" neg:"<<pos_neg.second<<" scans:"<<scans.size()<<endl;
    //random_shuffle(scans.begin(), scans.end());
    int nscans_train = scans.size() / 2;
    for (int idx = 0; idx < nscans_train;idx++) {
      train_scans.insert(scans[idx]);
    }
    for (int idx = nscans_train; idx < scans.size();idx++) {
      test_scans.insert(scans[idx]);
    }

  } 

  for (int idx = 0 ; idx < d.get_num_psms();idx++) {
    int scan = d.psmind2scan(idx);
    PSMScoreHolder psm;
    psm.psmind = idx;
    psm.label = d.psmind2label(idx);
    if (train_scans.find(scan) != train_scans.end()) {
      train.add_psm(psm);
    } else if (test_scans.find(scan) != test_scans.end()) {
      test.add_psm(psm);
    }
  }

  train.calc_factor();
  test.calc_factor();

  cerr <<" train pos:"<<train.pos<<" train neg:"<<train.neg<<" test pos:"<<test.pos<<" test neg"<<test.neg<<endl;

  //exit(-1);

}
/*

void PSMScores::fillFeaturesSplitScan(PSMScores& train, PSMScores& test, Dataset& d) {

  set<int> scans;

  for (int i = 0; i < d.get_num_psms(); i++) {
    scans.insert(d.psmind2scan(i));
  }

  vector<int> scans_vec(scans.begin(), scans.end());
  //random_shuffle(scans_vec.begin(), scans_vec.end());


  int nscans = scans_vec.size();

  set<int> train_scans;

  int nscans_train = nscans/2;

  for (int idx = 0; idx < nscans_train;idx++) {
    train_scans.insert(scans_vec[idx]);
  }

  for (int idx = 0 ; idx < d.get_num_psms();idx++) {
    int scan = d.psmind2scan(idx);
    PSMScoreHolder psm;
    psm.psmind = idx;
    psm.label = d.psmind2label(idx);
    if (train_scans.find(scan) != train_scans.end()) {
      train.add_psm(psm);
    } else {
      test.add_psm(psm);
    }

  }

  train.calc_factor();
  test.calc_factor();

}
*/

/*
void PSMScores::fillFeaturesSplitScan(
  PSMScores& in, 
  Dataset& d, 
  PSMScores& train, 
  PSMScores& test) {

  set<int> scans;

  for (int i = 0; i < in.size(); i++) {
    scans.insert(d.psmind2scan(in[i].psmind));
  }

  vector<int> scans_vec(scans.begin(), scans.end());
  //random_shuffle(scans_vec.begin(), scans_vec.end());


  int nscans = scans_vec.size();

  set<int> train_scans;

  int nscans_train = nscans/2;

  for (int idx = 0; idx < nscans_train;idx++) {
    train_scans.insert(scans_vec[idx]);
  }

  for (int idx = 0 ; idx < in.size();idx++) {
    int scan = d.psmind2scan(in[idx].psmind);
    PSMScoreHolder psm;
    psm.psmind = idx;
    psm.label = d.psmind2label(in[idx].psmind);
    if (train_scans.find(scan) != train_scans.end()) {
      train.add_psm(psm);
    } else {
      test.add_psm(psm);
    }

  }

  train.calc_factor();
  test.calc_factor();

}
*/

void PSMScores::fillFeaturesSplit(PSMScores& train,PSMScores& test, Dataset& d, double ratio) {
 
  int num_psms = d.get_num_psms();
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  int n = num_psms;
  int k = (int)(n*ratio);
  int l = n-k;
 
  PSMScoreHolder s;
  
  //collect everybody 
  vector<PSMScoreHolder> all_examples;
  all_examples.resize(n,s);
  for (int i = 0; i < num_psms; i++)
    {
      all_examples[i].psmind = i;
      all_examples[i].label = d.psmind2label(i);
    }
  
  
  //mix up the examples
  for(int i = 0; i < n; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(n-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(n-1)); 
      s = all_examples[p1];
      all_examples[p1] = all_examples[p2];
      all_examples[p2] = s;
    }
  
  train.scores.resize(k,s);
  test.scores.resize(l,s);
  

  int num_pos_train = 0;
  int num_neg_train = 0;
  int num_pos_test = 0;
  int num_neg_test = 0;
  //distribute the normal set between train and test
  for (int i = 0; i < k; i++)
    {
      train.scores[i] = all_examples[i];
      if (train.scores[i].label == 1)
	num_pos_train++;
      else
	num_neg_train++;
    }
  for(int i=0; i < l; i++)
    {
      test.scores[i] = all_examples[i+k];
      if (test.scores[i].label == 1)
	num_pos_test++;
      else
	num_neg_test++;
    }
  cout << num_pos_train << " " << num_neg_train << " " << num_pos_test << " " << num_neg_test << "\n";
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=num_neg_train;
  test.neg=num_neg_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}


void PSMScores::fillFeaturesFull(PSMScores& full, Dataset& d) {
 
  int n = d.get_num_psms();
 
  PSMScoreHolder s;
  full.scores.resize(n,s);
  int num_pos = 0;
  int num_neg = 0;

  for (int i = 0; i < n; i++)
    {
      full.scores[i].psmind = i;
      full.scores[i].label = d.psmind2label(i);
      if (full.scores[i].label == 1)
	num_pos++;
      else
	num_neg++;
  }
  
  full.pos=num_pos;
  full.neg=num_neg;
  full.factor = full.pos/(double)full.neg;
}


void PSMScores::fillFeaturesBootstrap(PSMScores& bootstrap, Dataset& d) {

  int n = d.get_num_psms();

  PSMScoreHolder s;
  bootstrap.scores.resize(n,s);
  int num_pos = 0;
  int num_neg = 0;


  for (int idx = 0;idx < n; idx++) {
    //sample with replacement.
    int psm_idx = (int)((double)rand()/RAND_MAX*(n-1)); 
    bootstrap.scores[idx].psmind = psm_idx;
    bootstrap.scores[idx].label = d.psmind2label(psm_idx);
    if (bootstrap.scores[idx].label == 1) {
      num_pos++;
    } else {
      num_neg++;
    }

  }

  bootstrap.pos = num_pos;
  bootstrap.neg = num_neg;

  bootstrap.factor = (double)num_pos / (double)num_neg; 
  cerr <<" bootstrap: num_pos:"<<num_pos<<" num_neg:"<<num_neg<<" factor:"<<bootstrap.factor<<endl;
}

void PSMScores::fillFeaturesBootstrap(PSMScores& in, PSMScores& bootstrap) {
  bootstrap.clear();

  int num_pos = 0;
  int num_neg = 0;

  int n = in.scores.size();

  for (int idx = 0;idx < n;idx++) {
    //sample with replacement.
    int idx2 = (int)((double)rand()/RAND_MAX*(n-1));
    bootstrap.add_psm(in[idx2]);
    if (bootstrap.scores[idx].label == 1) {
      num_pos++;
    } else {
      num_neg++;
    }

  }

  bootstrap.pos = num_pos;
  bootstrap.neg = num_neg;

  bootstrap.factor = (double)num_pos / (double)num_neg; 
  cerr <<" bootstrap: num_pos:"<<num_pos<<" num_neg:"<<num_neg<<" factor:"<<bootstrap.factor<<endl;

}


bool compareRank(const PSMScoreHolder& s1, const PSMScoreHolder &s2) {
  return s1.rank < s2.rank;
}

bool compareScore(const PSMScoreHolder& s1, const PSMScoreHolder& s2) {

  return s1.score > s2.score;

}


void PSMScores::sortByRank() {
  sort(scores.begin(), scores.end(), compareRank);
}

void PSMScores::calcMinRank(PSMScores& out, Dataset& d) {

  out.clear();

  sortByRank();

  int neg=0;
  int pos=0;

  set<int> target_scans;
  set<int> decoy_scans;

  vector<PSMScoreHolder>::iterator iter;

  for (iter = begin();iter != end();++iter) {
    int label = iter->label;
    int psmind = iter->psmind;
    int scan = d.psmind2scan(psmind);
    if (label == 1 && target_scans.find(scan) == target_scans.end()) {
      out.add_psm(*iter);
      target_scans.insert(scan);
      pos++;
    } else if (label == -1 && decoy_scans.find(scan) == decoy_scans.end()) {
      out.add_psm(*iter);
      decoy_scans.insert(scan);
      neg++;
   
    }
  }
  
  out.pos = pos;
  out.neg = neg;
  out.factor = (double)pos / (double)neg;

}

int PSMScores::calculateOverFDRRank(double fdr) {

  sortByRank();
  int positives=0,nulls=0;
  double efp=0.0,q;
  posNow = 0;
  register unsigned int ix=0;
  for(vector<PSMScoreHolder>::iterator it=scores.begin();it!=scores.end();it++) {
    //cerr << "psmind:"<<it->psmind<<" rank:"<<it->rank<<" label:"<<it->label<<endl;
    if (it->label!=-1)
      positives++;
    if (it->label==-1) {
      nulls++;
      efp=pi0*nulls*factor;
    }
    if (positives)
      q=efp/(double)positives;
    else
      q=pi0;
    if (q>pi0)
      q=pi0;
    it->q=q;
    if (fdr>=q)
      posNow = positives;
  }
  for (ix=scores.size();--ix;) {
    if (scores[ix-1].q > scores[ix].q)
      scores[ix-1].q = scores[ix].q;  
  }
  return posNow;
}








