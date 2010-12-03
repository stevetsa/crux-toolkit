#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <math.h>
using namespace std;
#include "Option.h"
#include "SanityCheck.h"
#include "SqtSanityCheck.h"
#include "DataSet.h"
#include "IntraSetRelation.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Net.h"
#include "ssl.h"
#include "Caller.h"
#include "Globals.h"
#include "MSReader.h"
#include "Spectrum.h"



namespace qranker {


string res_prefix="yeast_trypsin";



/********************* writing out results functions***************************************/

int Caller :: getOverFDR(Scores &set, NeuralNet &n, double fdr) {
  return getOverFDR(set, n, fdr, do_max_psm_);
}

int Caller :: getOverFDR(Scores &set, NeuralNet &n, double fdr, bool do_max_psm)
{

  calcScores(set, n);
  return set.calcOverFDR(fdr, do_max_psm);
}

void Caller :: calcScores(Scores &set, NeuralNet &n)
{
  double r = 0.0;
  const double* featVec;
  int label = 0;
  double loss = 0;

  for (unsigned int i = 0; i < set.size(); i++) 
    {
      featVec = set[i].pPSM->features;
      label = set[i].label;
      r = n.classify(featVec);
      //double err = n.get_err(label);
      //double err = 1-label*r;
      //if (err < 0)
      //	err = 0;
      //loss += err;
      set[i].score = r;
      set[i].pPSM->sc=r;

    }
  //cout << "loss " << loss << "\n";
}

void Caller :: calcPValues(Scores &set, NeuralNet &n, int &max_pos)
{
  calcScores(set, n);

  set.calcPValues(do_max_psm_, max_pos);
}

 

void Caller :: getMultiFDR(Scores &set, NeuralNet &n, vector<double> &qvalues, bool do_classify)
{
  double r = 0.0;
  const double* featVec;
  int label = 0;
  loss = 0;
 
  if (do_classify) {
    calcScores(set, n);
  }

  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti, do_max_psm_, do_classify);

}

void Caller :: setSigmoidZero(Scores &set, NeuralNet &n)
{
  
  int label;
  int q = overFDRmulti[0];
  //cout << "q " << q << "\n";
  int count = 0;
  unsigned int i;
  for (i = 0; i < set.size(); i++)
    {
      if (set[i].label == 1)
	count++;
      if (count>q)
	break;
    }
  
  label = set[i].label;
  double r = n.classify(set[i].pPSM->features);
  double err = n.get_err(label);
  //cout << "err " << err << "\n";
  n.adjust_bias_top_layer(r);
  r = n.classify(set[i].pPSM->features);
  err = n.get_err(label);
  //cout << "err " << err << "\n";

}


void Caller :: printNetResults(vector<int> &scores)
{
  double qv; int fdr;
  cerr << "QVALS SCORES:: ";
  for(unsigned int count = 0; count < qvals.size();count++)
    {  
      qv=qvals[count]; 
      fdr = scores[count];
      cerr << qv << ":" << fdr << " ";
    }
  //cerr << "loss " << loss << " ";
  cerr << endl;
}

void Caller :: write_max_nets(string filename, NeuralNet* max_net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  
  for(int count = 0; count < num_qvals; count++)
    {
      net = max_net[count];
      int r = getOverFDR(testset_,net, qvals[count]);
      int r1 = getOverFDR(trainset_,net, qvals[count]);
      f1 << qvals[count] << " " << r1 << " " << r << "\n";
      
      double qn;
      if(qvals[count] < 0.01)
	qn = 0.0012;
      else
	qn = 0.005;
      r = getOverFDR(testset_,net, qvals[count]+qn);
      r1 = getOverFDR(trainset_,net, qvals[count]+qn);
      f1 << qvals[count]+qn << " " << r1 << " " << r << "\n";
      
      //write out the net for a certain qvals
      ostringstream s1_net;
      s1_net << filename <<  "_net_" << qvals[count] << ".txt";
      ofstream f1_net(s1_net.str().c_str());
      net.write_to_file(f1_net);
      f1_net.close();
    }
  f1.close();
}


/********** xvalidation functions********************************************/

void Caller :: train_general_net(Scores &train, Scores &thresh, double qv)
{
  int overFDR = 0;
  int best_overFDR = 0;
  int best_iteration = 0;
  overFDR = getOverFDR(thresh,net, qv);
  ind_low = -1;

  NeuralNet best_net;
  best_net = net;
  best_overFDR = overFDR;

  
  for(int i=0;i<switch_iter;i++) {
    train_net_two(train);
    
    //record the best result
    overFDR = getOverFDR(thresh,net, qv);
    if(overFDR > best_overFDR)
      {
	best_overFDR = overFDR;
	best_net = net;
        best_iteration = i;
      }
      
    /*
    if((i % 100) == 0)
      {
	if(VERB>1) cerr << "Iteration " << i << " : \n";
	cerr << "train: ";
	cerr << net.get_mu() << " " << qv << " " << getOverFDR(train,net,qv) << "\n";
	cerr << "thresh: ";
	cerr << net.get_mu() << " " << qv << " " << getOverFDR(thresh,net,qv) << "\n";
	cerr << "testset: ";
	cerr << net.get_mu() << " " << qv << " " << getOverFDR(testset,net,qv) << "\n";
      }
    */
    
  }

  cerr<<"train_general_net(): best iter:"<<best_iteration<<" bestFDR("<<qv<<")="<<best_overFDR<<endl;

  net = best_net;
  best_net.clear();
}


void Caller :: train_target_net(Scores &train, Scores &thresh, double qv)
{
  
  int overFDR = 0;
  int best_overFDR = 0;
  
  overFDR = getOverFDR(train,net, qv);
  ind_low = overFDR;

  if (do_max_psm_) {
    ind_low = getOverFDR(train,net, qv, false);
  }

  NeuralNet best_net;
  best_net = net;
  best_overFDR = overFDR;
  //cerr << "target_q " << qv << " starting " << overFDR << "\n";
 
  for(int i = switch_iter; i < niter; i++)
    {
      overFDR = getOverFDR(train,net,qv);
      train_net_two(train);
      
      overFDR = getOverFDR(thresh,net,qv);
      if(overFDR > best_overFDR)
	{
	  best_overFDR = overFDR;
	  best_net = net;
	}
      /*
      if(i%100 == 0)
	{
	  cerr << i << " " << interval << "\n";
	  //cerr << "trainset: ";
	  //cerr << net.get_mu() << " " << qv << " " << getOverFDR(train,net,qv) << "\n";
	  cerr << "thresh: ";
	  cerr << net.get_mu() << " " << qv << " " << getOverFDR(thresh,net,qv) << "\n";
	  cerr << "testset: ";
	  cerr << net.get_mu() << " " << qv << " " << getOverFDR(testset,net,qv) << "\n";
	}
      */
    }
  //cerr << "end on interval " << interval << "\n";
  net = best_net;
  best_net.clear();
}

 
/**
 * Do cross-validation to select two hyperparameters: the learning
 * rate and the weight decay.  The optimization criterion is the
 * number of target scores below a specified q-value threshold.
 */
void Caller :: xvalidate_net(
  vector<Scores>& xv_train,
  vector<Scores>& xv_test,
  double qv ////< The q-value threshold -in
)
{
  
  vector <double> xv_wt;
  xv_wt.push_back(0.00001);xv_wt.push_back(0.00005);
  
  vector <double> xv_mu;
  //xv_mu.push_back(0.005);xv_mu.push_back(0.007);xv_mu.push_back(0.01);
  xv_mu.push_back(0.005);xv_mu.push_back(0.01);

  

  //split the trainset
  //for(unsigned int i = 0; i < xval_fold; i++)
   //cerr << "size " << trainset.size() << " " << xv_train[i].size() << " " << xv_test[i].size() << "\n";
  
  //set the linear flag: 1 if linear, 0 otherwise
  int lf = 0;
  //set whether there is bias in the linear units: 1 if yes, 0 otherwise
  int bs = 1;
  //cost linear flag indicating whether to use the sigmoid(0) or linear loss(1)
  int clf = 0;
 
  double best_sum = 0;
  int best_wt = 0;
  int best_mu = 0;
  cerr << "xvalidating \n";
  for(unsigned int h = 0; h < xv_wt.size();h++)
    {
      for(unsigned int m = 0; m < xv_mu.size(); m++)
	{
	  mu = xv_mu[m];
	  int overFDR = 0;
	  double sm=0;
	  for(unsigned int i = 0; i < xval_fold; i++)
	    {
	      net.initialize(FeatureNames::getNumFeatures(),num_hu,mu,clf,lf,bs);
	      net.set_weightDecay(xv_wt[h]);
	      cerr << "learning rate " << mu << " weight decay "<<xv_wt[h] << " xval set " << i << " size xval set " << xv_train[i].size() << "\n";
	      train_general_net(xv_train[i], xv_train[i],qv);

	      net.set_cost_flag(1);;
	      net.remove_bias();
	     
	      train_target_net(xv_train[i],xv_train[i],qv);
	      overFDR = getOverFDR(xv_test[i],net, qv);
	      sm +=overFDR;
	      net.clear();
	    }
	  sm /=xval_fold;
	  //cerr << "best_sum " << best_sum << " sum " << sm << "\n";
	  if(sm > best_sum)
	    {
	      best_sum = sm;
	      best_mu = m;
	      best_wt = h;
	    }
	}
    }
  mu = xv_mu[best_mu];
  weightDecay = xv_wt[best_wt];    
  
  cerr << "training target: mu " << best_mu << " " << mu << " wt_decay " << weightDecay << "\n";
   
}




/*********************** training net functions*******************************************/

void Caller :: train_net_two(Scores &set)
{
  double r1 = 0;
  double r2 = 0;
  double diff = 0;
  int label = 1;
  const double *featVec;

  double total_choose_time = 0;
  double total_classify_time = 0;
  double total_train_time = 0;
  double total_cn = 0;

  if (ind_low > 0) {
  /*
    vector<int> positive_indices;
    vector<int> negative_indices;
    int interval = min(set.size()-1, (unsigned int)ind_low*2);
    for (int i = 0;i<interval;i++) {
      if (set[i].label == 1) {
        positive_indices.push_back(i);
      } else {
        negative_indices.push_back(i);
      }
    }
    cerr<<"Looking at "<<interval<<" of "<<set.size()<<endl;
    cerr<<"Number of positives:"<<positive_indices.size()<<endl;
    cerr<<"Number of negatives:"<<negative_indices.size()<<endl;
  */
  }

  

  for(unsigned int i = 0; i < set.size(); i++)
    { 
      if(ind_low<0)
	{
	  unsigned int ind;
	  ind = rand()%set.size();
	  featVec = set[ind].pPSM->features;
	  label = set[ind].label;
	  //if (label == 1)
	    net.train(featVec,label);
	    //else {
	    //r1 = net.classify(set[ind].pPSM->features);
	    //if ((1-r1*label) > 0)
	    //  net.train1(set[ind].pPSM->features, set[ind].label); }
	}
      else
	{

          clock_t start_clock = clock();
	  interval = ind_low*2;
	  unsigned int ind1, ind2;
	  int label_flag = 1;
	  //get the first index
	  if(interval == 0)
	    ind1 = 0;
	  else
	    ind1 = rand()%(interval);
	  if(ind1<0 || ind1>set.size()-1) continue;
	  if(set[ind1].label == 1)
	    label_flag = -1;
	  else
	    label_flag = 1;
	  
	  int cn = 0;
	  while(1)
	    {
	      if(interval == 0)
		ind2 = rand()%set.size();
	      else
		ind2 = rand()%(interval);
	      if(ind2<0 || ind2>set.size()-1) continue;
	      if(set[ind2].label == label_flag) break;
	      if(cn > 1000)
		{
		  ind2 = rand()%set.size();
		  break;
	      	}
	      cn++;
	    }

          total_cn += cn;

	  clock_t choose_clock = clock();
	  
	  //pass both through the net
	  r1 = net.classify(set[ind1].pPSM->features);
	  r2 = net.classify(set[ind2].pPSM->features);
	  diff = r1-r2;

          clock_t classify_clock = clock();

	  label=0;

	  if(  set[ind1].label==1 && set[ind2].label==-1)
	    label=1;
	  
	  if( set[ind1].label==-1 && set[ind2].label==1)
	    label=-1;

	  
	  if(label != 0)
	    {
	      if(label*diff<1)
		{
		  net.train1(set[ind1].pPSM->features, set[ind1].label);
		  net.train1(set[ind2].pPSM->features, set[ind2].label);
		}
	    }
          clock_t train_clock = clock();
          total_choose_time += (double)(choose_clock - start_clock) / CLOCKS_PER_SEC;
          total_classify_time += (double)(classify_clock - choose_clock)/CLOCKS_PER_SEC;
          total_train_time += (double)(train_clock-classify_clock)/CLOCKS_PER_SEC;
	}
    }
    if (ind_low >= 0) {  
    /*
      cerr<<"Average choose:"<<(total_choose_time / (double)set.size())<<endl;
      cerr<<"Average classify:"<<(total_classify_time / (double)set.size()) << endl;
      cerr<<"Average train:"<<(total_train_time / (double)set.size())<<endl;
      cerr<<"Average cn:"<<total_cn / (double)set.size()<<endl;
      cerr<<"ind_low:"<<ind_low<<endl;
      double total = total_choose_time + total_classify_time + total_train_time;

      cerr<<"Ratio choose:"<<(total_choose_time/total)<<endl;
      cerr<<"Ratio classify:"<<(total_classify_time/total)<<endl;
      cerr<<"Ratio train:"<<(total_train_time/total)<<endl;
    */
    }
}




void Caller :: train_many_general_nets(
  Scores& trainset, 
  Scores& testset,
  Scores& thresholdset)
{
  ind_low = -1;
  for(int i=0;i<switch_iter;i++) {
    train_net_two(trainset);
    
    //record the best result
    getMultiFDR(thresholdset,net,qvals);
    for(int count = 0; count < num_qvals;count++)
      {
	if(overFDRmulti[count] > max_overFDR[count])
	  {
	    max_overFDR[count] = overFDRmulti[count];
	    max_net_gen[count] = net;
	  }
      }
    if((i % 30) == 0)
      {
	
	cerr << "Iteration " << i << " : \n";
        printResults(trainset, testset);
      }
  }

}


void Caller :: printResults(Scores& trainset, Scores& testset) {


  cerr << "trainset: ";
  getMultiFDR(trainset,net,qvals);
  printNetResults(overFDRmulti);
  cerr << "\n";
  cerr << "testset: ";
  getMultiFDR(testset,net,qvals);
  printNetResults(overFDRmulti);
  cerr << "\n";

  if (do_max_psm_) {
    do_max_psm_ = false;
    cerr << "trainset(orig): ";
    getMultiFDR(trainset,net,qvals);
    printNetResults(overFDRmulti);
    cerr << "\n";
    cerr << "testset(orig): ";
    getMultiFDR(testset,net,qvals);
    printNetResults(overFDRmulti);
    cerr << "\n";
    do_max_psm_ = true;
  }

}

void Caller :: printParameters() {
  cerr << "Hidden Units:"<<num_hu<<endl;
  cerr << "Learning Rate:"<<mu<<endl;
  cerr << "Weight Decay:"<<weightDecay<<endl;
  cerr << "Switch Iter:"<<switch_iter<<endl;
  cerr << "Total Iter:"<<niter<<endl;
  cerr << "do_xval:"<<do_xval<<endl;
  cerr << "do_max_psm_:"<<do_max_psm_<<endl;
  cerr << "do_pvalue:"<<do_pvalue<<endl;
}


void Caller :: train_many_target_nets_ave(
  Scores& trainset,
  Scores& testset,
  Scores& thresholdset) 
{

  int  thr_count = num_qvals-1;
  while (thr_count > 0)
    {
      net = max_net_gen[thr_count];
      //net = max_net_targ[thr_count];
      net.set_cost_flag(1);;
      net.remove_bias();
           
      //cerr << "training thresh " << thr_count << " mu " << net.get_mu() << "\n";
      cerr << "training thresh " << thr_count  << "\n";

      ind_low = getOverFDR(trainset, net, qvals[thr_count], false);
      cerr << "ind low:"<<ind_low<<endl;

      cerr <<" Before Iterating:"<<endl;
      printResults(trainset, testset);

      for(int i=switch_iter;i<niter;i++) {
		
	//sorts the examples in the training set according to the current net scores
	getMultiFDR(trainset,net,qvals);
	train_net_two(trainset);
	
	//record best result according to the average around q-value of interest
	for(int count = 0; count < num_qvals; count++)
	  ave_overFDR[count] = 0;
	//record the best result
	getMultiFDR(thresholdset,net, qvals);
	for(int count = 0; count < num_qvals; count++)
	  ave_overFDR[count] += overFDRmulti[count];
		
	//getMultiFDR(thresholdset,net, qvals1, false);
	//for(int count = 0; count < num_qvals; count++)
	//  ave_overFDR[count] += overFDRmulti[count];
	
	//getMultiFDR(thresholdset,net, qvals2, false);
	//for(int count = 0; count < num_qvals; count++)
	//  ave_overFDR[count] += overFDRmulti[count];
	
	//for(int count = 0; count < num_qvals; count++)
	//  ave_overFDR[count] /=3;
	  
	for(int count = 0; count < num_qvals;count++)
	  {
	    if(ave_overFDR[count] > max_overFDR[count])
	      {
                //cerr<<"Updating net:"<<count<<endl;
		max_overFDR[count] = ave_overFDR[count];
		max_net_targ[count] = net;
	      }
	  }

	if((i % 10) == 0)
        {
          
          cerr << "Iteration " << i << " : \n";
          printResults(trainset, testset);
        }

      }
      thr_count -= 2;
    }
}

void Caller::train_many_nets() 
{
  int max_pos = 0;
  int full_max_pos = 0;
  train_many_nets(trainset_, testset_, trainset_xv_train, trainset_xv_test, max_pos);
  full_max_pos = max_pos;
  
  if (do_pvalue) {
    train_many_nets(testset_, trainset_, testset_xv_train, testset_xv_test, max_pos);
    full_max_pos += max_pos;
    fullset.calcQValues(full_max_pos);
  } else {
    fullset.calcFDR_Decoy(do_max_psm_);
  }
}

void Caller::train_many_nets(
  Scores& trainset,
  Scores& testset,
  vector<Scores>& xv_train,
  vector<Scores>& xv_test,
  int& max_pos) 
{

  printParameters();
  
  thresholdset_ = trainset;

  switch_iter = 61;
  niter = 91;
   
  num_qvals = 14;
  qvals.clear();
  qvals.resize(num_qvals,0.0);
  qvals1.clear();
  qvals1.resize(num_qvals,0.0);
  qvals2.clear();
  qvals2.resize(num_qvals,0.0);
 
  overFDRmulti.clear();
  overFDRmulti.resize(num_qvals,0);
  ave_overFDR.clear();
  ave_overFDR.resize(num_qvals,0);
  max_overFDR.clear();
  max_overFDR.resize(num_qvals,0);

 
  max_net_gen = new NeuralNet[num_qvals];
  max_net_targ = new NeuralNet[num_qvals];
  
  double q = 0.0;
  for(int count = 0; count < num_qvals; count++)
    {
      qvals[count] = q;
      if (count < 2)
	qvals1[count] = q;
      else
	qvals1[count] = q-0.005;
      qvals2[count] = q+0.005;

      if(q < 0.01)
	q+=0.0025;
      else
	q+=0.01;
    }
  
  //set the linear flag: 1 if linear, 0 otherwise
  int lf = 0;
  if(num_hu == 1)
    lf = 1;
  //set whether there is bias in the linear units: 1 if yes, 0 otherwise
  int bs = 1;
  //cost linear flag indicating whether to use the sigmoid(0) or linear loss(1)
  int clf = 0;

  if (do_xval) {
    xvalidate_net(xv_train, xv_test, selectionfdr);
  }

  net.initialize(FeatureNames::getNumFeatures(),num_hu,mu,clf,lf,bs);
  net.set_weightDecay(weightDecay);
  for(int count = 0; count < num_qvals; count++){
    max_net_gen[count] = net;
  }
  initial_net = net;
  
  cerr << "Before iterating\n";
  cerr << "trainset: ";
  getMultiFDR(trainset,net,qvals);
  printNetResults(overFDRmulti);
  //setSigmoidZero(trainset,net);
  ofstream f_net_i("net_i.txt");
  net.write_to_file(f_net_i);
  f_net_i.close();

  cerr << "\n";
  cerr << "testset: ";
  getMultiFDR(testset,net,qvals);
  printNetResults(overFDRmulti);
  cerr << "\n";
  
  for(int count = 0; count < num_qvals;count++)
    max_overFDR[count] = overFDRmulti[count];


  train_many_general_nets(trainset, testset, thresholdset_);
  ofstream f_net("net.txt");
  net.write_to_file(f_net);
  f_net.close();

  //write out the results of the general net
  if (0){
    cerr << "general net results: ";
    ostringstream filename;
    filename << res_prefix << "_hu" << num_hu;
    write_max_nets(filename.str(), max_net_gen);
  }

  //copy the general net into target nets;
  for(int count = 0; count < num_qvals; count++)
    {
      if (getOverFDR(trainset, max_net_gen[count], qvals[count]) > getOverFDR(trainset, initial_net, qvals[count])) 
	max_net_targ[count] = max_net_gen[count];
      else
	{
	  //cout << count << " selecting init net\n";
	  max_net_targ[count] = initial_net;
	}
    }

  //calculate average of target q-values
  //and q-value +0.05 and -0.05
  getMultiFDR(thresholdset_,net, qvals);
  for(int count = 0; count < num_qvals; count++)
    ave_overFDR[count] += overFDRmulti[count];
  printNetResults(ave_overFDR);
  //getMultiFDR(thresholdset_,net, qvals1, false);
  //for(int count = 0; count < num_qvals; count++)
  //  ave_overFDR[count] += overFDRmulti[count];
  //printNetResults(ave_overFDR);
  //getMultiFDR(thresholdset_,net, qvals2, false);
  //for(int count = 0; count < num_qvals; count++)
  //  ave_overFDR[count] += overFDRmulti[count];
  //for(int count = 0; count < num_qvals; count++)
  //  ave_overFDR[count] /=3;
  //cerr<<"ave:"<<endl;
  //printNetResults(ave_overFDR);
  //for(int count = 0; count < num_qvals;count++)
  //  max_overFDR[count] = ave_overFDR[count];
  

  train_many_target_nets_ave(trainset, testset, thresholdset_);  
   
  //write out the results of the target net
  //cerr << "target net results: ";
  //ostringstream s2;
  //s2 << res_prefix << "_hu" << num_hu  << "_targ";
  //write_max_nets(s2.str(), max_net_targ);

  //choose the best net for the selectionfdr
  int max_fdr = 0;
  int fdr = 0;
  int ind = 0;
  cerr <<"Choosing best net"<<endl;
  for(unsigned int count = 0; count < qvals.size();count++)
    {
      fdr = getOverFDR(thresholdset_, max_net_targ[count], selectionfdr);

      if(fdr > max_fdr)
	{
	  max_fdr = fdr;
	  ind = count;
	}
    }
  net = max_net_targ[ind];
  cerr<<"Best :"<<ind<<":"<<max_fdr<<endl;


  //calculate pvalues on testset.
  if (do_pvalue) {
    cerr << " Found " << getOverFDR(testset, net, selectionfdr) << " over q<" << selectionfdr << "\n";
    cerr <<"Calculating pvalues"<<endl;
    calcPValues(testset, net, max_pos);
  } else {
    cerr << " Found " << getOverFDR(fullset, net, selectionfdr) << " over q<" << selectionfdr << "\n";
  }
  delete [] max_net_gen;
  delete [] max_net_targ;

}

} // qranker namspace

