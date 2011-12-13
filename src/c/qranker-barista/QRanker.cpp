#include "QRanker.h"

double pvalue_m_estimate = 0.0;

bool QRanker::no_delta_cn;

QRanker::QRanker() :  
  seed(0),
  selectionfdr(0.01),
  num_hu(5),mu(0.005),
  weightDecay(0.0000),
  max_net_gen(0),
  max_net_targ(0),
  nets(0),
  trained_nets(0),
  extra_fnames_exist(0)
{
}

QRanker::~QRanker()
{
  delete [] max_net_gen;
  delete [] max_net_targ;
  delete [] nets;
  delete [] trained_nets;
}

int QRanker :: getOverFDR(PSMScores &set, NeuralNet &n, double fdr)
{
  double *r;
  double* featVec;

  for(unsigned int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = n.fprop(featVec);
      set[i].score = r[0];

    }
  return set.calcOverFDR(fdr);
}

int QRanker :: getOverFDRPvalue(PSMScores &set, double fdr)
{

  set.sortByPValue();

  double total_targets = set.pos;

  double target_count = 0;

  for (unsigned int i=0; i < set.size();i++) {
      if (set[i].label == 1) {
        target_count++;
        //cerr << total_targets << "/" << 
        //        target_count << "*" << 
        //        PSMScores::pi0 <<"*"<<set[i].p;  
        set[i].q = total_targets / target_count * PSMScores::pi0 * set[i].p;
        //cerr << "="<<set[i].q<<endl;
      } else {
        set[i].q = PSMScores::pi0;
      }
  }

  int posNow = 0;

  double minFDR=PSMScores::pi0;

  for (int ix=set.size()-1;ix>=0;ix--) {
    if (set[ix].label == 1) {  
      if (set[ix].q >= minFDR) {
        set[ix].q = minFDR;
      } else {
        minFDR = set[ix].q; 
      }
      if (set[ix].q < fdr) {
        posNow++;
      }
    }
  }
  return posNow;
}

void QRanker :: calcPvalue(
  PSMScores &target_set, 
  PSMScores& target_decoy_set
  ) {
  //treat every example in the target set as a target, regardless of label.
  //use only the decoys in the target_decoy_set as decoys.
  cerr <<"Sort target set"<<endl;
  target_set.sortByScore();
  cerr <<"Sort target-decoy set"<<endl;
  target_decoy_set.sortByScore();

  cerr <<"Extracting decoy scores"<<endl;
  vector<double> decoy_scores;
  for (unsigned int idx=0;idx < target_decoy_set.size();idx++) {
    if (target_decoy_set[idx].label == -1) {
      decoy_scores.push_back(target_decoy_set[idx].score);
    }
  }

  cerr <<"Assigning pvalues"<<endl;
  unsigned int decoy_idx = 0;
 

  for (unsigned int target_idx = 0;target_idx < target_set.size();target_idx++) {
    
    double target_score = target_set[target_idx].score;
    while ((decoy_idx < decoy_scores.size()) && 
      (decoy_scores[decoy_idx] >= target_score)) {
      decoy_idx++;
    }

    target_set[target_idx].p = ((double)(decoy_idx + pvalue_m_estimate)) / 
      ((double)decoy_scores.size() + pvalue_m_estimate);
  }

}


void QRanker :: calcPvalue(PSMScores &set)
{

  //sort(set.begin(),set.end());
  //reverse(set.begin(),set.end());

  int decoy_count = 0;

  for(unsigned int i = 0; i < set.size(); i++) {
    if (set[i].label == -1) {
      decoy_count++;
    }
    set[i].p = (double)(decoy_count+pvalue_m_estimate)/(double)(set.neg+pvalue_m_estimate);
  }
}

void QRanker :: getMultiFDR(PSMScores &set, NeuralNet &n, vector<double> &qvalues)
{
  double *r;
  double* featVec;
   
  for(unsigned int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = n.fprop(featVec);
      set[i].score = r[0];
    }
 
  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti);
}


void QRanker :: printNetResults(vector<int> &scores)
{
  double qv; int fdr;
  cerr << "QVALS SCORES:: ";
  for(unsigned int count = 0; count < qvals.size();count++)
    {  
      qv=qvals[count]; 
      fdr = scores[count];
      cerr << qv << ":" << fdr << " ";
    }
  cerr << endl;
}

void QRanker :: write_results_pvalue(string prefix) {
  
  d.load_psm_data_for_reporting_results();
  
  net.copy(trained_nets[1]);
  getOverFDR(trainset, trained_nets[1], qvals[4]);
  calcPvalue(trainset);
  write_results(prefix, "train", trainset);
  
  net.copy(trained_nets[0]);
  getOverFDR(testset, net, qvals[4]);
  calcPvalue(testset);
  write_results(prefix, "test", testset);
    
  //PSMScores::fillFeaturesFull(fullset, d);
  fullset.combineSets(trainset,testset);
  //trainset.clear();
  //testset.clear();
  //thresholdset.clear();
  int r = getOverFDRPvalue(fullset, qvals[4]);
  write_results(prefix, "full", fullset);
  write_num_psm_per_spectrum(fullset, r);

  //treat the extra file as another decoy, we are going to use
  //this to calculate p-values/q-values on the peptide-level
  if(extra_fnames_exist) {
    PSMScores extra_trainset;
    PSMScores extra_testset;
    PSMScores extra_fullset;
    
    string summary_fn = "summary.txt";
    string psm_fn = "psm.txt";

    pars.clean_up(out_dir);
    
    if(!pars.run(extra_fnames))
      return;
    pars.clear();
    
    //Dataset d2;
    d.clear();
    d.load_psm_data_for_training(summary_fn, psm_fn);
    d.normalize_psms();
    
    cerr<<"Splitting extra files:"<<d.get_num_psms()<<endl;
    //PSMScores::fillFeaturesFull(extra_trainset, d2);
    PSMScores::fillFeaturesSplit(extra_trainset, extra_testset, d, 0.5);

    cerr <<"Scoring extra trainset:"<<extra_trainset.size()<<endl;
    getOverFDR(extra_trainset, trained_nets[1], qvals[4]);
    cerr <<"Calculating pvalues"<<endl;
    calcPvalue(extra_trainset, trainset);
    
    cerr <<"Scoring extra testset:"<<extra_testset.size()<<endl;
    getOverFDR(extra_testset, trained_nets[0], qvals[4]);
    
    cerr <<"Calculating test pvalues"<<endl;
    calcPvalue(extra_testset, testset);

    cerr <<"Writing results train"<<endl;
    d.load_psm_data_for_reporting_results();
    write_results(prefix, "extra.train", extra_trainset, d);
    
    cerr <<"Writing results test"<<endl;
    write_results(prefix, "extra.test", extra_testset, d);
    
    cerr <<"Writing extra full"<<endl;
    extra_fullset.combineSets(extra_trainset, extra_testset);

    getOverFDRPvalue(extra_fullset, qvals[4]);
    write_results(prefix, "extra.full", extra_fullset, d);

    extra_trainset.clear();
    extra_testset.clear();
    extra_fullset.clear();
  }

  cerr <<"Done write_results_pvalue"<<endl;
}


void QRanker :: write_results_extra(string prefix, NeuralNet& net)
{
  d.load_psm_data_for_reporting_results();
  getOverFDR(fullset, net, qvals[4]);
  write_results(prefix, "extra", fullset);
}    


void QRanker :: write_results(string prefix, NeuralNet& net) {
  
  d.load_psm_data_for_reporting_results();
  getOverFDR(trainset, net, qvals[4]);
  write_results(prefix, "train", trainset);
  getOverFDR(testset, net, qvals[4]);
  write_results(prefix, "test", testset);
    
  //PSMScores::fillFeaturesFull(fullset, d);
  fullset.combineSets(trainset,testset);
  trainset.clear();
  testset.clear();
  thresholdset.clear();
  int r = getOverFDR(fullset, net, qvals[4]);
  write_results(prefix, "full", fullset);
  write_num_psm_per_spectrum(fullset, r);
  
}

void QRanker :: write_results(
  string prefix,
  string filename,
  PSMScores& set) {

  write_results(prefix, filename, set, d);
}


void QRanker :: write_results(
  string prefix, 
  string filename, 
  PSMScores& set,
  Dataset& dataset
  ) 
{

  ostringstream s1;
  s1 << prefix << "." << filename << ".txt";
  ofstream f1(s1.str().c_str());

  f1 << "scan" << "\t" 
     << "charge" << "\t" 
     << "spectrum neutral mass" << "\t" 
     << "peptide mass" << "\t" 
     << "q-ranker score" << "\t" 
     << "q-ranker p-value" << "\t"
     << "q-ranker q-value" << "\t" 
     << "sequence" << "\t" 
     << "rtime max diff" << "\t"
     << "peptides/spectrum" << "\t"
     << "nzstates" << "\t"
     << "target/decoy"<<endl;

  for(unsigned int i = 0; i < set.size(); i++)
    {
      //cout<<"Writing "<<i<<endl;
      int psmind = set[i].psmind;
      //cout<<"psmind:"<<psmind<<endl;
      int num_pep = dataset.psmind2num_pep(psmind);
      //cout<<"num_pep:"<<num_pep<<endl;
      //write scan
      f1 << dataset.psmind2scan(psmind) << "\t";
      //write charges
      int *charges = dataset.psmind2charges(psmind);
      f1 << charges[0];
      for(int k = 1; k < num_pep; k++)
	{
	  f1 << ",";
	  int ch = charges[k];
	  f1 << ch; 
	}
      f1 << "\t";
      //write neutral_mass
      double* neutral_mass = dataset.psmind2neutral_mass(psmind);
      f1 << neutral_mass[0];
      for(int k = 1; k < num_pep; k++)
	{
	  f1 << ",";
	  double m = neutral_mass[k];
	  f1 << m; 
	}
      f1 << "\t";

      //write peptide mass
      double* peptide_mass = dataset.psmind2peptide_mass(psmind);
      f1 << peptide_mass[0];
      for(int k = 1; k < num_pep; k++)
	{
	  f1 << ",";
	  double m = peptide_mass[k];
	  f1 << m; 
	}
      f1 << "\t";
      f1 << set[i].score << "\t";
      f1 << set[i].p << "\t";
      f1 << set[i].q << "\t";

      //write peptide sequence
      int *pepinds = dataset.psmind2pepinds(psmind);
      f1 << dataset.ind2pep(pepinds[0]);
      for(int k = 1; k < num_pep; k++)
	{
	  f1 << ",";
	  f1 << dataset.ind2pep(pepinds[k]); 
	}

      //write rtime_max_diff
      double rtime_max_diff = dataset.psmind2rtime_max_diff(psmind);
      f1 << "\t" << rtime_max_diff;
            
      //write peptides/spectrum
      double peptides_spectrum = dataset.psmind2num_pep(psmind);
      f1 << "\t" << peptides_spectrum;

      //write nzstates
      int nzstates = dataset.psmind2nzstates(psmind);
      f1 << "\t" << nzstates;

      //write label
      f1 << "\t" << set[i].label;
      f1 << endl;


      
    }
  //cout<<"Done writing"<<endl;
  f1.close();
  //cout<<"Done..."<<endl;  

    


} 

void QRanker :: write_results_max(string filename, NeuralNet &net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
 
  set<int> scans_target;
  set<int> scans_decoy;
  
  for(unsigned int i = 0; i < fullset.size(); i++)
    {
      int scan = d.psmind2scan(fullset[i].psmind);
      if(fullset[i].label == 1)
	{
	  if(scans_target.find(scan) == scans_target.end())
	    {
	      fullset_max.add_psm(fullset[i]);
	      scans_target.insert(scan);
	    }
	}
      else
	{
	  if(scans_decoy.find(scan) == scans_decoy.end())
	    {
	      fullset_max.add_psm(fullset[i]);
	      scans_decoy.insert(scan);
	    }
	}
    }
  fullset_max.calc_factor();
  getOverFDR(fullset_max,net, qvals[4]);

  f1 << "scan" << "\t" << "charge" << "\t" << "spectrum neutral mass" << "\t" << "peptide mass" << "\t" << "q-ranker score" << "\t" << "q-ranker q-value" << "\t" << "sequence" << endl;
  for(unsigned int i = 0; i < fullset_max.size(); i++)
    {
      if(fullset_max[i].label == 1)
	{
	  int psmind = fullset_max[i].psmind;
	  int num_pep = d.psmind2num_pep(psmind);
	  
	  //write scan
	  f1 << d.psmind2scan(psmind) << "\t";
	  //write charges
	  int *charges = d.psmind2charges(psmind);
	  f1 << charges[0];
	  for(int k = 1; k < num_pep; k++)
	    {
	      f1 << ",";
	      int ch = charges[k];
	      f1 << ch; 
	    }
	  f1 << "\t";
	  //write neutral_mass
	  double* neutral_mass = d.psmind2neutral_mass(psmind);
	  f1 << neutral_mass[0];
	  for(int k = 1; k < num_pep; k++)
	    {
	      f1 << ",";
	      double m = neutral_mass[k];
	      f1 << m; 
	    }
	  f1 << "\t";
	  
	  //write peptide mass
	  double* peptide_mass = d.psmind2peptide_mass(psmind);
	  f1 << peptide_mass[0];
	  for(int k = 1; k < num_pep; k++)
	    {
	      f1 << ",";
	      double m = peptide_mass[k];
	      f1 << m; 
	    }
	  f1 << "\t";
	  f1 << fullset_max[i].score << "\t" << fullset_max[i].q << "\t";
	  
	  //write peptide sequence
	  int *pepinds = d.psmind2pepinds(psmind);
	  f1 << d.ind2pep(pepinds[0]);
	  for(int k = 1; k < num_pep; k++)
	    {
	      f1 << ",";
	      f1 << d.ind2pep(pepinds[k]); 
	    }
	  f1 << endl;
	}
    }
  f1.close();
}



void QRanker :: write_max_nets(string filename, NeuralNet* max_net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  
  for(int count = 0; count < num_qvals; count++)
    {
      net = max_net[count];
      int r = getOverFDR(testset,net, qvals[count]);
      int r1 = getOverFDR(trainset,net, qvals[count]);
      f1 << qvals[count] << " " << r1 << " " << r << "\n";
      
      double qn;
      if(qvals[count] < 0.01)
	qn = 0.0012;
      else
	qn = 0.005;
      r = getOverFDR(testset,net, qvals[count]+qn);
      r1 = getOverFDR(trainset,net, qvals[count]+qn);
      f1 << qvals[count]+qn << " " << r1 << " " << r << "\n";
      
    }
  f1.close();
}




void QRanker :: write_unique_peptides(string filename, NeuralNet* max_net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  
  for(int count = 0; count < num_qvals; count++)
    {
      net = max_net[count];
      int r = getOverFDR(testset,net, qvals[count]);
      set<int> peps;
      int cn = 0;
      for(unsigned int i = 0; i < testset.size();i++)
	{
	  if(testset[i].label == 1)
	    {
	      cn++;
	      int pepind = d.psmind2pepind(testset[i].psmind);
	      peps.insert(pepind);
	    }
	  if(cn > r) break;
	}
      int num_tst = 0;
      for(set<int>::iterator it = peps.begin(); it != peps.end(); it++)
	num_tst++;
      peps.clear();

      int r1 = getOverFDR(trainset,net, qvals[count]);
      cn = 0;
      for(unsigned int i = 0; i < trainset.size();i++)
	{
	  if(trainset[i].label == 1)
	    {
	      cn++;
	      int pepind = d.psmind2pepind(trainset[i].psmind);
	      peps.insert(pepind);
	    }
	  if(cn > r1) break;
	}
      int num_trn = 0;
      for(set<int>::iterator it = peps.begin(); it != peps.end(); it++)
	num_trn++;
      peps.clear();
      f1 << qvals[count] << " " << num_trn << " " << num_tst << "\n";
    }

  f1.close();
}


void QRanker :: write_num_psm_per_spectrum(PSMScores &fullset, int r)
{
  map<int,set<int> > scan_to_pepinds;
  int cn = 0;
  for(unsigned int i = 0; i < fullset.size();i++)
    {
      if(fullset[i].label == 1)
	{
	  cn++;
	  int psmind = fullset[i].psmind;
	  int *pepinds = d.psmind2pepinds(psmind);
	  int num_pep = d.psmind2num_pep(psmind);
	  for(int k = 0; k < num_pep; k++)
	    {
	      int pepind = pepinds[k];
	      int scan = d.psmind2scan(psmind);
	      if(scan_to_pepinds.find(scan) == scan_to_pepinds.end())
		{
		  set<int> peps;
		  scan_to_pepinds[scan] = peps;
		}
	      (scan_to_pepinds[scan]).insert(pepind);
	    }
	}
      if(cn > r) break;
    }
    
  vector<int> counts;

  for(map<int,set<int> >::iterator it = scan_to_pepinds.begin();
      it != scan_to_pepinds.end(); it++)
    {
      unsigned int cnt =  (it->second).size();
      while (counts.size() < cnt+1) {
        counts.push_back(0);
      }
      counts[cnt]++;
    }
  for(unsigned int i = 0; i < counts.size();i++)
    cout << i << " " << counts[i] << endl;

}



/*********************** training net functions*******************************************/
void QRanker :: train_net_hinge(PSMScores &set, int interval)
{
  double *r;
  int label;
  double *gc = new double[1];
  for(unsigned int i = 0; i < set.size(); i++)
    { 
      unsigned int ind;
      ind = rand()%(interval);
      //pass both through the net
      r = net.fprop(d.psmind2features(set[ind].psmind));
      label = d.psmind2label(set[ind].psmind);
      if(label*r[0]<1)
	{
	  net.clear_gradients();
	  gc[0] = -1.0*label;
	  net.bprop(gc);
	  net.update(mu,weightDecay);
	}
    }
  delete[] gc;
}

void QRanker :: train_net_ranking(PSMScores &set, int interval)
{
  double *r1;
  double *r2;
  double diff = 0;
  int label = 1;
  double *gc = new double[1];

  for(unsigned int i = 0; i < set.size(); i++)
    { 
      int ind1, ind2;
      int label_flag = 1;
      //get the first index
      if(interval == 0)
	ind1 = 0;
      else
	ind1 = rand()%(interval);
      if(ind1>set.size()-1) continue;
      if(set[ind1].label == 1)
	label_flag = -1;
      else
	label_flag = 1;
      
      int cn = 0;
      while(1)
	{
	  ind2 = rand()%(interval);
	  if(ind2>set.size()-1) continue;
	  if(set[ind2].label == label_flag) break;
	  if(cn > 1000)
	    {
	      ind2 = rand()%set.size();
	      break;
	    }
	  cn++;
	}
      
      //pass both through the net
      r1 = nets[0].fprop(d.psmind2features(set[ind1].psmind));
      r2 = nets[1].fprop(d.psmind2features(set[ind2].psmind));
      diff = r1[0]-r2[0];
            
      label=0;
      if(  set[ind1].label==1 && set[ind2].label==-1)
	label=1;
      if( set[ind1].label==-1 && set[ind2].label==1)
	    label=-1;
      
      if(label != 0)
	{
	  if(label*diff<1)
	    {
	      net.clear_gradients();
	      gc[0] = -1.0*label;
	      nets[0].bprop(gc);
	      gc[0] = 1.0*label;
	      nets[1].bprop(gc);
	      net.update(mu,weightDecay);
	      
	    }
	}
    }
  delete[] gc;
}

void QRanker :: train_many_general_nets()
{
  interval = trainset.size();
  for(int i=0;i<switch_iter;i++) {
    train_net_ranking(trainset, interval);
    
       
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
    if((i % 1) == 0)
      {
	cerr << "Iteration " << i << " : \n";
	cerr << "trainset: ";
	getMultiFDR(trainset,net,qvals);
	printNetResults(overFDRmulti);
	cerr << "\n";
	cerr << "testset: ";
	getMultiFDR(testset,net,qvals);
	printNetResults(overFDRmulti);
	cerr << "\n";
      }
  }

}

void QRanker :: train_many_target_nets()
{

  int  thr_count = num_qvals-1;
  while (thr_count > 0)
    {
      net.copy(max_net_gen[thr_count]);
          
      cerr << "training thresh " << thr_count  << "\n";
      //interval = getOverFDR(trainset, net, qvals[thr_count]);
      interval = max_overFDR[thr_count];
      for(int i=switch_iter;i<niter;i++) {
		
	//sorts the examples in the training set according to the current net scores
	getMultiFDR(trainset,net,qvals);
	train_net_ranking(trainset, interval);
			
	for(int count = 0; count < num_qvals;count++)
	  {
	    if(overFDRmulti[count] > max_overFDR[count])
	      {
		max_overFDR[count] = overFDRmulti[count];
		max_net_targ[count] = net;
	      }
	  }

	if((i % 1) == 0)
	  {
	    cerr << "Iteration " << i << " : \n";
	    getMultiFDR(testset,net,qvals);
	    cerr << "testset: ";
	    printNetResults(overFDRmulti);
	    cerr << "\n";
	  }

      }
      thr_count -= 3;
    }
}



void QRanker::train_many_nets()
{

  cerr << "Before iterating\n";
  cerr << "trainset: ";
  getMultiFDR(trainset,net,qvals);
  printNetResults(overFDRmulti);
  cerr << "\n";
  for(int count = 0; count < num_qvals;count++)
    if(overFDRmulti[count] > max_overFDR[count])
      max_overFDR[count] = overFDRmulti[count];
  cerr << "testset: ";
  getMultiFDR(testset,net,qvals);
  printNetResults(overFDRmulti);
  cerr << "\n";

  train_many_general_nets();
  
  //copy the general net into target nets;
  for(int count = 0; count < num_qvals; count++)
    max_net_targ[count] = max_net_gen[count];

  train_many_target_nets();  
 
  //choose the best net for the selectionfdr
  int max_fdr = 0;
  int fdr = 0;
  int ind = 0;
  for(unsigned int count = 0; count < qvals.size();count++)
    {
      fdr = getOverFDR(thresholdset, max_net_targ[count], selectionfdr);
      if(fdr > max_fdr)
	{
	  max_fdr = fdr;
	  ind = count;
	}
    }
  net.copy(max_net_targ[ind]);
  
}

int QRanker::run() {

  //set the seed
  srand(seed);
  
  //run the parser
  if(!pars.run(fnames))
    return 0;
  pars.clear();
    
  cout << "reading data\n";
  
  string summary_fn = "summary.txt";
  string psm_fn = "psm.txt";
  
  d.load_psm_data_for_training(summary_fn, psm_fn);
  d.normalize_psms();
  PSMScores::fillFeaturesSplit(trainset, testset, d, 0.5);
  thresholdset = trainset;
  
  switch_iter =30;
  niter = 60;
  
  /*
  switch_iter=1;
  niter=2;
  */
  num_qvals = 14;
  qvals.resize(num_qvals,0.0);
  qvals1.resize(num_qvals,0.0);
  qvals2.resize(num_qvals,0.0);

  overFDRmulti.resize(num_qvals,0);
  ave_overFDR.resize(num_qvals,0);
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
  int lf = 0; //num_hu = 1;//5;
  if(num_hu == 1)
    lf = 1;
  //set whether there is bias in the linear units: 1 if yes, 0 otherwise
  int bs = 0;
  
  net.initialize(d.get_num_features(),num_hu,lf,bs);
  for(int count = 0; count < num_qvals; count++){
    max_net_gen[count] = net;
    max_net_targ[count] = net;
  }

  nets = new NeuralNet[2];
  nets[0].clone(net);
  nets[1].clone(net);

  trained_nets = new NeuralNet[2];
  train_many_nets();
  trained_nets[0] = net;

  trainset = testset;
  testset = thresholdset;
  thresholdset = trainset;

  net.make_random();
  for(int count = 0; count < num_qvals; count++){
    max_net_gen[count] = net;
    max_net_targ[count] = net;
  }
  overFDRmulti.assign(overFDRmulti.size(),0);
  max_overFDR.assign(max_overFDR.size(),0);
  ave_overFDR.assign(ave_overFDR.size(),0);
  
  train_many_nets();
  trained_nets[1] = net;

  //write out main the results of the target net
  ostringstream res;
  res << out_dir << "/qranker_output";
  res_prefix = res.str();
  
  cerr << "target net results: ";
  ostringstream s2;
  s2 << res_prefix;
  cout << s2.str() << endl;
  //write_results(s2.str(),net);
  write_results_pvalue(s2.str());
  
  cerr <<"clean up"<<endl;
  pars.clean_up(out_dir);
  cerr <<"Done!"<<endl;
  return 1;
}


int QRanker::set_command_line_options(int argc, char **argv)
{
  int arg = 1;

  no_delta_cn = false;

  while(arg < argc)
    {

      string sarg(argv[arg]);
      if (sarg == "--nodeltacn") {
        cerr << "delta cn turned off"<<endl;
        no_delta_cn = true;
      } else if (sarg == "--num_hu") {
        num_hu = atoi(argv[arg+1]);
        arg++;
        cerr << "num_hu set to:"<<num_hu<<endl;
      } else if (sarg == "--mu") {
        mu = atof(argv[arg+1]);
        arg++;
        cerr << "mu set to:"<<mu << endl;
      } else if (sarg == "--decay") {
        weightDecay = atof(argv[arg+1]);
        arg++;
        cerr << "decay set to:"<<weightDecay<<endl; 
      }
      else if (sarg == "--extra-files") {
        string fn = argv[arg+1];
	extra_fnames.push_back(fn);
	extra_fnames_exist = 1;
	arg++;
	cerr << "--extra-files "<<fn<<endl; 
      }
      else {
        fnames.push_back(argv[arg]);
      }
      arg++;
    }
  string out_dir = "crux-output";
  pars.set_output_dir(out_dir);
  
  return 1;
}


int QRanker::main(int argc, char **argv) {

  if(!set_command_line_options(argc, argv))
    return -1;
  string dir = pars.get_output_dir();

  set_input_dir(dir);
  set_output_dir(dir);
  if(!run())
    return -1;
  cerr <<"Done qranker"<<endl;  
  return 0;
}   



string QRanker::getName() {
  return "q-ranker";
}

string QRanker::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the q-ranker algorithm (marina's new q-ranker)";

}
