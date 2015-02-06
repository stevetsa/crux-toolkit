#include "QRanker.h"

#include "analyze_psms.h"
#include "carp.h"
#include "objects.h"
#include "XLinkMatch.h"

QRanker::QRanker() :  seed(1),
                      selectionfdr(0.02),
                      num_hu(4),
                      mu(0.01),
                      weightDecay(0.0000),
                      xlink_mass(0.000),bootstrap_iters(5)
{
}

QRanker::~QRanker()
{
}

double QRanker:: getObjectiveErrorFDR(
  PSMScores& set,
  NeuralNet& n,
  double fdr) {

  int interval = get_interval(set, n, fdr, 0, 1);
  cerr << "interval is:"<<interval<<endl;
  return getObjectiveError(set, n, interval);

}

double QRanker:: getObjectiveError(
  PSMScores& set, 
  NeuralNet&n, 
  int interval) {

  int n_examples = set.size();

  calcScores(set, n);
  if (interval != -1) {
    set.calcOverFDR(0.1);
    n_examples = interval;
  }

  double sum = 0;
  double count = 0;
  for (int idx1 = 0;idx1 < n_examples;idx1++) {
    if (set[idx1].label == 1) {
      for (int idx2 = 0; idx2 < n_examples; idx2++) {
        if (set[idx2].label == -1) {
          double score1 = set[idx1].score;
          double score2 = set[idx2].score;
          double temp = max(0.0, 1.0 - (score1 - score2));
          sum += temp;
          count++;
        }
      }
    }
  }
  if (count == 0) {
    return 0.0;
  } else {
    return sum / count;
  }

}

int QRanker :: getOverFDR(PSMScores &set, NeuralNet &n, double fdr)
{
  double *r;
  double* featVec;

  for(int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = n.fprop(featVec);
      set[i].score = r[0];

    }
  return set.calcOverFDR(fdr);
}

void QRanker :: getMultiFDR(PSMScores &set, NeuralNet &n, vector<double> &qvalues)
{
  double *r;
  double* featVec;
   
  for(int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = n.fprop(featVec);
      set[i].score = r[0];
    }
 
  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti);
}


void QRanker :: getMultiFDRXCorr(PSMScores &set, vector<double> &qvalues)
{
  double* featVec;
   
  for(int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      set[i].score = featVec[0];
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


void QRanker :: write_results(string filename, NeuralNet &net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  trainset.clear();
  testset.clear();
  thresholdset.clear();
  PSMScores::fillFeaturesFull(fullset, d);
  getOverFDR(fullset,net, qvals[4]);
  
  int rank = 0;
  double last_score = -5000;
  //cerr << "Loading data for reporting results"<<endl;
  d.load_psm_data_for_reporting_results();
  //cerr << "Writing results" <<endl;
  f1 << "psm ind" << "\t" << 
        "q-value" << "\t" <<
        "score" << "\t" << 
        "scan" << "\t" << 
        "charge" << "\t" << 
        "label" << "\t" << 
        "rank" << "\t" << 
        "sequence" << endl;
  for(int i = 0; i < fullset.size(); i++)
    {
      //cerr << "Wrinting "<<i<<endl;
      int psmind = fullset[i].psmind;
      int scan = d.psmind2scan(psmind);
      int charge = d.psmind2charge(psmind);  
      double score = fullset[i].score;
      int label = fullset[i].label;
  
      if (score != last_score) {
        rank = (i+1);
        last_score = score;
      }

      ostringstream oss;
      if (label == 1) {
        oss << d.psmind2peptide1(psmind);
        if(string(d.psmind2peptide2(psmind)) != "_") {
          //cerr << "peptide2:"<<d.psmind2peptide2(psmind) <<" length:"<<d.psmind2peptide2(psmind).length() << endl;
          oss << ", " << d.psmind2peptide2(psmind);
        }
        if(string(d.psmind2loc(psmind)) != "") {
          oss << " " << d.psmind2loc(psmind);    
        }
     } else {
       oss << "null";
     }
  
      f1 << psmind << "\t" << fullset[i].q << "\t" << score << "\t" << scan << "\t" << charge << "\t" << label << "\t" << rank << "\t" << oss.str() << endl;
      

    }
  f1.close();
}

void QRanker :: write_results(string filename, PSMScores& set, bool decoy)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  
  cerr << "Loading data for reporting results"<<endl;
  d.load_psm_data_for_reporting_results();

  //use the best rank/scan to calculate q-values
  cerr <<"calculating FDR"<<endl;
  set.calcOverFDRBH(selectionfdr);

  cerr <<"calculating PEP"<<endl;
//  computePEP(set);


  cerr << "Writing results for " << set.size() <<" matches" << endl;
  f1 << "psmind" << "\t"; 
  f1 << "q-ranker q-value" << "\t" <<
        "q-ranker score" << "\t" <<
        "q-ranker rank" << "\t" <<
        "q-ranker p-value" << "\t" <<
        "PEP" << "\t" <<
        "scan" << "\t" << 
        "charge" << "\t" << 
        "spectrum neutral mass" << "\t"
        "sequence" << "\t" << 
        "protein id" << "\t" <<
        "product type" << "\t" <<
        "xcorr score" << "\t" <<
        "xcorr 1" << "\t" <<
        "xcorr 2" ; 
  //if (decoy) {
  f1 << "\tlabel";
  //}
  f1 << endl;
  for(int i = 0; i < set.size(); i++)
    {
      //cerr << "Wrinting "<<i<<endl;
      int psmind = set[i].psmind;
      int scan = d.psmind2scan(psmind);
      int charge = d.psmind2charge(psmind);
      int rank = set[i].rank;
      double spectrum_neutral_mass = d.psmind2neutral_mass(psmind);

  
      double score = set[i].score;
      int label = set[i].label;
      double qvalue = set[i].q;
      double pvalue = set[i].p;
      double pep = set[i].PEP;
      string proteins = d.psmind2protein1(psmind);
      double xcorr = d.psmind2xcorr(psmind);
      double xcorr1 = d.psmind2xcorr1(psmind);
      double xcorr2 = d.psmind2xcorr2(psmind);
      ostringstream oss;
      if (label == 1 || decoy) {
        oss << d.psmind2peptide1(psmind);
        if(string(d.psmind2peptide2(psmind)) != "_") {
          //cerr << "peptide2:"<<d.psmind2peptide2(psmind) <<" length:"<<d.psmind2peptide2(psmind).length() << endl;
          oss << ", " << d.psmind2peptide2(psmind);
        }
        if(string(d.psmind2loc(psmind)) != "") {
          oss << " " << d.psmind2loc(psmind);    
        }
  
        string sequence = oss.str();
        string product_type = "";

        product_type = XLinkMatch::getCandidateTypeString(d.psmind2product_type(psmind));

        f1 << psmind   << "\t" <<
              qvalue   << "\t" << 
              score    << "\t" <<
              rank     << "\t" <<
              pvalue   << "\t" <<
              pep      << "\t" <<
              scan     << "\t" << 
              charge   << "\t" << 
              spectrum_neutral_mass << "\t" << 
              sequence << "\t" << 
              proteins << "\t" << 
              product_type << "\t" <<
              xcorr << "\t" <<
              xcorr1 << "\t" <<
              xcorr2 ;
        //if (decoy) {
          f1 << "\t" << label;
        //}
        f1 << endl;
      }

    }
  f1.close();
}






void QRanker :: write_results_max(string filename, NeuralNet &net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  trainset.clear();
  testset.clear();
  thresholdset.clear();
  PSMScores::fillFeaturesFull(fullset, d);
  getOverFDR(fullset,net, qvals[4]);
  
  d.load_psm_data_for_reporting_results();

  set<int> scans_target;
  set<int> scans_decoy;
  
  for(int i = 0; i < fullset.size(); i++)
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

  f1 << "scan" << "\t" << "charge" << "\t" << "q-ranker score" << "\t" << "q-ranker q-value" << "\t" << "peptide sequences (loc)" <<"\t" << "protein loc1" << "\t" << "protein loc2" << endl;
  for(int i = 0; i < fullset_max.size(); i++)
    {
      if(fullset_max[i].label == 1)
	{
	  int psmind = fullset_max[i].psmind;
	  //write scan
	  f1 << d.psmind2scan(psmind) << "\t";
	  //write charges
	  f1 << d.psmind2charge(psmind) << "\t";
	  f1 << fullset_max[i].score << "\t" << fullset_max[i].q << "\t";
	  //write peptide sequence
	  f1 << d.psmind2peptide1(psmind);
	  if(d.psmind2peptide2(psmind) != "") 
	    f1 << "," << d.psmind2peptide2(psmind);
	  if(d.psmind2loc(psmind) != "")
	    f1 << " " << d.psmind2loc(psmind);
	  f1 << "\t";
	  f1 << d.psmind2protein1(psmind) << "\t" << d.psmind2protein2(psmind) << endl;
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
      for(int i = 0; i < testset.size();i++)
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
      for(int i = 0; i < trainset.size();i++)
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


void QRanker :: write_num_psm_per_spectrum(NeuralNet* max_net)
{
  //write out the results of the general net
  //ostringstream s1;
  //s1 << filename << ".txt";
  //ofstream f1(s1.str().c_str());
  
  int count = 5; 
  net = max_net[count];
  int r = getOverFDR(trainset,net, qvals[count]);
      
  map<int,set<int> > scan_to_pepinds;
  int cn = 0;
  for(int i = 0; i < trainset.size();i++)
    {
      if(trainset[i].label == 1)
	{
	  cn++;
	  int pepind = d.psmind2pepind(trainset[i].psmind);
	  int scan = d.psmind2scan(trainset[i].psmind);
	  if(scan_to_pepinds.find(scan) == scan_to_pepinds.end())
	    {
	      set<int> peps;
	      scan_to_pepinds[scan] = peps;
	    }
	  (scan_to_pepinds[scan]).insert(pepind);
	  
	}
      if(cn > r) break;
    }
    
  vector<int> counts;
  counts.resize(11,0);

  for(map<int,set<int> >::iterator it = scan_to_pepinds.begin();
      it != scan_to_pepinds.end(); it++)
    {
      int cnt =  (it->second).size();
      counts[cnt]++;
    }
  for(unsigned int i = 0; i < counts.size();i++)
    cout << i << " " << counts[i] << endl;

  //f1.close();
}



/*********************** training net functions*******************************************/
void QRanker :: train_net_sigmoid(PSMScores &set, int interval)
{
  double *r;
  int label;
  double *gc = new double[1];
  for(int i = 0; i < set.size(); i++)
    { 
      unsigned int ind;
      ind = rand()%(interval);
      //pass both through the net
      r = net.fprop(d.psmind2features(set[ind].psmind));
      label = d.psmind2label(set[ind].psmind);
      double a = exp(label*r[0]);
      net.clear_gradients();
      gc[0] = -a/((1+a)*(1+a))*label;
      net.bprop(gc);
      net.update(mu,weightDecay);

    }
  delete[] gc;
}


void QRanker :: train_net_hinge(PSMScores &set, int interval)
{
  double *r;
  int label;
  double *gc = new double[1];
  for(int i = 0; i < set.size(); i++)
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

  for(int i = 0; i < set.size(); i++)
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



void QRanker :: train_net_hybrid(PSMScores &set, int interval)
{
  double *r;
  int label;
  double *gc = new double[1];
  for(int i = 0; i < set.size(); i++)
    { 
      unsigned int ind;
      ind = rand()%(interval);
      //pass both through the net
      r = net.fprop(d.psmind2features(set[ind].psmind));
      label = d.psmind2label(set[ind].psmind);
      if(label == 1)
	{
	  double a = exp(label*r[0]);
	  net.clear_gradients();
	  gc[0] = -a/((1+a)*(1+a))*label;
	  net.bprop(gc);
	  net.update(mu,weightDecay);
	}
      else
	{
	  if(label*r[0]<1)
	    {
	      net.clear_gradients();
	      gc[0] = -1.0*label;
	      net.bprop(gc);
	      net.update(mu,weightDecay);
	    }
	}
    }
  delete[] gc;
}


void QRanker :: train_xcorr(PSMScores& set) {

  double* gc = new double[1];

  for (int idx = 0; idx < switch_iter ;idx++) {
    bool converged = true;
    getOverFDR(set, net, 1.0);
    double *r1;
    double *r2;
    for (int idx = 0; idx < set.size()-1; idx++) {
      if (d.psmind2xcorr(set[idx].psmind) < d.psmind2xcorr(set[idx+1].psmind)) {
        double xcorr1 = d.psmind2xcorr(set[idx].psmind);
        double xcorr2 = d.psmind2xcorr(set[idx+1].psmind);
        r1 = nets[0].fprop(d.psmind2features(set[idx].psmind));
        r2 = nets[1].fprop(d.psmind2features(set[idx+1].psmind));
        //cerr << idx << " " << r1[0] << " " << r2[0] << " " << d.psmind2xcorr(set[idx].psmind) << " " << d.psmind2xcorr(set[idx+1].psmind)<<endl; 
        double xdiff = xcorr2-xcorr1;
        if (xdiff > 0.01 && r2[0] < r1[0] ) {
          net.clear_gradients();
          gc[0] = 1;
          nets[0].bprop(gc);
          gc[0] = -1;
          nets[1].bprop(gc);
          net.update(mu, weightDecay);
          r1 = nets[0].fprop(d.psmind2features(set[idx].psmind));
          r2 = nets[1].fprop(d.psmind2features(set[idx+1].psmind));
          //cerr << "now "<<r1[0] << " " << r2[0] <<endl; 
          
          //cerr << "======================="<<endl;

          converged = false;
        }
      }
    }
    if (converged) break;
  }
  delete []gc;
}

void QRanker :: train_many_general_nets()
{
  interval = trainset.size();
  for(int i=0;i<switch_iter;i++) {
    train_net_ranking(trainset, interval);
    //train_net_hybrid(trainset, interval); 
           
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
    if (epoch_fout_ != NULL) {
/*
    (*epoch_fout_) << i << "\t" 
                   << getObjectiveError(trainset, net) << "\t"
                   << getObjectiveError(thresholdset, net) << "\t"
                   << getObjectiveError(testset, net) << "\t"
                   << getOverFDR(trainset, net, 0.01) << "\t" 
                   << getOverFDR(thresholdset, net, 0.01) << "\t" 
                   << getOverFDR(testset, net, 0.01) << endl;
*/
    }
    if((i % 10) == 0)
      {
    /*    
	cerr << "Iteration " << i << " : \n";
	cerr << "trainset: ";
	getMultiFDR(trainset,net,qvals);
	printNetResults(overFDRmulti);
	cerr << "\n";
	cerr << "testset: ";
	getMultiFDR(testset,net,qvals);
	printNetResults(overFDRmulti);
	cerr << "\n";
    */  
      }
  }

}


int QRanker :: get_interval(PSMScores& set, NeuralNet& net, double qval, int min_targets, int min_decoys) {

  int current_ntargets = 0;
  int current_ndecoys = 0;

  int interval = 0;
  int ntargets = 0;
  int ndecoys = 0;
  int idx = 0;

  getOverFDR(set, net, qval);

  for (int idx = 0; idx < set.size();idx++) {
    int label = set[idx].label;
    if (label == 1) {
      current_ntargets++;
    } else {
      current_ndecoys++;
    }

    if (current_ntargets > 0) {
      double fdr = (double)current_ndecoys / (double)current_ntargets * set.factor; 
      if (fdr <=qval) {
        interval = idx;
        ntargets = current_ntargets;
        ndecoys = current_ndecoys;
      }
    }
  }
/*
  cerr <<"fdr:"<<qval<<
         " interval:" << interval << 
         " ntargets:"<<ntargets<<
         " ndecoys:"<<ndecoys<<
         " size:"<<set.size() << endl;
*/
  while ((ntargets < min_targets || ndecoys < min_decoys) & interval < set.size()-1) {
    interval++;
    if (set[interval].label == 1) {
      ntargets++;
    } else {
      ndecoys++;
    }
  }

//  cerr <<"final interval:"<<interval<<" ntargets:"<<ntargets<<" ndecoys:"<<ndecoys<<endl;
  return interval;

}


void QRanker :: train_many_target_nets()
{

  //PSMScores trainset_max;
  //PSMScores testset_max;
  int  thr_count = num_qvals-1;
  while (thr_count > 0)
    {
      net.copy(max_net_gen[thr_count]);
          
      //calcScores(trainset, net);
      //trainset.getMaxPerScan(d, trainset_max);

      //cerr << "training thresh " << thr_count  << "\n";
      //interval = getOverFDR(trainset, net, qvals[thr_count]);
      //cerr << interval << " "<< get_interval(trainset, net, qvals[thr_count]);
      interval = get_interval(trainset, net, qvals[thr_count]);
      //cerr << "interval: "<<interval<<endl;
      //interval = max_overFDR[thr_count];
      if(interval < 10)
	interval = (int)trainset.size()/50;
      //interval = trainset.size()/thr_count;
      for(int i=switch_iter;i<niter;i++) {
		
	//sorts the examples in the training set according to the current net scores
	getMultiFDR(trainset,net,qvals);
	train_net_ranking(trainset, interval);
	//train_net_hybrid(trainset, interval);
        //calcScores(trainset, net);
        //trainset.getMaxPerScan(d, trainset_max);
        int new_interval = get_interval(trainset, net, qvals[thr_count]);
        if (new_interval > interval) {
          //cerr << "interval_old:"<<interval<<" new:"<<new_interval<<endl;
          interval = new_interval;
          
        }
	getMultiFDR(trainset,net,qvals);			
	for(int count = 0; count < num_qvals;count++)
	  {
	    if(overFDRmulti[count] > max_overFDR[count])
	      {
		max_overFDR[count] = overFDRmulti[count];
		max_net_targ[count] = net;
	      }
	  }

	if((i % 10) == 0)
	  {
/*            
	    cerr << "Iteration " << i << " : \n";
            cerr << "trainset: ";
            getMultiFDR(trainset_max,net,qvals);
            printNetResults(overFDRmulti);
            cerr << "\n";
            calcScores(testset,net);
            testset.getMaxPerScan(d, testset_max);
	    getMultiFDR(testset_max,net,qvals);
	    cerr << "testset: ";
	    printNetResults(overFDRmulti);
	    cerr << "\n";
*/            
	  }

      }
      thr_count -= 3;
    }
}



void QRanker::train_many_nets()
{

  int max_fdr = 0;
  NeuralNet best_net;

  //for (int idx = 0; idx < 10; idx++) {

   
    num_qvals = 14;
    qvals.resize(num_qvals,0.0);
    qvals1.resize(num_qvals,0.0);
    qvals2.resize(num_qvals,0.0);
  
    overFDRmulti.resize(num_qvals,0);
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
    int bs = 0;
    
    net.initialize(d.get_num_features(),num_hu,lf,bs);

  /*
    for(int count = 0; count < num_qvals; count++) {
      cerr <<"max_overFDR["<<count<<"]="<<max_overFDR[count]<<endl;
    }
  */
  
    nets = new NeuralNet[2];
      nets[0].clone(net);
    nets[1].clone(net);

    //train_xcorr(trainset);

    for(int count = 0; count < num_qvals; count++){
      max_net_gen[count] = net;
    }
    for(int count = 0; count < num_qvals; count++){
      max_net_targ[count] = net;
    }

/*
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
*/
    //switch_iter =0;niter=20;

    train_many_general_nets();
    
    //copy the general net into target nets;
    for(int count = 0; count < num_qvals; count++)
      max_net_targ[count] = max_net_gen[count];
  
    train_many_target_nets();  
   
        
    //choose the best net for the selectionfdr
    //int max_fdr = 0;
    int fdr = 0;
    int ind = -1;
    for(unsigned int count = 0; count < qvals.size();count++)
        {
        fdr = getOverFDR(thresholdset, max_net_targ[count], selectionfdr);
        //cerr << "t:"<< count<<" "<<fdr<<" "<<max_fdr<<endl;
        if(fdr > max_fdr)
          {
            max_fdr = fdr;
            ind = count;
          }
      }

    //cerr <<"max count:"<<max_fdr<<endl;
  //cerr << idx <<" "<<max_fdr<<endl;
  if (ind != -1) {
    
    best_net = max_net_targ[ind];
//  }
  //print out results, just to see
  /*
  for(unsigned int count = 0; count < qvals.size();count++)
    cout << qvals[count] << " " <<  getOverFDR(trainset, net, qvals[count]) << " " << getOverFDR(testset, net, qvals[count]) << endl;
  */

    delete [] max_net_gen;
    delete [] max_net_targ;
    delete [] nets;
  }

  net = best_net;
}


void QRanker::calcScores(PSMScores& set, NeuralNet& net) {

  double *r;
  double* featVec;

  for(int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = net.fprop(featVec);
      set[i].score = r[0];

    }
}

void QRanker::calcRanks(PSMScores& set, NeuralNet& net) {
  getOverFDR(set, net, 0);

  vector<PSMScoreHolder>::iterator iter;

  int current_rank = 1;
  int count = 1;
  double current_score = set.begin()->score;

  for (iter = set.begin(); iter != set.end(); ++iter) {
    if (iter->score != current_score) {
      current_score = iter->score;
      current_rank = count;
    }

    //cerr << "psmind:"<<iter->psmind<<" label:"<<iter->label<<" rank:"<<current_rank<<endl;

    //sum them
    iter->rank += current_rank;
    //cerr << "psmind:"<<iter->psmind<<" label:"<<iter->label<<" rank:"<<current_rank<<" sum:"<<iter->rank<<endl;
    count++;
  }
}

void QRanker::avgRanks(PSMScores& set, int n) {
  
  vector<PSMScoreHolder>::iterator iter;

  for (iter = set.begin(); iter != set.end(); ++iter) {
    iter->rank = iter->rank / (double)n;
  }



}
/*
void QRanker::getMinRank(PSMScores&in, PSMScores& out) {


}
*/
void QRanker::selectHyperParameters() {

  PSMScores previous_trainset = trainset;
  PSMScores previous_testset = testset;
  PSMScores previous_thresholdset = thresholdset;

  PSMScores cv_trainset;
  PSMScores cv_testset;
  PSMScores cv_fullset;
  PSMScores max_per_scan;

  PSMScores::fillFeaturesSplitScan(trainset, d, cv_trainset, cv_testset);

  
  vector<int> general_iters;
  general_iters.push_back(5);
  general_iters.push_back(10);
  general_iters.push_back(20);
  general_iters.push_back(40);

  vector<int> ranking_iters;
  ranking_iters.push_back(5);
  ranking_iters.push_back(10);
  ranking_iters.push_back(20);
  ranking_iters.push_back(40);

  vector<double> mus;
  mus.push_back(0.0005);
  mus.push_back(0.001);
  mus.push_back(0.005);
  mus.push_back(0.01);
//  mus.push_back(0.05);
//  mus.push_back(0.1);

  vector<double> wds;
  wds.push_back(0.000);
  wds.push_back(1e-7);
  wds.push_back(1e-6);
//  wds.push_back(1e-5);

  vector<int> num_hus;
  num_hus.push_back(1);
  num_hus.push_back(2);
  num_hus.push_back(3);
//num_hus.push_back(5);
  //num_hus.push_back(7);
 
  int best_general_iter_idx = -1;
  int best_ranking_iter_idx = -1;
  int best_mu_idx=-1;
  int best_wd_idx=-1;
  int best_num_hu_idx=-1;

  int best_score = -1;
  double best_obj = -1;

  for (size_t general_iter_idx = 0; general_iter_idx < general_iters.size();general_iter_idx++) {

    for (size_t ranking_iter_idx = 0; ranking_iter_idx < ranking_iters.size(); ranking_iter_idx++) {

      for (size_t mu_idx = 0; mu_idx < mus.size(); mu_idx++) {

        for (size_t wd_idx = 0; wd_idx < wds.size(); wd_idx++) {

          for (size_t num_hu_idx = 0; num_hu_idx < num_hus.size(); num_hu_idx++) {
            this->switch_iter = general_iters[general_iter_idx];
            this->niter = general_iters[general_iter_idx]+ranking_iters[ranking_iter_idx];
            if (this->switch_iter == 0 && this->niter == 0) {
              continue;
            }

            this->mu = mus[mu_idx];
            this->weightDecay = wds[wd_idx];
            this->num_hu = num_hus[num_hu_idx];
            cv_fullset.clear();
            trainset = cv_trainset;
            thresholdset = cv_trainset;
            testset = cv_testset;

            train_many_nets();

            calcScores(testset, net);
            testset.getMaxPerScan(d, max_per_scan);

            int score = max_per_scan.calcOverFDR(selectionfdr);
            double obj = getObjectiveError(max_per_scan, net);

//            int score = trainset.calcOverFDR(selectionfdr);
//            double obj = getObjectiveError(trainset, net);

            //swap trainset and testset.
            trainset = cv_testset;
            thresholdset = cv_testset;
            testset = cv_trainset;
  
            train_many_nets();
  
            calcScores(testset, net);
            testset.getMaxPerScan(d, max_per_scan);
            
            score+= max_per_scan.calcOverFDR(selectionfdr);
            obj += getObjectiveError(max_per_scan, net);
    

//            score+= testset.calcOverFDR(selectionfdr);
//            obj += getObjectiveError(testset, net);
            obj = obj / 2.0;


            if (score > best_score || (score == best_score && best_obj > obj)) {
              cout <<"switch:"<<switch_iter<<
                   " niter:"<<niter<<
                   " hu:"<< num_hu << 
                   " mu:"<<mu<<
                   " wd:"<<weightDecay<<
                   " score:"<<score<<
                   " obj: " << obj << endl;
              best_score = score;
              best_obj = obj;
              best_general_iter_idx = general_iter_idx;
              best_ranking_iter_idx = ranking_iter_idx;
              best_num_hu_idx = num_hu_idx;
              best_wd_idx = wd_idx;
              best_mu_idx = mu_idx;
            }
          }
        }
      }
    }
  }

  if (best_general_iter_idx != -1) {
    if ((best_general_iter_idx == 0) || (best_general_iter_idx == general_iters.size()-1)) {
      cerr <<"edge with general iters! :"<<general_iters[best_general_iter_idx]<<endl;
    }
    this->switch_iter = general_iters[best_general_iter_idx];
    this->niter = this->switch_iter;
  }

  if (best_ranking_iter_idx != -1) {
    if (best_ranking_iter_idx == 0 || best_ranking_iter_idx == ranking_iters.size() -1) {
      cerr <<"edge with ranking iters! :"<<ranking_iters[best_ranking_iter_idx] << endl;
    }
    this->niter = this->switch_iter + ranking_iters[best_ranking_iter_idx];
  }

  if (best_mu_idx != -1) {
    if (best_mu_idx == 0 || best_mu_idx == mus.size()-1) {
      cerr <<"edge with mu! :"<<mus[best_mu_idx]<<endl;
    }
    this->mu = mus[best_mu_idx];
  }
  if (best_wd_idx != -1) {
    if (best_wd_idx == 0 || best_wd_idx == wds.size()-1) {
      cerr << "egde with wd! :"<<wds[best_wd_idx] << endl;
    }
    this->weightDecay = wds[best_wd_idx];
  }

  if (best_num_hu_idx != -1) {
    if (best_num_hu_idx == 0 || best_num_hu_idx == num_hus.size()-1) {
      cerr << "edge with num_hus! :"<<num_hus[best_num_hu_idx]<<endl;
    }
    this->num_hu = num_hus[best_num_hu_idx];
  }

  cerr << "best switch:"<<switch_iter<<endl;
  cerr << "best iter:"<<niter<<endl;
  cerr << "best num_hu:"<<num_hu<<endl;;
  cerr << "best mu:"<<mu<<endl;
  cerr << "best wd:"<<weightDecay<<endl;
  cerr << "best score:"<<best_score<<endl;
  cerr << "best obj:"<<best_obj<<endl;

  trainset = previous_trainset;
  testset  = previous_testset;
  thresholdset = previous_thresholdset;



}


void QRanker::writeFeatures() {


  

  ofstream fout("features.txt");

  fout << "psmind" << "\t";
  fout << "label" << "\t";
  fout << "scan" << "\t";
  fout << "sequence" << "\t";
  fout << "xcorr score" << "\t";
  fout << "xcorr short" << "\t";
  fout << "xcorr long" << "\t";
  fout << "sp score" << "\t";
  fout << "log(sp rank)"<< "\t";
  fout << "frac"<<"\t";
  fout << "absMass"<<"\t";
  fout << "linear?"<<"\t";
  fout << "selfloop?"<<"\t";
  fout << "xlink?"<<"\t";
  fout << "dead-link?"<<"\t";
  fout << "length short"<<"\t";
  fout << "length long" << "\t";
  fout << "length sum" << "\t";
  fout << "Nenz" << "\t";
  fout << "nterm1" << "\t";
  fout << "cterm1" << "\t";
  fout << "cterm2" << "\t";
  fout << "cterm2" << "\t";
  fout << "nmissed" << "\t";
  fout << "obsmass" << "\t";
  fout << "ncmp"<<endl;

  d.load_psm_data_for_training();
  PSMScores scores;
  PSMScores::fillFeaturesFull(scores, d);
  
  for (int idx = 0; idx < scores.size();idx++) {
    int psmind = scores[idx].psmind;
    int scan = d.psmind2scan(psmind);
    //int charge = d.psmind2charge(psmind);
    int label = scores[idx].label;

    ostringstream oss;
    if (label == 1) {
      oss << d.psmind2peptide1(psmind);
      if(string(d.psmind2peptide2(psmind)) != "_") {
      //cerr << "peptide2:"<<d.psmind2peptide2(psmind) <<" length:"<<d.psmind2peptide2(psmind).length() << endl;
        oss << ", " << d.psmind2peptide2(psmind);
      }
      if(string(d.psmind2loc(psmind)) != "") {
        oss << " " << d.psmind2loc(psmind);    
      }
    } else {
      oss << "null";
    }
  
    double* x = d.psmind2features(psmind);

    fout << psmind << "\t";
    fout << label << "\t";
    fout << scan << "\t";
    fout << oss.str() << "\t";
    fout << x[0] << "\t"; // xcorr score
    fout << x[1] << "\t"; // xcorr short
    fout << x[2] << "\t"; // xcorr long
    fout << x[7] << "\t"; // sp score
    fout << x[8] << "\t"; // log(sp_rank)
    fout << x[9] << "\t"; // frac
    fout << x[10] << "\t"; // absMass
    fout << x[11] << "\t"; // linear?
    fout << x[12] << "\t"; // selfloop?
    fout << x[13] << "\t"; // xlink?
    fout << x[14] << "\t"; // dead-link?
    fout << x[15] << "\t"; // length short
    fout << x[16] << "\t"; // length long
    fout << x[17] << "\t"; // length sum
    fout << x[18] << "\t"; // nenz
    fout << x[19] << "\t"; // nterm1
    fout << x[20] << "\t"; // cterm1
    fout << x[21] << "\t"; // nterm2
    fout << x[22] << "\t"; // cterm2
    fout << x[23] << "\t"; // nmissed
    fout << x[24] << "\t"; // obsmass
    fout << x[25] << "\t"; // ncmp
    fout << x[26] << "\t"; // chargeX
    fout << x[27] << "\t"; // chargeX+1
    fout << x[28] << endl; // chargeX+2
  }
}




int QRanker::run( ) {
  srand(seed);
  cout << "reading data\n";
  
  writeFeatures();

  PSMScores max;
  PSMScores fullset_max;
  d.load_psm_data_for_training();
  d.normalize_psms();
  //PSMScores::fillFeaturesSplit(trainset, testset, d, 0.5);

  fullset.clear();

  PSMScores::fillFeaturesSplitScan(trainset, testset, d);
  thresholdset = trainset;
  fullset.clear();

  epoch_fout_ = NULL;
  selectHyperParameters();


  epoch_fout_ = new ofstream("qranker.1.obj.txt");
 
  train_many_nets();

  epoch_fout_ -> close();
  delete epoch_fout_;
  epoch_fout_ == NULL;
  double train1 = getObjectiveError(trainset, net);
  double test1 = getObjectiveError(testset, net);

  cerr << "train objective error:"<<getObjectiveError(trainset, net) << endl;
  cerr << "test objective error:"<<getObjectiveError(testset, net) << endl;

 

  calcScores(testset, net);
  calcRanks(testset, net);
  testset.getMaxPerScan(d, max);
  max.calcPValues();
  testset.calcPValues();

  // write out the testset results
  ostringstream res;
  res << out_dir << "/qranker.test";
  write_results(res.str(), testset, true);
  
  //fullset.add_psms(max);
  fullset.add_psms(testset);
  fullset_max.add_psms(max);
  //swap trainset and test set.
  thresholdset = testset;
  testset = trainset;
  trainset = thresholdset;

  selectHyperParameters();

  epoch_fout_ = new ofstream("qranker.2.obj.txt");

  train_many_nets();

  epoch_fout_ -> close();
  delete epoch_fout_;

  cerr << "train objective error:"<<getObjectiveError(trainset, net) << endl;
  cerr << "test objective error:"<<getObjectiveError(testset, net) << endl;

  double train2 = getObjectiveError(trainset, net);
  double test2 = getObjectiveError(testset, net);

  calcScores(testset, net);
  calcRanks(testset, net);
  testset.getMaxPerScan(d, max);
  max.calcPValues();
  testset.calcPValues();

  //write out the trainset results
  res.str("");
  res << out_dir << "/qranker.train";
  write_results(res.str(), testset, true);

  fullset.add_psms(testset);
  fullset_max.add_psms(max);
  fullset.calc_factor();
  fullset_max.calc_factor();

  //write out the fullset results
  res.str("");
  cerr << "writing all matches"<<endl;
  res << out_dir << "/qranker.all";
  write_results(res.str(), fullset, true);

  res.str("");
  res << out_dir << "/qranker.target";
  write_results(res.str(), fullset_max, true);

  int q = fullset_max.calcOverFDRBH(selectionfdr);
  
  ofstream fout("objective.txt");
  fout << "train1\ttest1\ttrain2\ttest2\tq"<<endl;
  fout << train1 << "\t" << test1 << "\t" << train2 << "\t" << test2 << "\t" << q << endl;
  fout.close();


  return 0;
}





int QRanker::set_command_line_options(int argc, char **argv)
{
  vector<string> fnames;
  string ms2fname;
  int ms2exists = 0;
  int arg = 1;

  switch_iter = 20;
  niter = 40;

  while(arg < argc)
    {
      string  str = argv[arg]; 
      size_t pos = str.find("=");
      if(pos != string::npos)
	{
          string tmp = str.substr(pos+1, str.size());
	  if(str.find("seed") != string::npos) {
	    cout << "found seed " << tmp << endl;
	    seed = atoi(tmp.c_str());
	  } else if (str.find("switch-iter") != string::npos) {
            cout << "found switch_iter "<< tmp << endl;
            switch_iter = atoi(tmp.c_str());
          } else if (str.find("niter") != string::npos) {
            cout << "found niter "<<tmp<<endl;
            niter = atoi(tmp.c_str());
          } else if (str.find("mu") != string::npos) {
            cout << "found mu "<<tmp<<endl;
            mu = atof(tmp.c_str());
          } else if (str.find("wd") != string::npos) {
            cout << "found wd "<<tmp<<endl;
            weightDecay = atof(tmp.c_str());
          } else if (str.find("num-hu") != string::npos) {
            cout << "found num hu" << tmp << endl;
            num_hu = atoi(tmp.c_str());
          } else if (str.find("xlink-mass") != string::npos) {
            cout << "found xlink-mass" << tmp << endl;
            xlink_mass = atof(tmp.c_str());
          } else if (str.find("ms2file") != string::npos) {
            cout << "found ms2file "<<tmp<<endl;
            ms2fname = tmp;
	    ms2exists = 1;
          } else if (str.find("quad") != string::npos) {
            cout << "found quad" << tmp << endl;
            int quad = atoi(tmp.c_str());
            pars.set_use_quadratic_features(quad);
          } else if (str.find("bootstrap") != string::npos) {
            cout << "found bootstrap "<<tmp<<endl;
            bootstrap_iters = atoi(tmp.c_str());
          }

	}
      else
	fnames.push_back(argv[arg]);
      arg++;
    }

  niter = max(switch_iter, niter);

  string out_dir = "crux-output";
  pars.set_output_dir(out_dir);
  if(ms2exists)
    {
      if(!pars.run_on_xlink(fnames, ms2fname, xlink_mass))
	return 0;
    }
  else
    {
      if(!pars.run_on_xlink(fnames))
	return 0;
    }
  writeFeatures();
  pars.clear();
  return 1;
}

/**
 * Uses the target and decoy scores to compute posterior error
 * probabilities. 
 */
void QRanker::computePEP(PSMScores& scores){
  carp(CARP_DEBUG, "Computing PEPs");
  vector<double> target_scores_vect;
  vector<double> decoy_scores_vect;

  // pull out the target and decoy scores
  for(int i = 0; i < scores.size(); i++){
    if( scores[i].label == 1 ){
      target_scores_vect.push_back(scores[i].score);
    } else { // == -1
      decoy_scores_vect.push_back(scores[i].score);
    }
  }

  int num_targets = target_scores_vect.size();
  int num_decoys = decoy_scores_vect.size();
  carp(CARP_DEBUG, "Found %d targets and %d decoys", num_targets, num_decoys); 

  // copy them to an array as required by the compute_PEP method
  double* target_scores = new double[num_targets];
  copy(target_scores_vect.begin(), target_scores_vect.end(), target_scores);
  double* decoy_scores = new double[num_decoys];
  copy(decoy_scores_vect.begin(), decoy_scores_vect.end(), decoy_scores);

  double* PEPs = compute_PEP(target_scores, num_targets, 
                             decoy_scores, num_decoys);

  // fill in the data set with the new scores for the targets
  int target_idx = 0;
  for(int full_idx = 0; full_idx < scores.size(); full_idx++){
    if( scores[full_idx].label == 1 ){
      scores[full_idx].PEP = PEPs[target_idx];
      target_idx++; 
    } // else, skip decoys
  }

  delete target_scores;
  delete decoy_scores;
  delete PEPs;
}


int QRanker::main(int argc, char **argv) {

  if(!set_command_line_options(argc, argv))
    return 0;
  string dir = pars.get_output_dir();

  set_input_dir(dir);
  set_output_dir(dir);
  run();
    
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
