#include "QRanker.h"

QRanker::QRanker() :  seed(0),selectionfdr(0.01),num_hu(5),mu(0.005),weightDecay(0.000)
{
}

QRanker::~QRanker()
{
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
 
  f1 << "psm ind" << "\t" << "q-value" << "\t" << "scan" << "\t" << "charge" << "\t" << "peptide" << endl;
  for(int i = 0; i < fullset.size(); i++)
    {
      int psmind = fullset[i].psmind;
      int pepind = d.psmind2pepind(psmind);
      f1 << psmind << "\t" << fullset[i].q << "\t" << d.psmind2scan(psmind) << "\t" << d.psmind2charge(psmind) << "\t" << d.ind2pep(pepind) << endl;
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
	      gc[0] = -1.0;
	      nets[0].bprop(gc);
	      gc[0] = 1.0;
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
  switch_iter =30;
  niter = 40;
   
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
  int lf = 0; num_hu = 5;
  if(num_hu == 1)
    lf = 1;
  //set whether there is bias in the linear units: 1 if yes, 0 otherwise
  int bs = 1;
  
  net.initialize(d.get_num_features(),num_hu,lf,bs);
  for(int count = 0; count < num_qvals; count++){
    max_net_gen[count] = net;
  }

  nets = new NeuralNet[2];
  nets[0].clone(net);
  nets[1].clone(net);

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
 
  /*
  //write out the results of the target net
  cerr << "target net results: ";
  ostringstream s2;
  s2 << res_prefix;
  cout << s2.str() << endl;
  //write_max_nets(s2.str(), max_net_targ);
  //write_unique_peptides(s2.str(), max_net_targ);

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
  net = max_net_targ[ind];

  //write_num_psm_per_spectrum("test.txt", max_net_targ);
  write_results(s2.str(),net);
  */
  delete [] max_net_gen;
  delete [] max_net_targ;
  delete [] nets;


}

int QRanker::run( ) {
  srand(seed);
  cout << "reading data\n";
  
  ostringstream res;
  res << out_dir << "/qranker_output";
  res_prefix = res.str();
    
  d.load_psm_data_for_training();
  d.normalize_psms();
  PSMScores::fillFeaturesSplit(trainset, testset, d, 0.5);
  thresholdset = trainset;
  train_many_nets();
  return 0;
}



int QRanker::set_command_line_options(int argc, char **argv)
{
  vector<string> fnames;
  int arg = 1;
  while(arg < argc)
    {
      fnames.push_back(argv[arg]);
      arg++;
    }
  string out_dir = "crux-output";
  pars.set_output_dir(out_dir);
  if(!pars.run_on_xlink(fnames))
    return 0;
  pars.clear();
  return 1;
}


int QRanker::main(int argc, char **argv) {

  if(!set_command_line_options(argc, argv))
    return 0;
  string dir = pars.get_output_dir();

  set_input_dir(dir);
  set_output_dir(dir);
  run();
    
  return 1;
}   



string QRanker::getName() {
  return "q-ranker";
}

string QRanker::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the q-ranker algorithm (marina's new q-ranker)";

}
