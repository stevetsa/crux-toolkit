#include "TabDelimParser.h"

/******************************/
/* tokens are
0. scan
1. charge
2. spectrum precursor mz
3. spectrum neutral mass   
4. peptide mass    
5. delta_cn        
6. sp score        
7. xcorr score     
8. xcorr rank      
9. matches/spectrum        
10. sequence        
11. cleavage type   
12. protein id      
13. flanking aa     
14. nzstates        
15. rtime max diff  
16. peptides/spectrum       
17. xcorr sum diff  
18. xcorr max diff
*/

const static int scan_col_idx = 0;
const static int charge_col_idx = 1;
const static int spectrum_mz_col_idx = 2;
const static int spectrum_mass_col_idx = 3;
const static int peptide_mass_col_idx = 4;
const static int sp_score_col_idx = 6;
const static int xcorr_score_col_idx = 7;
const static int matches_spectrum_col_idx = 9;
const static int sequence_col_idx = 10;

const static int max_charge=6;

TabDelimParser :: TabDelimParser() 
  : 	num_mixed_labels(0),     
	psmind_to_scan(0),
	psmind_to_charge(0),
	psmind_to_label(0), 
	psmind_to_num_pep(0),
	psmind_to_ofst(0),
	psmind_to_pepind(0),
	psmind_to_neutral_mass(0),
	psmind_to_peptide_mass(0),
	num_features(0),
	num_psm(0),
	num_pos_psm(0),
	num_neg_psm(0),
	num_pep(0),
	num_pep_in_all_psms(0),
	curr_ofst(0),
	psmind(0),
	x(0)
{
  
  //num_psm_features
  num_features = 19;
  
  //final_hits_per_spectrum
  fhps = 3;
  //decoy prefix
  decoy_prefix = "random_";
  //max peptide length to be considered
  max_len = 50;
  //min peptide length to be considered
   min_len = 7;
  
}

void TabDelimParser :: clear()
{
  delete[] x; x = (double*)0;
  delete[] psmind_to_scan; psmind_to_scan = (int*)0;
  delete[] psmind_to_charge; psmind_to_charge = (int*)0;
  delete[] psmind_to_label; psmind_to_label = (int*)0;
  delete[] psmind_to_num_pep; psmind_to_num_pep = (int*)0;
  delete[] psmind_to_ofst; psmind_to_ofst = (int*)0;
  delete[] psmind_to_pepind; psmind_to_pepind = (int*)0;
  delete[] psmind_to_neutral_mass; psmind_to_neutral_mass = (double*)0;
  delete[] psmind_to_peptide_mass; psmind_to_peptide_mass = (double*)0;
}



TabDelimParser :: ~TabDelimParser()
{
  clear();
}




void TabDelimParser :: get_tokens(string &line, vector<string>&tokens, string &delim)
{
  tokens.erase(tokens.begin(), tokens.end());
  string tmp = line;
  string tok;
  size_t pos = tmp.find(delim);
  while(pos != string::npos)
    {
      tok = tmp.substr(0,pos);
      tokens.push_back(tok);
      tmp = tmp.substr(pos+1, tmp.size());
      pos = tmp.find(delim);
    }
  //last token
  tokens.push_back(tmp);
}

void TabDelimParser :: first_pass(ifstream &fin)
{
  string line;
  vector<string> tokens;
  vector<string> peptides;
  string delim1 = "\t";
  string delim2 = ",";
  while(!fin.eof())
    {
      getline(fin,line);
      get_tokens(line, tokens, delim1);
  
      if(tokens.size() > 1)
	{
	  //get all peptides and fill pep_to_ind tables
	  string peps = tokens[sequence_col_idx];
	  get_tokens(peps,peptides,delim2);
	  
	  for(unsigned int i = 0; i < peptides.size(); i++)
	    {
	      string pep = peptides[i];
	      int pep_ind = -1;
	      //add peptide to pep_to_ind and ind_to_pep maps
	      if(pep_to_ind.find(pep) == pep_to_ind.end())
		{
		  pep_ind = num_pep;
		  //add a new peptide
		  pep_to_ind[pep] = pep_ind;
		  ind_to_pep[pep_ind] = pep;
		  num_pep++;
	  	}
	      else
		{
		  pep_ind = pep_to_ind[pep];
		  string p = ind_to_pep[pep_ind];
		  if(pep.compare(p) != 0)
		    cout << "warning : did not find peptide in ind_to_pep_table\n"; 
		}
	    }
	  num_pep_in_all_psms += peptides.size();
	  num_psm++;
	}
      
    }
}


void TabDelimParser :: allocate_feature_space()
{
  //space for feature vector
  x = new double[num_features];
  memset(x,0,sizeof(double)*num_features);

  //space for psminfo
  psmind_to_pepind = new int[num_pep_in_all_psms];
  memset(psmind_to_pepind,0,sizeof(int)*num_pep_in_all_psms);
  psmind_to_num_pep = new int[num_psm];
  memset(psmind_to_num_pep,0,sizeof(int)*num_psm);
  psmind_to_ofst = new int[num_psm];
  memset(psmind_to_ofst,0,sizeof(int)*num_psm);
  psmind_to_charge = new int[num_pep_in_all_psms];
  memset(psmind_to_charge,0,sizeof(int)*num_pep_in_all_psms);

  psmind_to_neutral_mass = new double[num_pep_in_all_psms];
  memset(psmind_to_neutral_mass,0,sizeof(double)*num_pep_in_all_psms);
  psmind_to_peptide_mass = new double[num_pep_in_all_psms];
  memset(psmind_to_peptide_mass,0,sizeof(double)*num_pep_in_all_psms);

  psmind_to_scan = new int[num_psm];
  memset(psmind_to_scan,0,sizeof(int)*num_psm);
  psmind_to_label = new int[num_psm];
  memset(psmind_to_label,0,sizeof(int)*num_psm);

}


void TabDelimParser :: extract_psm_features(
  int psmind, 
  vector<string> & tokens, 
  double *x
  ) {
  //cerr <<"extracting features for psm:"<<psmind<<endl;
  memset(x,0,sizeof(double)*num_features);

  int num_sequences = psmind_to_num_pep[psmind];
  int psm_offset = psmind_to_ofst[psmind];

  double xcorr = atof(tokens[xcorr_score_col_idx].c_str());
  double sp = atof(tokens[sp_score_col_idx].c_str());
  double log_sp_rank = 0; //TODO
  double by_ion_fraction_matched = 0;

  double avg_neutral_mass = 0.0;
  double max_neutral_mass = -1.0;

  double avg_weight_diff = 0.0;
  double max_weight_diff = -1.0;

  double avg_sequence_length = 0;
  double max_sequence_length = -1.0;

  double avg_missed_cleavages = 0;
  double max_missed_cleavages = 0;

  int charge_count[max_charge];
  for (int idx = 0; idx < max_charge;idx++) {
    charge_count[idx] = 0;
  }

  for (int idx = 0; idx < num_sequences;idx++) {

    double neutral_mass = psmind_to_neutral_mass[psm_offset+idx];
    double peptide_mass = psmind_to_peptide_mass[psm_offset+idx];
    int charge = min(psmind_to_charge[psm_offset+idx], max_charge);

    string& sequence = ind_to_pep[psmind_to_pepind[psm_offset+idx]];

    double weight_diff = fabs(neutral_mass - peptide_mass);
    
    avg_weight_diff += weight_diff;
    max_weight_diff = max(max_weight_diff, weight_diff);

    avg_neutral_mass += neutral_mass;
    max_neutral_mass = max(max_neutral_mass, neutral_mass);
    
    avg_sequence_length += sequence.length();
    max_sequence_length = max(max_sequence_length, (double)sequence.length());

    charge_count[min(max_charge,charge)-1]++;

  }

  avg_weight_diff = avg_weight_diff / (double)num_sequences;
  avg_neutral_mass = avg_neutral_mass / (double)num_sequences;
  avg_sequence_length = avg_sequence_length / (double)num_sequences;

  double ave_artd = 0; //TODO
  double max_artd = 0; //TODO

  double ln_experiment_size = logf(atof(tokens[matches_spectrum_col_idx].c_str()));
  
  //cerr << "Setting feature values"<<endl;

  x[0]  = xcorr;
  x[1]  = sp;
  x[2]  = log_sp_rank;
  x[3]  = by_ion_fraction_matched;
  x[4]  = avg_weight_diff;
  x[5]  = max_weight_diff;
  x[6]  = avg_neutral_mass;
  x[7]  = max_neutral_mass;
  x[8]  = ln_experiment_size;
  x[9]  = avg_missed_cleavages;
  x[10] = max_missed_cleavages;
  x[11] = avg_sequence_length;
  x[12] = max_sequence_length;
  x[13] = num_sequences;
  x[14] = ave_artd;
  x[15] = max_artd;
  x[13] = charge_count[0];
  x[14] = charge_count[1];
  x[15] = charge_count[2];
  x[16] = charge_count[3];
  x[17] = charge_count[4];
  x[18] = charge_count[5];

  /*
  cerr << psmind;
  for (int idx=0;idx < num_features;idx++) {
    cerr << " " << x[idx];
  }
  cerr<<endl;
  */
  //cerr <<"Done extracting features"<<endl;

}



void TabDelimParser :: second_pass(ifstream &fin, int label)
{
  string line;
  vector<string> tokens;
  
  vector<string> peptides;
  vector<string> charge_str;
  vector<string> spectrum_neutral_mass_str;
  vector<string> peptide_mass_str;
  string delim1 = "\t";
  string delim2 = ",";
  while(!fin.eof())
    {
      getline(fin,line);
      get_tokens(line, tokens, delim1);
  
      if(tokens.size() > 1)
	{
	  
	  //fill in tables

	  //get the scan
	  int scan = atoi(tokens[0].c_str());
	  psmind_to_scan[psmind] = scan;

	  //get charges
	  string ch = tokens[1];
	  get_tokens(ch, charge_str,delim2);
	  for(unsigned int i = 0; i < charge_str.size(); i++)
	    psmind_to_charge[curr_ofst+i] = atoi(charge_str[i].c_str());
	  
	  //get spectrum neutral mass
	  string neut_mass = tokens[3];
	  get_tokens(neut_mass, spectrum_neutral_mass_str,delim2);
	  for(unsigned int i = 0; i < spectrum_neutral_mass_str.size(); i++)
	    psmind_to_neutral_mass[curr_ofst+i] = atof(spectrum_neutral_mass_str[i].c_str());
	  
	  //get peptide mass
	  string pep_mass = tokens[4];
	  get_tokens(pep_mass, peptide_mass_str,delim2);
	  vector<double> peptide_mass;
	  peptide_mass.resize(peptide_mass_str.size(),0);
	  for(unsigned int i = 0; i < peptide_mass_str.size(); i++)
	    psmind_to_peptide_mass[curr_ofst+i] = atof(peptide_mass_str[i].c_str());
		  
	  //get all peptides
	  string peps = tokens[10];
	  get_tokens(peps,peptides,delim2);
	  for(unsigned int i = 0; i < peptides.size(); i++)
	    {
	      string pep = peptides[i];
	      int pep_ind = -1;
	      //add peptide to pep_to_ind and ind_to_pep maps
	      if(pep_to_ind.find(pep) == pep_to_ind.end())
		cout << "warning : did not find peptide in ind_to_pep_table\n"; 
	      else
		pep_ind = pep_to_ind[pep];
	      psmind_to_pepind[curr_ofst+i] = pep_ind;
	    }
	  psmind_to_num_pep[psmind] = peptides.size();
	  psmind_to_ofst[psmind] = curr_ofst;
	  
	  psmind_to_label[psmind] = label;

          //extract features
	  extract_psm_features(psmind, tokens, x);
	  f_psm.write((char*)x, sizeof(double)*num_features);

	  if(label == 1)
	    num_pos_psm++;
	  else
	    num_neg_psm++;

	  //augment counters
	  curr_ofst+= peptides.size();
	  psmind++;

  

  
	}
    }
}

/*****************************************************************************************/

void TabDelimParser :: save_data_in_binary(string out_dir)
{
  ostringstream fname;
  //write out data summary
  fname << out_dir << "/summary.txt";
  ofstream f_summary(fname.str().c_str());
  //psm info
  f_summary << num_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << " " << num_pep_in_all_psms << endl;
  //peptide info
  //f_summary << num_pep << endl;
  f_summary.close();
  fname.str("");
  
  //psmind_to_pepind
  fname << out_dir << "/psmind_to_pepind.txt";
  ofstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  f_psmind_to_pepind.write((char*)psmind_to_pepind,sizeof(int)*num_pep_in_all_psms);
  f_psmind_to_pepind.close();
  fname.str("");

  //psmind_to_neutral_mass
  fname << out_dir << "/psmind_to_neutral_mass.txt";
  ofstream f_psmind_to_neutral_mass(fname.str().c_str(),ios::binary);
  f_psmind_to_neutral_mass.write((char*)psmind_to_neutral_mass,sizeof(double)*num_pep_in_all_psms);
  f_psmind_to_neutral_mass.close();
  fname.str("");

  //psmind_to_peptide_mass
  fname << out_dir << "/psmind_to_peptide_mass.txt";
  ofstream f_psmind_to_peptide_mass(fname.str().c_str(),ios::binary);
  f_psmind_to_peptide_mass.write((char*)psmind_to_peptide_mass,sizeof(double)*num_pep_in_all_psms);
  f_psmind_to_peptide_mass.close();
  fname.str("");

  //psmind_to_num_pep
  fname << out_dir << "/psmind_to_num_pep.txt";
  ofstream f_psmind_to_num_pep(fname.str().c_str(),ios::binary);
  f_psmind_to_num_pep.write((char*)psmind_to_num_pep,sizeof(int)*num_psm);
  f_psmind_to_num_pep.close();
  fname.str("");

  //psmind_to_ofst
  fname << out_dir << "/psmind_to_ofst.txt";
  ofstream f_psmind_to_ofst(fname.str().c_str(),ios::binary);
  f_psmind_to_ofst.write((char*)psmind_to_ofst,sizeof(int)*num_psm);
  f_psmind_to_ofst.close();
  fname.str("");

  //psmind_to_scan
  fname << out_dir << "/psmind_to_scan.txt";
  ofstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  f_psmind_to_scan.write((char*)psmind_to_scan,sizeof(int)*num_psm);
  f_psmind_to_scan.close();
  fname.str("");

  //psmind_to_charge
  fname << out_dir << "/psmind_to_charge.txt";
  ofstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  f_psmind_to_charge.write((char*)psmind_to_charge,sizeof(int)*num_pep_in_all_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  ofstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  f_psmind_to_label.write((char*)psmind_to_label,sizeof(int)*num_psm);
  f_psmind_to_label.close();
  fname.str("");
  
  //ind_to_pep
  fname << out_dir << "/ind_to_pep.txt";
  ofstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = ind_to_pep.begin(); it != ind_to_pep.end(); it++)
    f_ind_to_pep << it->first << " " << it->second << "\n";
  f_ind_to_pep.close();
  fname.str("");

  //pep_to_ind
  fname << out_dir << "/pep_to_ind.txt";
  ofstream f_pep_to_ind(fname.str().c_str(),ios::binary);
  for(map<string,int>::iterator it = pep_to_ind.begin(); it != pep_to_ind.end(); it++)
    f_pep_to_ind << it->first << " " << it->second << "\n";
  f_pep_to_ind.close();
  fname.str("");

}

void TabDelimParser :: clean_up(string dir)
{
  
  cerr <<"Inside clean_up"<<endl;
  ostringstream fname;
      
  //fname << out_dir << "/summary.txt";
  //ofstream f_summary(fname.str().c_str());
  
  fname << dir << "/psm.txt";
  remove(fname.str().c_str());
  fname.str("");
  
  //psmind_to_pepind
  fname << out_dir << "/psmind_to_pepind.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_num_pep
  fname << out_dir << "/psmind_to_num_pep.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_num_pep
  fname << out_dir << "/psmind_to_ofst.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_scan
  fname << out_dir << "/psmind_to_scan.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_charge
  fname << out_dir << "/psmind_to_charge.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  remove(fname.str().c_str());
  fname.str("");

  //ind_to_pep
  fname << out_dir << "/ind_to_pep.txt";
  remove(fname.str().c_str());
  fname.str("");

  //pep_to_ind
  fname << out_dir << "/pep_to_ind.txt";
  remove(fname.str().c_str());
  fname.str("");

}



/******************************************************************************************/


int TabDelimParser :: run(vector<string> &filenames)
{
  string line;
  for(unsigned int i = 0; i < filenames.size(); i++)
    {
      string fname = filenames[i];
      cout << fname << endl;
      ifstream fin(fname.c_str());
      if(!fin.is_open())
	{
	  cout << "could not open " << fname << " for reading" << endl;
	  return 0;
	}
      getline(fin,line);
      first_pass(fin);
      fin.close();
    }
  cout << num_psm << " " << num_pep << " " << num_pep_in_all_psms << endl;
  allocate_feature_space();

  ostringstream fname;
  fname << out_dir << "/psm.txt";
  f_psm.open(fname.str().c_str());
  fname.str("");

  int label = 0;
  for(unsigned int i = 0; i < filenames.size(); i++)
    {
      string fname = filenames[i];
      cout << fname << endl;
      ifstream fin(fname.c_str());
      if(!fin.is_open())
	{
	  cout << "could not open " << fname << " for reading" << endl;
	  return 0;
	}
      getline(fin,line);
      if(i == 0)
	label = 1;
      else
	label = -1;
      second_pass(fin,label);
      fin.close();
    }
  f_psm.close();
  cout << psmind << " " << num_pos_psm << " " << num_neg_psm << endl;
  cerr <<"Saving data "<<out_dir<<endl;
  save_data_in_binary(out_dir);
  cerr <<"Returning"<<endl;
  return 1;
}


