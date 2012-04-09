#include "TabDelimParser.h"

/******************************/

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
	x(0),
	xs(0),
	num_xlink_features(0)
{
  
  //num_psm_features
  num_features = 17;
  //num_xlink_features
  num_xlink_features = 10;
  //num_spec_features 
  num_spec_features = 0;
  
  //final_hits_per_spectrum
  fhps = 1;
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
  delete[] xs; xs = (double*)0;
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
  tokens.erase(tokens.begin(),tokens.end());
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


/* xlink tokens
 * 0. scan    
 * 1. charge  
 * 2. spectrum precursor m/z  
 * 3. spectrum neutral mass   
 * 4. peptide mass
 * 5. delta_cn
 * 6. sp score
 * 7. sp rank
 * 8. xcorr score
 * 9. xcorr rank
 * 10. p-value
 * 11. b/y ions matched
 * 12. b/y ions total
 * 13. matches/spectrum
 * 14. sequence
 * 15. cleavage type 
 * 16. protein id
 * 17. flanking aa
 */

const static int scan_idx=0;
const static int charge_idx=1;
const static int spectrum_mz_idx=2;
const static int spectrum_mass_idx=3;
const static int peptide_mass_idx=4;
const static int delta_cn_idx=5;
const static int sp_score_idx=6;
const static int sp_rank_idx=7;
const static int xcorr_idx=8;
const static int xcorr_rank=9;
const static int pvalue_idx=10;
const static int by_matched_idx=11;
const static int by_total_idx=12;
const static int matches_idx=13;
const static int sequence_idx=14;



void TabDelimParser :: first_pass_xlink(ifstream &fin)
{
  string line;
  vector<string> tokens;
  vector<string> subtokens;
  vector<string> peptides;
  string delim1 = "\t";
  string delim2 = ",";
  string delim3 = " ";
  while(!fin.eof())
    {
      getline(fin,line);
      get_tokens(line, tokens, delim1);
      if(tokens.size() > 1)
	{
	  //get all peptides and fill tables
	  string pep_and_loc = tokens[sequence_idx];
	  get_tokens(pep_and_loc,subtokens,delim3);
	  //cout << pep_and_loc << endl;;
	  //for(size_t i = 0; i < subtokens.size(); i++)
          //  cout<< i << " " <<subtokens[i] << endl;
	  
          if(subtokens.size() > 0) {
            if (subtokens.size() == 2) {
              //it it either a linear peptide or a self loop.
              psmind_to_peptide1[num_psm] = subtokens[0];
              psmind_to_peptide2[num_psm] = "_";
              psmind_to_loc[num_psm] = subtokens[1];
              //cout << num_psm << endl;
              //cout << psmind_to_peptide1[num_psm] << endl;
              //cout << psmind_to_peptide2[num_psm] << endl;
              //cout << psmind_to_loc[num_psm] << endl;
            } else if (subtokens.size() == 3) {
              //it is a crosslinked peptide.
              //cerr << "parsing cross linked peptide:"<<subtokens[0]<<endl;
              string pep1 = subtokens.at(0).substr(0,subtokens[0].length()-1);
              //cerr << "pep1:"<<pep1<<endl;
              psmind_to_peptide1[num_psm] = pep1;
              psmind_to_peptide2[num_psm] = subtokens.at(1);
              psmind_to_loc[num_psm] = subtokens.at(2);
            } else {
              cerr << "Error!"<<endl;
              exit(-1);
            }
          
            
          }
	  //proteins
	  if(tokens[16].size() > 0)
	    psmind_to_protein1[num_psm] = tokens[16];
	  //if(tokens[17].size()>0)
	  //psmind_to_protein2[num_psm] = tokens[17];
	  num_psm++;
	}
    }
}


void TabDelimParser :: allocate_feature_space_xlink()
{
  //space for feature vector
  x = new double[num_xlink_features];
  memset(x,0,sizeof(double)*num_xlink_features);
  //space for spec feature vector
  if (num_spec_features > 0)
    {
      xs = new double[num_spec_features];
      memset(xs,0,sizeof(double)*num_spec_features);
    }
  psmind_to_label = new int[num_psm];
  memset(psmind_to_label,0,sizeof(int)*num_psm);
  psmind_to_scan = new int[num_psm];
  memset(psmind_to_scan,0,sizeof(int)*num_psm);
  psmind_to_charge = new int[num_psm];
  memset(psmind_to_charge,0,sizeof(int)*num_psm);


}

int TabDelimParser::get_peptide_length_sum(string& sequence) {

  string delim3 = " ";
  vector<string> subtokens;
  get_tokens(sequence,subtokens,delim3);
/*
  cout << sequence << endl;;
  for(size_t i = 0; i < subtokens.size(); i++) {
    cout<< i << " " <<subtokens[i] << endl;
  }
*/
  int ans = -1;

  if(subtokens.size() > 0) {
    if (subtokens.size() == 2) {
      //it it either a linear peptide or a self loop.
      cout << subtokens[0] <<" " << subtokens[0].length()<<endl;
      ans = subtokens[0].length();
      return subtokens[0].length();
    } else if (subtokens.size() == 3) {
      //it is a crosslinked peptide.
      cerr << "parsing cross linked peptide:" << subtokens[0] << endl;
      string pep1 = subtokens.at(0).substr(0,subtokens[0].length()-1);
      string pep2 = subtokens.at(1);
      cerr << pep1 <<" "<<pep1.length()<<endl;
      cerr << pep2 <<" "<<pep2.length()<<endl;
      ans = pep1.length() + pep2.length();
    } else {
      cerr << "getPeptideLengthSum:error:" << sequence << endl;
      exit(-1);
    }
  }
  return ans;

}



int TabDelimParser::get_peptide_type(string& sequence) {
  string delim3 = " ";
  vector<string> subtokens;
  get_tokens(sequence,subtokens,delim3);
/*
  cout << sequence << endl;;
  for(size_t i = 0; i < subtokens.size(); i++) {
    cout<< i << " " <<subtokens[i] << endl;
  }
*/
  int ans = 0;

  if(subtokens.size() > 0) {
    if (subtokens.size() == 2) {
      //it it either a linear peptide or a self loop.
      cout << subtokens[0] <<" " << subtokens[0].length()<<endl;
      ans = subtokens[0].length();
      if (subtokens[1].find(',') == string::npos) {
        //it is linear.
        cerr << sequence << " is linear"<<endl;
        ans = 0;
      } else {
        //it is a self-loop
        ans = 1;
        cerr << sequence << " is self loop"<<endl;
      }
    } else if (subtokens.size() == 3) {
      //it is a crosslinked peptide.
      ans = 2;
      cerr << sequence << " is cross-linked"<<endl;
    } else {
      cerr <<"get_peptide_type:error:"<<sequence<<endl;
      exit(-1);
    }
  }

  return ans;
}


/*
 0. XCorr Score
 1. Peptide Length
 2. product type (0 - linear peptide, 1 - self-loop peptide, 2- cross-linked peptide)
 3. mass diff
 4. abs mass diff
 5. sp score (unused)
 6. frac matched
 7. obs mass
 8. charge
 9. matches/spectrum (unused) 
*/
void TabDelimParser :: extract_xlink_features(vector<string> & tokens, double *x)
{
  memset(x,0,sizeof(double)*num_xlink_features);
  /*
  cerr << "scan:" << tokens[scan_idx] << 
          " charge:" << tokens[charge_idx] << 
          " sequence:" << tokens[sequence_idx] << 
          " xcorr:" << tokens[xcorr_idx] << 
          " sp:" << tokens[sp_score_idx] <<
          " sp rank:" << tokens[sp_rank_idx] << 
          " spectrum mass:" << tokens[spectrum_mass_idx] <<
          " peptide mass:" << tokens[peptide_mass_idx] << endl;
  */

  //xcorr score
  x[0] = atof(tokens[xcorr_idx].c_str());
  //x[0] = atof(tokens[pvalue_idx].c_str());

  /* peptide length */
  x[1] = get_peptide_length_sum(tokens[sequence_idx]);

  // peptide type (linear, self-loop, cross-link)
  x[2] = get_peptide_type(tokens[sequence_idx]);

/*     
  //log rank by Sp
  x[1] = 0;
  if(atof(tokens[sp_rank_idx].c_str()) > 0)
    x[1]=log(atof(tokens[sp_rank_idx].c_str()));
  //deltaCN
  //  x[2] = -log(atof(tokens[pvalue_idx].c_str()));
  */
  
  //difference between measured and calculated mass
  x[3] = atof(tokens[spectrum_mass_idx].c_str())-atof(tokens[peptide_mass_idx].c_str());
  
  // absolute value of difference between measured and calculated mass
  x[4] = fabs(atof(tokens[spectrum_mass_idx].c_str())-atof(tokens[peptide_mass_idx].c_str()));

  //sp score
  //x[5] = atof(tokens[sp_score_idx].c_str());

  //matched ions/predicted ions
  x[6] = 0;
  if(atof(tokens[by_total_idx].c_str()) != 0)
    x[6] = atof(tokens[by_matched_idx].c_str())/atof(tokens[by_total_idx].c_str());

  //observed mass
  x[7] = atof(tokens[spectrum_mass_idx].c_str());

  //charge
  x[8] = atof(tokens[charge_idx].c_str());
  // number of sequence_comparisons
  //x[9] = log(atof(tokens[matches_idx].c_str()));
  //whether n-terminus and c-terminus have proper cleavage sites
  // missed cleavages

}


void TabDelimParser :: second_pass_xlink(ifstream &fin, int label)
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
	  //get the scan and the charge
	  int scan = atoi(tokens[0].c_str());
	  int charge = atoi(tokens[1].c_str());
	  //extract features
	  extract_xlink_features(tokens, x);
	  f_psm.write((char*)x, sizeof(double)*num_xlink_features);
	  if(num_spec_features > 0)
	    {
	      
	      //clear out the spec features
	      memset(xs,0,sizeof(double)*num_spec_features);
	      
	      string pept;
	      
	      if(psmind_to_peptide1.find(psmind) != psmind_to_peptide1.end())
		{
		  pept = psmind_to_peptide1[psmind];
		  //cout << pept << endl;
		  if(num_spec_features == 7)
		    sfg.get_spec_features_m3(scan, charge,pept,xs,1);
		}
	      
	      if(psmind_to_peptide2.find(psmind) != psmind_to_peptide2.end())
		{
		  pept = psmind_to_peptide2[psmind];
		  //cout << psmind << " " << pept << endl;
		  if(num_spec_features == 7)
		    sfg.get_spec_features_m3(scan, charge,pept,xs,0);
		  }
	      f_psm.write((char*)xs, sizeof(double)*num_spec_features);
	      
	    }
  
	  //record charge and scan
	  psmind_to_scan[psmind] = scan;
	  psmind_to_charge[psmind] = charge;

	  psmind_to_label[psmind] = label;
	  if(label == 1)
	    num_pos_psm++;
	  else
	    num_neg_psm++;

	  //augment counters
	  psmind++;
	  if(psmind % 100 == 0)
	    cout << "psmind " << psmind << endl;
	}
    }
}


void TabDelimParser :: save_data_in_binary_xlink(string out_dir)
{
  ostringstream fname;
  //write out data summary
  fname << out_dir << "/summary.txt";
  ofstream f_summary(fname.str().c_str());
  //psm info
  f_summary << num_xlink_features+num_spec_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  cout << num_xlink_features+num_spec_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  f_summary.close();
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  ofstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  f_psmind_to_label.write((char*)psmind_to_label,sizeof(int)*num_psm);
  f_psmind_to_label.close();
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
  f_psmind_to_charge.write((char*)psmind_to_charge,sizeof(int)*num_psm);
  f_psmind_to_charge.close();
  fname.str("");
  
  //psmind_to_peptide1
  fname << out_dir << "/psmind_to_peptide1.txt";
  ofstream f_psmind_to_peptide1(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_peptide1.begin(); it != psmind_to_peptide1.end(); it++)
    f_psmind_to_peptide1 << it->first << " " << it->second << "\n";
  f_psmind_to_peptide1.close();
  fname.str("");
  
  //psmind_to_peptide2
  fname << out_dir << "/psmind_to_peptide2.txt";
  ofstream f_psmind_to_peptide2(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_peptide2.begin(); it != psmind_to_peptide2.end(); it++)
    f_psmind_to_peptide2 << it->first << " " << it->second << "\n";
  f_psmind_to_peptide2.close();
  fname.str("");

  //psmind_to_loc
  fname << out_dir << "/psmind_to_loc.txt";
  ofstream f_psmind_to_loc(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_loc.begin(); it != psmind_to_loc.end(); it++)
    f_psmind_to_loc << it->first << " " << it->second << "\n";
  f_psmind_to_loc.close();
  fname.str("");

  //psmind_to_protein1
  fname << out_dir << "/psmind_to_protein1.txt";
  ofstream f_psmind_to_protein1(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_protein1.begin(); it != psmind_to_protein1.end(); it++)
    f_psmind_to_protein1 << it->first << " " << it->second << "\n";
  f_psmind_to_protein1.close();
  fname.str("");
  
  //psmind_to_protein2
  fname << out_dir << "/psmind_to_protein2.txt";
  ofstream f_psmind_to_protein2(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_protein2.begin(); it != psmind_to_protein2.end(); it++)
    f_psmind_to_protein2 << it->first << " " << it->second << "\n";
  f_psmind_to_protein2.close();
  fname.str("");

}


void TabDelimParser :: clean_up_xlink(string dir)
{

  ostringstream fname;
      
  //fname << out_dir << "/summary.txt";
  //ofstream f_summary(fname.str().c_str());

  fname << dir << "/psm.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
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

  //psmind_to_peptide1
  fname << out_dir << "/psmind_to_peptide1.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_peptide2
  fname << out_dir << "/psmind_to_peptide2.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_loc
  fname << out_dir << "/psmind_to_loc.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_protein1
  fname << out_dir << "/psmind_to_protein1.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_protein2
  fname << out_dir << "/psmind_to_protein2.txt";
  remove(fname.str().c_str());
  fname.str("");

}



int TabDelimParser :: run_on_xlink(vector<string> &filenames)
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
      first_pass_xlink(fin);
      fin.close();
    
    }

  cout << num_psm  << endl;
  allocate_feature_space_xlink();
    
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
      second_pass_xlink(fin,label);
      fin.close();
    }
  f_psm.close();
  cout << psmind << " " << num_pos_psm << " " << num_neg_psm << endl;
  save_data_in_binary_xlink(out_dir);
  return 1;
}

int TabDelimParser :: run_on_xlink(vector<string> &filenames, string &ms2filename)
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
      first_pass_xlink(fin);
      fin.close();
    
    }

  cout << num_psm  << endl;
  num_spec_features = 7;
  allocate_feature_space_xlink();
    
  //prepare to generate spectrum features
  sfg.clear();
  if(!sfg.open_ms2_file_for_reading(ms2filename))
    return 0;
  sfg.read_ms2_file();
  sfg.initialize_aa_tables();
  
  
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
      second_pass_xlink(fin,label);
      fin.close();
    }
  f_psm.close();
  cout << psmind << " " << num_pos_psm << " " << num_neg_psm << endl;
  save_data_in_binary_xlink(out_dir);
  return 1;
}


