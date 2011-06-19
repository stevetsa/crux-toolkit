#include "SQTParser.h"

/******************************/

SQTParser :: SQTParser() 
  : num_mixed_labels(0),
    num_features(0),
    num_spec_features(0),
    num_spectra(0),
    num_psm(0),
    num_pos_psm(0),
    num_neg_psm(0),
    num_pep(0),
    num_pos_pep(0),
    num_neg_pep(0),
    num_prot(0),
    num_pos_prot(0),
    num_neg_prot(0),
    psmind(0),
    x(0), 
    xs(0), 
    psmind_to_pepind(0), 
    psmind_to_scan(0),
    psmind_to_charge(0), 
    psmind_to_label(0), 
    pepind_to_label(0),
    protind_to_label(0), 
    protind_to_num_all_pep(0)
{
  int capacity = 10;
  m.xcorr_rank.reserve(capacity);
  m.sp_rank.reserve(capacity);
  m.calc_mass.reserve(capacity);
  m.delta_cn.reserve(capacity);
  m.xcorr_score.reserve(capacity);
  m.sp_score.reserve(capacity);
  m.num_ions_matched.reserve(capacity);
  m.num_total_ions.reserve(capacity);
  m.peptides.reserve(capacity);
  m.num_proteins_in_match.reserve(capacity);
  m.proteins.reserve(capacity);

  //num_psm_features
  num_features = 17;
  //num_spec_features 
  num_spec_features = 7;
  
  //final_hits_per_spectrum
  fhps = 3;
  //enzyme
  e = TRYPSIN_ENZ;
  //decoy prefix
  decoy_prefix = "decoy_";
  //max peptide length to be considered
  max_len = 50;
  //min peptide length to be considered
   min_len = 7;
  
}


SQTParser :: ~SQTParser()
{
  clear();  
}

void SQTParser :: clear()
{
  clear_matches();
  delete[] x; x = (double*)0;
  delete[] xs; xs = (double*)0;
  delete[] psmind_to_pepind; psmind_to_pepind = (int*)0;
  delete[] psmind_to_scan; psmind_to_scan = (int*)0;
  delete[] psmind_to_charge; psmind_to_charge = (int*)0;
  delete[] psmind_to_label; psmind_to_label = (int*)0;
  delete[] pepind_to_label; pepind_to_label = (int*)0;
  delete[] protind_to_label; protind_to_label = (int*)0;
  delete[] protind_to_num_all_pep; protind_to_num_all_pep = (int*)0;
}

void SQTParser :: set_enzyme(string &enz)
{
  if(enz.find("elastase") != string::npos)
    {
      e = ELASTASE_ENZ;
    }
  else if (enz.find("chymotrypsin") != string::npos)
    {
      e = CHYMOTRYPSIN_ENZ;
    }
  else if (enz.find("trypsin") != string::npos)
    {
      e = TRYPSIN_ENZ;
    }
  else
    carp(CARP_WARNING, "could not determine enzyme, will assume trypsin");
}


void SQTParser :: clear_matches()
{
  m.xcorr_rank.clear();
  m.sp_rank.clear();
  m.calc_mass.clear();
  m.delta_cn.clear();
  m.xcorr_score.clear();
  m.sp_score.clear();
  m.num_ions_matched.clear();
  m.num_total_ions.clear();
  m.peptides.clear();
  m.num_proteins_in_match.clear();
  m.proteins.clear();
}


void SQTParser :: erase_matches()
{
  m.xcorr_rank.erase(m.xcorr_rank.begin(),m.xcorr_rank.end());
  m.sp_rank.erase(m.sp_rank.begin(),m.sp_rank.end());
  m.calc_mass.erase(m.calc_mass.begin(),m.calc_mass.end());
  m.delta_cn.erase(m.delta_cn.begin(),m.delta_cn.end());
  m.xcorr_score.erase(m.xcorr_score.begin(),m.xcorr_score.end());
  m.sp_score.erase(m.sp_score.begin(),m.sp_score.end());
  m.num_ions_matched.erase(m.num_ions_matched.begin(),m.num_ions_matched.end());
  m.num_total_ions.erase(m.num_total_ions.begin(),m.num_total_ions.end());
  m.peptides.erase(m.peptides.begin(),m.peptides.end());
  m.num_proteins_in_match.erase(m.num_proteins_in_match.begin(),m.num_proteins_in_match.end());
  m.proteins.erase(m.proteins.begin(),m.proteins.end());
}


void SQTParser :: add_matches_to_tables(sqt_match &m, string &decoy_prefix, int hits_read, int final_hits)
{
  int protein_pos = 0;
  for (int i = 0; i < min(hits_read,final_hits); i++)
    {
      set<string> proteins;
      int label = 0;
      // go through the proteins of the match
      for (int j = 0; j < m.num_proteins_in_match[i]; j++)
	{
	  string prot = m.proteins[protein_pos];
	  proteins.insert(prot);
	  if(prot.find(decoy_prefix) != string::npos)
	    label = -1;
	  else
	    label = 1;
	  protein_pos++;
	}
      //record the psm label
      psmind_to_label_map[num_psm] = label;

      string pep = m.peptides[i];
      int pep_ind = -1;
      //add peptide to pep_to_ind and ind_to_pep maps
      if(pep_to_ind.find(pep) == pep_to_ind.end())
	{
	  pep_ind = num_pep;
	  //add a new peptide
	  pep_to_ind[pep] = pep_ind;
	  ind_to_pep[pep_ind] = pep;
	  //start a pepind_to_psminds entry for the new peptide
	  set<int> tpsm;
	  pepind_to_psminds_map[pep_ind] = tpsm;
	  //start a pepind_to_protinds entry for the new peptide
	  set<int> t;
	  pepind_to_protinds_map[pep_ind] = t;
	  //set the label of the new peptide
	  pepind_to_label_map[pep_ind] = label;
 	  //augment num_pep count
	  num_pep++;
	  if(label == 1)
	    num_pos_pep++;
	  else
	    num_neg_pep++;
	}
      else
	{
	  pep_ind = pep_to_ind[pep];
	  string p = ind_to_pep[pep_ind];
	  //if(pep.compare(p) != 0)
	  //   cout << "warning : did not find peptide index in ind_to_pep_table\n"; 
	}
      
      for(set<string>::iterator it = proteins.begin(); it != proteins.end();it++)
	{
	  string prot = *it;
	  //if the label is -1 but the protein name does not contain decoy_prefix_,
	  // we don't include it
	  if((prot.find(decoy_prefix) == string::npos) && (label == -1))
	    num_mixed_labels++;
	  else
	    {
	      //add prot to tables
	      int prot_ind = -1;
	      if (prot_to_ind.find(prot) == prot_to_ind.end())
		{
		  prot_ind = num_prot;
		  //add new protein
		  prot_to_ind[prot] = prot_ind;
		  ind_to_prot[prot_ind] = prot;
		  //start a protind_to_pepinds entry for the new protein
		  set<int> t;
		  protind_to_pepinds_map[prot_ind] = t;
		  //set the prot label
		  protind_to_label_map[prot_ind] = label;
  		  num_prot++;
		  if(label == 1)
		    num_pos_prot++;
		  else
		    num_neg_prot++;

		  //find the prot in the prot_to_num_all_pep
		  int cnt = 0;
		  if(protein_to_num_all_pep_map.find(prot) == protein_to_num_all_pep_map.end())
		    {
		      for(map<string,int>::iterator itt = protein_to_num_all_pep_map.begin();
			  itt != protein_to_num_all_pep_map.end(); itt++)
			{
			  string protein = itt->first;
			  if(protein.find(prot) != string :: npos)
			    cnt = itt->second;
			}
		    }
		  else
		    cnt = protein_to_num_all_pep_map[prot];
		  //add the cnt to protind_to_num_all_pep_map
		  if(cnt == 0)
		    {
		      carp(CARP_WARNING, "did not find protein %s from sqt file in the database ", prot.c_str());
		    }
		  else
		    protind_to_num_all_pep_map[prot_ind] = cnt;
		}
	      else
		{
		  prot_ind = prot_to_ind[prot];
		  string p = ind_to_prot[prot_ind];
		  if(prot.compare(p) != 0)
		    {
		      cout << "did not find protein in the ind_to_prot_table\n";
		    }
		}
	      //augment the pepinds_to_psminds table
	      (pepind_to_psminds_map[pep_ind]).insert(num_psm);
	      //augment the pepind_to_protinds table
	      (pepind_to_protinds_map[pep_ind]).insert(prot_ind);
	      //augment to protind_to_pepinds table
	      protind_to_pepinds_map[prot_ind].insert(pep_ind);
	    }
	  
	}
      //remember the file a psm belongs to
      psmind_to_fname[num_psm] = cur_fname;

       //augment num psms
      num_psm++;
      if (label == 1)
	num_pos_psm++;
      else
	num_neg_psm++;
    }
  num_spectra++;
}


void SQTParser :: allocate_feature_space()
{
  //space for feature vector
  x = new double[num_features];
  memset(x,0,sizeof(double)*num_features);
  //space for spec feature vector
  if (num_spec_features > 0)
    {
      xs = new double[num_spec_features];
      memset(xs,0,sizeof(double)*num_spec_features);
    }
  //space for psminfo
  psmind_to_pepind = new int[num_psm];
  memset(psmind_to_pepind,0,sizeof(int)*num_psm);
  psmind_to_scan = new int[num_psm];
  memset(psmind_to_scan,0,sizeof(int)*num_psm);
  psmind_to_charge = new int[num_psm];
  memset(psmind_to_charge,0,sizeof(int)*num_psm);
  psmind_to_label = new int[num_psm];
  memset(psmind_to_label,0,sizeof(int)*num_psm);
  //space for pep info
  pepind_to_label = new int[num_pep];
  memset(pepind_to_label,0,sizeof(int)*num_pep);
  //space for prot info
  protind_to_label = new int[num_prot];
  memset(protind_to_label,0,sizeof(int)*num_prot);
  //space for prot to number of all peptides info
  protind_to_num_all_pep = new int[num_prot];
  memset(protind_to_num_all_pep,0,sizeof(int)*num_prot);
}


void SQTParser :: fill_graphs()
{
  pepind_to_psminds.create_bipartite_graph(pepind_to_psminds_map);
  pepind_to_protinds.create_bipartite_graph(pepind_to_protinds_map);
  protind_to_pepinds.create_bipartite_graph(protind_to_pepinds_map);
  for(map<int,int>::iterator it = psmind_to_label_map.begin(); it != psmind_to_label_map.end(); it++)
    psmind_to_label[it->first] = it->second;
  for(map<int,int>::iterator it = pepind_to_label_map.begin(); it != pepind_to_label_map.end(); it++)
    pepind_to_label[it->first] = it->second;
  for(map<int,int>::iterator it = protind_to_label_map.begin(); it != protind_to_label_map.end(); it++)
    {
      protind_to_label[it->first] = it->second;
      int protind = it->first;
      //if did not find the protein in the count of all proteins, then just get the cound of observed proteins
      if(protind_to_num_all_pep_map.find(protind) == protind_to_num_all_pep_map.end())
	{
	  int cnt = (protind_to_pepinds_map[protind]).size();
	  protind_to_num_all_pep[protind] = cnt;
	}
      else
	{
	  int cnt = protind_to_num_all_pep_map[protind];
	  protind_to_num_all_pep[protind] = cnt;
	}
    }
  pepind_to_protinds_map.clear();
  protind_to_pepinds_map.clear();
  psmind_to_label_map.clear();
  pepind_to_label_map.clear();
  protind_to_label_map.clear();
  protein_to_num_all_pep_map.clear();
  protind_to_num_all_pep_map.clear();
} 


  
/********* extracting features **********************************************************/
int SQTParser::cntEnz(const string& peptide,enzyme enz) {
    unsigned int pos=2, cnt=0;
    char n = peptide.at(pos++);
    while (pos<peptide.size()-2) {
      char c = peptide.at(pos++);
      if (isEnz(n,c,enz))
        cnt++;
      n=c;
    }
    return cnt;
}

double SQTParser::isTryptic(const char n,const char c) {
  return (
  (((n=='K' || n=='R') && c != 'P') ||
  n=='-' || c=='-')
  ?1.0:0.0);
}
// [FHWYLM].[^P]
double SQTParser::isChymoTryptic(const char n,const char c) {
  return (
  (((n=='F' || n=='H' || n=='W' || n=='Y' || n=='L' || n=='M') && c!= 'P') ||
  n=='-' || c=='-')
  ?1.0:0.0);
}

// [LVAG].[^P]
double SQTParser::isElastasic(const char n,const char c) {
  return (
  (((n=='L' || n=='V' || n=='A' || n=='G' ) && c!= 'P') ||
  n=='-' || c=='-')
  ?1.0:0.0);
}

double SQTParser::isEnz(const char n,const char c, enzyme enz) {
    switch(enz) {
      case TRYPSIN_ENZ:
        return isTryptic(n,c);
      case CHYMOTRYPSIN_ENZ:
        return isChymoTryptic(n,c);
      case ELASTASE_ENZ:
        return isElastasic(n,c);
    default:
        return 0;
    }
}

void SQTParser :: extract_psm_features(sqt_match &m, enzyme enz, double *x, int i)
{
  string pep = m.peptides[i];
  memset(x,0,sizeof(double)*num_features);
  //log rank by Sp
  if (m.sp_rank[i] > 0)
    x[0] = log(m.sp_rank[i]);
  else
    x[0] = 0.0;
  //deltaLCN
  x[1] = 0.0;
  //deltaCN
  x[2] = m.delta_cn[i];
  //xcorr score
  x[3] = m.xcorr_score[i];
  //sp score
  x[4] = m.sp_score[i];
  //matched ions/predicted ions
  if(m.num_total_ions[i] > 0)
    x[5] = m.num_ions_matched[i]/m.num_total_ions[i];
  else
    x[5] = 0.0;
  //observed mass
  x[6] = m.precursor_mass;
  //peptide length
  x[7] = pep.size();
  //charge
  x[8] = 0.0; x[9] = 0.0; x[10] = 0.0;
  if (m.charge == 1)
    x[8] = 1.0;
  if(m.charge == 2)
    x[9] = 1.0;
  if(m.charge == 3)
    x[10] = 1.0;
  //whether n-terminus and c-terminus have proper cleavage sites
  x[11]=isEnz(pep.at(0),pep.at(2),enz);        
  x[12]=isEnz(pep.at(pep.size()-3),pep.at(pep.size()-1),enz);
  // missed cleavages
  x[13]=(double)cntEnz(pep,enz);
  // number of sequence_comparisons
  //x[14] = log(m.num_sequence_comparisons);
  x[14] = m.num_sequence_comparisons;
  //difference between measured and calculated mass
  x[15] = m.precursor_mass-m.calc_mass[i];
  // absolute value of difference between measured and calculated mass
  x[16] = m.precursor_mass-m.calc_mass[i];
}


void SQTParser :: extract_features(sqt_match &m, string &decoy_prefix, int hits_read, int final_hits,enzyme enz)
{
  int protein_pos = 0;
  for (int i = 0; i < min(hits_read,final_hits); i++)
    {
      //write the feature vector out to file
      extract_psm_features(m, enz, x, i);
      
      if (num_spec_features > 0)
	{
	  ostringstream scan_stream;
	  scan_stream << m.scan << "." << m.charge;
	  string scan_str = scan_stream.str();
	  //if((m.scan %10000) == 0 && i == 0)
	  //cout << scan_str << endl;
	  if(psmind % 5000 == 0)
	    cout << "PMS number " << psmind << endl;
	  string peptide = m.peptides[i];
	  int pos = peptide.find(".");
	  string pept = peptide.substr(pos+1,peptide.size());
	  pos = pept.find(".");
	  pept = pept.substr(0,pos);
	  
	  if(num_spec_features == 3)
	      sfg.get_spec_features_m3(m.scan, m.charge,pept,xs);
	  if(num_spec_features == 6)
	    sfg.get_spec_features_m6(m.scan, m.charge,pept,xs);
	  if(num_spec_features == 7)
	    sfg.get_spec_features_m7(m.scan, m.charge,pept,xs);
	  
	  if(num_spec_features == 59)
	    {
	      sfg.get_spec_features_m7(m.scan, m.charge,pept,xs);
	      sfg.get_spec_features_aa_end(m.scan, m.charge,pept,xs+7);
	    }

	  //write out features
	  f_psm.write((char*)x, sizeof(double)*num_features);
	  f_psm.write((char*)xs, sizeof(double)*num_spec_features);

	}
      else
	f_psm.write((char*)x, sizeof(double)*num_features);

      //write psm tables
      psmind_to_scan[psmind] = m.scan;
      psmind_to_charge[psmind] = m.charge;
      //get the pepind of the peptide
      string pep = m.peptides[i];
      int pepind = pep_to_ind[pep];
      psmind_to_pepind[psmind] = pepind;
          
      //check the protein to peptide and peptide to protein tables
      //collect proteins
      set<string> proteins;
      int label = 0;
      // go through the proteins of the match
      for (int j = 0; j < m.num_proteins_in_match[i]; j++)
	{
	  string prot = m.proteins[protein_pos];
	  proteins.insert(prot);
	  if(prot.find(decoy_prefix) != string::npos)
	    label = -1;
	  else
	    label = 1;
	  protein_pos++;
	}
      assert(psmind_to_label[psmind] == label);
      //assert(pepind_to_label[pepind] == label);

      //go through the proteins
      for(set<string>::iterator it = proteins.begin(); it != proteins.end();it++)
	{
	  string prot = *it;
	  //if the label is -1 but the protein name does not contain decoy_prefix_,
	  // we don't include it
	  if((prot.find(decoy_prefix) == string::npos) && (label == -1))
	    ;
	  else
	    {
	      int protind = prot_to_ind[prot];
	      if(!protind_to_pepinds.is_index_in_range(pepind,protind))
#ifndef CRUX
		cout << "warning: peptide was not found in protein to peptide table\n";
#endif
	      if(!pepind_to_protinds.is_index_in_range(protind,pepind))
#ifndef CRUX
		cout << "warning: peptide was not found in peptide to protein table\n";
#endif
	      assert(protind_to_label[protind] == label);
	    }
	}
      psmind++;
    }
}



/************ parsing sqt file*******************************/

void SQTParser :: read_M_line(ifstream &is, sqt_match &m)
{
  //rank by scorr
  int xcorr_rank;
  is >> xcorr_rank;
  m.xcorr_rank.push_back(xcorr_rank);
  
  //rank by Sp
  int sp_rank;
  is >> sp_rank;
  m.sp_rank.push_back(sp_rank);
  
  //calculated mass
  double calc_mass;
  is >> calc_mass;
  m.calc_mass.push_back(calc_mass);

  //delta cn
  double delta_cn;
  is >> delta_cn;
  m.delta_cn.push_back(delta_cn);

  //xcorr
  double xcorr;
  is >> xcorr;
  m.xcorr_score.push_back(xcorr);

  //sp_score
  double sp;
  is >> sp;
  m.sp_score.push_back(sp);
  
  //number of matched ions
  double num_ions_matched;
  is >> num_ions_matched;
  m.num_ions_matched.push_back(num_ions_matched);

  //number of total ions;
  double num_total_ions;
  is >> num_total_ions;
  m.num_total_ions.push_back(num_total_ions);

  //peptide
  string peptide;
  is >> peptide;
  m.peptides.push_back(peptide);
}


void SQTParser :: read_S_line(ifstream &is, sqt_match &m)
{
  string tempstr;
  //scan begin
  is >> m.scan;
  //scan end
  is >> tempstr;
  //charge
  is >> m.charge;
  is >> tempstr;
  is >> tempstr;
  //precursor ion mass
  is >> m.precursor_mass;
  is >> tempstr;
  is >> tempstr;
  //number of matches considered
  is >> m.num_sequence_comparisons;
}

int SQTParser :: parse_sqt_spectrum_matches(ifstream &is, sqt_match &m)
{
  string tempstr;
  read_S_line(is,m);
  erase_matches();
  int num_hits = 0;
  int num_proteins_in_match = 0;
  while (!is.eof())
    {
      is >> tempstr;
      if (tempstr.compare("M") == 0)
	{
	  read_M_line(is,m);
	  if (num_hits > 0)
	    m.num_proteins_in_match.push_back(num_proteins_in_match);
	  num_proteins_in_match = 0;
	  num_hits++;
	}
      if(tempstr.compare("L") == 0)
	{
	  string prot;
	  is >> prot;
	  m.proteins.push_back(prot);
	  num_proteins_in_match++;
	}
      if(tempstr.compare("S") == 0)
	{
	  if (num_hits > 0)
	    m.num_proteins_in_match.push_back(num_proteins_in_match);
	  break;
	}
	
    }
  return num_hits;

  
}

void SQTParser :: read_sqt_file(ifstream &is, string &decoy_prefix, int final_hits, enzyme enz, int pass)
{

  string line;
  string tempstr;
  is >> tempstr;
  while(!is.eof())
    {
      if (tempstr.compare("H") != 0)
	break;
      getline(is, line);
      is >> tempstr;
    }
  int num_hits;
  while(!is.eof())
    {
      assert(tempstr.compare("S") == 0);
      num_hits = parse_sqt_spectrum_matches(is,m);
      if (pass == 1)
	add_matches_to_tables(m, decoy_prefix, num_hits, final_hits);
      else
	extract_features(m, decoy_prefix, num_hits, final_hits,enz);
    }
  

}

void SQTParser :: save_data_in_binary(string out_dir)
{

  ostringstream fname;
  //write out data summary
  fname << out_dir << "/summary.txt";
  ofstream f_summary(fname.str().c_str());
  //psm info
  //psm info
  f_summary << num_features+num_spec_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  //peptide info
  f_summary << num_pep << " " << num_pos_pep << " " << num_neg_pep << endl;
  //protein info
  f_summary << num_prot << " " << num_pos_prot << " " << num_neg_prot << endl;
  f_summary.close();
  fname.str("");
  
  //psmind_to_pepind
  fname << out_dir << "/psmind_to_pepind.txt";
  ofstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  f_psmind_to_pepind.write((char*)psmind_to_pepind,sizeof(int)*num_psm);
  f_psmind_to_pepind.close();
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

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  ofstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  f_psmind_to_label.write((char*)psmind_to_label,sizeof(int)*num_psm);
  f_psmind_to_label.close();
  fname.str("");

  //psmind_to_fname
  fname << out_dir << "/psmind_to_fname.txt";
  ofstream f_psmind_to_fname(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_fname.begin(); it != psmind_to_fname.end(); it++)
    f_psmind_to_fname << it->first << " " << it->second << "\n";
  f_psmind_to_fname.close();
  fname.str("");

  //pepind_to_label
  fname << out_dir << "/pepind_to_label.txt";
  ofstream f_pepind_to_label(fname.str().c_str(),ios::binary);
  f_pepind_to_label.write((char*)pepind_to_label,sizeof(int)*num_pep);
  f_pepind_to_label.close();
  fname.str("");

  //pepind_to_psminds
  fname << out_dir << "/pepind_to_psminds.txt";
  ofstream f_pepind_to_psminds(fname.str().c_str(),ios::binary);
  pepind_to_psminds.save(f_pepind_to_psminds);
  f_pepind_to_psminds.close();
  fname.str("");

  //pepind_to_protinds
  fname << out_dir << "/pepind_to_protinds.txt";
  ofstream f_pepind_to_protinds(fname.str().c_str(),ios::binary);
  pepind_to_protinds.save(f_pepind_to_protinds);
  f_pepind_to_protinds.close();
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

  //protind_to_label
  fname << out_dir << "/protind_to_label.txt";
  ofstream f_protind_to_label(fname.str().c_str(),ios::binary);
  f_protind_to_label.write((char*)protind_to_label,sizeof(int)*num_prot);
  f_protind_to_label.close();
  fname.str("");
  
  //protind_to_num_all_pep
  fname << out_dir << "/protind_to_num_all_pep.txt";
  ofstream f_protind_to_num_all_pep(fname.str().c_str(),ios::binary);
  f_protind_to_num_all_pep.write((char*)protind_to_num_all_pep,sizeof(int)*num_prot);
  f_protind_to_num_all_pep.close();
  fname.str("");


  //protind_to_pepinds
  fname << out_dir << "/protind_to_pepinds.txt";
  ofstream f_protind_to_pepinds(fname.str().c_str(),ios::binary);
  protind_to_pepinds.save(f_protind_to_pepinds);
  f_protind_to_pepinds.close();
  fname.str("");

  //ind_to_prot
  fname << out_dir << "/ind_to_prot.txt";
  ofstream f_ind_to_prot(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = ind_to_prot.begin(); it != ind_to_prot.end(); it++)
    f_ind_to_prot << it->first << " " << it->second << "\n";
  f_ind_to_prot.close();
  fname.str("");

  //prot_to_ind
  fname << out_dir << "/prot_to_ind.txt";
  ofstream f_prot_to_ind(fname.str().c_str(),ios::binary);
  for(map<string,int>::iterator it = prot_to_ind.begin(); it != prot_to_ind.end(); it++)
    f_prot_to_ind << it->first << " " << it->second << "\n";
  f_prot_to_ind.close();
  fname.str("");


}



void SQTParser :: clean_up(string dir)
{

  ostringstream fname;
      
  fname << out_dir << "/summary.txt";
  remove(fname.str().c_str());
  fname.str("");

  fname << dir << "/psm.txt";
  remove(fname.str().c_str());
  fname.str("");
  
  //psmind_to_pepind
  fname << out_dir << "/psmind_to_pepind.txt";
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

  //psmind_to_fname
  fname << out_dir << "/psmind_to_fname.txt";
  remove(fname.str().c_str());
  fname.str("");

  //pepind_to_label
  fname << out_dir << "/pepind_to_label.txt";
  remove(fname.str().c_str());
  fname.str("");

  //pepind_to_psminds
  fname << out_dir << "/pepind_to_psminds.txt";
  remove(fname.str().c_str());
  fname.str("");

  //pepind_to_protinds
  fname << out_dir << "/pepind_to_protinds.txt";
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

  //protind_to_label
  fname << out_dir << "/protind_to_label.txt";
  remove(fname.str().c_str());
  fname.str("");
  
  //protind_to_num_all_pep
  fname << out_dir << "/protind_to_num_all_pep.txt";
  remove(fname.str().c_str());
  fname.str("");


  //protind_to_pepinds
  fname << out_dir << "/protind_to_pepinds.txt";
  remove(fname.str().c_str());
  fname.str("");

  //ind_to_prot
  fname << out_dir << "/ind_to_prot.txt";
  remove(fname.str().c_str());
  fname.str("");

  //prot_to_ind
  fname << out_dir << "/prot_to_ind.txt";
  remove(fname.str().c_str());
  fname.str("");

}


int SQTParser::cntEnzConstraints(string& seq,enzyme enz) {
  int cnt = 0;
  unsigned int pos=0;
  unsigned int pos1 = pos;
  char n = seq.at(pos);
  pos++;
  while (pos<seq.size()-1) {
    char c = seq.at(pos);
    if (isEnz(n,c,enz))
      {
	int pep_len = pos-pos1;
	if(pep_len <= max_len && pep_len >= min_len)
	  {
	    cnt++;
	    pos1 = pos;
	  }
      }
    n=c;
    pos++;
  }
  return cnt;
}

void SQTParser :: digest_database(ifstream &f_db, enzyme e)
{
  string tempstr;
  string prot;
  ostringstream seq;
  int num_prot_read = 0;
  while(!f_db.eof())
    {
      f_db >> tempstr;
      if(tempstr[0] == '>')
	{
	  if(num_prot_read > 0)
	    {
	      string sequence = seq.str();
	      int cnt = cntEnzConstraints(sequence,e);
	      protein_to_num_all_pep_map[prot] = cnt+1;
	      seq.str("");
	    }
	  prot = tempstr.substr(1,tempstr.size());
	  num_prot_read++;
	  
	  getline(f_db, tempstr);
	}
      else
	{
	  seq << tempstr;
	}
    }
  if(num_prot_read > 0)
    {
      string sequence = seq.str();
      int cnt = cntEnzConstraints(sequence,e);
      protein_to_num_all_pep_map[prot] = cnt;
      seq.str("");
    }
  
}


int SQTParser :: run()
{
  //parse database
  for(unsigned int i = 0; i < db_file_names.size(); i++)
    {
      db_name = db_file_names[i];
      ifstream f_db(db_name.c_str());
      if(!f_db.is_open())
	{
	  carp(CARP_WARNING, "could not open database file: %s", db_name.c_str());
	  return 0;
	}
      carp(CARP_INFO,"digesting database %s", db_name.c_str());
      digest_database(f_db, e);
      f_db.close();
    }

  int pass;
  //input sqt file
  int num_files_read = 0;
  for(unsigned int i = 0; i < sqt_file_names.size(); i++)
    {
      cur_fname = sqt_file_names[i];
      carp(CARP_INFO, "parsing file %s", cur_fname.c_str());
      ifstream f_sqt(cur_fname.c_str());
      if(!f_sqt.is_open())
	{
	  carp(CARP_WARNING, "could not open sqt file: %s", cur_fname.c_str());
	  continue;
	}
      //first pass
      pass = 1;
      read_sqt_file(f_sqt, decoy_prefix, fhps,e,pass);
      f_sqt.close();
      num_files_read++;
    }
  if(num_files_read < 1)
    {
      carp(CARP_WARNING, "could not parse any sqt files");
      return 0;
    }

  allocate_feature_space();
  fill_graphs();
  
  ostringstream fname;
  fname << out_dir << "/psm.txt";
  f_psm.open(fname.str().c_str());
  fname.str("");
    
  for(unsigned int i = 0; i < sqt_file_names.size(); i++)
    {
      if(num_spec_features>0)
	{
	  string ms2_fn = ms2_file_names[i];
	  //prepare to generate spectrum features
	  sfg.clear();
	  if(!sfg.open_ms2_file_for_reading(ms2_fn))
	    {
	      carp(CARP_WARNING, "could not open ms2 file %s for reading", ms2_fn.c_str());
	      return 0;
	    }
	  carp(CARP_INFO, "reading file %s", ms2_fn.c_str());
	  sfg.read_ms2_file();
	  sfg.initialize_aa_tables();
	}
      //second pass
      pass = 2;
      cur_fname = sqt_file_names[i];
      carp(CARP_INFO, "extracting features from file %s", cur_fname.c_str()); 
      ifstream f_sqt(cur_fname.c_str());
      read_sqt_file(f_sqt, decoy_prefix, fhps,e,pass);
      f_sqt.close();
      
    }
  f_psm.close();
  
  carp(CARP_INFO, "Number of spectra: %d", num_spectra);
  carp(CARP_INFO, "Number of PSMs: total %d positives %d negatives %d", num_psm, num_pos_psm, num_neg_psm);
  carp(CARP_INFO, "Number of peptides: total %d positives %d negatives %d", num_pep, num_pos_pep, num_neg_pep);
  carp(CARP_INFO, "Number of proteins: total %d positives %d negatives %d", num_prot, num_pos_prot, num_neg_prot);

  //save the data
  save_data_in_binary(out_dir);

  return 1;
}

void SQTParser :: read_list_of_files(string &list, vector<string> &fnames)
{
  ifstream f(list.c_str());
  string str;
  f >> str;
  while(!f.eof())
    {
      fnames.push_back(str);
      f >> str;
    }
  f.close();
}


int SQTParser :: check_files(vector <string> &filenames)
{
  ifstream ftry;
  for(unsigned int i = 0; i < filenames.size(); i++)
    {
      string fname = filenames[i];
      ftry.open(fname.c_str());
      if(!ftry.is_open())
	{
	  cout << "could not open file " << fname << endl;
	  return 0;
	}
      ftry.close();
    }
  return 1;
}

int SQTParser :: set_output_dir(string &output_dir, int overwrite_flag)
{
  int intStat;
  struct stat stFileInfo;
  intStat = stat(output_dir.c_str(), &stFileInfo);

  if (intStat == 0)
    {
      //is this a directory?
      if(!((stFileInfo.st_mode & S_IFMT) == S_IFDIR))
	{
	  //it is not a directory
      	  cout << "WARNING: File " << output_dir << " already exists, but it is not a directory" << endl;
	  if(overwrite_flag == 1)
	    {
	      cout << "INFO: Creating output directory " << output_dir << endl;
	      remove(output_dir.c_str());
	      int dir_access = S_IRWXU + S_IRWXG + S_IRWXO;
	      if (mkdir(output_dir.c_str(), dir_access)) {
		// mkdir failed
		cout << "FATAL: Unable to create output directory " << output_dir << endl;
		return 0;
	      }
	    }
	  else
	    {
	      cout << "FATAL: File " << output_dir << " cannot be overwritten. Please use --overwrite T to replace or specify a different output directory." << endl;
	      return 0;
	    }
	}
    }
  else
    {
      cout << "INFO: creating output directory " << output_dir << endl;
      int dir_access = S_IRWXU + S_IRWXG + S_IRWXO;
      if (mkdir(output_dir.c_str(), dir_access)) {
	// mkdir failed
	cout << "FATAL: unable to create output directory " << output_dir << endl;
	return 0;
      }
    }

  out_dir = output_dir;
  return 1;
}

int SQTParser :: is_fasta(string &fname)
{
  if(fname.find(".fasta") != string::npos)
    return 1;
  if(fname.find(".fsa") != string::npos)
    return 1;
  if(fname.find(".fa") != string::npos)
    return 1;
  return 0;
}

int SQTParser :: set_database_source(string &db_source)
{
  DIR *dp;
  struct dirent *dirp;

  int intStat;
  struct stat stFileInfo;
  intStat = stat(db_source.c_str(), &stFileInfo);
  if(intStat != 0)
    {
      carp(CARP_WARNING, "%s does not exist", db_source.c_str());
      return 0;
    }
  //is this a directory?
  if((stFileInfo.st_mode & S_IFMT) == S_IFDIR)
    {
      //try to open it
      if((dp  = opendir(db_source.c_str())) == NULL)
	{
	  carp(CARP_WARNING, "openning directory %s failed ", db_source.c_str());
	  return 0;
	}
      int cn = 0;
      while ((dirp = readdir(dp)) != NULL) 
	{
	  string fname = string(dirp->d_name);
	  if(is_fasta(fname))
	  //if(fname.find(".fasta") != string::npos)
	    {
	      ostringstream fstr;
	      fstr << db_source;
	      if(db_source.at(db_source.size()-1) != '/')
		fstr <<"/";
	      fstr << fname;
	      string dbname = fstr.str();
	      db_file_names.push_back(dbname);
	      cn++;
	    }
	}
      closedir(dp);
      if(cn<1)
	{
	  carp(CARP_WARNING, "did not find any .fasta files in %s directory", db_source.c_str());
	  return 0;
	}
    }
  else
    {
      if(is_fasta(db_source))
      //if(db_source.find(".fasta") != string :: npos)
	db_file_names.push_back(db_source);
      else 
	read_list_of_files(db_source, db_file_names);
    }
     
  return 1;
}

/*
int SQTParser :: set_input_sources(string &sqt_source, string &ms2_source)
{
  DIR *dp;
  struct dirent *dirp;

  int intStat;
  struct stat stFileInfo;
  intStat = stat(sqt_source.c_str(), &stFileInfo);
  if(intStat != 0)
    {
      carp(CARP_WARNING, "%s does not exist", sqt_source.c_str());
      return 0;
    }
  if((stFileInfo.st_mode & S_IFMT) == S_IFDIR)
    {
      //try openning the directory
      if((dp  = opendir(sqt_source.c_str())) == NULL)
	{
	  carp(CARP_WARNING, "openning directory %s failed ", sqt_source.c_str());
	  return 0;
	}
      int cn = 0;
      //read sqt files in the directory and match them by name to the .ms2 files
      while ((dirp = readdir(dp)) != NULL) 
	{
	  string fname = string(dirp->d_name);
	  size_t pos = fname.find(".sqt"); 
	  if(pos != string::npos)
	    {
	      string prefix = fname.substr(0,pos);
	      ostringstream fstr;
	      fstr << sqt_source <<"/" << fname;
	      string sqtname = fstr.str();
	      fstr.str("");
	      fstr << ms2_source <<"/" << prefix << ".ms2";
	      string ms2name = fstr.str();
	      sqt_file_names.push_back(sqtname);
	      ms2_file_names.push_back(ms2name);
	      //cout << "found " << sqtname <<  " " << ms2name << endl;
	      cn++;
	    }
	}
      closedir(dp);
      if(cn<1)
	{
	  carp(CARP_WARNING, "did not find any .sqt files in %s directory", sqt_source.c_str());
	  return 0;
	}
    }
  else
    {
      if(sqt_source.find(".sqt") != string :: npos)
	{
	  sqt_file_names.push_back(sqt_source);
	  if(ms2_source.find(".ms2") == string::npos)
	    {
	      carp(CARP_WARNING,  "expecting ms2 file to accompany the sqt file");
	      return 0;
	    }
	  ms2_file_names.push_back(ms2_source);
	}
      else 
	{
	  read_list_of_files(sqt_source, sqt_file_names);
	  read_list_of_files(ms2_source, ms2_file_names);
	  if(ms2_file_names.size() != sqt_file_names.size())
	    {
	      carp(CARP_WARNING, " the number of sqt and ms2 files does not match: each sqt file should be accompaned by ms2 file");
	      return 0;
	    }
	}
    }  
    
  return 1;
}
*/

int SQTParser :: match_sqt_to_ms2(string &sqt_source, string &prefix)
{
  DIR *dp;
  struct dirent *dirp;
  //try openning the directory
  if((dp  = opendir(sqt_source.c_str())) == NULL)
    {
      carp(CARP_WARNING, "openning directory %s failed ", sqt_source.c_str());
      return 0;
    }
  int cn = 0;
  //read sqt files in the directory 
  while ((dirp = readdir(dp)) != NULL) 
    {
      string fname = string(dirp->d_name);
      size_t pos = fname.find(".sqt");
      if((pos != string::npos) && (fname.find(prefix) != string::npos))
	{
	  ostringstream fstr;
	  fstr << sqt_source;
	  if(sqt_source.at(sqt_source.size()-1) != '/')
	    fstr << "/";
	  fstr << fname;
	  string sqtname = fstr.str();
	  fstr.str("");
	  sqt_file_names.push_back(sqtname);
	  cout << "found " << sqtname << endl;
	  cn++;
	}
    }
  closedir(dp);
  return cn;
}


int SQTParser :: collect_ms2_files(string &ms2_source, string & sqt_source)
{
  DIR *dp;
  struct dirent *dirp;
  //try openning the directory
  if((dp  = opendir(ms2_source.c_str())) == NULL)
    {
      carp(CARP_WARNING, "openning directory %s failed ", ms2_source.c_str());
      return 0;
    }
  int cn = 0;
  //read ms2 files in the directory 
  while ((dirp = readdir(dp)) != NULL) 
    {
      string fname = string(dirp->d_name);
      size_t pos = fname.find(".ms2");
      if(pos != string::npos)
	{
	  string prefix = fname.substr(0,pos);
	  //collect the file
	  ostringstream fstr;
	  fstr << ms2_source;
	  if(ms2_source.at(ms2_source.size()-1) != '/')
	    fstr << "/";
	  fstr << fname;
	  string ms2name = fstr.str();
	  fstr.str("");
	  int num_matched = match_sqt_to_ms2(sqt_source, prefix); 
	  if(!num_matched)
	    carp(CARP_WARNING, "could not find %s*.sqt in directory %s to match %s, skipping", prefix.c_str(), sqt_source.c_str(), ms2name.c_str());
	  else
	    {
	      for(int i = 0; i < num_matched; i++)
		ms2_file_names.push_back(ms2name);
	    }
	  cn++;
	}
    }
  closedir(dp);

  if(cn<1)
    {
      carp(CARP_WARNING, "did not find any .ms2 files in %s directory", ms2_source.c_str());
      return 0;
    }
  
  return 1;
}


int SQTParser :: set_input_sources(string &ms2_source, string &sqt_source)
{
  
  int intStat;
  struct stat stFileInfo;
  intStat = stat(ms2_source.c_str(), &stFileInfo);
  if(intStat != 0)
    {
      carp(CARP_WARNING, "%s does not exist", ms2_source.c_str());
      return 0;
    }
  if((stFileInfo.st_mode & S_IFMT) == S_IFDIR)
    {
      if(!collect_ms2_files(ms2_source, sqt_source))
      return 0;
    }
  else
    {
      if(ms2_source.find(".ms2") != string :: npos)
	{
	  if(sqt_source.find(".sqt") == string::npos)
	    {
	      carp(CARP_WARNING,  "expecting sqt file to accompany the ms2 file");
	      return 0;
	    }
	  ms2_file_names.push_back(ms2_source);
	  sqt_file_names.push_back(sqt_source);
	}
      else 
	{
	  read_list_of_files(ms2_source, ms2_file_names);
	  read_list_of_files(sqt_source, sqt_file_names);
	  if(ms2_file_names.size() != sqt_file_names.size())
	    {
	      carp(CARP_WARNING, " the number of sqt and ms2 files does not match: each sqt file should be accompaned by ms2 file");
	      return 0;
	    }
	}
    }

  return 1;
}

/*************** for separate searches ********************************************************/

int SQTParser :: match_target_sqt_to_ms2(string &sqt_source, string &prefix)
{
  DIR *dp;
  struct dirent *dirp;
  //try openning the directory
  if((dp  = opendir(sqt_source.c_str())) == NULL)
    {
      carp(CARP_WARNING, "openning directory %s failed ", sqt_source.c_str());
      return 0;
    }
  int cn = 0;
  //read sqt files in the directory 
  while ((dirp = readdir(dp)) != NULL) 
    {
      string fname = string(dirp->d_name);
      size_t pos = fname.find(".target.sqt");
      if((pos != string::npos) && (fname.find(prefix) != string::npos))
	{
	  ostringstream fstr;
	  fstr << sqt_source;
	  if(sqt_source.at(sqt_source.size()-1) != '/')
	    fstr << "/";
	  fstr << fname;
	  string sqtname = fstr.str();
	  fstr.str("");
	  sqt_file_names.push_back(sqtname);
	  cout << "found " << sqtname << endl;
	  cn++;
	}
    }
  closedir(dp);
  return cn;
}


int SQTParser :: match_decoy_sqt_to_ms2(string &sqt_source, string &prefix)
{
  DIR *dp;
  struct dirent *dirp;
  //try openning the directory
  if((dp  = opendir(sqt_source.c_str())) == NULL)
    {
      carp(CARP_WARNING, "openning directory %s failed ", sqt_source.c_str());
      return 0;
    }
  int cn = 0;
  //read sqt files in the directory 
  while ((dirp = readdir(dp)) != NULL) 
    {
      string fname = string(dirp->d_name);
      size_t pos = fname.find(".decoy.sqt");
      if((pos != string::npos) && (fname.find(prefix) != string::npos))
	{
	  ostringstream fstr;
	  fstr << sqt_source;
	  if(sqt_source.at(sqt_source.size()-1) != '/')
	    fstr << "/";
	  fstr << fname;
	  string sqtname = fstr.str();
	  fstr.str("");
	  sqt_file_names.push_back(sqtname);
	  cout << "found " << sqtname << endl;
	  cn++;
	}
    }
  closedir(dp);
  return cn;
}


int SQTParser :: collect_ms2_files(string &ms2_source, string &sqt_target_source, string &sqt_decoy_source)
{
  DIR *dp;
  struct dirent *dirp;
  //try openning the directory
  if((dp  = opendir(ms2_source.c_str())) == NULL)
    {
      carp(CARP_WARNING, "openning directory %s failed ", ms2_source.c_str());
      return 0;
    }
  int cn = 0;
  //read ms2 files in the directory 
  while ((dirp = readdir(dp)) != NULL) 
    {
      string fname = string(dirp->d_name);
      size_t pos = fname.find(".ms2");
      if(pos != string::npos)
	{
	  string prefix = fname.substr(0,pos);
	  //collect the file
	  ostringstream fstr;
	  fstr << ms2_source;
	  if(ms2_source.at(ms2_source.size()-1) != '/')
	    fstr << "/";
	  fstr << fname;
	  string ms2name = fstr.str();
	  fstr.str("");
	  int num_matched_targets = match_target_sqt_to_ms2(sqt_target_source, prefix); 
	  int num_matched_decoys = match_decoy_sqt_to_ms2(sqt_decoy_source, prefix); 
	  	  
	  if(!num_matched_targets)
	    {
	      carp(CARP_WARNING, "could not find %s*.target.sqt in directory %s to match %s, skipping", prefix.c_str(), sqt_target_source.c_str(), ms2name.c_str());
	      continue;
	    }
	  if(!num_matched_decoys)
	    {
	      carp(CARP_WARNING, "could not find %s*.decoy.sqt in directory %s to match %s, skipping", prefix.c_str(), sqt_decoy_source.c_str(), ms2name.c_str());
	      continue;
	    }
	  for(int i = 0; i < (num_matched_targets+num_matched_decoys); i++)
	    ms2_file_names.push_back(ms2name);
	  	  
	  cn++;
	}
    }
  closedir(dp);

  if(cn<1)
    {
      carp(CARP_WARNING, "did not find any .ms2 files in %s directory or did not find any sqt files matching the ms2 files by name", ms2_source.c_str());
      return 0;
    }
  
  return 1;
}


int SQTParser :: set_input_sources(string &ms2_source, string &sqt_target_source, string &sqt_decoy_source)
{
  
  int intStat;
  struct stat stFileInfo;
  intStat = stat(ms2_source.c_str(), &stFileInfo);
  if(intStat != 0)
    {
      carp(CARP_WARNING, "%s does not exist", ms2_source.c_str());
      return 0;
    }
  if((stFileInfo.st_mode & S_IFMT) == S_IFDIR)
    {
      if(!collect_ms2_files(ms2_source, sqt_target_source, sqt_decoy_source))
      return 0;
    }
  else
    {
      if(ms2_source.find(".ms2") != string :: npos)
	{
	  if(sqt_target_source.find(".sqt") == string::npos)
	    {
	      carp(CARP_WARNING,  "expecting target sqt file to accompany the ms2 file");
	      return 0;
	    }
	  if(sqt_decoy_source.find(".sqt") == string::npos)
	    {
	      carp(CARP_WARNING,  "expecting decoy sqt file to accompany the ms2 file and the target sqt file for the separate searches");
	      return 0;
	    }
	  sqt_file_names.push_back(sqt_target_source);
	  ms2_file_names.push_back(ms2_source);
	  sqt_file_names.push_back(sqt_decoy_source);
	  ms2_file_names.push_back(ms2_source);
	}
      else 
	{
	  read_list_of_files(sqt_target_source, sqt_file_names);
	  read_list_of_files(ms2_source, ms2_file_names);
	  read_list_of_files(sqt_decoy_source, sqt_file_names);
	  read_list_of_files(ms2_source, ms2_file_names);
	  if(ms2_file_names.size() != sqt_file_names.size())
	    {
	      carp(CARP_WARNING, " the number of sqt and ms2 files does not match: each sqt file should be accompaned by ms2 file");
	      return 0;
	    }
	}
    }

  return 1;
}



