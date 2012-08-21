#include "TabDelimParser.h"
#include "carp.h"
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
	num_total_features(0),
	use_quadratic_features(0),
	num_psm(0),
	num_pos_psm(0),
	num_neg_psm(0),
	num_pep(0),
	num_pep_in_all_psms(0),
	curr_ofst(0),
	psmind(0),
	x(0),
	xs(0),
	num_xlink_features(0),
        have_spectra_(false)
{
  
  //num_psm_features
  num_features = 17;
  //num_xlink_features
  num_base_features = 26;
  num_charge_features = 0;


  num_xlink_features = 0;
  //num_spec_features 
  num_spec_features = 0;

  //number of charge features
  charges.clear();

  num_charge_features = 0;
  
  //final_hits_per_spectrum
  fhps = 1;
  //decoy prefix
  decoy_prefix = "rand_";
  //max peptide length to be considered
  max_len = 50;
  //min peptide length to be considered
   min_len = 7;
  
}

void TabDelimParser :: clear()
{

  cerr << "TabDelimParser::clear() - start"<<endl;
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

  cerr << "TabDelimParser::clear() - done"<<endl;
}



TabDelimParser :: ~TabDelimParser()
{
  clear();
}

void TabDelimParser :: set_use_quadratic_features(int use) {
  use_quadratic_features = use;
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
const static int cleavage_type_idx=15;
const static int protein_id_idx=16;
const static int flanking_aa_idx=17;


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
              cerr << "Error parsing sequence field:"<<pep_and_loc<<endl;
              cerr << "Error!"<<endl;
              exit(-1);
            }
          
            
          }
          int loc1,loc2;
          calc_xlink_locations(num_psm, loc1, loc2);
          psmind_to_loc1[num_psm] = loc1;
          psmind_to_loc2[num_psm] = loc2;
	  //proteins
	  if(tokens[16].size() > 0)
	    psmind_to_protein1[num_psm] = tokens[protein_id_idx];
	  //if(tokens[17].size()>0)
	  //psmind_to_protein2[num_psm] = tokens[17];
	  psmind_to_flankingaas[num_psm] = tokens[flanking_aa_idx];
          num_psm++;
	}
      int charge = atoi(tokens[charge_idx].c_str());
      charges.insert(charge);

    }



}

bool TabDelimParser :: isMissedTryptic(string& sequence, int idx) {

  if (idx == (int)(sequence.length()-1)) {
    return false;
  }

  if ((sequence.at(idx) == 'K' || sequence.at(idx) == 'R') && 
      (sequence.at(idx+1) != 'P')) {
    return true;
  }

  return false;

}

int TabDelimParser :: cntMissedCleavagesLinear(int psmind) {

  string sequence = psmind_to_peptide1[psmind];

  int cnt = 0;

  for (size_t idx = 0;idx < sequence.length();idx++) {
    if (isMissedTryptic(sequence,idx)) {
      cnt++;
    }
  }

  return cnt;
}

int TabDelimParser :: cntMissedCleavagesSelfLoop(int psmind) {
  int loc1,loc2;
  int cnt = 0;
  get_xlink_locations(psmind, loc1, loc2);
  string sequence = psmind_to_peptide1[psmind];
  
  if (loc1 == -1 || loc2 == -1) {
    carp(CARP_FATAL, "location %d %d is -1 for sequence %s!",loc1, loc2, sequence.c_str());
  }

  for (size_t idx=0;idx<sequence.length();idx++) {
    int lidx = idx+1;
    if (lidx != loc1 && lidx != loc2 && isMissedTryptic(sequence, idx)) {
      cnt++;
    }
  }

  return cnt;

}

int TabDelimParser :: cntMissedCleavagesDeadLink(int psmind) {
  int cnt = 0;
  int loc1 = psmind_to_loc1[psmind];
  string seq1 = psmind_to_peptide1[psmind];

  for (size_t idx = 0; idx < seq1.length();idx++) {
    int lidx = idx + 1;
    if (lidx != loc1 && isMissedTryptic(seq1, idx)) {
      cnt++;
    }
  }

  //cerr << "dead-link:"<<seq1<<" "<<cnt<<endl;

  return cnt;

}

int TabDelimParser :: cntMissedCleavagesCrossLink(int psmind) {
  int loc1, loc2;
  int cnt = 0;
  get_xlink_locations(psmind, loc1, loc2);
  string seq1 = psmind_to_peptide1[psmind];
  string seq2 = psmind_to_peptide2[psmind];
  if (loc1 == -1 || loc2 == -1) {
    carp(CARP_FATAL, "location %d %d is -1 for sequence %s %s!",loc1, loc2, seq1.c_str(), seq2.c_str());
  }

  for (size_t idx = 0; idx < seq1.length(); idx++) {
    int lidx = idx+1;
    if (lidx != loc1 && isMissedTryptic(seq1, idx)) {
      cnt++;
    }
  }

  for (size_t idx = 0; idx < seq2.length();idx++) {
    int lidx = idx+1;
    if (lidx != loc2 && isMissedTryptic(seq2, idx)) {
      cnt++;
    }
  }

  return cnt;
}


int TabDelimParser :: cntMissedCleavages(int psmind) {

  

  int type = get_peptide_type(psmind);
  int ans = 0;
  switch (type) {
    case XLINKPRODUCT_LINEAR:
      //linear peptide
      ans = cntMissedCleavagesLinear(psmind);
      break;
    case XLINKPRODUCT_SELFLOOP:
      ans = cntMissedCleavagesSelfLoop(psmind);
      break;
    case XLINKPRODUCT_XLINK:
      ans = cntMissedCleavagesCrossLink(psmind);
      break;
    case XLINKPRODUCT_DEADLINK:
      ans = cntMissedCleavagesDeadLink(psmind); 
      break;
    default:
      carp(CARP_FATAL, "Unknown peptide type:%d", type);
  }

  return ans;


}

void TabDelimParser :: get_xlink_locations(int psmind, int &loc1, int &loc2) {

  loc1 = psmind_to_loc1[psmind];
  loc2 = psmind_to_loc2[psmind];
}

void TabDelimParser :: calc_xlink_locations(int psmind, int &loc1, int & loc2) {

  loc1 = -1;
  loc2 = -1;

  if (psmind_to_loc.find(psmind) != psmind_to_loc.end()) {
    string location_string = psmind_to_loc[psmind];
    //cerr << "locating string:"<<location_string<<endl;
    if (location_string.find("()") == string::npos) {
      //it has a xlinker on it.
  
      size_t comma_pos = location_string.find(',');
      size_t rp_pos = location_string.find(')');
      
      
      if (comma_pos != string::npos) {
        //has two locations.
        string sloc1 = location_string.substr(1, comma_pos-1);
        //cerr << "loc1:"<<sloc1<<endl;
        loc1 = atoi(sloc1.c_str());

        string sloc2 = location_string.substr(comma_pos+1, rp_pos-1);
        //cerr << "loc2:"<<sloc2<<endl;
        loc2 = atoi(sloc2.c_str());
        //cerr << "loc2:"<<loc2<<endl;
      } else {
        //has one location
        string sloc1 = location_string.substr(1, rp_pos-1);
        //cerr << "loc1:"<<sloc1<<endl;
        loc1 = atoi(sloc1.c_str());
      }
    } else {
      //check to see if there is a modification indicating a deadlink
      string sequence = psmind_to_peptide1[psmind];
      size_t lb_pos = sequence.find("[");
      if (lb_pos != string::npos) {
        //cerr << "sequence:"<<sequence<<" "<<lb_pos<<endl;
        loc1 = lb_pos;
      }

    }
  }
}

void TabDelimParser :: allocate_feature_space_xlink()
{
  //space for feature vector

  for (set<int>::iterator iter = charges.begin();iter != charges.end();iter++) {
    cerr << "charge:"<<*iter<<endl;
    charge_vec.push_back(*iter);
  }
  num_charge_features = charge_vec.size();


  cerr << "num_base_features:"<<num_base_features<<endl;
  cerr << "num_charge_features:"<<num_charge_features<<endl;
  num_xlink_features = num_base_features + num_charge_features;
  cerr << "num_xlink_features:"<<num_xlink_features<<endl;
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
  psmind_to_neutral_mass = new double[num_psm];


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
      //cout << subtokens[0] <<" " << subtokens[0].length()<<endl;
      ans = subtokens[0].length();
      return subtokens[0].length();
    } else if (subtokens.size() == 3) {
      //it is a crosslinked peptide.
      //cerr << "parsing cross linked peptide:" << subtokens[0] << endl;
      string pep1 = subtokens.at(0).substr(0,subtokens[0].length()-1);
      string pep2 = subtokens.at(1);
      //cerr << pep1 <<" "<<pep1.length()<<endl;
      //cerr << pep2 <<" "<<pep2.length()<<endl;
      ans = pep1.length() + pep2.length();
    } else {
      cerr << "getPeptideLengthSum:error:" << sequence << endl;
      exit(-1);
    }
  }
  return ans;

}

void TabDelimParser::enzTerm(
  int psmind,
  bool& nterm1,
  bool& cterm1,
  bool& nterm2,
  bool& cterm2
  ) {

  nterm1 = false;
  cterm1 = false;
  nterm2 = false;
  cterm2 = false;

  string flankingaas_str = psmind_to_flankingaas[psmind];
  vector<string> flankingaas_vec1;
  string token1 = ";";
  get_tokens(flankingaas_str, flankingaas_vec1, token1);
  
  string sequence = psmind_to_peptide1[psmind];

  //cerr << "sequence:" << sequence;

  vector<string> flankingaas_vec2;

  
  string token2 = ",";
  get_tokens(flankingaas_vec1[0], flankingaas_vec2, token2);
  nterm1 = enzNTerm(sequence, flankingaas_vec2);
  cterm1 = enzCTerm(sequence, flankingaas_vec2);
  if (flankingaas_vec1.size() == 2) {
    sequence = psmind_to_peptide2[psmind];
    //cerr << "," << sequence;

    flankingaas_vec2.clear();
    get_tokens(flankingaas_vec1[1], flankingaas_vec2, token2);
    nterm2 = enzNTerm(sequence, flankingaas_vec2);
    cterm2 = enzCTerm(sequence, flankingaas_vec2);
      
  }

  //cerr << " flankingaas:"<<flankingaas_str;

  //cerr << " nterm1:" << nterm1 <<
  //        " cterm1:" << cterm1 << 
  //        " nterm2:" << nterm2 << 
  //        " cterm2:" << cterm2 <<endl;


}

bool TabDelimParser::enzNTerm(
  string& sequence, 
  vector<string>& flankingaas
  ) {

  if (sequence == "_") {
    return false;
  }

  for (size_t idx = 0 ; idx < flankingaas.size();idx++) {
    char flankingaa = flankingaas[idx][0];
    if (flankingaa == '-') {
      return true;
    }
    if ((flankingaa == 'R' || flankingaa == 'K') && sequence[0] != 'P') {
      return true;
    }
  }
  return false;

}

bool TabDelimParser::enzCTerm(
  string& sequence, 
  vector<string>& flankingaas
  ) {

  if (sequence == "_") {
    return false;
  }

  char last = sequence[sequence.length()-1];

  for (size_t idx = 0; idx < flankingaas.size();idx++) {
    char flankingaa = flankingaas[idx][1];

    if (flankingaa == '-') {
      return true;
    }

    if ((last == 'R' || last == 'K') && flankingaa != 'P') {
      return true;
    }
  }
  return false;

}

int TabDelimParser::get_peptide_length_short(int psmind) {

  string seq1 = psmind_to_peptide1[psmind];
  string seq2 = psmind_to_peptide2[psmind];

  if (seq2 == "_") {
    return seq1.length();
  } else {
    return min(seq1.length(), seq2.length());
  }

}

int TabDelimParser::get_peptide_length_long(int psmind) {

  string seq1 = psmind_to_peptide1[psmind];
  string seq2 = psmind_to_peptide2[psmind];

  if (seq2 == "_") {
    return seq1.length();
  } else {
    return max(seq1.length(), seq2.length());
  }

}

int TabDelimParser::get_peptide_length_sum(int psmind) {

  string seq1 = psmind_to_peptide1[psmind];
  string seq2 = psmind_to_peptide2[psmind];

  if (seq2 == "_") {
    return seq1.length();
  } else {
    return seq1.length() + seq2.length();
  }
}

void TabDelimParser::get_xcorr_short_long(
  int psmind, 
  double& short_xcorr, 
  double& long_xcorr
  ) {

  assert(get_peptide_type(psmind) == XLINKPRODUCT_XLINK);

  short_xcorr = 0;
  long_xcorr = 0;

  if (have_spectra_) {
    string seq1 = psmind_to_peptide1[psmind];
    string seq2 = psmind_to_peptide2[psmind];

    int loc1 = psmind_to_loc1[psmind];
    int loc2 = psmind_to_loc2[psmind];
    int scan = psmind_to_scan[psmind];
    int charge = psmind_to_charge[psmind];
    double xcorr1=0;
    double xcorr2=0;

    sfg.get_xlink_features(scan, charge, seq1, seq2, loc1, loc2, xcorr1, xcorr2);

    if (seq1.length() <= seq2.length()) {
      short_xcorr = xcorr1;
      long_xcorr = xcorr2;
    } else {
      short_xcorr = xcorr2;
      long_xcorr = xcorr1;
    }
  } else {
    carp(CARP_DEBUG, "Don't have spectra");
  }

}

void TabDelimParser::get_flanking_aas(int psmind, bool second, vector<string>& flanking_aas) {

  flanking_aas.clear();
  string flanking_aas_str = psmind_to_flankingaas[psmind];

  vector<string> temp;
  string delim = ";";
  get_tokens(flanking_aas_str, temp, delim);

  string current;

  if (second) {
    if (temp.size() > 1) {
      current = temp[1];
    } else {
      return;
    }         
  } else {
    current = temp[0];
  }

  delim = ",";
  get_tokens(current, flanking_aas, delim);

}


int TabDelimParser::get_num_enzymatic_ends(int psmind) {

  int ans = 0;

  string seq1 = psmind_to_peptide1[psmind];
  vector<string> flanking_aas;
  get_flanking_aas(psmind, false, flanking_aas);

  if (enzNTerm(seq1, flanking_aas)) {
    ans++;
  }
  if (enzCTerm(seq1, flanking_aas)) {
    ans++;
  }

  string seq2 = psmind_to_peptide2[psmind];

  if (seq2 != "_") {
    get_flanking_aas(psmind, true, flanking_aas);
    if (enzNTerm(seq2, flanking_aas)) {
      ans++;
    }
    if (enzCTerm(seq2, flanking_aas)) {
      ans++;
    }
  }

  return ans;
}

bool TabDelimParser::get_nterm1(int psmind) {

  string seq1 = psmind_to_peptide1[psmind];
  vector<string> flanking_aas;
  get_flanking_aas(psmind, false, flanking_aas);

  if (enzNTerm(seq1, flanking_aas)) {
    return true;
  }
  string seq2 = psmind_to_peptide2[psmind];

  if (seq2 != "_") {
    get_flanking_aas(psmind, true, flanking_aas);
    if (enzNTerm(seq2, flanking_aas)) {
      return true;
    }
  }
  return false;
}

bool TabDelimParser::get_cterm1(int psmind) {

  string seq1 = psmind_to_peptide1[psmind];
  vector<string> flanking_aas;
  get_flanking_aas(psmind, false, flanking_aas);

  if (enzCTerm(seq1, flanking_aas)) {
    return true;
  }
  string seq2 = psmind_to_peptide2[psmind];

  if (seq2 != "_") {
    get_flanking_aas(psmind, true, flanking_aas);
    if (enzCTerm(seq2, flanking_aas)) {
      return true;
    }
  }
  return false;

}


bool TabDelimParser::get_nterm2(int psmind) {
  string seq1 = psmind_to_peptide1[psmind];
  vector<string> flanking_aas;
  get_flanking_aas(psmind, false, flanking_aas);

  if (enzNTerm(seq1, flanking_aas)) {
    string seq2 = psmind_to_peptide2[psmind];

    if (seq2 != "_") {
      get_flanking_aas(psmind, true, flanking_aas);
      if (enzNTerm(seq2, flanking_aas)) {
        return true;
      }
    }
  }

  return false;


}

bool TabDelimParser::get_cterm2(int psmind) {
  string seq1 = psmind_to_peptide1[psmind];
  vector<string> flanking_aas;
  get_flanking_aas(psmind, false, flanking_aas);

  if (enzCTerm(seq1, flanking_aas)) {
    string seq2 = psmind_to_peptide2[psmind];

    if (seq2 != "_") {
      get_flanking_aas(psmind, true, flanking_aas);
      if (enzCTerm(seq2, flanking_aas)) {
        return true;
      }
    }
  }

  return false;


}


XLINK_PRODUCT_T TabDelimParser::get_peptide_type(int psmind) {
  if (psmind_to_peptide2[psmind] == "_") {
    //cerr << psmind_to_peptide1[psmind] << " : ";
    if (psmind_to_loc1[psmind] == -1 ) {
      //its a linear
      //cerr << "linear"<<endl;
      return XLINKPRODUCT_LINEAR;
    } else if (psmind_to_loc2[psmind] == -1) {
      //its a dead link
      //cerr <<"dead-link"<<endl;
      return XLINKPRODUCT_DEADLINK;
    } else {
      //its a self loop
      //cerr << "self-loop"<<endl;
      return XLINKPRODUCT_SELFLOOP;
    }

  } else {
    //its a cross-link
    return XLINKPRODUCT_XLINK;

  }
}



XLINK_PRODUCT_T TabDelimParser::get_peptide_type(string& sequence) {
  string delim3 = " ";
  vector<string> subtokens;
  get_tokens(sequence,subtokens,delim3);
/*
  cout << sequence << endl;;
  for(size_t i = 0; i < subtokens.size(); i++) {
    cout<< i << " " <<subtokens[i] << endl;
  }
*/
  XLINK_PRODUCT_T ans = XLINKPRODUCT_UNKNOWN;

  if(subtokens.size() > 0) {
    if (subtokens.size() == 2) {
      //it it either a linear peptide or a self loop.
      //cout << subtokens[0] <<" " << subtokens[0].length()<<endl;
      //ans = subtokens[0].length();
      if (subtokens[1].find(',') == string::npos) {
        if (subtokens[0].find('[') == string::npos) {
          //it is linear.
          //cerr << sequence << " is linear"<<endl;
          ans = XLINKPRODUCT_LINEAR;
        } else {
          //it is a dead-link
          //cerr << sequence << " is dead-link"<<endl;
          ans = XLINKPRODUCT_DEADLINK;
        }
      } else {
        //it is a self-loop
        ans = XLINKPRODUCT_SELFLOOP;
        //cerr << sequence << " is self loop"<<endl;
      }
    } else if (subtokens.size() == 3) {
      //it is a crosslinked peptide.
      ans = XLINKPRODUCT_XLINK;
      //cerr << sequence << " is cross-linked"<<endl;
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
void TabDelimParser :: extract_xlink_features(int psmind, vector<string> & tokens, double *x)
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
  
  XLINK_PRODUCT_T peptide_type = get_peptide_type(tokens[sequence_idx]);

  //xcorr score
  x[0] = atof(tokens[xcorr_idx].c_str());
  //x[0] = atof(tokens[pvalue_idx].c_str());

  x[1] = x[0];
  x[2] = x[0];

  if (peptide_type == XLINKPRODUCT_XLINK) {
    get_xcorr_short_long(psmind, x[1], x[2]);
  }

  x[7] = 0;
  //sp score
  //x[7] = atof(tokens[sp_score_idx].c_str());

 //log rank by Sp
  x[8] = 0;
  //if(atof(tokens[sp_rank_idx].c_str()) > 0)
  //  x[8]=log(atof(tokens[sp_rank_idx].c_str()));

  //matched ions/predicted ions
  x[9] = 0;
  if(atof(tokens[by_total_idx].c_str()) != 0)
    x[9] = atof(tokens[by_matched_idx].c_str())/atof(tokens[by_total_idx].c_str());

  // absolute value of difference between measured and calculated mass
  x[10] = fabs(atof(tokens[spectrum_mass_idx].c_str())-atof(tokens[peptide_mass_idx].c_str()));

  x[11] = peptide_type == XLINKPRODUCT_LINEAR; // linear peptide
  x[12] = peptide_type == XLINKPRODUCT_SELFLOOP; // self-loop
  x[13] = peptide_type == XLINKPRODUCT_XLINK; // cross-link
  x[14] = peptide_type == XLINKPRODUCT_DEADLINK; // dead-link

  /* short peptide length */
  x[15] = get_peptide_length_short(psmind);

  /* long peptide length */
  x[16] = get_peptide_length_long(psmind);

  /* peptide length */
  x[17] = get_peptide_length_sum(psmind);

  /* number of enzymatic ends */
  x[18] = get_num_enzymatic_ends(psmind);

  /* Does at least one of the peptides have an enzymatic cleavage in the N-terminus? */
  x[19] = get_nterm1(psmind);

  /* Does at least one of the peptides have an enzymatic cleavage in the C-terminus? */
  x[20] = get_cterm1(psmind);

  /* Do both peptides have an enzymatic cleavage in the N-terminus? */
  x[21] = get_nterm2(psmind);

  /* Do both peptides have an enzymatic cleavage in the C-terminus? */
  x[22] = get_cterm2(psmind);

  //number of missed-cleavages
  x[23] = cntMissedCleavages(psmind);

  //observed mass
  x[24] = atof(tokens[spectrum_mass_idx].c_str());

  // number of sequence_comparisons
  x[25] = log(atof(tokens[matches_idx].c_str()));

  //charge
  int charge = atoi(tokens[charge_idx].c_str());

  //cerr << "charge:"<<charge<<endl;

  for (size_t idx = 0; idx < charge_vec.size();idx++) {
    if (charge == charge_vec[idx]) {
      x[num_base_features + idx] = 1.0;
    } else {
      x[num_base_features + idx] = 0.0;
    }
  }
/*
  for (int idx = 0; idx < num_xlink_features;idx++) {
    cerr << idx << ":" << x[idx] << " " ;
  }
  cerr << endl;
*/
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

      if (line.find(decoy_prefix) != string::npos) {
        label = -1;
      } else {
        label = 1;
      }


      get_tokens(line, tokens, delim1);
      if(tokens.size() > 1)
	{
	  //get the scan and the charge
	  int scan = atoi(tokens[0].c_str());
	  int charge = atoi(tokens[1].c_str());

          //record charge and scan
	  psmind_to_scan[psmind] = scan;
	  psmind_to_charge[psmind] = charge;

          double neutral_mass = atof(tokens[spectrum_mass_idx].c_str());
	  //extract features
	  extract_xlink_features(psmind, tokens, x);
	  f_psm.write((char*)x, sizeof(double)*num_xlink_features);
	  if(num_spec_features > 0)
	    {
	      
	      //clear out the spec features
	      memset(xs,0,sizeof(double)*num_spec_features);
	      

              string pept1 = "_";
              string pept2 = "_";

              int loc1 = -1;
              int loc2 = -1;

	      //cerr <<" generating spec features for psmind:"<<psmind<<endl;

	      if(psmind_to_peptide1.find(psmind) != psmind_to_peptide1.end())
		{
		  pept1 = psmind_to_peptide1[psmind];
		  //cout << "psmind:" << psmind << "peptide1:" << pept1 << endl;
		}
	      
	      if(psmind_to_peptide2.find(psmind) != psmind_to_peptide2.end())
		{
		  pept2 = psmind_to_peptide2[psmind];
                  //cout << "psmind:" << psmind << "peptide2:" << pept2 << endl;
	        }

              get_xlink_locations(psmind, loc1, loc2);

              
              
              if (have_spectra_ && num_spec_features == 3) {
                sfg.get_spec_features_m3( scan, charge, pept1, pept2, loc1, loc2, xs);
              }

	      f_psm.write((char*)xs, sizeof(double)*num_spec_features);
	      
	    }
	  
	  num_total_features = num_xlink_features+num_spec_features;

	  if(use_quadratic_features)
	    {
	      double *b = new double[1];
	      for(int i = 0; i < num_xlink_features; i++)
		{
		  for(int j = i; j < num_xlink_features; j++)
		    {
		      b[0] = x[i]*x[j];
		      f_psm.write((char*)b, sizeof(double));
		      num_total_features++;
		    }
		}
	      
	      if(num_spec_features > 0)
		{
		  for(int i = 0; i < num_spec_features; i++)
		    {
		      for(int j = i; j < num_spec_features; j++)
			{
			  b[0] = xs[i]*xs[j];
			  f_psm.write((char*)b, sizeof(double));
			  num_total_features++;
			}
		    }
		}
	      
	    }


	  //record charge and scan
	  psmind_to_scan[psmind] = scan;
	  psmind_to_charge[psmind] = charge;
          psmind_to_neutral_mass[psmind] = neutral_mass;
	  psmind_to_label[psmind] = label;
	  if(label == 1)
	    num_pos_psm++;
	  else
	    num_neg_psm++;

	  //augment counters
	  psmind++;
	  if(psmind % 1000 == 0)
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
  //f_summary << num_xlink_features+num_spec_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  f_summary << num_total_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  //cout << num_xlink_features+num_spec_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  cout << num_total_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
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

  //psmind_to_neutralmass
  /*
  fname << out_dir << "/psmind_to_precursor_mass";
  ofstream f_psmind_to_precursor_mass(fname.str().c_str(),ios::binary);
  */


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
      cerr << line << endl;
      if (line.find(decoy_prefix) != string::npos) {
        label = -1;
      } else {
        label = 1;
      }
      /*
        Assumes second file is all decoys
      if(i == 0)
	label = 1;
      else
	label = -1;
      */
      second_pass_xlink(fin,label);
      fin.close();
    }
  f_psm.close();
  cout << psmind << " " << num_pos_psm << " " << num_neg_psm << endl;
  save_data_in_binary_xlink(out_dir);
  return 1;
}

int TabDelimParser :: run_on_xlink(vector<string> &filenames, string &ms2filename, double xlink_mass)
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
    
  //prepare to generate spectrum features
  sfg.clear();
  sfg.set_xlink_mass(xlink_mass);

  have_spectra_ = sfg.open_ms2_file_for_reading(ms2filename);

  if (have_spectra_) {
    sfg.read_ms2_file();
    sfg.initialize_aa_tables();
    num_spec_features = 3;
  }
  
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
      cerr << line << endl;
      if (line.find(decoy_prefix) != string::npos) {
        label = -1;
      } else {
        label = 1;
      }
/*
      if(i == 0)
	label = 1;
      else
	label = -1;
*/
      second_pass_xlink(fin,label);
      fin.close();
    }
  f_psm.close();
  cout << psmind << " " << num_pos_psm << " " << num_neg_psm << endl;
  save_data_in_binary_xlink(out_dir);
  return 1;
}


