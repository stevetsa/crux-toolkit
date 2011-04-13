#include "DataSet.h"

Dataset::Dataset() 
  : num_psms(0), num_pos_psms(0), num_neg_psms(0), 
    psmind_to_features((double*)0), 
    psmind_to_label((int*)0),
    psmind_to_pepind((int*)0),
    psmind_to_scan(0),
    psmind_to_charge(0),
    protind_to_label(0),
    protind_to_num_all_pep(0)
{
}

Dataset::~Dataset()
{
  delete[] psmind_to_features;
  delete[] psmind_to_label;
  delete[] psmind_to_pepind;
  delete[] psmind_to_scan;
  delete[] psmind_to_charge;
  delete[] protind_to_label;
  delete[] protind_to_num_all_pep;
}

void Dataset :: load_prot_data()
{
  ostringstream fname;
  fname << in_dir << "/summary.txt";
  ifstream f_summary(fname.str().c_str());
  if(!f_summary.is_open())
    {
      cout << "could not open files for reading data\n";
      return;
    }
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary >> num_pep;
  f_summary >> num_pos_pep;
  f_summary >> num_neg_pep;
  f_summary >> num_prot;
  f_summary >> num_pos_prot;
  f_summary >> num_neg_prot;
  f_summary.close();
  fname.str("");

  //cout << "num proteins in dataset " << num_prot << " pos " << num_pos_prot << " neg " << num_neg_prot << endl;

  //pepind_to_protinds
  fname << in_dir << "/pepind_to_protinds.txt";
  ifstream f_pepind_to_protinds(fname.str().c_str(),ios::binary);
  pepind_to_protinds.load(f_pepind_to_protinds);
  f_pepind_to_protinds.close();
  fname.str("");


  //pepind_to_psminds
  fname << in_dir << "/pepind_to_psminds.txt";
  ifstream f_pepind_to_psminds(fname.str().c_str(),ios::binary);
  pepind_to_psminds.load(f_pepind_to_psminds);
  f_pepind_to_psminds.close();
  fname.str("");


  //protind_to_label
  fname << in_dir << "/protind_to_label.txt";
  ifstream f_protind_to_label(fname.str().c_str(),ios::binary);
  protind_to_label = new int[num_prot];
  f_protind_to_label.read((char*)protind_to_label,sizeof(int)*num_prot);
  f_protind_to_label.close();
  fname.str("");

  //protind_to_num_all_pep
  fname << in_dir << "/protind_to_num_all_pep.txt";
  ifstream f_protind_to_num_all_pep(fname.str().c_str(),ios::binary);
  protind_to_num_all_pep = new int[num_prot];
  f_protind_to_num_all_pep.read((char*)protind_to_num_all_pep,sizeof(int)*num_prot);
  f_protind_to_num_all_pep.close();
  fname.str("");


  //protind_to_pepinds
  fname << in_dir << "/protind_to_pepinds.txt";
  ifstream f_protind_to_pepinds(fname.str().c_str(),ios::binary);
  protind_to_pepinds.load(f_protind_to_pepinds);
  f_protind_to_pepinds.close();
  fname.str("");

  //ind_to_prot
  fname << in_dir << "/ind_to_prot.txt";
  ifstream f_ind_to_prot(fname.str().c_str(),ios::binary);
  int ind;
  
  while(!f_ind_to_prot.eof())
    {
      string prot;
      f_ind_to_prot >> ind;
      f_ind_to_prot >> prot;
      ind_to_prot[ind] = prot;
    }
  f_ind_to_prot.close();
  fname.str("");

}




void Dataset :: load_data()
{
  ostringstream fname;
  fname << in_dir << "/summary.txt";
  ifstream f_summary(fname.str().c_str());
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary.close();
  fname.str("");

  //psm features
  long begin,end;
  fname << in_dir << "/psm.txt";
  ifstream f_psm_feat(fname.str().c_str(),ios::binary);
  if (0) {
  begin = f_psm_feat.tellg();
  f_psm_feat.seekg(0,ios::end);
  end = f_psm_feat.tellg();
  int psm_count = (end-begin)/(num_features*sizeof(double));
  assert(psm_count == num_psms);
  f_psm_feat.seekg(0,ios::beg);
  }
  psmind_to_features = new double[num_psms*num_features];
  f_psm_feat.read((char*)psmind_to_features,sizeof(double)*num_psms*num_features);
  f_psm_feat.close();
  fname.str("");

  //psmind_to_label
  fname << in_dir << "/psmind_to_label.txt";
  ifstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  psmind_to_label = new int[num_psms];
  f_psmind_to_label.read((char*)psmind_to_label,sizeof(int)*num_psms);
  f_psmind_to_label.close();
  fname.str("");
  //psmind_to_pepind
  fname << in_dir << "/psmind_to_pepind.txt";
  ifstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  psmind_to_pepind = new int[num_psms];
  f_psmind_to_pepind.read((char*)psmind_to_pepind,sizeof(int)*num_psms);
  f_psmind_to_pepind.close();
  fname.str("");
  //psmind_to_scan
  fname << in_dir << "/psmind_to_scan.txt";
  ifstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  psmind_to_scan = new int[num_psms];
  f_psmind_to_scan.read((char*)psmind_to_scan,sizeof(int)*num_psms);
  f_psmind_to_scan.close();
  fname.str("");
  //psmind_to_charge
  fname << in_dir << "/psmind_to_charge.txt";
  ifstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  psmind_to_charge = new int[num_psms];
  f_psmind_to_charge.read((char*)psmind_to_charge,sizeof(int)*num_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //ind_to_pep
  fname << in_dir << "/ind_to_pep.txt";
  ifstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  int ind;
  string pep;
  while(!f_ind_to_pep.eof())
    {
      f_ind_to_pep >> ind;
      f_ind_to_pep >> pep;
      ind_to_pep[ind] = pep;
    }
  f_ind_to_pep.close();
  fname.str("");
}


void Dataset :: load_data(string &summary_fn, string &psm_fn)
{
  ostringstream fname;
  fname << in_dir << "/" << summary_fn;
  ifstream f_summary(fname.str().c_str());
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary.close();
  fname.str("");

  //psm features
  long begin,end;
  fname << in_dir << "/" << psm_fn;
  ifstream f_psm_feat(fname.str().c_str(),ios::binary);
  if (0) {
    begin = f_psm_feat.tellg();
    f_psm_feat.seekg(0,ios::end);
    end = f_psm_feat.tellg();
    int psm_count = (end-begin)/(num_features*sizeof(double));
    assert(psm_count == num_psms);
    f_psm_feat.seekg(0,ios::beg);
  }
  psmind_to_features = new double[num_psms*num_features];
  f_psm_feat.read((char*)psmind_to_features,sizeof(double)*num_psms*num_features);
  f_psm_feat.close();
  fname.str("");

  //psmind_to_label
  fname << in_dir << "/psmind_to_label.txt";
  ifstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  psmind_to_label = new int[num_psms];
  f_psmind_to_label.read((char*)psmind_to_label,sizeof(int)*num_psms);
  f_psmind_to_label.close();
  fname.str("");
  //psmind_to_pepind
  fname << in_dir << "/psmind_to_pepind.txt";
  ifstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  psmind_to_pepind = new int[num_psms];
  f_psmind_to_pepind.read((char*)psmind_to_pepind,sizeof(int)*num_psms);
  f_psmind_to_pepind.close();
  fname.str("");
  //psmind_to_scan
  fname << in_dir << "/psmind_to_scan.txt";
  ifstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  psmind_to_scan = new int[num_psms];
  f_psmind_to_scan.read((char*)psmind_to_scan,sizeof(int)*num_psms);
  f_psmind_to_scan.close();
  fname.str("");
  //psmind_to_charge
  fname << in_dir << "/psmind_to_charge.txt";
  ifstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  psmind_to_charge = new int[num_psms];
  f_psmind_to_charge.read((char*)psmind_to_charge,sizeof(int)*num_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //ind_to_pep
  fname << in_dir << "/ind_to_pep.txt";
  ifstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  int ind;
  while(!f_ind_to_pep.eof())
    {
      string pep;
      f_ind_to_pep >> ind;
      f_ind_to_pep >> pep;
      ind_to_pep[ind] = pep;
    }
  f_ind_to_pep.close();
  fname.str("");
}

/****************************************************************************/

void Dataset :: load_psm_data_for_training()
{
  ostringstream fname;
  fname << in_dir << "/summary.txt";
  ifstream f_summary(fname.str().c_str());
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary.close();
  fname.str("");

  //psm features
  long begin,end;
  fname << in_dir << "/psm.txt";
  ifstream f_psm_feat(fname.str().c_str(),ios::binary);
  if (0) {
    begin = f_psm_feat.tellg();
    f_psm_feat.seekg(0,ios::end);
    end = f_psm_feat.tellg();
    int psm_count = (end-begin)/(num_features*sizeof(double));
    assert(psm_count == num_psms);
    f_psm_feat.seekg(0,ios::beg);
  }
  psmind_to_features = new double[num_psms*num_features];
  f_psm_feat.read((char*)psmind_to_features,sizeof(double)*num_psms*num_features);
  f_psm_feat.close();
  fname.str("");

  //psmind_to_label
  fname << in_dir << "/psmind_to_label.txt";
  ifstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  psmind_to_label = new int[num_psms];
  f_psmind_to_label.read((char*)psmind_to_label,sizeof(int)*num_psms);
  f_psmind_to_label.close();
  fname.str("");

}


void Dataset :: load_psm_data_for_reporting_results()
{
  ostringstream fname;
  fname << in_dir << "/summary.txt";
  ifstream f_summary(fname.str().c_str());
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary.close();
  fname.str("");

  //psmind_to_pepind
  fname << in_dir << "/psmind_to_pepind.txt";
  ifstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  psmind_to_pepind = new int[num_psms];
  f_psmind_to_pepind.read((char*)psmind_to_pepind,sizeof(int)*num_psms);
  f_psmind_to_pepind.close();
  fname.str("");
  //psmind_to_scan
  fname << in_dir << "/psmind_to_scan.txt";
  ifstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  psmind_to_scan = new int[num_psms];
  f_psmind_to_scan.read((char*)psmind_to_scan,sizeof(int)*num_psms);
  f_psmind_to_scan.close();
  fname.str("");
  //psmind_to_charge
  fname << in_dir << "/psmind_to_charge.txt";
  ifstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  psmind_to_charge = new int[num_psms];
  f_psmind_to_charge.read((char*)psmind_to_charge,sizeof(int)*num_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //ind_to_pep
  fname << in_dir << "/ind_to_pep.txt";
  ifstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  int ind;
  while(!f_ind_to_pep.eof())
    {
      string pep;
      f_ind_to_pep >> ind;
      f_ind_to_pep >> pep;
      ind_to_pep[ind] = pep;
    }
  f_ind_to_pep.close();
  fname.str("");

  //psmind_to_fname
  fname << in_dir << "/psmind_to_fname.txt";
  ifstream f_psmind_to_fname(fname.str().c_str(),ios::binary);
  int psmind;
  while(!f_psmind_to_fname.eof())
    {
      string filename;
      f_psmind_to_fname >> psmind;
      f_psmind_to_fname >> filename;
      psmind_to_fname[psmind] = filename;
    }
  f_psmind_to_fname.close();
  fname.str("");
}



void Dataset :: load_prot_data_for_training()
{
  ostringstream fname;
  fname << in_dir << "/summary.txt";
  ifstream f_summary(fname.str().c_str());
  if(!f_summary.is_open())
    {
      cout << "could not open files for reading data\n";
      return;
    }
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary >> num_pep;
  f_summary >> num_pos_pep;
  f_summary >> num_neg_pep;
  f_summary >> num_prot;
  f_summary >> num_pos_prot;
  f_summary >> num_neg_prot;
  f_summary.close();
  fname.str("");

    //pepind_to_protinds
  fname << in_dir << "/pepind_to_protinds.txt";
  ifstream f_pepind_to_protinds(fname.str().c_str(),ios::binary);
  pepind_to_protinds.load(f_pepind_to_protinds);
  f_pepind_to_protinds.close();
  fname.str("");


  //pepind_to_psminds
  fname << in_dir << "/pepind_to_psminds.txt";
  ifstream f_pepind_to_psminds(fname.str().c_str(),ios::binary);
  pepind_to_psminds.load(f_pepind_to_psminds);
  f_pepind_to_psminds.close();
  fname.str("");


  //protind_to_label
  fname << in_dir << "/protind_to_label.txt";
  ifstream f_protind_to_label(fname.str().c_str(),ios::binary);
  protind_to_label = new int[num_prot];
  f_protind_to_label.read((char*)protind_to_label,sizeof(int)*num_prot);
  f_protind_to_label.close();
  fname.str("");

  //protind_to_num_all_pep
  fname << in_dir << "/protind_to_num_all_pep.txt";
  ifstream f_protind_to_num_all_pep(fname.str().c_str(),ios::binary);
  protind_to_num_all_pep = new int[num_prot];
  f_protind_to_num_all_pep.read((char*)protind_to_num_all_pep,sizeof(int)*num_prot);
  f_protind_to_num_all_pep.close();
  fname.str("");


  //protind_to_pepinds
  fname << in_dir << "/protind_to_pepinds.txt";
  ifstream f_protind_to_pepinds(fname.str().c_str(),ios::binary);
  protind_to_pepinds.load(f_protind_to_pepinds);
  f_protind_to_pepinds.close();
  fname.str("");

}

void Dataset :: load_prot_data_for_reporting_results()
{
  ostringstream fname;
  fname << in_dir << "/summary.txt";
  ifstream f_summary(fname.str().c_str());
  if(!f_summary.is_open())
    {
      cout << "could not open files for reading data\n";
      return;
    }
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary >> num_pep;
  f_summary >> num_pos_pep;
  f_summary >> num_neg_pep;
  f_summary >> num_prot;
  f_summary >> num_pos_prot;
  f_summary >> num_neg_prot;
  f_summary.close();
  fname.str("");


  //ind_to_prot
  fname << in_dir << "/ind_to_prot.txt";
  ifstream f_ind_to_prot(fname.str().c_str(),ios::binary);
  int ind;
  
  while(!f_ind_to_prot.eof())
    {
      string prot;
      f_ind_to_prot >> ind;
      f_ind_to_prot >> prot;
      ind_to_prot[ind] = prot;
    }
  f_ind_to_prot.close();
  fname.str("");

}

void Dataset :: normalize_psms()
{
  for (int i = 0; i < num_features; i++)
    {
      double mean = 0;
      for (int j = 0; j < num_psms; j++)
	mean += psmind_to_features[num_features*j+i];
      mean /= num_psms;
      
      double std = 0;
      for (int j = 0; j < num_psms; j++)
	{
	  psmind_to_features[num_features*j+i] -= mean;
	  std += psmind_to_features[num_features*j+i]*psmind_to_features[num_features*j+i];
	}
      std = sqrt(std/num_psms);

      for (int j = 0; j < num_psms; j++)
	{
	  if(std > 0)
	    psmind_to_features[num_features*j+i] /= std;
	}

      double sm = 0;
      for (int j = 0; j < num_psms; j++)
	{
	  sm += psmind_to_features[num_features*j+i]*psmind_to_features[num_features*j+i];
	}
      //cout << i << " " << sm/num_psms << endl;
    }
}

/*
int main()
{
  Dataset* d = new Dataset();
  
  d->set_input_dir("yeast");
  d->load_prot_data();
  //d->normalize_psms();
  return 0;
}
*/
