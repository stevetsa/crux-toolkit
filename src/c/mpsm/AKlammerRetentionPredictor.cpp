#include "AKlammerRetentionPredictor.h"

#include <iostream>
#include <fstream>

#include "svm.h"
#include "DelimitedFile.h"

static const char* canonical_aminos[] = 
{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
 "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", 
 "X", "Y"};

static const int num_aminos = 25;



using namespace std;


AKlammerRetentionPredictor::AKlammerRetentionPredictor(bool training) {

  cerr<<"AKlammerRetentionPredictor("<<training<<")"<<endl;
  //cerr<<"Getting parameters"<<endl;
  nterm_ = get_boolean_parameter("aklammer-rtime-nterm");
  cterm_ = get_boolean_parameter("aklammer-rtime-cterm");
  tryptic_ = get_boolean_parameter("aklammer-rtime-tryptic");
  spectrum_mass_ = true;
  /*
  twomer_ = get_boolean_parameter("aklammer-rtime-twomer");
  spectrum_ = get_boolean_parameter("aklammer-rtime-spectrum");
  */
  //cerr<<"Initializing"<<endl;
  init(training);
}

AKlammerRetentionPredictor::AKlammerRetentionPredictor() {
  nterm_ = get_boolean_parameter("aklammer-rtime-nterm");
  cterm_ = get_boolean_parameter("aklammer-rtime-cterm");
  tryptic_ = get_boolean_parameter("aklammer-rtime-tryptic");
  spectrum_mass_ = true;
  init(false);
}

void AKlammerRetentionPredictor::init(bool training) {

  if (!training) {
    loadModel();
  }
  feature_names_.clear();
  for (int i=0;i<num_aminos;i++) {
    feature_names_.push_back(canonical_aminos[i]);
  }
  
  if (nterm_) {
    for (int i=0;i<num_aminos;i++) {
      feature_names_.push_back(string("n") + string(canonical_aminos[i]));
    }
  }

  if (cterm_) {
    if (tryptic_) {
      feature_names_.push_back("trypticK");
    }
    for (int i=0;i<num_aminos;i++) {
      feature_names_.push_back(string("c") + string(canonical_aminos[i]));
    }
  }

  if (spectrum_mass_) {
    feature_names_.push_back("spectrum_mass");
  }

  feature_names_.push_back("peptide_length");
  feature_names_.push_back("peptide_mass");
  feature_names_.push_back("charge");  
  data_ = (struct svm_node*)mymalloc(sizeof(struct svm_node)*(feature_names_.size()+1));
  
  for (int idx=0;idx<feature_names_.size();idx++) {
    data_[idx].index = idx;
  }
  data_[feature_names_.size()].index = -1;
  data_[feature_names_.size()].value = 0;
  /*
  if (!training) {
    loadModel();
  }
  */
}

AKlammerRetentionPredictor::~AKlammerRetentionPredictor() {
  
  if (model_ != NULL) {
    cerr<<"freeing model"<<endl;
    //svm_free_and_destroy_model(&model_);
  }
  
  if (data_ != NULL) {
    cerr<<"freeing data"<<endl;
    //free(data_);
  }
  cerr<<"Done ~AKlammerRetentionPredictor()"<<endl;
}

static const char* model_file = "aklammer_rtime.model";
static const char* scale_file = "aklammer_rtime.scales";
void AKlammerRetentionPredictor::loadModel() {
    cerr<<"AKlammerRetentionPredictor.loadModel(): start"<<endl;
    cerr<<"Loading model"<<endl;
    model_ = svm_load_model(model_file);
    if (model_ == NULL) {
      carp(CARP_FATAL,"AKlammerRetentionPredictor can't open model file!");
    }

    cerr<<"Loading scaling parameters"<<endl;
    DelimitedFileReader scaling_parameters(scale_file, true);
    while (scaling_parameters.hasNext()) {
    
      mins.push_back(scaling_parameters.getDouble("min"));
      maxs.push_back(scaling_parameters.getDouble("max"));
      scaling_parameters.next();
    }
    cerr<<"Loaded "<<mins.size()<<"scaling parameters"<<endl;
}



void AKlammerRetentionPredictor::saveModel() {

  svm_save_model(model_file, model_);
  
  DelimitedFile scaling_parameters;
  scaling_parameters.addColumn("feature");
  scaling_parameters.addColumn("min");
  scaling_parameters.addColumn("max");

  for (int idx=0;idx<feature_names_.size();idx++) {
    scaling_parameters.addRow();
    scaling_parameters.setValue<string>("feature",idx,feature_names_[idx]);
    scaling_parameters.setValue<double>("min",idx,mins[idx]);
    scaling_parameters.setValue<double>("max",idx,maxs[idx]);
  }

  scaling_parameters.addRow();
  
  scaling_parameters.setValue<string>("feature",feature_names_.size(),"rtime");
  scaling_parameters.setValue<double>("min",feature_names_.size(),mins.back());
  scaling_parameters.setValue<double>("max",feature_names_.size(),maxs.back());

  scaling_parameters.saveData(scale_file);

}


void AKlammerRetentionPredictor::trainScale(
  struct svm_node** data, 
  double* rtime, 
  int num_examples
) {

  mins.clear();
  maxs.clear();

  for (int idx=0;idx<feature_names_.size()+1;idx++) {
    mins.push_back(1e7);
    maxs.push_back(-1e7);
  }
  
  for (int example_idx = 0;example_idx < num_examples;example_idx++) {
    for(int feature_idx = 0;feature_idx < feature_names_.size();feature_idx++) {

      double current = data[example_idx][feature_idx].value;
      if (current < mins[feature_idx]) {
        mins[feature_idx] = current;
      }
      if (current > maxs[feature_idx]) {
        maxs[feature_idx] = current;
      }

    }
    double current = rtime[example_idx];
    if (current < mins[feature_names_.size()])
      mins[feature_names_.size()] = current;
    if (current > maxs[feature_names_.size()])
      maxs[feature_names_.size()] = current;
  }



}

void AKlammerRetentionPredictor::scale(struct svm_node** features, double* rtimes, int num_examples) {
  for (int example_idx = 0;example_idx < num_examples;example_idx++) {
    scale(features[example_idx],rtimes[example_idx]);
  }
}

void AKlammerRetentionPredictor::scale(struct svm_node* features, double& rtime) {
    

  for (int idx=0;idx<feature_names_.size();idx++) {
    //cerr<<"Scaling "<<idx<<" of "<<feature_names_.size()<<endl;
    double min = mins.at(idx);
    double max = maxs.at(idx);
    double div = max-min;
    if (div > 0) {
      features[idx].value = (features[idx].value - min)/div;
    }
  }

  //cerr<<"scaling rtime"<<endl;
  double min = mins.back();
  double max = maxs.back();
  double div = max-min;
  if (div > 0) {
    rtime = (rtime - min)/div;
  }
  //cerr<<"Done scaling"<<endl;
}

void AKlammerRetentionPredictor::unscale(double& rtime) {

  double min = mins.back();
  double max = maxs.back();
  double div = max-min;
  if (div > 0) {
    rtime = rtime * div + min;
  }

}


double getMSE(double* a, double* b, int n) {

  double ans = 0;
  for (int i = 0;i<n;i++) {

    double temp = b[i] - a[i];
    ans += (temp * temp);

  }

  ans = ans / (double)n;

  return ans;

}

void AKlammerRetentionPredictor::train(DelimitedFileReader& result_data) {
  cerr<<"AKlammerRetentionPredictor.train(): start"<<endl;

  cerr<<"Creating predictor"<<endl;
  AKlammerRetentionPredictor predictor(true);
  cerr<<"Creating svm problem"<<endl;
  struct svm_problem prob;
  struct svm_parameter param;
  
  memset(&prob, 0, sizeof(struct svm_problem));
  memset(&param, 0, sizeof(struct svm_parameter));
  
  param.svm_type = EPSILON_SVR;
  param.kernel_type = RBF;
  param.degree = 3;
  param.gamma = 1/(double)predictor.feature_names_.size();
  param.coef0 = 0;
  param.nu = 0.5;
  param.cache_size = 100;
  param.C = 2;
  param.eps = 1e-3;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;


  int num_examples = result_data.numRows();
  int num_features = predictor.feature_names_.size();
  int data_size = num_examples * num_features;

  cerr<<"Num examples:"<<num_examples<<endl;
  cerr<<"Num features:"<<num_features<<endl;
  cerr<<"total size of data"<<data_size<<endl;

  cerr<<"Allocating memory"<<endl;
  prob.l = num_examples;
  prob.x = (struct svm_node**)mymalloc(sizeof(struct svn_node*) * num_examples);
  prob.y = (double*)mymalloc(sizeof(double) * num_examples);
  
  for (int idx=0;idx<num_examples;idx++) {
    prob.x[idx] = (struct svm_node*)mymalloc(sizeof(struct svm_node)*(num_features+1));
  }

  cerr<<"filling features"<<endl;
  predictor.printFeatureNames(cerr);
  cerr<<endl;
  int example_idx = 0;
  while (result_data.hasNext()) {
    predictor.fillFeatures(prob.x[example_idx], result_data);
    //predictor.printFeatures(cerr,prob.x[example_idx]);
    //cerr<<endl;
    prob.y[example_idx] = result_data.getDouble("eta");

    example_idx++;
    result_data.next();
  }
  cerr<<"training model"<<endl;

  predictor.trainScale(prob.x, prob.y, num_examples);
  predictor.scale(prob.x, prob.y, num_examples);

  svm_check_parameter(&prob, &param);

  double* target = (double*)mymalloc(sizeof(double) * num_examples);
  
  double best_mse = 100000;
  double best_C = -1;
  double best_gamma = -1;

  for (int G_idx=-8;G_idx<3;G_idx++) {

    param.gamma = pow(10,G_idx);

    for (int C_idx=-3;C_idx<6;C_idx++) {

      param.C = pow(10,C_idx);
      svm_cross_validation(&prob, &param, 3, target);

      double mse = getMSE(prob.y,target, num_examples);
      cerr<<"=================="<<endl;
      cerr<<"C:"<<param.C<<" G:"<<param.gamma<<" mse:"<<mse<<endl;
    
 

      if (best_C == -1 || mse < best_mse) {
        best_C = param.C;
        best_gamma = param.gamma;
        best_mse = mse;
      }
    }
  }

    cerr <<"best C:"<<best_C<<" "<<best_mse<<endl;
  cerr<<"best G:"<<best_gamma<<endl;
  param.C = best_C;
  param.gamma = best_gamma;
  free(target);

  predictor.model_ = svm_train(&prob, &param);

  predictor.saveModel();


  AKlammerRetentionPredictor predictor2;

  result_data.reset();

  double sumsq = 0;
  cerr <<"Getting MSE from training data"<<endl;
  ofstream result("obs.v.pred.txt");
  result << "obs\tpred\tdiff\tadiff\tcharge"<<endl;
  while (result_data.hasNext()) {

    double rtime = result_data.getDouble("eta");
    double predict_rtime = predictor2.predict(result_data);

    //cerr<<"obs:" << rtime<<" pred:"<<predict_rtime<<endl;

    double diff = rtime - predict_rtime;
    double adiff = fabs(diff);

    result << rtime << "\t" 
           << predict_rtime << "\t" 
           << diff << "\t" 
           << adiff << "\t" 
           << result_data.getInteger("charge") << endl;

    sumsq += diff * diff;
    result_data.next();
  }

  result.flush();


  cerr<<"sumsq:"<<sumsq<<endl;
  cerr<<"MSE:"<<sumsq/(double)num_examples<<endl;

  cerr<<"saving model"<<endl;
  
  cerr<<"Cleanup"<<endl;
  
  for (int idx=0;idx<num_examples;idx++) {
    free(prob.x[idx]);
  }
  cerr<<"prob.x"<<endl;
  free(prob.x);
  cerr<<"prob.y"<<endl;
  free(prob.y);
  cerr<<"done"<<endl;
}

void AKlammerRetentionPredictor::printFeatureNames(ostream& os) {

  os << feature_names_[0];
  for (int idx=1;idx<feature_names_.size();idx++) {
    os << "\t" << feature_names_[idx]; 
  }
}

void AKlammerRetentionPredictor::printFeatures(
  ostream& os,
  struct svm_node* data
) {

  os << data[0].value;
  for (int idx=1;idx<feature_names_.size();idx++) {
    os << "\t" << data[idx].value;
  }


}

void AKlammerRetentionPredictor::fillFeatures(
  struct svm_node* data,
  Match* match) {

  char* sequence = match->getSequence();
  int charge = match->getCharge();
  int length = strlen(sequence);
  double spectrum_mz = match->getSpectrum()->getPrecursorMz();
 
  fillFeatures(data, sequence, length, charge, spectrum_mz);

  free(sequence);

}

void AKlammerRetentionPredictor::fillFeatures(
  struct svm_node* data,
  DelimitedFileReader& result_data) {

  string sequence = result_data.getString("sequence");
  int charge = result_data.getInteger("charge");
  int N = sequence.length();

  double spectrum_mz = result_data.getDouble("spectrum precursor m/z");
  
  fillFeatures(data, (char*)sequence.c_str(), N, charge, spectrum_mz);


}

void AKlammerRetentionPredictor::fillFeatures(
  struct svm_node* data, 
  char* sequence, 
  int N, 
  int charge,
  double spectrum_mz
) {

    //Set everything to 0.
    for (int idx=0;idx<feature_names_.size();idx++) {
      data[idx].index = idx;
      data[idx].value = 0.0;
    }
    data[feature_names_.size()].index = -1;
    data[feature_names_.size()].value = 0;
    
    //Set counts for amino acids
    int start_index = 0;
    //amino acid counts.
    for (int idx=0;idx<N;idx++) {
      int findex = (int) (sequence[idx] - 'A');
      //cerr <<sequence[idx]<<" "<<findex<<" "<<(findex+start_index)<<endl;
      data[start_index+findex].value++;
      
    }
    start_index += num_aminos;
    //Set count for nterminus amino acid
    if (nterm_) {
      int findex = (int) (sequence[0] - 'A');
      //cerr<<"nterm:"<<sequence[0]<<" "<<findex<<" "<<(start_index+findex)<<endl;
      data[start_index+findex].value = 1;
      start_index += num_aminos;
    }

    //Set count for cterminus amino acid
    if (cterm_) {
      if (tryptic_) {
        //Set boolean if last amino acid is K or R. (K=1,R=0).
        if (sequence[N-1] == 'K') {
          data[start_index].value = 1;
        }
        start_index++;
        //set penultimate amino acid count.
        int findex = (int)(sequence[N-2] - 'A');
        data[start_index+findex].value = 1;
        start_index += num_aminos;
      }  else {
        //not tryptic, set 
        int findex = (int)(sequence[N-1] - 'A');
        data[start_index+findex].value = 1;
        start_index += num_aminos;
      }
    }

    if (spectrum_mass_) {
      data[start_index].value = spectrum_mz;//(spectrum_mz - MASS_PROTON) * (double)charge;
      start_index++;
    }


    data[start_index].value = N; //peptide length
    data[start_index+1].value = Peptide::calcSequenceMass(sequence, MONO); //peptide mono mass
    data[start_index+2].value = charge; //peptide charge
}



FLOAT_T AKlammerRetentionPredictor::predict(struct svm_node* data) {

  double ans = 0.0;
  scale(data, ans);
  ans = svm_predict(model_, data);
  unscale(ans);
  return ans;

}

FLOAT_T AKlammerRetentionPredictor::predict(DelimitedFileReader& result_data) {
  fillFeatures(data_, result_data);
  return predict(data_);
}

FLOAT_T AKlammerRetentionPredictor::predictRTime(Match* match) {
  fillFeatures(data_, match);
  return predict(data_);
}
