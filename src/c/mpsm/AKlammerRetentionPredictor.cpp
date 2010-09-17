#include "AKlammerRetentionPredictor.h"

using namespace std;


AKlammerRetentionPredictor::AKlammerRetentionPredictor(bool training) {
  nterm_ = get_boolean_parameter("aklammer-rtime-nterm");
  cterm_ = get_boolean_parameter("aklammer-rtime-cterm");
  tryptic_ = get_boolean_parameter("aklammer-rtime-tryptic");
  /*
  twomer_ = get_boolean_parameter("aklammer-rtime-twomer");
  spectrum_ = get_boolean_parameter("aklammer-rtime-spectrum");
  */

  init();
  if (!training) {
    model = svm_load_model("aklammer_rtime.model");
  }

}

AKlammerRetentionPredictor::AKlammerRetentionPredictor() : AKlammerRetentionPredictor(false) {
}

AKlammerRetentionPredictor::init() {

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
      feature_names.push_back("trypticK");
    }
    for (int i=0;i<num_aminos;i++) {
      feature_names_.push_back(string("c") + string(canonical_aminos[i]));
    }
  }

  features_names_.push_back("peptide_length");
  feature_names_.push_back("peptide_mass");
  feature_names_.push_back("charge");  
  data_ = mymalloc(sizeof(struct svn_node)*feature_names_.size());
  
  for (int idx=0;idx<feature_names_.size();idx++) {
    data_[idx].index = idx;
  }

}

AKlammerRetentionPredictor::~AKlammerRetentionPredictor() {
  
  if (model_ != NULL) {
    svm_free_and_destroy_model(&model_);
  }
  
  if (data_ != NULL) {
    free(data_);
  }

}


void AKlammerRetentionPredictor::train(DelimitedFileReader& result_data) {


  AKlammerRetentionPredictor predictor(true);
  
  struct svm_problem prob;
  struct svm_parameter param;
  
  memset(&prob, 0, sizeof(struct svm_problem));
  memset(&param, 0, sizeof(struct svm_parameter));
  
  param.svm_type = EPSILON_SVR;
  param.kernel_type = RBF;
  param.degree = 3;
  param.gamma = 1e-7;
  param.coef0 = 0;
  param.nu = 0.5
  param.cache_size = 100;
  param.C = 1e6;
  param.eps = 0.1;
  param.p = 0.1;
  param.shrinking = 1;
  param.probability = 0;
  param.nr_weight = 0;
  param.weight_label = NULL;
  param.weight = NULL;


  int num_examples = result_data.numRows();
  int num_features = feature_names_.size();
  int data_size = num_examples * num_features;

  prob.l = num_examples;
  prob.x = (struct svn_node**)mymalloc(sizeof(struct svn_node*) * num_examples);
  prob.y = mymalloc(sizeof(double) * num_examples);
  
  for (int idx=0;idx<num_examples;idx++) {
    prob.x[idx] = mymalloc(sizeof(struct svn_node) * num_features);
  }

  int example_index = 0;
  while (result_data.hasNext()) {
      
    string sequence = result_data.getString("sequence");
    int charge = result_data.getInteger("charge");
    int N = sequence.length();

    fillFeatures(
      prob.x[example_index],
      sequence.c_str(),
      N,
      charge);

    prob.y[example_idx] = 
      result_data.getDouble("eta"); //This is the spectrum retention time.

    example_index++;
    result_data.next();
  }

  predictor.model = svm_train(&prob, &param);
  
  svm_save_model("aklammer_rtime.model", predictor.model);

  for (int idx=0;idx<num_examples;idx++) {
    free(prob.x[idx]);
  }
  free(prob.x);
  free(prob.y);

}

void AKlammerRetentionPredictor::printFeatureNames(ostream& os) {

  os << feature_names_[0];
  for (int idx=1;idx<feature_names_.size();idx++) {
    os << "\t" << feature_names_[idx]; 
  }
}

void AKlammerRetentionPredictor::printFeatures(
  ostream& os,
  struct svn_node* data
) {

  os << data[0].value;
  for (int idx=1;idx<feature_names_.size();idx++) {
    os << "\t" << data[idx].value;
  }


}

void AKlammerRetentionPredictor::fillFeatures(
  struct svn_node* data, 
  char* sequence, 
  int N, 
  int charge
) {


    //Set everything to 0.
    for (int idx=0;i<features_names_.size();idx++) {
      data[idx].value = 0.0;
    }

    //Set counts for amino acids
    int start_index = 0;
    //amino acid counts.
    for (int idx=0;idx<N;idx++) {
      int findex = (int) (sequence[idx] - 'A');
      data_[start_index+findex].value++;
      start_index += num_aminos;
    }
  
    //Set count for nterminus amino acid
    if (nterm_) {
      int findex = (int) (sequence[idx] - 'A');
      data_[start_index+findex].value = 1;
      start_index += num_aminos;
    }

    //Set count for cterminus amino acid
    if (cterm_) {
      if (tryptic_) {
        //Set boolean if last amino acid is K or R. (K=1,R=0).
        if (sequence[N-1] == 'K') {
          data_[start_index+idx].value = 1;
        }
        start_index++;
        //set penultimate amino acid count.
        int findex = (int)(sequence[N-2] - 'A');
        data_[start_index+findex].value = 1;
        start_index += num_aminos;
      }  else {
        //not tryptic, set 
        int findex = (int)(sequence[N-1] - 'A');
        data_[start_index+findex].value = 1;
        start_index += num_aminos;
    }
    
    data_[start_index] = N; //peptide length
    data_[start_index+1] = calc_peptide_mass(sequence, MONO); //peptide mono mass
    data_[start_index+2] = charge; //peptide charge
    

}

FLOAT_T AKlammerRetentionPredictor::predict(char* sequence, int length, int charge) {

  fillFeatures(data_, sequence, length, charge);
  cout <<sequence<<"\t"<<length<<"\t"<<charge<<"\t";
  printFeatures(cout, data_);

  double ans = svm_predict(model, data_);
  cout <<"\t"<<ans<<endl;

  return ans;
  


}


FLOAT_T AKlammerRetentionPredictor::predictRTime(MATCH_T* match) {


  char* sequence = get_match_sequence(match);
  int charge = get_match_charge(match);
  int length = strlen(sequence);

  FLOAT_T ans = predict(sequence, length, charge);
  
  free(sequence);

  return ans;
}
