#ifndef AKLAMMERRETENTIONPREDICTOR_H
#define AKLAMMERRETENTIONPREDICTOR_H

#include "RetentionPredictor.h"
#include "DelimitedFileReader.h"
#include "mass.h"
#include "Match.h"

#include <string>
#include <vector>



class AKlammerRetentionPredictor: public RetentionPredictor {
  protected:
    bool nterm_;
    bool cterm_;
    bool tryptic_;
    bool spectrum_mass_;
    /*
    bool twomer_;
    bool spectrum_;
    */

    struct svm_model* model_;
    struct svm_node* data_;

    std::vector<std::string> feature_names_;
    std::vector<double> mins;
    std::vector<double> maxs;

    

    void fillFeatures(
      struct svm_node* data,
      char* sequence, 
      int N, 
      int charge,
      double spectrum_mz);

    void fillFeatures(struct svm_node* data, Match* match);
    void fillFeatures(struct svm_node* data, DelimitedFileReader& result_data);


    FLOAT_T predict(struct svm_node* data);
    FLOAT_T predict(DelimitedFileReader& result_data);

    void printFeatureNames(std::ostream& os);
    void printFeatures(std::ostream& os, struct svm_node* data);

    void trainScale(struct svm_node** data, double* rtimes, int num_examples);
    void scale(struct svm_node** data, double* rtimes, int num_examples);
    void scale(struct svm_node* data, double& rtime);
    void unscale(double& rtime);
    
    void saveModel();
    void loadModel();


  public:
    AKlammerRetentionPredictor(bool training);
    AKlammerRetentionPredictor();
    void init(bool training=false);

    virtual ~AKlammerRetentionPredictor();

    virtual FLOAT_T predictRTime(Match* match);

    static void train(DelimitedFileReader& result_data);
};



#endif
