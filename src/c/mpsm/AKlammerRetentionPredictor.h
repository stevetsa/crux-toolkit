#ifndef AKLAMMERRETENTIONPREDICTOR_H
#define AKLAMMERRETENTIONPREDICTOR_H

#include "RetentionPredictor.h"

#include "match.h"

#include <map>


static const char* canonical_aminos = 
{"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
 "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", 
 "X", "Y"};

static const int num_aminos = 25;

class AKlammerRetentionPredictor: public RetentionPredictor {
  protected:
    bool nterm_;
    bool cterm_;
    bool tryptic_;

    /*
    bool twomer_;
    bool spectrum_;
    */

    struct svm_model* model_;
    struct svm_node* data_;

    vector<string> feature_names_;
        

    void fillFeatures(char* sequence, int N, int charge) {
    }

  public:
    AKlammerRetentionPredictor();
    void init();

    virtual ~AKlammerRetentionPredictor();

    virtual FLOAT_T predictRTime(MATCH_T* match);

    static void train(DelimitedFileReader& result_data);
};



#endif
