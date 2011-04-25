#ifndef AKLAMMERSTATICRETENTIONPREDICTOR_H
#define AKLAMMERSTATICRETENTIONPREDICTOR_H

#include "objects.h"
#include "RetentionPredictor.h"

class AKlammerStaticRetentionPredictor: public RetentionPredictor {

 protected:
  double* weights_composition;
  double* weights_cterm;
  double* weights_nterm;
  double weights_KorR;
  double weights_mass;
  double weights_length;

  double slope;
  double intercept;

 public:
    
  AKlammerStaticRetentionPredictor();

  virtual ~AKlammerStaticRetentionPredictor();

  virtual FLOAT_T predictRTime(MATCH_T* match);
  
  FLOAT_T predictRTimeS(const char* sequence);
  FLOAT_T predictRTimeS(const char* sequence, int N);
  FLOAT_T predictRTimeS(const char* sequence, int N, FLOAT_T mass);

};




#endif
