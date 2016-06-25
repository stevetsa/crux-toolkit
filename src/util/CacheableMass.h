#ifndef CACHEABLEMASS_H_
#define CACHEABLEMASS_H_

#include "objects.h"

class CacheableMass {

 protected:
  bool mass_calculated_[NUMBER_MASS_TYPES]; ///< is mass calculated?
  FLOAT_T mass_[NUMBER_MASS_TYPES]; ///< calculated mass

 public:
  void init();
  CacheableMass();
  virtual ~CacheableMass();

  // Override this
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type
    )=0;

  FLOAT_T getMass(
    MASS_TYPE_T mass_type = MONO
    );

  FLOAT_T getMassConst(
    MASS_TYPE_T mass_type = MONO
  ) const;

  static void copy(CacheableMass* src, CacheableMass* dest);
 
};

#endif
