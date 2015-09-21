#include "CacheableMass.h"
#include "io/carp.h"

void CacheableMass::init() {
  for (size_t idx=0;idx < NUMBER_MASS_TYPES;idx++) {
    mass_calculated_[idx] = false;
  }
}

CacheableMass::CacheableMass() {
  init();
}

CacheableMass::~CacheableMass() {
}


FLOAT_T CacheableMass::getMass(
  MASS_TYPE_T mass_type
  ) {
 
  if (!mass_calculated_[mass_type]) {
    //carp(CARP_DEBUG, "Calculating mass");
    mass_[mass_type] = calcMass(mass_type);
    mass_calculated_[mass_type] = true;
  }

  return mass_[mass_type];

}

FLOAT_T CacheableMass::getMassConst(
  MASS_TYPE_T mass_type
  ) const {

  if (!mass_calculated_[mass_type]) {
    carp(CARP_FATAL, "Internal Error: mass not calculated yet! %d", mass_type);
  }
  return (mass_[mass_type]);
}

void CacheableMass::copy(
  CacheableMass* src,
  CacheableMass* dest) {

  for (size_t idx=0;idx < NUMBER_MASS_TYPES;idx++) {
    if (src->mass_calculated_[idx]) {
      dest->mass_calculated_[idx] = true;
      dest->mass_[idx] = src->mass_[idx];
    } else {
      dest->mass_calculated_[idx] = false;
    }
  }

}

