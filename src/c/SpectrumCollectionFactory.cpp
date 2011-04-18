/**
 * \file SpectrumCollectionFactory.cpp 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 June 2011
 * \brief Return a SpectrumCollection object of the appropriate
 * derived class.
 */
#include "SpectrumCollectionFactory.h"

/**
 * Instantiates a SpectrumCollection based on the extension of the
 * given file and the use-mstoolkit option.
 */
SpectrumCollection* new_spectrum_collection(const char* filename){

  // for now, just produces one kind of collection
  Ms2SpectrumCollection* collection = new Ms2SpectrumCollection(filename);
  return collection;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
