/**
 * \file SpectrumCollectionFactory.h 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 June 2011
 * \brief Return a SpectrumCollection object of the appropriate
 * derived class.
 */
#ifndef SPECTRUM_COLLECTION_FACTORY_H
#define SPECTRUM_COLLECTION_FACTORY_H

#include "SpectrumCollection.h"
#include "Ms2SpectrumCollection.h"

/**
 * Instantiates a SpectrumCollection based on the extension of the
 * given file and the use-mstoolkit option.
 */
SpectrumCollection* new_spectrum_collection(const char* filename);


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif 
