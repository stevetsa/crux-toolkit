/**
 * \file Enzyme.h
 * \brief Decides which pairs of amino acids the digestion enzyme cleaves.
 */

/*
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 19 December 2012
 * $Revision: 1.22 $
 *****************************************************************************/
#ifndef ENZYME_H
#define ENZYME_H

#include <stdio.h>
#include <string>

using namespace Crux;

class Enzyme {

 private:
  /**
   * Initializes an enzyme of the specified type.
   */
  void init(
    const char* enzymeName
  );

 public:

  /**
   * \returns An enzyme of the specified type.
   */
  Enzyme(
    const char* enzymeName
  );

  /**
   * Decides whether the enzyme cleaves between the two specified amino acids.
   * Use 
   */
  bool cleaves(
    char preceding,
    char following
  );

  /**
   * Frees an allocated Enzyme object.
   */
  ~Scorer();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
