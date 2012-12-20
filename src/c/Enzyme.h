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
#include <map>

//using namespace Crux;

class Enzyme {

 private:

  /* Keep track of whether a cleavage if the given amino acid appears before
   * or after the position. */
  std::map<char,bool> precedingCleavage_;
  std::map<char,bool> followingCleavage_;

  /* List of all amino acids. */
  std::string allAminos_;

  /**
   * Remove ambiguous characters from a set of amino acids.  Note that
   * the resulting string may contain the same amino acid twice.
   */
  void removeAmbiguousAminos(
    std::string& aminos
  );

  /**
   * Take the complement of a set of amino acids.
   */
  void complementAminos(
    std::string& aminos
  );

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
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif