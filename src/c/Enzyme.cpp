/**
 * \file Enzyme.cpp
 * \brief Decides which pairs of amino acids the digestion enzyme cleaves.
 */

/*
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 19 December 2012
 * $Revision: 1.22 $
 *****************************************************************************/
#include "Enzyme.h"
#include "carp.h"
#include <string.h> // For the strcmp function.

//using namespace Crux;

/**
 * Remove ambiguous characters from a set of amino acids.  Note that
 * the resulting string may contain the same amino acid twice.
 */
void Enzyme::removeAmbiguousAminos(
  std::string& aminos
) {

  std::string nonAmbiguousAminos;
  for (int idx = 0; idx < aminos.size(); idx++) {
    char myChar = aminos[myChar];

    // B = D or N
    if (myChar == 'B') {
      nonAmbiguousAminos.append("DN");
    }
    // Z = E or Q
    else if (myChar == 'Z') {
      nonAmbiguousAminos.append("EQ");
    }
    // O, J, U = illegal
    else if ((myChar == 'O') || (myChar == 'U') || (myChar == 'J')) {
      carp(CARP_ERROR, "Illegal amino acid (%c) in enzyme rule.\n", myChar);
      exit(1);
    }
    // X = any amino acid
    else if (myChar == 'X') {
      nonAmbiguousAminos.append(allAminos_);
    } else {
      // FIXME: This doesn't work!!
      nonAmbiguousAminos.append(&myChar);
    }
  }

  carp(CARP_WARNING, "Remove ambiguous: %s -> %s", aminos.c_str(),
       nonAmbiguousAminos.c_str());
  aminos = nonAmbiguousAminos;
}

/**
 * Take the complement of a set of amino acids.
 */
void Enzyme::complementAminos(
  std::string& aminos
) {

  std::string complementAminos;
  for (int idx = 0; idx < allAminos_.size(); idx++) {
    char myChar = allAminos_[myChar];
    if (aminos.find(myChar) == std::string::npos) {
      complementAminos.append(&myChar);
    }
  }
  carp(CARP_WARNING, "Complement: %s -> %s.", aminos.c_str(), 
       complementAminos.c_str());
  aminos = complementAminos;
}

/**
 * Given a Boolean mapping without ambiguity codes, compute values
 * for ambiguity codes.
 */
void Enzyme::addAmbiguityCodes(
  std::map<char,bool> cleavageMap
) {

  if (cleavageMap['D'] && cleavageMap['N']) {
    cleavageMap['B'] = true;
  }

  if (cleavageMap['E'] && cleavageMap['Q']) {
    cleavageMap['Z'] = true;
  }

  cleavageMap['X'] = true;
  for (int idx = 0; idx < allAminos_.size(); idx++) {
    char myChar = allAminos_[myChar];
    if (not cleavageMap[myChar]) {
      cleavageMap['X'] = false;
    }
  }
}


/**
 * Initialize an enzyme object and precompute which pairs of amino
 * acids it cleaves.
 */
void Enzyme::init(
  const char* enzymeName
) {
  allAminos_ = "ACDEFGHIKLMNPQRSTVWY";

  // If the enzyme is named, convert it to the corresponding custom string.
  std::string enzymeRule;
  if (strcmp(enzymeName, "trypsin") == 0){
    enzymeRule = "[KR]|{P}";
  } else if (strcmp(enzymeName, "chymotrypsin") == 0){
    enzymeRule = "[FWY]|{P}";
  } else if (strcmp(enzymeName, "elastase") == 0){
    enzymeRule = "[ALIV]|{P}";
  } else if (strcmp(enzymeName, "clostripain") == 0){
    enzymeRule = "[R]|[]";
  } else if (strcmp(enzymeName, "cyanogen-bromide") == 0){
    enzymeRule = "[M]|[]";
  } else if (strcmp(enzymeName, "iodosobenzoate") == 0){
    enzymeRule = "[W]|[]";
  } else if (strcmp(enzymeName, "proline-endopeptidase") == 0){
    enzymeRule = "[P]|[]";
  } else if (strcmp(enzymeName, "staph-protease") == 0){
    enzymeRule = "[E]|[]";
  } else if (strcmp(enzymeName, "modified-chymotrypsin") == 0){
    enzymeRule = "[FWYL]|{P}";
  } else if (strcmp(enzymeName, "elastase-trypsin-chymotrypsin") == 0){
    enzymeRule = "[ALIVKRWFY]|{P}";
  } else if (strcmp(enzymeName, "aspn") == 0){
    enzymeRule = "[]|[D]";
  } else if (strcmp(enzymeName, "no-enzyme") == 0){
    enzymeRule = "[X]|[X]";
  } else if (strcmp(enzymeName, "no-cleavage") == 0){
    enzymeRule = "[]|[]";
  } else {
    enzymeRule = enzymeName;
  }

  // Find the vertical bar.
  int barPosition = -1;
  for (int idx = 0; idx < enzymeRule.length(); idx++) {
    if (enzymeRule[idx] == '|') {
      if (barPosition == -1) {
        barPosition = idx;
      } else {
	carp(CARP_ERROR, "Two vertical bars in enzyme rule (%s).\n",
	     enzymeRule.c_str());
	exit(1);
      }
    }
  }
  if (barPosition == -1) {
    carp(CARP_ERROR, "Unrecognized enzyme name (%s).\n", enzymeRule.c_str());
    exit(1);
  }
  carp(CARP_ERROR, "Vertical bar at %d.", barPosition);

  // Extract the two strings.
  std::string precedingAminos = enzymeRule.substr(1, barPosition - 2);
  std::string followingAminos = 
    enzymeRule.substr(barPosition + 2,
		      enzymeRule.length() - (barPosition + 3));
  carp(CARP_WARNING, "Rule=%s|%s", precedingAminos.c_str(), 
       followingAminos.c_str());

  // Remove ambiguous characters.
  removeAmbiguousAminos(precedingAminos);
  removeAmbiguousAminos(followingAminos);

  // Complement the two strings, if necessary.
  if ((enzymeRule[0] == '[') && (enzymeRule[barPosition - 1] == ']')) {
    complementAminos(precedingAminos);
  } else if (!((enzymeRule[0] == '{') && 
	       (enzymeRule[barPosition - 1] == '}'))) {
    carp(CARP_ERROR, "Failure to parse first half of enzyme rule (%s).\n",
	 enzymeRule.c_str());
    exit(1);
  }
  if ((enzymeRule[barPosition + 1] == '[') && 
      (enzymeRule[enzymeRule.length() - 1] == ']')) {
    complementAminos(followingAminos);
  } else if (!((enzymeRule[barPosition + 1] == '{') && 
	       (enzymeRule[enzymeRule.length() - 1] == '}'))) {
    carp(CARP_ERROR, "Failure to parse second half of enzyme rule (%s).\n",
	 enzymeRule.c_str());
    exit(1);
  }

  // Fill in the mappings.
  for (char myChar = 'A'; myChar <= 'Z'; myChar++) {
    precedingCleavage_[myChar] = false;
    followingCleavage_[myChar] = false;
  }
  for (int idx = 0; idx < precedingAminos.size(); idx++) {
    char myChar = precedingAminos[myChar];
    precedingCleavage_[myChar] = true;
  }
  for (int idx = 0; idx < followingAminos.size(); idx++) {
    char myChar = followingAminos[myChar];
    followingCleavage_[myChar] = true;
  }

  // Cope with ambiguity codes in input.
  addAmbiguityCodes(precedingCleavage_);
  addAmbiguityCodes(followingCleavage_);

  // Give the user an update on what we cleave.
  std::string message = "Enzyme rule: ";
  message.append(enzymeRule);
  message.append("\nCleavage sites: ");
  for (int first_idx = 0; first_idx < allAminos_.size(); first_idx++) {
    char firstChar = allAminos_[first_idx];
    for (int second_idx = 0; second_idx < allAminos_.size(); second_idx++) {
      char secondChar = allAminos_[second_idx];
      if (this->cleaves(firstChar, secondChar)) {
	message.append(&firstChar);
	message.append(&secondChar);
	message.append(" ");
      }
    }
  }
  carp(CARP_WARNING, message);
}

/**
 * \returns An enzyme of the specified type.
 */
Enzyme::Enzyme(
  const char* enzymeName
) {
  init(enzymeName);
}

/**  
 * Decides whether the enzyme cleaves between the two specified amino acids.
 * Use 
 */
bool Enzyme::cleaves(
  char preceding,
  char following
) {
  return(precedingCleavage_[preceding] && followingCleavage_[following]);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
