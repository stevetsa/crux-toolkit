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

/**
 * Remove ambiguous characters from a set of amino acids.  Note that
 * the resulting string may contain the same amino acid twice.
 */
void Enzyme::removeAmbiguousAminos(
  std::string& aminos
) {

  std::string nonAmbiguousAminos;
  for (std::string::iterator myChar = aminos.begin(); myChar < aminos.end(); 
       myChar++) {

    // B = D or N
    if (*myChar == 'B') {
      nonAmbiguousAminos.append("DN");
    }
    // Z = E or Q
    else if (*myChar == 'Z') {
      nonAmbiguousAminos.append("EQ");
    }
    // O, J, U = illegal
    else if ((*myChar == 'O') || (*myChar == 'U') || (*myChar == 'J')) {
      carp(CARP_ERROR, "Illegal amino acid (%c) in enzyme rule.\n", *myChar);
      exit(1);
    }
    // X = any amino acid
    else if (*myChar == 'X') {
      nonAmbiguousAminos.append(allAminos_);
    } else {
      nonAmbiguousAminos.append(1, *myChar);
    }
  }

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
    char myChar = allAminos_[idx];
    if (aminos.find(myChar) == std::string::npos) {
      complementAminos.append(1, myChar);
    }

  }
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
  for (std::string::iterator myChar = allAminos_.begin(); 
       myChar < allAminos_.end();
       myChar++) {
    if (not cleavageMap[*myChar]) {
      cleavageMap['X'] = false;
    }
  }
}


/**
 * Initialize an enzyme object and precompute which pairs of amino
 * acids it cleaves.
 */
void Enzyme::init(
  const char* name
) {
  allAminos_ = "ACDEFGHIKLMNPQRSTVWY";

  std::string enzymeName(name);

  // If the enzyme is named, convert it to the corresponding custom string.
  std::string enzymeRule;
  if (enzymeName == "trypsin"){
    enzymeRule = "[KR]|{P}";
  } else if (enzymeName == "chymotrypsin"){
    enzymeRule = "[FWY]|{P}";
  } else if (enzymeName == "elastase"){
    enzymeRule = "[ALIV]|{P}";
  } else if (enzymeName == "clostripain"){
    enzymeRule = "[R]|[]";
  } else if (enzymeName == "cyanogen-bromide"){
    enzymeRule = "[M]|[]";
  } else if (enzymeName == "iodosobenzoate"){
    enzymeRule = "[W]|[]";
  } else if (enzymeName == "proline-endopeptidase"){
    enzymeRule = "[P]|[]";
  } else if (enzymeName == "staph-protease"){
    enzymeRule = "[E]|[]";
  } else if (enzymeName == "modified-chymotrypsin"){
    enzymeRule = "[FWYL]|{P}";
  } else if (enzymeName == "elastase-trypsin-chymotrypsin"){
    enzymeRule = "[ALIVKRWFY]|{P}";
  } else if (enzymeName == "aspn"){
    enzymeRule = "[]|[D]";
  } else if (enzymeName == "no-enzyme"){
    enzymeRule = "[X]|[X]";
  } else if (enzymeName == "no-cleavage"){
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

  // Extract the two strings.
  std::string precedingAminos = enzymeRule.substr(1, barPosition - 2);
  std::string followingAminos = 
    enzymeRule.substr(barPosition + 2,
		      enzymeRule.length() - (barPosition + 3));

  // Remove ambiguous characters.
  removeAmbiguousAminos(precedingAminos);
  removeAmbiguousAminos(followingAminos);

  // Complement the two strings, if necessary.
  if ((enzymeRule[0] == '{') && (enzymeRule[barPosition - 1] == '}')) {
    complementAminos(precedingAminos);
  } else if (!((enzymeRule[0] == '[') && 
	       (enzymeRule[barPosition - 1] == ']'))) {
    carp(CARP_ERROR, "Failure to parse first half of enzyme rule (%s).\n",
	 enzymeRule.c_str());
    exit(1);
  }
  if ((enzymeRule[barPosition + 1] == '{') && 
      (enzymeRule[enzymeRule.length() - 1] == '}')) {
    complementAminos(followingAminos);
  } else if (!((enzymeRule[barPosition + 1] == '[') && 
	       (enzymeRule[enzymeRule.length() - 1] == ']'))) {
    carp(CARP_ERROR, "Failure to parse second half of enzyme rule (%s).\n",
	 enzymeRule.c_str());
    exit(1);
  }

  // Fill in the mappings.
  for (char myChar = 'A'; myChar <= 'Z'; myChar++) {
    precedingCleavage_[myChar] = false;
    followingCleavage_[myChar] = false;
  }
  for (std::string::iterator myChar = precedingAminos.begin(); 
       myChar < precedingAminos.end();
       myChar++) {
    precedingCleavage_[*myChar] = true;
  }
  for (std::string::iterator myChar = followingAminos.begin(); 
       myChar < followingAminos.end();
       myChar++) {
    followingCleavage_[*myChar] = true;
  }

  // Cope with ambiguity codes in input.
  addAmbiguityCodes(precedingCleavage_);
  addAmbiguityCodes(followingCleavage_);

  // Give the user an update on what we cleave.
  carp(CARP_INFO, "Expanded enzyme rule: %s|%s", 
       precedingAminos.c_str(), followingAminos.c_str());
  std::string siteList = "Cleavage sites:";
  int numSites = 0;
  for (std::string::iterator myChar1 = allAminos_.begin(); 
       myChar1 < allAminos_.end();
       myChar1++) {
    for (std::string::iterator myChar2 = allAminos_.begin(); 
	 myChar2 < allAminos_.end();
	 myChar2++) {
      if (this->cleaves(*myChar1, *myChar2)) {
	siteList.append(" ");
	siteList.append(1, *myChar1);
	siteList.append(1, *myChar2);
	numSites++;
      }
    }
  }
  carp(CARP_DETAILED_INFO, siteList);
  carp(CARP_INFO, "%d possible enymatic cleavage sites.", numSites);
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
