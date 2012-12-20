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

void Enzyme::init(
  const char* enzymeName
) {

  // If the enzyme is named, convert it to the corresponding custom string.
  // FIXME: Change rules to be appropriate per enzyme.
  string enzymeRule;
  if (strcmp(enzymeName, "trypsin") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "chymotrypsin") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "elastase") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "clostripain") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "cyanogen-bromide") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "iodosobenzoate") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "proline-endopeptidase") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "staph-protease") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "modified-chymotrypsin") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "elastase-trypsin-chymotrypsin") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "no-enzyme") == 0){
    enzymeRule = "[KR]|[X]";
  } else if (strcmp(enzymeName, "no-cleavage") == 0){
    enzymeRule = "[KR]|[X]";
  } else {
    enzymeRule = enzymeName;
  }

  // Find the vertical bar.
  barPosition = -1;
  for (int idx = 0; idx < enzymeRule.len(); idx++) {
    if (enzymeRule[idx] == "|") {
      if (barPosition == -1) {
        barPosition = idx;
      } else {
	carp(CARP_ERROR, "Two vertical bars in enzyme rule (%s).\n",
	     enzymeRule);
	exit(1);
      }
    }
  }
  if (barPosition == -1) {
    carp(CARP_ERROR, "Unrecognized enzyme name (%s).\n", enzymeRule);
  }

  // Extract the two strings.
  bool precedingComplement;
  if ((enzymeRule[0] == "[") && (enzymeRule[barPosition - 1] == "]")) {
    precedingComplement = false;
  } else if ((enzymeRule[0] == "{") && (enzymeRule[barPosition - 1] == "}")) {
    precedingComplement = true;
  } else {
    carp(CARP_ERROR, "Failure to parse enzyme rule (%s).\n", enzymeRule);
  }
  string precedingAminos = enzymeRule.substr(1, barPosition - 2);

  bool followingComplement;
  if ((enzymeRule[barPosition + 1] == "[") && 
      (enzymeRule[enzymeRule.length() - 1] == "]")) {
    followingComplement = false;
  } else if ((enzymeRule[barPosition + 1] == "{") && 
	     (enzymeRule[enzymeRule.length() - 1] == "}")) {
    followingComplement = true;
  } else {
    carp(CARP_ERROR, "Failure to parse enzyme rule (%s).\n", enzymeRule);
  }
  string followingAminos 
    = enzymeRule.substr(barPosition + 1,
			enzymeRule.length() - barPosition - 2);
  

  // Complement the two strings, if necessary.

  // Initialize two mappings from amino acid to false's.
  map<char,bool> precedingCleavage;
  map<char,bool> followingCleavage;
  for (myChar = 'A'; myChar <= 'Z'; myChar++) {
    precedingCleavage(myChar) = false;
    followingCleavage(myChar) = false;
  }

  // Fill in the mappings.
  for (myChar = precedingString.begin(); 
       myChar < precedingString.end(); myChar++) {
    precedingCleavage(myChar) = true;
  }
  for (myChar = followingString.begin(); 
       myChar < followingString.end(); myChar++) {
    followingCleavage(myChar) = true;
  }

  /*
   * Read the value enzyme cleavage string and enter values into local
   * variables.  Correct syntax is [A-Z]|[A-Z] or {A-Z}|{A-Z}.  An X
   * indicates that any residue is legal. Sets pre/post_list size and
   * allocates memory for pre/post_cleavage_list.  Sets
   * pre/post_for_inclusion as true if [] encloses list or false if {}
   * encloses list.  For special case of [X], set p_cleavage_list as
   * empty and inclusion as false.
   */
void parse_custom_enzyme(const char* enzymeRule){

  bool success = true;
  int len = strlen(enzymeRule);
  int idx = 0;
  int pipe_idx = 0;

  // 1. find the |
  for(int idx = 0; idx < len; idx++){
    if( enzymeRule[idx] == '|' ){
      pipe_idx = idx;
      break;
    }
  }
  // check that there isn't a second
  for(int idx = idx+1; idx < len; idx++){
    if( enzymeRule[idx] == '|' ){
      success = false;      
      break;
    }
  }


  // 2. set beginning and end of strings relative to pipe, start, end
  //    0 1    p-1 p p+1 p+2     len-1 len
  //    [ X    ]   | [   X       ]     '0'
  int pre_first_idx = 1;
  int pre_end_idx = pipe_idx - 1;
  int post_first_idx = pipe_idx + 2;
  int post_end_idx = len -1;

  // 3. check that braces match and set inclusion
  // pre-list
  if(pipe_idx < 1){
    success = false;
  }else if(enzymeRule[pre_first_idx-1] == '[' && 
           enzymeRule[pre_end_idx] == ']'){
    pre_for_inclusion = true;
  }else if(enzymeRule[pre_first_idx-1] == '{' && 
           enzymeRule[pre_end_idx] == '}'){
    pre_for_inclusion = false;
  }else{
    success = false;
  }

  // post list
  if(pipe_idx + 2 >= len ){
    success = false;
  }else if(enzymeRule[post_first_idx-1] == '[' && 
           enzymeRule[post_end_idx] == ']'){
    post_for_inclusion = true;
  }else if(enzymeRule[post_first_idx-1] == '{' && 
           enzymeRule[post_end_idx] == '}'){
    post_for_inclusion = false;
  }else{
    success = false;
  }

  // check that braces aren't empty 
  if(pre_first_idx >= pre_end_idx || post_first_idx >= post_end_idx ){
    success = false;
  }

  if( success == false ){
    carp(CARP_FATAL, "Custom enzyme syntax '%s' is incorrect.  "
         "Must be of the form [AZ]|[AZ] or with [] replaced by {}. "
         "AZ is a list of residues (letters A-Z) required [] or prohibited {}. "
         "Use [X] to indicate any reside is legal.",
         enzymeRule);
  }

  // 4. allocate lists and fill
  pre_list_size = pre_end_idx - pre_first_idx;
  pre_cleavage_list = (char*)mycalloc(pre_list_size, sizeof(char));
  for(idx = 0; idx < pre_list_size; idx++){
    pre_cleavage_list[idx] = enzymeRule[pre_first_idx+idx];
  }

  post_list_size = post_end_idx - post_first_idx;
  post_cleavage_list = (char*)mycalloc(post_list_size, sizeof(char));
  for(idx = 0; idx < post_list_size; idx++){
    post_cleavage_list[idx] = enzymeRule[post_first_idx+idx];
  }


  // 5. check special case of [X]
  if(strncmp( enzymeRule, "[X]", pre_list_size+2) == 0){
    free(pre_cleavage_list);
    pre_cleavage_list = NULL;
    pre_list_size = 0;
    pre_for_inclusion = false;
  }

  if(strncmp( enzymeRule+post_first_idx-1, "[X]", post_list_size+2) == 0){
    free(post_cleavage_list);
    post_cleavage_list = NULL;
    post_list_size = 0;
    post_for_inclusion = false;
  }

  string precedingAminos;
  string followingAminos;


  // Allocate a hash table that maps from pairs of amino acids to Booleans.

}

/**
 * \returns An enzyme of the specified type.
 */
Enzyme::Enzyme(
  const char* enzymeName
) {
  return(init(enzymeName));
}

/**  
 * Decides whether the enzyme cleaves between the two specified amino acids.
 * Use 
 */
bool Enzyme::cleaves(
  char preceding,
  char following
) {
  return(precedingCleavage(preceding) && followingCleavage(following));
}

/**
 * Frees an allocated Enzyme object.
 */
Enzyme::~Scorer() {
  // FIXME
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
