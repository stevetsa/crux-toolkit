/**
 * \file enzyme-main.cpp
 * AUTHOR: William Stafford Noble
 * CREATE DATE: December 20, 2012
 * \brief Tiny little program for testing the Enzyme object.
 *
 * Usage is "enzyme <name>|<rule>"
 * Output is only to stderr.
 **/

#include "Enzyme.h"
int main(int argc, char** argv){

  Enzyme myEnzyme(argv[1]);
  return(0)

}
