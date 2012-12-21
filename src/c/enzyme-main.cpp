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
#include "carp.h"
int main(int argc, char** argv){

  if (argc != 2) {
    fprintf(stderr, "enzyme <name>|<rule>\n");
    return(1);
  }
  set_verbosity_level(40);
  Enzyme myEnzyme(argv[1]);
  return(0);

}
