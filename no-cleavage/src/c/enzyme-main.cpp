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

  if (argc != 3) {
    fprintf(stderr, "enzyme <name>|<rule> <name>|<rule>\n");
    return(1);
  }
  set_verbosity_level(40);
  Enzyme myEnzyme1(argv[1]);
  Enzyme myEnzyme2(argv[2]);

  if (myEnzyme1.compare(&myEnzyme2)) {
    fprintf(stdout, "%s cleaves everywhere that %s cleaves.\n",
	    myEnzyme2.getName()->c_str(), myEnzyme1.getName()->c_str());
  }
  return(0);

}
