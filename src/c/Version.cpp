#include "Version.h"

#include "search.h"

using namespace std;

Version::Version() {

}

Version::~Version() {
}


int Version::main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return 0;
}

string Version::getName() {
  return "version";
}

string Version::getDescription() {
  return 
    "Print the Crux version number to standard output, "
    "then exit";
}
