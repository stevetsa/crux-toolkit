#include "SearchForXLinksMods.h"

#include "xlink_search.h"

using namespace std;

SearchForXLinksMods::SearchForXLinksMods() {

}

SearchForXLinksMods::~SearchForXLinksMods() {
}


int SearchForXLinksMods::main(int argc, char** argv) {
  return xlink_search_main(argc, argv);
}

string SearchForXLinksMods::getName() {
  return "search-for-xlinks-mods";
}

string SearchForXLinksMods::getDescription() {
  return 
    "Search a collection of spectra against a sequence "
    "database returning a collection of matches "
    "corresponding to linear and cross-lnked peptides "
    "scored by XCorr. This code also searches with "
    "post-translational modifications";
}
/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T SearchForXLinksMods::getCommand() {
  return XLINK_SEARCH_MODS_COMMAND;
}

bool SearchForXLinksMods::needsOutputDirectory() {
  return true;
}

