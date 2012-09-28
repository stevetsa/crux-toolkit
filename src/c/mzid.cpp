#include "objects.h"
#include "MzIdentMLReader.h"
#include "parameter.h"
#include "Database.h"
#include "Match.h"
#include "MatchCollectionParser.h"
#include "MatchCollection.h"
#include "MatchIterator.h"
#include "MzIdentMLReader.h"

#include <iostream>

using namespace Crux;


int main(int argc, char** argv) {

  initialize_parameters();
  set_verbosity_level(CARP_INFO);
  char* file_path = argv[1];
  const char* database_path = get_string_parameter_pointer("protein-database");

  Database* database = NULL;
  Database* decoy_database = NULL;
 
  MatchCollectionParser::loadDatabase(database_path, database, decoy_database);


  cerr << "creating reader"<<endl;
  MzIdentMLReader* reader = new MzIdentMLReader(file_path);
  reader->setDatabase(database);
  reader->setDecoyDatabase(decoy_database);
  cerr << "calling parse"<<endl;
  MatchCollection* match_collection = reader->parse();


  cerr << "there are "<<match_collection->getMatchTotal()<<" matches read"<<endl;

  MatchIterator* match_iterator = new MatchIterator(match_collection, XCORR, true);

  while(match_iterator->hasNext()) {
    Match* match = match_iterator->next();

    cout << "xcorr:"<<match->getScore(XCORR);
    if (match_collection->getScoredType(SP)) {
      cout <<" sp:"<<match->getScore(SP);
    }
    cout <<" rank:"<<match->getRank(XCORR);
    cout <<" sequence:"<<match->getPeptide()->getSequence();
//    cout <<" protein:"<< match->getPeptide()->getProteinIdsLocations()<<endl;
    cout << endl;
  }



  return 0;
}
