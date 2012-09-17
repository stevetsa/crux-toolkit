/*************************************************************************//**
 * \file MzIdentMLReader.cpp
 * \brief Object for parsing pepxml files
 ****************************************************************************/

#include "MzIdentMLReader.h"
#include "mass.h"
#include "expat.h"


#include "Protein.h"
#include "Peptide.h"


#include <cstdio>
#include <cstring>

#include <iostream>

#include "DelimitedFile.h"
#include "parameter.h"
#include "MatchCollectionParser.h"


using namespace std;
using namespace Crux;
using namespace pwiz;
using namespace identdata;


/**
 * Initializes the object
 */
void MzIdentMLReader::init() {
  current_match_collection_ = NULL;
}

/**
 * \returns an initialized object
 */
MzIdentMLReader::MzIdentMLReader() {
  init();
}

/**
 * \returns an object initialized with the file_path
 */
MzIdentMLReader::MzIdentMLReader(
  const string& file_path ///< the path of the pep.xml file
  ) {
  
  init();
  file_path_ = file_path;
}

/**
 * \returns an object initialized with the xml path, and the target,decoy databases
 */
MzIdentMLReader::MzIdentMLReader(
  const string& file_path, ///< the path of the pep.xml
  Database* database, ///< the protein database
  Database* decoy_database ///< the decoy protein database (can be null)
  ) {

  file_path_ = file_path;
  database_ = database;
  decoy_database_ = decoy_database;

}

void MzIdentMLReader::setDatabase(Database* database) {
  database_ = database;

}

void MzIdentMLReader::setDecoyDatabase(Database* decoy_database) {

  decoy_database_ = decoy_database;

}

/**
 * default destructor
 */
MzIdentMLReader::~MzIdentMLReader() {

}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */

void printScores(const SpectrumIdentificationItem& item) {
  vector<CVParam>::const_iterator iter = item.cvParams.begin();
  for (; iter!=item.cvParams.end(); ++iter) {
    string name = cvTermInfo((*iter).cvid).name;

    if (iter->cvid == MS_Sequest_xcorr) { cerr<< "XCORR!"<<endl;}

    cerr << name <<" = " << iter->value << endl;
  }

}


MatchCollection* MzIdentMLReader::parse() {

  MatchCollection* match_collection = new MatchCollection();
  //match_collection->preparePostProcess(database_->getNumProteins());
  cerr << "MzIdentMLReader::opening file:"<<file_path_<<endl;
  pwiz_reader_ = new IdentDataFile(file_path_);

  vector<SpectrumIdentificationListPtr>::const_iterator sil_iter;
  vector<SpectrumIdentificationListPtr>::const_iterator sil_end;
  vector<SpectrumIdentificationResultPtr>::const_iterator sir_iter;
  vector<SpectrumIdentificationItemPtr>::const_iterator sii_iter;

  sil_iter = pwiz_reader_->dataCollection.analysisData.spectrumIdentificationList.begin();
  sil_end = pwiz_reader_->dataCollection.analysisData.spectrumIdentificationList.end();

  int count = 0;

  for (; sil_iter != sil_end; ++sil_iter) {
    for (sir_iter = (**sil_iter).spectrumIdentificationResult.begin();
      sir_iter != (**sil_iter).spectrumIdentificationResult.end();
      ++sir_iter) {
      SpectrumIdentificationResult& result = **sir_iter;
      string idStr = result.spectrumID;
      string filename = result.spectraDataPtr->location;
      cerr << idStr <<" "<<filename << endl;

      for (sii_iter = result.spectrumIdentificationItem.begin();
        sii_iter != result.spectrumIdentificationItem.end();
        ++sii_iter) {
        SpectrumIdentificationItem& item = **sii_iter;



        int charge = item.chargeState;
        FLOAT_T obs_mz = item.experimentalMassToCharge;
  
        SpectrumZState zstate;
        zstate.setMZ(obs_mz, charge);
        vector<int> charge_vec;
        charge_vec.push_back(charge);

        Spectrum* spectrum = new Spectrum(0,0,obs_mz, charge_vec, "");

        


        FLOAT_T calc_mz = item.calculatedMassToCharge;
        int rank = item.rank;


        PeptidePtr peptide_ptr = item.peptidePtr;
        string sequence = peptide_ptr->peptideSequence;
        
        vector<PeptideEvidencePtr>& peptide_evidences = item.peptideEvidencePtr;
        PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences.front();
        string protein_id = peptide_evidence_ptr->dbSequencePtr->accession;
        int start_idx = peptide_evidence_ptr->start;
        bool is_decoy = peptide_evidence_ptr->isDecoy;
/*
        Protein* protein = database_->getProteinByIdString(protein_id.c_str());

        Peptide* peptide = new Peptide(length, obs_mass, protein, start_idx);

        Match* match = new Match(peptide, spectrum, zstate, is_decoy);

        match_collection->addMatchToPostMatchCollection(match);
*/
        
        for (int pe_idx = 1; pe_idx < peptide_evidences.size();pe_idx++) {
          PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences[pe_idx];
          int start = peptide_evidence_ptr->start;
          int end = peptide_evidence_ptr->end;
          string id = peptide_evidence_ptr->dbSequencePtr->accession; 
          bool decoy = peptide_evidence_ptr->isDecoy;
          cerr <<"id:"<<id<<" start:"<<start<<" end:"<<end<<" decoy:"<<decoy<<endl;
        }
        
        //cerr << "charge: "<<charge<<" obs mass:"<<obs_mass<<" calc mass:"<<calc_mass<<" sequence"<<sequence<<endl;

        printScores(item);
      }


    }


    cerr << "count:"<<count << endl;
    count++;
  }
  
  

  cerr << "done!"<<endl;
  exit(-1);
  



  delete pwiz_reader_;

  return match_collection;
}

//TODO - remove this code after some time of debugging.
#ifdef MAIN
int main(int argc, char** argv) {

  initialize_parameters();

  char* file_path = argv[1];
//  char* database_path = argv[2];

  Database* database = NULL;
  Database* decoy_database = NULL;
 
//  MatchCollectionParser::loadDatabase(database_path, database, decoy_database);


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
    cout <<" rank:"<<match->getRank(XCORR);
    cout <<" sequence:"<<match->getPeptide()->getSequence();
    cout <<" protein:"<< match->getPeptide()->getProteinIdsLocations()<<endl;

  }



  return 0;
}
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
