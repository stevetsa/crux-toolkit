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
#include "Peptide.h"


using namespace std;
using namespace Crux;
using namespace pwiz;
using namespace identdata;


/**
 * Initializes the object
 */
void MzIdentMLReader::init() {
  match_collection_ = NULL;
  use_pass_threshold_ = get_boolean_parameter("mzid-use-pass-threshold");
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

void MzIdentMLReader::addScores(
  const SpectrumIdentificationItem& item, 
  Match* match
) {
  vector<CVParam>::const_iterator iter = item.cvParams.begin();
  
  FLOAT_T fvalue;
  int ivalue;

  for (; iter != item.cvParams.end(); ++iter) {
    switch (iter->cvid) {
      case MS_Sequest_xcorr:
        from_string(fvalue, iter->value);
        match->setScore(XCORR, fvalue);
        match_collection_->setScoredType(XCORR, true);
        break;
      case MS_Sequest_PeptideSp:
        from_string(fvalue, iter->value);
        match->setScore(SP, fvalue);
        match_collection_->setScoredType(SP, true);
        break;
      case MS_Sequest_PeptideRankSp:
        from_string(ivalue, iter->value);
        match->setRank(SP, ivalue);
        break;
      case MS_Sequest_deltacn:
        from_string(fvalue, iter->value);
        match->setDeltaCn(fvalue);
        break;
      case MS_Sequest_matched_ions:
        from_string(ivalue, iter->value);
        match->setBYIonMatched(ivalue);
        break;
      case MS_Sequest_total_ions:
        from_string(ivalue, iter->value);
        match->setBYIonPossible(ivalue);
        break;
      default:
        carp(CARP_DEBUG, "Unknown score type, will be set in custom scores");
    }
    //go ahead and set all custom scores to the cvParam names.
    string name = cvTermInfo((*iter).cvid).name;
    from_string(fvalue, iter->value);
    match->setCustomScore(name, fvalue);
  }


  vector<UserParam>::const_iterator iter2 = item.userParams.begin();

  for (; iter2 != item.userParams.end(); ++iter2) {
    string name = iter2->name;

    bool success = from_string(fvalue, iter2->value);
    if (success) {
      match->setCustomScore(name, fvalue);
    }
  }




}

void printScores(const SpectrumIdentificationItem& item) {
  vector<CVParam>::const_iterator iter = item.cvParams.begin();
  for (; iter!=item.cvParams.end(); ++iter) {
    string name = cvTermInfo((*iter).cvid).name;

    if (iter->cvid == MS_Sequest_xcorr) { cerr<< "XCORR!"<<endl;}

    cerr << name <<" = " << iter->value << endl;
  }

}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* MzIdentMLReader::parse(
    const char* path, ///< path of the xml file
    Database* database, ///< target protein database
    Database* decoy_database ///< decoy protein database (can be null)
  ) {

  MzIdentMLReader* reader = new MzIdentMLReader(path);
  reader->setDatabase(database);
  reader->setDecoyDatabase(decoy_database);

  MatchCollection* collection = reader->parse();

  delete reader;

  return collection;

}
 

MatchCollection* MzIdentMLReader::parse() {

  match_collection_ = new MatchCollection();
  match_collection_ -> preparePostProcess();
  cerr << "MzIdentMLReader::opening file:"<<file_path_<<endl;
  pwiz_reader_ = new IdentDataFile(file_path_);

  parseDatabaseSequences();

  parsePSMs();

  return match_collection_;
}

void MzIdentMLReader::parseDatabaseSequences() {



}

void MzIdentMLReader::parsePSMs() {

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
        if (!use_pass_threshold_ || item.passThreshold) {
          int charge = item.chargeState;
          FLOAT_T obs_mz = item.experimentalMassToCharge;

          SpectrumZState zstate;
          zstate.setMZ(obs_mz, charge);
          vector<int> charge_vec;
          charge_vec.push_back(charge);

          //TODO crux requires scan numbers to be integer, where mzid can have
          //them be strings. Update crux to handle string type scan numbers.  

          Spectrum* spectrum = new Spectrum(0,0,obs_mz, charge_vec, "");

          FLOAT_T calc_mz = item.calculatedMassToCharge;
          FLOAT_T calc_mass = (calc_mz - MASS_PROTON ) * (FLOAT_T)charge;
          int rank = item.rank;


          PeptidePtr peptide_ptr = item.peptidePtr;
          string sequence = peptide_ptr->peptideSequence;
        
          vector<PeptideEvidencePtr>& peptide_evidences = item.peptideEvidencePtr;
          PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences.front();
          string protein_id = peptide_evidence_ptr->dbSequencePtr->accession;
          int start_idx = peptide_evidence_ptr->start;
          bool is_decoy = peptide_evidence_ptr->isDecoy;
          bool is_decoy_test;

          carp(CARP_DEBUG,"getting protein %s",protein_id.c_str());

          Protein* protein = MatchCollectionParser::getProtein(
            database_, decoy_database_, protein_id, is_decoy_test);

          if (is_decoy != is_decoy_test) {
            carp(CARP_WARNING, "mzid says %d, but database says %d", is_decoy, is_decoy_test);
          }
  
          start_idx = protein->findStart(sequence, "", "");
          int length = sequence.length();
        
          carp(CARP_DEBUG, "creating peptide %s %f %i",sequence.c_str(), calc_mass, start_idx);

          Crux::Peptide* peptide = 
            new Crux::Peptide(length, calc_mass, protein, start_idx);
        
          for (int pe_idx = 1; pe_idx < peptide_evidences.size();pe_idx++) {
            PeptideEvidencePtr peptide_evidence_ptr = peptide_evidences[pe_idx];
            int start = peptide_evidence_ptr->start;
            int end = peptide_evidence_ptr->end;
            protein_id = peptide_evidence_ptr->dbSequencePtr->accession; 
            bool decoy = peptide_evidence_ptr->isDecoy;
            carp(CARP_DEBUG, "id: %s start:%i end: %i decoy: %i", protein_id.c_str(),
             start, end, decoy);

            protein = MatchCollectionParser::getProtein(
              database_, decoy_database_, protein_id, is_decoy_test);
            start_idx = protein->findStart(sequence, "", "");
            PeptideSrc* src = new PeptideSrc((DIGEST_T)0, protein, start_idx);
            peptide->addPeptideSrc(src);
          }

          Match* match = new Match(peptide, spectrum, zstate, is_decoy);  
          match_collection_->addMatchToPostMatchCollection(match);

          match->setRank(XCORR, rank); // Is it safe to assume this?

          //cerr << "charge: "<<charge<<" obs mass:"<<obs_mass<<" calc mass:"<<calc_mass<<" sequence"<<sequence<<endl;
          addScores(item, match);
        }
      }
    }


    count++;
  }
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
