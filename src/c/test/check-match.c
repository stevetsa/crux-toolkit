#include <stdlib.h>
#include "check-match.h"
#include "../match.h"

// declare things to set up
MATCH_T *m1, *mdecoy;
PROTEIN_T *prot;
PEPTIDE_T *pep;
SPECTRUM_T *spec;
char* protseq = "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFVDQEFAQIKHVLHGISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLLAVGHYLESTVDIAESTRRIAASRIPADHMVLMAGFTAGNEKGELVVLGRNGSDYSAAVLAACLRADCCEIWTDVDGVYTCDPRQVPDARLLKSMSYQEAMELSYFGAKVLHPRTITPIAQFQIPCLIKNTGNPQAPGTLIGASRDEDELPVKGISNLNNMAMFSVSGPGMKGMVGMAARVFAAMSRARISVVLITQSSSEYSISFCVPQSDCVRAERAMQEEFYLELKEGLLEPLAVTERLAIISVVGDGMRTLRGISAKFFAALARANINIVAIAQGSSERSISVVVNNDDATTGVRVTHQMLFNTDQVIEVFVIGVGGVGGALLEQLKRQQSW";

void match_setup(){
  // create inputs
  prot = new_protein( "Protein1", protseq, strlen(protseq), 
                          NULL, 0, 0, NULL);//description, offset, idx, dbase
  pep = new_peptide( 10, 1087.20, prot, 20, TRYPTIC);//VADILESNAR

  // create match
  m1 = new_match();
  set_match_peptide(m1, pep);
  set_match_spectrum(m1, NULL);
  set_match_charge(m1, 2);
  set_match_null_peptide(m1, FALSE);

  // create match to null (decoy) peptide
  mdecoy = new_match();
  set_match_peptide(mdecoy, pep);
  set_match_spectrum(mdecoy, NULL);
  set_match_charge(mdecoy, 2);
  set_match_null_peptide(mdecoy, TRUE);
}

void match_teardown(){
  free_match(m1);
  free_match(mdecoy);
  free_protein(prot);
  free_peptide(pep);

}

START_TEST(test_create){
  fail_unless( m1 != NULL, "Failed to allocate a match.");

  // test getters
  char* seq = get_match_sequence(m1);
  fail_unless( strcmp(seq, "VADILESNAR") == 0,
               "Match returns %s as seq instead of %s", seq, "VADILESNAR");
  MODIFIED_AA_T* mod_seq = get_match_mod_sequence(m1);
  char* mod_str = modified_aa_string_to_string(mod_seq);
  fail_unless( strcmp(mod_str, seq) == 0,
               "MOD_AA string should be %s but is %s", seq, mod_str);

  // test getting seqs for null (decoy) peptide
  free(seq);
  seq = get_match_sequence(mdecoy);
  fail_unless( strcmp(seq, "VASDLINEAR") == 0,
               "Match decoy seq should be %s but is %s", "VASDLINEAR", seq );
  free(mod_seq);
  mod_seq = get_match_mod_sequence(mdecoy);
  free(mod_str);
  mod_str = modified_aa_string_to_string(mod_seq);
  fail_unless( strcmp(mod_str, seq) != 0,
               "For peptide with seq %s, shuffled MOD_AA string should " \
               "be different but is %s", seq, mod_str);
}
END_TEST

/*
START_TEST(test_set){
}
END_TEST

START_TEST(test_create){
}
END_TEST

START_TEST(test_create){
}
END_TEST
*/

Suite* match_suite(){
  Suite* s = suite_create("Match");
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  //  tcase_add_test(tc_core, ????);

  tcase_add_checked_fixture(tc_core, match_setup, match_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  /*
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, ????);
  tcase_add_checked_fixture(tc_limits, mod_setup, mod_teardown);
  suite_ad_tcase(s, tc_limits);
   */

  return s;
}

/*
#define scan_num 16
#define ms2_file "test.ms2"
#define parameter_file "test_parameter_file"

//THE parameter system will not work if set CK_FORK=no
START_TEST (test_create){
  SPECTRUM_T* spectrum = NULL;
  SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_T* match = NULL;
    
  // comment this parameter section out, when using CK_FORK=no, valgrind
  // Parameters have been already confirmed in check_scorer.
  
  //set verbbosity level
  int  verbosity = CARP_INFO;
  set_verbosity_level(verbosity);

  //parse paramter file
  //  parse_update_parameters(parameter_file);
  initialize_parameters();

  //add fasta file parameter_file fasta-file
  //add_parameter("fasta-file", "fasta_file");
  
  //parameters has been confirmed
  //parameters_confirmed();
  // ****************************************** end comment out

  //read ms2 file
  collection = new_spectrum_collection(ms2_file);
  spectrum = allocate_spectrum();

  DATABASE_T* database = new_database("fasta-file", FALSE);

  
  //search for spectrum with correct scan number
  fail_unless(get_spectrum_collection_spectrum(collection, scan_num, spectrum), "failed to find scan_num in ms3 file");
  
  //get match collection with perliminary score of SP, and main score of XCORR
  match_collection = new_match_collection_from_spectrum(spectrum, 1, 500, SP, XCORR, 0, FALSE, NULL, database);
  
  fail_unless(get_match_collection_scored_type(match_collection, SP), "failed to set match_collection scored type, SP");
  fail_unless(get_match_collection_scored_type(match_collection, XCORR), "failed to set match_collection scored type, SP");

  //LOGP_EXP_SP should not be scored yet
  fail_unless(!get_match_collection_scored_type(match_collection, LOGP_EXP_SP), "failed to set match_collection scored type, xcorr");
  fail_unless(!get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 
  
  //create match iterator
  match_iterator = new_match_iterator(match_collection, SP, TRUE);
  
  //match_collection should be locked now..
  fail_unless(get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 
  
  //iterate over all matches
  while(match_iterator_has_next(match_iterator)){
    match = match_iterator_next(match_iterator);
    print_match(match, stdout, TRUE, SP);
  }

  //free match iterator
  free_match_iterator(match_iterator);
  
  //should be unlocked
  fail_unless(!get_match_collection_iterator_lock(match_collection), "match_collection lock is not set correctly"); 

  free_match_collection(match_collection);
  free_spectrum_collection(collection);
  free_spectrum(spectrum);
  free_parameters();
}
END_TEST


Suite *match_suite(void){
  Suite *s = suite_create("match");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
*/
