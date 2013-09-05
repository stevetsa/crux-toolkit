/**
 * \file CometApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "Common.h"
#include "CometSearchManager.h"
#include "ModifiedPeptidesIterator.h"
#include "CarpStreamBuf.h"
#include "CometApplication.h"
#include "DelimitedFileWriter.h"
#include "DelimitedFile.h"

using namespace std;

/**
 * \returns a blank CometApplication object
 */
CometApplication::CometApplication() {

}

/**
 * Destructor
 */
CometApplication::~CometApplication() {
}

/**
 * main method for CometApplication
 */
int CometApplication::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };

  
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"input spectra","database_name"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  /* Re-route stderr to log file */
  CarpStreamBuf buffer;
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  /* set Parmeters */
  vector<InputFileInfo*> pv_input_files;
  CometSearchManager comet_search_mgr;
  setCometParameters(pv_input_files, comet_search_mgr);
  comet_search_mgr.AddInputFiles(pv_input_files);
  
  /* Run search */
  comet_search_mgr.DoSearch();

  /* Recover stderr */
  std::cerr.rdbuf(old);

  return 0;
}

void calcVarMods(const char* var_str, VarMods& varmods) {
  
  string temp = var_str;
  vector<string> tokens;
  DelimitedFile::tokenize(temp, tokens, ' ');
  
  from_string<double>(varmods.dVarModMass, tokens[0]);
  strcpy(varmods.szVarModChar, tokens[1].c_str());
  from_string<int>(varmods.bBinaryMod, tokens[2]);
  from_string<int>(varmods.iMaxNumVarModAAPerMod, tokens[3]);
}

void getDoubleRange(const char* str, DoubleRange& doubleRangeParam) {
  
  string temp = str;
  vector<string> tokens;
  DelimitedFile::tokenize(temp, tokens, ' ');
  
  from_string<double>(doubleRangeParam.dStart, tokens[0]);
  from_string<double>(doubleRangeParam.dEnd, tokens[1]);
}

void getIntRange(const char* str, IntRange& intRangeParam) {
  
  string temp = str;
  vector<string> tokens;
  DelimitedFile::tokenize(temp, tokens, ' ');
  
  from_string<int>(intRangeParam.iStart, tokens[0]);
  from_string<int>(intRangeParam.iEnd, tokens[1]);
}

void getEnzymeInfo(int search_enzyme_number, int sample_enzyme_number, EnzymeInfo& enzymeInformation) {
  
  switch(search_enzyme_number) {
      case 0:
        strcpy(enzymeInformation.szSearchEnzymeName, "No_enzyme");
        enzymeInformation.iSearchEnzymeOffSet = 0;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "-");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "-");
        break;
      case 1:
        strcpy(enzymeInformation.szSearchEnzymeName, "Trypsin");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "KR");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "P");
        break;
      case 2:
        strcpy(enzymeInformation.szSearchEnzymeName, "Trypsin/P");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "KR");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "-");
        break;
      case 3:
        strcpy(enzymeInformation.szSearchEnzymeName, "Lys_C");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "K");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "P");
        break;
      case 4:
        strcpy(enzymeInformation.szSearchEnzymeName, "Lys_N");
        enzymeInformation.iSearchEnzymeOffSet = 0;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "K");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "-");
        break;
      case 5:
        strcpy(enzymeInformation.szSearchEnzymeName, "Arg_C");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "R");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "P");
        break;
      case 6:
        strcpy(enzymeInformation.szSearchEnzymeName, "Asp_N");
        enzymeInformation.iSearchEnzymeOffSet = 0;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "D");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "-");
        break;
      case 7:
        strcpy(enzymeInformation.szSearchEnzymeName, "CNBr");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "M");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "-");
        break;
      case 8:
        strcpy(enzymeInformation.szSearchEnzymeName, "Glu_C");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "DE");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "P");
        break;
      case 9:
        strcpy(enzymeInformation.szSearchEnzymeName, "PepsinA");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "FL");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "P");
        break;
      case 10:
        strcpy(enzymeInformation.szSearchEnzymeName, "Chymotrypsin");
        enzymeInformation.iSearchEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSearchEnzymeBreakAA, "FWYL");
        strcpy(enzymeInformation.szSearchEnzymeNoBreakAA, "P");
        break;
      default:
        carp(CARP_FATAL, "Unknown enzyme number!");
  }

    
  switch(sample_enzyme_number) {
      case 0:
        strcpy(enzymeInformation.szSampleEnzymeName, "No_enzyme");
        enzymeInformation.iSampleEnzymeOffSet = 0;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "-");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "-");
        break;
      case 1:
        strcpy(enzymeInformation.szSampleEnzymeName, "Trypsin");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "KR");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "P");
        break;
      case 2:
        strcpy(enzymeInformation.szSampleEnzymeName, "Trypsin/P");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "KR");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "-");
        break;
      case 3:
        strcpy(enzymeInformation.szSampleEnzymeName, "Lys_C");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "K");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "P");
        break;
      case 4:
        strcpy(enzymeInformation.szSampleEnzymeName, "Lys_N");
        enzymeInformation.iSampleEnzymeOffSet = 0;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "K");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "-");
        break;
      case 5:
        strcpy(enzymeInformation.szSampleEnzymeName, "Arg_C");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "R");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "P");
        break;
      case 6:
        strcpy(enzymeInformation.szSampleEnzymeName, "Asp_N");
        enzymeInformation.iSampleEnzymeOffSet = 0;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "D");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "-");
        break;
      case 7:
        strcpy(enzymeInformation.szSampleEnzymeName, "CNBr");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "M");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "-");
        break;
      case 8:
        strcpy(enzymeInformation.szSampleEnzymeName, "Glu_C");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "DE");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "P");
        break;
      case 9:
        strcpy(enzymeInformation.szSampleEnzymeName, "PepsinA");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "FL");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "P");
        break;
      case 10:
        strcpy(enzymeInformation.szSampleEnzymeName, "Chymotrypsin");
        enzymeInformation.iSampleEnzymeOffSet = 1;
        strcpy(enzymeInformation.szSampleEnzymeBreakAA, "FWYL");
        strcpy(enzymeInformation.szSampleEnzymeNoBreakAA, "P");
        break;
      default:
        carp(CARP_FATAL, "Unknown enzyme number!");
  }


}

void CometApplication::setCometParameters(
  vector<InputFileInfo*> &pvInputFiles,
  CometSearchManager& searchMgr
  ) {
  
  VarMods varModsParam;
   IntRange intRangeParam;
   DoubleRange doubleRangeParam;

  
  InputFileInfo *pInputFile = new InputFileInfo();
  pInputFile->iAnalysisType = 0;
  strcpy(pInputFile->szFileName, get_string_parameter_pointer("input spectra"));
  
  //VALIDATE
  FILE *fp;
  if ((fp = fopen(pInputFile->szFileName, "r")) == NULL) {
      carp(CARP_FATAL, "Spectra File Not Found:%s", get_string_parameter_pointer("input spectra"));
  }
  fclose(fp);
  
  string scan_range_str = get_string_parameter_pointer("scan_range");
  if (scan_range_str == "0 0") {
    pInputFile->iAnalysisType = AnalysisType_EntireFile;
  } else {
    pInputFile->iAnalysisType = AnalysisType_SpecificScanRange;
    vector<string> tokens;
    DelimitedFile::tokenize(scan_range_str, tokens, ' ');
    from_string<int>(pInputFile->iFirstScan, tokens[0]);
    from_string<int>(pInputFile->iLastScan, tokens[1]);
  }

  pvInputFiles.push_back(pInputFile);
  //TODO - Comet allows multiple spectra to be searched, add this to crux.
  string basename = make_file_path(getName());
  searchMgr.SetOutputFileBaseName(basename.c_str());
  
  searchMgr.SetParam("database_name", get_string_parameter_pointer("database_name"), get_string_parameter_pointer("database_name"));
  searchMgr.SetParam("decoy_prefix", get_string_parameter_pointer("decoy_prefix"), get_string_parameter_pointer("decoy_prefix"));
  searchMgr.SetParam("nucleotide_reading_frame", get_string_parameter_pointer("nucleotide_reading_frame"), get_int_parameter("nucleotide_reading_frame"));
  searchMgr.SetParam("mass_type_parent", get_string_parameter_pointer("mass_type_parent"), get_int_parameter("mass_type_parent"));
  searchMgr.SetParam("mass_type_fragment", get_string_parameter_pointer("mass_type_fragment"), get_int_parameter("mass_type_fragment"));
  searchMgr.SetParam("show_fragment_ions", get_string_parameter_pointer("show_fragment_ions"), get_int_parameter("show_fragment_ions"));
  searchMgr.SetParam("num_threads", get_string_parameter_pointer("num_threads"), get_int_parameter("num_threads"));
  searchMgr.SetParam("clip_nterm_methionine", get_string_parameter_pointer("clip_nterm_methionine"), get_int_parameter("clip_nterm_methionine"));
  searchMgr.SetParam("theoretical_fragment_ions", get_string_parameter_pointer("theoretical_fragment_ions"), get_int_parameter("theoretical_fragment_ions"));
  searchMgr.SetParam("use_A_ions", get_string_parameter_pointer("use_A_ions"), get_int_parameter("use_A_ions"));
  searchMgr.SetParam("use_B_ions", get_string_parameter_pointer("use_B_ions"), get_int_parameter("use_B_ions"));
  searchMgr.SetParam("use_C_ions", get_string_parameter_pointer("use_C_ions"), get_int_parameter("use_C_ions"));
  searchMgr.SetParam("use_X_ions", get_string_parameter_pointer("use_X_ions"), get_int_parameter("use_X_ions"));
  searchMgr.SetParam("use_Y_ions", get_string_parameter_pointer("use_Y_ions"), get_int_parameter("use_Z_ions"));
  searchMgr.SetParam("use_Z_ions", get_string_parameter_pointer("use_Z_ions"), get_int_parameter("use_Z_ions"));
  searchMgr.SetParam("use_NL_ions", get_string_parameter_pointer("use_NL_ions"), get_int_parameter("use_NL_ions"));
  searchMgr.SetParam("use_sparse_matrix", get_string_parameter_pointer("use_sparse_matrix"), get_int_parameter("use_sparse_matrix"));
  
  calcVarMods(get_string_parameter_pointer("variable_mod1"), varModsParam);
  searchMgr.SetParam("variable_mod1", get_string_parameter_pointer("variable_mod1"), varModsParam );
  calcVarMods(get_string_parameter_pointer("variable_mod2"), varModsParam);
  searchMgr.SetParam("variable_mod2", get_string_parameter_pointer("variable_mod2"), varModsParam );
  calcVarMods(get_string_parameter_pointer("variable_mod3"), varModsParam);
  searchMgr.SetParam("variable_mod3", get_string_parameter_pointer("variable_mod3"), varModsParam );
  calcVarMods(get_string_parameter_pointer("variable_mod4"), varModsParam);
  searchMgr.SetParam("variable_mod4", get_string_parameter_pointer("variable_mod4"), varModsParam );
  calcVarMods(get_string_parameter_pointer("variable_mod5"), varModsParam);
  searchMgr.SetParam("variable_mod5", get_string_parameter_pointer("variable_mod5"), varModsParam );
  calcVarMods(get_string_parameter_pointer("variable_mod6"), varModsParam);
  searchMgr.SetParam("variable_mod6", get_string_parameter_pointer("variable_mod6"), varModsParam );
  
  searchMgr.SetParam("max_variable_mods_in_peptide", get_string_parameter_pointer("max_variable_mods_in_peptide"), get_int_parameter("max_variable_mods_in_peptide"));
  searchMgr.SetParam("fragment_bin_tol", get_string_parameter_pointer("fragment_bin_tol"), get_double_parameter("fragment_bin_tol"));
  searchMgr.SetParam("fragment_bin_offset", get_string_parameter_pointer("fragment_bin_offset"), get_double_parameter("fragment_bin_offset"));
  searchMgr.SetParam("peptide_mass_tolerance", get_string_parameter_pointer("peptide_mass_tolerance"), get_double_parameter("peptide_mass_tolerance"));
  searchMgr.SetParam("precursor_tolerance_type", get_string_parameter_pointer("precursor_tolerance_type"), get_int_parameter("precursor_tolerance_type"));
  searchMgr.SetParam("peptide_mass_units", get_string_parameter_pointer("peptide_mass_units"), get_int_parameter("peptide_mass_units"));
  searchMgr.SetParam("isotope_error", get_string_parameter_pointer("isotope_error"), get_int_parameter("isotope_error"));
  searchMgr.SetParam("num_output_lines", get_string_parameter_pointer("num_output_lines"), get_int_parameter("num_output_lines"));
  searchMgr.SetParam("num_results", get_string_parameter_pointer("num_results"), get_int_parameter("num_results"));
  searchMgr.SetParam("remove_precursor_peak", get_string_parameter_pointer("remove_precursor_peak"), get_int_parameter("remove_precursor_peak"));
  searchMgr.SetParam("remove_precursor_tolerance", get_string_parameter_pointer("remove_precursor_tolerance"), get_double_parameter("remove_precursor_tolerance"));
  
  getDoubleRange(get_string_parameter_pointer("clear_mz_range"), doubleRangeParam );
  searchMgr.SetParam("clear_mz_range", get_string_parameter_pointer("clear_mz_range"), doubleRangeParam );

  searchMgr.SetParam("print_expect_score", get_string_parameter_pointer("print_expect_score"), get_int_parameter("print_expect_score"));
  searchMgr.SetParam("output_sqtstream", "0", 0);
  
  searchMgr.SetParam("output_sqtfile", get_string_parameter_pointer("output_sqtfile"), get_int_parameter("output_sqtfile"));
  searchMgr.SetParam("output_txtfile", get_string_parameter_pointer("output_txtfile"), get_int_parameter("output_txtfile"));
  searchMgr.SetParam("output_pepxmlfile", get_string_parameter_pointer("output_pepxmlfile"), get_int_parameter("output_pepxmlfile"));
  //searchMgr.SetParam("output_pinxmlfile", get_string_parameter_pointer("output_pinxmlfile"), get_int_parameter("output_pinxmlfile"));
  searchMgr.SetParam("output_outfiles", "0", "0");
  searchMgr.SetParam("skip_researching", get_string_parameter_pointer("skip_researching"), get_int_parameter("skip_researching"));
  searchMgr.SetParam("variable_C_terminus", get_string_parameter_pointer("variable_C_terminus"), get_double_parameter("variable_C_terminus"));
  searchMgr.SetParam("variable_N_terminus", get_string_parameter_pointer("variable_N_terminus"), get_double_parameter("variable_N_terminus"));
  searchMgr.SetParam("variable_C_terminus_distance", get_string_parameter_pointer("variable_C_terminus_distance"), get_int_parameter("variable_C_terminus_distance"));
  searchMgr.SetParam("variable_N_terminus_distance", get_string_parameter_pointer("variable_N_terminus_distance"), get_int_parameter("variable_N_terminus_distance"));
  searchMgr.SetParam("add_Cterm_peptide", get_string_parameter_pointer("add_Cterm_peptide"), get_double_parameter("add_Cterm_peptide"));
  searchMgr.SetParam("add_Nterm_peptide", get_string_parameter_pointer("add_Nterm_peptide"), get_double_parameter("add_Nterm_peptide"));
  searchMgr.SetParam("add_Cterm_protein", get_string_parameter_pointer("add_Cterm_protein"), get_double_parameter("add_Cterm_protein"));
  searchMgr.SetParam("add_Nterm_protein", get_string_parameter_pointer("add_Nterm_protein"), get_double_parameter("add_Nterm_protein"));
  searchMgr.SetParam("add_G_glycine", get_string_parameter_pointer("add_G_glycine"), get_double_parameter("add_G_glycine"));
  searchMgr.SetParam("add_A_alanine", get_string_parameter_pointer("add_A_alanine"), get_double_parameter("add_A_alanine"));
  searchMgr.SetParam("add_S_serine", get_string_parameter_pointer("add_S_serine"), get_double_parameter("add_S_serine"));
  searchMgr.SetParam("add_P_proline", get_string_parameter_pointer("add_P_proline"), get_double_parameter("add_P_proline"));
  searchMgr.SetParam("add_V_valine", get_string_parameter_pointer("add_V_valine"), get_double_parameter("add_V_valine"));
  searchMgr.SetParam("add_T_threonine", get_string_parameter_pointer("add_T_threonine"), get_double_parameter("add_T_threonine"));
  searchMgr.SetParam("add_C_cysteine", get_string_parameter_pointer("add_C_cysteine"), get_double_parameter("add_C_cysteine"));
  searchMgr.SetParam("add_L_leucine", get_string_parameter_pointer("add_L_leucine"), get_double_parameter("add_L_leucine"));
  searchMgr.SetParam("add_I_isoleucine", get_string_parameter_pointer("add_I_isoleucine"), get_double_parameter("add_I_isoleucine"));
  searchMgr.SetParam("add_N_asparagine", get_string_parameter_pointer("add_N_asparagine"), get_double_parameter("add_N_asparagine"));
  searchMgr.SetParam("add_O_ornithine", get_string_parameter_pointer("add_O_ornithine"), get_double_parameter("add_O_ornithine"));
  searchMgr.SetParam("add_D_aspartic_acid", get_string_parameter_pointer("add_D_aspartic_acid"), get_double_parameter("add_D_aspartic_acid"));
  searchMgr.SetParam("add_Q_glutamine", get_string_parameter_pointer("add_Q_glutamine"), get_double_parameter("add_Q_glutamine"));
  searchMgr.SetParam("add_K_lysine", get_string_parameter_pointer("add_K_lysine"), get_double_parameter("add_K_lysine"));
  searchMgr.SetParam("add_E_glutamic_acid", get_string_parameter_pointer("add_E_glutamic_acid"), get_double_parameter("add_E_glutamic_acid"));
  searchMgr.SetParam("add_M_methionine", get_string_parameter_pointer("add_M_methionine"), get_double_parameter("add_M_methionine"));
  searchMgr.SetParam("add_H_histidine", get_string_parameter_pointer("add_H_histidine"), get_double_parameter("add_H_histidine"));
  searchMgr.SetParam("add_F_phenylalanine", get_string_parameter_pointer("add_F_phenylalanine"), get_double_parameter("add_F_phenylalanine"));
  searchMgr.SetParam("add_R_arginine", get_string_parameter_pointer("add_R_arginine"), get_double_parameter("add_R_arginine"));
  searchMgr.SetParam("add_Y_tyrosine", get_string_parameter_pointer("add_Y_tyrosine"), get_double_parameter("add_Y_tyrosine"));
  searchMgr.SetParam("add_W_tryptophan", get_string_parameter_pointer("add_W_tryptophan"), get_double_parameter("add_W_tryptophan"));
  searchMgr.SetParam("add_B_user_amino_acid", get_string_parameter_pointer("add_B_user_amino_acid"), get_double_parameter("add_B_user_amino_acid"));
  searchMgr.SetParam("add_J_user_amino_acid", get_string_parameter_pointer("add_J_user_amino_acid"), get_double_parameter("add_J_user_amino_acid"));
  searchMgr.SetParam("add_U_user_amino_acid", get_string_parameter_pointer("add_U_user_amino_acid"), get_double_parameter("add_U_user_amino_acid"));
  searchMgr.SetParam("add_X_user_amino_acid", get_string_parameter_pointer("add_X_user_amino_acid"), get_double_parameter("add_X_user_amino_acid"));
  searchMgr.SetParam("add_Z_user_amino_acid", get_string_parameter_pointer("add_Z_user_amino_acid"), get_double_parameter("add_Z_user_amino_acid"));
  searchMgr.SetParam("search_enzyme_number", get_string_parameter_pointer("search_enzyme_number"), get_int_parameter("search_enzyme_number"));
  searchMgr.SetParam("sample_enzyme_number", get_string_parameter_pointer("sample_enzyme_number"), get_int_parameter("sample_enzyme_number"));
  searchMgr.SetParam("num_enzyme_termini", get_string_parameter_pointer("num_enzyme_termini"), get_int_parameter("num_enzyme_termini"));
  searchMgr.SetParam("allowed_missed_cleavage", get_string_parameter_pointer("allowed_missed_cleavage"), get_int_parameter("allowed_missed_cleavage"));

  getIntRange(get_string_parameter_pointer("scan_range"), intRangeParam );
  searchMgr.SetParam("scan_range", get_string_parameter_pointer("scan_range"), intRangeParam );

  searchMgr.SetParam("spectrum_batch_size", get_string_parameter_pointer("spectrum_batch_size"), get_int_parameter("spectrum_batch_size"));
  searchMgr.SetParam("minimum_peaks", get_string_parameter_pointer("minimum_peaks"), get_int_parameter("minimum_peaks"));

  getIntRange(get_string_parameter_pointer("precursor_charge"), intRangeParam);
  searchMgr.SetParam("precursor_charge", get_string_parameter_pointer("precursor_charge"), intRangeParam);
  
  searchMgr.SetParam("max_fragment_charge", get_string_parameter_pointer("max_fragment_charge"), get_int_parameter("max_fragment_charge"));
  searchMgr.SetParam("max_precursor_charge", get_string_parameter_pointer("max_precursor_charge"), get_int_parameter("max_precursor_charge"));

  getDoubleRange(get_string_parameter_pointer("digest_mass_range"), doubleRangeParam);
  searchMgr.SetParam("digest_mass_range", get_string_parameter_pointer("digest_mass_range"), doubleRangeParam);
  
  searchMgr.SetParam("ms_level", get_string_parameter_pointer("ms_level"), get_int_parameter("ms_level"));
  searchMgr.SetParam("activation_method", get_string_parameter_pointer("activation_method"), get_string_parameter_pointer("activation_method"));
  searchMgr.SetParam("minimum_intensity", get_string_parameter_pointer("minimum_intensity"), get_double_parameter("minimum_intensity"));
  searchMgr.SetParam("decoy_search", get_string_parameter_pointer("decoy_search"), get_int_parameter("decoy_search"));
  
  EnzymeInfo enzymeInformation;
  
  getEnzymeInfo(get_int_parameter("search_enzyme_number"), get_int_parameter("sample_enzyme_number"), enzymeInformation);
  searchMgr.SetParam("[COMET_ENZYME_INFO]", "TODO", enzymeInformation);
  
}


/**
 * \write parameters 
 */
void CometApplication:: writeParams(ofstream &fout){
  string space="\t\t\t\t";
  string header_file= "# comet_version 2013.01 rev. 0"
 "\n# Comet MS/MS search engine parameters file."
  "\n# Everything following the \'#\' symbol is treated as a comment.\n";
  fout<<header_file<<endl;

  fout << "database_name="<<get_string_parameter("database_name")<<endl;
  fout << "decoy_search="<<get_int_parameter("decoy_search")<<endl;
  fout << "num_threads="<<get_int_parameter("num_threads")<<endl;
  fout << "peptide_mass_tolerance="<<get_double_parameter("peptide_mass_tolerance")<<endl;
  fout << "peptide_mass_units="<<get_int_parameter("peptide_mass_units")<<endl;
  fout << "mass_type_parent="<<get_int_parameter("mass_type_parent")<<endl;
  fout << "mass_type_fragment="<<get_int_parameter("mass_type_fragment")<<endl;
  fout << "precursor_tolerance_type="<<get_int_parameter("precursor_tolerance_type")<<endl;
  fout << "isotope_error="<<get_int_parameter("isotope_error")<<endl;
  fout << "search_enzyme_number="<<get_int_parameter("search_enzyme_number")<<endl;
  fout << "num_enzyme_termini="<<get_int_parameter("num_enzyme_termini")<<endl;
  fout << "allowed_missed_cleavage="<<get_int_parameter("allowed_missed_cleavage")<<endl;
  fout << "variable_mod1="<<get_string_parameter("variable_mod1")<<endl;
  fout << "variable_mod2="<<get_string_parameter("variable_mod2")<<endl;
  fout << "variable_mod3="<<get_string_parameter("variable_mod3")<<endl;
  fout << "variable_mod4="<<get_string_parameter("variable_mod4")<<endl;
  fout << "variable_mod5="<<get_string_parameter("variable_mod5")<<endl;
  fout << "variable_mod6="<<get_string_parameter("variable_mod6")<<endl;
  fout << "max_variable_mods_in_peptide="<<get_int_parameter("max_variable_mods_in_peptide")<<endl;
  fout << "fragment_bin_tol="<<get_double_parameter("fragment_bin_tol")<<endl;
  fout << "fragment_bin_offset="<<get_double_parameter("fragment_bin_offset")<<endl;
  fout << "theoretical_fragment_ions="<<get_int_parameter("theoretical_fragment_ions")<<endl;
  fout << "use_A_ions="<<get_int_parameter("use_A_ions")<<endl;
  fout << "use_B_ions="<<get_int_parameter("use_B_ions")<<endl;
  fout << "use_C_ions="<<get_int_parameter("use_C_ions")<<endl;
  fout << "use_X_ions="<<get_int_parameter("use_X_ions")<<endl;
  fout << "use_Y_ions="<<get_int_parameter("use_Y_ions")<<endl;
  fout << "use_Z_ions="<<get_int_parameter("use_Z_ions")<<endl;
  fout << "use_NL_ions="<<get_int_parameter("use_NL_ions")<<endl;
  fout << "use_sparse_matrix="<<get_int_parameter("use_sparse_matrix")<<endl;
  fout << "output_sqtfile="<<get_int_parameter("output_sqtfile")<<endl;
  //  fout << "output_txtfile="<<get_int_parameter("output_txtfile")<<endl;
  fout << "output_pepxmlfile="<<get_int_parameter("output_pepxmlfile")<<endl;
  fout << "output_outfiles="<<get_int_parameter("output_outfiles")<<endl;
  fout << "print_expect_score="<<get_int_parameter("print_expect_score")<<endl;
  fout << "num_output_lines="<<get_int_parameter("num_output_lines")<<endl;
  fout << "show_fragment_ions="<<get_int_parameter("show_fragment_ions")<<endl;
  fout << "sample_enzyme_number="<<get_int_parameter("sample_enzyme_number")<<endl;
  fout << "scan_range="<<get_string_parameter("scan_range")<<endl;
  fout << "precursor_charge="<<get_string_parameter("precursor_charge")<<endl;
  fout << "ms_level="<<get_int_parameter("ms_level")<<endl;
  fout << "activation_method="<<get_string_parameter("activation_method")<<endl;
  fout << "digest_mass_range="<<get_string_parameter("digest_mass_range")<<endl;
  fout << "num_results="<<get_int_parameter("num_results")<<endl;
  fout << "skip_researching="<<get_int_parameter("skip_researching")<<endl;
  fout << "max_fragment_charge="<<get_int_parameter("max_fragment_charge")<<endl;
  fout << "max_precursor_charge="<<get_int_parameter("max_precursor_charge")<<endl;
  fout << "nucleotide_reading_frame="<<get_int_parameter("nucleotide_reading_frame")<<endl;
  fout << "clip_nterm_methionine="<<get_int_parameter("clip_nterm_methionine")<<endl;
  fout << "spectrum_batch_size="<<get_int_parameter("spectrum_batch_size")<<endl;
  fout << "minimum_peaks="<<get_int_parameter("minimum_peaks")<<endl;
  fout << "minimum_intensity="<<get_double_parameter("minimum_intensity")<<endl;
  fout << "remove_precursor_peak="<<get_int_parameter("remove_precursor_peak")<<endl;
  fout << "remove_precursor_tolerance="<<get_double_parameter("remove_precursor_tolerance")<<endl;
  fout << "clear_mz_range="<<get_string_parameter("clear_mz_range")<<endl;
  fout << "variable_C_terminus="<<get_double_parameter("variable_C_terminus")<<endl;
  fout << "variable_N_terminus="<<get_double_parameter("variable_N_terminus")<<endl;
  fout << "variable_C_terminus_distance="<<get_int_parameter("variable_C_terminus_distance")<<endl;
  fout << "variable_N_terminus_distance="<<get_int_parameter("variable_N_terminus_distance")<<endl;
  fout << "add_Cterm_peptide="<<get_double_parameter("add_Cterm_peptide")<<endl;
  fout << "add_Nterm_peptide="<<get_double_parameter("add_Nterm_peptide")<<endl;
  fout << "add_Cterm_protein="<<get_double_parameter("add_Cterm_protein")<<endl;
  fout << "add_Nterm_protein="<<get_double_parameter("add_Nterm_protein")<<endl;
  
  fout << "add_G_glycine="<<get_double_parameter("add_G_glycine")<<endl;
  fout << "add_A_alanine="<<get_double_parameter("add_A_alanine")<<endl;
  fout << "add_S_serine="<<get_double_parameter("add_S_serine")<<endl;
  fout << "add_P_proline="<<get_double_parameter("add_P_proline")<<endl;
  fout << "add_V_valine="<<get_double_parameter("add_V_valine")<<endl;
  fout << "add_T_threonine="<<get_double_parameter("add_T_threonine")<<endl;
  fout << "add_C_cysteine="<<get_double_parameter("add_C_cysteine")<<endl;
  fout << "add_L_leucine="<<get_double_parameter("add_L_leucine")<<endl;
  fout << "add_I_isoleucine="<<get_double_parameter("add_I_isoleucine")<<endl;
  fout << "add_N_asparagine="<<get_double_parameter("add_N_asparagine")<<endl;
  fout << "add_D_aspartic_acid="<<get_double_parameter("add_D_aspartic_acid")<<endl;
  fout << "add_Q_glutamine="<<get_double_parameter("add_Q_glutamine")<<endl;
  fout << "add_K_lysine="<<get_double_parameter("add_K_lysine")<<endl;
  fout << "add_E_glutamic_acid="<<get_double_parameter("add_E_glutamic_acid")<<endl;
  fout << "add_M_methionine="<<get_double_parameter("add_M_methionine")<<endl;
  fout << "add_O_ornithine="<<get_double_parameter("add_O_ornithine")<<endl;
  fout << "add_H_histidine="<<get_double_parameter("add_H_histidine")<<endl;
  fout << "add_F_phenylalanine="<<get_double_parameter("add_F_phenylalanine")<<endl;
  fout << "add_R_arginine="<<get_double_parameter("add_R_arginine")<<endl;
  fout << "add_Y_tyrosine="<<get_double_parameter("add_Y_tyrosine")<<endl;
  fout << "add_W_tryptophan="<<get_double_parameter("add_W_tryptophan")<<endl;
  fout << "add_B_user_amino_acid="<<get_double_parameter("add_B_user_amino_acid")<<endl;
  fout << "add_J_user_amino_acid="<<get_double_parameter("add_J_user_amino_acid")<<endl;
  fout << "add_U_user_amino_acid="<<get_double_parameter("add_U_user_amino_acid")<<endl;
  fout << "add_X_user_amino_acid="<<get_double_parameter("add_X_user_amino_acid")<<endl;
  fout << "add_Z_user_amino_acid="<<get_double_parameter("add_Z_user_amino_acid")<<endl;
     
     
  /*enzyme information*/
    
  fout << "#"<<endl;
  fout << "# COMET_ENZYME_INFO _must_ be at the end of this parameters file" << endl;
  fout << "#" << endl;
  fout << "[COMET_ENZYME_INFO]" << endl;
  
  fout << "0.  No_enzyme"<<space;
  fout << 0;      
  fout << "       -           -" << endl;
  
  fout << "1.  Trypsin\t\t\t\t";
  fout << 1;
  fout << "      KR           P" << endl;
  
  fout << "2.  Trypsin/P\t\t\t\t";
  fout << 1;
  fout << "      KR           -" << endl;
  
  fout << "3.  Lys_C\t\t\t\t";
  fout << 1;
  fout << "      K            P" << endl;
  
  fout << "4.  Lys_N\t\t\t\t";
  fout << 0;
  fout << "      K            -" << endl;

  fout << "5.  Arg_C\t\t\t\t";
  fout << 1;
  fout << "      R            P" << endl;
  
  fout << "6.  Asp_N\t\t\t\t";
  fout << 0;
  fout << "      D            -" << endl;
  
  fout << "7.  CNBr\t\t\t\t";
  fout << 1;
  fout << "      M            -" << endl;

  fout << "8.  Glu_C"<<space;
  fout << 1;
  fout << "      DE           P" << endl;
  
  fout << "9.  PepsinA"<<space;
  fout << 1;
  fout << "      FL           P" << endl;

  fout << "10. Chymotrypsin"<<"\t\t\t";
  fout << 1;
  fout << "      FWYL         P" << endl;
   
  fout << "11. Elastase \t\t\t\t";
  fout << 1;
  fout << "      ALIV         P" << endl;

  fout << "12. Clostripai\t\t\t\t";
  fout << 1;
  fout << "      R            -" << endl;

  fout << "13. Iodosobenzoate\t\t\t";
  fout << 1;
  fout << "      W            -" << endl;

  fout << "14. Proline_Endopeptidase\t\t";
  fout << 1;
  fout << "      P            -" << endl;

  fout << "15. Staph_Protease\t\t\t";
  fout << 1;
  fout << "      E            -" << endl;

  fout << "16. Modified_Chymotrypsin\t\t";
  fout << 1;
  fout << "      FWYL         P" << endl;

  fout << "17. Elastase_Trypisn_Chymotrypsin\t";
  fout << 1;
  fout << "      ALIVKRWFY    P" << endl;

}
/**
 * \returns the command name for CometApplication
 */
string CometApplication::getName() {
  return "comet";
}

/**
 * \returns the description for CometApplication
 */
string CometApplication::getDescription() {

  return "Runs comet";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CometApplication::needsOutputDirectory() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
