/**
 * \file CometApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "Common.h"
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

  //get ms2 file 
  string spectra_file=get_string_parameter_pointer("input spectra");
                        
  // build argument list 
  vector<string> cmt_args_vec; 
  cmt_args_vec.push_back("comet");
  
  //pritn out comet params
  string param_file=make_file_path("comet.params");
  string param_file2=make_file_path("comet.params.txt");

  //Alterate ouput base name
  string basename =make_file_path(getName()); 
   
  
  ofstream fout;
  fout.open(param_file.c_str());
  if(fout.is_open()){
    carp(CARP_DEBUG,"start writing in file %s",param_file.c_str());
    writeParams(fout);
    fout.close();
  }else
    carp(CARP_FATAL,"Cann't open % file",param_file.c_str());
  //push param files 
  
  ostringstream param_oss;
  param_oss << "-P" << param_file2;
  cmt_args_vec.push_back(param_oss.str());

  ostringstream bsname_oss; 
  bsname_oss<<"-N"<<basename;  
  cmt_args_vec.push_back(bsname_oss.str());

  cmt_args_vec.push_back(spectra_file);

  /*build argv line*/
  int cmt_argc=cmt_args_vec.size();
  char** cmt_argv= new char*[cmt_argc];
  cmt_argv[0]= (char*)cmt_args_vec[0].c_str(); 
  
  for(unsigned idx=1; idx<cmt_argc;idx++){
    cmt_argv[idx]=(char*)cmt_args_vec[idx].c_str(); 
    
  }

  /* Re-route stderr to log file */
  CarpStreamBuf buffer;
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  /* Call comet_main */
  int retVal = -1;
  //retVal = comet_main(cmt_argc, cmt_argv);

  /* Recover stderr */
  std::cerr.rdbuf(old);

  delete []cmt_argv;

  return retVal;
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
