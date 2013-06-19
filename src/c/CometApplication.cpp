/**
 * \file CometApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "CometApplication.h"
#include "DelimitedFileWriter.h"
#include "DelimitedFile.h"
#include "Common.h"
#include "ModifiedPeptidesIterator.h"

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
  const char* argument_list[] = {"input spectra","protein-database"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);
  /*Get input file*/
  //1.protein database 
  string protein_dbase="/some/path/db.fasta";
  if( get_string_parameter("protein database")!= "_NULL_STR")
    protein_dbase= get_string_parameter("protein database");

  //2.get ms2 file 
  string spectra_file=get_string_parameter_pointer("input spectra");
                        
  // build argument list 
  vector<string> cmt_args_vec; 
  cmt_args_vec.push_back("comet");
  
  //pritn out comet params
  string param_file=make_file_path("comet.params");
 

  //Alterate ouput base name
  string basename =make_file_path(getName()); 
   
  
  ofstream fout;
  fout.open(param_file.c_str());
  if(fout.is_open()){
    carp(CARP_DEBUG,"start writing in file %s",param_file.c_str());
    writeParams(fout, protein_dbase);
    fout.close();
  }else
    carp(CARP_FATAL,"Cann't open % file",param_file.c_str());
  //push param files 
  
  ostringstream param_oss;
  param_oss << "-P" << param_file;
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
 
  return comet_main(cmt_argc, cmt_argv);
 
}


/**
 * \write parameters 
 */
void CometApplication:: writeParams(ofstream &fout, string protein_dbase){
  string space="\t\t\t\t";
  string header_file="# Comet MS/MS search engine parameters file."
  "\n# Everything following the \'#\' symbol is treated as a comment.\n";
  fout<<header_file<<endl;

  //database 
  fout<<"database_name = "<<protein_dbase<<"\n";
  //decoy_search 
  fout<<"decoy_search = " << 1;
  /*
  if(get_int_parameter("num-decoys-per-target")== 0)
    fout<<"decoy_search = " << 0;
  else if(get_int_parameter("num-decoys-per-target")== 1){
    string decoy_location= get_string_parameter("decoy-location");
    if(decoy_location == "target-file")
       fout<<"decoy_search = "<< 1;
    else if( decoy_location== "one-decoy-search" || decoy_location == "separate-decoy-files")
       fout<<"decoy_search = "<< 2;
    
  }else {
    carp(CARP_FATAL,"Invalid num-decoys-per-target %s", get_int_parameter("num-decoys-per-target"));
  }*/
  fout<<"                       # 0=no (default), 1=concatenated search, 2=separate search"<<endl;
    
  fout<<"num_threads = "<<get_int_parameter("num-threads");
  fout<<"                        # 0=poll CPU to set num threads; else specify num threads directly (max 32)"<<endl;

  /*mass*/
  fout <<""<<endl;
  fout <<"#"<<endl;
  fout <<"# masses"<<endl;
  fout <<"#"<<endl;
  fout <<"peptide_mass_tolerance = "<< get_double_parameter("precursor-window")<<endl;
  WINDOW_TYPE_T window_type = get_window_type_parameter("precursor-window-type");
  fout <<"peptide_mass_units = ";
  
  switch (window_type) {
    case WINDOW_MASS:
      fout << "0";
      break;
    case WINDOW_PPM:
      fout << "2";
      break;
    case WINDOW_INVALID:
    case WINDOW_MZ:
    case NUMBER_WINDOW_TYPES:
      carp(CARP_FATAL, "Unsupported precursor window type : %s", get_string_parameter_pointer("precursor-window-type"));
  }
  fout<<"\t\t"<<"       "<<"# 0=amu, 1=mmu, 2=ppm "<<endl;

    fout <<"mass_type_parent = "<<get_mass_type_parameter("isotopic-mass");
    fout<<"                   # 0=average masses, 1=monoisotopic masses "<<endl;

    fout <<"mass_type_fragment = "<<get_mass_type_parameter("fragment-mass");
    fout<<"                 # 0=average masses, 1=monoisotopic masses"<<endl;
    fout <<"precursor_tolerance_type = 0           # 0=MH+ (default), 1=precursor m/z"<<endl;
    fout <<"isotope_error = 0                      # 0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)"<<endl;

    /*search enzyme*/
    fout<<endl;
    fout <<"#"<<endl;
    fout <<"# search enzyme"<<endl;
    fout <<"#"<<endl;
   
    //print enzyme 
    ENZYME_T enzyme; //2.question 
    enzyme= get_enzyme_type_parameter("enzyme");

    fout<<"search_enzyme_number = ";
    map<ENZYME_T,int> enzyme_map;

    enzyme_map [NO_ENZYME ] = 0; 
    enzyme_map [TRYPSIN] = 1;
    enzyme_map [TRYPSINP] = 2;
    enzyme_map [LYSC]= 3 ;
    enzyme_map [LYSN]= 4 ;
    enzyme_map [ARGC]= 5 ;
    enzyme_map [ASPN]= 6 ;
    enzyme_map [CYANOGEN_BROMIDE]= 7 ;
    enzyme_map [GLUC]= 8 ;
    enzyme_map [PEPSINA]= 9 ; 
    enzyme_map [MODIFIED_CHYMOTRYPSIN]= 10 ;
    enzyme_map [ELASTASE]= 11 ;
    enzyme_map [CLOSTRIPAIN]= 12 ;
    enzyme_map [IODOSOBENZOATE]= 13 ;
    enzyme_map [PROLINE_ENDOPEPTIDASE]= 14 ;
    enzyme_map [STAPH_PROTEASE]= 15 ;
    enzyme_map [CHYMOTRYPSIN]= 16 ;
    enzyme_map [ELASTASE_TRYPSIN_CHYMOTRYPSIN]= 17 ;
    enzyme_map [CUSTOM_ENZYME]= 18 ;
    enzyme_map [INVALID_ENZYME]= 19 ;
    
    map<ENZYME_T,int>::iterator iter; 
    for(iter =enzyme_map.begin(); iter != enzyme_map.end(); iter++ ){
      if((*iter).first == enzyme){
        if((*iter).first == INVALID_ENZYME)
          carp(CARP_FATAL, "Invalid enzyme : %s", get_string_parameter("enzyme"));
          
        fout<< (*iter).second; 
        break; 
      }
    }
    
    fout<<"               "<<"# choose from list at end of this params file"<<endl;
    //digestion
    if(get_digest_type_parameter("digestion")==FULL_DIGEST)    
      fout <<"num_enzyme_termini = "<<"2";
    else if(get_digest_type_parameter("digestion")==PARTIAL_DIGEST)
      fout <<"num_enzyme_termini = "<<"1";
    else
      carp(CARP_FATAL,"Digestion value is not accaptable");
    fout<<"                 # valid values are 1 (semi-digested), 2 (fully digested, default), 8 N-term, 9 C-term"<<endl;
    //missed-cleavages
    fout <<"allowed_missed_cleavage = "<<get_int_parameter("missed-cleavages");
    fout<<"            # maximum value is 5; for enzyme search"<<endl;
    

    /*variable mod*/
    
    string filename = get_string_parameter("protein database");
    bool use_index = is_directory(filename.c_str());
    Index* index = NULL;
    Database* database = NULL;

    if( use_index == true ){
      index = new Index(filename.c_str()); 
    }else{
      database = new Database(filename.c_str(), false); // not memmapped
    }

    AA_MOD_T** aa_mods =NULL; 
    int num_aa_mods = get_all_aa_mod_list(&aa_mods);
   
        
    fout << endl;
    fout << "#" << endl;
    fout << "# Up to 6 variable modifications are supported" << endl;
    fout << "# format:  <mass> <residues> <0=variable/1=binary> <max mods per a peptide>" << endl;
    fout << "#     e.g. 79.966331 STY 0 3" << endl;
    fout << "#" << endl;

    int count = 1;
    for (int mod_idx = 0; mod_idx < num_aa_mods; mod_idx++) {
      AA_MOD_T* mod = aa_mods[mod_idx];
      if (aa_mod_get_position(mod) == ANY_POSITION && mod_idx < 7) {
        int max_per_peptide = aa_mod_get_max_per_peptide(mod);
        FLOAT_T mass_change = aa_mod_get_mass_change(mod);
        char* residues = aa_mod_get_aa_list_string(mod);
        fout <<"variable_mod"<<(count)<<" = " <<
               mass_change << " " << residues << " 0 " << max_per_peptide << endl;
        count++;
      }
    }
    fout <<"max_variable_mods_in_peptide = "<<get_int_parameter("max-mods")<<endl;

    /*fragment ions*/
    fout << endl;
    fout << "#"<< endl;
    fout << "# fragment ions"<< endl;
    fout << "#"<< endl;
    fout << "# ion trap ms/ms:  0.36 tolerance, 0.11 offset (mono masses)"<<endl;
    fout << "# high res ms/ms:  0.01 tolerance, 0.00 offset (mono masses)"<<endl;
    fout << "#"<< endl;
   
    fout << "fragment_bin_tol = "<<get_double_parameter("mz-bin-width");
    fout<< "             # binning to use on fragment ions"<<endl;

    fout << "fragment_bin_offset = "<<get_double_parameter("mz-bin-offset");
    fout<< "             # offset position to start the binning"<<endl;

    fout << "theoretical_fragment_ions = 0          # 0=default peak shape, 1=M peak only"<<endl;
    fout << "use_A_ions = " <<get_boolean_parameter("use-a-ions")<<endl;
    fout << "use_B_ions = " <<get_boolean_parameter("use-b-ions")<<endl;
    fout << "use_C_ions = " <<get_boolean_parameter("use-c-ions")<<endl;
    fout << "use_X_ions = " <<get_boolean_parameter("use-x-ions")<<endl;
    fout << "use_Y_ions = " <<get_boolean_parameter("use-y-ions")<<endl;
    fout << "use_Z_ions = " <<get_boolean_parameter("use-z-ions")<<endl;
    fout << "use_NL_ions = "<<get_boolean_parameter("use-nl-ions");
    fout << "                        # 0=no, 1=yes to consider NH3/H2O neutral loss peaks"<<endl;

    /*output*/
    fout << endl;
    fout << "#"<<endl;
    fout << "# output"<<endl;
    fout << "#"<<endl;

    fout << "output_format = 2                      # 0=sqt stdout (default), 1=sqt file, 2=out files, 3=pepXML file"<<endl;
    fout << "output_sqtfile = 1                     # 0=no, 1=yes write sqt file"<<endl;//always true 
    fout << "output_pepxmlfile = 1                  # 0=no, 1=yes write pep.xml file"<<endl;//always true 
    fout << "output_outfiles = 0                    # 0=no, 1=yes write pep.xml file"<<endl;//always false
 
    fout << "print_expect_score = "<<get_boolean_parameter("print-expect-score");
    fout<<"                 # 0=no, 1=yes to replace Sp with expect"<<endl;

    fout << "num_output_lines = "<<get_int_parameter("top-match");
    fout << "                   # num peptide results to show" << endl;

    fout << "show_fragment_ions = 0                 # 0=no, 1=yes"<<endl; //always false
    fout << ""<<endl;
    fout << "sample_enzyme_number = " << get_int_parameter("sample-enzyme");
    fout <<"               # Sample enzyme which is possibly different than the one applied to the search."<<endl;
    fout << "                                       # Used to calculate NTT & NMC in pepXML output (default=1 for trypsin)."<<endl;
    fout << endl;

    /*mzXML*/
    fout <<""<<endl;
    fout <<"#"<<endl;
    fout <<"# mzXML parameters"<<endl;
    fout <<"#"<<endl;
 
    fout <<"scan_range = " << get_int_parameter("start-scan") << " " ;
    fout << get_int_parameter("end-scan");
    fout <<"                       # start and scan scan range to search; 0 as 1st entry ignores parameter"<<endl;

    fout <<"precursor_charge = 0 0";
    fout <<"                 # precursor charge range to analyze; does not override mzXML charge;"; 
    fout <<"0 as 1st entry ignores parameter"<<endl;

    fout <<"ms_level = "<<get_int_parameter("ms-level");
    fout <<"                           # MS level to analyze, valid are levels 2 (default) or 3"<<endl;

    fout <<"activation_method = "<<get_string_parameter("activation-method");
    fout <<"                # activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD"<<endl;
    
    fout <<""<<endl;
    fout <<"#"<<endl;
    fout <<"# misc parameters"<<endl;
    fout <<"#"<<endl;

    fout <<"digest_mass_range = "<< get_double_parameter("min-mass");;
    fout<<" "<< get_double_parameter("max-mass");;
    fout <<"           # MH+ peptide mass range to analyze"<<endl;

    fout <<"num_results = 50"; //default it to 50 
    fout <<"                       # number of search hits to store internally"<<endl;

    fout <<"skip_researching = 0 "; //default it to 0 
    fout <<"                  # for '.out' file output only, 0=search everything again (default), 1=don't search if .out exists"<<endl;
    
    int max_ion_charge = 3; 
    if(max_ion_charge>5)
      carp(CARP_FATAL,"max_fragment_charge cann't be more than 5.");
    else 
      fout <<"max_fragment_charge = "<< max_ion_charge; 
    fout <<"                # set maximum fragment charge state to analyze (allowed max 5)"<<endl;
    
    int max_precursor_charge =  get_int_parameter("max-precursor-charge");
    if(max_precursor_charge > 9)
      carp(CARP_FATAL, "max_precusor_charge cann't be more than 9.");
    else 
      fout <<"max_precursor_charge = "<< max_precursor_charge;
    fout <<"               # set maximum precursor charge state to analyze (allowed max 9)"<<endl; //make sure not more than 9

    fout <<"nucleotide_reading_frame = 0"; //sean will fix it 
    fout <<"           # 0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six"<<endl;

    fout <<"clip_nterm_methionine = 0";
    fout <<"              # 0=leave sequences as-is; 1=also consider sequence w/o N-term methionine"<<endl;
    
    fout <<""<<endl;
    fout <<"#"<<endl;
    fout <<"# spectral processing"<<endl;
    fout <<"#"<<endl;
    fout <<"minimum_peaks = "<< get_int_parameter("min-peaks");
    fout <<"                     # minimum num. of peaks in spectrum to search (default 5)"<<endl;

    fout <<"remove_precursor_peak = "<<get_int_parameter("remove-precursor-peak");
    fout <<"              # 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)"<<endl;

    fout <<"remove_precursor_tolerance = "<<get_double_parameter("remove-precursor-tolerance");
    fout <<"       # +- Da tolerance for precursor removal"<<endl;
    fout<<endl;

    /*additional modifications*/ 
    fout <<"#"<<endl;
    fout <<"# additional modifications"<<endl;
    fout <<"#"<<endl;
    fout <<""<<endl;
    //question : send email to => ask Jimmy 
   
    AA_MOD_T** c_mods = NULL; 
    AA_MOD_T** n_mods = NULL; 
    int num_c_mods = get_c_mod_list(&c_mods);
    int num_n_mods = get_n_mod_list(&n_mods);
    
    //print c_mod
   
    for (int cmod_idx = 0; cmod_idx < num_c_mods; cmod_idx++) {
      AA_MOD_T* cmod = c_mods[cmod_idx];
     
      if (aa_mod_get_position(cmod) == C_TERM) {
        fout <<"variable_C_terminus = "<< aa_mod_get_mass_change(cmod) << endl;
        fout <<"variable_C_terminus_distance = " << aa_mod_get_max_distance(cmod);         
        fout <<"      # -1=all peptides, 0=protein terminus, 1-N = maximum offset from C-terminus"<<endl;
        
      }
     
    }
    
    //print n_mod 
    for (int nmod_idx = 0; nmod_idx < num_c_mods; nmod_idx++) {
      AA_MOD_T* nmod = n_mods[nmod_idx];
      if (aa_mod_get_position(nmod) == N_TERM) {
        fout <<"variable_N_terminus = "<< aa_mod_get_mass_change(nmod)<< endl;
        int n_distance = aa_mod_get_max_distance(nmod);
        fout <<"variable_N_terminus_distance = " << aa_mod_get_max_per_peptide(nmod);         
        fout <<"       # -1=all peptides, 0=protein terminus, 1-N = maximum offset from N-terminus"<<endl;        
      } 
    }   
    fout <<endl;


    fout <<"add_Cterm_peptide = 0.0"<<endl;
    fout <<"add_Nterm_peptide = 0.0"<<endl;
    fout <<"add_Cterm_protein = 0.0"<<endl;
    fout <<"add_Nterm_protein = 0.0"<<endl;
    fout <<endl;

    fout <<"add_G_glycine ="<< get_double_parameter("G");
    fout <<"                       # added to G - avg.  57.0513, mono.  57.02146"<<endl;

    fout <<"add_A_alanine = "<<get_double_parameter("A");
    fout <<"                      # added to A - avg.  71.0779, mono.  71.03711"<<endl;

    fout <<"add_S_serine = "<<get_double_parameter("S");
    fout <<"                       # added to S - avg.  87.0773, mono.  87.02303"<<endl;

    fout <<"add_P_proline = "<<get_double_parameter("P");
    fout <<"                      # added to P - avg.  97.1152, mono.  97.05276"<<endl;

    fout <<"add_V_valine = "<<get_double_parameter("V");
    fout <<"                       # added to V - avg.  99.1311, mono.  99.06841"<<endl;

    fout <<"add_T_threonine = "<<get_double_parameter("T");
    fout <<"                    # added to T - avg. 101.1038, mono. 101.04768"<<endl;

    fout <<"add_C_cysteine = "<<get_double_parameter("C");
    fout <<"               # added to C - avg. 103.1429, mono. 103.00918"<<endl;

    fout <<"add_L_leucine = "<<get_double_parameter("L");
    fout <<"                      # added to L - avg. 113.1576, mono. 113.08406"<<endl;

    fout <<"add_I_isoleucine = "<<get_double_parameter("I");
    fout <<"                   # added to I - avg. 113.1576, mono. 113.08406"<<endl;

    fout <<"add_N_asparagine = "<<get_double_parameter("N");
    fout <<"                   # added to N - avg. 114.1026, mono. 114.04293"<<endl;

    fout <<"add_D_aspartic_acid = "<<get_double_parameter("D");
    fout <<"                # added to D - avg. 115.0874, mono. 115.02694"<<endl;

    fout <<"add_Q_glutamine = "<<get_double_parameter("Q");
    fout <<"                    # added to Q - avg. 128.1292, mono. 128.05858"<<endl;

    fout <<"add_K_lysine =  "<<get_double_parameter("K");
    fout <<"                      # added to K - avg. 128.1723, mono. 128.09496"<<endl;

    fout <<"add_E_glutamic_acid =  "<<get_double_parameter("E");
    fout <<"               # added to E - avg. 129.1140, mono. 129.04259"<<endl;

    fout <<"add_M_methionine =  "<<get_double_parameter("M");
    fout <<"                  # added to M - avg. 131.1961, mono. 131.04048"<<endl;

    fout <<"add_O_ornithine =  "<<get_double_parameter("O");
    fout <<"                   # added to O - avg. 132.1610, mono  132.08988"<<endl;
          
    fout <<"add_H_histidine = "<<get_double_parameter("H");
    fout <<"                    # added to H - avg. 137.1393, mono. 137.05891"<<endl;

    fout <<"add_F_phenylalanine = "<<get_double_parameter("F");
    fout <<"                # added to F - avg. 147.1739, mono. 147.06841"<<endl; 

    fout <<"add_R_arginine = "<<get_double_parameter("R");
    fout <<"                     # added to R - avg. 156.1857, mono. 156.10111"<<endl;

    fout <<"add_Y_tyrosine = "<<get_double_parameter("Y");
    fout <<"                     # added to Y - avg. 163.0633, mono. 163.06333"<<endl;

    fout <<"add_W_tryptophan = "<<get_double_parameter("W");
    fout <<"                   # added to W - avg. 186.0793, mono. 186.07931"<<endl;

    fout <<"add_B_user_amino_acid = "<<get_double_parameter("B");
    fout <<"              # added to B - avg.   0.0000, mono.   0.00000"<<endl;

    fout <<"add_J_user_amino_acid = "<<get_double_parameter("J");
    fout <<"              # added to J - avg.   0.0000, mono.   0.00000"<<endl;
    fout <<"add_U_user_amino_acid = "<<get_double_parameter("U");
    fout <<"              # added to U - avg.   0.0000, mono.   0.00000"<<endl;

    fout <<"add_X_user_amino_acid = "<<get_double_parameter("X");
    fout <<"              # added to X - avg.   0.0000, mono.   0.00000"<<endl;

    fout <<"add_Z_user_amino_acid = "<<get_double_parameter("Z");
    fout <<"              # added to Z - avg.   0.0000, mono.   0.00000"<<endl;
    fout <<endl;
     
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

    fout << "18. Custom_Enzyme"<<"\t\t\t";
    fout << 1;
    fout << "      x            x" << endl;

    fout << "19. Invalid_Enzyme"<<"\t\t\t";
    fout << 1;
    fout << endl;

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
