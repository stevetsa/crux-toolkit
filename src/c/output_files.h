class OutputFiles {
 public:
  OutputFiles(int num_decoys, int num_proteins)
    : num_decoys_(num_decoys),
    psm_file_array_(NULL),
    sqt_file_(NULL),
    decoy_sqt_file_(NULL),
    tab_file_(NULL),
    decoy_tab_file_(NULL) {
    Open(num_proteins);
  }

  void Close(int num_spectra);

  void PrintMatches(MATCH_COLLECTION_T* match_collection, SPECTRUM_T* spectrum,
		    bool is_decoy = false, int psm_file_index = 0,
		    bool send_to_sqt_tab = true);
 private:
  void Open(int num_proteins);

  int num_decoys_;

  FILE** psm_file_array_; ///< binary psm filehandles
  FILE* sqt_file_;
  FILE* decoy_sqt_file_;
  FILE* tab_file_;
  FILE* decoy_tab_file_;
};
