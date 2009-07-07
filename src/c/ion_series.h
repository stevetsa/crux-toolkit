/**
 * \file ion_series.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * $Revision: 1.23 $
 * \brief Object for a series of ions.
 *****************************************************************************/
#ifndef ION_SERIES_H
#define ION_SERIES_H

#ifdef __cplusplus
#include <vector>
#include <stdio.h>
#include "objects.h"
#include "peptide.h"
#include "ion.h"
#include "ion_series.h"
#include "spectrum.h"


#define BINARY_GMTK 1
#define PRINT_NULL_IONS 1
#define MIN_FRAMES 3
#define MAX_IONS 10000
#define MAX_NUM_ION_TYPE 8 // number of different ion_types
/**
 * \class ion_series
 * \brief An object to represent a series of ions, and organize them.
 * For which additional data structures will be created as needed 
 * loss_limit can be equal to NULL, thus if need to use should always
 * check that it is not NULL.
 */

typedef std::vector<ION_T*>::iterator ION_ITERATOR_T;





class ION_SERIES_T {
  /********************************************************
   *PRIVATE
   ********************************************************/
 private:
  // TODO change name to unmodified_char_seq

  std::vector<ION_T*> ions;//< The ions in this series
  //< the number of ions of a specific ion_type
  std::vector<ION_T*> specific_ions[MAX_NUM_ION_TYPE]; 
  

  char* peptide; //< The peptide sequence for this ion series
  MODIFIED_AA_T* modified_aa_seq; //< sequence of the peptide
  FLOAT_T peptide_mass; //< The peptide neutral mass. For efficiency. 
  int charge; //< /<The charge state of the peptide for this ion series
  ION_CONSTRAINT_T* constraint; //< The constraints which these ions obey

  BOOLEAN_T is_predicted; //< has this ion_series been predicted already?

  //< specific ions in the series, reference to master array of ions
  LOSS_LIMIT_T* loss_limit; 
  //< nh3, h2o loss limit for a given cleavage index, 
  //< before using this array should always sheck if not NULL
  int peptide_length;   //< the length of the peptide
  // ??? what is the difference between peptide_length and num_ions

  /**
   * \brief Find instances of amino acid which can incur neutral
   * losses: H2O (S|T|E|D), NH3(R|K|Q|N).  
   * Set the count of those observed so far for each cleavage index.
   * If no instance of amino acid, the count is assigned to 0
   * The information is used to determine if how many nh3 or h2o neutral
   * losses are possible. 
   */
  void scan_for_aa_for_neutral_loss();

  /**
   * \brief Creates an array in which element i is the sum of the masses
   * of amino acids 0 to (i-1).  At i=0 is stored the length of the
   * peptide.  
   * \returns an array of ion masses for all sub sequences
   */
  FLOAT_T* create_ion_mass_matrix(
    //char* peptide, ///< The peptide for this ion series. -in
    MODIFIED_AA_T* modified_seq, ///< the sequence
    MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
    int peptide_length ///< the length of the peptide
    );
  
  /**
   * user must ensure that there is enough space for this ion
   * adds ion to ion_series' master ion_array and if B|Y ion to the specific ion_array
   */
  void add_ion_to_ion_series(
    ION_T* ion ///< ion to add -in
    );

  /**
   * helper function: add_ions
   * add all the ions to ion_series up to the max charge
   *\returns TRUE if successfully adds all ions, else FALSE
   */

  BOOLEAN_T add_ions_by_charge(
    FLOAT_T mass, ///< the base mass of the ion to add
    int cleavage_idx, ///< the absolute cleavage index (A,B,C from left X,Y,Z from right)
    ION_TYPE_T ion_type ///< the ion type of the ions to be added
  );

  /**
   * Creates all the ions with no modifications up to the max charge
   * Adds each ion to ion_series
   *\returns TRUE if successfully generates all the ions, else FALSE
   */
  BOOLEAN_T generate_ions_no_modification(
    FLOAT_T* mass_matrix ///< the mass matrix that stores the mass
  );

  /**
   * The modification depends on the loss/add && if the ion contains RKQN or STED
   * The number of losses possible cannot exceed the number of RKQN or STED in the ion
   * The loss_limit array in the ion_series must be populated prior to this method call
   *\returns TRUE if the ion can lose the mod_type modification, else FALSE
   */
  BOOLEAN_T can_ion_lose_modification(
    ION_T* ion, ///< the ion to check if can lose nh3 -in
    ION_MODIFICATION_T mod_type, ///< generate ions of this modification_type -in/out
    int increment  ///< the add/loss of the modification
    );
  /**
   * creates all the ions with specific modifications up to the max charge
   * copies all the existing ions that can be modified,
   * then applies the different modifications then adds the new modified ions to ion_series
   *\returns TRUE if successfully generates all the ions with modifications, else FALSE
   */
  BOOLEAN_T generate_ions(
    ION_MODIFICATION_T mod_type ///< generate ions of this modification_type -in/out
    );

  /**
   * creates all the flanking ions up to the max charge
   * can only create flanking ions that are B|Y ions and don't have modification
   * assumes the ions with no modification all are at the begining of the ion[] in ion_series
   * copies all the existing ions that can be modified,
   * then applies the different modifications then adds the new modified ions to ion_series
   *\returns TRUE if successfully generates all the ions with modifications, else FALSE
   */
  BOOLEAN_T generate_ions_flank();

  /********************************************************
   *PUBLIC
   ********************************************************/
 public:

  int numIons();
  int numSpecificIons(ION_TYPE_T ion_type);
  
  ION_ITERATOR_T begin();
  ION_ITERATOR_T end();

  void init();
  /**
   * \returns An (empty) ion_series object.
   */
  ION_SERIES_T();
  
  /**
   * new_ion_series()
   * copies in the peptide sequence
   * Use this method to create ion_series only when few are needed,
   * because the memory allocation process is expensive.
   * If need a repeated new ion-series for different peptides, 
   * use "new_ion_series_generic" & "update_ion_series" combination, thus only allocate one 
   * ion_seires object.
   *\returns Instantiates a new ion_series object from the given peptide sequence and charge
   */
  ION_SERIES_T(char* peptide, ///< The peptide for this ion series. -in
    int charge, ///< The charge for this ion series -in
    ION_CONSTRAINT_T* constraint ///< The constraints which the ions in this series obey.
    );
  /**
   * new_ion_series_generic()
   * Creates a heap allocated generic ion_series object that must be updated by "update_ion_series" method
   * to transform the object into a ion-series for a specific instance of a peptide sequence and charge.
   *\returns Instantiates a new generic ion_series object that must be updated for each peptide instance
   */
  ION_SERIES_T(
    ION_CONSTRAINT_T* constraint, ///< The constraints which the ions in this series obey.
    int charge ///< The charge for this ion series -in
    );

  /**
   * free_ion_series()
   * Frees an allocated ion_series object.
   */
  virtual ~ION_SERIES_T();


  /**
   * Updates an ion_series to a specific instance of a peptide sequence.
   * If the ion_series has been already generated its ions, will free ions up.
   * Copies in the peptide sequence.
   * and re-initialize for the new peptide sequence.
   */
  void update_ion_series(
    char* peptide, ///< The peptide sequence for this ion series. -in
    MODIFIED_AA_T* mod_seq ///< modified version of seq -in
    );
  
  /**
   * Prints a ion_series object to file.
   */
  void print_ion_series(
    FILE* file ///< file for output -out
  );

  /**
   * Prints a ion_series object to file, in GMTK single-ion format.
   */
  void print_ion_series_single_gmtk(
    ION_CONSTRAINT_T* ion_constraint, ///< ion_constraint to obey -in 
    FILE* file, ///< file -out
    int sentence_idx
  );

  /**
   * Prints a ion_series object to file, in GMTK paired-ion format.
   */
  void print_ion_series_paired_gmtk(
    ION_CONSTRAINT_T* first_ion_constraint, ///< ion_constraint to obey -in 
    ION_CONSTRAINT_T* second_ion_constraint, ///< ion_constraint to obey -in 
    FILE* file, ///< file output
    int sentence_idx
  );

  /**
   * Predict ion series
   */
  void predict_ions();

  /**
   * Assign peaks to the nearest ions, within a tolerance (set in param file)
   */
  void ion_series_assign_nearest_peaks(SPECTRUM_T* spectrum);

  /**
   * Copies ion_series object from src to dest.
   *  must pass in a memory allocated ION_SERIES_T dest
   */
  static void copy_ion_series(
    ION_SERIES_T* src,///< ion to copy from -in
    ION_SERIES_T* dest///< ion to copy to -out
    );

  /** 
   * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
   */

  /*************************************
   * ION_SERIES_T: get and set methods
   ************************************/

  /**
   * \returns the ion that meets the constraint or NULL
   */
  ION_T* get_ion_series_ion(
    ION_CONSTRAINT_T* ion_constraint,
    int cleavage_idx
  );

  /**
   * User should not free the peptide sequence seperate from the ion_series
   *\returns a pointer to the original parent peptide sequence of the ion_series object
   */
  char* get_ion_series_peptide();

  /**
   *\returns the peptide length of which the ions are made
   */
  int get_ion_series_peptide_length();

  /**
   * copies in the peptide sequence to heap allocated sequence.
   * set the parent peptide sequence of the ion_series object
   */
  void set_ion_series_peptide(
    char* peptide///< the peptide sequence to set -in
    );

  /**
   *\returns the charge of the ion_series object
   */
  int get_ion_series_charge();

  /**
   * set the charge of the ion_series object
   */
  void set_ion_series_charge(
    int charge///< the charge of the ion -in
    );

  /**
   *\returns the constraint of the ion_series object
   */
  ION_CONSTRAINT_T* get_ion_series_ion_constraint();

  /**
   * set the of the ion_series object
   */
  void set_ion_series_ion_constraint(
    ION_CONSTRAINT_T* constraint///<  -in
    );

  /**
   *\returns the total number of ions in the ion_series object
   */
  int get_ion_series_num_ions();

  /**
   *\returns the total number of ion_type in the ion_series object
   */
  int get_ion_series_num_ions_one_type(
    ION_TYPE_T ion_type ///< the type of ions -in
    );


  /*FilteredIterator for ion series*/
  class FilteredIterator :
  public std::iterator<std::forward_iterator_tag, ION_T*> {
  protected:
    ION_SERIES_T* ion_series;
    ION_CONSTRAINT_T* constraint;
    void satisfyConstraint();
  public:
    ION_ITERATOR_T current;
    ION_ITERATOR_T end_iter;
    FilteredIterator();
    FilteredIterator(ION_SERIES_T* ion_series, ION_CONSTRAINT_T* constraint);
    ~FilteredIterator();
    
    FilteredIterator& operator=(const FilteredIterator& other);
    bool operator==(const FilteredIterator& other);
    bool operator!=(const FilteredIterator& other);
    bool operator!=(const ION_ITERATOR_T& other);
    FilteredIterator& operator++();
    FilteredIterator& operator++(int);
    
    //return a reference to the current ion pointer.
    ION_T*& operator*();
    ION_T* operator->();
    

    
    
  };


  FilteredIterator begin(ION_CONSTRAINT_T* constraint);

};

typedef ION_SERIES_T::FilteredIterator ION_FILTERED_ITERATOR_T;



/******************************/

/**
 * \class ion_constraint
 * \brief An object to represent the contraints which the ions in this
 * series obey.
 */

class ION_CONSTRAINT_T{
 private:
  BOOLEAN_T use_neutral_losses; //< Should ions include neutral losses
  int modifications[MAX_MODIFICATIONS]; 
  ///< an array to indicate which modifications to perform
  MASS_TYPE_T mass_type; 
    ///< the mass_type to use MONO|AVERAGE
  int max_charge; 
    ///< maximum charge of the ions, cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type; 
    ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion; 
    ///< does a precursor-ion satisfy this constraint
  int min_charge; 
  BOOLEAN_T exact_modifications; 
    ///< TRUE  = ints in modfications array indicate exact number of mods
    ///< FALSE = ints in modfications array indicate maximum number of mods
 public:
  unsigned int pointer_count;///< Number of pointers referencing me 
  void init();

  void init(
    MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
    int max_charge, ///< max charge of the ions <= the parent peptide's charge
    ION_TYPE_T ion_type, ///< the ion types the peptide series should include
    BOOLEAN_T precursor_ion  ///< should include precursor ion?
    ); 

  /**
   * allocate_ion_constraint()
   *\returns an empty heap allocated ion_constraint
   */
  ION_CONSTRAINT_T();
  
  /**
   * new_ion_constraint()
   * modification, all modifications 0
   * add more modifications as needed using the set_ion_constraint_modification
   *\returns a new heap allocated ion_constraint
   */
  ION_CONSTRAINT_T(
    MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
    int max_charge, ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
    ION_TYPE_T ion_type, ///< the ion types the peptide series should include
    BOOLEAN_T precursor_ion  ///< should include precursor ion?
  );

  /**
   *new_ion_constrain_smart()
   * \brief Create a new ion constraint based on the score type and the
   * charge of the peptide to be modeled.  Uses other
   * new_ion_constraint_ methods for some types.
   *
   * \returns A newly allocated ion constraint.
   */
  ION_CONSTRAINT_T(
    SCORER_TYPE_T score_type,
    int charge
  );
  
  /**
   *new_ion_constraint_sequest()
   * modification, sets all fields for sequest settings
   *\returns a new heap allocated ion_constraint
   */
  void init_ion_constraint_sequest(
    MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
    int max_charge, ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
    ION_TYPE_T ion_type, ///< the ion types the peptide series should include
    BOOLEAN_T precursor_ion ///< should include precursor ion?
    );

  /**
   * new_ion_constrain_gmtk
   * modification, sets all fields for GMTK settings
   *\returns a new heap allocated ion_constraint
   */
  void init_ion_constraint_gmtk(
    int charge ///< the charge of the peptide for which to predict ions
  );

  static ION_CONSTRAINT_T* new_gmtk(
    int charge ///< the charge of the peptide for which to predict ions
  );


  /**
   * modification, sets all fields for sequest Sp scoring settings
   *\returns a new heap allocated ion_constraint
   */
  void init_ion_constraint_sequest_sp(
    int max_charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  );

  /**
   * modification, sets all fields for Sequest Xcorr scoring settings
   * make B, Y, A type ions
   *\returns a new heap allocated ion_constraint
   */
  void init_ion_constraint_sequest_xcorr(
    int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
    );

  virtual ~ION_CONSTRAINT_T();
  static void free(ION_CONSTRAINT_T*);
  /**
   * copies the ion_constraint pointer
   */
  static ION_CONSTRAINT_T* copy_ion_constraint_ptr(ION_CONSTRAINT_T* constraint);


  /**
   * copies ion_constraint object from src to dest
   * must pass in a memory allocated ION_CONSTRAINT_T dest
   */
  static void copy_ion_constraint(
    ION_CONSTRAINT_T* src,///< ion_constraint to copy from -in
    ION_CONSTRAINT_T* dest///< ion_constraint to copy to -out
    );

  /** 
   * Determines if a ion satisfies a ion_constraint.
   * \returns TRUE if the constraint is satisified. FALSE if not.
   */
  BOOLEAN_T is_satisfied(
   ION_T* ion ///< query ion -in
   );

  /**
   * Sets the modification count
   * can only add isotopes
   */
  void set_ion_constraint_modification(
    ION_MODIFICATION_T mod_type, ///< ion modification type -in
    int count  ///< the count of the modification -in  
    );

  /**
   * Sets the exact modification boolean to exactness criteria
   * and if exactness is true sets min_charge = max_charge.
   * In other words, the constraint is now exact, in that it refers to a
   * particular ion series, charge states, and modification state, as opposed
   * to e.g. b-ions of charge state +1 or +2, or with or without NH3 loss
   */
  void set_ion_constraint_exactness(
    BOOLEAN_T exactness ///< whether to be exact or not -in
    );

  /**
   * gets the modification count for specific mod_type
   */
  int get_modification(
    ION_MODIFICATION_T mod_type ///< ion modification type -in
    );

  /**
   * gets the modification count array for specific mod_types
   */
  int* get_modifications();

  /**
   * gets the mass type of the ion_constraint
   */
  MASS_TYPE_T get_mass_type();

  /**
   * get maximum charge of the ions, cannot exceed the parent peptide's charge
   */
  int get_max_charge();

  /**
   * get the ion types the peptide series should include
   */
  ION_TYPE_T get_ion_type(); 

  /*
   * get whether a precursor-ion satisfies this constraint
   */
  BOOLEAN_T get_precursor_ion();
  
  /*
   * get if ions should include neutral losses
   */
  BOOLEAN_T get_use_neutral_losses();



};


#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */



#endif
