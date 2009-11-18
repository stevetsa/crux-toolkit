#ifndef CRUX_PROTEIN_H_
#define CRUX_PROTEIN_H_


extern "C" {
#include "protein.h"

/**
 * Creates the data structures in the protein_peptide_iterator object necessary
 * for creating peptide objects.
 * - mass_array - cumulative distribution of masses. used to determine 
 *     the mass of any peptide subsequence.
 * - nterm_cleavage_positions - the nterm cleavage positions of the 
 *     peptides that satisfy the protein_peptide_iterator contraints
 * - peptide_lengths - the lengths of the peptides that satisfy the constraints
 * - cumulative_cleavages - cumulative distribution of cleavage positions
 *    used to determine if a cleavage location has been skipped
 */
void prepare_protein_peptide_iterator(
    PROTEIN_PEPTIDE_ITERATOR_T* iterator
    );
}

#endif
