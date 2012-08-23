#ifndef SPECTRUMSCOREDISTRIBUTION_H
#define SPECTRUMSCOREDISTRIBUTION_H

#include <stdio.h>
#ifndef _MSC_VER
#include <dirent.h>
#endif
#include <string>
#ifdef _MSC_VER
#include "windirent.h"
#endif
#include "objects.h"
#include "Spectrum.h"
#include "Peptide.h"
#include "Ion.h"

#include <limits.h>

using namespace Crux;

class SpectrumScoreDistribution {

private:

	/**
	 * Quantize range [lower, upper] in increments of delta.
	 *
	 * \side ovalues allocated and initialized.
	 * \side len initialized to length of ovalues.
	 */
	static void quantize(
		FLOAT_T lower, ///< lower value of range -in
		FLOAT_T upper, ///< upper value of range -in
		FLOAT_T delta, ///< increment of values in range -in
		FLOAT_T* ovalues, ///<
		int& olen
		);

protected:
	FLOAT_T** table_; ///< table_[s][m] = num of peptides of mass m with score m
	int nrows_; ///< number of rows (quantized scores) in table_ 
	int ncols_; ///< number of columns (masses) in table_
	int lower_; ///< index of column whose mass is neutral mass - tolerance
	int upper_; ///< index of column whose mass is neutral mass + tolerance

	FLOAT_T* mass_; ///< mass_[m] is mass represented by table_[...][m]
	FLOAT_T* score_; ///< scores_[s] is the score represented by table_[s][...]

	FLOAT_T* observed_; ///< observed_[m] FastSEQUEST processed spctrum (real-valued)
	int maxbin_; ///< length of observed_

	bool failed_;
	FLOAT_T neutralmass_;

	static const double epsilon_;
	int* aamass_;

	/**
	 * Deallocate all memory safely.
	 */
	void deallocate();

	/**
	 * Compute the distribution of peptides of a particular mass using dynamic
	 * programming.
	 *
	 * \side table_ is populated, and the class is either ready to convert xcorr -> pvalue
	 * or is initialized() is false.
	 */
	void computeScores();

	/**
	 * Count the fraction of peptides which have a score higher than xcorr.
	 */
	void countHigherScoring(FLOAT_T xcorr, FLOAT_T& nBetter, FLOAT_T& nPeptides) const;


public:

	SpectrumScoreDistribution(
		Spectrum* spectrum, ///< spectrum to compute score distribution against -in
		const SpectrumZState& zstate ///< charge state of the spectrum -in
		);

	~SpectrumScoreDistribution();

	/**
	 * Is the object correctly initialized, and ready to convert xcorr -> pvalue?
	 */
	bool initialized() const;

	/**
	 * Did the dynamic program fail for some reason?
	 */
	bool failed() const;

	/**
	 * Largest score considered by the dynamic program.
	 */
	FLOAT_T maxScore() const;

	/**
	 * Largest neutral mass of a peptide considered by the dynamic program.
	 */
	FLOAT_T maxMass() const;

	/**
	 * Convert xcorr to p-value.
	 *
	 * \returns A p-value, a floating point number in [0,1] where lower means
	 *   more confident about the match.
	 */
	FLOAT_T pvalue(FLOAT_T xcorr) const;

	// log base-e.
	FLOAT_T logpvalue(FLOAT_T xcorr) const;
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif