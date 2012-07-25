/*************************************************************************//**
 * \file SpectrumScoreDistribution.cpp
 * AUTHOR: Ajit P. Singh
 * CREATE DATE: 1 Jul 2012
 * \brief Object for converting XCorr to a per-spectrum p-value, using
 * dynamic programming.
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifndef _MSC_VER
#include <dirent.h>
#include <unistd.h>
#endif
#include <ctype.h>
#include <sys/stat.h>
#ifdef _MSC_VER
#include "windirent.h"
#endif
#include "objects.h"
#include "IonConstraint.h"
#include "IonFilteredIterator.h"
#include "IonSeries.h"
#include "crux-utils.h"
#include "Spectrum.h"
#include "Scorer.h"
#include "parameter.h"

#include "SpectrumScoreDistribution.h"

void SpectrumScoreDistribution::quantize(
	FLOAT_T lower,
	FLOAT_T upper,
	FLOAT_T delta,
	FLOAT_T* ovalues,
	int& olen
	)
{
	assert(lower < upper);
	assert(delta > 0.0);
	assert(ovalues == NULL);

	olen = ceil((upper-lower)/delta);
	ovalues = (FLOAT_T*)mycalloc(olen, sizeof(FLOAT_T));
	for (int i = 0; i < olen; ++i)
		ovalues[i] = i * delta;
}


SpectrumScoreDistribution::SpectrumScoreDistribution(
	Spectrum* spectrum,
	const SpectrumZState& zstate
	) : table_(NULL), nrows_(-1), ncols_(-1), mass_(NULL), score_(NULL), 
			observed_(NULL), maxbin_(-1), failed_(false)
{
	if ( spectrum == NULL )
		carp(CARP_FATAL, "Cannot compute score distribution, spectrum is NULL.");

	int charge = zstate.getCharge();
	double tol = get_double_parameter("precursor-window");
	FLOAT_T maxMass = ceil(charge * zstate.getNeutralMass() + tol);
	FLOAT_T massDelta = 1.0;
	quantize(0.0, maxMass, massDelta, mass_, ncols_);

	FLOAT_T maxScore = 10.0; // TODO(ajit): Do something smarter than fixing upper bound.
	FLOAT_T scoreDelta = 0.1;
	quantize(0.0, maxScore, scoreDelta, score_, nrows_);

	// Allocate memoization table.
	table_ = (FLOAT_T**)mycalloc(nrows_, sizeof(FLOAT_T*));
	for (int i = 0; i < nrows_; ++i)
		table_[i] = (FLOAT_T*)mycalloc(ncols_, sizeof(FLOAT_T));

	// Apply FastSEQUEST to the observed spectrum, storing into observed_.
	Scorer* scorer = new Scorer(XCORR);
  if(!scorer->createIntensityArrayObserved(spectrum, charge)) {
  	deallocate();
  	free(scorer);
    carp(CARP_ERROR, "Failed to preprocess observed spectrum for Xcorr.");
    return;
  }
  maxbin_ = scorer->getMaxBin();
  assert(maxbin_ > 0);
  observed_ = (FLOAT_T*)mycalloc(maxbin_, sizeof(FLOAT_T));
  memcpy(observed_, scorer->getIntensityArrayObserved(), maxbin_*sizeof(FLOAT_T));
  free(scorer);

  assert(nrows_ > 0);
  assert(ncols_ > 0);
  assert(mass_);
  assert(score_);
  assert(observed_);
  assert(maxbin_ > 0);
  assert(maxbin_ <= ncols_); // Max mass in DP table > mass of spectrum.

  computeScores();
}

void SpectrumScoreDistribution::deallocate() 
{
	for (int i = 0; i < nrows_; ++i)
		free(table_[i]);
	free(table_); table_ = NULL;
	free(mass_); mass_ = NULL;
	free(score_); score_ = NULL;
	free(observed_); observed_ = NULL;
}

/**
 * Compute all the values of table_[s][m] using dynamic programming from left
 * to right in table.
 *
 * table_[s][m] &=
 *   \begin{cases}
 *   1, if s = 0 \& m = 0,
 *   0, if s < 0,
 *   0, if s = 0 \& m > 0,
 *   \sum_{a in amino acids} table_[s - contrib[m]][m - mass(a)]
 *   \end{cases}
 *
 * where contrib[m] is the contribution of an ion at mass m. We only assume that
 * we are dealing with b-ions and y-ions: i.e., that "neutral-losses=none" in the
 * Crux parameters.
 *
 * contrib[m] is defined as observed_[m]/10000.0 (see Scorer::crossCorrelation).
 * We opt not to distinguish between b- and y-ions.
 */
void SpectrumScoreDistribution::computeScores()
{
	static const int mass[20] = { 71, 160, 129, 115, 57, 147, 113, 137, 128, 131, 113, 
																114, 128, 97, 87, 156, 101, 186, 99, 163 };

	// The only non-zero value in table_[0][*] is the first one.
	table_[0][0] = 1;
	for (int m = 0; m < ncols_; ++m) {
		for (int s = 0; s < nrows_; ++s) {
			if (s == 0) {
				continue;
			}
			// Recurrence
			int v = 0;
			for (int a = 0; a < 20; ++a) {
				int new_m = m - mass[a];
				int new_s = static_cast<int>(s - (50.0 * observed_[new_m])/10000.0);
				//int new_s = static_cast<int>(s - (observed_[new_m]/10000.0));
				if (new_s < 0 || (new_s == 0 && new_m > 0)) {
					// table_[new_s][new_m] will be zero, so add nothing to v.
				}
				else {
					v += table_[new_s][new_m];
				}
			}
			assert(v >= 0);
			table_[s][m] = v;		
		}
	}
}

SpectrumScoreDistribution::~SpectrumScoreDistribution()
{
	deallocate();
}

bool SpectrumScoreDistribution::initialized() const
{
	return (table_ != NULL);
}

bool SpectrumScoreDistribution::failed() const
{
	return failed_;
}

FLOAT_T SpectrumScoreDistribution::maxScore() const
{
	return score_[nrows_-1];
}

FLOAT_T SpectrumScoreDistribution::maxMass() const
{
	return mass_[ncols_-1];
}

FLOAT_T SpectrumScoreDistribution::pvalue(FLOAT_T xcorr) const
{
	int m = ncols_ - 1;
	FLOAT_T nPeptides = 0;
	FLOAT_T nBetter = 0;
	for (int s = 0; s < nrows_; ++s) {
		nPeptides += (FLOAT_T) table_[s][m];
		if (score_[s] >= xcorr)
			nBetter += (FLOAT_T) table_[s][m];
	}
	return nBetter/nPeptides;
}

FLOAT_T SpectrumScoreDistribution::logpvalue(FLOAT_T xcorr) const
{
	int m = ncols_ - 1;
	FLOAT_T nPeptides = 0;
	FLOAT_T nBetter = 0;
	for (int s = 0; s < nrows_; ++s) {
		nPeptides += (FLOAT_T) table_[s][m];
		if (score_[s] >= xcorr)
			nBetter += (FLOAT_T) table_[s][m];
	}
	return log(nBetter) - log(nPeptides);
}
