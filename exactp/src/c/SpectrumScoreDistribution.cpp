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
#include "Alphabet.h"

#include "SpectrumScoreDistribution.h"
#include <limits.h>
#include "mass.h"
 #include <algorithm>

#define ROUND(x) (floor(0.5+x))

const double SpectrumScoreDistribution::epsilon_ = numeric_limits<double>::epsilon();

template <class Iter, class T>
Iter closest(Iter begin, Iter end, T const& value)
{
	Iter result = std::upper_bound(begin, end, value);
	if (result != begin)
	{
		Iter lower_result = result;
		--lower_result;
		if (result == end || (value - *lower_result) < (*result - value))
		{
			result = lower_result;
		}
	}
	return result;
}

void SpectrumScoreDistribution::quantize(
	FLOAT_T lower,
	FLOAT_T upper,
	FLOAT_T delta,
	FLOAT_T** ovalues,
	int& olen
	)
{
	assert(lower < upper);
	assert(delta > 0.0);

	olen = ceil((upper-lower)/delta);
	*ovalues = (FLOAT_T*)mycalloc(olen, sizeof(FLOAT_T));
	for (int i = 0; i < olen; ++i)
		(*ovalues)[i] = i * delta;
}


SpectrumScoreDistribution::SpectrumScoreDistribution(
	Spectrum* spectrum,
	const SpectrumZState& zstate
	) : table_(NULL), nrows_(-1), ncols_(-1), lower_(-1), upper_(-1), 
			mass_(NULL), score_(NULL), observed_(NULL), maxbin_(-1), failed_(false), 
			neutralmass_(0.0), aamass_(NULL), naa_(-1)
{
	if ( spectrum == NULL )
		carp(CARP_FATAL, "Cannot compute score distribution, spectrum is NULL.");

	// Initialize the table of amino acids masses
	MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");
	naa_ = Alphabet::numAminoAcids;
	aamass_ = (int*)mycalloc(naa_, sizeof(int));
	for (int a = 0; a < Alphabet::numAminoAcids; ++a) {
		char residue = *(Alphabet::aminoAcids[a]);
		aamass_[a] = (int) ROUND(get_mass_amino_acid(residue, mass_type));
	}

	// Quantize the mass range.
	int charge = zstate.getCharge();
	double tol = get_double_parameter("precursor-window");
	FLOAT_T maxMass = ceil(zstate.getNeutralMass() + tol);
	neutralmass_ = zstate.getNeutralMass();
	FLOAT_T massDelta = get_double_parameter("mass-delta");  // 1.0
 	quantize(0.0, maxMass, massDelta, &mass_, ncols_);
 	assert(mass_);

#if 0
	// Determine the columns to sum up to form the null distribution, since we do
	// not assume that each column corresponds to 1.0Da.
	int l = max(0, (int) ceil(neutralmass_ - tol));
	//int h = min((int) floor(neutralmass_ + tol), ncols_-1);
	upper_ = ncols_ - 1;
	for (int m = 0; m < ncols_; ++m) {
		if (mass_[m] <= l) {
			lower_ = m;
			break;
		}
	}
#endif
	int l = max(0, (int)ROUND(neutralmass_));
	for (int m = 0; m < ncols_; ++m) {
		if (mass_[m] <= l) {
			lower_ = m;
			upper_ = m;
			break;
		}
	}
	assert(lower_ >= 0);
	assert(upper_ < ncols_);
	assert(lower_ <= upper_);

	FLOAT_T maxScore = get_double_parameter("max-xcorr"); // 10.0
	FLOAT_T scoreDelta = get_double_parameter("score-delta"); // 0.01
	//carp(CARP_INFO, "max-xcorr = %g, score-delta = %g", maxScore, scoreDelta);
	quantize(0.0, maxScore, scoreDelta, &score_, nrows_);
	assert(score_);	

	// Allocate memoization table.
	table_ = (double**)mycalloc(nrows_, sizeof(double*));
	for (int i = 0; i < nrows_; ++i)
		table_[i] = (double*)mycalloc(ncols_, sizeof(double));


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

  for (int b = 0; b < maxbin_; ++b)
  	assert(observed_[b] > -epsilon_);

  assert(nrows_ > 0);
  assert(ncols_ > 0);
  assert(mass_);
  assert(score_);
  assert(observed_);
  assert(maxbin_ > 0);
  //assert(maxbin_ <= ncols_); // Max mass in DP table > mass of spectrum. // NO, as teh bins are overalloc.

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

int SpectrumScoreDistribution::mass_to_index(FLOAT_T mass) const
{
	const FLOAT_T* ptr = closest(mass_, mass_ + ncols_, mass);
	return (int)(ptr - mass_);
}

int SpectrumScoreDistribution::score_to_index(FLOAT_T score) const
{
	const FLOAT_T* ptr = closest(score_, score_ + nrows_, score);
	return (int)(ptr - score_);
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
	//static const int mass[20] = { 71, 160, 129, 115, 57, 147, 113, 137, 128, 131, 113, 
	//															114, 128, 97, 87, 156, 101, 186, 99, 163 };
	//static double epsilon = numeric_limits<double>::epsilon();

	// The only non-zero value in table_[0][*] is the first one.
	int min_mass = *min_element(aamass_, aamass_ + naa_);
	int num_amino_acids = naa_;
	table_[0][0] = 1;
	for (int m = 0; m < ncols_; ++m) {
		for (int s = 0; s < nrows_; ++s) {
			if (s == 0) {
				continue;
			}
			// Recurrence
			double v = 0;
			for (int a = 0; a < num_amino_acids; ++a) {
				int new_mass = mass_[m] - aamass_[a];
				if (new_mass + min_mass < 0) {
					// new_mass is negative, so much so it cannot be explained by rounding
					// errors. No contribution to v, as this recursive case is impossible.
					continue;
				}
				new_mass = max(0, new_mass);
				//fprintf(stderr, "new_mass(m = %d, s = %d, a = %d) = %d\n", m, s, a, new_mass);

				// Get the contribution of the spectrum at new_mass to the score.
				// Check that the spectrum has been scaled correctly, i.e., that the
				// contribution is non-negative.
				assert(new_mass >= 0 && new_mass < maxbin_);
				FLOAT_T contrib = (50.0 * observed_[new_mass])/10000.0;
				assert(contrib > -epsilon_);
				//fprintf(stderr, "contrib = %g\n", contrib);

				FLOAT_T new_score = score_[s] - contrib;
				if (ROUND(new_score) < epsilon_ && new_mass < epsilon_) {
					// The score rounds to zero and the mass is close to zero
					v += table_[0][0];
					continue;
				}	else if (new_score < epsilon_) {
					// Score is too low for the recursive case to be possible.
					continue;
				}
				//fprintf(stderr, "new_score(m = %d, s = %d, a = %d) = %g\n", m, s, a, new_score);

				// Add in the recursive cases.
				int s_idx = score_to_index(new_score);
				int m_idx = mass_to_index(new_mass);
				assert(s_idx >= 0 && s_idx < nrows_);
				assert(m_idx >= 0 && m_idx < ncols_);
				//fprintf(stderr, "A[score = %g, idx = %d][mass = %d, idx = %d] = %g\n",
				//	new_score, s_idx, new_mass, m_idx, table_[s_idx][m_idx]);
				assert(table_[s_idx][m_idx] > -epsilon_);
				v += table_[s_idx][m_idx];				
			}
			//fprintf(stderr, "v(score=%g,mass=%g) = %g\n", score_[s],mass_[m],v);
			assert(v >= -epsilon_);
			table_[s][m] = v;		
		}
	}
}

/**
 * To compute the distribution over peptide scores, sum up the columns which
 * represent masses in the candidate peptide mass interval: [nm - tol, nm + tol],
 * where nm is the neutral mass of the precursor ion, observed during MS1.
 */
void SpectrumScoreDistribution::countHigherScoring(
	FLOAT_T xcorr,
	double& nBetter,
	double& nPeptides
	) const
{
	if (xcorr >= maxScore())
		carp(CARP_WARNING, "Cannot accurately score xcorr %g (max %g)", xcorr, maxScore());

	nPeptides = 0;
	nBetter = 0;
	for (int m = lower_; m <= upper_; ++m) {
		for (int s = 0; s < nrows_; ++s) {
			nPeptides += (FLOAT_T) table_[s][m];
			if (score_[s] >= xcorr)
				nBetter += (FLOAT_T) table_[s][m];
		}
	}

#if 0
	nPeptides = 0;
	nBetter = 0;
	//int m = ncols_ - 1;
	int m = (int) floor(neutralmass_);
	if (m >= ncols_)
		m -= 1;
	if (m >= ncols_)
		carp(CARP_FATAL, "Bad neutral mass");

	for (int s = 0; s < nrows_; ++s) {
		nPeptides += (FLOAT_T) table_[s][m];
		if (score_[s] >= xcorr)
			nBetter += (FLOAT_T) table_[s][m];
	}
#endif
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
	double nBetter, nPeptides;
	countHigherScoring(xcorr, nBetter, nPeptides);
	if (fabs(nBetter) < epsilon_) {
		return 10*epsilon_;
	}
	if (nPeptides < -epsilon_) {
		carp(CARP_WARNING, "Overflow, return worst p-value");
		return 1.0;
	}
	return nBetter/nPeptides;
}

FLOAT_T SpectrumScoreDistribution::logpvalue(FLOAT_T xcorr) const
{
	double nBetter, nPeptides;
	countHigherScoring(xcorr, nBetter, nPeptides);
	if (fabs(nBetter) < epsilon_) {
		return 10*epsilon_;
	}
	if (nPeptides < -epsilon_) {
		carp(CARP_WARNING, "Overflow, return worst p-value");
		return 0.0;
	}	
	return log(nBetter) - log(nPeptides);
}
