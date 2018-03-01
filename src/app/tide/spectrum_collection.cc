// Benjamin Diament
//
// See .h file.

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include "spectrum.pb.h"
#include "spectrum_collection.h"
#include "mass_constants.h"
#include "max_mz.h"
#include "records.h"
#include "records_to_vector-inl.h"
#include "util/mass.h"
#include "util/Params.h"
#if defined ( _MSC_VER ) || defined ( DARWIN )
#include <unordered_set>
#else
#include <tr1/unordered_set>
#endif

using namespace std;
using google::protobuf::uint64;

#define CHECK(x) GOOGLE_CHECK((x))

// Integerization constant for the XCorr p-value calculation.
#define EVIDENCE_INT_SCALE 500.0

Spectrum::Spectrum(const pb::Spectrum& spec) {
  spectrum_number_ = spec.spectrum_number();
  precursor_m_z_ = spec.precursor_m_z();
  rtime_ = spec.rtime();
  for (int i = 0; i < spec.charge_state_size(); ++i)
    charge_states_.push_back(spec.charge_state(i));
  int size = spec.peak_m_z_size();
  CHECK(size == spec.peak_intensity_size());
  uint64 total = 0;
  double m_z_denom = spec.peak_m_z_denominator();
  double intensity_denom = spec.peak_intensity_denominator();
  peaks_ = vector< pair<double, double> >(size);
  for (int i = 0; i < size; ++i) {
    CHECK(spec.peak_m_z(i) > 0);
    total += spec.peak_m_z(i); // deltas of m/z are stored
    peaks_[i] = make_pair(total / m_z_denom, spec.peak_intensity(i) / intensity_denom);
  }
  sorted_ = true;
}

// A spectrum can have multiple precursor charges assigned.  This
// reports the maximum such charge state.
int Spectrum::MaxCharge() const {
  vector<int>::const_iterator i = max_element(charge_states_.begin(), charge_states_.end());
  return i != charge_states_.end() ? *i : 1;
}

// Report maximum intensity peak in the given m/z range.
double Spectrum::MaxPeakInRange( double min_range, double max_range ) const {
  double maxIntensity = 0;
  for (vector< pair<double, double> >::const_iterator i = peaks_.begin(); i != peaks_.end(); i++) {
    double mz = i->first;
    if (min_range <= mz && mz <= max_range) {
      double intensity = i->second;
      if (intensity > maxIntensity) {
        maxIntensity = intensity;
      }
    }
  }
  return maxIntensity;
}

static inline bool IsInt(double x) {
  // See whether x is quite close to an integer (within 0.001).
  return fabs(x - uint64(x+0.5)) < 0.001;
}

static inline bool CheckDenom(const vector<double>& vals, int denom) {
  // See whether all vals can be accommodated by denom when rendered as a
  // fraction.
  double d_denom = denom;
  for (int i = 0; i < vals.size(); ++ i) {
    if (!IsInt(vals[i] * d_denom))
      return false;
  }
  return true;
}

static inline int GetDenom(const vector<double>& vals) {
  // See how much precision is given in the vals array. Not especially fast,
  // but fast enough. Used only for converting spectrum input format.
  const int kMaxPrecision = 10000; // store at most 3 digits of precision
  for (int precision = 1; precision < kMaxPrecision; precision *= 10)
    if (CheckDenom(vals, precision))
      return precision;
  return kMaxPrecision;
}

// Do Morpheus-style simple(-istic?) deisotoping.  "For each
// peak, lower m/z peaks are considered. If the reference peak
// lies where an expected peak would lie for a charge state from
// one to the charge state of the precursor, within mass
// tolerance, and is of lower abundance, the reference peak is
// considered to be an isotopic peak and removed."
bool Spectrum::Deisotope(int index, double deisotope_threshold) const {
  if (deisotope_threshold == 0.0) {
    return false;
  }
  double location, intensity;
  GetPeak(index, &location, &intensity);
  int maxCharge = MaxCharge();
  for (int fragCharge = 1; fragCharge < maxCharge; fragCharge++) {
    double isotopic_peak = location - (ISOTOPE_SPACING / fragCharge);
    double ppm_difference = (location * deisotope_threshold) / 1e6;
    double isotopic_intensity = MaxPeakInRange(isotopic_peak - ppm_difference,
                                               isotopic_peak + ppm_difference);

    if (intensity < isotopic_intensity) {
      carp(CARP_DETAILED_DEBUG,
           "Removing isotopic peak (%g, %g) because of peak in [%g, %g] with intensity %g.",
           location, intensity, isotopic_peak - ppm_difference,
           isotopic_peak + ppm_difference, isotopic_intensity);
      return true;
    }
  }
  return false;
}

/* Calculates vector of cleavage evidence for an observed spectrum, using XCorr
 * b/y/neutral peak sets and heights.
 *
 * Written by Jeff Howbert, May, 2013 (as function createEvidenceArrayObserved).
 * Extended and modified by Jeff Howbert, October, 2013.
 * Ported to and integrated with Tide by Jeff Howbert, November, 2013.
 */
vector<double> Spectrum::CreateEvidenceVector(
  double binWidth,
  double binOffset,
  int charge,
  double pepMassMonoMean,
  int maxPrecurMass,
  long int* num_range_skipped,
  long int* num_precursors_skipped,
  long int* num_isotopes_skipped,
  long int* num_retained
) const {
  // TODO need to review these constants, decide which can be moved to parameter file
  const double maxIntensPerRegion = 50.0;
  const double BYHeight = 50.0;
  const double NH3LossHeight = 10.0;
  const double COLossHeight = 10.0;    // for creating a ions on the fly from b ions
  const double H2OLossHeight = 10.0;
  const double FlankingHeight = BYHeight / 2;;
  // TODO end need to review
  int numPeaks = Size();
  double experimentalMassCutoff = PrecursorMZ() * charge + 50.0;
  double maxIonMass = 0.0;
  double maxIonIntens = 0.0;

  if (!sorted_) {
    carp(CARP_FATAL, "spectrum unsorted error");
  }

  // Find max ion mass and max ion intensity
  bool skipPreprocess = Params::GetBool("skip-preprocessing");
  bool remove_precursor = !skipPreprocess && Params::GetBool("remove-precursor-peak");
  double precursorMZExclude = Params::GetDouble("remove-precursor-tolerance");
  double deisotope_threshold = Params::GetDouble("deisotope");
  tr1::unordered_set<int> peakSkip;
  for (int ion = 0; ion < numPeaks; ion++) {
    double ionMass, ionIntens;
    GetPeak(ion, &ionMass, &ionIntens);
    if (ionMass >= experimentalMassCutoff) {
      peakSkip.insert(ion);
      if (num_range_skipped) {
        (*num_range_skipped)++;
      }
      continue;
    } else if (remove_precursor && ionMass > PrecursorMZ() - precursorMZExclude &&
               ionMass < PrecursorMZ() + precursorMZExclude) {
      peakSkip.insert(ion);
      if (num_precursors_skipped) {
        (*num_precursors_skipped)++;
      }
      continue;
    } else if (deisotope_threshold != 0.0 && Deisotope(ion, deisotope_threshold)) {
      peakSkip.insert(ion);
      if (num_isotopes_skipped) {
        (*num_isotopes_skipped)++;
      }
      continue;
    }
    if (num_retained) {
      (*num_retained)++;
    }
    if (maxIonMass < ionMass) {
      maxIonMass = ionMass;
    }
    if (maxIonIntens < ionIntens) {
      maxIonIntens = ionIntens;
    }
  }

  // 10 bin intensity normalization 
  double regionSelector = (int)floor(MassConstants::mass2bin(maxIonMass) / (double)NUM_SPECTRUM_REGIONS);
  vector<double> intensObs(maxPrecurMass, 0);
  vector<int> intensRegion(maxPrecurMass, -1);
  for (int ion = 0; ion < numPeaks; ion++) {
    if (peakSkip.find(ion) != peakSkip.end()) {
      continue;
    }
    double ionMass, ionIntens;
    GetPeak(ion, &ionMass, &ionIntens);
    int ionBin = MassConstants::mass2bin(ionMass);
    int region = (int)floor((double)(ionBin) / regionSelector);
    if (region >= NUM_SPECTRUM_REGIONS) {
      region = NUM_SPECTRUM_REGIONS - 1;
    }
    intensRegion[ionBin] = region;
    if (intensObs[ionBin] < ionIntens) {
      intensObs[ionBin] = ionIntens;
    }
  }

  const double intensThreshold = 0.05 * sqrt(maxIonIntens);
  for (vector<double>::iterator i = intensObs.begin(); i != intensObs.end(); i++) {
    double newIntens = sqrt(*i);
    if (newIntens <= intensThreshold) {
      newIntens = 0.0;
    }
    *i = newIntens;
  }

  vector<double> maxRegion(NUM_SPECTRUM_REGIONS, 0);
  for (int i = 0; i < maxPrecurMass; i++) {
    int reg = intensRegion[i];
    double ionIntens = intensObs[i];
    if (reg >= 0 && ionIntens > maxRegion[reg]) {
      maxRegion[reg] = ionIntens;
    }
  }

  vector<double> partial_sums(maxPrecurMass, 0);
  double total = 0.0;
  for (int i = 0; i < maxPrecurMass; i++) {
    int reg = intensRegion[i];
    if (reg >= 0 && maxRegion[reg] > 0.0) {
      intensObs[i] *= (maxIntensPerRegion / maxRegion[reg]);
    }
    partial_sums[i] = (total += intensObs[i]);
  }

  // ***** Adapted from tide/spectrum_preprocess2.cc.
  // TODO replace, if possible, with call to
  // static void SubtractBackground(double* observed, int end).
  // Note numerous small changes from Tide code.
  const double multiplier = 1.0 / (MAX_XCORR_OFFSET * 2.0 + 1.0);
  for (int i = 0; i < maxPrecurMass; ++i) {
    int right = i + MAX_XCORR_OFFSET;
    if (right >= maxPrecurMass) {
      right = maxPrecurMass - 1;
    }
    int left = i - MAX_XCORR_OFFSET - 1;
    if (left < 0) {
      left = 0;
    }
    intensObs[i] -= multiplier * (partial_sums[right] - partial_sums[left]);
  }

  bool flankingPeaks = Params::GetBool("use-flanking-peaks");
  bool nlPeaks = Params::GetBool("use-neutral-loss-peaks");
  int binFirst = MassConstants::mass2bin(30);
  int binLast = MassConstants::mass2bin(pepMassMonoMean - 47);
  vector<double> evidence(maxPrecurMass, 0);
  for (int i = binFirst; i <= binLast; i++) {
    // b ion
    double bIonMass = (i - 0.5 + binOffset) * binWidth;
    double current = intensObs[MassConstants::mass2bin(bIonMass)] * BYHeight;
    for (int j = 2; j < charge; j++) {
      current += intensObs[MassConstants::mass2bin(bIonMass, j)] * BYHeight;
    }
    // y ion
    double yIonMass = pepMassMonoMean + 2 * MASS_H_MONO - bIonMass;
    current += intensObs[MassConstants::mass2bin(yIonMass)] * BYHeight;
    for (int j = 2; j < charge; j++) {
      current += intensObs[MassConstants::mass2bin(yIonMass, j)] * BYHeight;
    }
    if (flankingPeaks) {
      // flanking peaks for b ions
      int ionBin = MassConstants::mass2bin(bIonMass, 1);
      current += intensObs[ionBin + 1] * FlankingHeight;
      current += intensObs[ionBin - 1] * FlankingHeight;
      for (int j = 2; j < charge; j++) {
        current += intensObs[MassConstants::mass2bin(bIonMass, j) + 1] * FlankingHeight;
        current += intensObs[MassConstants::mass2bin(bIonMass, j) - 1] * FlankingHeight;
      }
      // flanking peaks for y ions
      ionBin = MassConstants::mass2bin(yIonMass, charge);
      current += intensObs[ionBin + 1] * FlankingHeight;
      current += intensObs[ionBin - 1] * FlankingHeight;
      for (int j = 2; j < charge; j++) {
        current += intensObs[MassConstants::mass2bin(yIonMass, j) + 1] * FlankingHeight;
        current += intensObs[MassConstants::mass2bin(yIonMass, j) - 1] * FlankingHeight;
      }
    }
    if (nlPeaks) {
      // NH3 loss from b ion
      double ionMassNH3Loss = bIonMass - MASS_NH3_MONO;
      current += intensObs[MassConstants::mass2bin(ionMassNH3Loss)] * NH3LossHeight;
      for (int j = 2; j < charge; j++) {
        current += intensObs[MassConstants::mass2bin(ionMassNH3Loss, j)] * NH3LossHeight;
      }
      // NH3 loss from y ion
      ionMassNH3Loss = yIonMass - MASS_NH3_MONO;
      current += intensObs[MassConstants::mass2bin(ionMassNH3Loss)] * NH3LossHeight;
      for (int j = 2; j < charge; j++) {
        current += intensObs[MassConstants::mass2bin(ionMassNH3Loss, j)] * NH3LossHeight;
      }
      // CO and H2O loss from b ion
      double ionMassCOLoss = bIonMass - MASS_CO_MONO;
      double ionMassH2OLoss = bIonMass - MASS_H2O_MONO;
      current += intensObs[MassConstants::mass2bin(ionMassCOLoss)] * COLossHeight;
      current += intensObs[MassConstants::mass2bin(ionMassH2OLoss)] * H2OLossHeight;
      for (int j = 2; j < charge; j++) {
        current += intensObs[MassConstants::mass2bin(ionMassCOLoss, j)] * COLossHeight;
        current += intensObs[MassConstants::mass2bin(ionMassH2OLoss, j)] * H2OLossHeight;
      }
      // H2O loss from y ion
      ionMassH2OLoss = yIonMass - MASS_H2O_MONO;
      current += intensObs[MassConstants::mass2bin(ionMassH2OLoss)] * H2OLossHeight;
      for (int j = 2; j < charge; j++) {
        current += intensObs[MassConstants::mass2bin(ionMassH2OLoss, j)] * H2OLossHeight;
      }
    }
    evidence[i] = current;
  }
  return evidence;
}

vector<int> Spectrum::DiscretizeEvidenceVector(const std::vector<double>& evidence) {
  vector<int> discretized;
  discretized.reserve(evidence.size());
  for (vector<double>::const_iterator i = evidence.begin(); i != evidence.end(); i++) {
    discretized.push_back((int)floor(*i / EVIDENCE_INT_SCALE + 0.5));
  }
  return discretized;
}

bool SpectrumCollection::ReadSpectrumRecords(const string& filename,
                                             pb::Header* header) {
  pb::Header tmp_header;
  if (header == NULL)
    header = &tmp_header;
  HeadedRecordReader reader(filename, header);
  if (header->file_type() != pb::Header::SPECTRA)
    return false;
  pb::Spectrum pb_spectrum;
  while (!reader.Done()) {
    reader.Read(&pb_spectrum);
    spectra_.push_back(new Spectrum(pb_spectrum));
  }
  if (!reader.OK()) {
    for (int i = 0; i < spectra_.size(); ++i)
      delete spectra_[i];
    spectra_.clear();
    return false;
  }
  return true;
}

void SpectrumCollection::MakeSpecCharges() {
  // Create one entry in the spec_charges_ array for each
  // (spectrum, charge) pair.
  int spectrum_index = 0;
  vector<Spectrum*>::iterator i = spectra_.begin();
  for (; i != spectra_.end(); ++i) {
    for (int j = 0; j < (*i)->NumChargeStates(); ++j) {
      int charge = (*i)->ChargeState(j);
      double neutral_mass = (((*i)->PrecursorMZ() - MASS_PROTON)
                             * charge);
      spec_charges_.push_back(SpecCharge(neutral_mass, charge, *i,
                                         spectrum_index));
    }
    spectrum_index++;
  }
}

double SpectrumCollection::FindHighestMZ() const {
  // Return the maximum MZ seen across all input spectra.
  double highest = 0;
  vector<Spectrum*>::const_iterator i = spectra_.begin();
  for (; i != spectra_.end(); ++i) {
    CHECK((*i)->Size() > 0) << "ERROR: spectrum " << (*i)->SpectrumNumber()
                            << " has no peaks.\n";
    double last_peak = (*i)->M_Z((*i)->Size() - 1);
    if (last_peak > highest)
      highest = last_peak;
  }
  return highest;
}

void SpectrumCollection::Sort() {
  MakeSpecCharges();
  sort(spec_charges_.begin(), spec_charges_.end());
}
