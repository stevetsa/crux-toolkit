#ifndef SPECTRAL_COUNTS_H
#define SPECTRAL_COUNTS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <math.h>
#include <sstream>
#include <algorithm>
#include "match_collection.h"
#include "utils.h"
#include "match.h"
#include "carp.h"
#include "peptide.h"
#include "protein.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "protein.h"
#include "database.h"
#include "objects.h"
#include "SpectrumCollection.h"
#include "Spectrum.h"
#include "FilteredSpectrumChargeIterator.h"

int spectral_counts_main(int argc, char** argv);

#endif
