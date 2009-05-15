#ifndef QRANKER_CMD_H
#define QRANKER_CMD_H
/**
 * \file match_analysis.c
 */
/*
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 25, 2008
 * DESCRIPTION: Header file for crux q-ranker command. 
 *
 * $Revision: 1.1.2.3 $
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "protein.h"
#include "peptide.h"
#include "spectrum.h"
#include "parse_arguments.h" 
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "scorer.h"
#include "match.h"
#include "match_collection.h"
#include "QRankerCInterface.h"


int qranker_main(int argc, char** argv);



#endif //QRANKER_CMD_H

