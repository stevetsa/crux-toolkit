#ifndef QRANKER_CMD_H
#define QRANKER_CMD_H
/**
 * \file q-ranker.h
 */ 
/*
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 25, 2008
 * DESCRIPTION: Header file for crux q-ranker command. 
 */
#include "output-files.h"

MATCH_COLLECTION_T* run_qranker(
  char* psm_result_folder, 
  char* fasta_file, 
  OutputFiles& output);

#endif //QRANKER_CMD_H

