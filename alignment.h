/*
 *
 *    alignment.h
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <pthread.h>
#include "definitions.h"
#include "DV-Kernel.h"
#include "CPUfunctions.h"
#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "AlgnResult.h"
#include "DV-SemiDP.h"
#include "DV-DPForBothUnalign.h"
#include "DV-DPForSingleReads.h"
#include "OutputDPResult.h"
#include "IndexHandler.h"

// COPY INDEX TO DEVICE MEMORY
void GPUINDEXUpload ( Soap3Index * index,
                      uint ** _bwt, uint ** _occ,
                      uint ** _revBwt, uint ** _revOcc );

void GPUINDEXFree ( uint * _bwt, uint * _occ, uint * _revBwt, uint * _revOcc );

// perform round1 alignment in GPU
void perform_round1_alignment ( uint * nextQuery, uint * nextReadLength, uint * answers[][MAX_NUM_CASES],
                                uint numMismatch, uint numCases, uint sa_range_allowed, uint wordPerQuery, uint word_per_ans,
                                bool isExactNumMismatch, int doubleBufferIdx, uint blocksNeeded, ullint batchSize,
                                Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc );

// perform round2 alignment in GPU
void perform_round2_alignment ( uint * queries, uint * readLengths, uint * answers[][MAX_NUM_CASES],
                                uint numMismatch, uint numCases, uint sa_range_allowed_2, uint wordPerQuery, uint word_per_ans, uint word_per_ans_2,
                                bool isExactNumMismatch, int doubleBufferIdx, uint blocksNeeded, ullint batchSize,
                                Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                                uint processedQuery, uint * badReadIndices[][MAX_NUM_CASES],
                                uint * badAnswers[][MAX_NUM_CASES] );

// paired-end alignment: for all-valid, all-best, unique-best
// single-end alignment: for all-valid
void all_valid_alignment ( uint * queries, uint * readLengths, uint * seedLengths, uint numMismatch, uint wordPerQuery,
                           ullint maxBatchSize, uint numQueries, uint accumReadNum,
                           Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                           IniParams ini_params, InputOptions input_options,
                           char * upkdQualities,
                           uint * unAlignedReads, uint & numOfUnPaired,
                           uint * readIDs, char * upkdQueryNames,
                           char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                           unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                           uint8_t isTerminalCase,
                           ReadInputForDP ** readInputForDP,
                           ReadInputForDP ** readInputForNewDP,
                           ReadInputForDP ** otherSoap3Result,
                           BothUnalignedPairs ** bothUnalignedPairs );

// pair-end alignment: for random-best
// 4-phases [0,1,2,4]
void four_phases_alignment ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                             ullint maxBatchSize, uint numQueries, uint accumReadNum,
                             Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                             IniParams ini_params, InputOptions input_options,
                             char * upkdQualities,
                             uint * unAlignedPair, uint & numOfUnPaired,
                             uint * readIDs, char * upkdQueryNames,
                             char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                             unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                             ReadInputForDP ** readInputForDP,
                             ReadInputForDP ** readInputForNewDP,
                             ReadInputForDP ** otherSoap3Result,
                             BothUnalignedPairs ** bothUnalignedPairs );

// pair-end all-best alignment
// 3-phases [1,2,4]
void all_best_alignment ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                          ullint maxBatchSize, uint numQueries, uint accumReadNum,
                          Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                          IniParams ini_params, InputOptions input_options,
                          char * upkdQualities,
                          uint * unAlignedPair, uint & numOfUnPaired,
                          uint * readIDs, char * upkdQueryNames,
                          char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                          unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                          ReadInputForDPArrays * readInputForDPall,
                          ReadInputForDPArrays * readInputForNewDPall,
                          ReadInputForDPArrays * otherSoap3Resultall,
                          BothUnalignedPairsArrays * bothUnalignedPairsArrays );

// single-end alignment: for all-best, unique-best and random-best
void best_single_alignment ( uint * queries, uint * readLengths, uint * seedLengths, uint numMismatch, uint wordPerQuery,
                             ullint maxBatchSize, uint numQueries, uint accumReadNum,
                             Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                             IniParams ini_params, InputOptions input_options,
                             char * upkdQualities,
                             uint * unAlignedPair, uint & numOfUnPaired,
                             uint * readIDs, char * upkdQueryNames,
                             char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                             unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                             ReadInputForDP ** readInputForDP,
                             BothUnalignedPairs ** bothUnalignedPairs );

// Perform single-read all-best 1-mismatch alignment for seeds
// This function would reset the SingleAlgnResultArray
void single_1_mismatch_alignment2 ( unsigned int * queries, unsigned int * readLengths,
                                    unsigned int maxReadLength, unsigned int wordPerQuery, unsigned int maxBatchSize,
                                    unsigned int numQueries, Soap3Index * index, uint * _bwt, unsigned int * _revBwt,
                                    unsigned int * _occ, unsigned int * _revOcc,
                                    unsigned long long & numOfAnswer,
                                    unsigned int & numOfAlignedRead, unsigned int cpuNumThreads,
                                    SingleAlgnResultArray * alignResultArray, unsigned int maxHitNum );


// Perform SOAP3-DP Paired-End Alignment
void soap3_dp_pair_align ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                           ullint maxBatchSize, uint numQueries, uint accumReadNum,
                           Soap3Index * index,
                           uint * _bwt, uint * _revBwt,
                           uint * _occ, uint * _revOcc,
                           IniParams ini_params, InputOptions input_options,
                           uint maxReadLength, uint detected_read_length, uint detected_read_length2,
                           char * upkdQualities,
                           uint * unAlignedPair, uint & numOfUnPaired,
                           uint * readIDs, char * upkdQueryNames,
                           char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                           char * outputDPFileName, samfile_t * samOutputDPFilePtr,
                           samfile_t * samOutputUnpairFilePtr,
                           unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                           ReadInputForDPArrays * readInputForDPall,
                           ReadInputForDPArrays * readInputForNewDPall,
                           ReadInputForDPArrays * otherSoap3Resultall,
                           BothUnalignedPairsArrays * bothUnalignedPairsArrays,
                           double startTime, double & lastEventTime, double & totalAlignmentTime,
                           uint & indexLoadedToGPU );

// Perform SOAP3-DP Single Alignment
void soap3_dp_single_align ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                             ullint maxBatchSize, uint numQueries, uint accumReadNum,
                             Soap3Index * index,
                             uint * _bwt, uint * _revBwt,
                             uint * _occ, uint * _revOcc,
                             IniParams ini_params, InputOptions input_options,
                             uint maxReadLength, uint detected_read_length,
                             char * upkdQualities,
                             uint * unAlignedPair, uint & numOfUnPaired,
                             uint * readIDs, char * upkdQueryNames,
                             char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                             char * outputDPFileName, samfile_t * samOutputDPFilePtr,
                             unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                             ReadInputForDPArrays * readInputForDPall,
                             UnalignedSinglesArrays * unalignedSinglesArrays,
                             double startTime, double & lastEventTime, double & totalAlignmentTime,
                             uint & indexLoadedToGPU );

#endif
