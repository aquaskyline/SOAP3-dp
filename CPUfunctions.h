/*
 *
 *    CPUfunctions.h
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

#ifndef _CPUFUNCTIONS_H_
#define _CPUFUNCTIONS_H_

#include <pthread.h>
#include <math.h>
#include <ctype.h>
#include <sys/mman.h>
#include "2bwt-lib/dictionary.h"
#include "2bwt-lib/iniparser.h"
#include "2bwt-lib/Timing.h"
#include "BGS-HostAlgnmtAlgo2.h"
#include "BGS-HostAlgnmtAlgoSingle.h"
#include "definitions.h"
#include "SAList.h"
#include "PE.h"
#include "PEAlgnmt.h"
#include "AlgnResult.h"
#include "global_arrays.h"
#include "SAM.h"
#include "IniParam.h"
#include "UsageInterface.h"
#include "IndexHandler.h"

// The offset mask to retrieve the least significant 24-bit from the 32-bit word.
#define  BGS_GPU_ANSWER_OFFSET_LENGTH   24
#define  BGS_GPU_ANSWER_OFFSET_MASK     ((1<<BGS_GPU_ANSWER_OFFSET_LENGTH)-1)


typedef struct HostKernelArguements
{
    char * upkdQualities;
    char * upkdQueryNames;
    unsigned int batchFirstReadId;
    unsigned int skipFirst;
    unsigned int numQueries;

    OCC * occ;

    SRAIndex sraIndex;
    SRASetting sraQuerySettings;
    SRAQueryResultCount sraAccumulatedResultCount;
    SRAModel ** SRAMismatchModel;
    SRAModel ** SRAMismatchModel_neg;
    unsigned int numCases;

    // for computing the second best alignments
    // they are only used when:
    // (1) alignment type = all valid, AND output = SAM format; OR
    // (2) alignment type = the best and the second best
    SRAModel ** SRAMismatchModel2;
    SRAModel ** SRAMismatchModel2_neg;
    unsigned int numCases2;

    char * outputFileName;
    unsigned int * queries;
    unsigned int word_per_query;
    unsigned int ** answers;
    unsigned int ** badReadIndices;
    unsigned int ** badAnswers;
    unsigned int * badStartOffset;
    unsigned int * badCountOffset;
    unsigned int maxReadLength;
    unsigned int * readLengths;
    unsigned int * seedLengths;
    unsigned int * readIDs;
    unsigned int * unAlignedIDs;
    unsigned int unAlignedOcc;

    unsigned long long int alignedOcc;
    unsigned int alignedReads;
    int threadId;

    char outputGoodReads;
    unsigned int sa_range_allowed_1;
    unsigned int sa_range_allowed_2;
    char skip_round_2;
    unsigned int accumReadNum;
    unsigned int word_per_ans;

    char readType; // 1: single read; 2: pair-end read
    char maxNumMismatch;
    int insert_low;
    int insert_high;
    int peStrandLeftLeg;
    int peStrandRightLeg;
    char reportType; // 1: all valid; 2: all best; 3: unique best; 4: random best

    unsigned int peMaxOutputPerPair;
    uint8_t isTerminalCase;
    ReadInputForDP * readInput; // the alignment result of the reads which have hit but their mates have not.
    // these results are needed for proceeding to semi-global DP
    ReadInputForDP * readInputForNewDefault; // the alignment result of the reads which have hit but their mates have not.
    // When the number of alignment result is too many (i.e.  > 30)
    // these results are needed for proceeding to new semi-global DP
    ReadInputForDP * otherSoap3Result; // the alignment result of the unpaired reads whose both ends have hits
    BothUnalignedPairs * bothUnalignedPairs; // it stores the first read IDs of the pairs of which both ends cannot be aligned
    unsigned int maxHitNum; // maximum number of hits allowed for proceeding DP for even read of paired-end read
    unsigned int maxHitNum2; // maximum number of hits allowed for proceeding DP for odd read of paired-end read
    char enableDP; // 1: enable; 0: disable dynamic programming for the unmapped mate

    SingleAlgnResult * alignResult; // alignment results of single-end reads

} HostKernelArguements;

// for user does not specify # of mismatches and dp is disabled
// if the read length < 50, then return DEFAULT_NUM_MISMATCH_FOR_SHORT_READ
// else return DEFAULT_NUM_MISMATCH_FOR_NORMAL_READ
int getDefaultMismatchNum ( uint read_length );

// get the max hit # for default DP
int getMaxHitNumForDefaultDP ( uint read_length );

// get the parameters for All DP
void getParameterForAllDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params );

// get the parameters for Default DP
void getParameterForNewDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 );

// get the parameters for Default DP
void getParameterForDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 );

// get the parameters for Deep DP
void getParameterForDeepDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 );

// get the parameters for Single-end DP
void getParameterForSingleDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length );

// A Thread Wrapper to perform the CPU task
void * hostKernelThreadWrapper ( void * arg );

// A Thread Wrapper to perform the CPU task
// for single-end alignment and the results are stored in the arrays
void * hostKernelThreadWrapperSingle ( void * arg );

// set the multi-threading arguments
void setHostKernelArguments ( HostKernelArguements * hostKernelArguments, pthread_t * threads, IniParams ini_params, Soap3Index * index,
                              uint maxReadLength, uint word_per_query, uint word_per_ans,
                              InputOptions * input_options );
//char readType, uint insert_high, uint insert_low);

// set the multi-threading arguments
// for 1-mismatch single-end alignment
void setHostKernelArguments2 ( HostKernelArguements * hostKernelArguments, pthread_t * threads, Soap3Index * index,
                               uint maxReadLength, uint word_per_query, uint word_per_ans,
                               uint cpuNumThreads );

// obtain the number of cases for the number of mismatch
int GetNumCases ( uint numMismatch );

// obtain the number of cases for this number of mismatch
// and obtain the number of SA ranges allowed
void getParametersForThisMismatch ( uint numMismatch, uint & numCases, uint & sa_range_allowed_1,
                                    uint & sa_range_allowed_2, char & skip_round_2,
                                    uint & word_per_ans, uint & word_per_ans_2 );

// obtain the number of cases for this number of mismatch
// and obtain the number of SA ranges allowed
// this is designed for seed alignments
void getParametersForThisMismatch2 ( uint numMismatch, uint & numCases, uint & sa_range_allowed_1,
                                     uint & sa_range_allowed_2, char & skip_round_2,
                                     uint & word_per_ans, uint & word_per_ans_2 );

// pack the reads which are unpaired
void packUnPairedReads ( uint * queries, uint * readIDs, uint * readLengths, uint * unAlignedPair,
                         uint wordPerQuery, uint numOfUnPaired, ullint maxBatchSize );


// pack the reads with no alignment together
// return number of reads with no alignment
uint packReads ( uint * queries, uint * readIDs, uint * readLengths, uint * seedLengths, unsigned char * noAlignment,
                 uint wordPerQuery, uint numQueries );

// pack the reads with no alignment together
// return readIDS of the unaligned reads
uint * packReads2 ( uint * queries, uint * readLengths, unsigned char * noAlignment,
                    uint wordPerQuery, uint numQueries, uint & numUnAligned );


// repack the reads
// no read will be removed, but
// the reads which need to be processed in next-round by soap3 will be duplicated
// to the front of the list. The readIDs are stored inside the array called "needProcessPair"
// the corresponding readIDs inside "readInputForDP", "readInputForNewDP" and
// "bothUnalignedPairs" need to be updated correspondingly.
void repackUnPairedReads ( uint ** queries, uint ** readIDs, uint ** readLengths, uint * needProcessPair,
                           uint wordPerQuery, uint numOfReadsToProcess, ullint numOfTotalReads,
                           ReadInputForDPArrays * readInputForDPall,
                           ReadInputForDPArrays * readInputForNewDPall,
                           BothUnalignedPairsArrays * bothUnalignedPairsArrays );

// show the merging file command at the end of the program
void show_merge_file_command ( InputOptions input_options, char * queryFileName, char * queryFileName2 );

void show_merge_file_command2 ( InputOptions input_options, char * outputFileName, int cpuNumThreads );

// convert 2-bit-per-character array to 1-byte-per-character array
void retrieve2BWTQueryFormat ( unsigned int * packedPatterns,
                               unsigned int queryIdx, unsigned int wordPerQuery, unsigned char * unpackedPattern );

// For initialization of g_log_n array. This array will be used for calculation of MAPQ of
void bwase_initialize ( int * g_log_n );

#endif
