/*
 *
 *    definitions.h
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

#ifndef __DEFINITIONS_H__
#define __DEFINITIONS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "2bwt-flex/OCC.h"
#include "2bwt-flex/SRAArguments.h"
#include "PEAlgnmt.h"
#include "AlgnResult.h"

#define BGS_DEVICE_WARP_SIZE 32

#define MAX_READ_LENGTH 1024
#define MAX_CHROMOSOME_NUM 65536
#define MAX_FIELD_LEN 256 // the maximum length of each field in the input file for pair-multi mode

// Note that when the round 1's SA Range is updated,
// the value of "WORD_PER_ANSWER" needs to update too.
// 1st round constants
// #define WORD_PER_ANSWER 4 // for all types of mismatch
// for exact-match
#define MAX_SA_RANGES_ALLOWED1_0 2
// for 1-mismatch
#define MAX_SA_RANGES_ALLOWED1_1 4
// for 2-mismatch
#define MAX_SA_RANGES_ALLOWED1_2 4
// for 3-mismatch
#define MAX_SA_RANGES_ALLOWED1_3 2
// for 4-mismatch
#define MAX_SA_RANGES_ALLOWED1_4 1


// 2nd round constants
// changing the following values may affect the memory control of the program
// for exact-match
#define MAX_SA_RANGES_ALLOWED2_0 16
// for 1-mismatch
#define MAX_SA_RANGES_ALLOWED2_1 512
// for 2-mismatch
#define MAX_SA_RANGES_ALLOWED2_2 32
// for 3-mismatch
#define MAX_SA_RANGES_ALLOWED2_3 16
// for 4-mismatch
#define MAX_SA_RANGES_ALLOWED2_4 16
// for seed alignment (no scheduling)
#define MAX_SA_RANGES_ALLOWED1_SMALL_READS 4
#define MAX_SA_RANGES_ALLOWED2_SMALL_READS 1024

// best 8192 128
#define NUM_BLOCKS 8192
#define THREADS_PER_BLOCK 128
#define QUERIES_PER_THREAD 1
#define INPUT_BUFFER_SIZE 102400
#define MAX_NUM_BATCH 12
// # of reads will be loaded into the memory each time = MAX_NUM_BATCH*NUM_BLOCKS*THREADS_PER_BLOCK

// whether the second round is skipped
// for exact-match
#define SKIP_ROUND2_0 1
// for 1-mismatch
#define SKIP_ROUND2_1 0
// for 2-mismatch
#define SKIP_ROUND2_2 0
// for 3-mismatch
#define SKIP_ROUND2_3 1
// for 4-mismatch
#define SKIP_ROUND2_4 1

#define GPU_OCC_INTERVAL 128
#define GPU_OCC_VALUE_PER_WORD 1

// partition for 2-mismatch
#define SIZE_X_RATIO .3
#define SIZE_Y_RATIO .3
// SIZE_Z = the remaining

// partition for 3-mismatch
#define SIZE_1_RATIO .25
#define SIZE_2_RATIO .25
#define SIZE_3_RATIO .25
// SIZE_4 = the remaining

// partition for 4-mismatch
#define SIZE_A_RATIO .2
#define SIZE_B_RATIO .2
#define SIZE_C_RATIO .2
#define SIZE_D_RATIO .2
// SIZE_E = the remaining


#define NUM_CASES_0M 1
#define NUM_CASES_1M 2
#define NUM_CASES_2M 4
#define NUM_CASES_3M 6
#define NUM_CASES_4M 10
#define MAX_NUM_CASES 10
#define MAX_FILEEXT_LEN 6
#define MAX_FILE_NAME_LEN 255
#define MAX_NUM_CPU_THREADS 16
#define MAX_NUM_MISMATCH 4

#define OUTPUT_ALL_VALID 1
#define OUTPUT_ALL_BEST 2
#define OUTPUT_UNIQUE_BEST 3
#define OUTPUT_RANDOM_BEST 4

#define SINGLE_READ 1
#define PAIR_END_READ 2

// For reading quality base values
#define DEFAULT_QUAL_CONST 33
#define ILLUMINA_QUAL_CONST 64

// Long read: read with length > 120
#define LONG_READ_LEN 120
#define SOAP3_SEED_LEN 100 // the seed length for the long read for SOAP3
#define MISMATCH_RATIO_FOR_LONG 0.02 // the mismatch ratio allowed for long read

///////////////////
//  Parameters   //
///////////////////

// default parameters
#define DEFAULT_MAX_READ_LEN 120

// default parameters for DP (in general)
#define DEFAULT_NUM_MISMATCH_DP 2 // when DP is enabled (by default), the mismatches # = 2
#define DP_SCORE_THRESHOLD_RATIO 0.3 // DP score threshold = DP_SCORE_THRESHOLD_RATIO * read-length
#define MIN_READ_LEN_FOR_DP 30 // if read length < 30, DP is disabled
#define NUM_CPU_THREADS_FOR_SEEDING 2 // cannot be changed

// default parameters for default DP
#define MAX_SEED_HITS_DEFAULT_DP_FOR_NORMAL_READ 50 // read length > 50
#define MAX_SEED_HITS_DEFAULT_DP_FOR_SHORT_READ 70 // read length <= 50
#define MAX_SEED_HITS_NEW_DEFAULT_DP_FOR_NORMAL_READ 150 // read length > 50
#define MAX_SEED_HITS_NEW_DEFAULT_DP_FOR_SHORT_READ 200 // read length <= 50
// #define MAX_LEN_FROM_HEAD_CLIPPING_DEFAULT_DP 3 // the max length from head of the read allowed for soft-clipping
// #define MAX_LEN_FROM_TAIL_CLIPPING_DEFAULT_DP 8 // the max length from tail of the read allowed for soft-clipping

// default parameters for deep DP (ROUND 1)
#define MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ 100 // read length > 50
#define MAX_SEED_HITS_DEEP_DP_FOR_SHORT_READ 150 // read length <= 50
#define SEED_LEN_DEEP_DP_FOR_V_LONG_READ 45 // read length > 150
#define SEED_LEN_DEEP_DP_FOR_LONG_READ 26 // read length > 80
#define SEED_LEN_DEEP_DP_FOR_MEDIAN_READ 24 // read length > 60 and <= 80
#define SEED_LEN_DEEP_DP_FOR_SHORT_READ 22 // read length > 40 and <= 60
#define SEED_LEN_DEEP_DP_FOR_V_SHORT_READ 20 // read length <= 40
#define SEED_SAMPLE_RATE_DEEP_DP 0.5 // 0.5 of the seed length
#define TAIL_TRIM_SEEDING_DEEP_DP 0 // trail trimmed for seeding
// #define MAX_LEN_FROM_HEAD_CLIPPING_DEEP_DP 0 // cannot be changed
// #define MAX_LEN_FROM_TAIL_CLIPPING_DEEP_DP 0 // cannot be changed

// default parameters for deep DP (ROUND 2)
#define MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ_2 1000 // read length > 50
#define SEED_LEN_DEEP_DP_FOR_V_LONG_READ_2 52 // read length > 150
#define SEED_LEN_DEEP_DP_FOR_LONG_READ_2 30 // read length > 80
#define SEED_LEN_DEEP_DP_FOR_MEDIAN_READ_2 28 // read length > 60 and <= 80
#define SEED_LEN_DEEP_DP_FOR_SHORT_READ_2 26 // read length > 40 and <= 60
#define SEED_LEN_DEEP_DP_FOR_V_SHORT_READ_2 24 // read length <= 40
#define SEED_SAMPLE_RATE_DEEP_DP_2 0.5 // 0.5 of the seed length
#define TAIL_TRIM_SEEDING_DEEP_DP_2 0 // trail trimmed for seeding
// #define MAX_LEN_FROM_HEAD_CLIPPING_DEEP_DP_2 0 // cannot be changed
// #define MAX_LEN_FROM_TAIL_CLIPPING_DEEP_DP_2 0 // cannot be changed


// default parameters for single-end DP
#define SEED_NUM_SINGLE_DP 3 // number of seeds
// #define MAX_LEN_FROM_HEAD_CLIPPING_SINGLE_DP 0 // cannot be changed
// #define MAX_LEN_FROM_TAIL_CLIPPING_SINGLE_DP 0 // cannot be changed
// for read length > 300
#define SEED_LEN_SINGLE_DP_FOR_V_LONG_READ 70 // seed length
#define MAX_SEED_HITS_SINGLE_DP_FOR_V_LONG_READ 4 // max # of hits for seeding
// for read length > 80
#define SEED_LEN_SINGLE_DP_FOR_LONG_READ 38 // seed length
#define MAX_SEED_HITS_SINGLE_DP_FOR_LONG_READ 10 // max # of hits for seeding
#define TAIL_TRIM_SEEDING_SINGLE_DP_FOR_LONG_READ 10 // trail trimmed for seeding
// for read length > 60 and <= 80
#define SEED_LEN_SINGLE_DP_FOR_MEDIAN_READ 32 // seed length
#define MAX_SEED_HITS_SINGLE_DP_FOR_MEDIAN_READ 20 // max # of hits for seeding
#define TAIL_TRIM_SEEDING_SINGLE_DP_FOR_MEDIAN_READ 4 // trail trimmed for seeding
// for read length > 40 and <= 60
#define SEED_LEN_SINGLE_DP_FOR_SHORT_READ 26 // seed length
#define MAX_SEED_HITS_SINGLE_DP_FOR_SHORT_READ 30 // max # of hits for seeding
#define TAIL_TRIM_SEEDING_SINGLE_DP_FOR_SHORT_READ 4 // trail trimmed for seeding
// for read length <= 40
#define SEED_LEN_SINGLE_DP_FOR_V_SHORT_READ 22 // seed length
#define MAX_SEED_HITS_SINGLE_DP_FOR_V_SHORT_READ 40 // max # of hits for seeding
#define TAIL_TRIM_SEEDING_SINGLE_DP_FOR_V_SHORT_READ 0 // trail trimmed for seeding

// default parameters when DP is disabled
#define DEFAULT_NUM_MISMATCH_NO_DP_NORMAL_READ 3 // read length >= 50 (valid only when DP is disabled)
#define DEFAULT_NUM_MISMATCH_NO_DP_SHORT_READ 2 // read length < 50 (valid only when DP is disabled)

// default sample name
#define DEFAULT_SAMPLE_NAME "default"

// Dynamics Programming Parameters..
#define PERFORM_NEW_DEFAULT_DP_FOR_SEMI_ALIGNED_PAIR
#define PERFORM_DEFAULT_DP_FOR_SEMI_ALIGNED_PAIR
#define PERFORM_DEEP_DP_FOR_UNALIGN_PAIR

// uncomment the below to perform deep-dp for unaligned paired-end reads
#define PERFORM_DEEP_DP_FOR_UNALIGN_READS
#define MAX_UNALIGN_READS_NUM_FOR_DEEP_DP 2097152 // 2M

// uncomment the below to enable the constraint
// on the number of hits for each end
// when finding a valid pair of hit
// #define NO_CONSTRAINT_SINGLE_READ_NUM_FOR_PAIRING


// uncomment the below to skip the default DP
// and all the half-aligned reads and unaligned reads
// will proceed to deep DP directly
// #define SKIP_DEFAULT_DP

// DEBUG SWITCH
// ============
//
// Uncomment the below to enable BGS to output
// found SA ranges to standard out.
// Note that sm_20 is needed for this to work.
// #define BGS_OUTPUT_SA_RANGE_TO_SCREEN
//----
//
// Uncomment the below to disable the translation from
// the positions on packed sequence to the positions on chromosomes
// #define DEBUG_2BWT_NO_TRANSLATION
//----
//
// Uncomment the below to show debug messages.
// IT IS HIGHLY RECOMMENDED that this functionality
// is used against only small amount of read input, like below 10,
//// as the amount of output generated is enormous.
//#define BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
//
//#define BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
//----
//
// Uncomment the below to disable the cases.
//#define BGS_DISABLE_CASE_A
//#define BGS_DISABLE_CASE_B
//#define BGS_DISABLE_CASE_C
//#define BGS_DISABLE_CASE_D
//#define BGS_DISABLE_CASE_E
//#define BGS_DISABLE_CASE_F
//#define BGS_DISABLE_CASE_G
//#define BGS_DISABLE_CASE_H
//#define BGS_DISABLE_CASE_I
//#define BGS_DISABLE_CASE_J
//----
//
// Uncomment the below to disable negative strand for both device and host (for single-end alignment only)
// #define BGS_DISABLE_NEGATIVE_STRAND
//----
//
// Uncomment the below to enable extra timing figures.
// #define BGS_GPU_CASE_BREAKDOWN_TIME
// #define BGS_CPU_CASE_BREAKDOWN_TIME
// #define BGS_ROUND_BREAKDOWN_TIME
// #define BGS_CPU_JOBS_BREAKDOWN_TIME
//----
//
// Uncomment the below to enable file name output.
// #define BGS_OUTPUT_FILENAMES_TO_SCREEN
//----
//
// Uncomment the below to enable breakdown in reported occurrences.
// #define BGS_OCC_RESULT_BREAKDOWN
//----
// Uncomment the below to enable memory usage output.
// #define BGS_OUTPUT_MEMORY_USAGE
// #define BGS_OUTPUT_GPU_MEMORY_USAGE
//----
// Uncomment the below to enable bad read information.
// #define BGS_BAD_READS_INFO
//----
// Uncomment the below to output the readID for second iteration (i.e. no-hit reads for trimming).
// #define BGS_OUTPUT_NO_HIT_READS
//----
// Uncomment the below to output the input parameters.
// #define BGS_OUTPUT_PARAMETERS
//----
// Uncomment the below to show DP for aligned reads message.
// #define BGS_OUTPUT_DP_MESSAGE
//----


typedef unsigned int uint;
typedef unsigned long long ullint;

#define STAGE_SINGLE_DP      1
#define STAGE_DEFAULT_DP     2
#define STAGE_NEW_DEFAULT_DP 3
#define STAGE_DEEP_DP_ROUND1 4
#define STAGE_DEEP_DP_ROUND2 5

static void getSeedPositions ( int stage, int readLength, int * seedLength, int * seedPositions, int * seedNum )
{
    if ( stage == STAGE_SINGLE_DP ||
            stage == STAGE_NEW_DEFAULT_DP )
    {
        if ( readLength > 300 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_V_LONG_READ; }
        else if ( readLength > 80 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_LONG_READ; }
        else if ( readLength > 60 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_MEDIAN_READ; }
        else if ( readLength > 40 )
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_SHORT_READ; }
        else
        { ( *seedLength ) = SEED_LEN_SINGLE_DP_FOR_V_SHORT_READ; }

        // update: for reads longer than 120, one more seed for every extra 100 bases
        if ( readLength > 120 )
        {
            ( *seedNum ) = SEED_NUM_SINGLE_DP + readLength / 100;
        }
        else
        {
            ( *seedNum ) = SEED_NUM_SINGLE_DP;
        }

        int H = 0;

        if ( readLength > 300 )
        { H = readLength * 0.15; }

        int X;

        if ( readLength > 300 )
        { X = readLength * 0.15; }
        else if ( readLength > 80 )
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_LONG_READ; }
        else if ( readLength > 60 )
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_MEDIAN_READ; }
        else if ( readLength > 40 )
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_SHORT_READ; }
        else
        { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_V_SHORT_READ; }

        int seed_apart = ( readLength - X - H ) / ( *seedNum );
        int i;

        for ( i = 0; i < ( *seedNum ); i++ )
        {
            seedPositions[i] = H + i * seed_apart;
        }

        if ( seedPositions[ ( *seedNum ) - 1] > readLength - ( *seedLength ) - X )
        { seedPositions[ ( *seedNum ) - 1] = readLength - ( *seedLength ) - X; }
    }
    else if ( stage == STAGE_DEEP_DP_ROUND1 )
    {
        if ( readLength > 150 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_LONG_READ; }
        else if ( readLength > 80 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_LONG_READ; }
        else if ( readLength > 60 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_MEDIAN_READ; }
        else if ( readLength > 40 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_SHORT_READ; }
        else
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_SHORT_READ; }

        int H = 0;
        int T = 0;

        if ( readLength > 150 )
        {
            H = readLength * 0.1;
            T = readLength * 0.2;
        }

        ( *seedNum ) = 0;

        for ( int i = readLength - ( *seedLength ) - T; i >= H; i -= ( *seedLength ) / 2 )
        {
            seedPositions[ ( *seedNum ) ++] = i;
        }

        if ( seedPositions[ ( *seedNum ) - 1] > H )
        { seedPositions[ ( *seedNum ) ++] = H; }
    }
    else if ( stage == STAGE_DEEP_DP_ROUND2 )
    {
        if ( readLength > 150 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_LONG_READ_2; }
        else if ( readLength > 80 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_LONG_READ_2; }
        else if ( readLength > 60 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_MEDIAN_READ_2; }
        else if ( readLength > 40 )
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_SHORT_READ_2; }
        else
        { ( *seedLength ) = SEED_LEN_DEEP_DP_FOR_V_SHORT_READ_2; }

        int H = 0;
        int T = 0;

        if ( readLength > 150 )
        {
            H = readLength * 0.1;
            T = readLength * 0.2;
        }

        ( *seedNum ) = 0;

        for ( int i = readLength - ( *seedLength ) - T; i >= H; i -= ( *seedLength ) / 2 )
        {
            seedPositions[ ( *seedNum ) ++] = i;
        }

        if ( seedPositions[ ( *seedNum ) - 1] > H )
        { seedPositions[ ( *seedNum ) ++] = H; }
    }
}

static unsigned int getWordPerQuery ( unsigned int maxReadLength )
{
    // determine number of words per query. rounded up to power of 2
    unsigned int wordPerQuery = 1;

    while ( wordPerQuery < maxReadLength )
    { wordPerQuery *= 2; }

    wordPerQuery = wordPerQuery / CHAR_PER_WORD;
    return wordPerQuery;
}

#endif
