/*
 *
 *    PEAlgnmt.h
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

/////////////////////////////////////////////////////
/*

            Each pair-end result is considered to have two legs
            which each of them is a position on the reference sequence.

            The left leg is the position of the left-most position of
            the match (which must be a %strandLeftLeg% strand alignment) ;
            the right leg is the position of its mate
            (which must be a %strandRightLeg% strand alignment);

                | left-leg                  | right-leg
                v                           v
            ___|XXXXXXXXXXXXX|_____________|XXXXXXXXXXXXX|_

                |<----------- insertion size ---------->|
*/
/////////////////////////////////////////////////////

#ifndef __PE_ALIGNMENT_H__
#define __PE_ALIGNMENT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "2bwt-flex/SRACore.h"
#include "HSPAux.h"

#define PE_MAX_BUCKET_SIZE 1024

#define PE_REPORT_ALL      0
#define PE_REPORT_ONE      1

#define PE_ALIGNMENT_COMPLETED    0
#define PE_ALIGNMENT_INITIALISED  1
#define PE_ALIGNMENT_INPUT_ERROR  2
#define PE_ALIGNMENT_FAILED       3

#define PE_NON_INPUT_VALUE 99999999

#define INITIAL_SIZE_SA_LIST_FOR_DP 1048576 // 1M
#define INITIAL_SIZE_OCC_LIST_FOR_DP 1048576 // 1M

#define MAX_SEED_NUMBER_FOR_SINGLE_DP 10

#define INITIAL_SIZE_CHAR_ARRAY 40000 // 40K    //1048576 // 1M
#define SIZE_1_M 1048576


// Uncomment the below two sentences to ensure the first read for forward direction only
// and the second read for reverse direction only if paired-end reads are inputted.
// this option is only valid for simple output format and all-best alignment
// #define BGS_FWD_FOR_FIRST_REV_FOR_SECOND
// #define BGS_DISABLE_NEGATIVE_STRAND
// NOTE: Remember to delete the file GPUfunctions.cpp before recompile the codes.
//----



/////////////////////////////////////////////////////
/*
    Structure SRAOccurrence / PESRAAlignmentResult (Input)
    They are the structures to hold the input to Pair-End alignment.

    SRAOccurrence is supposed to hold a list of occurrences.
    The function PEMappingOccurrences() takes two lists of SRAOccurrence as input.

    PESRAAlignmentResult is supposed to hold a list of SA ranges.
    The function PEMapping() takes two lists of PESRAAlignmentResult as input.
*/
/////////////////////////////////////////////////////
/*
// SRAOccurrence is part of SRA now
typedef struct SRAOccurrence
{
    unsigned int readID; // read ID
    unsigned int ambPosition;
    char strand; // (i.e. 1: +ve; 2: -ve)
    char mismatchCount;
} SRAOccurrence;
*/

typedef struct PESRAAlignmentResult
{
    unsigned int readID; // read ID
    unsigned int saIndexLeft;
    unsigned int saIndexRight;
    unsigned char strand; // (i.e. 1: +ve; 2: -ve)
    unsigned char mismatchCount;
} PESRAAlignmentResult;

/////////////////////////////////////////////////////
/*
    Structure PEInput (Parameters)
    It is the structure to hold the query parameters to Pair-End alignment.
    Some of the parameters are compulsory; some of them are optional.
    PEInput SHOULD ALWAYS BE constructed by PEInputConstruct() function.
*/
/////////////////////////////////////////////////////
typedef struct PEInput
{
    BWT * bwt;                  //Compulsory
    HSP * hsp;                  //Compulsory
    int OutputType;             //Compulsory. Valid value being PE_REPORT_*.

    int patternLength;          //Compulsory

    int strandLeftLeg;          //Compulsory
    int strandRightLeg;         //Compulsory

    // Compulsory
    ///////////////////////////////////////////////////
    // Either option A or B should be defined and it is compulsory
    // Option B takes precedence over option A.
    //
    // Option A: Insertion size defined by Mean and StdDev.
    int insertSizeMean;         //Optional.
    int insertSizeStdDev;       //Optional.
    // Option B: Insertion size defined by Lower bound and Upper bound.
    int insertLbound;           //Optional.
    int insertUbound;           //Optional.
    ///////////////////////////////////////////////////
} PEInput;


/////////////////////////////////////////////////////
/*
    Structure PEPairs (Output)
    It is the structure to hold 1 Pair-End alignment result.
    The result being held is w.r.t to the read set source.
    The alignment/strand from the first read set is stored in
    algnmt_1/strand_1 while those from the second read set is
    stored in algnmt_2/strand_2.
*/
/////////////////////////////////////////////////////
typedef struct PEPairs
{

    unsigned int algnmt_1;
    char strand_1;
    char mismatch_1;

    unsigned int algnmt_2;
    char strand_2;
    char mismatch_2;

    char totalMismatchCount;
    int insertion;

} PEPairs;


/////////////////////////////////////////////////////
/*
    Structure PEPairList (Output)
    It is the structure to hold a bucket of alignment result(PEPairs).
*/
/////////////////////////////////////////////////////
typedef struct PEPairList
{
    // The bucket of pair-end results
    PEPairs pairs[PE_MAX_BUCKET_SIZE];
    unsigned int pairsCount;

    // Link to the next bucket
    struct PEPairList * next;
} PEPairList;


/////////////////////////////////////////////////////
/*
    Structure PEOutput (Output)
    It is the structure to enclose a linked-list of bucket(PEPairList).
*/
/////////////////////////////////////////////////////
typedef struct PEOutput
{
    int flag;
    struct PEPairList * root;
    struct PEPairList * tail;
} PEOutput;


/////////////////////////////////////////////////////
/*
    Function SRAEnrichSARanges

    This function fills the gap between SOAP3-GPU and the PE Aligner
    by enriching each reported SA range with the information of the
    strand and mismatchCount of the alignment.
*/
/////////////////////////////////////////////////////
/*
unsigned int SRAEnrichSARanges(BWT * bwt, HSP * hsp,
                                unsigned int saIndexLeft, unsigned int saIndexRight,
                                char * strand, int * mismatchCount);

*/
/////////////////////////////////////////////////////
/*
    Function PEMapping
    This function takes in 2 lists of SA ranges, which
    are assumed to be enriched by SRAEnrichSARanges.

    This function will perform the following,
    1. Initialise output collector
    2. Retrieve the occurrences of all SA ranges
    3. Allocate enough memory for all occurrences
    4. Sort the 2 lists of occurrences
    5. Merge the occurrences to generate a list of PE
       alignment

    Parameters:
        peInput - The input parameter for PE
        peOutput - The output of the PE alignment, which is
                   expected to be constructed by the caller.
                   It is initialised again by PEMapping to allow
                   reusing of constructed PEOutput.
        resultListA,
        resultCountA - The first list of SA ranges being
                       passed into PEMapping.
        resultListB,
        resultCountB - The second list of SA ranges being
                       passed into PEMapping.

*/
/////////////////////////////////////////////////////
void PEMapping ( PEInput * peInput, PEOutput * peOutput,
                 PESRAAlignmentResult * resultListA, unsigned int resultCountA,
                 PESRAAlignmentResult * resultListB, unsigned int resultCountB );


/////////////////////////////////////////////////////
/*
    Function PEMappingOccurrences
    This function is basically behaves the same as the above function,
    however it takes in 2 lists of Occurrences, instead of two lists of
    SA Ranges.

    This function will perform the following,
    1. Allocate enough memory for sorting all occurrences
    2. Sort the 2 lists of occurrences
    3. Merge the occurrences to generate a list of PE
       alignment

    Parameters:
        peInput  - The input parameter for PE
        peOutput - The output of the PE alignment, which is
                   expected to be constructed by the caller.
                   It is initialised again by PEMapping to allow
                   reusing of constructed PEOutput.
        occList_1,
        occCount_1 - The first list of occurrences
                   being passed into PEMapping.
        occList_2,
        occCount_2 - The second list of occurrences
                   being passed into PEMapping.

*/
/////////////////////////////////////////////////////
void PEMappingOccurrences ( PEInput * peInput, PEOutput * peOutput,
                            SRAOccurrence * occList_1, unsigned int occCount_1,
                            SRAOccurrence * occList_2, unsigned int occCount_2 );


// unsigned int PEMappingToDisk(PEInput * peInput, PEPairList * pePairList);


/////////////////////////////////////////////////////
// Constructor and Destructor
/////////////////////////////////////////////////////
PEInput * PEInputConstruct ( BWT * bwt, HSP * hsp );
void PEInputFree ( PEInput * peInput );
PEOutput * PEOutputConstruct ();
void PEOutputFree ( PEOutput * peOutput );

/////////////////////////////////////////////////////
// Utilities
/////////////////////////////////////////////////////
unsigned int PECountPEOutput ( PEOutput * peOutput );
unsigned int PEStatsPEOutput ( PEOutput * peOutput, PEPairs ** optimal, PEPairs ** suboptimal, unsigned int * mismatchStats );

///////////////////////////////////////////////////////////
// Structures for semi-global DP                         //
///////////////////////////////////////////////////////////


// To hold SOAP3 alignment results of the reads which can be mapped but their mates cannot be mapped
// these reads are required to be processed by semi-global DP
// this array is also used for the input of gap alignment
typedef struct ReadInputForDP
{
    unsigned int readNum; // number of reads
    PESRAAlignmentResult * sa_list; // the list of SA ranges for all reads
    SRAOccurrence * occ_list; // the list of occurrences for all reads
    unsigned int saRangeTotalNum; // the number of SA Ranges in sa_list
    unsigned int occTotalNum; // the number of occurrences in occ_list
    unsigned int saListSize; // the size of sa_list
    unsigned int occListSize; // the size of occ_list
    unsigned int lastReadID; // the ID (i.e. even id) of the latest pair inputted in the array
} ReadInputForDP;

// To hold a set of arrays of ReadInputForDP
typedef struct ReadInputForDPArrays
{
    ReadInputForDP ** inputArrays;
    unsigned int numArrays; // number of arrays of ReadInputForDP
} ReadInputForDPArrays;


// To contain the parameters for DP

typedef struct DPParam
{
    int cutoffThreshold; // cutoff threshold when performing DP
    int maxHitNum; // maximum number of hits allowed for soap3 seeding
    int sampleDist; // sample distane of seeding (i.e. sample rate = 1/sampleDist)
    int seedLength; // length of each seed for seeding
} DPParam;

typedef struct DPParameters
{
    int matchScore; // score for match
    int mismatchScore; // score for mismatch
    int openGapScore; // score for gap opening
    int extendGapScore; // score for gap extension
    int numOfCPUThreads; // number of CPU threads used for matrix computation
    int numOfCPUForSeeding; // number of CPU threads used for seeding
    int softClipLeft; // left
    int softClipRight; // right
    int tailTrimLen; // the length of tail trimmed for seeding
    int singleDPSeedNum; // the number of seeds for single-end DP
    int singleDPSeedPos[MAX_SEED_NUMBER_FOR_SINGLE_DP]; // the seed positions for single-end DP
    DPParam paramRead[2]; // the specific paramters for the first and the second reads.
    // for single-end DP, use the first one.
} DPParameters;

// Construct the structure ReadInputForDP for each CPU thread
ReadInputForDP * constructReadInputForDP ( int cpuNum );

// Release the memory for the structure ReadInputForDP
void freeReadInputForDP ( ReadInputForDP * readInput );

// Reset the structure ReadInputForDP
void resetReadInputForDP ( ReadInputForDP * readInput );

// To add the alignment results of a read to ReadInputForDP
void addToReadInputForDP ( ReadInputForDP * readInput, unsigned int readid, PESRAAlignmentResult * saList, unsigned int saNum,
                           SRAOccurrence * occList, unsigned int occNum );

//////////////////////////////////////////////////////////
//// The followings structure is for DEFAULT DP    ///////
//////////////////////////////////////////////////////////

// To hold one pair-end alignment result
typedef struct AlgnmtDPResult
{
    unsigned int readID; // read ID of first read (i.e. even)
    char whichFromDP; // 0: first read; 1: second read; 2: none
    char * cigarString; // alignment pattern of the hit from DP
    int editdist; // edit distance of the hit from DP
    int insertSize; // insert size
    int num_sameScore;

    // alignment result of the first read
    unsigned int algnmt_1; // 0xFFFFFFFF if unaligned
    char strand_1; // 1: +ve; 2: -ve
    int score_1; // # of mismatches or score from DP

    // alignment result of the second read
    unsigned int algnmt_2; // 0xFFFFFFFF if unaligned
    char strand_2; // 1: +ve; 2: -ve
    int score_2; // # of mismatches or score from DP
} AlgnmtDPResult;

///////////////////////////////////////////////////////
//// The followings structure is for DEEP DP    ///////
///////////////////////////////////////////////////////

// To hold one pair-end alignment result of Deep DP for unaligned pairs
typedef struct DeepDPAlignResult
{
    unsigned int readID; // read ID of first read (i.e. even)
    int insertSize; // insert size

    // alignment result of the first read
    unsigned int algnmt_1; // 0xFFFFFFFF if unaligned
    char strand_1; // 1: +ve; 2: -ve
    int score_1; // score from DP
    int editdist_1;
    char * cigarString_1; // alignment pattern of the hit from DP
    int num_sameScore_1;

    // alignment result of the second read
    unsigned int algnmt_2; // 0xFFFFFFFF if unaligned
    char strand_2; // 1: +ve; 2: -ve
    int score_2; // score from DP
    int editdist_2;
    char * cigarString_2; // alignment pattern of the hit from DP
    int num_sameScore_2;
} DeepDPAlignResult;



// To hold one single-end alignment result
typedef struct SingleAlgnmtResult
{
    unsigned int readID; // read ID
    char * cigarString; // alignment pattern of the hit from DP

    // alignment result of the read
    unsigned int algnmt; // alignment pos; 0xFFFFFFFF if unaligned
    char strand; // 1: +ve; 2: -ve
    int score; // score from DP
    int editdist;
    int num_sameScore;
} SingleAlgnmtResult;



// dynamic-size character array
typedef struct DynamicUint8Array
{
    uint8_t * charStr; // the character array
    unsigned long long length; // length of the array
    unsigned long long size; // available size of the char array
} DynamicUint8Array;

DynamicUint8Array * DynamicUint8ArrayConstruct ();
void DynamicUint8ArrayFree ( DynamicUint8Array * uint8Array );
void DynamicUint8ArrayReset ( DynamicUint8Array * uint8Array );
void appendStringToUint8Array ( DynamicUint8Array * uint8Array, char * charString, int len );


// ================================================
// To hold a set of single-end alignment results
// ================================================

// each algnment
typedef struct Algnmt
{
    char * cigarString; // alignment pattern of the hit from DP

    // an alignment of the read
    unsigned int algnmt; // alignment pos; 0xFFFFFFFF if unaligned
    char strand; // 1: +ve; 2: -ve
    int score; // score from DP
    int editdist;
    int num_sameScore;
    int isFromDP; // 1: from DP; 0: from SOAP3
} Algnmt;

// read pointer
typedef struct ReadPtr
{
    unsigned int readID;
    unsigned int startIndex; // start index of the array of algnmt
    unsigned int numAlgnmt; // number of alignments
} ReadPtr;

// the structure to hold a set of alignment results
typedef struct AllHits
{
    Algnmt * hitArray;
    int hitNum;
    int hitArrayAvailSize;
    ReadPtr * readPtrArray;
    int readNum;
    int readArrayAvailSize;
} AllHits;

AllHits * constructAllHits ();
void inputAlgnmtsToArray ( AllHits * allHits, SingleAlgnmtResult * algnResults, int algnNum );
void sortReadPtrs ( AllHits * allHits );
void resetAllHits ( AllHits * allHits );
void releaseAllHits ( AllHits * allHits );

void printOutHits ( AllHits * allHits ); // for debugging

void inputSoap3AnsToArray ( AllHits * allHits, unsigned int readID, HSPAux * hspaux, BWT * bwt, int readLength );
// input soap3 answer to array
// soap3 answer should be stored in HSP->soap3AnsArray

// ================================================


#endif

