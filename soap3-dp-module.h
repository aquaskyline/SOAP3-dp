#ifndef __SOAP3_DP_MODULE_H__
#define __SOAP3_DP_MODULE_H__

#define MAX_READ_NAME_LENGTH 256

#include "global_arrays.h"
#include "CPUfunctions.h"
#include "alignment.h"
#include "QueryParser.h"
#include "IndexHandler.h"

// scoring system for DP
typedef struct DPScoring
{
    int cutoffThreshold; // cutoff threshold when performing DP
    int matchScore; // score for match
    int mismatchScore; // score for mismatch
    int openGapScore; // score for gap opening
    int extendGapScore; // score for gap extension
} DPScoring;

// paramaters for paired alignment
typedef struct PairAlignParam
{
    int maxReadLength1;
    int maxReadLength2;
    int numMismatch; // the maximum number of mismatch allowed for the alignment of input reads (default: 2)
    int outputOption; // 1: all valid; 2: all best; 3: unique best; 4: random best
    int insert_low; // minimum value of insert size
    int insert_high; // maximum value of insert size
    int cpuNumThreads; // number of CPU threads
    int enableDP; // 1: DP is enable; 2: DP is disable

    // DP scoring
    DPScoring scoring[2]; // for DP is enabled
} PairAlignParam;

// parameters for single alignment
typedef struct SingleAlignParam
{
    int maxReadLength;
    int numMismatch; // the maximum number of mismatch allowed for the alignment of input reads
    int outputOption; // 1: all valid; 2: all best; 3: unique best; 4: random best
    int cpuNumThreads; // number of CPU threads
    int maxHitNum; // if there are too many hits, then only outputs the first "maxHitNum" hits
    int enableDP; // 1: DP is enable; 2: DP is disable

    // DP scoring
    DPScoring scoring; // for DP is enabled
} SingleAlignParam;



//-------------------------------------------//
// The functions for CX to use               //
//-------------------------------------------//

// To perform the paired alignment
// The resulting alignments are stored inside "algnResultArrays"
void alignPairR ( unsigned int * queries, unsigned int * readLengths, unsigned int * readIDs,
                  unsigned int wordPerQuery,
                  unsigned int numQueries, Soap3Index * index,
                  PairAlignParam * param,
                  unsigned long long & numOfAnswer,
                  unsigned int & numOfAlignedRead );

// To perform the single alignment
// The resulting alignments are stored inside "algnResultArrays"
void alignSingleR ( unsigned int * queries, unsigned int * readLengths, unsigned int * readIDs,
                    unsigned int wordPerQuery,
                    unsigned int numQueries, Soap3Index * index,
                    SingleAlignParam * param,
                    unsigned long long & numOfAnswer,
                    unsigned int & numOfAlignedRead, AlgnResultArrays * algnResultArrays );




#endif
