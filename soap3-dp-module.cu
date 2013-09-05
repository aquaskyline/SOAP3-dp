/*
 *  soap3-dp-module.cpp
 *  soap3_ext_r5
 *
 *  Created by kfwong on 14/3/12.
 *  Copyright 2012 HKU. All rights reserved.
 *
 */


#include "soap3-dp-module.h"

// set Parameters
void setParam ( SingleAlignParam * param, IniParams & ini_params, InputOptions & input_options )
{
    ini_params.Ini_NumOfCpuThreads = param->cpuNumThreads; // number of CPU threads
    // ini_params.Ini_HostAlignmentModelStr = "16G";
    ini_params.Ini_HostAlignmentModel = 1; // 16G
    ini_params.Ini_GPUMemory = 6; // 6G
    ini_params.Ini_PEStrandLeftLeg = 1; // pos strand
    ini_params.Ini_PEStrandRightLeg = 2; // neg strand
    ini_params.Ini_MaxOutputPerRead = param->maxHitNum;
    // if there are too many hits, then only outputs the first "maxHitNum" hits
    ini_params.Ini_PEMaxOutputPerPair = 1000; // not applicable

    ini_params.Ini_shareIndex = 0; // index to be shared among multiple copies of soap3-dp
    ini_params.Ini_maxReadNameLen = 64; // max length allowed for read name
    // max length from the front of the read for clipping
    ini_params.Ini_maxFrontLenClipped = 3;
    // max length from the end of the read for clipping
    ini_params.Ini_maxEndLenClipped = 8;
    // whether the seed will proceed to perform DP if there are too many hits
    ini_params.Ini_proceedDPForTooManyHits = 0;
    // whether the read will perform SOAP3 module
    ini_params.Ini_skipSOAP3Alignment = 0;

    if ( param->enableDP == 1 )
    {
        ini_params.Ini_MatchScore = param->scoring.matchScore;
        ini_params.Ini_MismatchScore = param->scoring.mismatchScore;
        ini_params.Ini_GapOpenScore = param->scoring.openGapScore;
        ini_params.Ini_GapExtendScore = param->scoring.extendGapScore;
        ini_params.Ini_DPScoreThreshold = param->scoring.cutoffThreshold;
        ini_params.Ini_isDefaultThreshold = 0;
    }

    input_options.maxReadLength = param->maxReadLength;
    input_options.outputFormat = 1; // plain format (fixed)
    input_options.enableDP = param->enableDP;
    input_options.numMismatch = param->numMismatch;
    input_options.alignmentType = param->outputOption;
    input_options.insert_low = -1; // not applicable
    input_options.insert_high = -1; // not applicable
    input_options.readType = 1; // single-end
    input_options.isReadList = 0; // not applicable
}


// To perform the single alignment
// The resulting alignments are stored inside "algnResultArrays"
// "algnResultArrays" has to be constructed before calling this function
void alignSingleR ( unsigned int * queries, unsigned int * readLengths, unsigned int * readIDs,
                    unsigned int wordPerQuery,
                    unsigned int numQueries, Soap3Index * index,
                    SingleAlignParam * param,
                    unsigned long long & numOfAnswer,
                    unsigned int & numOfAlignedRead,
                    AlgnResultArrays * algnResultArrays )
{
    uint * _bwt, *_occ;
    uint * _revBwt, *_revOcc;
    double startTime, copyTime;
    double lastEventTime;
    double totalAlignmentTime = 0.0;
    //Start measuring runtime..
    startTime = setStartTime ();
    lastEventTime = 0;
    // Indicate whether the index has been loaded to GPU
    uint indexLoadedToGPU = 0;

    if ( indexLoadedToGPU == 0 )
    {
        GPUINDEXUpload ( index, &_bwt, &_occ,
                         &_revBwt, &_revOcc );
        copyTime = getElapsedTime ( startTime );
        printf ( "[Main] Finished copying index into device (GPU).\n" );
        printf ( "[Main] Loading time : %9.4f seconds\n\n", copyTime - lastEventTime );
        lastEventTime = copyTime;
        indexLoadedToGPU = 1;
    }

    // ======================================================================================
    // | ALLOCATE MEMORY FOR THE ARRAYS                                                     |
    // ======================================================================================
    ullint roundUp = ( numQueries + 31 ) / 32 * 32;
    uint maxReadLength = param->maxReadLength;
    char * upkdQualities = ( char * ) malloc ( roundUp * maxReadLength * sizeof ( char ) ); // may not need
    memset ( upkdQualities, 0, roundUp * maxReadLength ); // may not need
    unsigned int * unAlignedPair = ( unsigned int * ) malloc ( roundUp * sizeof ( unsigned int ) );
    char * upkdQueryNames = ( char * ) malloc ( roundUp * MAX_READ_NAME_LENGTH * sizeof ( char ) );
    uint maxBatchSize = NUM_BLOCKS * THREADS_PER_BLOCK * QUERIES_PER_THREAD; // queries processed in one kernel call
    // maxBatchSize has to be divible by 2
    maxBatchSize = maxBatchSize / 2 * 2;
    // set parameters
    InputOptions input_options;
    IniParams ini_params;
    setParam ( param, ini_params, input_options );
    index->sraIndex->hspaux->algnResultArrays = algnResultArrays;
    // Declare the structure for storing the alignment results
    // for the pairs of reads with one end has no hit but another has.
    // The structure will be used for proceeding semi-global DP
    ReadInputForDPArrays readInputForDPall;
    readInputForDPall.inputArrays = ( ReadInputForDP ** ) malloc ( ini_params.Ini_NumOfCpuThreads * sizeof ( ReadInputForDP * ) );
    readInputForDPall.numArrays = ini_params.Ini_NumOfCpuThreads;

    for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        readInputForDPall.inputArrays[threadId] = constructReadInputForDP ( ini_params.Ini_NumOfCpuThreads );
    }

    // Declare the structure for storing the read IDs
    // for those read cannot be aligned
    UnalignedSinglesArrays * unalignedSinglesArrays = constructBothUnalignedPairsArrays ( ini_params.Ini_NumOfCpuThreads + 1 );
    // for single-end reads
    // get the max read length for the first ten reads (i.e. 0, 1, ..., 9)
    // ==================================================================
    // | DETECT THE READ LENGTH                                         |
    // ==================================================================
    uint detected_read_length = GetReadLength ( readLengths, numQueries, 1 );
    numOfAnswer = 0;
    numOfAlignedRead = 0;
    unsigned int numOfUnAlignedPairs = 0;
    unsigned int accumReadNum = 0;
    // ======================================================================================
    // | Configuration on GPU functions                                                     |
    // ======================================================================================
    // cudaFuncSetCacheConfig(kernel, cudaFuncCachePreferShared);
    // cudaFuncSetCacheConfig(kernel_4mismatch_1, cudaFuncCachePreferShared);
    // cudaFuncSetCacheConfig(kernel_4mismatch_2, cudaFuncCachePreferShared);
    // ======================================================================================
    // | Perform alignment                                                                  |
    // ======================================================================================
    soap3_dp_single_align ( queries, readLengths, param->numMismatch, wordPerQuery,
                            maxBatchSize, numQueries, accumReadNum,
                            index, _bwt, _revBwt,
                            _occ, _revOcc,
                            ini_params, input_options,
                            param->maxReadLength, detected_read_length,
                            upkdQualities,
                            unAlignedPair, numOfUnAlignedPairs,
                            readIDs, upkdQueryNames,
                            NULL, NULL,
                            NULL, NULL,
                            numOfAnswer, numOfAlignedRead,
                            &readInputForDPall,
                            unalignedSinglesArrays,
                            startTime, lastEventTime, totalAlignmentTime,
                            indexLoadedToGPU );
    // ======================================================================================
    // | CLEAN UP                                                                           |
    // ======================================================================================
    free ( upkdQualities );
    free ( upkdQueryNames );
    free ( unAlignedPair );

    // free device memory
    if ( indexLoadedToGPU == 1 )
    {
        printf ( "[Main] Free device memory..\n" );
        GPUINDEXFree ( _bwt, _occ, _revBwt, _revOcc );
    }

    for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        freeReadInputForDP ( readInputForDPall.inputArrays[threadId] );
    }

    free ( readInputForDPall.inputArrays ); // for half-aligned pairs
    freeBothUnalignedPairsArrays ( unalignedSinglesArrays ); // for single-unaligned reads
}


// pending.... NOT YET IMPLEMENTED
// To perform the paired alignment
// The resulting alignments are stored inside "algnResultArrays"
void alignPairR ( unsigned int * queries, unsigned int * readLengths, unsigned int * readIDs,
                  unsigned int wordPerQuery,
                  unsigned int numQueries, Soap3Index * index,
                  PairAlignParam * param,
                  unsigned long long & numOfAnswer,
                  unsigned int & numOfAlignedRead )
{
}



