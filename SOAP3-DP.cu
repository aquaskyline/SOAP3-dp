/*
 *
 *    SOAP3-DP.cu
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <zlib.h>


#include "Release.h"
#include "DV-Kernel.h"
#include "CPUfunctions.h"
#include "alignment.h"
#include "SAM.h"
#include "samtools-0.1.18/bam.h"

#include "BGS-IO.h"
#include "2bwt-flex/SRAArguments.h"
#include "PEAlgnmt.h"
#include "AlgnResult.h"
#include "OutputDPResult.h"
#include "DV-SemiDP.h"
#include "DV-DPForBothUnalign.h"
#include "DV-DPForSingleReads.h"

#include "aio_thread.h"

int main ( int argc, char ** argv )
{
    // ======================================================================================
    // | VARIABLE DECLARATION                                                               |
    // ======================================================================================
    // local variables used in main, like BWT indexes and such.
    int i;
    double startTime, indexLoadTime, readLoadTime;
    double lastEventTime;
    double totalReadLoadTime = 0.0;
    double totalAlignmentTime = 0.0;
    double totalTrimAlignmentTime = 0.0;
    char * queryFileName = "";
    char * queryFileName2 = "";
    FILE * outputFile;
    FILE * outputDPFile;
    FILE * outputDoneFile;
    // FILE *queryFile;
    // FILE *queryFile2;
    //gzFile * gzQueryFile;
    //gzFile * gzQueryFile2;
    gzFile gzQueryFile;
    gzFile gzQueryFile2;
    bamFile bamQueryFile;
    bam_header_t * bamHeader;
    bam1_t * bam;
    char isFastq = 0;
    // user-specified maximum read length
    // and number of words per query
    uint maxReadLength;
    uint wordPerQuery;
    // uint numQueries;
    ullint roundUp;
    ullint totalQueryLength;
    uint * queries;
    uint * readLengths;
    uint * readIDs;
    char * upkdQueryNames;
    char * upkdQualities;
    uint * queries0;
    uint * readLengths0;
    uint * readIDs0;
    char * upkdQueryNames0;
    char * upkdQualities0;
    uint * queries1;
    uint * readLengths1;
    uint * readIDs1;
    char * upkdQueryNames1;
    char * upkdQualities1;
    unsigned int * unAlignedPair;
    char queryFileBuffer[INPUT_BUFFER_SIZE];
    char queryFileBuffer2[INPUT_BUFFER_SIZE];
    unsigned long long numOfAnswer;
    unsigned int numOfAlignedRead;
    unsigned int numOfUnAlignedPairs;
    uint detected_read_length = 0;
    // for single-end reads
    // it represents the max read length for the first ten reads (i.e. 0, 1, ..., 9)
    // for paired-end reads
    // the max read length for the first ten reads with even readIDs (i.e. 0,2,...,18)
    uint detected_read_length2 = 0;
    // for paired-end reads only
    // the max read length for the first ten reads with odd readIDs (i.e. 1,3,...,19)
    // Declare variables and set up preference for device.
    // #ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    // cudaEvent_t start, stop;
    // float time1, totalDeviceTime;
    // #endif
    uint * _bwt, *_occ;
    uint * _revBwt, *_revOcc;
    // Accumulated number of reads aligned and alignments
    ullint totalReadsAlignedForInputReads = 0;
    ullint totalAnsForInputReads = 0;
    // Input parameters
    InputOptions input_options;
    // Indicate whether the index has been loaded to GPU
    uint indexLoadedToGPU = 0;
    // The inputs for multi option
    MultiInputItem * multiInput = NULL;
    int expNum = 0;
    int currExp = 0;
    // ======================================================================================
    // | Configuration on GPU functions                                                     |
    // ======================================================================================
    // not much effect, also this is not supported by old version of cuda library
    // like the machines in TJ
    // thus these are depreciated
    cudaDeviceSetCacheConfig ( cudaFuncCachePreferL1 );
    // cudaFuncSetCacheConfig(kernel, cudaFuncCachePreferShared);
    // cudaFuncSetCacheConfig(kernel_4mismatch_1, cudaFuncCachePreferShared);
    // cudaFuncSetCacheConfig(kernel_4mismatch_2, cudaFuncCachePreferShared);
    // ======================================================================================
    // | PARSING CONFIGURATION FILE                                                         |
    // ======================================================================================
    char * iniFileName;
    iniFileName = ( char * ) malloc ( strlen ( argv[0] ) + 5 );
    strcpy ( iniFileName, argv[0] );
    strcpy ( iniFileName + strlen ( argv[0] ), ".ini" );
    IniParams ini_params;

    if ( ParseIniFile ( iniFileName, ini_params ) != 0 )
    {
        fprintf ( stderr, "Failed to parse config file ... %s\n", iniFileName );
        return 1;
    }

    printf ( "\n[Main] %s v%d.%d.%d (%s)\n", PROJECT_NAME, PROJECT_MAJOR, PROJECT_MINOR, PROJECT_REV, PROJECT_SPECIAL );
    // printf("[Main] Finished parsing ini file %s.\n\n", iniFileName);
    free ( iniFileName );
    // ======================================================================================
    // | CHECK THE INPUT ARGUMENTS                                                          |
    // ======================================================================================
    bool inputValid = parseInputArgs ( argc, argv, input_options );

    if ( !inputValid )
    { exit ( 1 ); }

    // for multi mode
    if ( input_options.isReadList == 1 )
    {
        multiInput = loadMultiInputFile ( input_options.queryFileName, ( input_options.readType == PAIR_END_READ ) ? 1 : 0,
                                          input_options.isReadBAM, expNum );
        updateInputOption ( ( &input_options ), multiInput, currExp++ );
    }

    // get the name of the query file (and the second query file for pair-ended reads)
    queryFileName = input_options.queryFileName;

    if ( input_options.readType == PAIR_END_READ && input_options.isReadBAM == 0 )
    {
        queryFileName2 = input_options.queryFileName2;
    }

    maxReadLength = input_options.maxReadLength;

    if ( input_options.readType == SINGLE_READ )
    { printf ( "[Main] Loading read file %s\n", queryFileName ); }
    else if ( input_options.isReadBAM == 0 )
    { printf ( "[Main] Loading read files %s and %s\n", queryFileName, queryFileName2 ); }

    // ======================================================================================
    // | Restriction on GPU card with 3G memory                                             |
    // ======================================================================================
    // restriction on the maxReadLength when using GPU card with 3G memory
    if ( ( ini_params.Ini_GPUMemory == 3 ) && ( maxReadLength > 128 ) )
    {
        printf ( "For GPU card with 3G memory, the program cannot support maximum read length more than 128.\n" );
        exit ( 1 );
    }

    // restriction on the number of hits of each end for pairing
#ifdef NO_CONSTRAINT_SINGLE_READ_NUM_FOR_PAIRING

    if ( input_options.readType == PAIR_END_READ )
    { ini_params.Ini_MaxOutputPerRead = 0xFFFFFFFF; }

#endif

    // set the maximum number of mismatches allowed for soap3 alignment
    // if DP is enabled.
    if ( input_options.enableDP == 1 )
    {
        input_options.numMismatch = ini_params.Ini_Soap3MisMatchAllow;
    }

    // ======================================================================================
    // | Selection of GPU device                                                            |
    // ======================================================================================

    if ( input_options.GPUDeviceID > -1 )
    { cudaSetDevice ( input_options.GPUDeviceID ); }

    // ======================================================================================
    // | VARIABLES SETTING AND INITIALISATION                                               |
    // ======================================================================================
    // determine number of words per query. rounded up to power of 2
    wordPerQuery = 1;

    while ( wordPerQuery < maxReadLength )
    { wordPerQuery *= 2; }

    wordPerQuery = wordPerQuery / CHAR_PER_WORD;
    uint maxNumQueries = MAX_NUM_BATCH * NUM_BLOCKS * THREADS_PER_BLOCK;
    // maxNumQueries has to be divible by 2
    maxNumQueries = maxNumQueries / 2 * 2;
    // For Output filenames
    char * outputFileName[MAX_NUM_CPU_THREADS];

    for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
    { outputFileName[i] = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 10 ); }

    // For Output of the DP results
    char * outputDPFileName;
    outputDPFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 9 );
    // For Output of the unpaired results
    char * outputUnpairFileName;
    outputUnpairFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 8 );
    // For Output the DONE file
    char * outputDoneFileName;
    outputDoneFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 6 );

    for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
    {
        sprintf ( outputFileName[i], "%s.gout.%d", input_options.outputPrefix, i + 1 );
    }

    sprintf ( outputDPFileName, "%s.dpout.1", input_options.outputPrefix );
    sprintf ( outputUnpairFileName, "%s.unpair", input_options.outputPrefix );
    sprintf ( outputDoneFileName, "%s.done", input_options.outputPrefix );
    // ======================================================================================
    // | STRUCTURES FOR SEMI-GLOBAL ALIGNMENT DP                                            |
    // ======================================================================================
    // Declare the structure for storing the alignment results
    // for the pairs of reads with one end has no hit but another has.
    // The structure will be used for proceeding semi-global DP
    ReadInputForDPArrays readInputForDPall;
    readInputForDPall.inputArrays = ( ReadInputForDP ** ) malloc ( MAX_NUM_CPU_THREADS * sizeof ( ReadInputForDP * ) );
    readInputForDPall.numArrays = ini_params.Ini_NumOfCpuThreads;
    // The following structure will be used for proceeding NEW semi-global DP
    ReadInputForDPArrays readInputForNewDPall;
    readInputForNewDPall.inputArrays = ( ReadInputForDP ** ) malloc ( MAX_NUM_CPU_THREADS * sizeof ( ReadInputForDP * ) );
    readInputForNewDPall.numArrays = ini_params.Ini_NumOfCpuThreads;
    // Declare the structure for storing the read IDs of the first end of the pairs
    // which both ends cannot be aligned (if the input read is paired-end)
    BothUnalignedPairsArrays * bothUnalignedPairsArrays = constructBothUnalignedPairsArrays ( ini_params.Ini_NumOfCpuThreads + 1 );
    // OR for single end cannot be aligned (if the input read is single-end)
    UnalignedSinglesArrays * unalignedSinglesArrays = bothUnalignedPairsArrays; // same array but just different name
    // Declare the structure for storing the alignment results
    // for the pairs of reads with both ends have hits but do not have valid insert size or proper strands.
    ReadInputForDPArrays readAnsForNonValidInsert;
    readAnsForNonValidInsert.inputArrays = ( ReadInputForDP ** ) malloc ( MAX_NUM_CPU_THREADS * sizeof ( ReadInputForDP * ) );
    readAnsForNonValidInsert.numArrays = ini_params.Ini_NumOfCpuThreads;
    // Parameters for DP
    DPParameters dpParameters;
    getParameterForAllDP ( dpParameters, input_options, ini_params );
    // ======================================================================================
    // | INDEX LOADING                                                                      |
    // ======================================================================================
    //Start measuring runtime..
    startTime = setStartTime ();
    lastEventTime = startTime;
    Soap3Index * index = INDEXLoad ( &ini_params, input_options.indexName, ini_params.Ini_shareIndex );
    HSP * hsp = index->sraIndex->hsp;
    HSPAux * hspaux = index->sraIndex->hspaux;
    indexLoadTime = getElapsedTime ( startTime );
    printf ( "[Main] Finished loading index into host.\n" );
    printf ( "[Main] Loading time : %9.4f seconds\n", indexLoadTime );
    printf ( "[Main] Reference sequence length : %u\n\n", index->sraIndex->bwt->textLength );
    lastEventTime = indexLoadTime;
    roundUp = ( maxNumQueries + 31 ) / 32 * 32;
    totalQueryLength = roundUp * wordPerQuery;
    // ==============================================================
    // | QUALITY CONSTANT and DP MATCH SCORE and alignment type     |
    // ==============================================================
    int quality_constant = DEFAULT_QUAL_CONST;

    if ( input_options.isIlluminaQual == 1 )
    { quality_constant = ILLUMINA_QUAL_CONST; }

    hspaux->dpMatchScore = ini_params.Ini_MatchScore;
    hspaux->dpMisMatchScore = ini_params.Ini_MismatchScore;
    hspaux->alignmentType = input_options.alignmentType;
    hspaux->readType = input_options.readType;
    hspaux->peMaxOutputPerRead = ini_params.Ini_PEMaxOutputPerPair;

    if ( input_options.readType == SINGLE_READ )
    { hspaux->ProceedDPForTooManyHits = ini_params.Ini_proceedDPForTooManyHits; }
    else
    { hspaux->ProceedDPForTooManyHits = 0; }

    // For Mapping Quality Score Calculation
    hspaux->x0_array = ( int * ) malloc ( roundUp * sizeof ( int ) );
    hspaux->x1_array = ( int * ) malloc ( roundUp * sizeof ( int ) );
    hspaux->mismatch_array = ( int * ) malloc ( roundUp * sizeof ( int ) );
    hspaux->minMAPQ = ini_params.Ini_minMAPQ;
    hspaux->maxMAPQ = ini_params.Ini_maxMAPQ;
    hspaux->bwaLikeScore = ini_params.Ini_bwaLikeScore;
    hspaux->maxLenReadName = ini_params.Ini_maxReadNameLen;

    if ( hspaux->bwaLikeScore )
    { bwase_initialize ( hspaux->g_log_n ); }

    // For storing the results of the unpaired reads
    hspaux->soap3AnsArray = ( ReadInputForDP ** ) malloc ( roundUp * sizeof ( ReadInputForDP * ) );
    hspaux->sa_start = ( unsigned int * ) malloc ( roundUp * sizeof ( unsigned int ) );
    hspaux->occ_start = ( unsigned int * ) malloc ( roundUp * sizeof ( unsigned int ) );
    hspaux->sa_num = ( unsigned int * ) malloc ( roundUp * sizeof ( unsigned int ) );
    hspaux->occ_num = ( unsigned int * ) malloc ( roundUp * sizeof ( unsigned int ) );

    // print MD string and NM tag?
    hspaux->isPrintMDNM = input_options.isPrintMDNM;

    // For the SAM output information
    hspaux->readGroup = input_options.readGroup;

    if ( strlen ( hspaux->readGroup ) == 0 )
    { hspaux->readGroup = input_options.queryFileName; }

    hspaux->sampleName = input_options.sampleName;

    if ( strlen ( hspaux->sampleName ) == 0 )
    { hspaux->sampleName = DEFAULT_SAMPLE_NAME; }

    hspaux->readGrpOption = input_options.readGrpOption;
    // ==================================================================
    // | For construction of arrays to store the unaligned read IDs
    // | for proceeding single deep-dp on them
    // ==================================================================
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS
    hspaux->readsIDForSingleDP = ( BothUnalignedPairsArrays * ) constructBothUnalignedPairsArrays ( 1 );
    hspaux->allHits = ( AllHits * ) constructAllHits (); // for storing the corresponding algnments
#endif
    // ======================================================================================
    // | ALLOCATE MEMORY FOR THE ARRAYS                                                     |
    // ======================================================================================
    queries0 = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) ); // a large array to store all queries
    readLengths0 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    readIDs0 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    upkdQualities0 = ( char * ) malloc ( roundUp * maxReadLength * sizeof ( char ) );
    memset ( upkdQualities0, 0, roundUp * maxReadLength );
    upkdQueryNames0 = ( char * ) malloc ( roundUp * ini_params.Ini_maxReadNameLen * sizeof ( char ) );
    queries1 = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) ); // a large array to store all queries
    readLengths1 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    readIDs1 = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    upkdQualities1 = ( char * ) malloc ( roundUp * maxReadLength * sizeof ( char ) );
    memset ( upkdQualities1, 0, roundUp * maxReadLength );
    upkdQueryNames1 = ( char * ) malloc ( roundUp * ini_params.Ini_maxReadNameLen * sizeof ( char ) );
    /*
    queries = (uint*) malloc(totalQueryLength * sizeof(uint)); // a large array to store all queries
    readLengths = (uint*) malloc(roundUp * sizeof(uint));
    readIDs = (uint*) malloc(roundUp * sizeof(uint));
    upkdQualities = (char*) malloc(roundUp*maxReadLength * sizeof(char));
    memset(upkdQualities,0,roundUp*maxReadLength);
    upkdQueryNames = (char*) malloc(roundUp*ini_params.Ini_maxReadNameLen * sizeof(char));
    */
    unAlignedPair = ( unsigned int * ) malloc ( roundUp * sizeof ( unsigned int ) );

    for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        readInputForDPall.inputArrays[threadId] = constructReadInputForDP ( ini_params.Ini_NumOfCpuThreads );
        readInputForNewDPall.inputArrays[threadId] = constructReadInputForDP ( ini_params.Ini_NumOfCpuThreads );
        readAnsForNonValidInsert.inputArrays[threadId] = constructReadInputForDP ( ini_params.Ini_NumOfCpuThreads );
    }

    uint maxBatchSize = NUM_BLOCKS * THREADS_PER_BLOCK * QUERIES_PER_THREAD; // queries processed in one kernel call
    // maxBatchSize has to be divible by 2
    maxBatchSize = maxBatchSize / 2 * 2;
    // #ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    // totalDeviceTime=0;
    // #endif
    // initialize the output files
    bam_header_t samOutputHeader;
    samfile_t * samOutputFilePtr[MAX_NUM_CPU_THREADS];
    samfile_t * samOutputDPFilePtr;
    samfile_t * samOutputUnpairFilePtr;

    switch ( input_options.outputFormat )
    {
        case SRA_OUTPUT_FORMAT_SAM_API:
            SAMOutputHeaderConstruct ( &samOutputHeader, hsp, hspaux, maxReadLength );

            for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
            {
                if ( input_options.isOutputBinary == 1 )
                { samOutputFilePtr[i] = samopen ( outputFileName[i], "wb", &samOutputHeader ); }
                else
                { samOutputFilePtr[i] = samopen ( outputFileName[i], "wh", &samOutputHeader ); }

                if ( samOutputFilePtr[i] == NULL )
                {
                    fprintf ( stderr, "Could not open the output file %s\n", outputFileName[i] );
                    exit ( 1 );
                }
            }

            break;

        default:
            for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
            {
                outputFile = ( FILE * ) fopen ( outputFileName[i], "w" );

                if ( outputFile == NULL ) { fprintf ( stderr, "Could not open the output file %s\n", outputFileName[i] ); exit ( 1 );}

                OCCWriteOutputHeader ( hsp, outputFile, maxReadLength, 1, input_options.outputFormat ); // will modify the number of reads later
                fclose ( outputFile );
            }

            break;
    }

    // For DP
    switch ( input_options.outputFormat )
    {
        case SRA_OUTPUT_FORMAT_SAM_API:
            if ( input_options.isOutputBinary == 1 )
            { samOutputDPFilePtr = samopen ( outputDPFileName, "wb", &samOutputHeader ); }
            else
            { samOutputDPFilePtr = samopen ( outputDPFileName, "wh", &samOutputHeader ); }

            if ( samOutputDPFilePtr == NULL )
            {
                fprintf ( stderr, "Could not open the output file %s\n", outputDPFileName );
                exit ( 1 );
            }

            break;

        default:
            outputDPFile = ( FILE * ) fopen ( outputDPFileName, "w" );

            if ( outputDPFile == NULL ) { fprintf ( stderr, "Could not open the output file %s\n", outputDPFileName ); exit ( 1 );}

            OCCWriteOutputHeader ( hsp, outputDPFile, maxReadLength, 1, input_options.outputFormat ); // will modify the number of reads later
            fclose ( outputDPFile );
            break;
    }

    // For Unpaired reads
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS

    if ( input_options.readType == PAIR_END_READ )
    {
        switch ( input_options.outputFormat )
        {
            case SRA_OUTPUT_FORMAT_SAM_API:
                if ( input_options.isOutputBinary == 1 )
                { samOutputUnpairFilePtr = samopen ( outputUnpairFileName, "wb", &samOutputHeader ); }
                else
                { samOutputUnpairFilePtr = samopen ( outputUnpairFileName, "wh", &samOutputHeader ); }

                if ( samOutputUnpairFilePtr == NULL )
                {
                    fprintf ( stderr, "Could not open the output file %s\n", outputUnpairFileName );
                    exit ( 1 );
                }

                break;
        }
    }

#endif
    // ======================================================================================
    // | LOADING INPUT SHORT READ FILE                                                      |
    // ======================================================================================
    size_t bufferSize;
    uint bufferIndex;
    char queryChar;
    size_t bufferSize2;
    uint bufferIndex2;
    char queryChar2;

    if ( input_options.isReadBAM )
    {
        // the query file is in BAM format
        bamQueryFile = bam_open ( queryFileName, "r" );
        bamHeader = bam_header_init ();
        bamHeader = bam_header_read ( bamQueryFile );
        bam = bam_init1 ();
    }
    else
    {
        //gzQueryFile = ( gzFile * ) gzopen ( queryFileName, "r" );
        gzQueryFile = gzopen ( queryFileName, "r" );

        if ( gzQueryFile == NULL ) { fprintf ( stderr, "Cannot open queryFile\n" ); exit ( 1 );}

        bufferSize = gzread ( gzQueryFile, queryFileBuffer, INPUT_BUFFER_SIZE );

        if ( bufferSize < INPUT_BUFFER_SIZE && ( !gzeof ( gzQueryFile ) ) )
        {
            const char * error_string;
            int err;
            error_string = gzerror ( gzQueryFile, & err );

            if ( err )
            {
                fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                exit ( EXIT_FAILURE );
            }
        }

        bufferIndex = 0;
        queryChar = queryFileBuffer[bufferIndex++];

        if ( input_options.readType == PAIR_END_READ )
        {
            // pair-ended reads
            //gzQueryFile2 = ( gzFile * ) gzopen ( queryFileName2, "r" );
            gzQueryFile2 = gzopen ( queryFileName2, "r" );

            if ( gzQueryFile2 == NULL ) { fprintf ( stderr, "Cannot open queryFile2\n" ); exit ( 1 );}

            bufferSize2 = gzread ( gzQueryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );

            if ( bufferSize2 < INPUT_BUFFER_SIZE && ( !gzeof ( gzQueryFile2 ) ) )
            {
                const char * error_string;
                int err;
                error_string = gzerror ( gzQueryFile2, & err );

                if ( err )
                {
                    fprintf ( stderr, "Error in reading the read file: %s.\n", error_string );
                    exit ( EXIT_FAILURE );
                }
            }

            bufferIndex2 = 0;
            queryChar2 = queryFileBuffer2[bufferIndex2++];
        }
    }

    uint accumReadNum = 0;
    uint numQueries;
    // create buffers
    InputReadsBuffer * buffer0 = InputReadsBufferFullCreate ( maxReadLength,
                                 maxNumQueries, wordPerQuery, quality_constant, queries0,
                                 readLengths0, readIDs0, upkdQualities0, upkdQueryNames0,
                                 isFastq, ini_params.Ini_maxReadNameLen );
    InputReadsBuffer * buffer1 = InputReadsBufferFullCreate ( maxReadLength,
                                 maxNumQueries, wordPerQuery, quality_constant, queries1,
                                 readLengths1, readIDs1, upkdQualities1, upkdQueryNames1,
                                 isFastq, ini_params.Ini_maxReadNameLen );
    AIOInputBuffer * aiob = AIOInputBufferCreate ( buffer0, buffer1 );
    InputFilePointers * ifp = InputFilePointersCreate ();

    if ( input_options.isReadBAM )
    {
        InputFilePointersSetBam ( ifp, bamQueryFile, bamHeader, bam );
    }
    else if ( input_options.readType == SINGLE_READ )
    {
        InputFilePointersSetSingle ( ifp, gzQueryFile );
    }
    else
    {
        InputFilePointersSetPair ( ifp, gzQueryFile, gzQueryFile2 );
    }

    aiob->reads = ifp;
    InputReadsBuffer * readyReadsBuffer;
    //create io thread
    AIOInputThreadCreate ( aiob, bufferSize, bufferIndex, index->charMap, queryChar,
                           queryFileBuffer, bufferSize2, bufferIndex2, queryChar2,
                           queryFileBuffer2 );
    //load reads
    readyReadsBuffer = LoadReadsFromAIOBuffer ( aiob );
    queries = readyReadsBuffer->queries;
    readLengths = readyReadsBuffer->readLengths;
    readIDs = readyReadsBuffer->readIDs;
    upkdQualities = readyReadsBuffer->upkdQualities;
    upkdQueryNames = readyReadsBuffer->upkdQueryNames;
    isFastq = readyReadsBuffer->isFastq;
    hspaux->isFastq = isFastq;
    numQueries = readyReadsBuffer->filledNum;

    while ( numQueries > 0 )
    {
        printf ( "[Main] Loaded %u short reads from the query file.\n", numQueries );
        readLoadTime = getElapsedTime ( startTime );
        printf ( "[Main] Elapsed time on host : %9.4f seconds\n\n", readLoadTime - lastEventTime );
        totalReadLoadTime += readLoadTime - lastEventTime;
        lastEventTime = readLoadTime;

        numOfAnswer = 0;
        numOfAlignedRead = 0;
        numOfUnAlignedPairs = 0;
        uint origNumQueries = numQueries;

        if ( detected_read_length == 0 )
        {
            // printParameters(input_options, ini_params);

            // ==================================================================
            // | DETECT THE READ LENGTH                                         |
            // ==================================================================
            if ( input_options.readType == PAIR_END_READ )
            {
                detected_read_length = GetReadLength ( readLengths, numQueries, 2 );
                detected_read_length2 = GetReadLength ( readLengths + 1, numQueries, 2 );

                // the minimum insert size cannot be smaller than detected_read_length2
                if ( input_options.insert_low < detected_read_length2 )
                { input_options.insert_low = detected_read_length2; }
            }
            else
            {
                detected_read_length = GetReadLength ( readLengths, numQueries, 1 );
            }

            // ==================================================================
            // | FOR DP IS ENABLED                                              |
            // | IF READ LENGTH < MIN_READ_LEN_FOR_DP (i.e. 30)                 |
            // |    THEN DP IS DISABLE.,                                        |
            // | IF READ LENGTH > 150, SKIP SOAP3 MODULE.                       |
            // | IF READ LENGTH <= 50,  ONLY ALLOW 1 MISMATCH IN SOAP3          |
            // ==================================================================
            if ( input_options.readType == SINGLE_READ && detected_read_length < MIN_READ_LEN_FOR_DP && input_options.enableDP == 1 )
            {
                input_options.enableDP = 0;
                printf ( "Dynamic programming is disabled because read length < %i\n", MIN_READ_LEN_FOR_DP );
            }
            else if ( input_options.readType == PAIR_END_READ && ( detected_read_length < MIN_READ_LEN_FOR_DP || detected_read_length2 < MIN_READ_LEN_FOR_DP ) && input_options.enableDP == 1 )
            {
                input_options.enableDP = 0;
                printf ( "Dynamic programming is disabled because read length < %i\n", MIN_READ_LEN_FOR_DP );
            }
            else if ( input_options.readType == PAIR_END_READ && ( detected_read_length > 150 || detected_read_length2 > 150 ) && input_options.enableDP == 1 )
            {
                ini_params.Ini_skipSOAP3Alignment = 1;
                printf ( "All reads are directly processed by DP, because the read length > 150\n" );
            }
            else if ( input_options.readType == PAIR_END_READ && ( detected_read_length <= 50 || detected_read_length2 <= 50 ) && input_options.enableDP == 1 )
            {
                input_options.numMismatch = 1;
            }

            // ==================================================================
            // | IF DP IS DISABLE AND USER DOES NOT SPECIFY # OF MISMATCHES     |
            // | THEN SET THE DEFAULT # OF MISMATCHES AS:                       |
            // |   - IF READ LENGTH < 50, DEFAULT_NUM_MISMATCH_FOR_SHORT_READ   |
            // |   - IF READ LENGTH >= 50, DEFAULT_NUM_MISMATCH_FOR_NORMAL_READ |
            // ==================================================================
            if ( input_options.enableDP == 0 && input_options.numMismatch == -1 )
            {
                // user does not specify # of mismatches
                input_options.numMismatch = getDefaultMismatchNum ( detected_read_length );
                printf ( "Maximum number of mismatches allowed: %i\n",  input_options.numMismatch );
            }

            // get the max hit # for default DP
            if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ )
            {
                input_options.maxHitNum = getMaxHitNumForDefaultDP ( detected_read_length );
                input_options.maxHitNum2 = getMaxHitNumForDefaultDP ( detected_read_length2 );
            }
        }

        // Reset the array for mapping quality score calculation
        memset ( hspaux->x0_array, 0, roundUp * sizeof ( int ) );
        memset ( hspaux->x1_array, 0, roundUp * sizeof ( int ) );
        memset ( hspaux->mismatch_array, 0, roundUp * sizeof ( int ) );
        memset ( hspaux->sa_start, 0, roundUp * sizeof ( unsigned int ) );
        memset ( hspaux->occ_start, 0, roundUp * sizeof ( unsigned int ) );
        memset ( hspaux->sa_num, 0, roundUp * sizeof ( unsigned int ) );
        memset ( hspaux->occ_num, 0, roundUp * sizeof ( unsigned int ) );
        memset ( ( ( ReadInputForDP ** ) hspaux->soap3AnsArray ), 0, roundUp * sizeof ( ReadInputForDP * ) );
        // Reset the array for storing the single alignment for those unaligned paired-end reads
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS
        resetAllHits ( ( AllHits * ) hspaux->allHits );
#endif

        // =====================================================================
        // | Reset the arrays for storing alignment results for semi-global DP |
        // =====================================================================

        if ( input_options.enableDP == 1 )
        {
            if ( input_options.readType == PAIR_END_READ )
                for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
                {
                    resetReadInputForDP ( readInputForDPall.inputArrays[threadId] );
                    resetReadInputForDP ( readInputForNewDPall.inputArrays[threadId] );
                    resetReadInputForDP ( readAnsForNonValidInsert.inputArrays[threadId] );
                }

            resetBothUnalignedPairsArrays ( bothUnalignedPairsArrays );
        }

        // =======================================
        // | BATCH PROCESS READ (GPU/GPU/CPU)    |
        // =======================================
        char ** currOutputFileName = outputFileName;
        samfile_t ** currSamOutputFilePtr = samOutputFilePtr;

        // ========================================
        // | IF THE INDEX IS NOT IN DEVICE, THEN  |
        // | COPY INDEX TO DEVICE MEMORY          |
        // ========================================

        if ( indexLoadedToGPU == 0 )
        {
            GPUINDEXUpload ( index, &_bwt, &_occ,
                             &_revBwt, &_revOcc );
            double copyTime = getElapsedTime ( startTime );
            printf ( "[Main] Finished copying index into device (GPU).\n" );
            printf ( "[Main] Loading time : %9.4f seconds\n\n", copyTime - lastEventTime );
            lastEventTime = copyTime;
            indexLoadedToGPU = 1;
        }

        //*******************************//
        // Perform Alignment             //
        //*******************************//

        if ( input_options.readType == PAIR_END_READ )
        {
            soap3_dp_pair_align ( queries, readLengths, input_options.numMismatch, wordPerQuery,
                                  maxBatchSize, numQueries, accumReadNum,
                                  index,
                                  _bwt, _revBwt,
                                  _occ, _revOcc,
                                  ini_params, input_options,
                                  maxReadLength, detected_read_length, detected_read_length2,
                                  upkdQualities,
                                  unAlignedPair, numOfUnAlignedPairs,
                                  readIDs, upkdQueryNames,
                                  currOutputFileName, currSamOutputFilePtr,
                                  outputDPFileName, samOutputDPFilePtr,
                                  samOutputUnpairFilePtr,
                                  numOfAnswer, numOfAlignedRead,
                                  &readInputForDPall, &readInputForNewDPall,
                                  &readAnsForNonValidInsert,
                                  bothUnalignedPairsArrays,
                                  startTime, lastEventTime, totalAlignmentTime,
                                  indexLoadedToGPU );
        }
        else
        {
            soap3_dp_single_align ( queries, readLengths, input_options.numMismatch, wordPerQuery,
                                    maxBatchSize, numQueries, accumReadNum,
                                    index, _bwt, _revBwt,
                                    _occ, _revOcc,
                                    ini_params, input_options,
                                    maxReadLength, detected_read_length,
                                    upkdQualities,
                                    unAlignedPair, numOfUnAlignedPairs,
                                    readIDs, upkdQueryNames,
                                    currOutputFileName, currSamOutputFilePtr,
                                    outputDPFileName, samOutputDPFilePtr,
                                    numOfAnswer, numOfAlignedRead,
                                    &readInputForDPall,
                                    unalignedSinglesArrays,
                                    startTime, lastEventTime, totalAlignmentTime,
                                    indexLoadedToGPU );
        }

        // ========================================
        // | GET NEXT BATCH OF READ               |
        // ========================================
        totalReadsAlignedForInputReads += numOfAlignedRead;
        totalAnsForInputReads += numOfAnswer;
        accumReadNum += origNumQueries;
        ResetBufferStatusToUnfilled ( aiob );
        readyReadsBuffer = LoadReadsFromAIOBuffer ( aiob );
        queries = readyReadsBuffer->queries;
        readLengths = readyReadsBuffer->readLengths;
        readIDs = readyReadsBuffer->readIDs;
        upkdQualities = readyReadsBuffer->upkdQualities;
        upkdQueryNames = readyReadsBuffer->upkdQueryNames;
        isFastq = readyReadsBuffer->isFastq;
        numQueries = readyReadsBuffer->filledNum;

        // If the current opened file still have queries returned
        // skip the code to process the result / open next file; and continue
        if ( numQueries > 0 )
        {
            // Skip
            continue;
        }

        // show the summary of the result and then load another pair of read files
        // ======================================================================================
        // | SHOW THE SUMMARY                                                                   |
        // ======================================================================================
        if ( input_options.readType == PAIR_END_READ )
        {
            // printf("[Main] Overall number of pairs of reads aligned: %llu (number of alignments: %llu)\n", totalReadsAlignedForInputReads/2, totalAnsForInputReads);
            printf ( "[Main] Overall number of pairs of reads aligned: %llu\n", totalReadsAlignedForInputReads / 2 );
        }
        else
        {
            // printf("[Main] Overall number of reads aligned: %llu (number of alignments: %llu)\n", totalReadsAlignedForInputReads, totalAnsForInputReads);
            // printf("[Main] Overall number of unaligned reads: %llu\n", accumReadNum-totalReadsAlignedForInputReads);
            printf ( "[Main] Overall number of reads aligned: %llu\n", totalReadsAlignedForInputReads );
            printf ( "[Main] Overall number of unaligned reads: %llu\n", accumReadNum - totalReadsAlignedForInputReads );
        }

        printf ( "[Main] Overall read load time : %9.4f seconds\n", totalReadLoadTime );
        printf ( "[Main] Overall alignment time (excl. read loading) : %9.4f seconds\n", totalAlignmentTime + totalTrimAlignmentTime );
        // ======================================================================================
        // | SHOW THE COMMAND FOR MERGING THE OUTPUT FILES INTO ONE                             |
        // ======================================================================================

        // update the output files
        switch ( input_options.outputFormat )
        {
            case SRA_OUTPUT_FORMAT_DEFAULT:
                // update the header of the output files
                outputDPFile = ( FILE * ) fopen ( outputDPFileName, "r+" );

                if ( outputDPFile == NULL ) { fprintf ( stderr, "Cannot open outputFile %s\n", outputDPFileName ); exit ( 1 );}

                fseek ( outputDPFile , 0 , SEEK_SET );
                OCCWriteOutputHeader ( hsp, outputDPFile, maxReadLength, accumReadNum, input_options.outputFormat ); // update the number of reads
                fclose ( outputDPFile );
                break;

            case SRA_OUTPUT_FORMAT_SAM_API:
                samclose ( samOutputDPFilePtr );
                break;
        }

        switch ( input_options.outputFormat )
        {
            case SRA_OUTPUT_FORMAT_DEFAULT:

                // update the header of the output files
                for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
                {
                    outputFile = ( FILE * ) fopen ( outputFileName[i], "r+" );

                    if ( outputFile == NULL ) { fprintf ( stderr, "Cannot open outputFile %s\n", outputFileName[i] ); exit ( 1 );}

                    fseek ( outputFile , 0 , SEEK_SET );
                    OCCWriteOutputHeader ( hsp, outputFile, maxReadLength, accumReadNum, input_options.outputFormat ); // update the number of reads
                    fclose ( outputFile );
                }

                break;

            case SRA_OUTPUT_FORMAT_SAM_API:
                for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
                {
                    samclose ( samOutputFilePtr[i] );
                }

                break;
        }

        // For Unpaired reads
        if ( input_options.readType == PAIR_END_READ )
        {
            switch ( input_options.outputFormat )
            {
                case SRA_OUTPUT_FORMAT_SAM_API:
                    samclose ( samOutputUnpairFilePtr );
                    break;
            }
        }

        // create the done file
        outputDoneFile = ( FILE * ) fopen ( outputDoneFileName, "w" );

        if ( outputDoneFile == NULL )
        {
            fprintf ( stderr, "Could not create the output file %s\n",
                      outputDoneFileName );
        }

        fclose ( outputDoneFile );
        free ( outputDoneFileName );

        // release memory and close files
        for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
        { free ( outputFileName[i] ); }

        free ( outputDPFileName );
        free ( outputUnpairFileName );

        if ( input_options.isReadBAM )
        {
            bam_close ( bamQueryFile );
            bam_destroy1 ( bam );
        }
        else
        {
            gzclose ( gzQueryFile );

            if ( input_options.readType == PAIR_END_READ )
            {
                gzclose ( gzQueryFile2 );
            }
        }

        // ======================================================================================
        // | Load the next set of files                                                         |
        // ======================================================================================

        if ( input_options.isReadList == 1 && currExp < expNum )
        {
            printf ( "\n[Main] Load the next set of files\n" );
            updateInputOption ( ( &input_options ), multiInput, currExp++ );
            queryFileName = input_options.queryFileName;

            if ( input_options.readType == PAIR_END_READ && input_options.isReadBAM == 0 )
            {
                queryFileName2 = input_options.queryFileName2;
            }

            // Update the HSP information
            hspaux->readGroup = input_options.readGroup;

            if ( strlen ( hspaux->readGroup ) == 0 )
            { hspaux->readGroup = input_options.queryFileName; }

            hspaux->sampleName = input_options.sampleName;

            if ( strlen ( hspaux->sampleName ) == 0 )
            { hspaux->sampleName = DEFAULT_SAMPLE_NAME; }

            hspaux->readGrpOption = input_options.readGrpOption;

            // update the output files
            for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
            { outputFileName[i] = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 10 ); }

            outputDPFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 9 );
            outputUnpairFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 8 );

            for ( i = 0; i < MAX_NUM_CPU_THREADS; i++ )
            {
                sprintf ( outputFileName[i], "%s.gout.%d", input_options.outputPrefix, i + 1 );
            }

            sprintf ( outputDPFileName, "%s.dpout.1", input_options.outputPrefix );
            sprintf ( outputUnpairFileName, "%s.unpair", input_options.outputPrefix );
            outputDoneFileName = ( char * ) malloc ( strlen ( input_options.outputPrefix ) + 6 );
            sprintf ( outputDoneFileName, "%s.done", input_options.outputPrefix );

            // initialize the output files
            switch ( input_options.outputFormat )
            {
                case SRA_OUTPUT_FORMAT_SAM_API:
                    SAMOutputHeaderConstruct ( &samOutputHeader, hsp, hspaux, maxReadLength );

                    for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
                    {
                        if ( input_options.isOutputBinary == 1 )
                        { samOutputFilePtr[i] = samopen ( outputFileName[i], "wb", &samOutputHeader ); }
                        else
                        { samOutputFilePtr[i] = samopen ( outputFileName[i], "wh", &samOutputHeader ); }

                        if ( samOutputFilePtr[i] == NULL )
                        {
                            fprintf ( stderr, "Could not open the output file %s\n", outputFileName[i] );
                            exit ( 1 );
                        }
                    }

                    break;

                default:
                    for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
                    {
                        outputFile = ( FILE * ) fopen ( outputFileName[i], "w" );

                        if ( outputFile == NULL ) { fprintf ( stderr, "Cannot open outputFile %s\n", outputFileName[i] ); exit ( 1 );}

                        OCCWriteOutputHeader ( hsp, outputFile, maxReadLength, 1, input_options.outputFormat ); // will modify the number of reads later
                        fclose ( outputFile );
                    }

                    break;
            }

            // For DP
            switch ( input_options.outputFormat )
            {
                case SRA_OUTPUT_FORMAT_SAM_API:
                    if ( input_options.isOutputBinary == 1 )
                    { samOutputDPFilePtr = samopen ( outputDPFileName, "wb", &samOutputHeader ); }
                    else
                    { samOutputDPFilePtr = samopen ( outputDPFileName, "wh", &samOutputHeader ); }

                    if ( samOutputDPFilePtr == NULL )
                    {
                        fprintf ( stderr, "Could not open the output file %s\n", outputDPFileName );
                        exit ( 1 );
                    }

                    break;

                default:
                    outputDPFile = ( FILE * ) fopen ( outputDPFileName, "w" );

                    if ( outputDPFile == NULL ) { fprintf ( stderr, "Cannot open outputFile %s\n", outputDPFileName ); exit ( 1 );}

                    OCCWriteOutputHeader ( hsp, outputDPFile, maxReadLength, 1, input_options.outputFormat ); // will modify the number of reads later
                    fclose ( outputDPFile );
                    break;
            }

            // For Unpaired reads
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS

            if ( input_options.readType == PAIR_END_READ )
            {
                switch ( input_options.outputFormat )
                {
                    case SRA_OUTPUT_FORMAT_SAM_API:
                        if ( input_options.isOutputBinary == 1 )
                        { samOutputUnpairFilePtr = samopen ( outputUnpairFileName, "wb", &samOutputHeader ); }
                        else
                        { samOutputUnpairFilePtr = samopen ( outputUnpairFileName, "wh", &samOutputHeader ); }

                        if ( samOutputUnpairFilePtr == NULL )
                        {
                            fprintf ( stderr, "Could not open the output file %s\n", outputUnpairFileName );
                            exit ( 1 );
                        }

                        break;
                }
            }

#endif
            // reset the variables
            accumReadNum = 0;
            totalReadsAlignedForInputReads = 0;
            totalAnsForInputReads = 0;
            detected_read_length = 0; // the max read length for the first ten reads
            totalAlignmentTime = 0;
            totalReadLoadTime = 0;

            // load reads
            if ( input_options.readType == SINGLE_READ )
            { printf ( "\n\n[Main] Loading read file %s\n", queryFileName ); }
            else if ( input_options.isReadBAM == 0 )
            { printf ( "\n\n[Main] Loading read files %s and %s\n", queryFileName, queryFileName2 ); }

            if ( input_options.isReadBAM )
            {
                // the query file is in BAM format
                bamQueryFile = bam_open ( queryFileName, "r" );
                bamHeader = bam_header_init ();
                bamHeader = bam_header_read ( bamQueryFile );
                bam = bam_init1 ();
            }
            else
            {
                //gzQueryFile = ( gzFile * ) gzopen ( queryFileName, "r" );
                gzQueryFile = gzopen ( queryFileName, "r" );

                if ( gzQueryFile == NULL ) { fprintf ( stderr, "Cannot open queryFile\n" ); exit ( 1 );}

                bufferSize = gzread ( gzQueryFile, queryFileBuffer, INPUT_BUFFER_SIZE );
                bufferIndex = 0;
                queryChar = queryFileBuffer[bufferIndex++];

                if ( input_options.readType == PAIR_END_READ )
                {
                    // pair-ended reads
                    //gzQueryFile2 = ( gzFile * ) gzopen ( queryFileName2, "r" );
                    gzQueryFile2 = gzopen ( queryFileName2, "r" );

                    if ( gzQueryFile2 == NULL ) { fprintf ( stderr, "Cannot open queryFile2\n" ); exit ( 1 );}

                    bufferSize2 = gzread ( gzQueryFile2, queryFileBuffer2, INPUT_BUFFER_SIZE );
                    bufferIndex2 = 0;
                    queryChar2 = queryFileBuffer2[bufferIndex2++];
                }
            }

            if ( input_options.isReadBAM )
            {
                InputFilePointersSetBam ( ifp, bamQueryFile, bamHeader, bam );
            }
            else if ( input_options.readType == SINGLE_READ )
            {
                InputFilePointersSetSingle ( ifp, gzQueryFile );
            }
            else
            {
                InputFilePointersSetPair ( ifp, gzQueryFile, gzQueryFile2 );
            }

            aiob->reads = ifp;
            //clear buffer' status
            AIOInputBufferClear ( aiob );
            //create io thread
            AIOInputThreadCreate ( aiob, bufferSize, bufferIndex, index->charMap, queryChar,
                                   queryFileBuffer, bufferSize2, bufferIndex2, queryChar2,
                                   queryFileBuffer2 );
            //load reads
            readyReadsBuffer = LoadReadsFromAIOBuffer ( aiob );
            queries = readyReadsBuffer->queries;
            readLengths = readyReadsBuffer->readLengths;
            readIDs = readyReadsBuffer->readIDs;
            upkdQualities = readyReadsBuffer->upkdQualities;
            upkdQueryNames = readyReadsBuffer->upkdQueryNames;
            isFastq = readyReadsBuffer->isFastq;
            numQueries = readyReadsBuffer->filledNum;

            // Update the isFastq variable inside HSP
            hspaux->isFastq = isFastq;
        }
    }

    cudaThreadExit ();
    // ======================================================================================
    // | CLEAN UP                                                                           |
    // ======================================================================================
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS
    freeBothUnalignedPairsArrays ( ( BothUnalignedPairsArrays * ) hspaux->readsIDForSingleDP );
    releaseAllHits ( ( AllHits * ) hspaux->allHits );
#endif
    free ( hspaux->x0_array );
    free ( hspaux->x1_array );
    free ( hspaux->mismatch_array );
    free ( hspaux->sa_start );
    free ( hspaux->occ_start );
    free ( hspaux->sa_num );
    free ( hspaux->occ_num );
    ReadInputForDP ** array = ( ReadInputForDP ** ) hspaux->soap3AnsArray;
    free ( array );

    // free device memory
    if ( indexLoadedToGPU == 1 )
    {
        printf ( "[Main] Free device memory..\n" );
        GPUINDEXFree ( _bwt, _occ, _revBwt, _revOcc );
    }

    printf ( "[Main] Free index from host memory..\n" );
    INDEXFree ( index, ini_params.Ini_shareIndex );
    printf ( "[Main] Free host memory..\n" );
    free ( queries0 );
    free ( readLengths0 );
    free ( readIDs0 );
    free ( upkdQualities0 );
    free ( upkdQueryNames0 );
    free ( queries1 );
    free ( readLengths1 );
    free ( readIDs1 );
    free ( upkdQualities1 );
    free ( upkdQueryNames1 );
    AIOInputBufferFree ( aiob );
    // free(queries);
    // free(readLengths);
    // free(readIDs);
    free ( unAlignedPair );
    // free(isBad);
    // free(upkdQueries);
    // free(upkdQualities);
    // free(upkdQueryNames);

    for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        freeReadInputForDP ( readInputForDPall.inputArrays[threadId] );
        freeReadInputForDP ( readInputForNewDPall.inputArrays[threadId] );
        freeReadInputForDP ( readAnsForNonValidInsert.inputArrays[threadId] );
    }

    free ( readInputForDPall.inputArrays ); // for half-aligned pairs
    free ( readInputForNewDPall.inputArrays ); // for half-aligned pairs
    free ( readAnsForNonValidInsert.inputArrays );
    freeBothUnalignedPairsArrays ( bothUnalignedPairsArrays ); // for both-unaligned pairs OR single-unaligned reads

    if ( multiInput != NULL )
    {
        free ( multiInput );
    }

    if ( input_options.outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    {
        SAMOutputHeaderDestruct ( &samOutputHeader );
    }

    return 0;
}
