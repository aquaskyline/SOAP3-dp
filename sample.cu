/*
 *  sample.cpp
 *  soap3_ext_r5
 *
 *  Created by kfwong on 14/3/12.
 *  Copyright 2012 HKU. All rights reserved.
 *
 */


#include "soap3-dp-module.h"


// This function is to load the single reads for at most "maxNumQueries" # of reads
int loadSingleR ( FILE * queryFile, char * queryFileBuffer, uint * queries, uint * readLengths, uint * readIDs,
                  uint maxReadLength, uint maxNumQueries, size_t & bufferSize, char & queryChar,
                  uint & bufferIndex, uint wordPerQuery )
{
    // for single reads in FASTA format
    // load reads from the read file
    // return how many reads are loaded
    // Convert ACGT to 0123
    uint charMap[256];

    for ( int i = 0; i < 256; i++ )
    { charMap[i] = 0; }

    charMap['A'] = 0;
    charMap['C'] = 1;
    charMap['G'] = 2;
    charMap['T'] = 3;
    charMap['U'] = 3;
    charMap['N'] = 2; // N -> G
    charMap['a'] = 0;
    charMap['c'] = 1;
    charMap['g'] = 2;
    charMap['t'] = 3;
    charMap['u'] = 3;
    charMap['n'] = 2; // N -> G
    ullint queriesRead = 0;
    uint currentWord = 0;
    uint offset = 0;
    uint bits = 0;
    uint * queryPtr = queries;

    while ( bufferSize != 0 )
    {
        //Read everything before the entry point of a read, character ">"
        while ( bufferSize != 0 && queryChar != '>' )
        {
            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }
        }

        //Read the header of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = fread ( queryFileBuffer, sizeof ( char ), INPUT_BUFFER_SIZE, queryFile );
            bufferIndex = 0;
        }

        while ( bufferSize != 0  && queryChar != '\n' )
        {
            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }
        }

        //Read the pattern body of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = fread ( queryFileBuffer, sizeof ( char ), INPUT_BUFFER_SIZE, queryFile );
            bufferIndex = 0;
        }

        uint nucleoId = 0;

        while ( bufferSize != 0  && queryChar != '>' && queryChar != '@' && queryChar != '+' )
        {
            if ( queryChar != '\n' )
            {
                bits = charMap[queryChar];

                if ( nucleoId < maxReadLength )
                {
                    currentWord |= ( bits << ( offset * BIT_PER_CHAR ) );
                    offset++;

                    if ( offset == CHAR_PER_WORD )
                    {
                        *queryPtr = currentWord;
                        queryPtr += 32;
                        offset = 0;
                        currentWord = 0;
                    }
                }

                nucleoId++;
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }
        }

        if ( offset > 0 )
        { *queryPtr = currentWord; }

        currentWord = 0;
        offset = 0;
        readLengths[queriesRead] = nucleoId;
        readIDs[queriesRead] = queriesRead + 1;

        if ( nucleoId > maxReadLength )
        {
            // printf("[WARNING] Read #%u is longer than %u! Read truncated.\n", queriesRead+1+accumReadNum, maxReadLength);
            readLengths[queriesRead] = maxReadLength;
        }

        queriesRead++;
        queryPtr = queries + ( queriesRead / 32 * 32 * wordPerQuery + queriesRead % 32 );

        if ( queriesRead >= maxNumQueries )
        { break; }
    }

    return queriesRead;
}


int main ( int argc, char ** argv )
{
    char * indexFileName = argv[1]; // index file name
    char * queryFileName = argv[2]; // read file name
    int output_option = 1;         // 1: all valid; 2: all best; 3: unique best; 4: random best
    uint maxReadLength = 120;      // maximum read length
    uint maxNumQueries = NUM_BLOCKS * THREADS_PER_BLOCK * QUERIES_PER_THREAD;  // 1 million reads
    int numCPUThreads = 3;
    // parameters
    SingleAlignParam param;
    param.maxReadLength = 120;
    param.numMismatch = 2;
    param.outputOption = output_option;
    param.cpuNumThreads = numCPUThreads;
    param.maxHitNum = 1000;       // maximum number of hits reported for each read
    param.enableDP = 1;
    param.scoring.cutoffThreshold = 30;
    param.scoring.matchScore = 1;
    param.scoring.mismatchScore = -2;
    param.scoring.openGapScore = -3;
    param.scoring.extendGapScore = -1;
    // first construct the global array for storing the alignment results
    AlgnResultArrays * algnResultArrays = resultArraysConstruct ( numCPUThreads );
    // load index
    Soap3Index * index = INDEXLoad ( NULL, indexFileName, FALSE );
    // variables
    uint wordPerQuery = 1;

    while ( wordPerQuery < maxReadLength )
    { wordPerQuery *= 2; }

    wordPerQuery = wordPerQuery / CHAR_PER_WORD;
    ullint roundUp = ( maxNumQueries + 31 ) / 32 * 32;
    ullint totalQueryLength = roundUp * wordPerQuery;
    uint * queries = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) );
    uint * readLengths = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    uint * readIDs = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
    unsigned long long numOfAnswer = 0;
    unsigned int numOfAlignedRead = 0;
    // load read file
    char queryFileBuffer[INPUT_BUFFER_SIZE];
    FILE * queryFile = ( FILE * ) fopen ( queryFileName, "r" );

    if ( queryFile == NULL ) { fprintf ( stderr, "Cannot open queryFile\n" ); exit ( 1 );}

    size_t bufferSize = fread ( queryFileBuffer, sizeof ( char ), INPUT_BUFFER_SIZE, queryFile );
    uint bufferIndex = 0;
    char queryChar = queryFileBuffer[bufferIndex++];
    uint numQueries = loadSingleR ( queryFile, queryFileBuffer, queries, readLengths, readIDs,
                                    maxReadLength, maxNumQueries, bufferSize, queryChar,
                                    bufferIndex, wordPerQuery );
    // perform alignment
    alignSingleR ( queries, readLengths, readIDs,
                   wordPerQuery, numQueries, index,
                   &param, numOfAnswer, numOfAlignedRead, algnResultArrays );

    // print the alignments
    for ( int i = 0; i < algnResultArrays->numArrays; i++ )
    {
        AlgnResult * algnArray = algnResultArrays->algnArrays[i];

        for ( int j = 0; j < algnArray->occTotalNum; j++ )
        {
            occRec occ = algnArray->occ_list[j];
            printf ( "%u 1 %u %c %i\n", occ.readID, occ.ambPosition, ( occ.strand == 1 ? '+' : '-' ), ( int ) occ.score );
        }
    }

    // to reset the array for storing alignment results
    // resultArraysReset(algnResultArrays);
    // clean the memory
    free ( queries );
    free ( readLengths );
    free ( readIDs );
    INDEXFree ( index, FALSE );
    resultArraysFree ( algnResultArrays );
}
