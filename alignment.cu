/*
 *
 *    alignment.cu
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

#include "alignment.h"

// COPY INDEX TO DEVICE MEMORY
void GPUINDEXUpload ( Soap3Index * index, uint ** _bwt, uint ** _occ,
                      uint ** _revBwt, uint ** _revOcc )
{
    BWT * bwt = index->sraIndex->bwt;
    BWT * revBwt = index->sraIndex->rev_bwt;
    unsigned int * revOccValue = index->gpu_revOccValue;
    unsigned int * occValue = index->gpu_occValue;
    unsigned int numOfOccValue = index->gpu_numOfOccValue;
    cudaError_t gpuErr;
    gpuErr = cudaMalloc ( ( void ** ) _bwt, bwt->bwtSizeInWord * sizeof ( uint ) );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MALLOC FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMalloc ( ( void ** ) _occ, numOfOccValue * ALPHABET_SIZE * sizeof ( uint ) );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MALLOC FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMemcpy ( ( *_bwt ), bwt->bwtCode, bwt->bwtSizeInWord * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMemcpy ( ( *_occ ), occValue, numOfOccValue * ALPHABET_SIZE * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMalloc ( ( void ** ) _revBwt, revBwt->bwtSizeInWord * sizeof ( uint ) );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MALLOC FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMalloc ( ( void ** ) _revOcc, numOfOccValue * ALPHABET_SIZE * sizeof ( uint ) );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MALLOC FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMemcpy ( ( *_revBwt ), revBwt->bwtCode, revBwt->bwtSizeInWord * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMemcpy ( ( *_revOcc ), revOccValue, numOfOccValue * ALPHABET_SIZE * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMemcpyToSymbol ( gpuCharMap, index->charMap, sizeof ( unsigned char ) * 256 );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }
}

void GPUINDEXFree ( uint * _bwt, uint * _occ, uint * _revBwt, uint * _revOcc )
{
    cudaFree ( _bwt );
    cudaFree ( _occ );
    cudaFree ( _revBwt );
    cudaFree ( _revOcc );
}

// perform round1 alignment in GPU
void perform_round1_alignment ( uint * nextQuery, uint * nextReadLength, uint * answers[][MAX_NUM_CASES],
                                uint numMismatch, uint numCases, uint sa_range_allowed, uint wordPerQuery, uint word_per_ans,
                                bool isExactNumMismatch, int doubleBufferIdx, uint blocksNeeded, ullint batchSize,
                                Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc )
{
    cudaError_t gpuErr;
    uint * _queries, *_readLengths, *_answers;
    bool * _isBad;
    ullint roundUp = ( batchSize + 31 ) / 32 * 32;
    // allocated device memory for bad read indicator
    bool * isBad;
    isBad = ( bool * ) malloc ( roundUp * sizeof ( bool ) ); // an array to store bad read indicator
    memset ( isBad, 0, roundUp );
    gpuErr = cudaMalloc ( ( void ** ) &_isBad, roundUp * sizeof ( bool ) );
    BWT * bwt = index->sraIndex->bwt;
    BWT * revBwt = index->sraIndex->rev_bwt;

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MALLOC FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    // to initialize the array _isBad
    gpuErr = cudaMemcpy ( _isBad, isBad,
                          roundUp * sizeof ( bool ),
                          cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    // to initialize the array answers
    // for (uint caseno=0; caseno < numCases; caseno++) {
    //       memset(answers[doubleBufferIdx][caseno], 0, roundUp * word_per_ans * sizeof(uint));
    // }
    // printf("[perform_round1_alignment] sa_range_allowed = %u; word_per_ans = %u\n", sa_range_allowed, word_per_ans);
    // allocate device memory for queries and answers
    cudaMalloc ( ( void ** ) &_queries, roundUp * wordPerQuery * sizeof ( uint ) );
    cudaMalloc ( ( void ** ) &_readLengths, roundUp * sizeof ( uint ) );
    cudaMalloc ( ( void ** ) &_answers, roundUp * word_per_ans * sizeof ( uint ) );
    gpuErr = cudaMemcpy ( _queries, nextQuery, roundUp * wordPerQuery * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMemcpy ( _readLengths, nextReadLength, roundUp * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    for ( uint caseno = 0; caseno < numCases; caseno++ )
    {
        // =======================================
        // | GPU-1: FOR EACH CASE                |
        // =======================================
        if ( numMismatch <= 3 )
            kernel <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( caseno, _queries, _readLengths, batchSize, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0, bwt->textLength,
              _answers, _isBad, 0 , numMismatch, sa_range_allowed, word_per_ans, isExactNumMismatch );
        else if ( caseno < 5 ) // 4 mismatch and case 0 - 4
            kernel_4mismatch_1 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( caseno, _queries, _readLengths, batchSize, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0, bwt->textLength,
              _answers, _isBad, 0 , sa_range_allowed, word_per_ans, isExactNumMismatch );
        else  // 4 mismatch and case 5 - 9
            kernel_4mismatch_2 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( caseno, _queries, _readLengths, batchSize, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0, bwt->textLength,
              _answers, _isBad, 0 , sa_range_allowed, word_per_ans, isExactNumMismatch );

        gpuErr = cudaMemcpy ( answers[doubleBufferIdx][caseno], _answers, roundUp * word_per_ans * sizeof ( uint ), cudaMemcpyDeviceToHost );

        if ( gpuErr != cudaSuccess )
        {
            printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
            exit ( 1 );
        }
    }

    // free the memories
    free ( isBad );
    cudaFree ( _isBad );
    cudaFree ( _queries );
    cudaFree ( _readLengths );
    cudaFree ( _answers );
}



// perform round2 alignment in GPU
void perform_round2_alignment ( uint * queries, uint * readLengths, uint * answers[][MAX_NUM_CASES],
                                uint numMismatch, uint numCases, uint sa_range_allowed_2, uint wordPerQuery, uint word_per_ans, uint word_per_ans_2,
                                bool isExactNumMismatch, int doubleBufferIdx, uint blocksNeeded, ullint batchSize,
                                Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                                uint processedQuery, uint * badReadIndices[][MAX_NUM_CASES],
                                uint * badAnswers[][MAX_NUM_CASES] )
{
    uint * badQueries;
    uint * badReadLengths;
    BWT * bwt = index->sraIndex->bwt;
    BWT * revBwt = index->sraIndex->rev_bwt;

    for ( int whichCase = 0; whichCase < numCases; ++whichCase )
    {
        ullint numBads = 0;

        // Count number of bad reads and
        for ( ullint readId = 0; readId < batchSize; readId++ )
        {
            ullint srcOffset = ( ( readId ) / 32 * 32 * word_per_ans + readId % 32 );
            numBads += ( answers[doubleBufferIdx][whichCase][srcOffset] > 0xFFFFFFFD );
        }

        if ( numBads == 0 )
        { continue; }

        // Allocate memory and copy bad reads to another array
        ullint roundUp = ( numBads + 31 ) / 32 * 32;
        badQueries = ( uint * ) malloc ( roundUp * wordPerQuery * sizeof ( uint ) );
        // printf("size of badQueries = %u\n", roundUp * wordPerQuery);
        badReadLengths = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
        // printf("size of badReadLengths = %u\n", roundUp);
        numBads = 0;

        for ( ullint readId = 0; readId < batchSize; ++readId )
        {
            ullint srcOffset = ( ( readId ) / 32 * 32 * word_per_ans + readId % 32 );
            ullint srcQueryOffset = ( processedQuery + readId ) / 32 * 32 * wordPerQuery + readId % 32;

            if ( answers[doubleBufferIdx][whichCase][srcOffset] > 0xFFFFFFFD ) // is bad
            {
                badReadIndices[doubleBufferIdx][whichCase][numBads] = readId;
                ullint targetOffset = numBads / 32 * 32 * wordPerQuery + numBads % 32;

                for ( ullint i = 0; i < wordPerQuery; ++i ) // copy each word of the read
                {
                    badQueries[targetOffset + i * 32] = * ( queries + ( srcQueryOffset + i * 32 ) );
                }

                badReadLengths[numBads] = readLengths[processedQuery + readId];
                numBads++;
            }
        }

        // Copy reads to device memory, call kernel, and copy results back to host
        uint * _queries, *_readLengths, *_answers;
        bool * _isBad = NULL;
        // bool *isBad = (bool*) malloc((batchSize + 31) / 32 * 32 * sizeof(bool)); // an array to store bad read indicator
        // memset(isBad,0, (batchSize + 31) / 32 * 32 * sizeof(bool));
        // cudaMalloc((void**)&_isBad, (batchSize + 31) / 32 * 32 * sizeof(bool));
        // cudaMemcpy(_isBad, isBad,(batchSize + 31) / 32 * 32 * sizeof(bool), cudaMemcpyHostToDevice);
        // free(isBad);
        cudaMalloc ( ( void ** ) &_queries, roundUp * wordPerQuery * sizeof ( uint ) );
        cudaMemcpy ( _queries, badQueries,
                     roundUp * wordPerQuery * sizeof ( uint ), cudaMemcpyHostToDevice );
        free ( badQueries );
        cudaMalloc ( ( void ** ) &_readLengths, roundUp * sizeof ( uint ) );
        cudaMemcpy ( _readLengths, badReadLengths,
                     roundUp * sizeof ( uint ), cudaMemcpyHostToDevice );
        free ( badReadLengths );
        // reset the values of the array badAnswers
        // memset(badAnswers[doubleBufferIdx][whichCase], 0, roundUp * word_per_ans_2 * sizeof(uint));
        cudaMalloc ( ( void ** ) &_answers, roundUp * word_per_ans_2 * sizeof ( uint ) );

        // to reset the values inside the array _answers
        // cudaMemcpy(_answers, badAnswers[doubleBufferIdx][whichCase],
        //          roundUp * word_per_ans_2 * sizeof(uint),
        //          cudaMemcpyHostToDevice);
        if ( numMismatch <= 3 ) // 3 mismatch
            kernel <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( whichCase, _queries, _readLengths, numBads, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0,
              bwt->textLength, _answers, _isBad, 1, numMismatch, sa_range_allowed_2, word_per_ans_2, isExactNumMismatch );
        else if ( whichCase < 5 ) // 4 mismatch and case no. 0 - 4
            kernel_4mismatch_1 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( whichCase, _queries, _readLengths, numBads, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0,
              bwt->textLength, _answers, _isBad, 1, sa_range_allowed_2, word_per_ans_2, isExactNumMismatch );
        else   // 4 mismatch and case no. 5 - 9
            kernel_4mismatch_2 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( whichCase, _queries, _readLengths, numBads, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0,
              bwt->textLength, _answers, _isBad, 1, sa_range_allowed_2, word_per_ans_2, isExactNumMismatch );

        cudaMemcpy ( badAnswers[doubleBufferIdx][whichCase], _answers,
                     roundUp * word_per_ans_2 * sizeof ( uint ), cudaMemcpyDeviceToHost );
        // free the memory in the device
        // cudaFree(_isBad);
        cudaFree ( _queries );
        cudaFree ( _answers );
        cudaFree ( _readLengths );
    }
}

// perform round1 alignment in GPU for 1 mismatch (no pipeline)
void perform_round1_alignment_no_pipeline ( uint * nextQuery, uint * nextReadLength, uint * answers[MAX_NUM_CASES],
        uint numMismatch, uint numCases, uint sa_range_allowed, uint wordPerQuery, uint word_per_ans,
        bool isExactNumMismatch, uint blocksNeeded, ullint batchSize,
        Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc )
{
    cudaError_t gpuErr;
    uint * _queries, *_readLengths, *_answers;
    bool * _isBad;
    ullint roundUp = ( batchSize + 31 ) / 32 * 32;
    // allocated device memory for bad read indicator
    bool * isBad;
    isBad = ( bool * ) malloc ( roundUp * sizeof ( bool ) ); // an array to store bad read indicator
    memset ( isBad, 0, roundUp );
    gpuErr = cudaMalloc ( ( void ** ) &_isBad, roundUp * sizeof ( bool ) );
    BWT * bwt = index->sraIndex->bwt;
    BWT * revBwt = index->sraIndex->rev_bwt;

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MALLOC FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    // to initialize the array _isBad
    gpuErr = cudaMemcpy ( _isBad, isBad,
                          roundUp * sizeof ( bool ),
                          cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    // printf("[perform_round1_alignment] sa_range_allowed = %u; word_per_ans = %u\n", sa_range_allowed, word_per_ans);
    // allocate device memory for queries and answers
    cudaMalloc ( ( void ** ) &_queries, roundUp * wordPerQuery * sizeof ( uint ) );
    cudaMalloc ( ( void ** ) &_readLengths, roundUp * sizeof ( uint ) );
    cudaMalloc ( ( void ** ) &_answers, roundUp * word_per_ans * sizeof ( uint ) );
    gpuErr = cudaMemcpy ( _queries, nextQuery, roundUp * wordPerQuery * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    gpuErr = cudaMemcpy ( _readLengths, nextReadLength, roundUp * sizeof ( uint ), cudaMemcpyHostToDevice );

    if ( gpuErr != cudaSuccess )
    {
        printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
        exit ( 1 );
    }

    for ( uint caseno = 0; caseno < numCases; caseno++ )
    {
        // =======================================
        // | GPU-1: FOR EACH CASE                |
        // =======================================
        if ( numMismatch <= 3 )
            kernel <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( caseno, _queries, _readLengths, batchSize, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0, bwt->textLength,
              _answers, _isBad, 0 , numMismatch, sa_range_allowed, word_per_ans, isExactNumMismatch );
        else if ( caseno < 5 ) // 4 mismatch and case 0 - 4
            kernel_4mismatch_1 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( caseno, _queries, _readLengths, batchSize, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0, bwt->textLength,
              _answers, _isBad, 0 , sa_range_allowed, word_per_ans, isExactNumMismatch );
        else  // 4 mismatch and case 5 - 9
            kernel_4mismatch_2 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( caseno, _queries, _readLengths, batchSize, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0, bwt->textLength,
              _answers, _isBad, 0 , sa_range_allowed, word_per_ans, isExactNumMismatch );

        gpuErr = cudaMemcpy ( answers[caseno], _answers, roundUp * word_per_ans * sizeof ( uint ), cudaMemcpyDeviceToHost );

        if ( gpuErr != cudaSuccess )
        {
            printf ( "CUDA MEMCOPY FAILED .. %s(%d)\n", cudaGetErrorString ( gpuErr ), gpuErr );
            exit ( 1 );
        }
    }

    // free the memories
    free ( isBad );
    cudaFree ( _isBad );
    cudaFree ( _queries );
    cudaFree ( _readLengths );
    cudaFree ( _answers );
}

// perform round2 alignment in GPU for 1 mismatch
void perform_round2_alignment_no_pipeline ( uint * queries, uint * readLengths, uint * answers[MAX_NUM_CASES],
        uint numMismatch, uint numCases, uint sa_range_allowed_2, uint wordPerQuery, uint word_per_ans, uint word_per_ans_2,
        bool isExactNumMismatch, uint blocksNeeded, ullint batchSize,
        Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
        uint processedQuery, uint * badReadIndices[MAX_NUM_CASES],
        uint * badAnswers[MAX_NUM_CASES] )
{
    uint * badQueries;
    uint * badReadLengths;
    BWT * bwt = index->sraIndex->bwt;
    BWT * revBwt = index->sraIndex->rev_bwt;

    for ( int whichCase = 0; whichCase < numCases; ++whichCase )
    {
        ullint numBads = 0;

        // Count number of bad reads and
        for ( ullint readId = 0; readId < batchSize; readId++ )
        {
            ullint srcOffset = ( ( readId ) / 32 * 32 * word_per_ans + readId % 32 );
            numBads += ( answers[whichCase][srcOffset] > 0xFFFFFFFD );
        }

        if ( numBads == 0 )
        { continue; }

        // Allocate memory and copy bad reads to another array
        ullint roundUp = ( numBads + 31 ) / 32 * 32;
        badQueries = ( uint * ) malloc ( roundUp * wordPerQuery * sizeof ( uint ) );
        // printf("size of badQueries = %u\n", roundUp * wordPerQuery);
        badReadLengths = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
        // printf("size of badReadLengths = %u\n", roundUp);
        numBads = 0;

        for ( ullint readId = 0; readId < batchSize; ++readId )
        {
            ullint srcOffset = ( ( readId ) / 32 * 32 * word_per_ans + readId % 32 );
            ullint srcQueryOffset = ( processedQuery + readId ) / 32 * 32 * wordPerQuery + readId % 32;

            if ( answers[whichCase][srcOffset] > 0xFFFFFFFD ) // is bad
            {
                badReadIndices[whichCase][numBads] = readId;
                ullint targetOffset = numBads / 32 * 32 * wordPerQuery + numBads % 32;

                for ( ullint i = 0; i < wordPerQuery; ++i ) // copy each word of the read
                {
                    badQueries[targetOffset + i * 32] = * ( queries + ( srcQueryOffset + i * 32 ) );
                }

                badReadLengths[numBads] = readLengths[processedQuery + readId];
                numBads++;
            }
        }

        // Copy reads to device memory, call kernel, and copy results back to host
        uint * _queries, *_readLengths, *_answers;
        bool * _isBad = NULL;
        // bool *isBad = (bool*) malloc((batchSize + 31) / 32 * 32 * sizeof(bool)); // an array to store bad read indicator
        // memset(isBad,0, (batchSize + 31) / 32 * 32 * sizeof(bool));
        // cudaMalloc((void**)&_isBad, (batchSize + 31) / 32 * 32 * sizeof(bool));
        // cudaMemcpy(_isBad, isBad,(batchSize + 31) / 32 * 32 * sizeof(bool), cudaMemcpyHostToDevice);
        // free(isBad);
        cudaMalloc ( ( void ** ) &_queries, roundUp * wordPerQuery * sizeof ( uint ) );
        cudaMemcpy ( _queries, badQueries,
                     roundUp * wordPerQuery * sizeof ( uint ), cudaMemcpyHostToDevice );
        free ( badQueries );
        cudaMalloc ( ( void ** ) &_readLengths, roundUp * sizeof ( uint ) );
        cudaMemcpy ( _readLengths, badReadLengths,
                     roundUp * sizeof ( uint ), cudaMemcpyHostToDevice );
        free ( badReadLengths );
        // reset the values of the array badAnswers
        // memset(badAnswers[whichCase], 0, roundUp * word_per_ans_2 * sizeof(uint));
        cudaMalloc ( ( void ** ) &_answers, roundUp * word_per_ans_2 * sizeof ( uint ) );

        // to reset the values inside the array _answers
        // cudaMemcpy(_answers, badAnswers[whichCase],
        //          roundUp * word_per_ans_2 * sizeof(uint),
        //          cudaMemcpyHostToDevice);
        if ( numMismatch <= 3 ) // 3 mismatch
            kernel <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( whichCase, _queries, _readLengths, numBads, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0,
              bwt->textLength, _answers, _isBad, 1, numMismatch, sa_range_allowed_2, word_per_ans_2, isExactNumMismatch );
        else if ( whichCase < 5 ) // 4 mismatch and case no. 0 - 4
            kernel_4mismatch_1 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( whichCase, _queries, _readLengths, numBads, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0,
              bwt->textLength, _answers, _isBad, 1, sa_range_allowed_2, word_per_ans_2, isExactNumMismatch );
        else   // 4 mismatch and case no. 5 - 9
            kernel_4mismatch_2 <<< blocksNeeded, THREADS_PER_BLOCK>>>
            ( whichCase, _queries, _readLengths, numBads, wordPerQuery,
              _bwt, _occ, bwt->inverseSa0,
              _revBwt, _revOcc, revBwt->inverseSa0,
              bwt->textLength, _answers, _isBad, 1, sa_range_allowed_2, word_per_ans_2, isExactNumMismatch );

        cudaMemcpy ( badAnswers[whichCase], _answers,
                     roundUp * word_per_ans_2 * sizeof ( uint ), cudaMemcpyDeviceToHost );
        // free the memory in the device
        // cudaFree(_isBad);
        cudaFree ( _queries );
        cudaFree ( _answers );
        cudaFree ( _readLengths );
    }
}


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
                           BothUnalignedPairs ** bothUnalignedPairs )
{
#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    double startTime2 = setStartTime ();
    double lastEventTime2 = 0;
    double currEventTime2;
#endif
    //Modified to be double buffered
    uint * answers[2][MAX_NUM_CASES];
    uint * badReadIndices[2][MAX_NUM_CASES]; // IDs of bad reads for each case
    uint * badAnswers[2][MAX_NUM_CASES]; // stores answers (and reads temporarily) for 2nd round
    uint threadBadStarts[2][MAX_NUM_CPU_THREADS][MAX_NUM_CASES];
    uint threadBadCounts[2][MAX_NUM_CPU_THREADS][MAX_NUM_CASES];
    uint * unAlignedId[MAX_NUM_CPU_THREADS];
    int doubleBufferIdx = 0; // expect to be either 0 or 1
    // Host alignment model
    SRAModel * SRAMismatchModel[MAX_READ_LENGTH + 1];
    SRAModel * SRAMismatchModel2[MAX_READ_LENGTH + 1];
    SRAModel * SRAMismatchModel_neg[MAX_READ_LENGTH + 1];
    SRAModel * SRAMismatchModel2_neg[MAX_READ_LENGTH + 1];
    HostKernelArguements hostKernelArguments[MAX_NUM_CPU_THREADS];
    // Host multi-threading variables
    int threadId;
    pthread_t threads[MAX_NUM_CPU_THREADS];
    uint threadBucketSize[MAX_NUM_CPU_THREADS];
    numOfAnswer = 0;
    numOfAlignedRead = 0;
    numOfUnPaired = 0;
    uint blocksNeeded;
    ullint batchSize;
    uint numCases;
    uint sa_range_allowed_1;
    uint sa_range_allowed_2;
    char skip_round_2;
    ullint queriesLeft = numQueries;
    uint maxReadLength = input_options.maxReadLength;
    uint word_per_ans;
    uint word_per_ans_2;
    int i;

    SRASetting _mapqTempSRASetting;

    SRAIndex * sraIndex = index->sraIndex;

    // initialization of arrays
    for ( i = 0; i < 2; i++ )
    {
        for ( int j = 0; j < MAX_NUM_CASES; j++ )
        {
            badReadIndices[i][j] = NULL;
            badAnswers[i][j] = NULL;
            answers[i][j] = NULL;

            for ( int k = 0; k < MAX_NUM_CPU_THREADS; k++ )
            {
                threadBadStarts[i][k][j] = 0;
                threadBadCounts[i][k][j] = 0;
            }
        }
    }

    // initialization of unAlignedId
    for ( int i = 0; i < ini_params.Ini_NumOfCpuThreads; i++ )
    {
        unAlignedId[i] = ( uint * ) malloc ( maxBatchSize * sizeof ( uint ) );
    }

    // obtain the number of cases for this number of mismatch
    // and obtain the number of SA ranges allowed
    getParametersForThisMismatch ( numMismatch, numCases, sa_range_allowed_1,
                                   sa_range_allowed_2, skip_round_2, word_per_ans, word_per_ans_2 );

    // For single-read alignment (and not a long-read mode),
    // unique best and random best requires only 1 answer
    if ( input_options.readType == SINGLE_READ && seedLengths == NULL )
    {
        if ( input_options.alignmentType == OUTPUT_UNIQUE_BEST )
        {
            sa_range_allowed_1 = 1;
            word_per_ans = 2;
            skip_round_2 = 1;
        }

        if ( input_options.alignmentType == OUTPUT_RANDOM_BEST )
        {
            sa_range_allowed_1 = 2;
            word_per_ans = 4;
            skip_round_2 = 1;
        }
    }

    // For 4-mismatch, if all-valid alignment
    // then sa_range_allowed_1 is changed from 1 to 2
    if ( numMismatch == 4 && input_options.alignmentType == OUTPUT_ALL_VALID )
    {
        sa_range_allowed_1 = 2;
        word_per_ans = 4;
    }

    // For paired-end read alignment, (all-best or all-valid) and SAM format
    bool needOutputMAPQ = ( input_options.alignmentType == OUTPUT_ALL_VALID ||
                            input_options.alignmentType == OUTPUT_ALL_BEST ) &&
                          ( input_options.outputFormat == SRA_OUTPUT_FORMAT_SAM_API );

    if ( needOutputMAPQ && ( numMismatch == 2 ) )
    {
        sa_range_allowed_2 = 256;
        word_per_ans_2 = 512;
    }

    // set the multi-threading arguments
    setHostKernelArguments ( hostKernelArguments, threads, ini_params, index, maxReadLength, wordPerQuery, word_per_ans, &input_options );

    for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        hostKernelArguments[threadId].upkdQualities = upkdQualities;
        hostKernelArguments[threadId].upkdQueryNames = upkdQueryNames;
        hostKernelArguments[threadId].SRAMismatchModel = SRAMismatchModel;
        hostKernelArguments[threadId].SRAMismatchModel2 = SRAMismatchModel2;
        hostKernelArguments[threadId].SRAMismatchModel_neg = SRAMismatchModel_neg;
        hostKernelArguments[threadId].SRAMismatchModel2_neg = SRAMismatchModel2_neg;

        if ( currOutputFileName != NULL )
        { hostKernelArguments[threadId].outputFileName = currOutputFileName[threadId]; }
        else
        { hostKernelArguments[threadId].outputFileName = NULL; }

        hostKernelArguments[threadId].readLengths = readLengths;
        hostKernelArguments[threadId].seedLengths = seedLengths;
        hostKernelArguments[threadId].readIDs = readIDs;
        hostKernelArguments[threadId].accumReadNum = accumReadNum;
        hostKernelArguments[threadId].unAlignedIDs = unAlignedId[threadId];
        hostKernelArguments[threadId].unAlignedOcc = 0;
    }

    ullint roundUp = ( maxBatchSize + 31 ) / 32 * 32;

    // allocate host memory for answers
    for ( i = 0; i < 2; i++ )
    {
        for ( int j = 0; j < numCases; j++ )
        {
            answers[i][j] = ( uint * ) malloc ( roundUp * word_per_ans * sizeof ( uint ) );
            //memset(answers[i][j],0xFFFFFFFF, roundUp * word_per_ans * sizeof ( uint ) );
        }
    }

    // if there is only 3G GPU memory, then skip round 2
    if ( ini_params.Ini_GPUMemory == 3 )
    { skip_round_2 = 1; }

    // printParameters(input_options, ini_params);
    // printf("numMismatch = %u; numCases = %u; sa_range_allowed_1 = %u; sa_range_allowed_2 = %u; skip_round_2 = %i; word_per_ans = %u; word_per_ans_2 = %u \n", numMismatch, numCases, sa_range_allowed_1, sa_range_allowed_2, skip_round_2, word_per_ans, word_per_ans_2);

    // update the num of mismatch of the model for host alignment
    for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        hostKernelArguments[threadId].sraQuerySettings.MaxError = numMismatch;
        hostKernelArguments[threadId].sraQuerySettings.MaxNBMismatch = 0;
        hostKernelArguments[threadId].numCases = numCases;
        hostKernelArguments[threadId].numCases2 = GetNumCases ( numMismatch + 1 );
        hostKernelArguments[threadId].sa_range_allowed_1 = sa_range_allowed_1;
        hostKernelArguments[threadId].sa_range_allowed_2 = sa_range_allowed_2;
        hostKernelArguments[threadId].maxNumMismatch = numMismatch;

        if ( readInputForDP != NULL )
        { hostKernelArguments[threadId].readInput = readInputForDP[threadId]; }

        if ( readInputForNewDP != NULL )
        { hostKernelArguments[threadId].readInputForNewDefault = readInputForNewDP[threadId]; }

        if ( otherSoap3Result != NULL )
        { hostKernelArguments[threadId].otherSoap3Result = otherSoap3Result[threadId]; }

        if ( bothUnalignedPairs != NULL )
        { hostKernelArguments[threadId].bothUnalignedPairs = bothUnalignedPairs[threadId]; }

        //Update the SAM output file ptr in QuerySetting
        if ( currSamOutputFilePtr != NULL )
        { hostKernelArguments[threadId].sraQuerySettings.SAMOutFilePtr = currSamOutputFilePtr[threadId]; }
        else
        { hostKernelArguments[threadId].sraQuerySettings.SAMOutFilePtr = NULL; }
    }

    for ( i = 0; i <= MAX_READ_LENGTH; i++ )
    {
        SRAMismatchModel[i] = NULL;
        SRAMismatchModel2[i] = NULL;
        SRAMismatchModel_neg[i] = NULL;
        SRAMismatchModel2_neg[i] = NULL;
    }

    if ( seedLengths == NULL )
    {
        for ( unsigned readId = 0; readId < numQueries; ++readId )
        {
            if ( SRAMismatchModel[readLengths[readId]] == NULL )
            {
                SRAMismatchModel[readLengths[readId]] = SRAModelConstruct ( readLengths[readId], QUERY_POS_STRAND, & ( hostKernelArguments[0].sraQuerySettings ), sraIndex,  ini_params.Ini_HostAlignmentModel );

            }

            if ( SRAMismatchModel_neg[readLengths[readId]] == NULL )
            {
                SRAMismatchModel_neg[readLengths[readId]] = SRAModelConstruct ( readLengths[readId], QUERY_NEG_STRAND, & ( hostKernelArguments[0].sraQuerySettings ), sraIndex,  ini_params.Ini_HostAlignmentModel );

            }
        }
    }
    else
    {
        for ( unsigned readId = 0; readId < numQueries; ++readId )
        {
            if ( SRAMismatchModel[seedLengths[readId]] == NULL )
            {
                SRAMismatchModel[seedLengths[readId]] = SRAModelConstruct ( seedLengths[readId], QUERY_POS_STRAND, & ( hostKernelArguments[0].sraQuerySettings ), sraIndex, ini_params.Ini_HostAlignmentModel );
            }

            if ( SRAMismatchModel_neg[seedLengths[readId]] == NULL )
            {
                SRAMismatchModel_neg[seedLengths[readId]] = SRAModelConstruct ( seedLengths[readId], QUERY_NEG_STRAND, & ( hostKernelArguments[0].sraQuerySettings ), sraIndex, ini_params.Ini_HostAlignmentModel );
            }
        }
    }

    // if (ALL-VALID or ALL-BEST) and SAM format, then need to update the model for mismatch+1
    if ( needOutputMAPQ && numMismatch < 4 )
    {
        memcpy ( &_mapqTempSRASetting, & ( hostKernelArguments[0].sraQuerySettings ), sizeof ( SRASetting ) );
        _mapqTempSRASetting.MaxError = numMismatch + 1;

        if ( seedLengths == NULL )
        {
            for ( unsigned readId = 0; readId < numQueries; ++readId )
            {
                if ( SRAMismatchModel2[readLengths[readId]] == NULL )
                {
                    SRAMismatchModel2[readLengths[readId]] = SRAModelConstruct ( readLengths[readId], QUERY_POS_STRAND, &_mapqTempSRASetting, sraIndex, ini_params.Ini_HostAlignmentModel );
                }

                if ( SRAMismatchModel2_neg[readLengths[readId]] == NULL )
                {
                    SRAMismatchModel2_neg[readLengths[readId]] = SRAModelConstruct ( readLengths[readId], QUERY_NEG_STRAND, &_mapqTempSRASetting, sraIndex, ini_params.Ini_HostAlignmentModel );
                }
            }
        }
        else
        {
            for ( unsigned readId = 0; readId < numQueries; ++readId )
            {
                if ( SRAMismatchModel2[seedLengths[readId]] == NULL )
                {
                    SRAMismatchModel2[seedLengths[readId]] = SRAModelConstruct ( seedLengths[readId], QUERY_POS_STRAND, &_mapqTempSRASetting, sraIndex, ini_params.Ini_HostAlignmentModel );
                }

                if ( SRAMismatchModel2_neg[seedLengths[readId]] == NULL )
                {
                    SRAMismatchModel2_neg[seedLengths[readId]] = SRAModelConstruct ( seedLengths[readId], QUERY_NEG_STRAND, &_mapqTempSRASetting, sraIndex, ini_params.Ini_HostAlignmentModel );
                }
            }
        }
    }

    uint * nextQuery = queries;
    uint * nextReadLength;

    if ( seedLengths != NULL )
    { nextReadLength = seedLengths; }
    else
    { nextReadLength = readLengths; }

    uint * nextUnAlignedReads = unAlignedReads;
#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    currEventTime2 = getElapsedTime ( startTime2 );
    printf ( "[Main] Time elapsed for initialization: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
    lastEventTime2 = currEventTime2;
#endif

    while ( queriesLeft > 0 )
    {
        uint processedQuery = numQueries - queriesLeft;

        if ( queriesLeft > maxBatchSize )
        {
            blocksNeeded = NUM_BLOCKS;
            batchSize = maxBatchSize;
        }
        else
        {
            blocksNeeded = ( queriesLeft + THREADS_PER_BLOCK * QUERIES_PER_THREAD - 1 ) /
                           ( THREADS_PER_BLOCK * QUERIES_PER_THREAD );
            batchSize = queriesLeft;
        }

        // allocate the reads to different CPU threads
        uint roughBucketSize = batchSize / ini_params.Ini_NumOfCpuThreads;
        // roughBucketSize has to be divible by 2
        roughBucketSize = roughBucketSize / 2 * 2;
        uint batchSizeUnalloc = batchSize;

        for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads - 1; threadId++ )
        {
            threadBucketSize[threadId] = roughBucketSize;
            batchSizeUnalloc -= roughBucketSize;
        }

        threadBucketSize[ini_params.Ini_NumOfCpuThreads - 1] = batchSizeUnalloc;
        // perform first round alignment in GPU
        perform_round1_alignment ( nextQuery, nextReadLength, answers,
                                   numMismatch, numCases, sa_range_allowed_1, wordPerQuery,
                                   word_per_ans, false, doubleBufferIdx, blocksNeeded, batchSize,
                                   index, _bwt, _revBwt, _occ, _revOcc );
#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
        printf ( "[Main] # of mismatches allowed: %u with # of sa ranges: %u\n", numMismatch, sa_range_allowed_1 );
        currEventTime2 = getElapsedTime ( startTime2 );
        printf ( "[Main] Time elapsed for first round alignment in GPU: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
        lastEventTime2 = currEventTime2;
#endif

        if ( skip_round_2 == 0 )
        {
            // =======================================
            // | GPU-2                               |
            // =======================================
            // perform second round alignment in GPU
            for ( int whichCase = 0; whichCase < numCases; ++whichCase )
            {
                ullint numBads = 0;

                // Count number of bad reads and
                for ( ullint readId = 0; readId < batchSize; readId++ )
                {
                    ullint srcOffset = ( ( readId ) / 32 * 32 * word_per_ans + readId % 32 );
                    numBads += ( answers[doubleBufferIdx][whichCase][srcOffset] > 0xFFFFFFFD );
                }

                if ( numBads > 0 )
                {
                    // Allocate memory and copy bad reads to another array
                    ullint roundUp = ( numBads + 31 ) / 32 * 32;
                    badReadIndices[doubleBufferIdx][whichCase] = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
                    // printf("size of badReadIndices[%i][%i] = %u\n",doubleBufferIdx, whichCase, roundUp);
                    badAnswers[doubleBufferIdx][whichCase] = ( uint * ) malloc ( roundUp * word_per_ans_2 * sizeof ( uint ) );
                    // printf("size of badAnswers[%i][%i] = %u\n",doubleBufferIdx, whichCase, roundUp * word_per_ans_2);
                }

                // printf("numBads == %u in case %i of mismatch %u\n", numBads, whichCase, numMismatch);
            }

            if ( seedLengths == NULL )
            {
                perform_round2_alignment ( queries, readLengths, answers,
                                           numMismatch, numCases, sa_range_allowed_2, wordPerQuery, word_per_ans,
                                           word_per_ans_2, false, doubleBufferIdx, blocksNeeded, batchSize,
                                           index, _bwt, _revBwt, _occ, _revOcc,
                                           processedQuery, badReadIndices, badAnswers );
            }
            else
            {
                perform_round2_alignment ( queries, seedLengths, answers,
                                           numMismatch, numCases, sa_range_allowed_2, wordPerQuery, word_per_ans,
                                           word_per_ans_2, false, doubleBufferIdx, blocksNeeded, batchSize,
                                           index, _bwt, _revBwt, _occ, _revOcc,
                                           processedQuery, badReadIndices, badAnswers );
            }

            for ( int whichCase = 0; whichCase < numCases; ++whichCase )
            {
                ullint readId = 0;

                for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; threadId ++ )
                {
                    uint threadNumBad = 0;

                    // printf("whichCase = %i; threadId = %i; threadBucketSize[threadId] = %u\n",
                    //        whichCase, threadId, threadBucketSize[threadId]);
                    for ( uint n = 0 ; n < threadBucketSize[threadId]; ++n )
                    {
                        ullint srcOffset = ( readId / 32 * 32 * word_per_ans + readId % 32 );
                        threadNumBad += ( answers[doubleBufferIdx][whichCase][srcOffset] > 0xFFFFFFFD );
                        readId++;
                    }

                    threadBadCounts[doubleBufferIdx][threadId][whichCase] = threadNumBad;
                    // printf("threadBadCounts[%i][%i][%i] = %u\n",
                    //    doubleBufferIdx, threadId, whichCase, threadBadCounts[doubleBufferIdx][threadId][whichCase]);
                }
            }

            for ( int whichCase = 0; whichCase < numCases; ++whichCase )
            {
                threadBadStarts[doubleBufferIdx][0][whichCase] = 0;

                for ( threadId = 1; threadId < ini_params.Ini_NumOfCpuThreads; threadId ++ )
                {
                    threadBadStarts[doubleBufferIdx][threadId][whichCase] = threadBadStarts[doubleBufferIdx][threadId - 1][whichCase] + threadBadCounts[doubleBufferIdx][threadId - 1][whichCase];
                    // printf("threadBadStarts[%i][%i][%i] = %u\n",
                    //  doubleBufferIdx, threadId, whichCase, threadBadStarts[doubleBufferIdx][threadId][whichCase]);
                }
            }

#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
            currEventTime2 = getElapsedTime ( startTime2 );
            printf ( "[Main] Time elapsed for second round alignment in GPU: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
            printf ( "[Main] # of sa ranges allowed: %u \n", sa_range_allowed_2 );
            lastEventTime2 = currEventTime2;
#endif
        } // skip the round 2

        for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads ; threadId++ )
        {
            if ( threads[threadId] != 0 )
            {
                if ( pthread_join ( threads[threadId], NULL ) )
                { fprintf ( stderr, "[Main:Thread%u] Crash!\n", threadId ), exit ( 1 ); }

                // pthread_detach(threads[threadId]);
                numOfAnswer += hostKernelArguments[threadId].alignedOcc;
                numOfAlignedRead += hostKernelArguments[threadId].alignedReads;
                threads[threadId] = 0;

                // consolidate the unAlignedId;
                if ( hostKernelArguments[threadId].unAlignedOcc > 0 )
                {
                    memcpy ( nextUnAlignedReads, unAlignedId[threadId],
                             ( ullint ) hostKernelArguments[threadId].unAlignedOcc * sizeof ( uint ) );
                    nextUnAlignedReads += hostKernelArguments[threadId].unAlignedOcc;
                    numOfUnPaired += hostKernelArguments[threadId].unAlignedOcc;
                }
            }
        }

#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
        currEventTime2 = getElapsedTime ( startTime2 );
        printf ( "[Main] Time elapsed for waiting for CPU threads: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
        lastEventTime2 = currEventTime2;
#endif

        for ( int whichCase = 0; whichCase < numCases; ++whichCase )
        {
            if ( badReadIndices[1 - doubleBufferIdx][whichCase] != NULL )
            {
                free ( badReadIndices[1 - doubleBufferIdx][whichCase] );
                badReadIndices[1 - doubleBufferIdx][whichCase] = NULL;
            }

            if ( badAnswers[1 - doubleBufferIdx][whichCase] != NULL )
            {
                free ( badAnswers[1 - doubleBufferIdx][whichCase] );
                badAnswers[1 - doubleBufferIdx][whichCase] = NULL;
            }
        }

        // =======================================
        // | CPU: Thread #0,1,2,3,...            |
        // =======================================
        unsigned int threadProcessedQuery = 0;

        for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads ; threadId++ )
        {
            hostKernelArguments[threadId].batchFirstReadId = processedQuery;
            hostKernelArguments[threadId].skipFirst = threadProcessedQuery;
            hostKernelArguments[threadId].numQueries = threadBucketSize[threadId];
            hostKernelArguments[threadId].word_per_query = wordPerQuery;
            hostKernelArguments[threadId].queries = queries;
            hostKernelArguments[threadId].answers = answers[doubleBufferIdx];
            hostKernelArguments[threadId].alignedOcc = 0;
            hostKernelArguments[threadId].alignedReads = 0;
            hostKernelArguments[threadId].badReadIndices = badReadIndices[doubleBufferIdx];
            hostKernelArguments[threadId].badAnswers = badAnswers[doubleBufferIdx];
            hostKernelArguments[threadId].badStartOffset = threadBadStarts[doubleBufferIdx][threadId];
            hostKernelArguments[threadId].badCountOffset = threadBadCounts[doubleBufferIdx][threadId];
            hostKernelArguments[threadId].outputGoodReads = TRUE;
            hostKernelArguments[threadId].skip_round_2 = skip_round_2;
            hostKernelArguments[threadId].isTerminalCase = isTerminalCase;

            if ( pthread_create ( & ( threads[threadId] ), NULL, hostKernelThreadWrapper, ( void * ) & ( hostKernelArguments[threadId] ) ) )
            { fprintf ( stderr, "[Main:Threads%u] Can't create hostKernelThreadWrapper\n", threadId ), exit ( 1 ); }

            threadProcessedQuery += threadBucketSize[threadId];
        }

        // Swap the double buffer
        doubleBufferIdx = 1 - doubleBufferIdx;
        queriesLeft -= batchSize;
        nextQuery += batchSize * wordPerQuery;
        nextReadLength += batchSize;
    }

    for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads ; threadId++ )
    {
        if ( threads[threadId] != 0 )
        {
            if ( pthread_join ( threads[threadId], NULL ) )
            { fprintf ( stderr, "[Main:Thread%u] Crash!\n", threadId ), exit ( 1 ); }

            // pthread_detach(threads[threadId]);
            numOfAnswer += hostKernelArguments[threadId].alignedOcc;
            numOfAlignedRead += hostKernelArguments[threadId].alignedReads;
            threads[threadId] = 0;

            // consolidate the unAlignedId;
            if ( hostKernelArguments[threadId].unAlignedOcc > 0 )
            {
                memcpy ( nextUnAlignedReads, unAlignedId[threadId],
                         ( ullint ) hostKernelArguments[threadId].unAlignedOcc * sizeof ( uint ) );
                nextUnAlignedReads += hostKernelArguments[threadId].unAlignedOcc;
                numOfUnPaired += hostKernelArguments[threadId].unAlignedOcc;
            }
        }
    }

    for ( int whichCase = 0; whichCase < numCases; ++whichCase )
    {
        if ( badReadIndices[1 - doubleBufferIdx][whichCase] != NULL )
        {
            free ( badReadIndices[1 - doubleBufferIdx][whichCase] );
            badReadIndices[1 - doubleBufferIdx][whichCase] = NULL;
        }

        if ( badAnswers[1 - doubleBufferIdx][whichCase] != NULL )
        {
            free ( badAnswers[1 - doubleBufferIdx][whichCase] );
            badAnswers[1 - doubleBufferIdx][whichCase] = NULL;
        }
    }

    // CLEAN UP for each MISMATCH iteration                                                                          |
    for ( i = 0; i <= MAX_READ_LENGTH; i++ )
    {
        if ( SRAMismatchModel[i] != NULL )
        {
            SRAModelFree ( SRAMismatchModel[i] );
            SRAMismatchModel[i] = NULL;
        }

        if ( SRAMismatchModel2[i] != NULL )
        {
            SRAModelFree ( SRAMismatchModel2[i] );
            SRAMismatchModel2[i] = NULL;
        }

        if ( SRAMismatchModel_neg[i] != NULL )
        {
            SRAModelFree ( SRAMismatchModel_neg[i] );
            SRAMismatchModel_neg[i] = NULL;
        }

        if ( SRAMismatchModel2_neg[i] != NULL )
        {
            SRAModelFree ( SRAMismatchModel2_neg[i] );
            SRAMismatchModel2_neg[i] = NULL;
        }
    }

    for ( threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        OCCFree ( hostKernelArguments[threadId].occ );
        free ( unAlignedId[threadId] );
    }

    for ( i = 0; i < 2; i++ )
    {
        for ( int j = 0; j < numCases; j++ )
        {
            free ( answers[i][j] );
        }
    }
}

// pair-end alignment: for random-best
// 4-phases [0,1,2,4]
void four_phases_alignment ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                             ullint maxBatchSize, uint numQueries, uint accumReadNum,
                             Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                             IniParams ini_params, InputOptions input_options,
                             char * upkdQualities,
                             uint * unAlignedReads, uint & numOfUnPaired,
                             uint * readIDs, char * upkdQueryNames,
                             char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                             unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                             ReadInputForDP ** readInputForDP,
                             ReadInputForDP ** readInputForNewDP,
                             ReadInputForDP ** otherSoap3Result,
                             BothUnalignedPairs ** bothUnalignedPairs )
{
    // for performing pair-end alignment with random-best output
    unsigned long long numOfAnswer0 = 0;
    uint numOfAlignedRead0 = 0;
    unsigned long long numOfAnswer1 = 0;
    uint numOfAlignedRead1 = 0;
    unsigned long long numOfAnswer2 = 0;
    uint numOfAlignedRead2 = 0;
    unsigned long long numOfAnswer3 = 0;
    uint numOfAlignedRead3 = 0;
    //********************************//
    // Phase 0: Perform 0 Alignment //
    //********************************//
    uint numMismatch_phase0 = 0;
    all_valid_alignment ( queries, readLengths, NULL, numMismatch_phase0, wordPerQuery,
                          maxBatchSize, numQueries, accumReadNum,
                          index, _bwt, _revBwt, _occ, _revOcc,
                          ini_params, input_options,
                          upkdQualities,
                          unAlignedReads, numOfUnPaired,
                          readIDs, upkdQueryNames,
                          currOutputFileName,  currSamOutputFilePtr,
                          numOfAnswer0, numOfAlignedRead0, numMismatch == 0,
                          readInputForDP, readInputForNewDP, otherSoap3Result, bothUnalignedPairs );

    if ( numMismatch > 0 && numOfAlignedRead0 < numQueries )
    {
        //********************************//
        // Phase 1: Perform 0/1 Alignment //
        //********************************//
        uint numMismatch_phase1 = 1;
        // pack the reads which are not paired
        packUnPairedReads ( queries, readIDs, readLengths, unAlignedReads,
                            wordPerQuery, numOfUnPaired, maxBatchSize );
        numQueries = numOfUnPaired;
        all_valid_alignment ( queries, readLengths, NULL, numMismatch_phase1, wordPerQuery,
                              maxBatchSize, numQueries, accumReadNum,
                              index, _bwt, _revBwt, _occ, _revOcc,
                              ini_params, input_options,
                              upkdQualities,
                              unAlignedReads, numOfUnPaired,
                              readIDs, upkdQueryNames,
                              currOutputFileName,  currSamOutputFilePtr,
                              numOfAnswer1, numOfAlignedRead1, numMismatch == 1,
                              readInputForDP, readInputForNewDP, otherSoap3Result, bothUnalignedPairs );

        // printf("numOfAlignedRead1 = %u\n", numOfAlignedRead1);
        // printf("numOfUnPaired = %u\n", numOfUnPaired);

        if ( numMismatch > 1 && numOfAlignedRead1 < numQueries )
        {
            //****************************************//
            // Phase 2: Perform 0/1/2 alignment     //
            //****************************************//
            uint numMismatch_phase2 = 2;
            // pack the reads which are not paired
            packUnPairedReads ( queries, readIDs, readLengths, unAlignedReads,
                                wordPerQuery, numOfUnPaired, maxBatchSize );
            numQueries = numOfUnPaired;
            all_valid_alignment ( queries, readLengths, NULL, numMismatch_phase2, wordPerQuery,
                                  maxBatchSize, numQueries, accumReadNum,
                                  index, _bwt, _revBwt, _occ, _revOcc,
                                  ini_params, input_options,
                                  upkdQualities,
                                  unAlignedReads, numOfUnPaired,
                                  readIDs, upkdQueryNames,
                                  currOutputFileName, currSamOutputFilePtr,
                                  numOfAnswer2, numOfAlignedRead2, numMismatch == 2,
                                  readInputForDP, readInputForNewDP, otherSoap3Result, bothUnalignedPairs );

            // printf("numOfAlignedRead2 = %u\n", numOfAlignedRead2);
            // printf("numOfUnPaired = %u\n", numOfUnPaired, 0);

            if ( numMismatch > 2 && numOfAlignedRead2 < numQueries )
            {
                //******************************************//
                // Phase 3: Perform 0/1/2/3/4 alignment     //
                //******************************************//
                // pack the reads which are not paired
                packUnPairedReads ( queries, readIDs, readLengths, unAlignedReads,
                                    wordPerQuery, numOfUnPaired, maxBatchSize );
                numQueries = numOfUnPaired;
                all_valid_alignment ( queries, readLengths, NULL, numMismatch, wordPerQuery,
                                      maxBatchSize, numQueries, accumReadNum,
                                      index, _bwt, _revBwt, _occ, _revOcc,
                                      ini_params, input_options,
                                      upkdQualities,
                                      unAlignedReads, numOfUnPaired,
                                      readIDs, upkdQueryNames,
                                      currOutputFileName, currSamOutputFilePtr,
                                      numOfAnswer3, numOfAlignedRead3, 1,
                                      readInputForDP, readInputForNewDP, otherSoap3Result, bothUnalignedPairs );
                // printf("numOfAlignedRead3 = %u\n", numOfAlignedRead3);
                // printf("numOfUnPaired = %u\n", numOfUnPaired)read;
            }
        }
    }

    numOfAnswer = numOfAnswer0 + numOfAnswer1 + numOfAnswer2 + numOfAnswer3;
    numOfAlignedRead = numOfAlignedRead0 + numOfAlignedRead1 + numOfAlignedRead2 + numOfAlignedRead3;
}

// pair-end all-best alignment
// 3-phases [1,2,4]
void all_best_alignment ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                          ullint maxBatchSize, uint numQueries, uint accumReadNum,
                          Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                          IniParams ini_params, InputOptions input_options,
                          char * upkdQualities,
                          uint * unAlignedReads, uint & numOfUnPaired,
                          uint * readIDs, char * upkdQueryNames,
                          char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                          unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                          ReadInputForDPArrays * readInputForDPall,
                          ReadInputForDPArrays * readInputForNewDPall,
                          ReadInputForDPArrays * otherSoap3Resultall,
                          BothUnalignedPairsArrays * bothUnalignedPairsArrays )
{
    // for performing pair-end alignment with the-best and the-second-best output
    unsigned long long numOfAnswer1 = 0;
    uint numOfAlignedRead1 = 0;
    unsigned long long numOfAnswer2 = 0;
    uint numOfAlignedRead2 = 0;
    unsigned long long numOfAnswer3 = 0;
    uint numOfAlignedRead3 = 0;
    unsigned long long numOfAnswer4 = 0;
    uint numOfAlignedRead4 = 0;
    uint numMismatchPerform = 1;
    //********************************//
    // Phase 1: Perform 0/1 Alignment //
    //********************************//
    all_valid_alignment ( queries, readLengths, NULL, numMismatchPerform, wordPerQuery,
                          maxBatchSize, numQueries, accumReadNum,
                          index, _bwt, _revBwt, _occ, _revOcc,
                          ini_params, input_options,
                          upkdQualities,
                          unAlignedReads, numOfUnPaired,
                          readIDs, upkdQueryNames,
                          currOutputFileName,  currSamOutputFilePtr,
                          numOfAnswer1, numOfAlignedRead1, numMismatch <= 1,
                          readInputForDPall->inputArrays, readInputForNewDPall->inputArrays,
                          otherSoap3Resultall->inputArrays,
                          bothUnalignedPairsArrays->array );

    // printf("numOfAlignedRead1 = %u\n", numOfAlignedRead1);
    // printf("numOfUnPaired = %u\n", numOfUnPaired);

    if ( numMismatch > 1 && numOfAlignedRead1 < numQueries )
    {
        //****************************************//
        // Phase 2: Perform 0/1/2 alignment     //
        //****************************************//
        numMismatchPerform = 2;
        // pack the reads which are not paired
        packUnPairedReads ( queries, readIDs, readLengths, unAlignedReads,
                            wordPerQuery, numOfUnPaired, maxBatchSize );
        numQueries = numOfUnPaired;
        all_valid_alignment ( queries, readLengths, NULL, numMismatchPerform, wordPerQuery,
                              maxBatchSize, numQueries, accumReadNum,
                              index, _bwt, _revBwt, _occ, _revOcc,
                              ini_params, input_options,
                              upkdQualities,
                              unAlignedReads, numOfUnPaired,
                              readIDs, upkdQueryNames,
                              currOutputFileName, currSamOutputFilePtr,
                              numOfAnswer2, numOfAlignedRead2, numMismatch == 2,
                              readInputForDPall->inputArrays, readInputForNewDPall->inputArrays,
                              otherSoap3Resultall->inputArrays,
                              bothUnalignedPairsArrays->array );

        // printf("numOfAlignedRead2 = %u\n", numOfAlignedRead2);
        // printf("numOfUnPaired = %u\n", numOfUnPaired, 0);

        if ( numMismatch > 2 && numOfAlignedRead2 < numQueries )
        {
            //******************************************//
            // Phase 3: Perform 0/1/2/3/4 alignment     //
            //******************************************//
            numMismatchPerform = numMismatch;
            // pack the reads which are not paired
            packUnPairedReads ( queries, readIDs, readLengths, unAlignedReads,
                                wordPerQuery, numOfUnPaired, maxBatchSize );
            numQueries = numOfUnPaired;
            all_valid_alignment ( queries, readLengths, NULL, numMismatchPerform, wordPerQuery,
                                  maxBatchSize, numQueries, accumReadNum,
                                  index, _bwt, _revBwt, _occ, _revOcc,
                                  ini_params, input_options,
                                  upkdQualities,
                                  unAlignedReads, numOfUnPaired,
                                  readIDs, upkdQueryNames,
                                  currOutputFileName, currSamOutputFilePtr,
                                  numOfAnswer3, numOfAlignedRead3, 1,
                                  readInputForDPall->inputArrays, readInputForNewDPall->inputArrays,
                                  otherSoap3Resultall->inputArrays,
                                  bothUnalignedPairsArrays->array );
            // printf("numOfAlignedRead3 = %u\n", numOfAlignedRead3);
            // printf("numOfUnPaired = %u\n", numOfUnPaired);
        }
    }

    /*
    // for those reads which needs to further process to get the second-best hits
    if (numOfUnPaired > 0 && numMismatchPerform < 4) {

          numMismatchPerform++;
          uint numOfReadToProcess = numOfUnPaired;

          // repack the reads
          // no read will be removed, but
          // the reads which need to be processed in next-round by soap3 will be duplicated
          // to the front of the list. The readIDs are stored inside the array called "needProcessPair"
          // the corresponding readIDs inside "readInputForDP", "readInputForNewDP" and
          // "bothUnalignedPairs" need to be updated correspondingly.

          printf("Start repacking the reads....\n");
          printf("numOfReadToProcess = %u\n", numOfReadToProcess);
          repackUnPairedReads(&queries, &readIDs, &readLengths, unAlignedReads,
                              wordPerQuery, numOfReadToProcess, numQueries,
                              readInputForDPall, readInputForNewDPall,
                              bothUnalignedPairsArrays);

          printf("Finish repacking the reads.\n");

          all_valid_alignment(queries, readLengths, NULL, numMismatchPerform, wordPerQuery,
                 maxBatchSize, numOfReadToProcess, accumReadNum,
                 index, _bwt, _revBwt, _occ, _revOcc,
                 ini_params, input_options,
                 upkdQualities,
                 unAlignedReads, numOfUnPaired,
                 readIDs, upkdQueryNames,
                 currOutputFileName, currSamOutputFilePtr,
                 numOfAnswer4, numOfAlignedRead4, 1,
                 readInputForDPall->inputArrays, readInputForNewDPall->inputArrays,
                 otherSoap3Resultall->inputArrays,
                 bothUnalignedPairsArrays->array);

          printf("numOfAlignedRead4 = %u\n", numOfAlignedRead4);
          printf("numOfUnPaired = %u\n", numOfUnPaired);
    }
    */
    numOfAnswer = numOfAnswer1 + numOfAnswer2 + numOfAnswer3 + numOfAnswer4;
    numOfAlignedRead = numOfAlignedRead1 + numOfAlignedRead2 + numOfAlignedRead3 + numOfAlignedRead4;
}


void best_single_alignment ( uint * queries, uint * readLengths, uint * seedLengths, uint numMismatch, uint wordPerQuery,
                             ullint maxBatchSize, uint numQueries, uint accumReadNum,
                             Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                             IniParams ini_params, InputOptions input_options,
                             char * upkdQualities,
                             uint * unAlignedReads, uint & numOfUnAligned,
                             uint * readIDs, char * upkdQueryNames,
                             char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                             unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                             ReadInputForDP ** readInputForDP,
                             BothUnalignedPairs ** bothUnalignedPairs )
{
    // for performing single-end alignment
    // random-best, unique-best
    numOfAnswer = 0;
    numOfAlignedRead = 0;

    for ( uint currNumMismatch = 0; currNumMismatch <= numMismatch; currNumMismatch++ )
    {
        unsigned long long currNumOfAnswer = 0;
        uint currNumOfAlignedRead = 0;
        all_valid_alignment ( queries, readLengths, seedLengths, currNumMismatch, wordPerQuery,
                              maxBatchSize, numQueries, accumReadNum,
                              index, _bwt, _revBwt, _occ, _revOcc,
                              ini_params, input_options,
                              upkdQualities,
                              unAlignedReads, numOfUnAligned,
                              readIDs, upkdQueryNames,
                              currOutputFileName,  currSamOutputFilePtr,
                              currNumOfAnswer, currNumOfAlignedRead, currNumMismatch == numMismatch,
                              readInputForDP, NULL, NULL, bothUnalignedPairs );
        numOfAnswer += currNumOfAnswer;
        numOfAlignedRead += currNumOfAlignedRead;

        if ( currNumMismatch < numMismatch )
        {
            // need to proceed the next round
            // pack the reads which have no hits
            packUnPairedReads ( queries, readIDs, readLengths, unAlignedReads,
                                wordPerQuery, numOfUnAligned, maxBatchSize );
            numQueries = numOfUnAligned;
        }
    }
}

void all_best_single_alignment ( uint * queries, uint * readLengths, uint * seedLengths, uint numMismatch, uint wordPerQuery,
                                 ullint maxBatchSize, uint numQueries, uint accumReadNum,
                                 Soap3Index * index, uint * _bwt, uint * _revBwt, uint * _occ, uint * _revOcc,
                                 IniParams ini_params, InputOptions input_options,
                                 char * upkdQualities,
                                 uint * unAlignedReads, uint & numOfUnAligned,
                                 uint * readIDs, char * upkdQueryNames,
                                 char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                                 unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                                 ReadInputForDP ** readInputForDP,
                                 BothUnalignedPairs ** bothUnalignedPairs )
{
    // for performing single-end alignment
    // all-best
    numOfAnswer = 0;
    numOfAlignedRead = 0;

    for ( uint currNumMismatch = 1; currNumMismatch <= ( numMismatch > 0 ? numMismatch : 1 ); currNumMismatch++ )
    {
        unsigned long long currNumOfAnswer = 0;
        uint currNumOfAlignedRead = 0;
        all_valid_alignment ( queries, readLengths, seedLengths, currNumMismatch, wordPerQuery,
                              maxBatchSize, numQueries, accumReadNum,
                              index, _bwt, _revBwt, _occ, _revOcc,
                              ini_params, input_options,
                              upkdQualities,
                              unAlignedReads, numOfUnAligned,
                              readIDs, upkdQueryNames,
                              currOutputFileName,  currSamOutputFilePtr,
                              currNumOfAnswer, currNumOfAlignedRead, currNumMismatch >= numMismatch,
                              readInputForDP, NULL, NULL, bothUnalignedPairs );
        numOfAnswer += currNumOfAnswer;
        numOfAlignedRead += currNumOfAlignedRead;

        if ( currNumMismatch < numMismatch )
        {
            // need to proceed the next round
            // pack the reads which have no hits
            packUnPairedReads ( queries, readIDs, readLengths, unAlignedReads,
                                wordPerQuery, numOfUnAligned, maxBatchSize );
            numQueries = numOfUnAligned;
        }
    }
}


// Perform single-read all-valid alignment for seeds
// This function would NOT reset the SingleAlgnResultArray
void single_all_valid_seed_alignment (
    unsigned int * queries, unsigned int * readLengths, uint numMismatch,
    unsigned int maxReadLength, unsigned int wordPerQuery, unsigned int maxBatchSize,
    unsigned int numQueries,
    Soap3Index * index,
    uint * _bwt, unsigned int * _revBwt,
    unsigned int * _occ, unsigned int * _revOcc,
    unsigned long long & numOfAnswer,
    unsigned int & numOfAlignedRead, unsigned int cpuNumThreads,
    SingleAlgnResultArray * alignResultArray, unsigned int maxHitNum,
    unsigned char outNoAlign, unsigned char * noAlignment,
    unsigned int * readIDs )
{
#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    double startTime2 = setStartTime ();
    double lastEventTime2 = 0;
    double currEventTime2;
#endif
    uint * answers[MAX_NUM_CASES];
    uint * badReadIndices[MAX_NUM_CASES]; // IDs of bad reads for each case
    uint * badAnswers[MAX_NUM_CASES]; // stores answers (and reads temporarily) for 2nd round
    uint threadBadStarts[MAX_NUM_CPU_THREADS][MAX_NUM_CASES];
    uint threadBadCounts[MAX_NUM_CPU_THREADS][MAX_NUM_CASES];
    // Host alignment model
    SRAModel * SRAMismatchModel[MAX_READ_LENGTH + 1];
    SRAModel * SRAMismatchModel_neg[MAX_READ_LENGTH + 1];
    HostKernelArguements hostKernelArguments[MAX_NUM_CPU_THREADS];
    // Host multi-threading variables
    int threadId;
    pthread_t threads[MAX_NUM_CPU_THREADS];
    uint threadBucketSize[MAX_NUM_CPU_THREADS];
    numOfAnswer = 0;
    numOfAlignedRead = 0;
    uint blocksNeeded;
    ullint batchSize;
    uint numCases;
    uint sa_range_allowed_1;
    uint sa_range_allowed_2;
    char skip_round_2;
    ullint queriesLeft = numQueries;
    uint word_per_ans;
    uint word_per_ans_2;
    int i;

    if ( numQueries == 0 )
    { return; }

    SRAIndex * sraIndex = index->sraIndex;

    // initialization of arrays
    for ( int j = 0; j < MAX_NUM_CASES; j++ )
    {
        badReadIndices[j] = NULL;
        badAnswers[j] = NULL;
        answers[j] = NULL;

        for ( int k = 0; k < MAX_NUM_CPU_THREADS; k++ )
        {
            threadBadStarts[k][j] = 0;
            threadBadCounts[k][j] = 0;
        }
    }

    // obtain the number of cases for this number of mismatch
    // and obtain the number of SA ranges allowed
    getParametersForThisMismatch2 ( numMismatch, numCases, sa_range_allowed_1,
                                    sa_range_allowed_2, skip_round_2, word_per_ans, word_per_ans_2 );
    // set the multi-threading arguments
    setHostKernelArguments2 ( hostKernelArguments, threads, index, maxReadLength, wordPerQuery, word_per_ans, cpuNumThreads );

    for ( threadId = 0; threadId < cpuNumThreads; ++threadId )
    {
        hostKernelArguments[threadId].SRAMismatchModel = SRAMismatchModel;
        hostKernelArguments[threadId].SRAMismatchModel_neg = SRAMismatchModel_neg;
        hostKernelArguments[threadId].readLengths = readLengths;
        hostKernelArguments[threadId].maxHitNum = maxHitNum;
        hostKernelArguments[threadId].readIDs = readIDs;
    }

    ullint roundUp = ( maxBatchSize + 31 ) / 32 * 32;

    // allocate host memory for answers
    for ( int j = 0; j < numCases; j++ )
    {
        answers[j] = ( uint * ) malloc ( roundUp * word_per_ans * sizeof ( uint ) );
    }

    // printParameters(input_options, ini_params);
    // printf("numMismatch = %u; numCases = %u; sa_range_allowed_1 = %u; sa_range_allowed_2 = %u; skip_round_2 = %i; word_per_ans = %u; word_per_ans_2 = %u \n\n\n\n\n", 1, numCases, sa_range_allowed_1, sa_range_allowed_2, skip_round_2, word_per_ans, word_per_ans_2);

    // update the num of mismatch of the model for host alignment
    for ( threadId = 0; threadId < cpuNumThreads; ++threadId )
    {
        hostKernelArguments[threadId].sraQuerySettings.MaxError = numMismatch;
        hostKernelArguments[threadId].sraQuerySettings.MaxNBMismatch = 0;
        hostKernelArguments[threadId].numCases = numCases;
        hostKernelArguments[threadId].sa_range_allowed_1 = sa_range_allowed_1;
        hostKernelArguments[threadId].sa_range_allowed_2 = sa_range_allowed_2;
        hostKernelArguments[threadId].maxNumMismatch = numMismatch;
        hostKernelArguments[threadId].alignResult = alignResultArray->array[threadId];
    }

    for ( i = 0; i <= MAX_READ_LENGTH; i++ )
    {
        SRAMismatchModel[i] = NULL;
        SRAMismatchModel_neg[i] = NULL;
    }

    for ( unsigned readId = 0; readId < numQueries; ++readId )
    {
        if ( SRAMismatchModel[readLengths[readId]] == NULL )
        {
            SRAMismatchModel[readLengths[readId]] = SRAModelConstruct ( readLengths[readId], QUERY_POS_STRAND, & ( hostKernelArguments[0].sraQuerySettings ), sraIndex, SRA_MODEL_16G );
        }

        if ( SRAMismatchModel_neg[readLengths[readId]] == NULL )
        {
            SRAMismatchModel_neg[readLengths[readId]] = SRAModelConstruct ( readLengths[readId], QUERY_NEG_STRAND, & ( hostKernelArguments[0].sraQuerySettings ), sraIndex, SRA_MODEL_16G );
        }
    }

    uint * nextQuery = queries;
    uint * nextReadLength = readLengths;
#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
    currEventTime2 = getElapsedTime ( startTime2 );
    printf ( "[Main] Time elapsed for initialization: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
    lastEventTime2 = currEventTime2;
#endif

    while ( queriesLeft > 0 )
    {
        uint processedQuery = numQueries - queriesLeft;

        if ( queriesLeft > maxBatchSize )
        {
            blocksNeeded = NUM_BLOCKS;
            batchSize = maxBatchSize;
        }
        else
        {
            blocksNeeded = ( queriesLeft + THREADS_PER_BLOCK * QUERIES_PER_THREAD - 1 ) /
                           ( THREADS_PER_BLOCK * QUERIES_PER_THREAD );
            batchSize = queriesLeft;
        }

        // allocate the reads to different CPU threads
        uint roughBucketSize = batchSize / cpuNumThreads;
        uint batchSizeUnalloc = batchSize;

        for ( threadId = 0; threadId < cpuNumThreads - 1; threadId++ )
        {
            threadBucketSize[threadId] = roughBucketSize;
            batchSizeUnalloc -= roughBucketSize;
        }

        threadBucketSize[cpuNumThreads - 1] = batchSizeUnalloc;
        // perform first round alignment in GPU
        perform_round1_alignment_no_pipeline ( nextQuery, nextReadLength, answers,
                                               numMismatch, numCases, sa_range_allowed_1, wordPerQuery,
                                               word_per_ans, false, blocksNeeded, batchSize,
                                               index, _bwt, _revBwt, _occ, _revOcc );
#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
        printf ( "[Main] # of mismatches allowed: %u with # of sa ranges: %u\n", 1, sa_range_allowed_1 );
        currEventTime2 = getElapsedTime ( startTime2 );
        printf ( "[Main] Time elapsed for first round alignment in GPU: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
        lastEventTime2 = currEventTime2;
#endif

        // identify which read has no hit
        if ( outNoAlign == 1 )
        {
            uint * ans;

            for ( ullint j = 0; j < batchSize; ++j )
            {
                ullint batchReadId = j;
                ullint readId = processedQuery + j;
                char isNoAlign = 1;

                for ( int whichCase = 0; whichCase < numCases && isNoAlign; whichCase++ )
                {
                    ans = answers[whichCase] + ( ( batchReadId ) / 32 * 32 * word_per_ans + ( batchReadId ) % 32 );

                    if ( * ( ans ) != 0xFFFFFFFD )
                    { isNoAlign = 0; }
                }

                if ( isNoAlign == 1 )
                {
                    noAlignment[readId] = 1;
                }
            }
        }

        if ( skip_round_2 == 0 )
        {
            // =======================================
            // | GPU-2                               |
            // =======================================
            // perform second round alignment in GPU
            for ( int whichCase = 0; whichCase < numCases; ++whichCase )
            {
                ullint numBads = 0;

                // Count number of bad reads and
                for ( ullint readId = 0; readId < batchSize; readId++ )
                {
                    ullint srcOffset = ( ( readId ) / 32 * 32 * word_per_ans + readId % 32 );
                    numBads += ( answers[whichCase][srcOffset] > 0xFFFFFFFD );
                }

                if ( numBads > 0 )
                {
                    // Allocate memory and copy bad reads to another array
                    ullint roundUp = ( numBads + 31 ) / 32 * 32;
                    badReadIndices[whichCase] = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
                    // printf("size of badReadIndices[%i] = %u\n",whichCase, roundUp);
                    badAnswers[whichCase] = ( uint * ) malloc ( roundUp * word_per_ans_2 * sizeof ( uint ) );
                    // printf("size of badAnswers[%i] = %u\n", whichCase, roundUp * word_per_ans_2);
                }

                // printf("numBads == %u in case %i of mismatch %u\n", numBads, whichCase, numMismatch);
            }

            perform_round2_alignment_no_pipeline ( queries, readLengths, answers,
                                                   numMismatch, numCases, sa_range_allowed_2, wordPerQuery, word_per_ans,
                                                   word_per_ans_2, false, blocksNeeded, batchSize,
                                                   index, _bwt, _revBwt, _occ, _revOcc,
                                                   processedQuery, badReadIndices, badAnswers );

            for ( int whichCase = 0; whichCase < numCases; ++whichCase )
            {
                ullint readId = 0;

                for ( threadId = 0; threadId < cpuNumThreads; threadId ++ )
                {
                    uint threadNumBad = 0;

                    // printf("whichCase = %i; threadId = %i; threadBucketSize[threadId] = %u\n",
                    //        whichCase, threadId, threadBucketSize[threadId]);
                    for ( uint n = 0 ; n < threadBucketSize[threadId]; ++n )
                    {
                        ullint srcOffset = ( readId / 32 * 32 * word_per_ans + readId % 32 );
                        threadNumBad += ( answers[whichCase][srcOffset] > 0xFFFFFFFD );
                        readId++;
                    }

                    threadBadCounts[threadId][whichCase] = threadNumBad;
                    // printf("threadBadCounts[%i][%i] = %u\n",
                    //    threadId, whichCase, threadBadCounts[threadId][whichCase]);
                }
            }

            for ( int whichCase = 0; whichCase < numCases; ++whichCase )
            {
                threadBadStarts[0][whichCase] = 0;

                for ( threadId = 1; threadId < cpuNumThreads; threadId ++ )
                {
                    threadBadStarts[threadId][whichCase] = threadBadStarts[threadId - 1][whichCase] + threadBadCounts[threadId - 1][whichCase];
                    // printf("threadBadStarts[%i][%i] = %u\n",
                    //  threadId, whichCase, threadBadStarts[threadId][whichCase]);
                }
            }

#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
            currEventTime2 = getElapsedTime ( startTime2 );
            printf ( "[Main] Time elapsed for second round alignment in GPU: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
            printf ( "[Main] # of sa ranges allowed: %u \n", sa_range_allowed_2 );
            lastEventTime2 = currEventTime2;
#endif
        } // skip the round 2

#ifdef BGS_GPU_CASE_BREAKDOWN_TIME
        currEventTime2 = getElapsedTime ( startTime2 );
        printf ( "[Main] Time elapsed for waiting for CPU threads: %9.4f seconds.\n", currEventTime2 - lastEventTime2 );
        lastEventTime2 = currEventTime2;
#endif
        // ================ PIPE SECTION FINAL========================
        // =======================================
        // | CPU: Thread #0,1,2,3,...            |
        // =======================================
        unsigned int threadProcessedQuery = 0;

        for ( threadId = 0; threadId < cpuNumThreads ; threadId++ )
        {
            hostKernelArguments[threadId].batchFirstReadId = processedQuery;
            hostKernelArguments[threadId].skipFirst = threadProcessedQuery;
            hostKernelArguments[threadId].numQueries = threadBucketSize[threadId];
            hostKernelArguments[threadId].queries = queries;
            hostKernelArguments[threadId].word_per_query = wordPerQuery;
            hostKernelArguments[threadId].answers = answers;
            hostKernelArguments[threadId].alignedOcc = 0;
            hostKernelArguments[threadId].alignedReads = 0;
            hostKernelArguments[threadId].badReadIndices = badReadIndices;
            hostKernelArguments[threadId].badAnswers = badAnswers;
            hostKernelArguments[threadId].badStartOffset = threadBadStarts[threadId];
            hostKernelArguments[threadId].badCountOffset = threadBadCounts[threadId];
            hostKernelArguments[threadId].skip_round_2 = skip_round_2;

            if ( pthread_create ( & ( threads[threadId] ), NULL, hostKernelThreadWrapperSingle, ( void * ) & ( hostKernelArguments[threadId] ) ) )
            { fprintf ( stderr, "[Main:Threads%u] Can't create hostKernelThreadWrapper\n", threadId ), exit ( 1 ); }

            threadProcessedQuery += threadBucketSize[threadId];
        }

        for ( threadId = 0; threadId < cpuNumThreads ; threadId++ )
        {
            if ( threads[threadId] != 0 )
            {
                if ( pthread_join ( threads[threadId], NULL ) )
                { fprintf ( stderr, "[Main:Thread%u] Crash!\n", threadId ), exit ( 1 ); }

                numOfAnswer += hostKernelArguments[threadId].alignedOcc;
                numOfAlignedRead += hostKernelArguments[threadId].alignedReads;
                threads[threadId] = 0;
            }
        }

        for ( int whichCase = 0; whichCase < numCases; ++whichCase )
        {
            if ( badReadIndices[whichCase] != NULL )
            {
                free ( badReadIndices[whichCase] );
                badReadIndices[whichCase] = NULL;
            }

            if ( badAnswers[whichCase] != NULL )
            {
                free ( badAnswers[whichCase] );
                badAnswers[whichCase] = NULL;
            }
        }

        // ================ PIPE SECTION FINAL========================
        queriesLeft -= batchSize;
        nextQuery += batchSize * wordPerQuery;
        nextReadLength += batchSize;
    }

    // CLEAN UP for each MISMATCH iteration                                                                          |
    for ( i = 0; i <= MAX_READ_LENGTH; i++ )
    {
        if ( SRAMismatchModel[i] != NULL )
        {
            SRAModelFree ( SRAMismatchModel[i] );
            SRAMismatchModel[i] = NULL;
        }

        if ( SRAMismatchModel_neg[i] != NULL )
        {
            SRAModelFree ( SRAMismatchModel_neg[i] );
            SRAMismatchModel_neg[i] = NULL;
        }
    }

    for ( threadId = 0; threadId < cpuNumThreads; ++threadId )
    {
        OCCFree ( hostKernelArguments[threadId].occ );
    }

    for ( int j = 0; j < numCases; j++ )
    {
        free ( answers[j] );
    }
}

// Perform single-read all-best 1-mismatch alignment for seeds
// This function would reset the SingleAlgnResultArray
void single_1_mismatch_alignment2 ( unsigned int * queries, unsigned int * readLengths,
                                    unsigned int maxReadLength, unsigned int wordPerQuery, unsigned int maxBatchSize,
                                    unsigned int numQueries, Soap3Index * index, uint * _bwt, unsigned int * _revBwt,
                                    unsigned int * _occ, unsigned int * _revOcc,
                                    unsigned long long & numOfAnswer,
                                    unsigned int & numOfAlignedRead, unsigned int cpuNumThreads,
                                    SingleAlgnResultArray * alignResultArray, unsigned int maxHitNum )
{
    // reset the array alignResult
    for ( unsigned int threadId = 0; threadId < cpuNumThreads; ++threadId )
    {
        resetSingleAlgnResult ( alignResultArray->array[threadId] );
    }

    // create the arrays for processing
    unsigned char * noAlignment = ( unsigned char * ) malloc ( sizeof ( unsigned char ) * numQueries );
    memset ( noAlignment, 0, numQueries );
    single_all_valid_seed_alignment (
        queries, readLengths, 0,
        maxReadLength, wordPerQuery, maxBatchSize,
        numQueries, index, _bwt, _revBwt,
        _occ, _revOcc,
        numOfAnswer,
        numOfAlignedRead, cpuNumThreads,
        alignResultArray, maxHitNum,
        1, noAlignment, NULL );
    // create the arrays for processing
    unsigned long long roundUp = ( numQueries + 31 ) / 32 * 32;
    unsigned int * queries2 = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * roundUp * wordPerQuery );
    unsigned int * readLengths2 = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * roundUp );
    memcpy ( queries2, queries, sizeof ( unsigned int ) * roundUp * wordPerQuery );
    memcpy ( readLengths2, readLengths, sizeof ( unsigned int ) * numQueries );
    // pack the unaligned reads
    unsigned int numUnAligned;
    unsigned int * readIDs2 = packReads2 ( queries2, readLengths2, noAlignment,
                                           wordPerQuery, numQueries, numUnAligned );
    unsigned int numOfAlignedRead2;
    unsigned long long numOfAnswer2;
    single_all_valid_seed_alignment (
        queries2, readLengths2, 1,
        maxReadLength, wordPerQuery, maxBatchSize,
        numUnAligned, index, _bwt, _revBwt,
        _occ, _revOcc,
        numOfAnswer2,
        numOfAlignedRead2, cpuNumThreads,
        alignResultArray, maxHitNum,
        0, noAlignment, readIDs2 );
    numOfAnswer += numOfAnswer2;
    numOfAlignedRead += numOfAlignedRead2;
    free ( queries2 );
    free ( readLengths2 );
    free ( readIDs2 );
    free ( noAlignment );
}


// Perform SOAP3-DP Paired-End Alignment
void soap3_dp_pair_align ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                           ullint maxBatchSize, uint numQueries, uint accumReadNum,
                           Soap3Index * index,
                           uint * _bwt, uint * _revBwt,
                           uint * _occ, uint * _revOcc,
                           IniParams ini_params, InputOptions input_options,
                           uint maxReadLength, uint detected_read_length, uint detected_read_length2,
                           char * upkdQualities,
                           uint * unAlignedReads, uint & numOfUnPaired,
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
                           uint & indexLoadedToGPU )
{
    double alignmentTime, copyTime;
    // ======================================================================================
    // | IF THE INDEX IS NOT IN DEVICE, THEN                                                |
    // | COPY INDEX TO DEVICE MEMORY                                                        |
    // ======================================================================================
    unsigned int numDPAlignedPair = 0;
    unsigned int numDPAlignment = 0;
    DPParameters dpParameters;
    int orig_align_type;
    HSPAux * hspaux = index->sraIndex->hspaux;

    if ( ini_params.Ini_skipSOAP3Alignment == 1 )
    {
        // Add all first reads of the pairs to BothUnalignedPairs (i.e. 0, 2, 4, ..., totalReadNum-2)
        addAllFirstReadIDToBothUnalignedPairs ( bothUnalignedPairsArrays->array[0], numQueries );
        // For DP module, if SAM format and all-best are both selected, then
        // output format is needed to set to all-valid.
        orig_align_type = input_options.alignmentType;

        if ( input_options.alignmentType == OUTPUT_ALL_BEST &&
                input_options.outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            input_options.alignmentType = OUTPUT_ALL_VALID;
        }

        // Parameters for DP
        getParameterForAllDP ( dpParameters, input_options, ini_params );
    }
    else
    {
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

        //*******************************//
        // Perform Alignment             //
        //*******************************//

        if ( input_options.alignmentType == OUTPUT_ALL_BEST &&
                input_options.outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            all_best_alignment ( queries, readLengths, input_options.numMismatch, wordPerQuery,
                                 maxBatchSize, numQueries, accumReadNum,
                                 index,
                                 _bwt, _revBwt, _occ, _revOcc,
                                 ini_params, input_options,
                                 upkdQualities,
                                 unAlignedReads, numOfUnPaired,
                                 readIDs, upkdQueryNames,
                                 currOutputFileName, currSamOutputFilePtr,
                                 numOfAnswer, numOfAlignedRead,
                                 readInputForDPall, readInputForNewDPall,
                                 otherSoap3Resultall,
                                 bothUnalignedPairsArrays );
        }
        else if ( input_options.alignmentType == OUTPUT_UNIQUE_BEST ||
                  input_options.alignmentType == OUTPUT_RANDOM_BEST ||
                  input_options.alignmentType == OUTPUT_ALL_BEST )
        {
            four_phases_alignment ( queries, readLengths, input_options.numMismatch, wordPerQuery,
                                    maxBatchSize, numQueries, accumReadNum,
                                    index,
                                    _bwt, _revBwt, _occ, _revOcc,
                                    ini_params, input_options,
                                    upkdQualities,
                                    unAlignedReads, numOfUnPaired,
                                    readIDs, upkdQueryNames,
                                    currOutputFileName, currSamOutputFilePtr,
                                    numOfAnswer, numOfAlignedRead,
                                    readInputForDPall->inputArrays, readInputForNewDPall->inputArrays,
                                    otherSoap3Resultall->inputArrays,
                                    bothUnalignedPairsArrays->array );
        }
        else
        {
            all_valid_alignment ( queries, readLengths, NULL, input_options.numMismatch, wordPerQuery,
                                  maxBatchSize, numQueries, accumReadNum,
                                  index,
                                  _bwt, _revBwt, _occ, _revOcc,
                                  ini_params, input_options,
                                  upkdQualities,
                                  unAlignedReads, numOfUnPaired,
                                  readIDs, upkdQueryNames,
                                  currOutputFileName, currSamOutputFilePtr,
                                  numOfAnswer, numOfAlignedRead, 1,
                                  readInputForDPall->inputArrays, readInputForNewDPall->inputArrays,
                                  otherSoap3Resultall->inputArrays,
                                  bothUnalignedPairsArrays->array );
        }

        printf ( "[Main] Finished alignment with <= %i mismatches\n", input_options.numMismatch );
        // printf("[Main] Number of pairs aligned: %u (number of alignments: %llu)\n", numOfAlignedRead/2, numOfAnswer);
        printf ( "[Main] Number of pairs aligned: %u\n", numOfAlignedRead / 2 );
        alignmentTime = getElapsedTime ( startTime );
        printf ( "[Main] Elapsed time : %9.4f seconds\n\n", alignmentTime - lastEventTime );
        totalAlignmentTime += alignmentTime - lastEventTime;
        lastEventTime = alignmentTime;
        // For DP module, if SAM format and all-best are both selected, then
        // output format is needed to set to all-valid.
        orig_align_type = input_options.alignmentType;

        if ( input_options.alignmentType == OUTPUT_ALL_BEST &&
                input_options.outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            input_options.alignmentType = OUTPUT_ALL_VALID;
        }

        numDPAlignedPair = 0;
        numDPAlignment = 0;
        // Parameters for DP
        getParameterForAllDP ( dpParameters, input_options, ini_params );
        
        numDPAlignedPair = 0;
        numDPAlignment = 0;

        ////////////////////////////////////////////////////////////
        // PERFORM NEW SEMI-GLOBAL DP IF NECESSARY                    //
        ////////////////////////////////////////////////////////////

#ifdef PERFORM_NEW_DEFAULT_DP_FOR_SEMI_ALIGNED_PAIR
        if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ )
        {
            // DP Parameters for half-aligned reads
            getParameterForNewDefaultDP ( dpParameters, input_options, ini_params, detected_read_length, detected_read_length2 );
            printDPParameters ( dpParameters );
            /*
            if (indexLoadedToGPU == 1) {
                  // free device memory
                  // printf("[Main] Free index from device memory..\n");
                  cudaFree(_bwt);
                  cudaFree(_occ);
                  cudaFree(_revBwt);
                  cudaFree(_revOcc);


                  indexLoadedToGPU=0;
            }
            */
#ifdef BGS_OUTPUT_DP_MESSAGE
            printf ( "*********************************************************\n" );
            printf ( "NEW SEMI-GLOBAL MESSAGE:\n" );
#endif
            unsigned int totalReadsProceedToDP = 0;
            unsigned int totalSARanges = 0;
            unsigned int totalOccs = 0;
            unsigned int totalHits = 0;

            for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
            {
                totalReadsProceedToDP += readInputForNewDPall->inputArrays[threadId]->readNum;
                totalSARanges += readInputForNewDPall->inputArrays[threadId]->saRangeTotalNum;
                totalOccs += readInputForNewDPall->inputArrays[threadId]->occTotalNum;

                for ( int h = 0; h < readInputForNewDPall->inputArrays[threadId]->saRangeTotalNum; h++ )
                { totalHits += readInputForNewDPall->inputArrays[threadId]->sa_list[h].saIndexRight - readInputForNewDPall->inputArrays[threadId]->sa_list[h].saIndexLeft + 1; }

                totalHits += readInputForNewDPall->inputArrays[threadId]->occTotalNum;
            }

            if ( totalReadsProceedToDP > 0 )
            {
                printf ( "[Main] %u half-aligned pairs of reads are proceeded to new default DP.\n", totalReadsProceedToDP );
#ifdef BGS_OUTPUT_HALF_ALIGNED_DP_MESSAGE
                printf ( "%u half-aligned pairs of reads (with %u SA Ranges and \n", totalReadsProceedToDP, totalSARanges );
                printf ( "%u occurrences) are proceeded to new default DP.\n", totalOccs );
                printf ( "total %u number of hits) are proceeded to new default DP.\n", totalHits );
#endif
                // prepare the file pointer for output
                FILE * outputDPFile = NULL;

                //Output file for SAM output is handled by SAM API
                //=> OutFilePtr is not used for SAM API output format.
                switch ( input_options.outputFormat )
                {
                    case SRA_OUTPUT_FORMAT_SAM_API:
                        break;

                    default:
                        outputDPFile = ( FILE * ) fopen ( outputDPFileName, "a" );
                }

                newSemiGlobalDP ( readInputForNewDPall, input_options.insert_high, input_options.insert_low,
                                  queries, readLengths, readIDs, upkdQueryNames, upkdQualities,
                                  maxReadLength, ini_params.Ini_PEStrandLeftLeg, ini_params.Ini_PEStrandRightLeg,
                                  index,
                                  _bwt,  _revBwt,
                                  _occ,  _revOcc,
                                  input_options.alignmentType, &dpParameters,
                                  numDPAlignedPair, numDPAlignment,
                                  accumReadNum, input_options.outputFormat,
                                  outputDPFile, samOutputDPFilePtr, bothUnalignedPairsArrays->array[ini_params.Ini_NumOfCpuThreads] );
                // the index is released in the function "newSemiGlobalDP"
                indexLoadedToGPU = 0;

                switch ( input_options.outputFormat )
                {
                    case SRA_OUTPUT_FORMAT_SAM_API:
                        break;

                    default:
                        fclose ( outputDPFile );
                }

                // printf("Finished semi-global DP\n");
                // printf("[Main] Number of pairs aligned by new default DP: %u (number of alignments: %u)\n", numDPAlignedPair, numDPAlignment);
                printf ( "[Main] Number of pairs aligned by new default DP: %u\n", numDPAlignedPair );
                alignmentTime = getElapsedTime ( startTime );
                printf ( "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
                totalAlignmentTime += alignmentTime - lastEventTime;
                lastEventTime = alignmentTime;
                numOfAlignedRead +=  numDPAlignedPair * 2;
                numOfAnswer += numDPAlignment;
                // printf("[Main] Total Number of pairs aligned: %u (number of alignments: %llu)\n", numOfAlignedRead/2, numOfAnswer);
                printf ( "[Main] Total Number of pairs aligned: %u\n", numOfAlignedRead / 2 );
#ifdef BGS_OUTPUT_DP_MESSAGE
                printf ( "*********************************************************\n" );
#endif
                printf ( "\n" );
            }
        }
#endif

        numDPAlignedPair = 0;
        numDPAlignment = 0;
    }

    ////////////////////////////////////////////////////////////
    // PERFORM SEMI-GLOBAL DP IF NECESSARY                    //
    ////////////////////////////////////////////////////////////

#ifdef PERFORM_DEFAULT_DP_FOR_SEMI_ALIGNED_PAIR
    if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ )
    {
        // DP Parameters for half-aligned reads
        getParameterForDefaultDP ( dpParameters, input_options, ini_params, detected_read_length, detected_read_length2 );
        printDPParameters ( dpParameters );

        if ( indexLoadedToGPU == 1 )
        {
            // free device memory
            // printf("[Main] Free index from device memory..\n");
            cudaFree ( _bwt );
            cudaFree ( _occ );
            cudaFree ( _revBwt );
            cudaFree ( _revOcc );
            indexLoadedToGPU = 0;
        }

#ifdef BGS_OUTPUT_DP_MESSAGE
        printf ( "*********************************************************\n" );
        printf ( "SEMI-GLOBAL MESSAGE:\n" );
#endif
        unsigned int totalReadsProceedToDP = 0;
        unsigned int totalSARanges = 0;
        unsigned int totalOccs = 0;

        for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
        {
            totalReadsProceedToDP += readInputForDPall->inputArrays[threadId]->readNum;
            totalSARanges += readInputForDPall->inputArrays[threadId]->saRangeTotalNum;
            totalOccs += readInputForDPall->inputArrays[threadId]->occTotalNum;
        }

        if ( totalReadsProceedToDP > 0 )
        {
            printf ( "[Main] %u half-aligned pairs of reads are proceeded to DP.\n", totalReadsProceedToDP );
#ifdef BGS_OUTPUT_HALF_ALIGNED_DP_MESSAGE
            printf ( "%u half-aligned pairs of reads (with %u SA Ranges and \n", totalReadsProceedToDP, totalSARanges );
            printf ( "%u occurrences) are proceeded to DP.\n", totalOccs );
#endif
            // prepare the file pointer for output
            FILE * outputDPFile = NULL;

            //Output file for SAM output is handled by SAM API
            //=> OutFilePtr is not used for SAM API output format.
            switch ( input_options.outputFormat )
            {
                case SRA_OUTPUT_FORMAT_SAM_API:
                    break;

                default:
                    outputDPFile = ( FILE * ) fopen ( outputDPFileName, "a" );
            }

            semiGlobalDP2 ( readInputForDPall, input_options.insert_high, input_options.insert_low,
                            queries, readLengths, readIDs, upkdQueryNames, upkdQualities,
                            maxReadLength, ini_params.Ini_PEStrandLeftLeg, ini_params.Ini_PEStrandRightLeg,
                            index,
                            _bwt,  _revBwt,
                            _occ,  _revOcc,
                            input_options.alignmentType, &dpParameters,
                            numDPAlignedPair, numDPAlignment,
                            accumReadNum, input_options.outputFormat,
                            outputDPFile, samOutputDPFilePtr, bothUnalignedPairsArrays->array[ini_params.Ini_NumOfCpuThreads] );

            switch ( input_options.outputFormat )
            {
                case SRA_OUTPUT_FORMAT_SAM_API:
                    break;

                default:
                    fclose ( outputDPFile );
            }

            // printf("Finished semi-global DP\n");
            // printf("[Main] Number of pairs aligned by DP: %u (number of alignments: %u)\n", numDPAlignedPair, numDPAlignment);
            printf ( "[Main] Number of pairs aligned by DP: %u\n", numDPAlignedPair );
            // printf(" Number of reads passed to deep-dp: %u\n", bothUnalignedPairsArrays->array[ini_params.Ini_NumOfCpuThreads]->totalNum);
            alignmentTime = getElapsedTime ( startTime );
            printf ( "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
            totalAlignmentTime += alignmentTime - lastEventTime;
            lastEventTime = alignmentTime;
            numOfAlignedRead +=  numDPAlignedPair * 2;
            numOfAnswer += numDPAlignment;
            // printf("[Main] Total Number of pairs aligned: %u (number of alignments: %llu)\n", numOfAlignedRead/2, numOfAnswer);
            printf ( "[Main] Total Number of pairs aligned: %u\n", numOfAlignedRead / 2 );
#ifdef BGS_OUTPUT_DP_MESSAGE
            printf ( "*********************************************************\n" );
#endif
            printf ( "\n" );
        }
    }

#endif
    ////////////////////////////////////////////////////////////
    // PERFORM DP FOR BOTH ENDS UNALIGNED                     //
    ////////////////////////////////////////////////////////////
    numDPAlignedPair = 0;
    numDPAlignment = 0;

#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_PAIR
    if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ )
    {
        // DP Parameters for Deep DP
        getParameterForDeepDP ( dpParameters, input_options, ini_params, detected_read_length, detected_read_length2 );
        printDPParameters ( dpParameters );

        // ======================================================================================
        // | IF THE INDEX IS NOT IN DEVICE, THEN                                                |
        // | COPY INDEX TO DEVICE MEMORY                                                        |
        // ======================================================================================

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

#ifdef BGS_OUTPUT_DP_MESSAGE
        printf ( "*********************************************************\n" );
        printf ( "DP FOR BOTH UNALIGNED READS\n" );
#endif
        unsigned int totalReadsProceedToDP = 0;

        for ( int threadId = 0; threadId <= ini_params.Ini_NumOfCpuThreads; ++threadId )
        {
            totalReadsProceedToDP += bothUnalignedPairsArrays->array[threadId]->totalNum;
        }

        if ( totalReadsProceedToDP > 0 )
        {
            printf ( "[Main] %u pairs of reads are proceeded to deep DP.\n", totalReadsProceedToDP );
            // print out the read IDs
            // printAllReadIDs(bothUnalignedPairsArrays, accumReadNum);
            // prepare the file pointer for output
            FILE * outputDPFile = NULL;

            //Output file for SAM output is handled by SAM API
            //=> OutFilePtr is not used for SAM API output format.
            switch ( input_options.outputFormat )
            {
                case SRA_OUTPUT_FORMAT_SAM_API:
                    break;

                default:
                    outputDPFile = ( FILE * ) fopen ( outputDPFileName, "a" );
            }

            DPForUnalignPairs2 ( bothUnalignedPairsArrays, input_options.insert_high, input_options.insert_low,
                                 queries, readLengths, readIDs, upkdQueryNames, upkdQualities,
                                 maxReadLength, ini_params.Ini_PEStrandLeftLeg, ini_params.Ini_PEStrandRightLeg,
                                 index,
                                 _bwt, _revBwt,
                                 _occ, _revOcc,
                                 input_options.alignmentType, &dpParameters,
                                 numDPAlignedPair, numDPAlignment,
                                 accumReadNum, input_options.outputFormat,
                                 outputDPFile, samOutputDPFilePtr );
            indexLoadedToGPU = 0;

            switch ( input_options.outputFormat )
            {
                case SRA_OUTPUT_FORMAT_SAM_API:
                    break;

                default:
                    fclose ( outputDPFile );
            }

            // printf("Finished DP for both-end unaligned reads\n");
            // printf("[Main] Number of pairs aligned by DP: %u (number of alignments: %u)\n", numDPAlignedPair, numDPAlignment);
            printf ( "[Main] Number of pairs aligned by DP: %u\n", numDPAlignedPair );
            alignmentTime = getElapsedTime ( startTime );
            printf ( "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
            totalAlignmentTime += alignmentTime - lastEventTime;
            lastEventTime = alignmentTime;
            numOfAlignedRead +=  numDPAlignedPair * 2;
            numOfAnswer += numDPAlignment;
            // printf("[Main] Total Number of pairs aligned: %u (number of alignments: %llu)\n", numOfAlignedRead/2, numOfAnswer);
            printf ( "[Main] Total Number of pairs aligned: %u\n", numOfAlignedRead / 2 );
#ifdef BGS_OUTPUT_DP_MESSAGE
            printf ( "*********************************************************\n" );
#endif
            printf ( "\n" );
        }
    }
#endif

    ////////////////////////////////////////////////////////////
    // TO PERFORM SINGLE DEEP-DP FOR UNALIGNED READS          //
    ////////////////////////////////////////////////////////////
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS
    unsigned int numSingleDPAligned = 0;
    unsigned int numSingleDPAlignment = 0;

    if ( input_options.enableDP == 1 && input_options.readType == PAIR_END_READ && input_options.outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    {
        // reset the readType to single read
        input_options.readType = SINGLE_READ;
        // DP Parameters for single-end alignment
        getParameterForSingleDP ( dpParameters, input_options, ini_params, detected_read_length );
        hspaux->singleDPcutoffThreshold = dpParameters.paramRead[0].cutoffThreshold;
        // printDPParameters(dpParameters);
#ifdef BGS_OUTPUT_DP_MESSAGE
        printf ( "*********************************************************\n" );
        printf ( "DP FOR UNALIGNED SINGLE READS\n" );
#endif
        unsigned int totalReadsProceedToDP = 0;
        UnalignedSinglesArrays * unalignedSingleEndArrays = ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP );
        totalReadsProceedToDP += unalignedSingleEndArrays->array[0]->totalNum;
        printf ( "[Main] %u unaligned reads are proceeded to DP.\n", totalReadsProceedToDP );

        if ( totalReadsProceedToDP > 0 )
        {
            // ======================================================================================
            // | IF THE INDEX IS NOT IN DEVICE, THEN                                                |
            // | COPY INDEX TO DEVICE MEMORY                                                        |
            // ======================================================================================
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

            DPForUnalignSingle2 ( unalignedSingleEndArrays,
                                  queries, readLengths, readIDs, upkdQueryNames, upkdQualities,
                                  maxReadLength,
                                  index,
                                  _bwt, _revBwt,
                                  _occ, _revOcc,
                                  input_options.alignmentType, &dpParameters,
                                  numSingleDPAligned, numSingleDPAlignment,
                                  accumReadNum, input_options.outputFormat,
                                  NULL, samOutputUnpairFilePtr );
            indexLoadedToGPU = 0;
            // sort the alignment results
            sortReadPtrs ( ( AllHits * ) hspaux->allHits );
            // output the alignment results
            outputSingleResultForPairEnds ( ( AllHits * ) hspaux->allHits,
                                            queries, readLengths, readIDs, upkdQueryNames, upkdQualities,
                                            maxReadLength, accumReadNum, input_options.outputFormat,
                                            samOutputUnpairFilePtr, index );
        }

        // printf("Finished DP for single unaligned reads\n");
        // printf("[Main] Number of reads aligned by DP: %u (number of alignments: %u)\n", numDPAlignedSingle, numDPAlignment);
        printf ( "[Main] Number of reads aligned by single-end DP: %u\n", numSingleDPAligned );
        alignmentTime = getElapsedTime ( startTime );
        printf ( "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
        totalAlignmentTime += alignmentTime - lastEventTime;
        lastEventTime = alignmentTime;
        // numOfAlignedRead +=  numDPAlignedSingle;
        // numOfAnswer += numDPAlignment;
        // printf("[Main] Total Number of reads aligned: %u (number of alignments: %llu)\n", numOfAlignedRead, numOfAnswer);
        // printf("[Main] Total Number of reads aligned: %u\n", numOfAlignedRead);
#ifdef BGS_OUTPUT_DP_MESSAGE
        printf ( "*********************************************************\n" );
#endif
        printf ( "\n" );
        // reset the unaligned-singles array
        unalignedSingleEndArrays->array[0]->totalNum = 0;
        input_options.readType = PAIR_END_READ;
    }

#endif
    input_options.alignmentType = orig_align_type;
}

// Perform SOAP3-DP Single Alignment
void soap3_dp_single_align ( uint * queries, uint * readLengths, uint numMismatch, uint wordPerQuery,
                             ullint maxBatchSize, uint numQueries, uint accumReadNum,
                             Soap3Index * index,
                             uint * _bwt, uint * _revBwt,
                             uint * _occ, uint * _revOcc,
                             IniParams ini_params, InputOptions input_options,
                             uint maxReadLength, uint detected_read_length,
                             char * upkdQualities,
                             uint * unAlignedReads, uint & numOfUnPaired,
                             uint * readIDs, char * upkdQueryNames,
                             char ** currOutputFileName, samfile_t ** currSamOutputFilePtr,
                             char * outputDPFileName, samfile_t * samOutputDPFilePtr,
                             unsigned long long & numOfAnswer, uint & numOfAlignedRead,
                             ReadInputForDPArrays * readInputForDPall,
                             UnalignedSinglesArrays * unalignedSinglesArrays,
                             double startTime, double & lastEventTime, double & totalAlignmentTime,
                             uint & indexLoadedToGPU )
{
    double alignmentTime, copyTime;
    HSPAux * hspaux = index->sraIndex->hspaux;

    // ======================================================================================
    // | IF THE INDEX IS NOT IN DEVICE, THEN                                                |
    // | COPY INDEX TO DEVICE MEMORY                                                        |
    // ======================================================================================

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
    // | FOR LONG READS (with length > 120), only align the first 100 bases
    // ======================================================================================
    uint * seedLengths = NULL;

    if ( maxReadLength > LONG_READ_LEN )
    {
        seedLengths = ( uint * ) malloc ( numQueries * sizeof ( uint ) );
        uint i;

        for ( i = 0; i < numQueries; i++ )
        {
            if ( readLengths[i] > LONG_READ_LEN )
            {
                seedLengths[i] = SOAP3_SEED_LEN;
            }
            else
            {
                seedLengths[i] = readLengths[i];
            }
        }
    }

    //*******************************//
    // Perform Alignment             //
    //*******************************//

    if ( input_options.alignmentType == OUTPUT_ALL_BEST )
    {
        all_best_single_alignment ( queries, readLengths, seedLengths, input_options.numMismatch, wordPerQuery,
                                    maxBatchSize, numQueries, accumReadNum,
                                    index,
                                    _bwt, _revBwt, _occ, _revOcc,
                                    ini_params, input_options,
                                    upkdQualities,
                                    unAlignedReads, numOfUnPaired,
                                    readIDs, upkdQueryNames,
                                    currOutputFileName, currSamOutputFilePtr,
                                    numOfAnswer, numOfAlignedRead, readInputForDPall->inputArrays,
                                    unalignedSinglesArrays->array );
    }
    else if ( input_options.alignmentType == OUTPUT_UNIQUE_BEST ||
              input_options.alignmentType == OUTPUT_RANDOM_BEST )
    {
        best_single_alignment ( queries, readLengths, seedLengths, input_options.numMismatch, wordPerQuery,
                                maxBatchSize, numQueries, accumReadNum,
                                index,
                                _bwt, _revBwt, _occ, _revOcc,
                                ini_params, input_options,
                                upkdQualities,
                                unAlignedReads, numOfUnPaired,
                                readIDs, upkdQueryNames,
                                currOutputFileName, currSamOutputFilePtr,
                                numOfAnswer, numOfAlignedRead, readInputForDPall->inputArrays,
                                unalignedSinglesArrays->array );
    }
    else
    {
        all_valid_alignment ( queries, readLengths, seedLengths, input_options.numMismatch, wordPerQuery,
                              maxBatchSize, numQueries, accumReadNum,
                              index,
                              _bwt, _revBwt, _occ, _revOcc,
                              ini_params, input_options,
                              upkdQualities,
                              unAlignedReads, numOfUnPaired,
                              readIDs, upkdQueryNames,
                              currOutputFileName, currSamOutputFilePtr,
                              numOfAnswer, numOfAlignedRead, 1, readInputForDPall->inputArrays,
                              readInputForDPall->inputArrays, readInputForDPall->inputArrays,
                              unalignedSinglesArrays->array );
    }

    printf ( "[Main] Finished alignment with <= %i mismatches\n", input_options.numMismatch );
    // printf("[Main] Number of reads aligned: %u (number of alignments: %llu)\n", numOfAlignedRead, numOfAnswer);
    printf ( "[Main] Number of reads aligned: %u\n", numOfAlignedRead );
    alignmentTime = getElapsedTime ( startTime );
    printf ( "[Main] Elapsed time : %9.4f seconds\n\n", alignmentTime - lastEventTime );
    totalAlignmentTime += alignmentTime - lastEventTime;
    lastEventTime = alignmentTime;

    // ======================================================================================
    // | FOR LONG READS (with length > 120), release the seedLengths
    // ======================================================================================

    if ( maxReadLength > LONG_READ_LEN )
    {
        delete ( seedLengths );
    }

    /////////////////////////////////////////////////
    // FOR SINGLE READS TO PERFORM DP IF NECESSARY //
    /////////////////////////////////////////////////
    // For DP module, if SAM format and all-best are both selected, then
    // output format is needed to set to all-valid.
    int orig_align_type = input_options.alignmentType;

    if ( input_options.alignmentType == OUTPUT_ALL_BEST &&
            input_options.outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    {
        input_options.alignmentType = OUTPUT_ALL_VALID;
    }

    unsigned int numDPAlignedSingle = 0;
    unsigned int numDPAlignment = 0;
    // Parameters for DP
    DPParameters dpParameters;
    getParameterForAllDP ( dpParameters, input_options, ini_params );

    if ( input_options.enableDP == 1 )
    {
        // DP Parameters for single-end alignment
        getParameterForSingleDP ( dpParameters, input_options, ini_params, detected_read_length );
        hspaux->singleDPcutoffThreshold = dpParameters.paramRead[0].cutoffThreshold;
        // printDPParameters(dpParameters);
#ifdef BGS_OUTPUT_DP_MESSAGE
        printf ( "*********************************************************\n" );
        printf ( "DP FOR UNALIGNED SINGLE READS\n" );
#endif
        unsigned int totalReadsProceedToDP = 0;

        for ( int threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
        {
            totalReadsProceedToDP += unalignedSinglesArrays->array[threadId]->totalNum;
        }

        printf ( "[Main] %u unaligned reads are proceeded to DP.\n", totalReadsProceedToDP );

        // print out the read IDs
        // printAllReadIDs(unalignedSinglesArrays, accumReadNum);

        if ( totalReadsProceedToDP > 0 )
        {
            // ======================================================================================
            // | IF THE INDEX IS NOT IN DEVICE, THEN                                                |
            // | COPY INDEX TO DEVICE MEMORY                                                        |
            // ======================================================================================
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

            // prepare the file pointer for output
            FILE * outputDPFile = NULL;

            //Output file for SAM output is handled by SAM API
            //=> OutFilePtr is not used for SAM API output format.
            if ( outputDPFileName != NULL )
            {
                switch ( input_options.outputFormat )
                {
                    case SRA_OUTPUT_FORMAT_SAM_API:
                        break;

                    default:
                        outputDPFile = ( FILE * ) fopen ( outputDPFileName, "a" );
                }
            }

            // printAllReadIDs(unalignedSinglesArrays, accumReadNum, readIDs);
            DPForUnalignSingle2 ( unalignedSinglesArrays,
                                  queries, readLengths, readIDs, upkdQueryNames, upkdQualities,
                                  maxReadLength,
                                  index,
                                  _bwt, _revBwt,
                                  _occ, _revOcc,
                                  input_options.alignmentType, &dpParameters,
                                  numDPAlignedSingle, numDPAlignment,
                                  accumReadNum, input_options.outputFormat,
                                  outputDPFile, samOutputDPFilePtr );
            indexLoadedToGPU = 0;

            if ( outputDPFileName != NULL )
            {
                switch ( input_options.outputFormat )
                {
                    case SRA_OUTPUT_FORMAT_SAM_API:
                        break;

                    default:
                        fclose ( outputDPFile );
                }
            }
        }

        // printf("Finished DP for single unaligned reads\n");
        // printf("[Main] Number of reads aligned by DP: %u (number of alignments: %u)\n", numDPAlignedSingle, numDPAlignment);
        printf ( "[Main] Number of reads aligned by DP: %u\n", numDPAlignedSingle );
        alignmentTime = getElapsedTime ( startTime );
        printf ( "[Main] Elapsed time : %9.4f seconds\n", alignmentTime - lastEventTime );
        totalAlignmentTime += alignmentTime - lastEventTime;
        lastEventTime = alignmentTime;
        numOfAlignedRead +=  numDPAlignedSingle;
        numOfAnswer += numDPAlignment;
        //printf("[Main] Total Number of reads aligned: %u (number of alignments: %llu)\n", numOfAlignedRead, numOfAnswer);
        printf ( "[Main] Total Number of reads aligned: %u\n", numOfAlignedRead );
#ifdef BGS_OUTPUT_DP_MESSAGE
        printf ( "*********************************************************\n" );
#endif
        printf ( "\n" );
    }

    input_options.alignmentType = orig_align_type;
}
