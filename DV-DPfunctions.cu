/*
 *
 *    DV-DPfunctions.cu
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

#include "DV-DPfunctions.h"
#include "OutputDPResult.h"
#include <assert.h>

#include <algorithm>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// GPU function definition //////////////////////////
///////////////////////////////////////////////////////////////////////////////

__forceinline__ __device__ int _MAX ( int x, int y )
{
    //return (x > y ? x : y);
    return max ( x, y );
}
__forceinline__ __device__ int _MIN ( int x, int y )
{
    //return (x < y ? x : y);
    return min ( x, y );
}
__forceinline__ __device__ short _LOW_THRESHOLD ( int x )
{
    return ( short ) max ( x, -32000 );
}

texture <uint> texPatterns;
texture <uint> texSequences;
#define DP_SCORE_NEG_INFINITY -32000
#define dist(x,y) (x==y ? MatchScore : MismatchScore)
#define GapInit (GapOpenScore - GapExtendScore)
#define PatternLength() (maxReadLength + maxDPTableLength)

#define MC_ScoreAddr(X,j,i) *(X + ((j)*LPARA + (((i)>>1)<<6) + TPARA + ((i)&0x1)))
#define MC_DnaUnpack(X,i) ((X[dnaTPARA + (((i)>>4)<<5)] >> ((15-((i)&0xF))<<1)) & 3)
#define MC_ReadUnpack(X,i) ((tex1Dfetch(texPatterns, readTPARA + (((i)>>4)<<5)) >> ((15-((i)&0xF))<<1)) & 3)

__device__ void DPScoreNHitPos ( uint * packedDNASequence, uint DNALength, uint maxDNALength, uint maxDPTableLength,
                                 uint * readSequence, uint readLength, uint maxReadLength,
                                 int MatchScore, int MismatchScore,
                                 int GapOpenScore, int GapExtendScore,
                                 int clipLtCheckLoc, int clipRtCheckLoc,
                                 int anchorLeftLoc, int anchorRightLoc,
                                 int & maxScore, uint & hitPos,
                                 void * DPTable, uint threadId )
{
    short * score = ( short * ) DPTable + ( threadId / 32 ) * ( maxDPTableLength * maxReadLength * 32 * 2 );
    short * scoreOpen = ( short * ) score + maxDPTableLength * maxReadLength * 32;
    uint LPARA = maxReadLength << 5;
    uint TPARA = ( threadId & 0x1F ) << 1;
    uint readTPARA = ( threadId >> 5 ) * ( MC_CeilDivide16 ( maxReadLength ) << 5 ) + ( threadId & 0x1F );
    uint dnaTPARA = ( threadId >> 5 ) * ( MC_CeilDivide16 ( maxDNALength ) << 5 ) + ( threadId & 0x1F );
    int i, j;
    maxScore = DP_SCORE_NEG_INFINITY;
    //Initialize the first column
    MC_ScoreAddr ( score, 0, 0 ) = _LOW_THRESHOLD ( 0 );
    MC_ScoreAddr ( scoreOpen, 0, 0 ) = _LOW_THRESHOLD ( GapInit );
    int upScore = GapInit;

    for ( i = 1; i <= readLength; i++ )
    {
        if ( i <= clipLtCheckLoc )
        {
            MC_ScoreAddr ( score, 0, i ) = _LOW_THRESHOLD ( GapOpenScore );
            MC_ScoreAddr ( scoreOpen, 0, i ) = _LOW_THRESHOLD ( GapOpenScore + GapInit );
        }
        else
        {
            upScore += GapExtendScore;
            MC_ScoreAddr ( score, 0, i ) = _LOW_THRESHOLD ( upScore );
            MC_ScoreAddr ( scoreOpen, 0, i ) = _LOW_THRESHOLD ( upScore + GapInit );
        }
    }

    //prepare for filling table
    //start!
    int prevInitScore = 0;

    for ( j = 1; j <= DNALength; j++ )
    {
        uchar refChar = MC_DnaUnpack ( packedDNASequence, j );
        int initScore = ( ( j >= anchorLeftLoc ) ? DP_SCORE_NEG_INFINITY : 0 );
        int upScore = initScore;
        int scoreOpenUp = initScore + GapInit;
        int prevScoreUp = prevInitScore;
        int prevScoreR;
        int gappedScore;
#define Mc_innerDPSO_UpdateScoreTable() { \
        prevScoreR = MC_ScoreAddr(score,0,i); \
        gappedScore = _MAX(GapOpenScore + prevScoreR, GapExtendScore + MC_ScoreAddr(scoreOpen,0,i)); \
        MC_ScoreAddr(scoreOpen,0,i) = _LOW_THRESHOLD( gappedScore ); \
        scoreOpenUp = _MAX(GapExtendScore + scoreOpenUp, GapOpenScore + upScore); \
        gappedScore = _MAX(scoreOpenUp, gappedScore); \
        upScore = _MAX(gappedScore, prevScoreUp + dist(refChar, MC_ReadUnpack(readSequence, i))); \
        MC_ScoreAddr(score,0,i) = _LOW_THRESHOLD( upScore ); \
        prevScoreUp = prevScoreR; \
    }

        for ( i = 1; i <= readLength; i++ )
        {
            Mc_innerDPSO_UpdateScoreTable ();

            if ( i <= clipLtCheckLoc )
            {
                scoreOpenUp = _MAX ( initScore + GapInit, scoreOpenUp );
                prevScoreUp = _MAX ( prevInitScore, prevScoreUp );
            }

            if ( i >= clipRtCheckLoc &&
                    j >= anchorRightLoc &&
                    upScore > maxScore )
            {
                // update max score
                maxScore = upScore;
                hitPos = j;
            }
        }

        prevInitScore = initScore;
    }
}

__device__ void GenerateDPTable ( uint * packedDNASequence, uint DNALength, uint maxDNALength, uint maxDPTableLength,
                                  uint * readSequence, uint readLength, uint maxReadLength,
                                  int MatchScore, int MismatchScore,
                                  int GapOpenScore, int GapExtendScore,
                                  int clipLtCheckLoc, int clipRtCheckLoc,
                                  int anchorLeftLoc, int anchorRightLoc,
                                  uint refOffset, int & maxScore, uint & hitPos, uint & scRight,
                                  uint & maxScoreCount,
                                  void * DPTable, uint threadId )
{
    short * score = ( short * ) DPTable + ( threadId / 32 ) * ( maxDPTableLength * maxReadLength * 32 * 2 );
    short * scoreOpen = ( short * ) score + maxDPTableLength * maxReadLength * 32;
    uint LPARA = maxReadLength << 5;
    uint TPARA = ( threadId & 0x1F ) << 1;
    uint readTPARA = ( threadId >> 5 ) * ( MC_CeilDivide16 ( maxReadLength ) << 5 ) + ( threadId & 0x1F );
    uint dnaTPARA = ( threadId >> 5 ) * ( MC_CeilDivide16 ( maxDNALength ) << 5 ) + ( threadId & 0x1F );
    int i, j;
    maxScore = DP_SCORE_NEG_INFINITY;
    maxScoreCount = 0;
    //Initialize the first column
    MC_ScoreAddr ( score, 0, 0 ) = _LOW_THRESHOLD ( 0 );
    MC_ScoreAddr ( scoreOpen, 0, 0 ) = _LOW_THRESHOLD ( GapInit );
    int upScore = GapInit;

    for ( i = 1; i <= readLength; i++ )
    {
        if ( i <= clipLtCheckLoc )
        {
            MC_ScoreAddr ( score, 0, i ) = _LOW_THRESHOLD ( GapOpenScore );
            MC_ScoreAddr ( scoreOpen, 0, i ) = _LOW_THRESHOLD ( GapOpenScore + GapInit );
        }
        else
        {
            upScore += GapExtendScore;
            MC_ScoreAddr ( score, 0, i ) = _LOW_THRESHOLD ( upScore );
            MC_ScoreAddr ( scoreOpen, 0, i ) = _LOW_THRESHOLD ( upScore + GapInit );
        }
    }

    //prepare for filling table
    //start!
    int prevInitScore = 0;

    for ( j = 1; j <= DNALength; j++ )
    {
        uchar refChar = MC_DnaUnpack ( packedDNASequence, j + refOffset );
        int initScore = ( ( j + refOffset >= anchorLeftLoc ) ? DP_SCORE_NEG_INFINITY : 0 );
        int upScore = initScore;
        int scoreOpenUp = initScore + GapInit;
        int prevScoreUp = prevInitScore;
        int prevScoreR;
        int gappedScore;
        MC_ScoreAddr ( score, j, 0 ) = _LOW_THRESHOLD ( upScore );
        MC_ScoreAddr ( scoreOpen, j, 0 ) = _LOW_THRESHOLD ( scoreOpenUp );
#define Mc_UpdateScoreTable() { \
        prevScoreR = MC_ScoreAddr(score,j-1,i); \
        gappedScore = _MAX(GapOpenScore + prevScoreR, GapExtendScore + MC_ScoreAddr(scoreOpen,j-1,i)); \
        MC_ScoreAddr(scoreOpen,j,i) = _LOW_THRESHOLD( gappedScore ); \
        scoreOpenUp = _MAX(GapExtendScore + scoreOpenUp, GapOpenScore + upScore); \
        gappedScore = _MAX(scoreOpenUp, gappedScore); \
        upScore = _MAX(gappedScore, prevScoreUp + dist(refChar, MC_ReadUnpack(readSequence, i))); \
        MC_ScoreAddr(score,j,i) = _LOW_THRESHOLD( upScore ); \
        prevScoreUp = prevScoreR; \
    }

        for ( i = 1; i <= readLength; i++ )
        {
            Mc_UpdateScoreTable ();

            if ( i <= clipLtCheckLoc )
            {
                scoreOpenUp = _MAX ( initScore + GapInit, scoreOpenUp );
                prevScoreUp = _MAX ( prevInitScore, prevScoreUp );
            }

            if ( i >= clipRtCheckLoc &&
                    j + refOffset >= anchorRightLoc )
            {
                // update max score
                if ( upScore > maxScore )
                {
                    maxScore = upScore;
                    hitPos = j;
                    scRight = readLength - i;
                    maxScoreCount = 1;
                }
                else if ( upScore == maxScore )
                {
                    ++maxScoreCount;
                }
            }
        }

        prevInitScore = initScore;
    }
}

__global__ void SemiGlobalAligntment ( uint * packedDNASequence, uint * DNALengths, uint maxDNALength, uint maxDPTableLength,
                                       uint * readSequence, uint * readLengths, uint maxReadLength,
                                       int * maxScores, uint * hitLocs, uint * startOffsets,
                                       uint * clipLtSizes, uint * clipRtSizes,
                                       uint * anchorLeftLocs, uint * anchorRightLocs, uint numOfThreads,
                                       int MatchScore, int MismatchScore,
                                       int GapOpenScore, int GapExtendScore,
                                       void * DPTable, uint * maxScoreCounts,
                                       int alignmentScheme = 1 )
{
    // precondition: MAX_READ_LENGTH should be multiple of 4
    uint threadId = blockIdx.x * DP_THREADS_PER_BLOCK + threadIdx.x;

    if ( threadId < numOfThreads )
    {
        uint clipLtSize = ( clipLtSizes == NULL ) ? 0 : clipLtSizes[threadId];
        uint clipRtSize = ( clipRtSizes == NULL ) ? 0 : clipRtSizes[threadId];
        uint anchorLeftLoc = ( anchorLeftLocs == NULL ) ?
                             maxDNALength : anchorLeftLocs[threadId];
        uint anchorRightLoc = ( anchorRightLocs == NULL ) ?
                              0 : anchorRightLocs[threadId];
        uint readLength = readLengths[threadId];
        uint DNALength = DNALengths[threadId];
        uint refOffset = 0, hitPos = 0, scRight = 0;
        uint maxScoreCount = 0;
        int maxScore;

        if ( alignmentScheme == 2 )
        {
            // get maxScore & hitPos
            DPScoreNHitPos ( packedDNASequence, DNALength, maxDNALength, maxDPTableLength,
                             readSequence, readLength, maxReadLength,
                             MatchScore, MismatchScore,
                             GapOpenScore, GapExtendScore,
                             clipLtSize, readLength - clipRtSize,
                             anchorLeftLoc, anchorRightLoc,
                             maxScore, hitPos,
                             DPTable, threadId );
            // decide offset
            DNALength = _MIN ( DNALength, maxDPTableLength - 1 );

            if ( hitPos > DNALength )
            {
                refOffset = hitPos - DNALength;

                if ( refOffset >= anchorLeftLoc )
                {
                    refOffset = ( anchorLeftLoc > 0 ? anchorLeftLoc - 1 : 0 );
                }
            }
        }

        GenerateDPTable ( packedDNASequence, DNALength, maxDNALength, maxDPTableLength,
                          readSequence, readLength, maxReadLength,
                          MatchScore, MismatchScore,
                          GapOpenScore, GapExtendScore,
                          clipLtSize, readLength - clipRtSize,
                          anchorLeftLoc, anchorRightLoc,
                          refOffset, maxScore, hitPos, scRight,
                          maxScoreCount,
                          DPTable, threadId );
        startOffsets[threadId] = refOffset;
        hitLocs[threadId] = hitPos;
        hitLocs[threadId] = hitLocs[threadId];

        if ( clipRtSizes != NULL )
        { clipRtSizes[threadId] = scRight; }

        maxScores[threadId] = maxScore;
        maxScoreCounts[threadId] = maxScoreCount;
    }
}

__global__ void GPUBacktrack ( uint * packedDNASequence, uint * DNALengths, uint maxDNALength, uint maxDPTableLength,
                               uint * readSequence, uint * readLengths, uint maxReadLength,
                               int * maxScores, uint * hitLocs, uint * startOffsets,
                               uint * clipLtSizes, uint * clipRtSizes, uint * anchorLeftLocs, uint numOfThreads,
                               int MatchScore, int MismatchScore,
                               int GapOpenScore, int GapExtendScore, int * cutoffThresholds,
                               void * DPTable, uchar * pattern )
{
    uint threadId = blockIdx.x * DP_THREADS_PER_BLOCK + threadIdx.x;

    if ( threadId < numOfThreads )
    {
        if ( maxScores[threadId] >= cutoffThresholds[threadId] )
        {
            short * score = ( short * ) DPTable + ( threadId / 32 ) * ( maxDPTableLength * maxReadLength * 32 * 2 );
            short * scoreOpen = ( short * ) score + maxDPTableLength * maxReadLength * 32;
            uint LPARA = maxReadLength << 5;
            uint TPARA = ( threadId & 0x1F ) << 1;
            uint readTPARA = ( threadId >> 5 ) * ( MC_CeilDivide16 ( maxReadLength ) << 5 ) + ( threadId & 0x1F );
            uint dnaTPARA = ( threadId >> 5 ) * ( MC_CeilDivide16 ( maxDNALength ) << 5 ) + ( threadId & 0x1F );
            uchar * curPattern = pattern + threadId * PatternLength ();
            uint pIndex = 0;
            uint readLength = readLengths[threadId];
            uint refOffset = startOffsets[threadId];
            uint anchorLeftLoc = ( anchorLeftLocs == NULL ) ?
                                 maxDNALength : anchorLeftLocs[threadId];
#define MC_PatternAppend(x) { curPattern[pIndex] = (x); ++pIndex; }
            uint clipLtCheckLoc = ( ( clipLtSizes == NULL ) ? 0 : clipLtSizes[threadId] );
            uint clipRtLength = ( clipRtSizes == NULL ) ? 0 : clipRtSizes[threadId];

            if ( clipRtLength > 0 )
            {
                MC_PatternAppend ( 'S' );
                MC_PatternAppend ( 'V' );
                MC_PatternAppend ( clipRtLength );
            }

            uint readPos = readLength - clipRtLength;
            uint refIndex = hitLocs[threadId];
            uchar readChar = MC_ReadUnpack ( readSequence, readPos );
            uchar refChar = MC_DnaUnpack ( packedDNASequence, refOffset + refIndex );
#define MC_NextRefCharNInitScore() { \
        --refIndex; \
        refChar = MC_DnaUnpack(packedDNASequence, refOffset+refIndex); \
        initScore = prevInitScore; \
        prevInitScore = ((refOffset+refIndex > anchorLeftLoc) ? \
                         DP_SCORE_NEG_INFINITY : 0); \
    }
#define MC_NextReadChar() { --readPos; readChar = MC_ReadUnpack(readSequence, readPos); }
            short curScore = MC_ScoreAddr ( score, refIndex, readPos );
            short nextScore;
            short initScore = ( ( refOffset + refIndex >= anchorLeftLoc ) ?
                                DP_SCORE_NEG_INFINITY : 0 );
            short prevInitScore = ( ( refOffset + refIndex > anchorLeftLoc ) ?
                                    DP_SCORE_NEG_INFINITY : 0 );
            enum DP_BacktrackState { DP_BT_NORMAL, DP_BT_I_EXT, DP_BT_D_EXT, \
                                     DP_BT_SM_EXIT, DP_BT_SI_EXIT
                                   };
            DP_BacktrackState state = DP_BT_NORMAL;

            while ( readPos > 0 && refIndex > 0 )
            {
                // check match
                if ( state == DP_BT_NORMAL )
                {
                    if ( curScore == dist ( refChar, readChar ) +
                            ( nextScore = MC_ScoreAddr ( score, refIndex - 1, readPos - 1 ) ) )
                    {
                        // Match/Mismatch --> reference X : read X
                        MC_PatternAppend ( 'm' - ( ( refChar == readChar ) << 5 ) );
                        MC_NextRefCharNInitScore ();
                        MC_NextReadChar ();
                        curScore = nextScore;
                    }
                    else if ( curScore == GapOpenScore +
                              ( nextScore = MC_ScoreAddr ( score, refIndex - 1, readPos ) ) )
                    {
                        // Deletion --> reference X : read -
                        MC_PatternAppend ( 'D' );
                        MC_NextRefCharNInitScore ();
                        curScore = nextScore;
                    }
                    else if ( curScore == GapExtendScore +
                              MC_ScoreAddr ( scoreOpen, refIndex - 1, readPos ) )
                    {
                        // (start extension) Deletion --> reference X : read -
                        MC_PatternAppend ( 'D' );
                        MC_NextRefCharNInitScore ();
                        curScore -= GapExtendScore;
                        state = DP_BT_D_EXT;
                    }
                    else
                    {
                        // check for left soft clip
                        if ( readPos <= clipLtCheckLoc + 1 )
                        {
                            if ( curScore == prevInitScore + dist ( refChar, readChar ) )
                            {
                                state = DP_BT_SM_EXIT;
                                break;
                            }
                            else if ( curScore == initScore + GapOpenScore )
                            {
                                state = DP_BT_SI_EXIT;
                                break;
                            }
                        }

                        if ( curScore == GapOpenScore +
                                ( nextScore = MC_ScoreAddr ( score, refIndex, readPos - 1 ) ) )
                        {
                            //  Insertion --> reference - : read X
                            MC_PatternAppend ( 'I' );
                            MC_NextReadChar ();
                            curScore = nextScore;
                        }
                        else
                        {
                            // (start extension) Insertion --> reference - : read X
                            MC_PatternAppend ( 'I' );
                            MC_NextReadChar ();
                            curScore -= GapExtendScore;
                            state = DP_BT_I_EXT;
                        }
                    }
                }
                else
                {
                    // extension state
                    if ( state == DP_BT_D_EXT )
                    {
                        // (extension) Deletion --> reference X : read -
                        MC_PatternAppend ( 'D' );
                        MC_NextRefCharNInitScore ();
                    }
                    else
                    {
                        // (extension) check for left soft clip
                        if ( readPos <= clipLtCheckLoc + 1 &&
                                curScore == initScore + GapOpenScore )
                        {
                            state = DP_BT_SI_EXIT;
                            break;
                        }

                        // (extension) Insertion --> reference - : read X
                        MC_PatternAppend ( 'I' );
                        MC_NextReadChar ();
                    }

                    if ( curScore == GapOpenScore + ( nextScore = MC_ScoreAddr ( score, refIndex, readPos ) ) )
                    {
                        state = DP_BT_NORMAL;
                        curScore = nextScore;
                    }
                    else
                    { curScore -= GapExtendScore; }
                }
            }

            //last proc
            if ( refIndex == 0 )
            {
                uint scNum = min ( clipLtCheckLoc, readPos );

                if ( scNum < readPos )
                {
                    MC_PatternAppend ( 'I' );
                    MC_PatternAppend ( 'V' );
                    MC_PatternAppend ( readPos - scNum );
                }

                MC_PatternAppend ( 'S' );
                MC_PatternAppend ( 'V' );
                MC_PatternAppend ( scNum );
            }
            else if ( state == DP_BT_SI_EXIT )
            {
                MC_PatternAppend ( 'I' );
                MC_PatternAppend ( 'S' );
                MC_PatternAppend ( 'V' );
                MC_PatternAppend ( readPos - 1 );
            }
            else if ( state == DP_BT_SM_EXIT )
            {
                MC_PatternAppend ( 'm' - ( ( refChar == readChar ) << 5 ) );
                MC_PatternAppend ( 'S' );
                MC_PatternAppend ( 'V' );
                MC_PatternAppend ( readPos - 1 );
                refIndex -= 1;
            }

            MC_PatternAppend ( 0 );
            hitLocs[threadId] = refOffset + refIndex;
        }
    }
}



///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// SemiGlobalAligner  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////

uint SemiGlobalAligner::estimateThreadSize ( int maxReadLength, int maxDNALength )
{
    uint tableSizeOfThread = 2 * maxDNALength * maxReadLength * sizeof ( short );
    uint otherSizeOfThread = MC_CeilDivide16 ( maxDNALength ) * sizeof ( uint ) + //_packedDNASequence
                             MC_CeilDivide16 ( maxReadLength ) * sizeof ( uint ) + //_packedReadSequence
                             sizeof ( uint ) + //_DNALengths
                             sizeof ( uint ) + //_readLengths
                             sizeof ( uint ) + //_startLocs
                             sizeof ( int ) * 2 + //_scores, _cutoffThresholds
                             sizeof ( uint ) + //_startOffsets
                             sizeof ( uint ) * 2 + //_clipLtSizes, Rt
                             sizeof ( uint ) * 2 + //_anchorLeftLocs, Right
                             sizeof ( uint ) + //_hitLocs
                             sizeof ( uint ) + //_maxScoreCounts
                             ( maxReadLength + maxDNALength ) * sizeof ( uchar ); //_pattern
    return tableSizeOfThread + otherSizeOfThread + 16; // 16 is padding size
}

SemiGlobalAligner::SemiGlobalAligner ()
{
    n_conf = 6;

    for ( int i = 0; i < 4; i++ )
    {
        blockConf[i] = 64 - 16 * i;
        coefConf[i] = 3;
    }

    blockConf[4] = 8;
    coefConf[4] = 2;
    blockConf[5] = 2;
    coefConf[5] = 1.5; // <- desperate
}

int SemiGlobalAligner::tryAlloc ( size_t estimatedThreadSize, size_t numOfBlocks )
{
    void * _testMalloc;
    cudaError_t err = cudaMalloc ( ( void ** ) &_testMalloc, estimatedThreadSize * numOfBlocks * DP_THREADS_PER_BLOCK );

    if ( err == cudaSuccess )
    {
        cudaFree ( _testMalloc );
        return 0;
    }

    return -1;
}

void SemiGlobalAligner::decideConfiguration (
    int maxReadLength, int maxDNALength,
    int & maxDPTableLength, int & numOfBlocks,
    int & patternLength, DPParameters & dpPara
)
{
#define DP_Effective_Region(l)  (l + l/2 + 8)
    size_t avail, total;
    cudaMemGetInfo ( &avail, &total );
    size_t availableMemory = avail * 88 / 100; // 88%
    size_t estimatedThreadSize      = estimateThreadSize ( maxReadLength, maxDNALength );
    size_t availableBlocks          = availableMemory / ( estimatedThreadSize * DP_THREADS_PER_BLOCK );
    size_t estimatedThreadSize_2    = estimateThreadSize ( maxReadLength, DP_Effective_Region ( maxReadLength ) );
    size_t availableBlocks_2        = availableMemory / ( estimatedThreadSize_2 * DP_THREADS_PER_BLOCK );
    int successFlag = 0;

    for ( int i = 0; i < n_conf; i++ )
    {
        if ( availableBlocks >= blockConf[i] )
        {
            //              printf("[%d] scheme 1, try allocate\n", i);
            if ( 0 == tryAlloc ( estimatedThreadSize, blockConf[i] ) )
            {
                successFlag = 1;
                alignmentScheme = 1;
                numOfBlocks = blockConf[i];
                maxDPTableLength = maxDNALength;
                break;
            }
        }

        if ( availableBlocks_2 >= blockConf[i] &&
                maxDNALength >= ( int ) ( maxReadLength * coefConf[i] ) )
        {
            //              printf("[%d] scheme 2, try allocate\n", i);
            if ( 0 == tryAlloc ( estimatedThreadSize_2, blockConf[i] ) )
            {
                successFlag = 1;
                alignmentScheme = 2;
                numOfBlocks = blockConf[i];
                maxDPTableLength = DP_Effective_Region ( maxReadLength );
                break;
            }
        }
    }

    if ( !successFlag )
    {
        // configuration failed
        printf ( "[DPfunc] error: insufficient GPU memory, cannot perform DP\n" );
        exit ( -1 );
    }

    // check invalid configuration
    if ( numOfBlocks < 32 )
    {
        printf ( "[DPfunc] warning: insufficient GPU memory, performance might degrade\n" );
    }

    if ( dpPara.matchScore > 30 )
    {
        printf ( "[DPfunc] warning: MatchScore (set to %d) should not exceed 30\n", dpPara.matchScore );
    }

    fflush ( stdout );
    patternLength = PatternLength ();
}

void SemiGlobalAligner::init (
    int batchSize,
    int maxReadLength, int maxDNALength, int maxDPTableLength,
    DPParameters & dpPara
)
{
    MC_MemberCopy5 ( this->, , batchSize, maxReadLength, maxDNALength, maxDPTableLength, dpPara );
    cudaFuncSetCacheConfig ( SemiGlobalAligntment, cudaFuncCachePreferL1 );
    //      cudaFuncSetCacheConfig( GPUBacktrack, cudaFuncCachePreferL1 );
    //      showGPUMemInfo("before");
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_DPTable, ( size_t ) 2 * maxDPTableLength * maxReadLength * batchSize * sizeof ( short ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_packedDNASequence, batchSize * MC_CeilDivide16 ( maxDNALength ) * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_packedReadSequence, batchSize * MC_CeilDivide16 ( maxReadLength ) * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_DNALengths, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_readLengths, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_startLocs, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_startOffsets, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_hitLocs, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_scores, batchSize * sizeof ( int ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_cutoffThresholds, batchSize * sizeof ( int ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_clipLtSizes, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_clipRtSizes, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_anchorLeftLocs, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_anchorRightLocs, batchSize * sizeof ( uint ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_pattern, batchSize * PatternLength () * sizeof ( uchar ) ) );
    DP_HANDLE_ERROR ( cudaMalloc ( ( void ** ) &_maxScoreCounts, batchSize * sizeof ( uint ) ) );
    //      showGPUMemInfo("alloc");
    cudaBindTexture(NULL, texPatterns, _packedReadSequence,
            batchSize * MC_CeilDivide16 ( maxReadLength ) * sizeof ( uint ));
    cudaBindTexture(NULL, texSequences, _packedDNASequence,
            batchSize * MC_CeilDivide16 ( maxDNALength ) * sizeof ( uint ));
}

void SemiGlobalAligner::performAlignment (
    uint * packedDNASequence, uint * DNALengths,
    uint * packedReadSequence, uint * readLengths,
    int * cutoffThresholds, int * scores, uint * hitLocs,
    uint * maxScoreCounts,
    uchar * pattern, int numOfThreads,
    uint * clipLtSizes, uint * clipRtSizes,
    uint * anchorLeftLocs, uint * anchorRightLocs )
{
    DP_HANDLE_ERROR ( cudaMemcpy ( _packedDNASequence, packedDNASequence, batchSize * MC_CeilDivide16 ( maxDNALength ) * sizeof ( uint ), cudaMemcpyHostToDevice ) );
    DP_HANDLE_ERROR ( cudaMemcpy ( _packedReadSequence, packedReadSequence, batchSize * MC_CeilDivide16 ( maxReadLength ) * sizeof ( uint ), cudaMemcpyHostToDevice ) );
    DP_HANDLE_ERROR ( cudaMemcpy ( _DNALengths, DNALengths, batchSize * sizeof ( uint ), cudaMemcpyHostToDevice ) );
    DP_HANDLE_ERROR ( cudaMemcpy ( _readLengths, readLengths, batchSize * sizeof ( uint ), cudaMemcpyHostToDevice ) );
    DP_HANDLE_ERROR ( cudaMemcpy ( _cutoffThresholds, cutoffThresholds, batchSize * sizeof ( int ), cudaMemcpyHostToDevice ) );
#define MC_CheckCopy_TypeUINT(CPU_para, GPU_para, tmp_para) { \
        if (CPU_para != NULL) { \
            DP_HANDLE_ERROR( cudaMemcpy(GPU_para, CPU_para, batchSize * sizeof(uint), cudaMemcpyHostToDevice) ); \
            tmp_para = GPU_para; \
        } \
    }
    uint * _clipLt = NULL, *_clipRt = NULL;
    MC_CheckCopy_TypeUINT ( clipLtSizes, _clipLtSizes, _clipLt );
    MC_CheckCopy_TypeUINT ( clipRtSizes, _clipRtSizes, _clipRt );
    uint * _anchorLtLocs = NULL, *_anchorRtLocs = NULL;
    MC_CheckCopy_TypeUINT ( anchorLeftLocs, _anchorLeftLocs, _anchorLtLocs );
    MC_CheckCopy_TypeUINT ( anchorRightLocs, _anchorRightLocs, _anchorRtLocs );
    int blocksNeeded = ( numOfThreads + DP_THREADS_PER_BLOCK - 1 ) / DP_THREADS_PER_BLOCK;
    SemiGlobalAligntment <<< blocksNeeded, DP_THREADS_PER_BLOCK>>> (
        _packedDNASequence, _DNALengths, maxDNALength, maxDPTableLength,
        _packedReadSequence, _readLengths, maxReadLength,
        _scores, _hitLocs, _startOffsets,
        _clipLt, _clipRt, _anchorLtLocs, _anchorRtLocs, numOfThreads,
        dpPara.matchScore, dpPara.mismatchScore,
        dpPara.openGapScore, dpPara.extendGapScore,
        _DPTable, _maxScoreCounts,
        alignmentScheme
    );
    GPUBacktrack <<< blocksNeeded, DP_THREADS_PER_BLOCK>>> (
        _packedDNASequence, _DNALengths, maxDNALength, maxDPTableLength,
        _packedReadSequence, _readLengths, maxReadLength,
        _scores, _hitLocs, _startOffsets,
        _clipLt, _clipRt, _anchorLtLocs, numOfThreads,
        dpPara.matchScore, dpPara.mismatchScore,
        dpPara.openGapScore, dpPara.extendGapScore, _cutoffThresholds,
        _DPTable, _pattern
    );
    //Fetch results
    DP_HANDLE_ERROR ( cudaMemcpy ( scores, _scores, batchSize * sizeof ( int ), cudaMemcpyDeviceToHost ) );
    DP_HANDLE_ERROR ( cudaMemcpy ( hitLocs, _hitLocs, batchSize * sizeof ( uint ), cudaMemcpyDeviceToHost ) );
    DP_HANDLE_ERROR ( cudaMemcpy ( pattern, _pattern, batchSize * PatternLength () * sizeof ( uchar ), cudaMemcpyDeviceToHost ) );
    DP_HANDLE_ERROR ( cudaMemcpy ( maxScoreCounts, _maxScoreCounts, batchSize * sizeof ( uint ), cudaMemcpyDeviceToHost ) );

}

void SemiGlobalAligner::freeMemory ()
{
    cudaFree ( _DPTable );
    cudaFree ( _packedDNASequence );
    cudaFree ( _packedReadSequence );
    cudaFree ( _DNALengths );
    cudaFree ( _readLengths );
    cudaFree ( _hitLocs );
    cudaFree ( _startLocs );
    cudaFree ( _startOffsets );
    cudaFree ( _scores );
    cudaFree ( _cutoffThresholds );
    cudaFree ( _clipLtSizes );
    cudaFree ( _clipRtSizes );
    cudaFree ( _anchorLeftLocs );
    cudaFree ( _anchorRightLocs );
    cudaFree ( _pattern );
    cudaFree ( _maxScoreCounts );
}



//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// For output ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

AlgnmtFlags::AlgnmtFlags ( uint range )
{
    size = ( range + 31 ) / 32;
    MC_CheckMalloc ( flags,   uint,   size );

    for ( int i = 0; i < 32; i++ )
    {
        MASK[i] = 1 << i;
    }

    pthread_mutex_init ( &occupy_mutex, NULL );
    clear ();
}

void AlgnmtFlags::clear ()
{
    memset ( flags, 0, size * sizeof ( uint ) );
}

void AlgnmtFlags::increaseSize ( uint newSize )
{
    uint * oldFlags = flags;
    uint oldSize = size;
    size  = newSize;
    MC_CheckMalloc ( flags,   uint,   size );
    memcpy ( flags, oldFlags, oldSize * sizeof ( uint ) );
    memset ( flags + oldSize, 0, ( newSize - oldSize ) * sizeof ( uint ) );
    free ( oldFlags );
}

void AlgnmtFlags::set ( int readID )
{
    pthread_mutex_lock ( &occupy_mutex );
    uint offset = readID >> 5;

    if ( offset >= size )
    {
        uint newSize = size * 2;

        while ( offset >= newSize )
        { newSize *= 2; }

        increaseSize ( newSize );
    }

    flags[offset] |= MASK[readID & 0x1F];
    pthread_mutex_unlock ( &occupy_mutex );
}

#define AlgnmtFlags_Get(in_flag, int32Offset, diff) { \
        uint flag = in_flag; \
        if (flag != 0) { \
            int offset = int32Offset << 5; \
            for (int j = 0; flag != 0; j++) { \
                if (flag & 1) { \
                    diff->push_back(offset + j); \
                } \
                flag >>= 1; \
            } \
        } \
    }

void AlgnmtFlags::get ( vector<int> * diff )
{
    for ( int i = 0; i < size; i++ )
    {
        AlgnmtFlags_Get ( flags[i], i, diff );
    }
}

inline void AlgnmtFlags::reserveSize ( AlgnmtFlags * algnFlags )
{
    if ( size < algnFlags->size )
    {
        this->increaseSize ( algnFlags->size );
    }
    else if ( size > algnFlags->size )
    {
        algnFlags->increaseSize ( size );
    }
}

void AlgnmtFlags::getXOR ( AlgnmtFlags * algnFlags, vector<int> * diff )
{
    reserveSize ( algnFlags );

    for ( int i = 0; i < size; i++ )
    {
        AlgnmtFlags_Get ( flags[i] ^ algnFlags->flags[i], i, diff );
    }
}

void AlgnmtFlags::AND ( AlgnmtFlags * algnFlags )
{
    reserveSize ( algnFlags );

    for ( int i = 0; i < size; i++ )
    {
        flags[i] &= algnFlags->flags[i];
    }
}

void AlgnmtFlags::XOR ( AlgnmtFlags * algnFlags )
{
    reserveSize ( algnFlags );

    for ( int i = 0; i < size; i++ )
    {
        flags[i] ^= algnFlags->flags[i];
    }
}

AlgnmtFlags::~AlgnmtFlags ()
{
    free ( flags );
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Alignment modules //////////////////////////////////////
/////////////////// The following code better be placed in a seperate file ///////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// standard space ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////


template <>
int isValid ( AlgnmtDPResult & a )
{
    return ( a.whichFromDP < 2 );
}

template <>
int ScoreCompare ( AlgnmtDPResult & a, AlgnmtDPResult & b )
{
#define MC_DPScoreCompare_SetValue(result, aligned, mismatch, score) { \
        if (result.whichFromDP == 0) { \
            aligned = 1; \
            mismatch = result.score_2; \
            score = result.score_1; \
        } \
        else \
            if (result.whichFromDP == 1) { \
                aligned = 1; \
                mismatch = result.score_1; \
                score = result.score_2; \
            } \
            else { \
                aligned = 0; \
                if (result.algnmt_1 != 0xFFFFFFFF) \
                    mismatch = result.score_1; \
                else \
                    mismatch = result.score_2; \
                score = 0; \
            } \
    }
    uint aligned_a, mismatch_a, score_a;
    uint aligned_b, mismatch_b, score_b;
    MC_DPScoreCompare_SetValue ( a, aligned_a, mismatch_a, score_a );
    MC_DPScoreCompare_SetValue ( b, aligned_b, mismatch_b, score_b );
    uint64 value_a = ( ( uint64 ) aligned_a << 63 ) | ( ( uint64 ) ( 0x1FFFFFFF - mismatch_a ) << 32 ) | ( score_a + 0x1FFFFFFF );
    uint64 value_b = ( ( uint64 ) aligned_b << 63 ) | ( ( uint64 ) ( 0x1FFFFFFF - mismatch_b ) << 32 ) | ( score_b + 0x1FFFFFFF );

    if ( value_a > value_b )
    { return 1; }
    else if ( value_a < value_b )
    { return -1; }
    else
    { return 0; }
}

template <>
int ScoreCompare ( SingleAlgnmtResult & a, SingleAlgnmtResult & b )
{
    if ( a.score > b.score )
    { return 1; }
    else if ( a.score < b.score )
    { return -1; }
    else
    { return 0; }
}

template <>
int ScoreCompare ( DeepDPAlignResult & a, DeepDPAlignResult & b )
{
    int score_a = a.score_1 + a.score_2;
    int score_b = b.score_1 + b.score_2;

    if ( score_a > score_b )
    { return 1; }
    else if ( score_a < score_b )
    { return -1; }
    else
    { return 0; }
}

template <>
bool ResultCompare ( const AlgnmtDPResult & a, const AlgnmtDPResult & b )
{
    return ( a.algnmt_1 < b.algnmt_1 );
}
template <>
bool ResultCompare ( const SingleAlgnmtResult & a, const SingleAlgnmtResult & b )
{
    return ( a.algnmt < b.algnmt );
}
template <>
bool ResultCompare ( const DeepDPAlignResult & a, const DeepDPAlignResult & b )
{
    return ( a.algnmt_1 < b.algnmt_1 );
}

QueryIDStream::QueryIDStream ()
{
    data = new vector<int>;
}
QueryIDStream::QueryIDStream ( BothUnalignedPairsArrays * input )
{
    data = new vector<int>;

    for ( int arrIndex = 0; arrIndex < input->arrayNum; arrIndex++ )
    {
        BothUnalignedPairs * array = input->array[arrIndex];

        for ( int pairIndex = 0; pairIndex < array->totalNum; pairIndex++ )
        {
            data->push_back ( array->readIDs[pairIndex] );
        }
    }
}
QueryIDStream::~QueryIDStream ()
{
    delete data;
}
void QueryIDStream::append ( QueryIDStream * stream )
{
    data->insert ( data->end (), stream->data->begin (), stream->data->end () );
}
void QueryIDStream::setBuffer ( vector<int> * input )
{
    delete data;
    data = input;
}



//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// single-dp space //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

using namespace SingleDP_Space;
#define DPS_SEEDING_BATCH_SIZE 256 * 1024
#define DPS_MARGIN(l) ((l>100) ? (l>>2) : 25)

CandidateStream::CandidateStream ()
{
    pthread_mutex_init ( &occupy_mutex, NULL );
}
void CandidateStream::append ( vector<CandidateInfo> * canInfo, AlgnmtFlags * alignFlags )
{
    pthread_mutex_lock ( &occupy_mutex );

    for ( vector<CandidateInfo>::iterator it = canInfo->begin ();
            it < canInfo->end (); ++it )
    {
        data.push_back ( *it );
        alignFlags->set ( it->readID );
    }

    pthread_mutex_unlock ( &occupy_mutex );
}

// ****
SingleEndSeedingEngine::SingleEndSeedingBatch::SingleEndSeedingBatch (
    uint batchSize, DPParameters * dpPara,
    uint * queries, uint * queryLengths, uint inputMaxReadLength
)
{
    MC_MemberCopy3 ( this->, , batchSize, queries, queryLengths );
    this->numOfCPUForSeeding    = dpPara->numOfCPUForSeeding;
    this->maxHitNum             = dpPara->paramRead[0].maxHitNum;
    this->maxSeedLength     = inputMaxReadLength;
    this->wordPerSeed       = MC_CeilDivide16 ( maxSeedLength );
    this->wordPerQuery      = getWordPerQuery ( inputMaxReadLength );
    this->algnResultArray = constructSingleAlgnResultArray ( numOfCPUForSeeding );
    MC_CheckMalloc ( readIDs, uint,   batchSize );
    MC_CheckMalloc ( lengths, uint,   batchSize );
    MC_CheckMalloc ( offsets, uint,   batchSize );
    MC_CheckMalloc ( seeds,   uint,   batchSize * wordPerSeed );
    MC_CheckMalloc ( seedPositions,   int,    inputMaxReadLength );
    clear ();
}
SingleEndSeedingEngine::SingleEndSeedingBatch::~SingleEndSeedingBatch ()
{
    freeSingleAlgnResultArray ( algnResultArray );
    free ( readIDs );
    free ( lengths );
    free ( offsets );
    free ( seeds );
    free ( seedPositions );
}
void SingleEndSeedingEngine::SingleEndSeedingBatch::clear ()
{
    numQueries = 0;
}
inline void SingleEndSeedingEngine::SingleEndSeedingBatch::pack ( uint readID, int off, int seedLength )
{
    readIDs[numQueries] = readID;
    lengths[numQueries] = seedLength;
    offsets[numQueries] = off;
#define MC_OldReadUnpackIn(X,i) ((X[oldReadTPARA + ((i>>4)<<5)] >> ((i & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerQuery + ( readID % 32 );
    uint seedTPARA = ( numQueries / 32 ) * 32 * wordPerSeed + ( numQueries % 32 );

    for ( int i = 0; i < wordPerSeed; i++ )
    {
        seeds[seedTPARA + ( i << 5 )] = 0;
    }

    for ( int i = 0; i < seedLength; i++ )
    {
        int pos = off + i;
        seeds[seedTPARA + ( ( i >> 4 ) << 5 )] |= ( uint ) MC_OldReadUnpackIn ( queries, pos ) << ( ( i & 0xF ) << 1 ) ;
    }

    ++numQueries;
}

int SingleEndSeedingEngine::SingleEndSeedingBatch::packSeeds (
    uint readID, int stage
)
{
    int seedNum, seedLength;
    int readLength = queryLengths[readID];
    getSeedPositions ( stage, readLength, &seedLength, seedPositions, &seedNum );

    if ( numQueries + seedNum > batchSize )
    {
        return 0;
    }

    for ( int i = 0; i < seedNum; i++ )
    {
        pack ( readID, seedPositions[i], seedLength );
    }

    return seedNum;
}

vector<CandidateInfo> * SingleEndSeedingEngine::SingleEndSeedingBatch::singleMerge (
    SeedPos * readPos
)
{
    vector<CandidateInfo> * canInfo = new vector<CandidateInfo> ();
    SeedPos * readIter = readPos;
    uint seedReadCnt = 0;

    while ( true )
    {
        uint readID = readIter->readID;

        if ( readID == 0x7FFFFFFF )
        { break; }

        SeedPos * readStart = readIter;

        while ( readIter->readID == readID )
        { ++readIter; }

        SeedPos * readEnd = readIter;
        canInfo->push_back ( *readStart );
        register uint prevLoc = readStart->pos;

        for ( SeedPos * p = readStart + 1; p < readEnd; p++ )
        {
            register uint curLoc = p->pos;

            if ( prevLoc + DPS_DIVIDE_GAP < curLoc )
            {
                canInfo->push_back ( *p );
                prevLoc = curLoc;
            }
        }

        ++seedReadCnt;
    }

    return canInfo;
}

SeedPos * SingleEndSeedingEngine::SingleEndSeedingBatch::decodePositions (
    BWT * bwt
)
{
#define MC_AppendPos(posIter, id, seedStrand, ePos, off) { \
        posIter->readID = id; \
        posIter->pos = ePos; \
        posIter->strand = seedStrand; \
        ++posIter; \
    }
    // printf("num of raw positions = %llu\n", numOfAnswer);
    SeedPos * pos, *auxPos;
    MC_CheckMalloc ( pos,     SeedPos,    numOfAnswer + 1 );
    MC_CheckMalloc ( auxPos,  SeedPos,    numOfAnswer + 1 );
    SeedPos * iter_pos = pos;

    for ( uint cpuThread = 0; cpuThread < algnResultArray->arrayNum; cpuThread++ )
    {
#define MC_EstimatedPos(x) ( strand == 1 ? \
                             x - offset : \
                             x + seedLength + offset - readLength )
        SingleAlgnResult * result = algnResultArray->array[cpuThread];
        SARecord * iter_sa = result->sa_list;
        OccRecord * iter_occ = result->occ_list;

        for ( uint i = 0; i < result->readNum; i++ )
        {
            uint seedID = result->readIDs[i];
            uint readID = readIDs[seedID];
            uint offset = offsets[seedID];
            uint seedLength = lengths[seedID];
            uint readLength = queryLengths[readID];

            if ( result->saEnds[i] < 0xFFFFFFFF )
            {
                SARecord * end_sa = result->sa_list + result->saEnds[i];

                while ( iter_sa <= end_sa )
                {
                    uint strand = iter_sa->strand;

                    for ( uint k = iter_sa->saLeft; k <= iter_sa->saRight; k++ )
                    {
                        uint estimatedPos = MC_EstimatedPos ( ( *bwt->_bwtSaValue ) ( bwt, k ) );
                        MC_AppendPos ( iter_pos, readID, strand,
                                       estimatedPos, offset );
                    }

                    ++iter_sa;
                }
            }

            if ( result->occEnds[i] < 0xFFFFFFFF )
            {
                OccRecord * end_occ = result->occ_list + result->occEnds[i];

                while ( iter_occ <= end_occ )
                {
                    uint strand = iter_occ->strand;
                    uint estimatedPos = MC_EstimatedPos ( iter_occ->pos );
                    MC_AppendPos ( iter_pos, readID, strand,
                                   estimatedPos, offset );
                    ++iter_occ;
                }
            }
        }
    }

    // array guard
    MC_AppendPos ( iter_pos, 0x7FFFFFFF, 0, 0xFFFFFFFF, 0 );
    uint len = iter_pos - pos;
    MC_RadixSort_32_16 ( pos, pos, auxPos, len );
    MC_RadixSort_32_16 ( pos, readID, auxPos, len );
    MC_RadixSort_8_8 ( pos, strand, auxPos, len );
    free ( auxPos );
    return pos;
}

vector<CandidateInfo> * SingleEndSeedingEngine::SingleEndSeedingBatch::decodeMergePositions (
    BWT * bwt
)
{
    SeedPos * pos = decodePositions ( bwt );
    vector<CandidateInfo> * canInfo = singleMerge ( pos );
    free ( pos );
    return canInfo;
}

// ****
void SingleEndSeedingEngine::SingleEndSeedingThreadContext::init ( SingleEndSeedingBatch * batch )
{
    sem_init ( &ACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    this->batch = batch;
    this->batch->clear ();
}

void SingleEndSeedingEngine::SingleEndSeedingThreadContext::freeMemory ()
{
    delete batch;
}

// ****
SingleEndSeedingEngine::SingleEndSeedingEngine () {}

void SingleEndSeedingEngine::performSeeding ()
{
    cuCtxPopCurrent ( & ( ctx ) );
    seedingSwapBatch =
        new SingleEndSeedingBatch ( DPS_SEEDING_BATCH_SIZE, dpPara,
                                    queries, queryLengths, inputMaxReadLength );
    seedingThreadContext =
        new SingleEndSeedingThreadContext[dpPara->numOfCPUForSeeding];

    for ( int i = 0; i < dpPara->numOfCPUForSeeding; i++ )
    {
        SingleEndSeedingBatch * batch =
            new SingleEndSeedingBatch ( DPS_SEEDING_BATCH_SIZE, dpPara,
                                        queries, queryLengths, inputMaxReadLength );
        seedingThreadContext[i].init ( batch );
    }

    seedingGPUThreadDelegator.init ( 1, SeedingGPUThread,
                                     SeedingGPUThreadInit, SeedingGPUThreadFinalize );
    seedingCPUThreadDelegator.init ( dpPara->numOfCPUForSeeding,
                                     SeedingCPUThread );
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    int threadId;
    void * empty;

    for ( uint i = 0; i < queryIDStream->data->size (); i++ )
    {
        int readID = ( * ( queryIDStream->data ) ) [i];
        inputFlags->set ( readID );

        if ( !seedingSwapBatch->packSeeds ( readID, STAGE_SINGLE_DP ) )
        {
            // launch one batch
            threadId = seedingCPUThreadDelegator.schedule ( empty );
            sem_wait ( & ( seedingThreadContext[threadId].ACKSem ) );
            seedingSwapBatch->clear ();
            seedingSwapBatch->packSeeds ( readID, STAGE_SINGLE_DP );
        }
    }

    // last batch
    if ( seedingSwapBatch->numQueries > 0 )
    {
        threadId = seedingCPUThreadDelegator.schedule ( empty );
        sem_wait ( & ( seedingThreadContext[threadId].ACKSem ) );
    }

    seedingCPUThreadDelegator.finalize ();
    seedingGPUThreadDelegator.finalize ();
    alignFlags->getXOR ( inputFlags, unseededIDStream->data );
    delete inputFlags;
    delete alignFlags;
    delete seedingSwapBatch;

    for ( int i = 0; i < dpPara->numOfCPUForSeeding; i++ )
    {
        seedingThreadContext[i].freeMemory ();
    }

    delete[] seedingThreadContext;
    cuCtxPushCurrent ( ctx );
}

void SingleEndSeedingEngine::performSeeding (
    /* input */
    QueryIDStream    *    queryIDStream,
    DPParameters     *    dpPara,
    uint * queries, uint * queryLengths, int inputMaxReadLength,
    /* soap3 seeding related */
    SOAP3Wrapper<void>  * soap3Wrapper,
    Soap3Index      *     index,
    /* output */
    CandidateStream   *   canStream,
    QueryIDStream    *    unseededIDStream
)
{
    engine = new SingleEndSeedingEngine ();
    MC_MemberCopy5 ( engine->, , queryIDStream, dpPara, queries, queryLengths, inputMaxReadLength );
    MC_MemberCopy4 ( engine->, , soap3Wrapper, index, canStream, unseededIDStream );
    engine->performSeeding ();
    delete engine;
}
SingleEndSeedingEngine * SingleEndSeedingEngine::engine;

void SingleDP_Space::SeedingCPUThread ( int threadId, void *& empty )
{
    SingleEndSeedingEngine * engine = SingleEndSeedingEngine::engine;
    SingleEndSeedingEngine::SingleEndSeedingBatch * batch = engine->seedingSwapBatch;
    engine->seedingSwapBatch = engine->seedingThreadContext[threadId].batch;
    sem_post ( & ( engine->seedingThreadContext[threadId].ACKSem ) );
    engine->seedingThreadContext[threadId].batch = batch;
    int * pThreadId = &threadId;
    engine->seedingGPUThreadDelegator.schedule ( pThreadId );
    sem_wait ( & ( engine->seedingThreadContext[threadId].GPUFinishSem ) );
    vector<CandidateInfo> * canInfo = batch->decodeMergePositions ( engine->index->sraIndex->bwt );
    engine->canStream->append ( canInfo, engine->alignFlags );
    delete canInfo;
}

void SingleDP_Space::SeedingGPUThreadInit ()
{
    cuCtxPushCurrent ( SingleEndSeedingEngine::engine->ctx );
}

void SingleDP_Space::SeedingGPUThread ( int threadId, int *& pCallThreadId )
{
    SingleEndSeedingEngine * engine = SingleEndSeedingEngine::engine;
    SingleEndSeedingEngine::SingleEndSeedingBatch * batch = engine->seedingThreadContext[*pCallThreadId].batch;
    engine->soap3Wrapper->seeding (
        batch->seeds, batch->lengths,
        batch->maxSeedLength, batch->wordPerSeed, batch->batchSize,
        batch->numQueries, batch->numOfAnswer, batch->numOfAlignedRead,
        batch->numOfCPUForSeeding,
        batch->algnResultArray, batch->maxHitNum
    );
    sem_post ( & ( engine->seedingThreadContext[*pCallThreadId].GPUFinishSem ) );
}

void SingleDP_Space::SeedingGPUThreadFinalize ()
{
    cuCtxPopCurrent ( & ( SingleEndSeedingEngine::engine->ctx ) );
}

// ****
SingleEndAlignmentEngine::SingleEndAlgnBatch::SingleEndAlgnBatch (
    int batchSize, DPParameters * dpPara,
    int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
    Soap3Index * index, uint * queries, uint inputMaxReadLength, uint * upkdLengths
)
{
    MC_MemberCopy5 ( this->, , batchSize, maxReadLength, maxDNALength, maxDPTableLength, patternLength );
    MC_MemberCopy3 ( this->, , queries, inputMaxReadLength, upkdLengths );
    MC_MemberCopy2 ( this->, dpPara->, softClipLeft, softClipRight );
    this->cutoffThreshold   = dpPara->paramRead[0].cutoffThreshold;
    this->wordPerOldQuery   = getWordPerQuery ( inputMaxReadLength );
    this->wordPerQuery      = MC_CeilDivide16 ( maxReadLength );
    this->wordPerDNA        = MC_CeilDivide16 ( maxDNALength );
    this->packedDNA         = index->sraIndex->hsp->packedDNA;
    this->fullDNALength     = index->sraIndex->hsp->dnaLength;
    this->index             = index;
    MC_CheckMalloc ( canInfos,            CandidateInfo,  batchSize );
    MC_CheckMalloc ( DNALengths,          uint,           batchSize );
    MC_CheckMalloc ( lengths,             uint,           batchSize );
    MC_CheckMalloc ( packedDNASeq,        uint,           batchSize * MC_CeilDivide16 ( maxDNALength ) );
    MC_CheckMalloc ( packedReadSeq,       uint,           batchSize * MC_CeilDivide16 ( maxReadLength ) );
    MC_CheckMalloc ( scores,              int,            batchSize );
    MC_CheckMalloc ( cutoffThresholds,    int,            batchSize );
    MC_CheckMalloc ( softClipLtSizes,     uint,           batchSize );
    MC_CheckMalloc ( softClipRtSizes,     uint,           batchSize );
    MC_CheckMalloc ( hitLocs,             uint,           batchSize );
    MC_CheckMalloc ( pattern,             uchar,          batchSize * patternLength );
    MC_CheckMalloc ( maxScoreCounts,      uint,           batchSize );
    clear ();
}

SingleEndAlignmentEngine::SingleEndAlgnBatch::~SingleEndAlgnBatch ()
{
    free ( canInfos );
    free ( DNALengths );
    free ( lengths );
    free ( packedDNASeq );
    free ( packedReadSeq );
    free ( scores );
    free ( cutoffThresholds );
    free ( softClipLtSizes );
    free ( softClipRtSizes );
    free ( hitLocs );
    free ( pattern );
    free ( maxScoreCounts );
}

void SingleEndAlignmentEngine::SingleEndAlgnBatch::clear ()
{
    numOfThreads = 0;
}

int SingleEndAlignmentEngine::SingleEndAlgnBatch::pack (
    CandidateInfo & canInfo
)
{
    if ( numOfThreads >= batchSize )
    {
        return 0;
    }

    uint readID = canInfo.readID;
    uint readLength = upkdLengths[readID];
    int margin = DPS_MARGIN ( readLength );
    uint DNAStart = canInfo.pos - margin;

    if ( DNAStart >= fullDNALength )
    {
        DNAStart = 0;
    }

    uint DNALength = readLength + margin * 2;

    if ( DNAStart + DNALength > fullDNALength )
    {
        DNALength = fullDNALength - DNAStart;
    }

    packRead ( packedReadSeq, numOfThreads,
               readID, readLength,
               canInfo.strand );
    repackDNA ( packedDNASeq, numOfThreads,
                packedDNA, DNAStart, DNALength );
    softClipLtSizes[numOfThreads] = ( canInfo.strand == 1 ) ?
                                    softClipLeft : softClipRight;
    softClipRtSizes[numOfThreads] = ( canInfo.strand == 1 ) ?
                                    softClipRight : softClipLeft;
    DNALengths[numOfThreads] = DNALength;
    lengths[numOfThreads] = readLength;
    cutoffThresholds[numOfThreads] = cutoffThreshold;
    canInfo.pos = DNAStart;
    canInfos[numOfThreads] = canInfo;
    ++numOfThreads;
    return 1;
}

inline void SingleEndAlignmentEngine::SingleEndAlgnBatch::packRead (
    uint * packedSeq, uint threadId,
    uint readID, uint length, int strand
)
{
#define MC_OldReadUnpack(X,i) ((X[oldReadTPARA + (((i)>>4)<<5)] >> (((i) & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerOldQuery + ( readID % 32 );
    uint readTPARA = ( threadId / 32 ) * 32 * wordPerQuery + ( threadId % 32 );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[readTPARA + ( i << 5 )] = 0;
    }

    if ( strand == 1 )
    {
        for ( int i = 1; i <= length; i++ )
        {
            int fwd_i = i - 1;
            register uint c_nucleotide = ( uint ) MC_OldReadUnpack ( queries, fwd_i );
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
    else   // strand == 2
    {
        for ( int i = 1; i <= length; i++ )
        {
            int rev_i = length - i;
            register uint c_nucleotide = soap3DnaComplement[ ( uint ) MC_OldReadUnpack ( queries, rev_i )];
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
}

inline void SingleEndAlignmentEngine::SingleEndAlgnBatch::repackDNA (
    uint * packedSeq, uint threadId,
    uint * seq, uint start, uint length
)
{
#define MC_OldDnaUnpack(X,i) ((X[(i)>>4] >> ((15-((i)&0xF))<<1)) & 3)
    uint dnaTPARA = ( threadId / 32 ) * 32 * wordPerDNA + ( threadId & 0x1F );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[dnaTPARA + ( i << 5 )] = 0;
    }

    for ( int i = 1; i <= length; i++ )
    { packedSeq[dnaTPARA + ( ( i >> 4 ) << 5 )] |= ( uint ) ( MC_OldDnaUnpack ( seq, start + i - 1 ) ) << ( ( 15 - ( i & 0xF ) ) << 1 ); }
}

// ****
void SingleEndAlignmentEngine::SingleEndAlgnThreadContext::init ( SingleEndAlgnBatch * batch )
{
    sem_init ( &ACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    sem_init ( &outputACKSem, 0, 0 );
    this->batch = batch;
}

void SingleEndAlignmentEngine::SingleEndAlgnThreadContext::freeMemory ()
{
    delete batch;
}

// ****
SingleEndAlignmentEngine::AlgnmtResultStream::AlgnmtResultStream ()
{
    numOut = 0;
    pthread_mutex_init ( &occupy_mutex, NULL );
}

SingleEndAlignmentEngine::AlgnmtResultStream::~AlgnmtResultStream ()
{
    for ( int i = 0; i < dpSResult.size (); i++ )
    {
        SingleDPResultBatch & resultBatch = * ( dpSResult[i] );

        for ( int j = 0; j < resultBatch.size (); j++ )
        {
            free ( resultBatch[j].cigarString );
        }

        delete dpSResult[i];
    }

    dpSResult.clear ();
}

// ****
void SingleEndAlignmentEngine::performAlignment (
    uint & numDPAlignedRead, uint & numDPAlignment
)
{
    /* initialize */
    cuCtxPopCurrent ( & ( ctx ) );
    algnBatchCount = 0;
    dpSAlignedRead = 0;
    dpSAlignment = 0;
    lastReadID = -1;
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    resultStream = new AlgnmtResultStream;
    outputBuf = new OutputBuffer<SingleAlgnmtResult> ();
    outputBuf->setAlignmentType ( alignmentType );
    maxReadLength = ( inputMaxReadLength / 4 + 1 ) * 4;
    maxDNALength = maxReadLength + 2 * DPS_MARGIN ( inputMaxReadLength ) + 8;
    semiGlobalAligner.decideConfiguration ( maxReadLength, maxDNALength,
                                            maxDPTableLength, DPS_ALGN_NUM_OF_BLOCKS,
                                            patternLength, *dpPara );
    algnSwapBatch =
        new SingleEndAlgnBatch ( DPS_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                 maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                 index, queries, inputMaxReadLength, upkdReadLengths );
    algnThreadContext = new SingleEndAlgnThreadContext[dpPara->numOfCPUThreads];

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        SingleEndAlgnBatch * batch =
            new SingleEndAlgnBatch ( DPS_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                     maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                     index, queries, inputMaxReadLength, upkdReadLengths );
        algnThreadContext[i].init ( batch );
    }

    algnmtGPUThreadDelegator.init ( 1, algnmtGPUThread,
                                    algnmtGPUThreadInit, algnmtGPUThreadFinalize );
    outputThreadDelegator.init ( 1, DPSOutputThread,
                                 NULL, DPSOutputThreadFinalize );
    algnmtCPUThreadDelegator.init ( dpPara->numOfCPUThreads, algnmtCPUThread );
    /* perform alignment */
    int threadId;
    void * empty;

    for ( uint i = 0; i < canStream->data.size (); i++ )
    {
        CandidateInfo & info = canStream->data[i];
        inputFlags->set ( info.readID );

        if ( !algnSwapBatch->pack ( info ) )
        {
            // launch one batch
            threadId = algnmtCPUThreadDelegator.schedule ( empty );
            sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
            algnSwapBatch->clear ();
            algnSwapBatch->pack ( info );
        }
    }

    // last batch
    if ( algnSwapBatch->numOfThreads > 0 )
    {
        threadId = algnmtCPUThreadDelegator.schedule ( empty );
        sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
    }

    /* finalize */
    algnmtCPUThreadDelegator.finalize ();
    algnmtGPUThreadDelegator.finalize ();
    outputThreadDelegator.finalize ();
    alignFlags->getXOR ( inputFlags, unalignedIDStream->data );
    delete inputFlags;
    delete alignFlags;
    delete algnSwapBatch;

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        algnThreadContext[i].freeMemory ();
    }

    delete[] algnThreadContext;
    delete outputBuf;
    delete resultStream;
    numDPAlignedRead = this->dpSAlignedRead;
    numDPAlignment = this->dpSAlignment;
    cuCtxPushCurrent ( ctx );
}

void SingleEndAlignmentEngine::performAlignment (
    /* input */
    CandidateStream   *   canStream,
    DPParameters     *    dpPara,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    char * upkdQueryNames, uint * origReadIDs, char * upkdQualities,
    Soap3Index * index,
    int alignmentType,
    uint accumReadNum, int outputFormat,
    FILE * outputFile, samfile_t * samOutputDPFilePtr,
    /* output */
    QueryIDStream    *    unalignedIDStream,
    uint         &        numDPAlignedRead,
    uint         &        numDPAlignment
)
{
    engine = new SingleEndAlignmentEngine ();
    MC_MemberCopy2 ( engine->, , canStream, dpPara );
    MC_MemberCopy4 ( engine->, , queries, upkdQueryNames, upkdReadLengths, inputMaxReadLength );
    MC_MemberCopy2 ( engine->, , origReadIDs, upkdQualities );
    MC_MemberCopy ( engine->, , index );
    MC_MemberCopy4 ( engine->, , accumReadNum, outputFormat, outputFile, samOutputDPFilePtr );
    MC_MemberCopy2 ( engine->, , alignmentType, unalignedIDStream );
    engine->performAlignment ( numDPAlignedRead, numDPAlignment );
    delete engine;
}

SingleEndAlignmentEngine * SingleEndAlignmentEngine::engine;

void SingleDP_Space::algnmtCPUThread ( int threadId, void *& empty )
{
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    SingleEndAlignmentEngine::SingleEndAlgnBatch * batch = engine->algnSwapBatch;
    engine->algnSwapBatch = engine->algnThreadContext[threadId].batch;
    engine->algnThreadContext[threadId].batchID = engine->algnBatchCount++;
    sem_post ( & ( engine->algnThreadContext[threadId].ACKSem ) );
    engine->algnThreadContext[threadId].batch = batch;
    int * pThreadId = &threadId;
    engine->algnmtGPUThreadDelegator.schedule ( pThreadId );
    sem_wait ( & ( engine->algnThreadContext[threadId].GPUFinishSem ) );
    MC_MemberCopy2 ( int, engine->dpPara->, matchScore, mismatchScore );
    MC_MemberCopy2 ( int, engine->dpPara->, openGapScore, extendGapScore );
    int cutoffThreshold = engine->dpPara->paramRead[0].cutoffThreshold;
    // Rearrange result and Output
    SingleDPResultBatch * resultBatch = new SingleDPResultBatch;

    for ( int i = 0; i < batch->numOfThreads; i++ )
    {
        if ( batch->scores[i] >= cutoffThreshold )
        {
            CigarStringEncoder<void> encoder;
            uchar lastType = 'N';

            for ( uchar * p = batch->pattern + i * engine->patternLength; *p != 0; p++ )
            {
                if ( *p == 'V' )
                {
                    encoder.append ( lastType, ( int ) ( * ( ++p ) ) - 1 );
                }
                else
                {
                    encoder.append ( *p, 1 );
                    lastType = *p;
                }
            }

            SingleAlgnmtResult result;
            result.readID = batch->canInfos[i].readID;
            result.strand = batch->canInfos[i].strand;
            result.algnmt = batch->canInfos[i].pos + batch->hitLocs[i];
            result.score = batch->scores[i];
            encoder.encodeCigarString ( openGapScore, extendGapScore );
            result.cigarString = encoder.cigarString;
            int L = batch->lengths[i] - encoder.charCount['I'] - encoder.charCount['S'];
            int numOfMismatch = ( L * matchScore + encoder.gapPenalty - batch->scores[i] ) /
                                ( matchScore - mismatchScore );
            result.editdist = encoder.charCount['I'] + encoder.charCount['D'] + numOfMismatch;
            result.num_sameScore = batch->maxScoreCounts[i]; // TODO
            resultBatch->push_back ( result );
        }
    }

    // Output
    engine->algnThreadContext[threadId].resultBatch = resultBatch;
    int * pid = &threadId;
    engine->outputThreadDelegator.schedule ( pid );
    // printf("ALgn CPU Thread done.\n");
    sem_wait ( & ( engine->algnThreadContext[threadId].outputACKSem ) );
}

void SingleDP_Space::algnmtGPUThreadInit ()
{
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    cuCtxPushCurrent ( engine->ctx );
    //  showGPUMemInfo("algn enter");
    int batchSize = engine->DPS_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK;
    engine->semiGlobalAligner.init ( batchSize, engine->maxReadLength,
                                     engine->maxDNALength, engine->maxDPTableLength, * ( engine->dpPara ) );
}

void SingleDP_Space::algnmtGPUThread ( int threadId, int *& pCallThreadId )
{
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    SingleEndAlignmentEngine::SingleEndAlgnBatch * batch =
        engine->algnThreadContext[*pCallThreadId].batch;
    engine->semiGlobalAligner.performAlignment (
        batch->packedDNASeq, batch->DNALengths,
        batch->packedReadSeq, batch->lengths,
        batch->cutoffThresholds, batch->scores, batch->hitLocs,
        batch->maxScoreCounts,
        batch->pattern, batch->numOfThreads,
        batch->softClipLtSizes, batch->softClipRtSizes );
    sem_post ( & ( engine->algnThreadContext[*pCallThreadId].GPUFinishSem ) );
}

void SingleDP_Space::algnmtGPUThreadFinalize ()
{
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    engine->semiGlobalAligner.freeMemory ();
    //  showGPUMemInfo("algn exit");
    cuCtxPopCurrent ( & ( engine->ctx ) );
}

void SingleDP_Space::DPSOutputThread ( int threadId, int *& pCallThreadId )
{
    int callThreadId = *pCallThreadId;
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    int batchID = engine->algnThreadContext[callThreadId].batchID;
    SingleDPResultBatch * resultBatch = engine->algnThreadContext[callThreadId].resultBatch;
    sem_post ( & ( engine->algnThreadContext[callThreadId].outputACKSem ) );
    vector<SingleDPResultBatch *> & dpResult = engine->resultStream->dpSResult;

    while ( dpResult.size () <= batchID )
    {
        dpResult.push_back ( NULL );
    }

    dpResult[batchID] = resultBatch;
#define MC_DPSOutputRead() { \
        engine->outputBuf->ready(); \
        if (engine->outputBuf->size > 0) { \
            outputDPSingleResult2( \
                                   engine->outputBuf->elements, engine->outputBuf->size, \
                                   engine->queries, engine->upkdReadLengths, engine->origReadIDs, \
                                   engine->upkdQueryNames, engine->upkdQualities, \
                                   engine->inputMaxReadLength, engine->accumReadNum, engine->outputFormat, \
                                   engine->outputFile, engine->samOutputDPFilePtr, engine->index); \
            engine->dpSAlignedRead += 1; \
            engine->dpSAlignment += engine->outputBuf->size; \
            engine->alignFlags->set(engine->lastReadID); \
        } \
    }
    uint numOut = engine->resultStream->numOut;

    while ( numOut < dpResult.size () && dpResult[numOut] != NULL )
    {
        //OUTPUT HERE
        SingleDPResultBatch & batch = *dpResult[numOut];

        for ( int i = 0; i < batch.size (); i++ )
        {
            SingleAlgnmtResult & result = batch[i];
            int readID = result.readID;

            if ( readID != engine->lastReadID )
            {
                MC_DPSOutputRead ();
                engine->outputBuf->clear ();
                engine->lastReadID = readID;
            }

            engine->outputBuf->add ( result );
        }

        ++numOut;
    }

    engine->resultStream->numOut = numOut;
}

void SingleDP_Space::DPSOutputThreadFinalize ()
{
    SingleEndAlignmentEngine * engine = SingleEndAlignmentEngine::engine;
    MC_DPSOutputRead ();
    engine->outputBuf->clear ();
}

void SingleDP_Space::DPSOutputUnalignedReads (
    QueryIDStream * unalignedIDStream,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    Soap3Index * index,
    char * upkdQueryNames, uint * origReadIDs, char * upkdQualities,
    uint accumReadNum, int outputFormat,
    FILE * outputFile, samfile_t * samOutputDPFilePtr
)
{
    // output unaligned result
#define MC_DPSOutputUnalgnRead() { \
        outputDPSingleResult2(buf, idx, \
                              queries, upkdReadLengths, origReadIDs, \
                              upkdQueryNames, upkdQualities, \
                              inputMaxReadLength, accumReadNum, outputFormat, \
                              outputFile, samOutputDPFilePtr, index); }
    SingleAlgnmtResult * buf;
    MC_CheckMalloc ( buf, SingleAlgnmtResult, 1024 );
    int idx = 0;

    for ( uint i = 0; i < unalignedIDStream->data->size (); i++ )
    {
        buf[idx].readID = ( * ( unalignedIDStream->data ) ) [i];
        buf[idx].algnmt = 0xFFFFFFFF;
        buf[idx].cigarString = NULL;
        ++idx;

        if ( idx >= 1024 )
        {
            MC_DPSOutputUnalgnRead ();
            idx = 0;
        }
    }

    if ( idx > 0 )
    { MC_DPSOutputUnalgnRead (); }

    free ( buf );
}



//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// default-dp space /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
using namespace DP_Space;

// ****
HalfEndOccStream::HalfEndOccStream ( ReadInputForDPArrays * input, BWT * bwt )
{
    this->data = input;
    this->bwt = bwt;
    arrayIndex = 0;
    iter_readInput = data->inputArrays[arrayIndex];
    iter_occ = iter_readInput->occ_list;
    end_occ = iter_occ + iter_readInput->occTotalNum;
    iter_sa = iter_readInput->sa_list;
    end_sa = iter_sa + iter_readInput->saRangeTotalNum;
    nextSAIndex = -1;
}

int HalfEndOccStream::fetchNextOcc ( SRAOccurrence & occ )
{
    while ( true )
    {
#define SA2OCC() { \
        occ.readID = iter_sa->readID; \
        occ.mismatchCount = iter_sa->mismatchCount; \
        occ.strand = iter_sa->strand; \
        occ.ambPosition = (*bwt->_bwtSaValue)(bwt, nextSAIndex++); \
        if (nextSAIndex > iter_sa->saIndexRight) { \
            nextSAIndex = -1; \
            ++iter_sa; \
        } \
    }

        if ( nextSAIndex != -1 )
        {
            SA2OCC ();
        }
        else if ( iter_occ == end_occ )
        {
            if ( iter_sa == end_sa )
            {
                ++arrayIndex;

                if ( arrayIndex >= data->numArrays )
                { return 0; }  // finished
                else
                {
                    iter_readInput = data->inputArrays[arrayIndex];
                    iter_occ = iter_readInput->occ_list;
                    end_occ = iter_occ + iter_readInput->occTotalNum;
                    iter_sa = iter_readInput->sa_list;
                    end_sa = iter_sa + iter_readInput->saRangeTotalNum;
                    continue;
                }
            }
            else
            {
                nextSAIndex = iter_sa->saIndexLeft;
                SA2OCC ();
            }
        }
        else
        {
            if ( iter_sa == end_sa )
            {
                occ = * ( iter_occ++ );
            }
            else
            {
                if ( ( iter_occ->readID >> 1 ) < ( iter_sa->readID >> 1 ) )
                { occ = * ( iter_occ++ ); }
                else
                {
                    nextSAIndex = iter_sa->saIndexLeft;
                    SA2OCC ();
                }
            }
        }

        return 1;
    }
}

// ****
HalfEndAlignmentEngine::HalfEndAlgnBatch::HalfEndAlgnBatch (
    int batchSize, DPParameters * dpPara,
    int peStrandLeftLeg, int peStrandRightLeg, int insert_high, int insert_low,
    int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
    Soap3Index * index, uint * queries, int inputMaxReadLength, uint * upkdReadLengths
)
{
    MC_MemberCopy5 ( this->, , batchSize, peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low );
    MC_MemberCopy4 ( this->, , maxReadLength, maxDNALength, maxDPTableLength, patternLength );
    MC_MemberCopy4 ( this->, , index, queries, inputMaxReadLength, upkdReadLengths );
    MC_MemberCopy2 ( this->, dpPara->, softClipLeft, softClipRight );
    this->cutoffThreshold[0]    = dpPara->paramRead[0].cutoffThreshold;
    this->cutoffThreshold[1]    = dpPara->paramRead[1].cutoffThreshold;
    this->isDoubleStrand        = ( peStrandLeftLeg == peStrandRightLeg );
    this->fullDNALength         = index->sraIndex->hsp->dnaLength;
    this->wordPerOldQuery       = getWordPerQuery ( inputMaxReadLength );
    this->wordPerQuery          = MC_CeilDivide16 ( maxReadLength );
    this->wordPerDNA            = MC_CeilDivide16 ( maxDNALength );
    MC_CheckMalloc ( canInfo,             CandidateInfo,  batchSize );
    MC_CheckMalloc ( DNALengths,          uint,           batchSize );
    MC_CheckMalloc ( lengths,             uint,           batchSize );
    MC_CheckMalloc ( packedDNASequence,   uint,           batchSize * wordPerDNA );
    MC_CheckMalloc ( packedReadSequence,  uint,           batchSize * wordPerQuery );
    MC_CheckMalloc ( startLocs,           uint,           batchSize );
    MC_CheckMalloc ( hitLocs,             uint,           batchSize );
    MC_CheckMalloc ( scores,              int,            batchSize );
    MC_CheckMalloc ( cutoffThresholds,    int,            batchSize );
    MC_CheckMalloc ( softClipLtSizes,     uint,           batchSize );
    MC_CheckMalloc ( softClipRtSizes,     uint,           batchSize );
    MC_CheckMalloc ( peLeftAnchorLocs,    uint,           batchSize );
    MC_CheckMalloc ( peRightAnchorLocs,   uint,           batchSize );
    MC_CheckMalloc ( pattern,             uchar,          batchSize * patternLength );
    MC_CheckMalloc ( maxScoreCounts,      uint,           batchSize );
    clear ();
}

HalfEndAlignmentEngine::HalfEndAlgnBatch::~HalfEndAlgnBatch ()
{
    free ( canInfo );
    free ( DNALengths );
    free ( packedDNASequence );
    free ( lengths );
    free ( packedReadSequence );
    free ( startLocs );
    free ( hitLocs );
    free ( scores );
    free ( cutoffThresholds );
    free ( softClipLtSizes );
    free ( softClipRtSizes );
    free ( peLeftAnchorLocs );
    free ( peRightAnchorLocs );
    free ( pattern );
    free ( maxScoreCounts );
}

void HalfEndAlignmentEngine::HalfEndAlgnBatch::clear ()
{
    numOfThreads = 0;
}

int HalfEndAlignmentEngine::HalfEndAlgnBatch::pack (
    SRAOccurrence & curOcc
)
{
    uint alignedStrand = curOcc.strand;

    if ( alignedStrand != peStrandLeftLeg && alignedStrand != peStrandRightLeg )
    { return numOfThreads; }
    else
    {
        if ( numOfThreads + 1 + isDoubleStrand > batchSize )
        { return -1; }
    }

    uint alignedReadID = curOcc.readID;
    uint alignedPos = curOcc.ambPosition;
    uint alignedReadLength = upkdReadLengths[alignedReadID];
    int  unalignedIsReadOrMate = 1 - ( alignedReadID & 1 );
    uint unalignedReadID = ( unalignedIsReadOrMate == 0 ?
                             alignedReadID - 1 : alignedReadID + 1 );
    uint unalignedReadLength = upkdReadLengths[unalignedReadID];
#define MC_SetRead(strand) { \
        packRead(packedReadSequence, numOfThreads, \
                 unalignedReadID, unalignedReadLength, strand); \
        cutoffThresholds[numOfThreads] = cutoffThreshold[unalignedIsReadOrMate]; \
        softClipLtSizes[numOfThreads] = (strand == 1) ? softClipLeft : softClipRight;  \
        softClipRtSizes[numOfThreads] = (strand == 1) ? softClipRight : softClipLeft; \
    }

    if ( peStrandLeftLeg == alignedStrand )
    {
        //aligned read: at left, unaligned read: at right
        uint rightEnd = alignedPos + insert_high;
        uint rightStart = alignedPos + insert_low - unalignedReadLength;

        // rightStart has to be >= alignedPos
        if ( rightStart < alignedPos )
        { rightStart = alignedPos; }

        if ( rightStart < fullDNALength && rightEnd <= fullDNALength )
        {
            canInfo[numOfThreads].refer = curOcc;
            canInfo[numOfThreads].leftOrRight = 1;
            lengths[numOfThreads] = unalignedReadLength;
            startLocs[numOfThreads] = rightStart;
            DNALengths[numOfThreads] = rightEnd - rightStart;
            peLeftAnchorLocs[numOfThreads] = maxDNALength;
            peRightAnchorLocs[numOfThreads] = unalignedReadLength;
            repackDNA ( packedDNASequence, numOfThreads,
                        ( uint * ) index->sraIndex->hsp->packedDNA, rightStart, DNALengths[numOfThreads] );
            MC_SetRead ( peStrandRightLeg );
            ++numOfThreads;
        }
    }

    if ( peStrandRightLeg == alignedStrand )
    {
        //aligned read: at right, unaligned read: at left
        uint leftStart = alignedPos + alignedReadLength - insert_high;
        uint leftEnd = alignedPos + alignedReadLength - insert_low + unalignedReadLength;

        // leftEnd has to be < alignedPos + alignedReadLength
        if ( leftEnd >= alignedPos + alignedReadLength )
        { leftEnd = alignedPos + alignedReadLength - 1; }

        if ( leftStart < fullDNALength && leftEnd <= fullDNALength )
        {
            canInfo[numOfThreads].refer = curOcc;
            canInfo[numOfThreads].leftOrRight = 0;
            lengths[numOfThreads] = unalignedReadLength;
            startLocs[numOfThreads] = leftStart;
            DNALengths[numOfThreads] = leftEnd - leftStart;
            peLeftAnchorLocs[numOfThreads] = insert_high - insert_low + 1;
            peRightAnchorLocs[numOfThreads] = 0;
            repackDNA ( packedDNASequence, numOfThreads,
                        ( uint * ) index->sraIndex->hsp->packedDNA, leftStart, DNALengths[numOfThreads] );
            MC_SetRead ( peStrandLeftLeg );
            ++numOfThreads;
        }
    }

    return numOfThreads;
}

inline void HalfEndAlignmentEngine::HalfEndAlgnBatch::packRead (
    uint * packedSeq, uint threadId,
    uint readID, uint length, int strand
)
{
#define MC_OldReadUnpack(X,i) ((X[oldReadTPARA + (((i)>>4)<<5)] >> (((i) & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerOldQuery + ( readID % 32 );
    uint readTPARA = ( threadId / 32 ) * 32 * wordPerQuery + ( threadId % 32 );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[readTPARA + ( i << 5 )] = 0;
    }

    if ( strand == 1 )
    {
        for ( int i = 1; i <= length; i++ )
        {
            int fwd_i = i - 1;
            register uint c_nucleotide = ( uint ) MC_OldReadUnpack ( queries, fwd_i );
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
    else   // strand == 2
    {
        for ( int i = 1; i <= length; i++ )
        {
            int rev_i = length - i;
            register uint c_nucleotide = soap3DnaComplement[ ( uint ) MC_OldReadUnpack ( queries, rev_i )];
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
}
inline void HalfEndAlignmentEngine::HalfEndAlgnBatch::repackDNA (
    uint * packedSeq, uint threadId,
    uint * seq, uint start, uint length
)
{
#define MC_OldDnaUnpack(X,i) ((X[(i)>>4] >> ((15-((i)&0xF))<<1)) & 3)
    uint dnaTPARA = ( threadId / 32 ) * 32 * wordPerDNA + ( threadId & 0x1F );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[dnaTPARA + ( i << 5 )] = 0;
    }

    for ( int i = 1; i <= length; i++ )
    {
        packedSeq[dnaTPARA + ( ( i >> 4 ) << 5 )] |= ( uint ) ( MC_OldDnaUnpack ( seq, start + i - 1 ) ) << ( ( 15 - ( i & 0xF ) ) << 1 );
    }
}

// ****
void HalfEndAlignmentEngine::HalfEndAlgnThreadContext::init (
    HalfEndAlgnBatch * batch
)
{
    sem_init ( &dispatchACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    sem_init ( &outputACKSem, 0, 0 );
    this->batch = batch;
}

void HalfEndAlignmentEngine::HalfEndAlgnThreadContext::freeMemory ()
{
    resultBatch = NULL;
    delete batch;
}

// ****
HalfEndAlignmentEngine::AlgnmtResultStream::AlgnmtResultStream ()
{
    numOut = 0;
    pthread_mutex_init ( &occupy_mutex, NULL );
}
HalfEndAlignmentEngine::AlgnmtResultStream::~AlgnmtResultStream ()
{
    for ( int i = 0; i < dpResult.size (); i++ )
    {
        DPResultBatch & resultBatch = * ( dpResult[i] );

        for ( int j = 0; j < resultBatch.size (); j++ )
        {
            free ( resultBatch[j].cigarString );
        }

        delete dpResult[i];
    }

    dpResult.clear ();
}

// ****
void HalfEndAlignmentEngine::performAlignment ( uint & numDPAlignedRead, uint & numDPAlignment )
{
    /* initialize */
    cuCtxPopCurrent ( & ( ctx ) );
    algnBatchCount = 0;
    dpAlignedRead = 0;
    dpAlignment = 0;
    lastReadID = -1;
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    resultStream = new AlgnmtResultStream;
    outputBuf = new OutputBuffer<AlgnmtDPResult> ();
    outputBuf->setAlignmentType ( alignmentType );
    maxReadLength = ( inputMaxReadLength / 4 + 1 ) * 4;
    maxDNALength = insert_high - insert_low + inputMaxReadLength + 1;
    semiGlobalAligner.decideConfiguration ( maxReadLength, maxDNALength,
                                            maxDPTableLength, DP_ALGN_NUM_OF_BLOCKS,
                                            patternLength, *dpPara );
    algnSwapBatch =
        new HalfEndAlgnBatch ( DP_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                               peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                               maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                               index, queries, inputMaxReadLength, upkdReadLengths );
    algnThreadContext = new HalfEndAlgnThreadContext[dpPara->numOfCPUThreads];

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        HalfEndAlgnBatch * batch =
            new HalfEndAlgnBatch ( DP_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                   peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                                   maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                   index, queries, inputMaxReadLength, upkdReadLengths );
        algnThreadContext[i].init ( batch );
    }

    algnmtGPUThreadDelegator.init ( 1, algnmtGPUThread,
                                    algnmtGPUThreadInit, algnmtGPUThreadFinalize );
    outputThreadDelegator.init ( 1, DPOutputThread,
                                 NULL, DPOutputThreadFinalize );
    algnmtCPUThreadDelegator.init ( dpPara->numOfCPUThreads, algnmtCPUThread );
    /* perform alignment */
    int threadId;
    void * empty;
    SRAOccurrence occ;

    while ( canStream->fetchNextOcc ( occ ) )
    {
        inputFlags->set ( ( occ.readID >> 1 ) << 1 );

        if ( algnSwapBatch->pack ( occ ) == -1 )
        {
            threadId = algnmtCPUThreadDelegator.schedule ( empty );
            sem_wait ( & ( algnThreadContext[threadId].dispatchACKSem ) );
            algnSwapBatch->clear ();
            algnSwapBatch->pack ( occ );
        }
    }

    // last batch
    if ( algnSwapBatch->numOfThreads > 0 )
    {
        threadId = algnmtCPUThreadDelegator.schedule ( empty );
        sem_wait ( & ( algnThreadContext[threadId].dispatchACKSem ) );
    }

    /* finalize */
    algnmtCPUThreadDelegator.finalize ();
    algnmtGPUThreadDelegator.finalize ();
    outputThreadDelegator.finalize ();
    alignFlags->getXOR ( inputFlags, unalignedIDStream->data );
    delete inputFlags;
    delete alignFlags;
    delete algnSwapBatch;

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        algnThreadContext[i].freeMemory ();
    }

    delete[] algnThreadContext;
    delete outputBuf;
    delete resultStream;
    numDPAlignedRead = this->dpAlignedRead;
    numDPAlignment = this->dpAlignment;
    cuCtxPushCurrent ( ctx );
}

void HalfEndAlignmentEngine::performAlignment (
    /* input */
    HalfEndOccStream   *  canStream,
    DPParameters     *    dpPara,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    char * upkdQueryNames, uint * origReadIDs, char * upkdQualities,
    Soap3Index * index,
    int alignmentType,
    uint accumReadNum, int outputFormat,
    FILE * outputFile, samfile_t * samOutputDPFilePtr,
    /* output */
    QueryIDStream    *    unalignedIDStream,
    uint         &        numDPAlignedRead,
    uint         &        numDPAlignment
)
{
    engine = new HalfEndAlignmentEngine ();
    MC_MemberCopy2 ( engine->, , canStream, dpPara );
    MC_MemberCopy4 ( engine->, , queries, upkdQueryNames, upkdReadLengths, inputMaxReadLength );
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy2 ( engine->, , origReadIDs, upkdQualities );
    MC_MemberCopy ( engine->, , index );
    MC_MemberCopy4 ( engine->, , accumReadNum, outputFormat, outputFile, samOutputDPFilePtr );
    MC_MemberCopy2 ( engine->, , alignmentType, unalignedIDStream );
    engine->performAlignment ( numDPAlignedRead, numDPAlignment );
    delete engine;
}
HalfEndAlignmentEngine::HalfEndAlignmentEngine * HalfEndAlignmentEngine::engine;

// ****
void DP_Space::algnmtCPUThread ( int threadId, void *& empty )
{
    // Copy data, then ACK to dispatching thread
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    engine->algnThreadContext[threadId].batchID = engine->algnBatchCount++;
    HalfEndAlignmentEngine::HalfEndAlgnBatch * batch = engine->algnSwapBatch;
    engine->algnSwapBatch = engine->algnThreadContext[threadId].batch;
    sem_post ( & ( engine->algnThreadContext[threadId].dispatchACKSem ) );
    engine->algnThreadContext[threadId].batch = batch;
    // launch kernel
    int * pThreadId = &threadId;
    engine->algnmtGPUThreadDelegator.schedule ( pThreadId );
    sem_wait ( & ( engine->algnThreadContext[threadId].GPUFinishSem ) );
    // rearrange result and Output
    MC_MemberCopy2 ( int, engine->, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy2 ( uint *, batch->, startLocs, hitLocs );
    MC_MemberCopy2 ( int, engine->dpPara->, matchScore, mismatchScore );
    MC_MemberCopy2 ( int, engine->dpPara->, openGapScore, extendGapScore );
    uchar * pattern = batch->pattern;
    CandidateInfo * canInfo =  batch->canInfo;
    DPResultBatch * resultBatch = new DPResultBatch;

    for ( int id = 0; id < batch->numOfThreads; id++ )
    {
        //Create record for AlgnmtDPResult;
        AlgnmtDPResult result;
        int alignedID = canInfo[id].refer.readID;
        int alignedIsReadOrMate = alignedID & 1;
        result.readID = alignedID - alignedIsReadOrMate;
        uint dpAlgnmtPos;

        if ( batch->scores[id] >= engine->dpPara->paramRead[1 - alignedIsReadOrMate].cutoffThreshold )
        {
//fprintf ( stderr, "%u %u %u\n", alignedID, batch->scores[id], canInfo[id].refer.strand );
            CigarStringEncoder<void> encoder;
            uchar lastType = 'N';

            for ( uchar * p = pattern + id * engine->patternLength; *p != 0; p++ )
            {
                if ( *p == 'V' )
                {
                    encoder.append ( lastType, ( int ) ( * ( ++p ) ) - 1 );
                }
                else
                {
                    encoder.append ( *p, 1 );
                    lastType = *p;
                }
            }

            encoder.encodeCigarString ( openGapScore, extendGapScore );
            result.cigarString = encoder.cigarString;
            // To get edit distance
            int L = batch->lengths[id] - encoder.charCount['I'] - encoder.charCount['S'];
            int numOfMismatch = ( L * matchScore + encoder.gapPenalty - batch->scores[id] ) /
                                ( matchScore - mismatchScore );
            result.editdist = encoder.charCount['I'] + encoder.charCount['D'] + numOfMismatch;
            result.whichFromDP = 1 - alignedIsReadOrMate;
            dpAlgnmtPos = startLocs[id] + hitLocs[id];

            if ( dpAlgnmtPos < canInfo[id].refer.ambPosition )
            {
                // dp is on left
                result.insertSize = canInfo[id].refer.ambPosition - dpAlgnmtPos +
                                    engine->upkdReadLengths[alignedID];
            }
            else
            {
                // dp is on right
                result.insertSize = dpAlgnmtPos - canInfo[id].refer.ambPosition +
                                    batch->lengths[id] + encoder.charCount['D'] -
                                    encoder.charCount['I'] - encoder.charCount['S'];
            }

            result.num_sameScore = batch->maxScoreCounts[id]; //TODO
        }
        else
        {
            result.cigarString = NULL;
            result.whichFromDP = 2;
            dpAlgnmtPos = 0xFFFFFFFF;
        }

        if ( alignedIsReadOrMate == 0 )
        {
            // aligned is read, unaligned is mate
            result.algnmt_1 = canInfo[id].refer.ambPosition;
            result.algnmt_2 = dpAlgnmtPos;
            result.score_1 = canInfo[id].refer.mismatchCount;
            result.score_2 = batch->scores[id];
            result.strand_1 = canInfo[id].refer.strand;
            result.strand_2 = ( canInfo[id].leftOrRight == 0 ? peStrandLeftLeg : peStrandRightLeg );
        }
        else
        {
            // aligned is mate, unaligned is read
            result.algnmt_1 = dpAlgnmtPos;
            result.algnmt_2 = canInfo[id].refer.ambPosition;
            result.score_1 = batch->scores[id];
            result.score_2 = canInfo[id].refer.mismatchCount;
            result.strand_1 = ( canInfo[id].leftOrRight == 0 ? peStrandLeftLeg : peStrandRightLeg );
            result.strand_2 = canInfo[id].refer.strand;
        }

        resultBatch->push_back ( result );
    }

    // output thread
    engine->algnThreadContext[threadId].resultBatch = resultBatch;
    pThreadId = &threadId;
    engine->outputThreadDelegator.schedule ( pThreadId );
    sem_wait ( & ( engine->algnThreadContext[threadId].outputACKSem ) );
}

void DP_Space::algnmtGPUThreadInit ()
{
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    cuCtxPushCurrent ( engine->ctx );
    int batchSize = engine->DP_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK;
    engine->semiGlobalAligner.init ( batchSize, engine->maxReadLength,
                                     engine->maxDNALength, engine->maxDPTableLength, * ( engine->dpPara ) );
}

void DP_Space::algnmtGPUThreadFinalize ()
{
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    engine->semiGlobalAligner.freeMemory ();
    cuCtxPopCurrent ( & ( engine->ctx ) );
}

void DP_Space::algnmtGPUThread ( int gpuThreadId, int *& pCallThreadId )
{
    int threadId = *pCallThreadId;
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    HalfEndAlignmentEngine::HalfEndAlgnBatch
    *batch = engine->algnThreadContext[threadId].batch;
    //  timeRecorder.appendStart("GPUTime");
    engine->semiGlobalAligner.performAlignment ( batch->packedDNASequence, batch->DNALengths,
            batch->packedReadSequence, batch->lengths,
            batch->cutoffThresholds, batch->scores, batch->hitLocs,
            batch->maxScoreCounts,
            batch->pattern, batch->numOfThreads,
            batch->softClipLtSizes, batch->softClipRtSizes,
            batch->peLeftAnchorLocs, batch->peRightAnchorLocs );
    //  timeRecorder.appendEnd("GPUTime");
    sem_post ( & ( engine->algnThreadContext[threadId].GPUFinishSem ) );
}

void DP_Space::DPOutputThread ( int outputThreadId, int *& pCallThreadId )
{
    int threadId = *pCallThreadId;
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    DPResultBatch * resultBatch = engine->algnThreadContext[threadId].resultBatch;
    int batchID = engine->algnThreadContext[threadId].batchID;
    sem_post ( & ( engine->algnThreadContext[threadId].outputACKSem ) );
    vector<DPResultBatch *> & dpResult = engine->resultStream->dpResult;

    while ( dpResult.size () <= batchID )
    {
        dpResult.push_back ( NULL );
    }

    dpResult[batchID] = resultBatch;
#define MC_OutputRead() { \
        engine->outputBuf->ready(3); \
        if (engine->outputBuf->size > 0 && engine->outputBuf->elements[0].whichFromDP < 2) { \
            outputRead2(engine->outputBuf->elements, engine->outputBuf->size, \
                        engine->queries, engine->upkdReadLengths, \
                        engine->origReadIDs, engine->upkdQueryNames, \
                        engine->upkdQualities, engine->inputMaxReadLength, \
                        engine->accumReadNum, engine->outputFormat, \
                        engine->outputFile, engine->samOutputDPFilePtr, engine->index, \
                        engine->peStrandLeftLeg, engine->peStrandRightLeg); \
            engine->dpAlignedRead += 1; \
            engine->dpAlignment += engine->outputBuf->size; \
            engine->alignFlags->set(engine->lastReadID); \
        } \
    }
    uint numOut = engine->resultStream->numOut;

    while ( numOut < dpResult.size () && dpResult[numOut] != NULL )
    {
        //OUTPUT HERE
        DPResultBatch & batch = *dpResult[numOut];

        for ( int i = 0; i < batch.size (); i++ )
        {
            AlgnmtDPResult & result = batch[i];

            if ( result.readID != engine->lastReadID )
            {
                MC_OutputRead ();
                engine->outputBuf->clear ();
                engine->lastReadID = result.readID;
            }

            engine->outputBuf->add ( result );
        }

        ++numOut;
    }

    engine->resultStream->numOut = numOut;
}

void DP_Space::DPOutputThreadFinalize ()
{
    // last read
    HalfEndAlignmentEngine * engine = HalfEndAlignmentEngine::engine;
    MC_OutputRead ();
    engine->outputBuf->clear ();
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// deep-dp space ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
using namespace DeepDP_Space;
#define DP2_SEEDING_BATCH_SIZE 128 * 1024
#define DP2_MARGIN(l) ((l>100) ? (l>>2) : 25)

// ****
DeepDP_Space::CandidateStream::CandidateStream ()
{
    pthread_mutex_init ( &occupy_mutex, NULL );
}
void DeepDP_Space::CandidateStream::append ( vector<CandidateInfo> * canInfo, AlgnmtFlags * alignFlags )
{
    pthread_mutex_lock ( &occupy_mutex );

    for ( vector<CandidateInfo>::iterator it = canInfo->begin ();
            it < canInfo->end (); ++it )
    {
        data.push_back ( *it );
        alignFlags->set ( ( it->readIDLeft >> 1 ) << 1 );
    }

    pthread_mutex_unlock ( &occupy_mutex );
}

// ****
PairEndSeedingEngine::PairEndSeedingBatch::PairEndSeedingBatch (
    uint batchSize, DPParameters * dpPara,
    uint * queries, uint * readLengths, uint inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg, BWT * bwt
)
{
    MC_MemberCopy4 ( this->, , batchSize, queries, readLengths, bwt );
    MC_MemberCopy4 ( this->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    this->numOfCPUForSeeding = dpPara->numOfCPUForSeeding;

    for ( int i = 0; i < 2; i++ )
    {
        this->maxHitNum[i]  = dpPara->paramRead[i].maxHitNum;
    }

    this->maxSeedLength = inputMaxReadLength;
    this->wordPerSeed   = MC_CeilDivide16 ( maxSeedLength );
    this->wordPerQuery  = getWordPerQuery ( inputMaxReadLength );

    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        algnResultArray[lOr] = constructSingleAlgnResultArray ( numOfCPUForSeeding );
        MC_CheckMalloc ( readIDs[lOr],    uint,   batchSize );
        MC_CheckMalloc ( lengths[lOr],    uint,   batchSize );
        MC_CheckMalloc ( offsets[lOr],    uint,   batchSize );
        MC_CheckMalloc ( seeds[lOr],      uint,   batchSize * wordPerSeed );
    }

    MC_CheckMalloc ( seedPositions,       int,    inputMaxReadLength );
    clear ();
}

PairEndSeedingEngine::PairEndSeedingBatch::~PairEndSeedingBatch ()
{
    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        freeSingleAlgnResultArray ( algnResultArray[lOr] );
        free ( readIDs[lOr] );
        free ( lengths[lOr] );
        free ( offsets[lOr] );
        free ( seeds[lOr] );
    }
}

void PairEndSeedingEngine::PairEndSeedingBatch::clear ()
{
    for ( int i = 0; i < 2; i++ )
    {
        numQueries[i] = 0;
        inPosArr[i].clear ();
    }

    lastPairID = -1;
}
uint PairEndSeedingEngine::PairEndSeedingBatch::findRevStart (
    SeedPos * arr, uint len
)
{
    if ( len == 0 || ! ( arr[len - 1].strand_readID >> 31 ) )
    { return len; }

    uint start = 0;
    uint end = len - 1;

    while ( start < end )
    {
        uint mid = ( start + end ) / 2;

        if ( ( arr[mid].strand_readID >> 31 ) )
        {
            // reverse
            end = mid;
        }
        else
        {
            // forward
            start = mid + 1;
        }
    }

    return start;
}

inline void PairEndSeedingEngine::PairEndSeedingBatch::pack (
    uint evenReadID, uint readID, int off, int seedLength, int readOrMate
)
{
    int seedID = numQueries[readOrMate];
    readIDs[readOrMate][seedID] = evenReadID;
    lengths[readOrMate][seedID] = seedLength;
    offsets[readOrMate][seedID] = off;
#define MC_OldReadUnpackIn(X,i) ((X[oldReadTPARA + ((i>>4)<<5)] >> ((i & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerQuery + ( readID % 32 );
    uint seedTPARA = ( seedID / 32 ) * 32 * wordPerSeed + ( seedID % 32 );

    for ( int i = 0; i < wordPerSeed; i++ )
    {
        seeds[readOrMate][seedTPARA + ( i << 5 )] = 0;
    }

    for ( int i = 0; i < seedLength; i++ )
    {
        int pos = off + i;
        seeds[readOrMate][seedTPARA + ( ( i >> 4 ) << 5 )] |= ( uint ) MC_OldReadUnpackIn ( queries, pos ) << ( ( i & 0xF ) << 1 ) ;
    }

    ++numQueries[readOrMate];
}

int PairEndSeedingEngine::PairEndSeedingBatch::packSeeds (
    uint evenReadID, int stage
)
{
    for ( int i = 0; i < 2; i++ )
    {
        uint readID = evenReadID + i;
        uint readLength = readLengths[readID];
        int seedNum, seedLength;
        getSeedPositions ( stage, readLength, &seedLength, seedPositions, &seedNum );

        if ( numQueries[i] + seedNum > batchSize )
        {
            return 0;
        }

        for ( int j = 0; j < seedNum; j++ )
        {
            pack ( evenReadID, readID, seedPositions[j], seedLength, i );
        }

        // pack(evenReadID, readID, 0, 40, i);
    }

    lastPairID = evenReadID >> 1;
    return 1;
}

int PairEndSeedingEngine::PairEndSeedingBatch::packSeeds (
    SRAOccurrence & occ, int stage
)
{
    uint readID = occ.readID;
    uint pairID = readID >> 1;
    int readOrMate = readID & 1;
    uint evenReadID = readID - readOrMate;

    if ( pairID != lastPairID )
    {
        if ( !packSeeds ( evenReadID, stage ) )
        {
            return 0;
        }
    }

    SeedPos pos;
    pos.pos = occ.ambPosition;
    pos.strand_readID = ( ( occ.strand - 1 ) << 31 ) | evenReadID;
    inPosArr[readOrMate].push_back ( pos );
    return 1;
}

inline int PairEndSeedingEngine::PairEndSeedingBatch::packSeedsOneSide (
    uint evenReadID, int readOrMate, int stage
)
{
    uint readID = evenReadID + readOrMate;
    uint readLength = readLengths[readID];
    int seedNum, seedLength;
    getSeedPositions ( stage, readLength, &seedLength, seedPositions, &seedNum );

    if ( numQueries[readOrMate] + seedNum > batchSize )
    {
        return 0;
    }

    for ( int j = 0; j < seedNum; j++ )
    {
        pack ( evenReadID, readID, seedPositions[j], seedLength, readOrMate );
    }

    lastPairID = evenReadID >> 1;
    return 1;
}

int PairEndSeedingEngine::PairEndSeedingBatch::packSeedsOneSide (
    SRAOccurrence & occ, int stage
)
{
    uint readID = occ.readID;
    uint pairID = readID >> 1;
    int readOrMate = readID & 1;
    uint evenReadID = readID - readOrMate;

    if ( pairID != lastPairID )
    {
        if ( !packSeedsOneSide ( evenReadID, 1 - readOrMate, stage ) )
        {
            return 0;
        }
    }

    SeedPos pos;
    pos.pos = occ.ambPosition;
    pos.strand_readID = ( ( occ.strand - 1 ) << 31 ) | evenReadID;
    inPosArr[readOrMate].push_back ( pos );
    return 1;
}

void PairEndSeedingEngine::PairEndSeedingBatch::pairEndMerge (
    vector<CandidateInfo> * pairEndPos,
    SeedPos * readPos, SeedPos * matePos,
    int leftReadOrMate
)
{
#define MC_DecodePos(x) ((x)->pos)
#define MC_DecodeID(x) ((x)->strand_readID & 0x7FFFFFFF)
#define MC_ReadID() MC_DecodeID(readIter)
#define MC_MateID() MC_DecodeID(mateIter)
    SeedPos * readIter = readPos;
    SeedPos * mateIter = matePos;
    uint readID;
    uint mateID;

    while ( true )
    {
        mateID = MC_MateID ();

        while ( MC_ReadID () < mateID )
        { ++readIter; }

        readID = MC_ReadID ();

        while ( MC_MateID () < readID )
        { ++mateIter; }

        mateID = MC_MateID ();

        if ( mateID == 0x7FFFFFFF )
        { break; }
        else if ( readID < mateID )
        { continue; }

        SeedPos * readStart = readIter;
        SeedPos * mateStart = mateIter;

        // assert : readID == mateID
        while ( MC_ReadID () == readID )
        { ++readIter; }

        while ( MC_MateID () == mateID )
        { ++mateIter; }

        SeedPos * readEnd = readIter;
        SeedPos * mateEnd = mateIter;
#define MC_Compress(start, end, divideGap) { \
        SeedPos *cmprReadIter = start; \
        register uint prevLoc = MC_DecodePos(cmprReadIter); \
        for (SeedPos* p = start+1; p < end; p++) { \
            register uint curLoc = MC_DecodePos(p); \
            if (prevLoc + divideGap < curLoc) { \
                *(++cmprReadIter) = *p; \
                prevLoc = curLoc; \
            } \
        } \
        end = cmprReadIter + 1; \
    }
        MC_Compress ( readStart, readEnd, DP2_DIVIDE_GAP );
        //          MC_Compress(mateStart, mateEnd, DP2_DIVIDE_GAP);
        int readLength = readLengths[readID];
        int margin = DP2_MARGIN ( readLength );
        int length_low = insert_low - readLength - margin;

        if ( length_low < 0 )
        { length_low = 0; }

        int length_high = insert_high - readLength + margin;
        SeedPos * readP = readStart;
        SeedPos * mateP = mateStart;
        register uint readLoc = MC_DecodePos ( readP );
        register uint mateLoc = MC_DecodePos ( mateP );

        while ( readP < readEnd && mateP < mateEnd )
        {
            if ( readLoc + length_low > mateLoc )
            {
                ++mateP;
                mateLoc = MC_DecodePos ( mateP );
            }
            else if ( readLoc + length_high < mateLoc )
            {
                ++readP;
                readLoc = MC_DecodePos ( readP );
            }
            else
            {
                CandidateInfo ci;
                ci.pos[0] = readLoc;
                ci.pos[1] = mateLoc;
                ci.readIDLeft = MC_DecodeID ( readP ) + leftReadOrMate;
                pairEndPos->push_back ( ci );
                //                  ++mateP;
                //                  mateLoc = MC_DecodePos(mateP);
                // TODO
                ++readP;
                readLoc = MC_DecodePos ( readP );
            }
        }
    }
}

int PairEndSeedingEngine::PairEndSeedingBatch::decodePositions (
    int readOrMate, SeedPos *& pos, AlgnmtFlags * tooManyHitFlags
)
{
#define MC_Inc2(x) (x+2)
#define MC_SingleDP_AppendPos(posIter, readID, strandIndex, ePos, off) { \
        posIter->strand_readID = readID | (strandIndex << 31); \
        posIter->pos = ePos; \
        ++posIter; \
    }
    int inPosSize = inPosArr[readOrMate].size ();
    int arrSize = inPosSize + MC_Inc2 ( numOfAnswer[readOrMate] );
    SeedPos * auxPos;
    MC_CheckMalloc ( pos,     SeedPos,    arrSize );
    MC_CheckMalloc ( auxPos,  SeedPos,    arrSize );
    // pre-input value
    copy ( inPosArr[readOrMate].begin (), inPosArr[readOrMate].end (), pos );
    inPosArr[readOrMate].clear ();
    SeedPos * iter_pos = pos + inPosSize;

    for ( uint cpuThread = 0; cpuThread < algnResultArray[readOrMate]->arrayNum; cpuThread++ )
    {
#define MC_SingleDP_EstimatedPos(x) ( \
                                      strandIndex == 0 ? \
                                      x - offset : \
                                      x + seedLength + offset - readLength \
                                    )
        SingleAlgnResult * result = algnResultArray[readOrMate]->array[cpuThread];
        SARecord * iter_sa = result->sa_list;
        OccRecord * iter_occ = result->occ_list;

        for ( uint i = 0; i < result->readNum; i++ )
        {
            uint seedID = result->readIDs[i];
            uint readID = readIDs[readOrMate][seedID];
            uint offset = offsets[readOrMate][seedID];
            uint seedLength = lengths[readOrMate][seedID];
            uint readLength = readLengths[readID + readOrMate];

            if ( result->saEnds[i] < 0xFFFFFFFF )
            {
                SARecord * end_sa = result->sa_list + result->saEnds[i];

                while ( iter_sa <= end_sa )
                {
                    uint strandIndex = iter_sa->strand - 1;

                    for ( uint k = iter_sa->saLeft; k <= iter_sa->saRight; k++ )
                    {
                        uint estimatedPos = MC_SingleDP_EstimatedPos ( ( *bwt->_bwtSaValue ) ( bwt, k ) );
                        MC_SingleDP_AppendPos ( iter_pos, readID, strandIndex,
                                                estimatedPos, offset );
                    }

                    ++iter_sa;
                }
            }

            if ( result->occEnds[i] < 0xFFFFFFFF )
            {
                OccRecord * end_occ = result->occ_list + result->occEnds[i];

                while ( iter_occ <= end_occ )
                {
                    uint strandIndex = iter_occ->strand - 1;
                    uint estimatedPos = MC_SingleDP_EstimatedPos ( iter_occ->pos );
                    MC_SingleDP_AppendPos ( iter_pos, readID, strandIndex,
                                            estimatedPos, offset );
                    ++iter_occ;
                }
            }

            // set tooManyHitflags
            if ( result->isTooManyHit[i] )
            {
                tooManyHitFlags->set ( readID );
            }
        }
    }

    // array guard
    MC_SingleDP_AppendPos ( iter_pos, 0x7FFFFFFF, 0, 0xFFFFFFFF, 0 );
    MC_SingleDP_AppendPos ( iter_pos, 0x7FFFFFFF, 1, 0xFFFFFFFF, 0 );
    uint len = iter_pos - pos;
    MC_RadixSort_32_16 ( pos, pos, auxPos, len );
    MC_RadixSort_32_16 ( pos, strand_readID, auxPos, len );
    free ( auxPos );
    return len;
}

vector<DeepDP_Space::CandidateInfo> * PairEndSeedingEngine::PairEndSeedingBatch::decodeMergePositions (
    AlgnmtFlags * tooManyHitFlags
)
{
    SeedPos * readPos, *matePos;
    uint readPosLen = decodePositions ( 0, readPos, tooManyHitFlags );
    uint matePosLen = decodePositions ( 1, matePos, tooManyHitFlags );
    SeedPos * readArr[2], *mateArr[2];
    // 0 -- forward, 1 -- reverse
    readArr[0] = readPos;
    mateArr[0] = matePos;
    readArr[1] = readPos + findRevStart ( readPos, readPosLen );
    mateArr[1] = matePos + findRevStart ( matePos, matePosLen );
    vector<CandidateInfo> * canInfo = new vector<CandidateInfo>;
    // read left, mate right
    pairEndMerge ( canInfo, readArr[peStrandLeftLeg - 1], mateArr[peStrandRightLeg - 1], 0 );
    pairEndMerge ( canInfo, mateArr[peStrandLeftLeg - 1], readArr[peStrandRightLeg - 1], 1 );
    free ( readPos );
    free ( matePos );
    // Sort the candidates so that readID will be in order
    // To be revised
    vector<CandidateInfo> & candArr = *canInfo;
    uint arrLength = candArr.size ();
    CandidateInfo * auxCandArr;
    MC_CheckMalloc ( auxCandArr, CandidateInfo, arrLength );
    MC_RadixSort_32_16 ( candArr, readIDLeft, auxCandArr, arrLength );
    free ( auxCandArr );
    return canInfo;
}

// ****
void PairEndSeedingEngine::PairEndSeedingThreadContext::init (
    PairEndSeedingBatch * batch
)
{
    sem_init ( &ACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    this->batch = batch;
    this->batch->clear ();
}
void PairEndSeedingEngine::PairEndSeedingThreadContext::freeMemory ()
{
    delete batch;
}

// ****
PairEndSeedingEngine::PairEndSeedingEngine ()
{
    queryIDStream = NULL;
    halfEndOccStream = NULL;
    tooManyHitIDStream = NULL;
}

void PairEndSeedingEngine::performSeeding ()
{
    cuCtxPopCurrent ( & ( ctx ) );
    seedingSwapBatch =
        new PairEndSeedingBatch ( DP2_SEEDING_BATCH_SIZE, dpPara,
                                  queries, queryLengths, inputMaxReadLength,
                                  insert_high, insert_low,
                                  peStrandLeftLeg, peStrandRightLeg, index->sraIndex->bwt );
    seedingThreadContext =
        new PairEndSeedingThreadContext[dpPara->numOfCPUForSeeding];

    for ( int i = 0; i < dpPara->numOfCPUForSeeding; i++ )
    {
        PairEndSeedingBatch * batch =
            new PairEndSeedingBatch ( DP2_SEEDING_BATCH_SIZE, dpPara,
                                      queries, queryLengths, inputMaxReadLength,
                                      insert_high, insert_low,
                                      peStrandLeftLeg, peStrandRightLeg, index->sraIndex->bwt );
        seedingThreadContext[i].init ( batch );
    }

    seedingGPUThreadDelegator.init ( 1, SeedingGPUThread,
                                     SeedingGPUThreadInit, SeedingGPUThreadFinalize );
    seedingCPUThreadDelegator.init ( dpPara->numOfCPUForSeeding,
                                     SeedingCPUThread );
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    tooManyHitFlags = new AlgnmtFlags;
    int threadId;
    void * empty;
    int lastPairID = -1;

    if ( halfEndOccStream != NULL )
    {
        SRAOccurrence occ;

        while ( halfEndOccStream->fetchNextOcc ( occ ) )
        {
            int pairID = occ.readID >> 1;

            if ( lastPairID != pairID )
            {
                inputFlags->set ( pairID << 1 );
                lastPairID = pairID;
            }

            if ( !seedingSwapBatch->packSeeds ( occ, seedingStage ) )
            {
                // launch one batch
                threadId = seedingCPUThreadDelegator.schedule ( empty );
                sem_wait ( & ( seedingThreadContext[threadId].ACKSem ) );
                seedingSwapBatch->clear ();
                seedingSwapBatch->packSeeds ( occ, seedingStage );
            }
        }
    }

    if ( queryIDStream != NULL )
    {
        for ( uint i = 0; i < queryIDStream->data->size (); i++ )
        {
            int readID = ( * ( queryIDStream->data ) ) [i];
            inputFlags->set ( readID );

            if ( !seedingSwapBatch->packSeeds ( readID, seedingStage ) )
            {
                // launch one batch
                threadId = seedingCPUThreadDelegator.schedule ( empty );
                sem_wait ( & ( seedingThreadContext[threadId].ACKSem ) );
                seedingSwapBatch->clear ();
                seedingSwapBatch->packSeeds ( readID, seedingStage );
            }
        }
    }

    // last batch
    threadId = seedingCPUThreadDelegator.schedule ( empty );
    sem_wait ( & ( seedingThreadContext[threadId].ACKSem ) );
    seedingCPUThreadDelegator.finalize ();
    seedingGPUThreadDelegator.finalize ();

    if ( tooManyHitIDStream == NULL )
    {
        // put together
        alignFlags->getXOR ( inputFlags, unseededIDStream->data );
    }
    else
    {
        // should separate
        // step 1, get unaligned reads
        alignFlags->XOR ( inputFlags );
        // step 2, get unaligned & tooManyHitReads
        tooManyHitFlags->AND ( alignFlags );
        // step 3
        alignFlags->getXOR ( tooManyHitFlags, unseededIDStream->data );
        tooManyHitFlags->get ( tooManyHitIDStream->data );
    }

    delete inputFlags;
    delete alignFlags;
    delete tooManyHitFlags;
    delete seedingSwapBatch;

    for ( int i = 0; i < dpPara->numOfCPUForSeeding; i++ )
    {
        seedingThreadContext[i].freeMemory ();
    }

    delete[] seedingThreadContext;
    cuCtxPushCurrent ( ctx );
}

void PairEndSeedingEngine::performSeeding (
    /* input */
    QueryIDStream    *    queryIDStream,
    DPParameters     *    dpPara,
    uint * queries, uint * queryLengths, int inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    int seedingStage,
    /* soap3 seeding related */
    SOAP3Wrapper<void>  * soap3Wrapper,
    Soap3Index      *     index,
    /* output */
    CandidateStream   *   canStream,
    QueryIDStream    *    unseededIDStream
)
{
    engine = new PairEndSeedingEngine ();
    MC_MemberCopy5 ( engine->, , queryIDStream, dpPara, queries, queryLengths, inputMaxReadLength );
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy4 ( engine->, , soap3Wrapper, index, canStream, unseededIDStream );
    engine->seedingStage = seedingStage;
    engine->performSeeding ();
    delete engine;
}

void PairEndSeedingEngine::performSeeding (
    /* input */
    QueryIDStream    *    queryIDStream,
    DPParameters     *    dpPara,
    uint * queries, uint * queryLengths, int inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    int seedingStage,
    /* soap3 seeding related */
    SOAP3Wrapper<void>  * soap3Wrapper,
    Soap3Index      *     index,
    /* output */
    CandidateStream   *   canStream,
    QueryIDStream    *    tooManyHitIDStream,
    QueryIDStream    *    unseededIDStream
)
{
    engine = new PairEndSeedingEngine ();
    MC_MemberCopy5 ( engine->, , queryIDStream, dpPara, queries, queryLengths, inputMaxReadLength );
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy5 ( engine->, , soap3Wrapper, index, canStream, tooManyHitIDStream, unseededIDStream );
    engine->seedingStage = seedingStage;
    engine->performSeeding ();
    delete engine;
}


void PairEndSeedingEngine::performSeeding (
    /* input */
    DP_Space::HalfEndOccStream
    * halfEndOccStream,
    DPParameters     *    dpPara,
    uint * queries, uint * queryLengths, int inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    int seedingStage,
    /* soap3 seeding related */
    SOAP3Wrapper<void>  * soap3Wrapper,
    Soap3Index      *     index,
    /* output */
    CandidateStream   *   canStream,
    QueryIDStream    *    unseededIDStream
)
{
    engine = new PairEndSeedingEngine ();
    MC_MemberCopy5 ( engine->, , halfEndOccStream, dpPara, queries, queryLengths, inputMaxReadLength );
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy4 ( engine->, , soap3Wrapper, index, canStream, unseededIDStream );
    engine->seedingStage = seedingStage;
    engine->performSeeding ();
    delete engine;
}

PairEndSeedingEngine * PairEndSeedingEngine::engine;

void DeepDP_Space::SeedingGPUThreadInit ()
{
    cuCtxPushCurrent ( PairEndSeedingEngine::engine->ctx );
    //  showGPUMemInfo("seeding enter");
}
void DeepDP_Space::SeedingGPUThread ( int threadId, int *& pCallThreadId )
{
    PairEndSeedingEngine * engine = PairEndSeedingEngine::engine;
    PairEndSeedingEngine::PairEndSeedingBatch * batch =
        engine->seedingThreadContext[*pCallThreadId].batch;

    for ( int r = 0; r < 2; r++ )
    {
        engine->soap3Wrapper->seeding (
            batch->seeds[r], batch->lengths[r],
            batch->maxSeedLength, batch->wordPerSeed, batch->batchSize,
            batch->numQueries[r], batch->numOfAnswer[r], batch->numOfAlignedRead[r],
            batch->numOfCPUForSeeding,
            batch->algnResultArray[r], batch->maxHitNum[r]
        );
    }

    sem_post ( & ( engine->seedingThreadContext[*pCallThreadId].GPUFinishSem ) );
}
void DeepDP_Space::SeedingGPUThreadFinalize ()
{
    //  showGPUMemInfo("seeding exit");
    cuCtxPopCurrent ( & ( PairEndSeedingEngine::engine->ctx ) );
}
void DeepDP_Space::SeedingCPUThread ( int threadId, void *& empty )
{
    PairEndSeedingEngine * engine = PairEndSeedingEngine::engine;
    PairEndSeedingEngine::PairEndSeedingBatch * batch = engine->seedingSwapBatch;
    engine->seedingSwapBatch = engine->seedingThreadContext[threadId].batch;
    // printf("[%u][%d] Launching a seeding batch... ", threadId, engine->seedingBatchCount++); fflush(stdout);
    sem_post ( & ( engine->seedingThreadContext[threadId].ACKSem ) );
    engine->seedingThreadContext[threadId].batch = batch;
    int * pThreadId = &threadId;
    engine->seedingGPUThreadDelegator.schedule ( pThreadId );
    sem_wait ( & ( engine->seedingThreadContext[threadId].GPUFinishSem ) );
    // printf("[%u] Seeding thread done.\n", threadId); fflush(stdout);
    vector<CandidateInfo> * candidates = batch->decodeMergePositions ( engine->tooManyHitFlags );
    engine->canStream->append ( candidates, engine->alignFlags );
    delete candidates;
}

void DeepDP_Space::DP2OutputUnalignedReads (
    QueryIDStream * unalignedIDStream,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg,
    char * upkdQueryNames, uint * origReadIDs, char * upkdQualities,
    uint accumReadNum, int outputFormat,
    FILE * outputFile, samfile_t * samOutputDPFilePtr
)
{
    // output unaligned result
#define MC_DP2OutputUnalgnRead() { \
        outputDeepDPResult2(buf, idx, \
                            queries, upkdReadLengths, \
                            origReadIDs, upkdQueryNames, upkdQualities, \
                            inputMaxReadLength, accumReadNum, outputFormat, \
                            outputFile, samOutputDPFilePtr, index, \
                            peStrandLeftLeg, peStrandRightLeg); }
    DeepDPAlignResult * buf;
    MC_CheckMalloc ( buf, DeepDPAlignResult, 1024 );
    int idx = 0;

    for ( uint i = 0; i < unalignedIDStream->data->size (); i++ )
    {
        buf[idx].readID = ( * ( unalignedIDStream->data ) ) [i];
        buf[idx].algnmt_1 = 0xFFFFFFFF;
        buf[idx].algnmt_2 = 0xFFFFFFFF;
        buf[idx].cigarString_1 = NULL;
        buf[idx].cigarString_2 = NULL;
        ++idx;

        if ( idx >= 1024 )
        {
            MC_DP2OutputUnalgnRead ();
            idx = 0;
        }
    }

    if ( idx > 0 )
    { MC_DP2OutputUnalgnRead (); }

    free ( buf );
}

// ****
PairEndAlignmentEngine::PairEndAlgnBatch::PairEndAlgnBatch (
    int batchSize, DPParameters * dpPara,
    int peStrandLeftLeg, int peStrandRightLeg, int insert_high, int insert_low,
    int maxReadLength, int maxDNALength, int maxDPTableLength, int patternLength,
    Soap3Index * index, uint * queries, uint inputMaxReadLength, uint * upkdLengths
)
{
    MC_MemberCopy5 ( this->, , batchSize, maxReadLength, maxDNALength, maxDPTableLength, patternLength );
    MC_MemberCopy4 ( this->, , peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low );
    MC_MemberCopy3 ( this->, , queries, inputMaxReadLength, upkdLengths );
    MC_MemberCopy2 ( this->, dpPara->, softClipLeft, softClipRight );
    this->cutoffThreshold[0]    = dpPara->paramRead[0].cutoffThreshold;
    this->cutoffThreshold[1]    = dpPara->paramRead[1].cutoffThreshold;
    this->wordPerOldQuery   = getWordPerQuery ( inputMaxReadLength );
    this->wordPerQuery      = MC_CeilDivide16 ( maxReadLength );
    this->wordPerDNA        = MC_CeilDivide16 ( maxDNALength );
    this->packedDNA         = index->sraIndex->hsp->packedDNA;
    this->fullDNALength     = index->sraIndex->hsp->dnaLength;
    this->index             = index;
    MC_CheckMalloc ( packedDNASeq,        uint,           batchSize * wordPerDNA );
    MC_CheckMalloc ( packedReadSeq,       uint,           batchSize * wordPerQuery );
    MC_CheckMalloc ( canInfos,            CandidateInfo,  batchSize );
    MC_CheckMalloc ( DNALengths,          uint,           batchSize );
    MC_CheckMalloc ( lengths,             uint,           batchSize );
    MC_CheckMalloc ( cutoffThresholds,    int,            batchSize );

    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        MC_CheckMalloc ( scores[lOr],     int,            batchSize );
        MC_CheckMalloc ( hitLocs[lOr],    uint,           batchSize );
        MC_CheckMalloc ( pattern[lOr],    uchar,          batchSize * patternLength );
        MC_CheckMalloc ( maxScoreCounts[lOr],     uint,           batchSize );
    }

    MC_CheckMalloc ( softClipLtSizes,     uint,           batchSize );
    MC_CheckMalloc ( softClipRtSizes,     uint,           batchSize );
    MC_CheckMalloc ( peLeftAnchorLocs,    uint,           batchSize );
    MC_CheckMalloc ( peRightAnchorLocs,   uint,           batchSize );
    clear ();
}
PairEndAlignmentEngine::PairEndAlgnBatch::~PairEndAlgnBatch ()
{
    free ( packedDNASeq );
    free ( packedReadSeq );
    free ( canInfos );
    free ( DNALengths );
    free ( lengths );
    free ( cutoffThresholds );

    for ( int lOr = 0; lOr < 2; lOr++ )
    {
        free ( scores[lOr] );
        free ( hitLocs[lOr] );
        free ( pattern[lOr] );
        free ( maxScoreCounts[lOr] );
    }

    free ( softClipLtSizes );
    free ( softClipRtSizes );
    free ( peLeftAnchorLocs );
    free ( peRightAnchorLocs );
}
void PairEndAlignmentEngine::PairEndAlgnBatch::clear ()
{
    numOfThreads = 0;
}

int PairEndAlignmentEngine::PairEndAlgnBatch::packLeft (
    CandidateInfo & canInfo
)
{
    if ( numOfThreads >= batchSize )
    {
        return 0;
    }

    uint readIDLeft = canInfo.readIDLeft;
    uint readLength = upkdLengths[readIDLeft];
    int margin = DP2_MARGIN ( readLength );
    uint DNAStartLeft = canInfo.pos[0] - margin;

    if ( DNAStartLeft >= fullDNALength )
    {
        DNAStartLeft = 0;
    }

    uint DNALength = readLength + margin * 2;

    if ( DNAStartLeft + DNALength > fullDNALength )
    {
        DNALength = fullDNALength - DNAStartLeft;
    }

    // no anchor requirement
    peLeftAnchorLocs[numOfThreads] = maxDNALength;
    peRightAnchorLocs[numOfThreads] = 0;
    packRead ( packedReadSeq, numOfThreads,
               readIDLeft, readLength, peStrandLeftLeg );
    repackDNA ( packedDNASeq, numOfThreads,
                packedDNA, DNAStartLeft, DNALength );
    softClipLtSizes[numOfThreads] = ( peStrandLeftLeg == 1 ) ?
                                    softClipLeft : softClipRight;
    softClipRtSizes[numOfThreads] = ( peStrandLeftLeg == 1 ) ?
                                    softClipRight : softClipLeft;
    DNALengths[numOfThreads] = DNALength;
    lengths[numOfThreads] = readLength;
    cutoffThresholds[numOfThreads] = cutoffThreshold[readIDLeft & 1];
    canInfo.pos[0] = DNAStartLeft;
    canInfos[numOfThreads] = canInfo;
    ++numOfThreads;
    return 1;
}

void PairEndAlignmentEngine::PairEndAlgnBatch::packRight ()
{
    for ( int i = 0; i < numOfThreads; i++ )
    {
        uint readIDLeft = canInfos[i].readIDLeft;
        uint leftIsOdd = readIDLeft & 1;

        if ( scores[0][i] >= cutoffThreshold[leftIsOdd] )
        {
            uint readIDLeft = canInfos[i].readIDLeft;
            uint readIDRight = ( leftIsOdd ) ? ( readIDLeft - 1 ) : ( readIDLeft + 1 );
            uint readLength = upkdLengths[readIDRight];
            uint margin = DP2_MARGIN ( readLength );
            uint DNAStartRight = canInfos[i].pos[1] - margin;

            if ( DNAStartRight >= fullDNALength )
            {
                DNAStartRight = 0;
            }

            uint DNALength = readLength + margin * 2;

            if ( DNAStartRight + DNALength > fullDNALength )
            {
                DNALength = fullDNALength - DNAStartRight;
            }

            uint hitPosLeft = canInfos[i].pos[0] + hitLocs[0][i];
            // restrict maximum insert size
            uint boundedLength = hitPosLeft + insert_high - DNAStartRight;

            if ( boundedLength < DNALength ) \
                DNALength = boundedLength;

            // set pair-end anchor boundary, restrict minimum insert size
            peLeftAnchorLocs[i] = maxDNALength;
            int rightAnchor = hitPosLeft + insert_low - DNAStartRight;
            peRightAnchorLocs[i] = rightAnchor > 0 ? rightAnchor : 0;
            packRead ( packedReadSeq, i,
                       readIDRight, readLength, peStrandRightLeg );
            repackDNA ( packedDNASeq, i,
                        packedDNA, DNAStartRight, DNALength );
            softClipLtSizes[i] = ( peStrandRightLeg == 1 ) ?
                                 softClipLeft : softClipRight;
            softClipRtSizes[i] = ( peStrandRightLeg == 1 ) ?
                                 softClipRight : softClipLeft;
            DNALengths[i] = DNALength;
            lengths[i] = readLength;
            cutoffThresholds[i] = cutoffThreshold[readIDRight & 1];
            canInfos[i].pos[1] = DNAStartRight;
        }
    }
}

inline void PairEndAlignmentEngine::PairEndAlgnBatch::packRead (
    uint * packedSeq, uint threadId,
    uint readID, uint length, int strand
)
{
#define MC_OldReadUnpack(X,i) ((X[oldReadTPARA + (((i)>>4)<<5)] >> (((i) & 0xF) << 1)) & 0x3)
    uint oldReadTPARA = ( readID / 32 ) * 32 * wordPerOldQuery + ( readID % 32 );
    uint readTPARA = ( threadId / 32 ) * 32 * wordPerQuery + ( threadId % 32 );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[readTPARA + ( i << 5 )] = 0;
    }

    if ( strand == 1 )
    {
        for ( int i = 1; i <= length; i++ )
        {
            int fwd_i = i - 1;
            register uint c_nucleotide = ( uint ) MC_OldReadUnpack ( queries, fwd_i );
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
    else   // strand == 2
    {
        for ( int i = 1; i <= length; i++ )
        {
            int rev_i = length - i;
            register uint c_nucleotide = soap3DnaComplement[ ( uint ) MC_OldReadUnpack ( queries, rev_i )];
#ifdef BS_MOD
            c_nucleotide = c_nucleotide ^ ( ( c_nucleotide == index->sraIndex->hsp->flag ) << 1 );
#endif
            packedSeq[readTPARA + ( ( i >> 4 ) << 5 )] |= c_nucleotide << ( ( 15 - ( i & 0xF ) ) << 1 );
        }
    }
}

inline void PairEndAlignmentEngine::PairEndAlgnBatch::repackDNA (
    uint * packedSeq, uint threadId,
    uint * seq, uint start, uint length
)
{
#define MC_OldDnaUnpack(X,i) ((X[(i)>>4] >> ((15-((i)&0xF))<<1)) & 3)
    uint dnaTPARA = ( threadId / 32 ) * 32 * wordPerDNA + ( threadId & 0x1F );

    for ( uint i = 0; i <= ( length / CHAR_PER_WORD ); i++ )
    {
        packedSeq[dnaTPARA + ( i << 5 )] = 0;
    }

    for ( int i = 1; i <= length; i++ )
    { packedSeq[dnaTPARA + ( ( i >> 4 ) << 5 )] |= ( uint ) ( MC_OldDnaUnpack ( seq, start + i - 1 ) ) << ( ( 15 - ( i & 0xF ) ) << 1 ); }
}

// ****
void PairEndAlignmentEngine::PairEndAlgnThreadContext::init (
    PairEndAlgnBatch * batch
)
{
    sem_init ( &ACKSem, 0, 0 );
    sem_init ( &GPUFinishSem, 0, 0 );
    sem_init ( &outputACKSem, 0, 0 );
    this->batch = batch;
}
void PairEndAlignmentEngine::PairEndAlgnThreadContext::freeMemory ()
{
    delete batch;
}

// ****
PairEndAlignmentEngine::AlgnmtResultStream::AlgnmtResultStream ()
{
    numOut = 0;
    pthread_mutex_init ( &occupy_mutex, NULL );
}

PairEndAlignmentEngine::AlgnmtResultStream::~AlgnmtResultStream ()
{
    for ( int i = 0; i < dp2Result.size (); i++ )
    {
        DP2ResultBatch & resultBatch = * ( dp2Result[i] );

        for ( int j = 0; j < resultBatch.size (); j++ )
        {
            free ( resultBatch[j].cigarString_1 );
            free ( resultBatch[j].cigarString_2 );
        }

        delete dp2Result[i];
    }

    dp2Result.clear ();
}

// ****
PairEndAlignmentEngine::PairEndAlignmentEngine () {}

void PairEndAlignmentEngine::performAlignment ( uint & numDPAlignedRead, uint & numDPAlignment )
{
    /* initialize */
    cuCtxPopCurrent ( & ( ctx ) );
    algnBatchCount = 0;
    dp2AlignedRead = 0;
    dp2Alignment = 0;
    lastReadID = -1;
    inputFlags = new AlgnmtFlags;
    alignFlags = new AlgnmtFlags;
    resultStream = new AlgnmtResultStream;
    outputBuf = new OutputBuffer<DeepDPAlignResult> ();
    outputBuf->setAlignmentType ( alignmentType );
    maxReadLength = ( inputMaxReadLength / 4 + 1 ) * 4;
    maxDNALength = maxReadLength + 2 * DP2_MARGIN ( inputMaxReadLength ) + 8;
    semiGlobalAligner.decideConfiguration ( maxReadLength, maxDNALength,
                                            maxDPTableLength, DP2_ALGN_NUM_OF_BLOCKS,
                                            patternLength, *dpPara );
    algnSwapBatch =
        new PairEndAlgnBatch ( DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                               peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                               maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                               index, queries, inputMaxReadLength, upkdReadLengths );
    algnThreadContext = new PairEndAlgnThreadContext[dpPara->numOfCPUThreads];

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        PairEndAlgnBatch * batch =
            new PairEndAlgnBatch ( DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK, dpPara,
                                   peStrandLeftLeg, peStrandRightLeg, insert_high, insert_low,
                                   maxReadLength, maxDNALength, maxDPTableLength, patternLength,
                                   index, queries, inputMaxReadLength, upkdReadLengths );
        algnThreadContext[i].init ( batch );
    }

    algnmtGPUThreadDelegator.init ( 1, DP2GPUAlgnThread,
                                    DP2GPUAlgnThreadInit, DP2GPUAlgnThreadFinalize );
    outputThreadDelegator.init ( 1, DP2OutputThread,
                                 NULL, DP2OutputThreadFinalize );
    algnmtCPUThreadDelegator.init ( dpPara->numOfCPUThreads, DP2CPUAlgnThread );
    /* perform alignment */
    int threadId;
    void * empty;

    for ( uint i = 0; i < canStream->data.size (); i++ )
    {
        CandidateInfo & info = canStream->data[i];
        inputFlags->set ( ( info.readIDLeft >> 1 ) << 1 );

        if ( !algnSwapBatch->packLeft ( info ) )
        {
            // launch one batch
            threadId = algnmtCPUThreadDelegator.schedule ( empty );
            sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
            algnSwapBatch->clear ();
            algnSwapBatch->packLeft ( info );
        }
    }

    // last batch
    if ( algnSwapBatch->numOfThreads > 0 )
    {
        threadId = algnmtCPUThreadDelegator.schedule ( empty );
        sem_wait ( & ( algnThreadContext[threadId].ACKSem ) );
    }

    /* finalize */
    algnmtCPUThreadDelegator.finalize ();
    algnmtGPUThreadDelegator.finalize ();
    outputThreadDelegator.finalize ();
    alignFlags->getXOR ( inputFlags, unalignedIDStream->data );
    delete inputFlags;
    delete alignFlags;
    delete algnSwapBatch;

    for ( int i = 0; i < dpPara->numOfCPUThreads; i++ )
    {
        algnThreadContext[i].freeMemory ();
    }

    delete[] algnThreadContext;
    delete outputBuf;
    delete resultStream;
    numDPAlignedRead = this->dp2AlignedRead;
    numDPAlignment = this->dp2Alignment;
    cuCtxPushCurrent ( ctx );
}

void PairEndAlignmentEngine::performAlignment (
    /* input */
    CandidateStream   *   canStream,
    DPParameters     *    dpPara,
    uint * queries, uint * upkdReadLengths, int inputMaxReadLength,
    int insert_high, int insert_low,
    int peStrandLeftLeg, int peStrandRightLeg,
    char * upkdQueryNames, uint * origReadIDs, char * upkdQualities,
    Soap3Index * index,
    int alignmentType,
    uint accumReadNum, int outputFormat,
    FILE * outputFile, samfile_t * samOutputDPFilePtr,
    /* output */
    QueryIDStream    *    unalignedIDStream,
    uint         &        numDPAlignedRead,
    uint         &        numDPAlignment
)
{
    engine = new PairEndAlignmentEngine ();
    MC_MemberCopy2 ( engine->, , canStream, dpPara );
    MC_MemberCopy4 ( engine->, , queries, upkdQueryNames, upkdReadLengths, inputMaxReadLength );
    MC_MemberCopy4 ( engine->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
    MC_MemberCopy2 ( engine->, , origReadIDs, upkdQualities );
    MC_MemberCopy ( engine->, , index );
    MC_MemberCopy4 ( engine->, , accumReadNum, outputFormat, outputFile, samOutputDPFilePtr );
    MC_MemberCopy2 ( engine->, , alignmentType, unalignedIDStream );
    engine->performAlignment ( numDPAlignedRead, numDPAlignment );
    delete engine;
}
PairEndAlignmentEngine * PairEndAlignmentEngine::engine;

// ****
void DeepDP_Space::DP2GPUAlgnThreadInit ()
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    cuCtxPushCurrent ( engine->ctx );
    //showGPUMemInfo("algn enter");
    int batchSize = engine->DP2_ALGN_NUM_OF_BLOCKS * DP_THREADS_PER_BLOCK;
    engine->semiGlobalAligner.init ( batchSize, engine->maxReadLength,
                                     engine->maxDNALength, engine->maxDPTableLength, * ( engine->dpPara ) );
}

void DeepDP_Space::DP2GPUAlgnThread ( int threadId, int *& pCallThreadId )
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    PairEndAlignmentEngine::PairEndAlgnBatch * batch =
        engine->algnThreadContext[*pCallThreadId].batch;
    //  timeRecorder.appendStart("DP2");
    int lOr = batch->leftOrRight;
    engine->semiGlobalAligner.performAlignment (
        batch->packedDNASeq, batch->DNALengths,
        batch->packedReadSeq, batch->lengths,
        batch->cutoffThresholds, batch->scores[lOr], batch->hitLocs[lOr],
        batch->maxScoreCounts[lOr],
        batch->pattern[lOr], batch->numOfThreads,
        batch->softClipLtSizes, batch->softClipRtSizes,
        batch->peLeftAnchorLocs, batch->peRightAnchorLocs
    );
    //  timeRecorder.appendEnd("DP2");
    sem_post ( & ( engine->algnThreadContext[*pCallThreadId].GPUFinishSem ) );
}

void DeepDP_Space::DP2GPUAlgnThreadFinalize ()
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    engine->semiGlobalAligner.freeMemory ();
    cuCtxPopCurrent ( & ( engine->ctx ) );
}

void DeepDP_Space::DP2CPUAlgnThread ( int threadId, void *& empty )
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    PairEndAlignmentEngine::PairEndAlgnBatch * batch = engine->algnSwapBatch;
    engine->algnSwapBatch = engine->algnThreadContext[threadId].batch;
    engine->algnThreadContext[threadId].batchID = engine->algnBatchCount++;
    sem_post ( & ( engine->algnThreadContext[threadId].ACKSem ) );
    engine->algnThreadContext[threadId].batch = batch;
    int * pThreadId = &threadId;
    // align left side
    batch->leftOrRight = 0;
    engine->algnmtGPUThreadDelegator.schedule ( pThreadId );
    sem_wait ( & ( engine->algnThreadContext[threadId].GPUFinishSem ) );
    // align right side
    batch->packRight ();
    batch->leftOrRight = 1;
    engine->algnmtGPUThreadDelegator.schedule ( pThreadId );
    sem_wait ( & ( engine->algnThreadContext[threadId].GPUFinishSem ) );
    MC_MemberCopy2 ( int, engine->dpPara->, matchScore, mismatchScore );
    MC_MemberCopy2 ( int, engine->dpPara->, openGapScore, extendGapScore );
    int cutoffThreshold[2];
    cutoffThreshold[0] = engine->dpPara->paramRead[0].cutoffThreshold;
    cutoffThreshold[1] = engine->dpPara->paramRead[1].cutoffThreshold;
    // rearrange result and Output
    vector<DeepDPAlignResult> * resultBatch = new vector<DeepDPAlignResult> ();

    for ( int i = 0; i < batch->numOfThreads; i++ )
    {
        int readSide = batch->canInfos[i].readIDLeft & 1;
        int mateSide = 1 - readSide;

        if ( batch->scores[0][i] >= cutoffThreshold[readSide] &&
                batch->scores[1][i] >= cutoffThreshold[mateSide] )
        {
            char * cigarString[2];
            int editdist[2], DIS[2];

            for ( int lOr = 0; lOr < 2; lOr++ )
            {
                CigarStringEncoder<void> encoder;
                uchar lastType = 'N';

                for ( uchar * p = batch->pattern[lOr] + i * engine->patternLength; *p != 0; p++ )
                {
                    if ( *p == 'V' )
                    {
                        encoder.append ( lastType, ( int ) ( * ( ++p ) ) - 1 );
                    }
                    else
                    {
                        encoder.append ( *p, 1 );
                        lastType = *p;
                    }
                }

                encoder.encodeCigarString ( openGapScore, extendGapScore );
                cigarString[lOr] = encoder.cigarString;
                // To get edit distance
                int L = batch->lengths[i] - encoder.charCount['I'] - encoder.charCount['S'];
                int numOfMismatch = ( L * matchScore + encoder.gapPenalty - batch->scores[lOr][i] ) /
                                    ( matchScore - mismatchScore );
                editdist[lOr] = encoder.charCount['I'] + encoder.charCount['D'] + numOfMismatch;
                DIS[lOr] = encoder.charCount['D'] - encoder.charCount['I'] - encoder.charCount['S'];
            }

            //#define MC_GetMateID(x) (((x)&1)?((x)-1):((x)+1))
            DeepDPAlignResult result;
            result.readID = batch->canInfos[i].readIDLeft - readSide;
            result.strand_1 = ( ( readSide == 0 ) ? engine->peStrandLeftLeg : engine->peStrandRightLeg );
            result.strand_2 = ( ( mateSide == 0 ) ? engine->peStrandLeftLeg : engine->peStrandRightLeg );
            result.algnmt_1 = batch->canInfos[i].pos[readSide] + batch->hitLocs[readSide][i];
            result.algnmt_2 = batch->canInfos[i].pos[mateSide] + batch->hitLocs[mateSide][i];
            result.score_1 = batch->scores[readSide][i];
            result.score_2 = batch->scores[mateSide][i];
            result.cigarString_1 = cigarString[readSide];
            result.cigarString_2 = cigarString[mateSide];
            result.editdist_1 = editdist[readSide];
            result.editdist_2 = editdist[mateSide];

            if ( result.algnmt_1 < result.algnmt_2 )
                result.insertSize = result.algnmt_2 - result.algnmt_1 +
                                    batch->lengths[i] + DIS[mateSide];
            else
                result.insertSize = result.algnmt_1 - result.algnmt_2 +
                                    batch->lengths[i] + DIS[readSide];

            result.num_sameScore_1 = batch->maxScoreCounts[readSide][i]; //TODO
            result.num_sameScore_2 = batch->maxScoreCounts[mateSide][i];
            resultBatch->push_back ( result );
        }
    }

    // Output
    engine->algnThreadContext[threadId].resultBatch = resultBatch;
    int * pid = &threadId;
    engine->outputThreadDelegator.schedule ( pid );
    sem_wait ( & ( engine->algnThreadContext[threadId].outputACKSem ) );
}

void DeepDP_Space::DP2OutputThread ( int threadId, int *& pCallThreadId )
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    int callThreadId = *pCallThreadId;
    int batchID = engine->algnThreadContext[callThreadId].batchID;
    DP2ResultBatch * resultBatch = engine->algnThreadContext[callThreadId].resultBatch;
    sem_post ( & ( engine->algnThreadContext[callThreadId].outputACKSem ) );
    vector<DP2ResultBatch *> & dpResult = engine->resultStream->dp2Result;

    while ( dpResult.size () <= batchID )
    {
        dpResult.push_back ( NULL );
    }

    dpResult[batchID] = resultBatch;
#define MC_DP2OutputRead() { \
        engine->outputBuf->ready(); \
        if (engine->outputBuf->size > 0) { \
            outputDeepDPResult2(engine->outputBuf->elements, engine->outputBuf->size, \
                                engine->queries, engine->upkdReadLengths, \
                                engine->origReadIDs, engine->upkdQueryNames, engine->upkdQualities, \
                                engine->inputMaxReadLength, engine->accumReadNum, engine->outputFormat, \
                                engine->outputFile, engine->samOutputDPFilePtr, engine->index, \
                                engine->peStrandLeftLeg, engine->peStrandRightLeg); \
            engine->dp2AlignedRead += 1; \
            engine->dp2Alignment += engine->outputBuf->size; \
            engine->alignFlags->set(engine->lastReadID << 1); \
        } \
    }
    uint numOut = engine->resultStream->numOut;

    while ( numOut < dpResult.size () && dpResult[numOut] != NULL )
    {
        //OUTPUT HERE
        DP2ResultBatch & batch = *dpResult[numOut];

        for ( int i = 0; i < batch.size (); i++ )
        {
            DeepDPAlignResult & result = batch[i];
            int pairID = result.readID >> 1;

            if ( pairID != engine->lastReadID )
            {
                MC_DP2OutputRead ();
                engine->outputBuf->clear ();
                engine->lastReadID = pairID;
            }

            engine->outputBuf->add ( result );
        }

        ++numOut;
    }

    engine->resultStream->numOut = numOut;
}

void DeepDP_Space::DP2OutputThreadFinalize ()
{
    PairEndAlignmentEngine * engine = PairEndAlignmentEngine::engine;
    // last read
    MC_DP2OutputRead ();
    engine->outputBuf->clear ();
}




// Temporary data for testing
// DPParameters Constants::deepDPPara_Len100;

Constants::Constants ()
{
    /*  deepDPPara_Len100.paramRead[0].cutoffThreshold = 30;
        deepDPPara_Len100.paramRead[0].maxHitNum = 100;
        deepDPPara_Len100.paramRead[0].seedLength = 26;
        deepDPPara_Len100.paramRead[0].sampleDist = 13;
    */
}

Constants constants;

