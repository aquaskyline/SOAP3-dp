/*
 *
 *    CPUfunctions.cpp
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

#include "CPUfunctions.h"

int getDefaultMismatchNum ( uint read_length )
{
    // for user does not specify # of mismatches and dp is disabled
    // if the read length < 50, then return DEFAULT_NUM_MISMATCH_FOR_SHORT_READ
    // else return DEFAULT_NUM_MISMATCH_FOR_NORMAL_READ
    if ( read_length > 50 )
    { return DEFAULT_NUM_MISMATCH_NO_DP_NORMAL_READ; }
    else
    { return DEFAULT_NUM_MISMATCH_NO_DP_SHORT_READ; }
}

int getMaxHitNumForDefaultDP ( uint read_length )
{
    // get the max hit # for default DP
    if ( read_length > 50 )
    { return MAX_SEED_HITS_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { return MAX_SEED_HITS_DEFAULT_DP_FOR_SHORT_READ; }
}

void getParameterForAllDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params )
{
    // get the parameters for All DP
    // not related to read length
    dp_params.matchScore = ini_params.Ini_MatchScore;
    dp_params.mismatchScore = ini_params.Ini_MismatchScore;
    dp_params.openGapScore = ini_params.Ini_GapOpenScore;
    dp_params.extendGapScore = ini_params.Ini_GapExtendScore;
    dp_params.numOfCPUThreads = ini_params.Ini_NumOfCpuThreads;
    dp_params.numOfCPUForSeeding = NUM_CPU_THREADS_FOR_SEEDING;
}


void getParameterForDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 )
{
    // get the parameters for Default DP
    // for paired-end reads
    if ( ini_params.Ini_isDefaultThreshold == 1 )
    {
        dp_params.paramRead[0].cutoffThreshold = ( int ) ceil ( DP_SCORE_THRESHOLD_RATIO * ( double ) read_length );
        dp_params.paramRead[1].cutoffThreshold = ( int ) ceil ( DP_SCORE_THRESHOLD_RATIO * ( double ) read_length2 );
    }
    else
    {
        dp_params.paramRead[0].cutoffThreshold = ini_params.Ini_DPScoreThreshold;
        dp_params.paramRead[1].cutoffThreshold = ini_params.Ini_DPScoreThreshold;
    }

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;

    if ( read_length > 50 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_SHORT_READ; }

    if ( read_length2 > 50 )
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_DEFAULT_DP_FOR_SHORT_READ; }

    // the following parameters are not needed for default DP
    // seedLength, sampleDist, tailTrimLen, singleDPSeedNum and singleDPSeedPos[]
}

void getParameterForNewDefaultDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 )
{
    // get the parameters for Default DP
    // for paired-end reads
    if ( ini_params.Ini_isDefaultThreshold == 1 )
    {
        dp_params.paramRead[0].cutoffThreshold = ( int ) ceil ( DP_SCORE_THRESHOLD_RATIO * ( double ) read_length );
        dp_params.paramRead[1].cutoffThreshold = ( int ) ceil ( DP_SCORE_THRESHOLD_RATIO * ( double ) read_length2 );
    }
    else
    {
        dp_params.paramRead[0].cutoffThreshold = ini_params.Ini_DPScoreThreshold;
        dp_params.paramRead[1].cutoffThreshold = ini_params.Ini_DPScoreThreshold;
    }

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;

    if ( read_length > 50 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_NEW_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_NEW_DEFAULT_DP_FOR_SHORT_READ; }

    if ( read_length2 > 50 )
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_NEW_DEFAULT_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_NEW_DEFAULT_DP_FOR_SHORT_READ; }

    if ( read_length > 75 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_LONG_READ; }
    else if ( read_length > 50 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_MEDIAN_READ; }
    else
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_SHORT_READ; }

    if ( read_length2 > 75 )
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_LONG_READ; }
    else if ( read_length2 > 50 )
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_MEDIAN_READ; }
    else
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_SHORT_READ; }
}


void getParameterForDeepDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length, uint read_length2 )
{
    // get the parameters for Deep DP
    if ( ini_params.Ini_isDefaultThreshold == 1 )
    {
        dp_params.paramRead[0].cutoffThreshold = ( int ) ceil ( DP_SCORE_THRESHOLD_RATIO * ( double ) read_length );
        dp_params.paramRead[1].cutoffThreshold = ( int ) ceil ( DP_SCORE_THRESHOLD_RATIO * ( double ) read_length2 );
    }
    else
    {
        dp_params.paramRead[0].cutoffThreshold = ini_params.Ini_DPScoreThreshold;
        dp_params.paramRead[1].cutoffThreshold = ini_params.Ini_DPScoreThreshold;
    }

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;

    if ( read_length > 50 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_SHORT_READ; }

    if ( read_length2 > 50 )
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ; }
    else
    { dp_params.paramRead[1].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_SHORT_READ; }

    if ( read_length > 150 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_V_LONG_READ; }
    else if ( read_length > 80 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_LONG_READ; }
    else if ( read_length > 60 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_MEDIAN_READ; }
    else if ( read_length > 40 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_SHORT_READ; }
    else
    { dp_params.paramRead[0].seedLength = SEED_LEN_DEEP_DP_FOR_V_SHORT_READ; }

    if ( read_length2 > 150 )
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_V_LONG_READ; }
    else if ( read_length2 > 80 )
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_LONG_READ; }
    else if ( read_length2 > 60 )
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_MEDIAN_READ; }
    else if ( read_length2 > 40 )
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_SHORT_READ; }
    else
    { dp_params.paramRead[1].seedLength = SEED_LEN_DEEP_DP_FOR_V_SHORT_READ; }

    dp_params.paramRead[0].sampleDist = dp_params.paramRead[0].seedLength * SEED_SAMPLE_RATE_DEEP_DP;
    dp_params.paramRead[1].sampleDist = dp_params.paramRead[1].seedLength * SEED_SAMPLE_RATE_DEEP_DP;
    dp_params.tailTrimLen = TAIL_TRIM_SEEDING_DEEP_DP;
    // the following parameters are not needed for deep DP
    // singleDPSeedNum and singleDPSeedPos[]
}

void getParameterForSingleDP ( DPParameters & dp_params, InputOptions & input_options, IniParams & ini_params, uint read_length )
{
    // get the parameters for Single-end DP
    if ( ini_params.Ini_isDefaultThreshold == 1 )
    {
        dp_params.paramRead[0].cutoffThreshold = ( int ) ceil ( DP_SCORE_THRESHOLD_RATIO * ( double ) read_length );
    }
    else
    {
        dp_params.paramRead[0].cutoffThreshold = ini_params.Ini_DPScoreThreshold;
    }

    dp_params.softClipLeft = ini_params.Ini_maxFrontLenClipped;
    dp_params.softClipRight = ini_params.Ini_maxEndLenClipped;

    if ( read_length > 300 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_V_LONG_READ; }
    else if ( read_length > 80 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_LONG_READ; }
    else if ( read_length > 60 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_MEDIAN_READ; }
    else if ( read_length > 40 )
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_SHORT_READ; }
    else
    { dp_params.paramRead[0].maxHitNum = MAX_SEED_HITS_SINGLE_DP_FOR_V_SHORT_READ; }

    if ( read_length > 300 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_SINGLE_DP_FOR_V_LONG_READ; }
    else if ( read_length > 80 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_SINGLE_DP_FOR_LONG_READ; }
    else if ( read_length > 60 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_SINGLE_DP_FOR_MEDIAN_READ; }
    else if ( read_length > 40 )
    { dp_params.paramRead[0].seedLength = SEED_LEN_SINGLE_DP_FOR_SHORT_READ; }
    else
    { dp_params.paramRead[0].seedLength = SEED_LEN_SINGLE_DP_FOR_V_SHORT_READ; }

    // update: for reads longer than 100, one more seed for every extra 100 bases
    if ( read_length > 100 )
    {
        dp_params.singleDPSeedNum = SEED_NUM_SINGLE_DP + read_length / 100;
    }
    else
    {
        dp_params.singleDPSeedNum = SEED_NUM_SINGLE_DP;
    }

    // the followings are obsoleted
    dp_params.singleDPSeedPos[0] = 0;
    int X;

    if ( read_length > 80 )
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_LONG_READ; }
    else if ( read_length > 60 )
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_MEDIAN_READ; }
    else if ( read_length > 40 )
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_SHORT_READ; }
    else
    { X = TAIL_TRIM_SEEDING_SINGLE_DP_FOR_V_SHORT_READ; }

    dp_params.singleDPSeedPos[2] = ( read_length - X ) * 0.5 - 1;

    if ( dp_params.singleDPSeedPos[2] > read_length - dp_params.paramRead[0].seedLength )
    { dp_params.singleDPSeedPos[2] = read_length - dp_params.paramRead[0].seedLength; }

    dp_params.singleDPSeedPos[1] = ( dp_params.singleDPSeedPos[0] + dp_params.singleDPSeedPos[2] ) / 2 - 1;
    dp_params.tailTrimLen = 0;
    // the following parameters are not needed for single DP
    // sampleDist and tailTrimLen
}



// pack the reads with no alignment together
// return number of reads with no alignment
uint packReads ( uint * queries, uint * readIDs, uint * readLengths, uint * seedLengths, unsigned char * noAlignment,
                 uint wordPerQuery, uint numQueries )
{
    ullint total_no_alignment = 0;

    for ( ullint readId = 0; readId < numQueries; ++readId )
    {
        if ( noAlignment[readId] == 1 )
        {
            // the read has no alignment
            if ( total_no_alignment < readId )
            {
                ullint srcQueryOffset = ( readId ) / 32 * 32 * wordPerQuery + readId % 32;
                ullint srcNewQueryOffset = ( total_no_alignment ) / 32 * 32 * wordPerQuery + total_no_alignment % 32;

                for ( ullint i = 0; i < wordPerQuery; ++i ) // copy each word of the read
                {
                    queries[srcNewQueryOffset + i * 32] = * ( queries + ( srcQueryOffset + i * 32 ) );
                }

                readLengths[total_no_alignment] = readLengths[readId];

                if ( seedLengths != NULL )
                { seedLengths[total_no_alignment] = seedLengths[readId]; }

                readIDs[total_no_alignment] = readIDs[readId];
            }

            total_no_alignment++;
        }
    }

    return total_no_alignment;
}

// pack the reads with no alignment together
// return readIDS of the unaligned reads
uint * packReads2 ( uint * queries, uint * readLengths, unsigned char * noAlignment,
                    uint wordPerQuery, uint numQueries, uint & numUnAligned )
{
    uint * readIDs = ( uint * ) malloc ( sizeof ( unsigned int ) * numQueries );
    ullint total_no_alignment = 0;

    for ( ullint readId = 0; readId < numQueries; ++readId )
    {
        if ( noAlignment[readId] == 1 )
        {
            // the read has no alignment
            if ( total_no_alignment < readId )
            {
                ullint srcQueryOffset = ( readId ) / 32 * 32 * wordPerQuery + readId % 32;
                ullint srcNewQueryOffset = ( total_no_alignment ) / 32 * 32 * wordPerQuery + total_no_alignment % 32;

                for ( ullint i = 0; i < wordPerQuery; ++i ) // copy each word of the read
                {
                    queries[srcNewQueryOffset + i * 32] = * ( queries + ( srcQueryOffset + i * 32 ) );
                }

                readLengths[total_no_alignment] = readLengths[readId];
            }

            readIDs[total_no_alignment] = readId + 1;
            total_no_alignment++;
            // printf("total_no_alignment = %u\n", total_no_alignment);
        }
    }

    numUnAligned = total_no_alignment;
    return readIDs;
}



// pack the reads which are unpaired
void packUnPairedReads ( uint * queries, uint * readIDs, uint * readLengths, uint * unAlignedPair,
                         uint wordPerQuery, uint numOfUnPaired, ullint maxBatchSize )
{
    ullint roundUp = ( numOfUnPaired + 31 ) / 32 * 32;
    ullint totalQueryLength = roundUp * wordPerQuery;
    uint * newqueries;
    uint * newreadLengths;
    uint * newreadIDs;

    if ( numOfUnPaired > 0 )
    {
        // copy the content to the new array
        newqueries = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) ); // a large array to store all queries
        newreadLengths = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
        newreadIDs = ( uint * ) malloc ( roundUp * sizeof ( uint ) );

        for ( ullint i = 0; i < numOfUnPaired; i++ )
        {
            ullint readId = unAlignedPair[i];
            ullint srcQueryOffset = readId / 32 * 32 * wordPerQuery + readId % 32;
            ullint srcNewQueryOffset = i / 32 * 32 * wordPerQuery + i % 32;

            for ( ullint j = 0; j < wordPerQuery; ++j ) // copy each word of the read
            {
                newqueries[srcNewQueryOffset + j * 32] = * ( queries + ( srcQueryOffset + j * 32 ) );
            }

            newreadLengths[i] = readLengths[readId];
            newreadIDs[i] = readIDs[readId];
        }
    }

    // clear the content of the old array
    ullint maxRoundUp = ( maxBatchSize + 31 ) / 32 * 32;
    ullint maxTotalQueryLength = maxRoundUp * wordPerQuery;
    memset ( queries, 0, maxTotalQueryLength * sizeof ( uint ) );
    memset ( readLengths, 0, maxRoundUp * sizeof ( uint ) );
    memset ( readIDs, 0, maxRoundUp * sizeof ( uint ) );

    if ( numOfUnPaired > 0 )
    {
        // copy the content of the new array to the old array
        memcpy ( queries, newqueries, totalQueryLength * sizeof ( uint ) );
        memcpy ( readLengths, newreadLengths, roundUp * sizeof ( uint ) );
        memcpy ( readIDs, newreadIDs, roundUp * sizeof ( uint ) );
        // free the memory
        free ( newqueries );
        free ( newreadLengths );
        free ( newreadIDs );
    }
}



// repack the reads
// no read will be removed, but
// the reads which need to be processed in next-round by soap3 will be duplicated
// to the front of the list. The readIDs are stored inside the array called "needProcessPair"
// the corresponding readIDs inside "readInputForDP", "readInputForNewDP" and
// "bothUnalignedPairs" need to be updated correspondingly.
void repackUnPairedReads ( uint ** orgqueries, uint ** orgreadIDs, uint ** orgreadLengths, uint * needProcessPair,
                           uint wordPerQuery, uint numOfReadsToProcess, ullint numOfTotalReads,
                           ReadInputForDPArrays * readInputForDPall,
                           ReadInputForDPArrays * readInputForNewDPall,
                           BothUnalignedPairsArrays * bothUnalignedPairsArrays )
{
    uint * queries = ( *orgqueries );
    uint * readIDs = ( *orgreadIDs );
    uint * readLengths = ( *orgreadLengths );

    if ( numOfReadsToProcess > 0 )
    {
        ullint i, j;
        ullint roundUp = ( numOfReadsToProcess + numOfTotalReads + 31 ) / 32 * 32;
        ullint totalQueryLength = roundUp * wordPerQuery;
        ullint readId, srcQueryOffset, srcNewQueryOffset, newReadId;
        // copy the content to the new array
        uint * newqueries = ( uint * ) malloc ( totalQueryLength * sizeof ( uint ) );
        uint * newreadLengths = ( uint * ) malloc ( roundUp * sizeof ( uint ) );
        uint * newreadIDs = ( uint * ) malloc ( roundUp * sizeof ( uint ) );

        for ( i = 0; i < numOfReadsToProcess; i++ )
        {
            readId = needProcessPair[i];
            srcQueryOffset = readId / 32 * 32 * wordPerQuery + readId % 32;
            srcNewQueryOffset = i / 32 * 32 * wordPerQuery + i % 32;

            for ( j = 0; j < wordPerQuery; ++j ) // copy each word of the read
            {
                newqueries[srcNewQueryOffset + j * 32] = * ( queries + ( srcQueryOffset + j * 32 ) );
            }

            newreadLengths[i] = readLengths[readId];
            newreadIDs[i] = readIDs[readId];
        }

        // append the content of the old array to the new array
        newReadId = numOfReadsToProcess;

        for ( readId = 0; readId < numOfTotalReads; readId++ )
        {
            srcQueryOffset = readId / 32 * 32 * wordPerQuery + readId % 32;
            srcNewQueryOffset = newReadId / 32 * 32 * wordPerQuery + newReadId % 32;

            for ( j = 0; j < wordPerQuery; ++j ) // copy each word of the read
            {
                newqueries[srcNewQueryOffset + j * 32] = * ( queries + ( srcQueryOffset + j * 32 ) );
            }

            newreadLengths[newReadId] = readLengths[readId];
            newreadIDs[newReadId] = readIDs[readId];
            newReadId++;
        }

        // update the read IDs in the readInputForDPall
        for ( i = 0; i < readInputForDPall->numArrays; i++ )
        {
            ReadInputForDP * currArray = readInputForDPall->inputArrays[i];

            // sa list
            for ( j = 0; j < currArray->saRangeTotalNum; j++ )
            {
                ( currArray->sa_list[j] ).readID += numOfReadsToProcess;
            }

            // occ list
            for ( j = 0; j < currArray->occTotalNum; j++ )
            {
                ( currArray->occ_list[j] ).readID += numOfReadsToProcess;
            }
        }

        // update the read IDs in the readInputForNewDPall
        for ( i = 0; i < readInputForNewDPall->numArrays; i++ )
        {
            ReadInputForDP * currArray = readInputForNewDPall->inputArrays[i];

            // sa list
            for ( j = 0; j < currArray->saRangeTotalNum; j++ )
            {
                ( currArray->sa_list[j] ).readID += numOfReadsToProcess;
            }

            // occ list
            for ( j = 0; j < currArray->occTotalNum; j++ )
            {
                ( currArray->occ_list[j] ).readID += numOfReadsToProcess;
            }
        }

        // update the read IDs in the bothUnalignedPairsArrays
        for ( i = 0; i < bothUnalignedPairsArrays->arrayNum; i++ )
        {
            BothUnalignedPairs * currArray = bothUnalignedPairsArrays->array[i];

            for ( j = 0; j < currArray->totalNum; j++ )
            {
                ( currArray->readIDs ) [j] += numOfReadsToProcess;
            }
        }

        // free the memory of the original array
        free ( queries );
        free ( readLengths );
        free ( readIDs );
        ( *orgqueries ) = newqueries;
        ( *orgreadLengths ) = newreadLengths;
        ( *orgreadIDs ) = newreadIDs;
    }
}

unsigned long long ProcessReadSingleStrand2 ( SRAQueryInput * qInput, SRAModel * aModel, int whichCase, SAList * sa_list, OCCList * occ_list )
{
    unsigned long long saCount = 0;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    BWT * bwt = aIndex->bwt;
    SRACase * thisCase = & ( aModel->cases[whichCase] );
    unsigned long long saRanges[4];
    rOutput->IsOpened = 0;
    rOutput->IsClosed = 0;
    rOutput->IsUnique = 0;
    rOutput->IsBest = 0;
    saRanges[0] = 0;
    saRanges[1] = bwt->textLength;
    saRanges[2] = 0;
    saRanges[3] = bwt->textLength;

    if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP )
    {
        saCount += BWTExactModelBackward_Lookup2 ( qInput, thisCase, 0, saRanges, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP )
    {
        saCount += BWTExactModelForward_Lookup2 ( qInput, thisCase, 0, saRanges, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
    {
        saCount += BWTModelSwitchBackward2 ( qInput, 0, 0,
                                             thisCase, 0,
                                             saRanges, 0, 0, 0, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT )
    {
        saCount += BWTModelSwitchAnyDirection2 ( qInput, 0, 0,
                   thisCase, 0,
                   saRanges, 0, 0, 0, sa_list, occ_list );
    }

    return saCount;
}


unsigned long long ProcessReadDoubleStrand2 ( SRAQueryInput * qInput, SRAQueryInput * qInput_neg, SRAModel * aModel, SRAModel * aModel_neg, int whichCase, SAList * sa_list, OCCList * occ_list )
{
    unsigned long long saCount = 0;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    BWT * bwt = aIndex->bwt;
    SRACase * thisCase =  & ( aModel->cases[whichCase] );
    unsigned long long saRanges[4];
    rOutput->IsOpened = 0;
    rOutput->IsClosed = 0;
    rOutput->IsUnique = 0;
    rOutput->IsBest = 0;
    saRanges[0] = 0;
    saRanges[1] = bwt->textLength;
    saRanges[2] = 0;
    saRanges[3] = bwt->textLength;

    if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP )
    {
        saCount += BWTExactModelBackward_Lookup2 ( qInput, thisCase, 0, saRanges, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP )
    {
        saCount += BWTExactModelForward_Lookup2 ( qInput, thisCase, 0, saRanges, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
    {
        saCount += BWTModelSwitchBackward2 ( qInput, 0, 0,
                                             thisCase, 0,
                                             saRanges, 0, 0, 0, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT )
    {
        saCount += BWTModelSwitchAnyDirection2 ( qInput, 0, 0,
                   thisCase, 0,
                   saRanges, 0, 0, 0, sa_list, occ_list );
    }

    saRanges[0] = 0;
    saRanges[1] = bwt->textLength;
    saRanges[2] = 0;
    saRanges[3] = bwt->textLength;
    thisCase =  & ( aModel_neg->cases[whichCase] );

    if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP )
    {
        saCount += BWTExactModelBackward_Lookup2 ( qInput_neg, thisCase, 0, saRanges, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP )
    {
        saCount += BWTExactModelForward_Lookup2 ( qInput_neg, thisCase, 0, saRanges, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
    {
        saCount += BWTModelSwitchBackward2 ( qInput_neg, 0, 0,
                                             thisCase, 0,
                                             saRanges, 0, 0, 0, sa_list, occ_list );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT )
    {
        saCount += BWTModelSwitchAnyDirection2 ( qInput_neg, 0, 0,
                   thisCase, 0,
                   saRanges, 0, 0, 0, sa_list, occ_list );
    }

    return saCount;
}

unsigned long long ProcessOneMoreMismatchAllCases ( unsigned char * charMap, SRAQueryInput * qInput_Positive, SRAQueryInput * qInput_Negative,
        SRAModel * aModel, SRAModel * aModel_neg, SAList * sa_list, OCCList * occ_list, int totalNumCases )
{
    unsigned long long saCount = 0;
    int i;
    SRAQueryInfo * qInfo_Positive = qInput_Positive->QueryInfo;
    SRAQueryInfo * qInfo_Negative = qInput_Negative->QueryInfo;
    unsigned int readLength = qInfo_Positive->ReadLength;
    SRASetting * qSetting = qInput_Positive->QuerySetting;
    SRAQueryResultCount * rOutput = qInput_Positive->QueryOutput;
    rOutput->TotalOccurrences = 0;
    memset ( rOutput->WithError, 0, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );
    qSetting->MaxError++;
#ifndef BGS_DISABLE_NEGATIVE_STRAND
    unsigned char oStrandQuery[MAX_READ_LENGTH];

    for ( i = 0; i < readLength; i++ )
    {
        oStrandQuery[i] = soap3DnaComplement[qInfo_Positive->ReadCode[readLength - i - 1]];
    }

    qInfo_Negative->ReadCode = oStrandQuery;
    qInfo_Positive->ReadStrand = QUERY_POS_STRAND;
    qInfo_Negative->ReadStrand = QUERY_NEG_STRAND;
    qInfo_Negative->ReadLength = readLength;
    qInput_Negative->QuerySetting = qInput_Positive->QuerySetting;
    //Reversing the read from positive strand to negative strand.
#endif

    for ( int whichCase = 0; whichCase < totalNumCases; whichCase++ )
    {
#ifdef BGS_DISABLE_NEGATIVE_STRAND
        saCount += ProcessReadSingleStrand2 ( qInput_Positive, aModel, whichCase, sa_list, occ_list );
#else
        saCount += ProcessReadDoubleStrand2 ( qInput_Positive, qInput_Negative, aModel, aModel_neg, whichCase, sa_list, occ_list );
#endif
    }

    qSetting->MaxError--;
    return saCount;
}



unsigned long long ProcessReadSingleStrand3 ( SRAQueryInput * qInput, SRAModel * aModel, int whichCase, SingleAlgnResult * alignResult )
{
    unsigned long long saCount = 0;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    BWT * bwt = aIndex->bwt;
    SRACase * thisCase = & ( aModel->cases[whichCase] );
    unsigned long long saRanges[4];
    rOutput->IsOpened = 0;
    rOutput->IsClosed = 0;
    rOutput->IsUnique = 0;
    rOutput->IsBest = 0;
    saRanges[0] = 0;
    saRanges[1] = bwt->textLength;
    saRanges[2] = 0;
    saRanges[3] = bwt->textLength;

    if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP )
    {
        saCount += BWTExactModelBackward_Lookup3 ( qInput, thisCase, 0, saRanges, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP )
    {
        saCount += BWTExactModelForward_Lookup3 ( qInput, thisCase, 0, saRanges, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
    {
        saCount += BWTModelSwitchBackward3 ( qInput, 0, 0,
                                             thisCase, 0,
                                             saRanges, 0, 0, 0, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT )
    {
        saCount += BWTModelSwitchAnyDirection3 ( qInput, 0, 0,
                   thisCase, 0,
                   saRanges, 0, 0, 0, alignResult );
    }

    return saCount;
}



unsigned long long ProcessReadDoubleStrand3 ( SRAQueryInput * qInput, SRAQueryInput * qInput_neg, SRAModel * aModel, SRAModel * aModel_neg, int whichCase, SingleAlgnResult * alignResult )
{
    unsigned long long saCount = 0;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    BWT * bwt = aIndex->bwt;
    SRACase * thisCase =  & ( aModel->cases[whichCase] );
    unsigned long long saRanges[4];
    rOutput->IsOpened = 0;
    rOutput->IsClosed = 0;
    rOutput->IsUnique = 0;
    rOutput->IsBest = 0;
    saRanges[0] = 0;
    saRanges[1] = bwt->textLength;
    saRanges[2] = 0;
    saRanges[3] = bwt->textLength;

    if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP )
    {
        saCount += BWTExactModelBackward_Lookup3 ( qInput, thisCase, 0, saRanges, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP )
    {
        saCount += BWTExactModelForward_Lookup3 ( qInput, thisCase, 0, saRanges, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
    {
        saCount += BWTModelSwitchBackward3 ( qInput, 0, 0,
                                             thisCase, 0,
                                             saRanges, 0, 0, 0, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT )
    {
        saCount += BWTModelSwitchAnyDirection3 ( qInput, 0, 0,
                   thisCase, 0,
                   saRanges, 0, 0, 0, alignResult );
    }

    saRanges[0] = 0;
    saRanges[1] = bwt->textLength;
    saRanges[2] = 0;
    saRanges[3] = bwt->textLength;
    thisCase =  & ( aModel_neg->cases[whichCase] );

    if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP )
    {
        saCount += BWTExactModelBackward_Lookup3 ( qInput_neg, thisCase, 0, saRanges, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP )
    {
        saCount += BWTExactModelForward_Lookup3 ( qInput_neg, thisCase, 0, saRanges, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
    {
        saCount += BWTModelSwitchBackward3 ( qInput_neg, 0, 0,
                                             thisCase, 0,
                                             saRanges, 0, 0, 0, alignResult );
    }
    else if ( thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT )
    {
        saCount += BWTModelSwitchAnyDirection3 ( qInput_neg, 0, 0,
                   thisCase, 0,
                   saRanges, 0, 0, 0, alignResult );
    }

    return saCount;
}


// set the multi-threading arguments
void setHostKernelArguments ( HostKernelArguements * hostKernelArguments, pthread_t * threads, IniParams ini_params, Soap3Index * index,
                              uint maxReadLength, uint word_per_query, uint word_per_ans,
                              InputOptions * input_options )
{
    char readType = input_options->readType;
    int insert_high = input_options->insert_high;
    int insert_low = input_options->insert_low;
    int alignmentType = input_options->alignmentType; // for CPU Alignment
    int reportType = input_options->alignmentType; // for reporting the alignments
    int outputFormat = input_options->outputFormat;

    // Convert to CPU Kernel data
    // This shouldn't have been needed...
    if ( readType == SINGLE_READ )
    {
        // for single-end alignment
        if ( alignmentType == OUTPUT_ALL_BEST )
        {
            alignmentType = SRA_REPORT_ALL_BEST;
        }
        else if ( alignmentType == OUTPUT_UNIQUE_BEST )
        {
            alignmentType = SRA_REPORT_UNIQUE_BEST;
        }
        else if ( alignmentType == OUTPUT_RANDOM_BEST )
        {
            alignmentType = SRA_REPORT_RANDOM_BEST;
        }
        else
        {
            alignmentType = SRA_REPORT_ALL;
        }
    }
    else
    {
        // for paired-end alignment, the alignment type always is REPORT_ALL
        alignmentType = SRA_REPORT_ALL;
    }

    // MULTI-THREADING arguments
    for ( uint threadId = 0; threadId < ini_params.Ini_NumOfCpuThreads; ++threadId )
    {
        threads[threadId] = 0;
        hostKernelArguments[threadId].occ = OCCConstruct ();

        for ( uint i = 0 ; i <= MAX_NUM_OF_MISMATCH ; i++ )
        {
            hostKernelArguments[threadId].sraAccumulatedResultCount.WithError[i] = 0;
        }

        // Copy SRAIndex
        memcpy ( & ( hostKernelArguments[threadId].sraIndex ), index->sraIndex, sizeof ( SRAIndex ) );

        hostKernelArguments[threadId].sraQuerySettings.occ = hostKernelArguments[threadId].occ;
        hostKernelArguments[threadId].sraQuerySettings.OutFilePtr = NULL;
        hostKernelArguments[threadId].sraQuerySettings.OutFileName = NULL;
        hostKernelArguments[threadId].sraQuerySettings.OutFilePart = 0;
        hostKernelArguments[threadId].sraQuerySettings.OutFileFormat = outputFormat;

        if ( readType == SINGLE_READ )
        { hostKernelArguments[threadId].sraQuerySettings.MaxOutputPerRead = ini_params.Ini_MaxOutputPerRead; }
        else
        { hostKernelArguments[threadId].sraQuerySettings.MaxOutputPerRead = ini_params.Ini_MaxHitsEachEndForPairing; }

        hostKernelArguments[threadId].sraQuerySettings.writtenCacheCount = 0;
        hostKernelArguments[threadId].sraQuerySettings.ReadStrand = QUERY_BOTH_STRAND;
        hostKernelArguments[threadId].sraQuerySettings.ErrorType = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        hostKernelArguments[threadId].sraQuerySettings.OutputType = alignmentType;
        hostKernelArguments[threadId].sraQuerySettings.SAMOutFilePtr = NULL;

        hostKernelArguments[threadId].sraAccumulatedResultCount.TotalOccurrences = 0;
        hostKernelArguments[threadId].sraAccumulatedResultCount.RetrievedBySa = 0;
        hostKernelArguments[threadId].sraAccumulatedResultCount.RetrievedByCe = 0;
        hostKernelArguments[threadId].sraAccumulatedResultCount.RetrievedByHocc = 0;
        hostKernelArguments[threadId].batchFirstReadId = 0;
        hostKernelArguments[threadId].skipFirst = 0;
        hostKernelArguments[threadId].numQueries = 0;
        hostKernelArguments[threadId].queries = NULL;
        hostKernelArguments[threadId].word_per_query = word_per_query;
        hostKernelArguments[threadId].answers = NULL;
        hostKernelArguments[threadId].badReadIndices = NULL;
        hostKernelArguments[threadId].badAnswers = NULL;
        hostKernelArguments[threadId].badStartOffset = NULL;
        hostKernelArguments[threadId].badCountOffset = NULL;
        hostKernelArguments[threadId].alignedOcc = 0;
        hostKernelArguments[threadId].alignedReads = 0;
        hostKernelArguments[threadId].threadId = threadId;
        hostKernelArguments[threadId].maxReadLength = maxReadLength;
        hostKernelArguments[threadId].word_per_ans = word_per_ans;
        hostKernelArguments[threadId].readType = readType;
        hostKernelArguments[threadId].reportType = reportType;
        hostKernelArguments[threadId].insert_low = insert_low;
        hostKernelArguments[threadId].insert_high = insert_high;
        hostKernelArguments[threadId].peStrandLeftLeg = ini_params.Ini_PEStrandLeftLeg;
        hostKernelArguments[threadId].peStrandRightLeg = ini_params.Ini_PEStrandRightLeg;
        hostKernelArguments[threadId].peMaxOutputPerPair = ini_params.Ini_PEMaxOutputPerPair;
        hostKernelArguments[threadId].maxHitNum = input_options->maxHitNum;
        hostKernelArguments[threadId].maxHitNum2 = input_options->maxHitNum2;
        hostKernelArguments[threadId].enableDP = input_options->enableDP;
    }
}


// obtain the number of cases for this number of mismatch
// and obtain the number of SA ranges allowed
void getParametersForThisMismatch ( uint numMismatch, uint & numCases, uint & sa_range_allowed_1,
                                    uint & sa_range_allowed_2, char & skip_round_2,
                                    uint & word_per_ans, uint & word_per_ans_2 )
{
    switch ( numMismatch )
    {
        case 0:
            numCases = NUM_CASES_0M;
            sa_range_allowed_1 = MAX_SA_RANGES_ALLOWED1_0;
            sa_range_allowed_2 = MAX_SA_RANGES_ALLOWED2_0;
            skip_round_2 = SKIP_ROUND2_0;
            break;

        case 1:
            numCases = NUM_CASES_1M;
            sa_range_allowed_1 = MAX_SA_RANGES_ALLOWED1_1;
            sa_range_allowed_2 = MAX_SA_RANGES_ALLOWED2_1;
            skip_round_2 = SKIP_ROUND2_1;
            break;

        case 2:
            numCases = NUM_CASES_2M;
            sa_range_allowed_1 = MAX_SA_RANGES_ALLOWED1_2;
            sa_range_allowed_2 = MAX_SA_RANGES_ALLOWED2_2;
            skip_round_2 = SKIP_ROUND2_2;
            break;

        case 3:
            numCases = NUM_CASES_3M;
            sa_range_allowed_1 = MAX_SA_RANGES_ALLOWED1_3;
            sa_range_allowed_2 = MAX_SA_RANGES_ALLOWED2_3;
            skip_round_2 = SKIP_ROUND2_3;
            break;

        case 4:
            numCases = NUM_CASES_4M;
            sa_range_allowed_1 = MAX_SA_RANGES_ALLOWED1_4;
            sa_range_allowed_2 = MAX_SA_RANGES_ALLOWED2_4;
            skip_round_2 = SKIP_ROUND2_4;
            break;

        default:
            numCases = MAX_NUM_CASES;
            sa_range_allowed_1 = 4;
            sa_range_allowed_2 = 1024;
            skip_round_2 = 0;
            break;
    }

    word_per_ans = sa_range_allowed_1 * 2;
    word_per_ans_2 = sa_range_allowed_2 * 2;
}

// obtain the number of cases for the number of mismatch
int GetNumCases ( uint numMismatch )
{
    uint numCases = MAX_NUM_CASES;

    switch ( numMismatch )
    {
        case 0:
            numCases = NUM_CASES_0M;
            break;

        case 1:
            numCases = NUM_CASES_1M;
            break;

        case 2:
            numCases = NUM_CASES_2M;
            break;

        case 3:
            numCases = NUM_CASES_3M;
            break;

        case 4:
            numCases = NUM_CASES_4M;
            break;
    }

    return numCases;
}


// obtain the number of cases for this number of mismatch
// and obtain the number of SA ranges allowed
// this is designed for seed alignments
void getParametersForThisMismatch2 ( uint numMismatch, uint & numCases, uint & sa_range_allowed_1,
                                     uint & sa_range_allowed_2, char & skip_round_2,
                                     uint & word_per_ans, uint & word_per_ans_2 )
{
    sa_range_allowed_1 = MAX_SA_RANGES_ALLOWED1_SMALL_READS;
    sa_range_allowed_2 = MAX_SA_RANGES_ALLOWED2_SMALL_READS;
    skip_round_2 = 0;

    switch ( numMismatch )
    {
        case 0:
            numCases = NUM_CASES_0M;
            sa_range_allowed_1 = MAX_SA_RANGES_ALLOWED1_0;
            skip_round_2 = 1;
            break;

        case 1:
            numCases = NUM_CASES_1M;
            break;

        case 2:
            numCases = NUM_CASES_2M;
            break;

        case 3:
            numCases = NUM_CASES_3M;
            break;

        case 4:
            numCases = NUM_CASES_4M;
            break;

        default:
            numCases = MAX_NUM_CASES;
            break;
    }

    word_per_ans = sa_range_allowed_1 * 2;
    word_per_ans_2 = sa_range_allowed_2 * 2;
}



// show the merging file command at the end of the program
void show_merge_file_command ( InputOptions input_options, char * queryFileName, char * queryFileName2 )
{
    printf ( "\n" );
    printf ( "Please type the following command to merge the alignment results into one file:\n" );

    if ( input_options.outputFormat == SRA_OUTPUT_FORMAT_DEFAULT )
    {
        // binary output format
        if ( input_options.readType == SINGLE_READ )
        {
            // single reads
            printf ( "./merge-binary.sh single %s\n", queryFileName );
        }
        else
        {
            // paired-end reads
            printf ( "./merge-binary.sh pair %s %s\n", queryFileName, queryFileName2 );
        }
    }
    else if ( input_options.outputFormat == SRA_OUTPUT_FORMAT_PLAIN )
    {
        // plain output format
        printf ( "./merge-succinct.sh %s\n", queryFileName );
    }
    else
    {
        // SAM output format
        printf ( "./merge-sam.sh %s\n", queryFileName );
    }
}

// show the merging file command at the end of the program
void show_merge_file_command2 ( InputOptions input_options, char * outputFileName, int cpuNumThreads )
{
    printf ( "\n" );
    printf ( "If necessary, " );

    if ( input_options.isOutputBinary == 1 )
    {
        // for merging the BAM output files
        printf ( "you may use 'samtools' to concatenate all the alignment results into one BAM file:\n" );
        printf ( "./samtools cat -o [output].bam" );
        int i;

        for ( i = 0; i < cpuNumThreads; i++ )
        {
            printf ( " %s.gout.%i", outputFileName, i + 1 );
        }

        printf ( " %s.dpout.1\n", outputFileName );
    }
    else
    {
        printf ( "please type the following command to merge the alignment results into one file:\n" );

        if ( input_options.outputFormat == SRA_OUTPUT_FORMAT_PLAIN )
        {
            // plain output format
            printf ( "./merge-succinct.sh %s\n", outputFileName );
        }
        else
        {
            // SAM output format
            printf ( "./merge-sam.sh %s\n", outputFileName );
        }
    }

    printf ( "\n" );
}



#define lastAlignedBoundary(offset, alignBoundary)              ( (offset) & (- (alignBoundary)) )
void retrieve2BWTQueryFormat ( unsigned int * packedPatterns,
                               unsigned int queryIdx, unsigned int wordPerQuery, unsigned char * unpackedPattern )
{
    unsigned int i;
    int j;
    int k = 0;
    unsigned int warpSkipSize = lastAlignedBoundary ( queryIdx, BGS_DEVICE_WARP_SIZE ) * wordPerQuery;
    packedPatterns += ( warpSkipSize + queryIdx % BGS_DEVICE_WARP_SIZE );

    for ( i = 0; i < wordPerQuery; i++ )
    {
        unsigned int currentWord = * ( packedPatterns );

        for ( j = 0; j < CHAR_PER_WORD; j++ )
        {
            char bit = currentWord & CHAR_MASK;
            currentWord >>= BIT_PER_CHAR;
            unpackedPattern[k++] = bit;
        }

        packedPatterns += BGS_DEVICE_WARP_SIZE;
    }
}

void keepSoap3Results ( unsigned int readID, SAList * sa_list, OCCList * occ_list,
                        ReadInputForDP * dpInput, HSPAux * hspaux )
{
    hspaux->sa_start[readID] = dpInput->saRangeTotalNum;
    hspaux->sa_num[readID] = sa_list->curr_size;
    hspaux->occ_start[readID] = dpInput->occTotalNum;
    hspaux->occ_num[readID] = occ_list->curr_size;
    ( ( ReadInputForDP ** ) hspaux->soap3AnsArray ) [readID] = dpInput;
}



void validateAlignments ( OCCList * occ_list, unsigned char * thisQuery,
                          unsigned int seed_len, unsigned int read_len, HSP * hsp, bool only_keep_best_ans,
                          int min_seed_mismatch, int max_mismatch, int max_hit_num, int threadID )
{
    // fprintf(stderr, "start validateAlignments\n");
    // fprintf(stderr, "occ_list->curr_size = %u\n", occ_list->curr_size);
    if ( read_len <= seed_len )
    { return; }

    unsigned int extension_len = read_len - seed_len;
    unsigned int packedQuery[ ( MAX_READ_LENGTH + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD];
    unsigned int revPackedQuery[ ( MAX_READ_LENGTH + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD];
    unsigned int targetPackedDNA[ ( MAX_READ_LENGTH + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD];
    thisQuery += seed_len;
    // make the packed DNA for the query
    createQueryPackedDNA ( thisQuery, extension_len, packedQuery );
    // make the packed DNA for the reverse strand of the query
    createRevQueryPackedDNA ( thisQuery, extension_len, revPackedQuery );
    int preTotMismatch = max_mismatch;
    unsigned int i;

    if ( occ_list->curr_size > 0 )
    {
        int new_occ_size = 0;

        for ( i = 0; i < occ_list->curr_size; i++ )
        {
            if ( ( occ_list->occ[i] ).mismatchCount < min_seed_mismatch )
            { continue; }

            int mismatch = max_mismatch + 1;

            // fprintf(stderr, "(occ_list->occ[i]).strand = %i\n", (occ_list->occ[i]).strand);
            // fprintf(stderr, "(occ_list->occ[i]).mismatchCount = %u\n", (occ_list->occ[i]).mismatchCount);
            // fprintf(stderr, "deduced ambpos = %u \n", ((occ_list->occ[i]).strand==1)?(occ_list->occ[i]).ambPosition:(occ_list->occ[i]).ambPosition-extension_len);
            if ( ( occ_list->occ[i] ).strand == 1 && ( occ_list->occ[i] ).ambPosition + read_len <= hsp->dnaLength )
            {
                // create the packedDNA for the target sequence;
                createTargetPackedDNA ( hsp->packedDNA, ( occ_list->occ[i] ).ambPosition + seed_len, extension_len, targetPackedDNA );
                mismatch = numMismatchNew ( packedQuery, targetPackedDNA, extension_len, threadID );
            }

            if ( ( occ_list->occ[i] ).strand == 2 && ( occ_list->occ[i] ).ambPosition >= extension_len )
            {
                // create the packedDNA for the target sequence;
                createTargetPackedDNA ( hsp->packedDNA, ( occ_list->occ[i] ).ambPosition - extension_len, extension_len, targetPackedDNA );
                mismatch = numMismatchNew ( revPackedQuery, targetPackedDNA, extension_len, threadID );
            }

            // fprintf(stderr, "mismatch = %i\n", (int) mismatch);
            int totMismatch = mismatch + ( occ_list->occ[i] ).mismatchCount;

            if ( totMismatch <= max_mismatch )
            {
                // the alignment is valid
                if ( only_keep_best_ans && totMismatch > preTotMismatch )
                { continue; }

                if ( only_keep_best_ans && totMismatch < preTotMismatch )
                {
                    new_occ_size = 0; // reset to zero, remove all previous results
                    preTotMismatch = totMismatch;
                }

                if ( new_occ_size < i )
                {
                    // move the alignment to the next valid entry
                    ( occ_list->occ[new_occ_size] ).readID = ( occ_list->occ[i] ).readID;
                    ( occ_list->occ[new_occ_size] ).strand = ( occ_list->occ[i] ).strand;
                    ( occ_list->occ[new_occ_size] ).mismatchCount = totMismatch;
                }
                else
                {
                    ( occ_list->occ[new_occ_size] ).mismatchCount = totMismatch;
                }

                // update the position
                unsigned int newAmbPos = ( ( occ_list->occ[i] ).strand == 1 ) ? ( occ_list->occ[i] ).ambPosition : ( occ_list->occ[i] ).ambPosition - extension_len;
                ( occ_list->occ[new_occ_size] ).ambPosition = newAmbPos;
                new_occ_size++;

                // stop when there are too many answers
                if ( preTotMismatch == min_seed_mismatch && new_occ_size >= max_hit_num )
                {
                    break;
                }
            }
        }

        occ_list->curr_size = new_occ_size;
    }

    // fprintf(stderr, "end validateAlignments\n");
}



void collect_all_answers ( SAList * sa_list, OCCList * occ_list, unsigned int ** answers, ullint batchReadId,
                           unsigned int word_per_ans, SRAQueryInput * qInput_Positive, SRAQueryInput * qInput_Negative,
                           unsigned int sa_range_allowed_1, unsigned int sa_range_allowed_2, char skip_round_2,
                           unsigned int numCases, int readLength, unsigned char * thisQuery, SRAModel ** SRAMismatchModel, SRAModel ** SRAMismatchModel_neg,
                           unsigned char * charMap,
                           unsigned int ** badReadIndices, unsigned int ** badAnswers,
                           unsigned int * badStartOffset, unsigned int * badCountOffset,
                           unsigned int * badReadPos, bool only_get_round1_ans, bool & isMoreThanSA1 )
{
    // fprintf(stderr, "readLength = %i\n", readLength);
    // collect all the answers
    // All the answers will be stored in sa_list and occ_list
    // output: isMoreThanSA1 - whether the number of SA ranges > SA_range_allowed_1
    isMoreThanSA1 = false;
    SRAQueryInfo * qInfo_Positive = qInput_Positive->QueryInfo;
    SRAQueryInfo * qInfo_Negative = qInput_Negative->QueryInfo;
    SRASetting * qSetting = qInput_Positive->QuerySetting;
    SRAIndex * aIndex = qInput_Positive->AlgnmtIndex;
    SRAQueryResultCount * rOutput = qInput_Positive->QueryOutput;
    rOutput->TotalOccurrences = 0;
    memset ( rOutput->WithError, 0, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );
    BWT * bwt = aIndex->bwt;
    uint * ans;
    unsigned int l, r;
    int i;
    int strand;
    char num_mis;
    int word_per_ans_2 = sa_range_allowed_2 * 2;
#ifndef BGS_DISABLE_NEGATIVE_STRAND
    unsigned char oStrandQuery[MAX_READ_LENGTH];
#endif

    for ( int whichCase = 0; whichCase < numCases; whichCase++ )
    {
        ans = answers[whichCase] + ( ( batchReadId ) / 32 * 32 * word_per_ans + ( batchReadId ) % 32 );

        if ( * ( ans ) > 0xFFFFFFFD )
        {
            // number of SA ranges > SA_range_allowed_1
            isMoreThanSA1 = true;
        }

        if ( * ( ans ) < 0xFFFFFFFD || only_get_round1_ans )
        {
            // Number of SA ranges reported within MAX_SA_RANGES_ALLOWED
            int start = ( * ( ans ) > 0xFFFFFFFD ) ? 1 : 0;

            for ( i = start; i < sa_range_allowed_1; i++ )
            {
                if ( ( * ( ans + ( i * 2 ) * 32 ) ) >= 0xFFFFFFFD || ( * ( ans + ( i * 2 + 1 ) * 32 ) ) >= 0xFFFFFFFD )
                {
                    break;
                }

                l = * ( ans + ( i * 2 ) * 32 );
                r = l + ( * ( ans + ( i * 2 + 1 ) * 32 ) &BGS_GPU_ANSWER_OFFSET_MASK );
                strand = ( ( ( * ( ans + ( i * 2 + 1 ) * 32 ) >> ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) & 1 ) ) + 1;
                num_mis = ( ( ( * ( ans + ( i * 2 + 1 ) * 32 ) >> ( BGS_GPU_ANSWER_OFFSET_LENGTH ) ) & 7 ) );

                if ( l <= r && r <= bwt->textLength ) // assume unique hits
                {
                    if ( rOutput->TotalOccurrences < qSetting->MaxOutputPerRead )
                    {
                        if ( ( rOutput->TotalOccurrences + ( r - l + 1 ) ) > qSetting->MaxOutputPerRead )
                        {
                            r = l + qSetting->MaxOutputPerRead - rOutput->TotalOccurrences - 1;
                        }

                        addSAToSAList ( sa_list, l, r, strand, num_mis );
                        rOutput->TotalOccurrences += r - l + 1;
                        rOutput->WithError[num_mis] += r - l + 1;
                    }

                    if ( rOutput->TotalOccurrences >= qSetting->MaxOutputPerRead )
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        else if ( ( skip_round_2 == 1 ) && ( * ( ans ) > 0xFFFFFFFD ) )
        {
            // for the case when (1) round 2 is skipped and (2) it is a bad read
            // this read will be handled by CPU
#ifndef BGS_DISABLE_NEGATIVE_STRAND
            for ( i = 0; i < readLength; i++ )
            {
                oStrandQuery[i] = soap3DnaComplement[thisQuery[readLength - i - 1]];
            }

            qInfo_Negative->ReadCode = oStrandQuery;
            qInfo_Positive->ReadStrand = QUERY_POS_STRAND;
            qInfo_Negative->ReadStrand = QUERY_NEG_STRAND;
            //Reversing the read from positive strand to negative strand.
#endif
#ifdef BGS_DISABLE_NEGATIVE_STRAND
            ProcessReadSingleStrand2 ( qInput_Positive, SRAMismatchModel[readLength], whichCase, sa_list, occ_list );
#else
            ProcessReadDoubleStrand2 ( qInput_Positive, qInput_Negative, SRAMismatchModel[readLength], SRAMismatchModel_neg[readLength], whichCase, sa_list, occ_list );
#endif
        }
        else if ( ( skip_round_2 == 0 ) && ( * ( ans ) > 0xFFFFFFFD ) )
        {
            // check the second round answer
            bool isSuperBadRead = true;

            while ( badReadPos[whichCase] < badCountOffset[whichCase] && badReadIndices[whichCase][badReadPos[whichCase] + badStartOffset[whichCase]] <= batchReadId )
            {
                if ( badReadIndices[whichCase][badReadPos[whichCase] + badStartOffset[whichCase]] == batchReadId )
                {
                    ullint pos = badReadPos[whichCase];
                    ans = badAnswers[whichCase] + ( ( pos + badStartOffset[whichCase] ) / 32 * 32 * word_per_ans_2 + ( pos + badStartOffset[whichCase] ) % 32 );

                    if ( * ( ans ) <= 0xFFFFFFFD )
                    {
                        // either no hit or number of SA ranges reported within sa_range_allowed_2
                        isSuperBadRead = false;

                        if ( * ( ans ) < 0xFFFFFFFD )
                        {
                            //Number of SA ranges reported within sa_range_allowed_2
                            for ( i = 0; i < sa_range_allowed_2; i++ )
                            {
                                if ( ( * ( ans + ( i * 2 ) * 32 ) ) >= 0xFFFFFFFD || ( * ( ans + ( i * 2 + 1 ) * 32 ) ) >= 0xFFFFFFFD )
                                {
                                    break;
                                }

                                l = * ( ans + ( i * 2 ) * 32 );
                                r = l + ( * ( ans + ( i * 2 + 1 ) * 32 ) &BGS_GPU_ANSWER_OFFSET_MASK );
                                strand = ( ( * ( ans + ( i * 2 + 1 ) * 32 ) >> ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) & 1 ) + 1;
                                num_mis = ( ( * ( ans + ( i * 2 + 1 ) * 32 ) >> ( BGS_GPU_ANSWER_OFFSET_LENGTH ) ) & 7 );

                                if ( l <= r && r <= bwt->textLength )
                                {
                                    if ( rOutput->TotalOccurrences < qSetting->MaxOutputPerRead )
                                    {
                                        if ( ( rOutput->TotalOccurrences + ( r - l + 1 ) ) > qSetting->MaxOutputPerRead )
                                        {
                                            r = l + qSetting->MaxOutputPerRead - rOutput->TotalOccurrences - 1;
                                        }

                                        addSAToSAList ( sa_list, l, r, strand, num_mis );
                                        rOutput->TotalOccurrences += r - l + 1;
                                        rOutput->WithError[num_mis] += r - l + 1;
                                    }

                                    if ( rOutput->TotalOccurrences >= qSetting->MaxOutputPerRead )
                                    {
                                        break;
                                    }
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }
                    }
                }

                badReadPos[whichCase]++;
            }

            if ( isSuperBadRead )
            {
                // super bad read
#ifndef BGS_DISABLE_NEGATIVE_STRAND
                //Reversing the read from positive strand to negative strand.
                for ( i = 0; i < readLength; i++ )
                {
                    oStrandQuery[i] = soap3DnaComplement[thisQuery[readLength - i - 1]];
                }

                qInfo_Negative->ReadCode = oStrandQuery;
                qInfo_Positive->ReadStrand = QUERY_POS_STRAND;
                qInfo_Negative->ReadStrand = QUERY_NEG_STRAND;
#endif
#ifdef BGS_DISABLE_NEGATIVE_STRAND
                ProcessReadSingleStrand2 ( qInput_Positive, SRAMismatchModel[readLength], whichCase, sa_list, occ_list );
#else
                ProcessReadDoubleStrand2 ( qInput_Positive, qInput_Negative, SRAMismatchModel[readLength], SRAMismatchModel_neg[readLength], whichCase, sa_list, occ_list );
#endif
            }
        }
    }
}


inline unsigned long long PEDumpAllOccurrence ( PEOutput * pe_out, SRAQueryInput * qInput, char reportType, int minMismatchCount, unsigned int maxOutputCount,
        int preReadLength, int readLength )
{
    unsigned long long numOfPairEndAlignment = 0;
    unsigned long long numOfAnswer = 0;
    PEPairList * pairList = pe_out->root;
    unsigned int i;
    SRASetting * qSetting = qInput->QuerySetting;
    OCC * occ = qSetting->occ;

    if ( pairList != NULL && pairList->pairsCount > 0 )
    {
        OCCReportDelimitor ( qInput );

        while ( pairList != NULL && pairList->pairsCount > 0 )
        {
            unsigned int pairsCount = pairList->pairsCount;

            for ( i = 0; i < pairsCount; i++ )
            {
                PEPairs * pePair = & ( pairList->pairs[i] );

                if ( reportType == OUTPUT_ALL_BEST && pePair->totalMismatchCount > minMismatchCount )
                { continue; }

#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
                pePair->strand_2 = 2;
#endif
                numOfAnswer++;
                numOfPairEndAlignment++;

                // output the alignment of the first read
                if ( 1 || occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );} // condition made always true by cx

                occ->occPositionCache[occ->occPositionCacheCount].tp = pePair->algnmt_1;
                occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = pePair->strand_1;
                occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                occ->occPositionCache[occ->occPositionCacheCount].occMismatch = pePair->mismatch_1;
                occ->occPositionCache[occ->occPositionCacheCount].len = preReadLength; // cx
                occ->occPositionCacheCount++;

                // output the alignment of the second read
                // OCCReportDelimitor(occ,hsp,outputFile,(unsigned long long)readId+accumReadNum);
                if ( 1 ||  occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );} // condition made always true by cx

                occ->occPositionCache[occ->occPositionCacheCount].tp = pePair->algnmt_2;
                occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = pePair->strand_2;
                occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                occ->occPositionCache[occ->occPositionCacheCount].occMismatch = pePair->mismatch_2;
                occ->occPositionCache[occ->occPositionCacheCount].len = readLength; // cx
                occ->occPositionCacheCount++;

                if ( numOfPairEndAlignment >= maxOutputCount )
                {
                    break;
                }
            }

            if ( numOfPairEndAlignment >= maxOutputCount )
            {
                break;
            }

            pairList = pairList->next;
        }
    }

    return numOfAnswer;
}

//Function hostKernel : core function of output report and host alignment.
//
//batchFirstReadId: is the readId of the first read in the batch
//                  w.r.t to the read input file. (same among all threads)
//skipFirst       : (Multi-threading) define the number of reads in the batch to skip for 1st batch-run.
//answers         : is the answer arrays used by GPU to store alignment result of the 1st batch-run.
//badAnswers      : is the answer arrays used by GPU to store alignment result of the 2st batch-run.
//badReadIndices  : is the array of the read id w.r.t to the batch.

inline uint hostKernel ( char * upkdQualities, char * upkdQueryNames, unsigned int batchFirstReadId, unsigned int skipFirst, unsigned int numQueries,
                         unsigned int * queries, unsigned int word_per_query,
                         SRASetting * qSetting, SRAIndex * aIndex, SRAQueryResultCount * rOutput,
                         SRAModel ** SRAMismatchModel, SRAModel ** SRAMismatchModel_neg,
                         char * outputFileName,
                         unsigned int ** answers,
                         unsigned long long int * alignedOcc, unsigned int * alignedReads,
                         unsigned int ** badReadIndices, unsigned int ** badAnswers,
                         unsigned int * badStartOffset, unsigned int * badCountOffset,
                         unsigned int maxReadLength, unsigned int * readLengths,
                         unsigned int * seedLengths, unsigned int numCases,
                         unsigned int sa_range_allowed_1, unsigned int sa_range_allowed_2,
                         unsigned int * readIDs, char skip_round_2, unsigned int accumReadNum,
                         unsigned int word_per_ans, char readType,
                         int insert_low, int insert_high,
                         int peStrandLeftLeg, int peStrandRightLeg,
                         int threadID, unsigned int * unAlignedIDs, unsigned int * unAlignedOcc,
                         char reportType, int peMaxOutputPerPair, uint8_t isTerminalCase,
                         ReadInputForDP * dpInput, ReadInputForDP * dpInputForNewDefault,
                         ReadInputForDP * otherSoap3Result,
                         BothUnalignedPairs * bothUnalignedPairs,
                         unsigned int maxHitNumForDP, unsigned int maxHitNumForDP2, char enableDP,
                         SRAModel ** SRAMismatchModel2, SRAModel ** SRAMismatchModel2_neg, unsigned int numCases2 )
{
    ullint numOfAnswer = 0;
    ullint pre_numOfAnswer = 0;
    uint numOfAlignedRead = 0;
    uint numOfUnAligned = 0;
    unsigned char pStrandQuery[2][MAX_READ_LENGTH];
    unsigned char * thisQuery = NULL;
    char * thisQueryName = NULL;
    unsigned char * preQuery = NULL;
    char * thisQualities = NULL;
    char * preQualities = NULL;
    // char oQualities[MAX_READ_LENGTH];
    int preFlag = 0;
    int i, j;
    ullint readId = 0;
    ullint preReadId = 0;
    ullint batchReadId;
    int readLength = MAX_READ_LENGTH;
    int preReadLength = MAX_READ_LENGTH;
    char * preQueryName = NULL;
    // for collecting the SA ranges from GPU alignment
    SAList * sa_list1 = SAListConstruct ();
    SAList * sa_list2 = SAListConstruct ();
    SAList * curr_sa_list = sa_list1;
    OCCList * occ_list1 = OCCListConstruct ();
    OCCList * occ_list2 = OCCListConstruct ();
    OCCList * curr_occ_list = occ_list1;
    PEOutput * pe_out = PEOutputConstruct ();
    BWT * bwt = aIndex->bwt;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    OCC * occ = qSetting->occ;
    PEInput * pe_in = PEInputConstruct ( bwt, hsp );
    AlgnResultArrays * algnResultArrays = hspaux->algnResultArrays;
    int outputFormat = qSetting->OutFileFormat;
    DynamicUint8Array * charArray = NULL;

    unsigned int previousTotalOccurrences;
    unsigned int previousWithError[ MAX_NUM_OF_ERROR + 1];
    int previousMinNumMismatch;
    unsigned int currentTotalOccurrences;
    unsigned int currentWithError[ MAX_NUM_OF_ERROR + 1];
    int currentMinNumMismatch;

    char dummyQuality[MAX_READ_LENGTH];
    memset ( dummyQuality, 1, sizeof ( char ) *MAX_READ_LENGTH );
    SRAQueryInput qInput_Positive;
    SRAQueryInfo qInfo_Positive;
    qInfo_Positive.ReadStrand = QUERY_POS_STRAND;
    qInfo_Positive.ReadQuality = dummyQuality;
    qInput_Positive.QueryInfo = &qInfo_Positive;
    qInput_Positive.QuerySetting = qSetting;
    qInput_Positive.AlgnmtIndex = aIndex;
    qInput_Positive.QueryOutput = rOutput;

    SRAQueryInput qInput_Negative;
    SRAQueryInfo qInfo_Negative;
    qInfo_Negative.ReadStrand = QUERY_NEG_STRAND;
    qInfo_Negative.ReadQuality = dummyQuality;
    qInput_Negative.QueryInfo = &qInfo_Negative;
    qInput_Negative.QuerySetting = qSetting;
    qInput_Negative.AlgnmtIndex = aIndex;
    qInput_Negative.QueryOutput = rOutput;

    /////////////////////////////////
    // Allocate and Fill CharMap
    unsigned char charMap[256];
    INDEXFillCharMap ( charMap );

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    {
        charArray = DynamicUint8ArrayConstruct ();
    }

#ifdef BGS_ROUND_BREAKDOWN_TIME
    double start_time = setStartTime ();
#endif
    FILE * outputFile = NULL;

    //Output file for SAM output is handled by SAM API
    //=> OutFilePtr is not used for SAM API output format.
    if ( outputFileName != NULL )
    {
        switch ( outputFormat )
        {
            case SRA_OUTPUT_FORMAT_SAM_API:
                break;

            default:
                outputFile = ( FILE * ) fopen ( outputFileName, "a" );
        }
    }

    qSetting->OutFilePtr = outputFile;
    qSetting->OutFileName = outputFileName;
    int currNumMismatch = qSetting->MaxError;
    pe_in->insertLbound = insert_low;
    pe_in->insertUbound = insert_high;
    pe_in->strandLeftLeg = peStrandLeftLeg;
    pe_in->strandRightLeg = peStrandRightLeg;
    unsigned int numOfPairEndAlignment = 0;
    // printf("[hostkernel] numCases = %u; numQueries = %u; accumReadNum = %u; skip_round_2 = %i\n", numCases, numQueries, accumReadNum, (int) skip_round_2);
    // printf("[hostkernel] sa_range_allowed_1 = %u; sa_range_allowed_2 = %u; word_per_ans = %u; word_per_ans_2 = %i\n", sa_range_allowed_1, sa_range_allowed_2, word_per_ans, word_per_ans_2);
    // printf("[hostkernel] batchFirstReadId = %u; skipFirst = %u; numQueries = %u\n", batchFirstReadId, skipFirst, numQueries);
    // ======================================================== //
    // If this is the terminal case (i.e. isTerminalCase = 1)   //
    // and PERFORM_DP = 1
    // THEN the unaligned reads may need to proceed to DP       //
    // ======================================================== //
    bool needProceedDP = ( isTerminalCase == 1 ) && ( enableDP == 1 );
    // ======================================================== //
    // If this is the terminal case (i.e. isTerminalCase = 1)   //
    // and (outputFormat == SRA_OUTPUT_FORMAT_SAM or            //
    //      outputFormat == SRA_OUTPUT_FORMAT_SAM_API)          //
    // THEN the unaligned reads are needed to output            //
    // ======================================================== //
    bool needOutputUnalignedReads = ( isTerminalCase == 1 ) && ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API );
    // ========================================================= //
    // If alignment type is ALL-VALID and output format is SAM,  //
    // then mapping quality score is required to compute.        //
    // ========================================================= //
    bool needOutputMAPQ = ( reportType == OUTPUT_ALL_VALID || reportType == OUTPUT_ALL_BEST ) && ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API );
    // ========================================================= //
    // For LONG READ (read length > 120),                        //
    // only the first 100 bases are aligned                      //
    // then for each alignment hit, extend the alignment and     //
    // check the total # of mismatches (by using popcount)       //
    // ========================================================= //
    bool longReadMode = ( seedLengths != NULL );
    unsigned int origHitLimit = -1;

    if ( longReadMode )
    {
        // remove the limit on the number of alignments for each read
        origHitLimit = qSetting->MaxOutputPerRead;
        qSetting->MaxOutputPerRead = -1;
    }

    // ========================================================================================= //
    // If OUTPUT_ANS_TO_ARRAY is defined, then all results have to be stored in the global array //
    // ========================================================================================= //
    uint * badReadPos = ( uint * ) malloc ( numCases * sizeof ( uint ) );

    for ( j = 0; j < numCases; j++ )
    { badReadPos[j] = 0; }

    for ( j = 0; j < numQueries; ++j )
    {
        numOfPairEndAlignment = 0;
        // OBTAIN THE ALIGNMENT RESULT
        batchReadId = skipFirst + j;
        preReadId = readId;
        readId = readIDs[batchFirstReadId + skipFirst + j];
        preQuery = thisQuery;
        thisQuery = pStrandQuery[preFlag];
        preQualities = thisQualities;
        thisQualities = upkdQualities + maxReadLength * ( readId - 1 );
        retrieve2BWTQueryFormat ( queries, batchFirstReadId + skipFirst + j, word_per_query, thisQuery );
        preFlag++;
        preFlag = preFlag % 2;
        preQueryName = thisQueryName;
        thisQueryName = upkdQueryNames + hspaux->maxLenReadName * ( readId - 1 );
        preReadLength = readLength;

        if ( longReadMode )
        { readLength = seedLengths[batchFirstReadId + batchReadId]; }
        else
        { readLength = readLengths[batchFirstReadId + batchReadId]; }

        pre_numOfAnswer = numOfAnswer;
        qInfo_Positive.ReadName = thisQueryName;
        qInfo_Negative.ReadName = thisQueryName;

        if ( ( readType == PAIR_END_READ ) && ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN ) )
        {
            qInfo_Positive.ReadId = ( unsigned long long ) ( readId + accumReadNum + 1 ) / 2;
            qInfo_Negative.ReadId = ( unsigned long long ) ( readId + accumReadNum + 1 ) / 2;
        }
        else
        {
            qInfo_Positive.ReadId = ( unsigned long long ) readId + accumReadNum;
            qInfo_Negative.ReadId = ( unsigned long long ) readId + accumReadNum;
        }

        qInfo_Positive.ReadCode = thisQuery;
        qInfo_Positive.ReadLength = readLength;
        qInfo_Positive.ReadQuality = thisQualities;
        qInfo_Negative.ReadLength = readLength;
        qInfo_Negative.ReadQuality = thisQualities;

        if ( readType == SINGLE_READ )
        {
            SAListReset ( sa_list1 );
            OCCListReset ( occ_list1 );

            if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API && outputFileName != NULL )
            {
                OCCReportDelimitor ( &qInput_Positive );
            }
        }
        else if ( readType == PAIR_END_READ )
        {
            if ( j % 2 == 0 )
            {
                SAListReset ( sa_list1 );
                SAListReset ( sa_list2 );
                curr_sa_list = sa_list1;
                OCCListReset ( occ_list1 );
                OCCListReset ( occ_list2 );
                curr_occ_list = occ_list1;
            }
            else
            {
                curr_sa_list = sa_list2;
                curr_occ_list = occ_list2;
            }
        }

        bool only_get_round1_ans = ( readType == SINGLE_READ && ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API ) &&
                                     ( !longReadMode ) && ( reportType == OUTPUT_RANDOM_BEST ||
                                             reportType == OUTPUT_UNIQUE_BEST ) );
        bool isMoreThanSA1 = false; // whether the number of SA ranges > sa_range_allowed_1
        // initialize to false
        collect_all_answers ( curr_sa_list, curr_occ_list, answers, batchReadId,
                              word_per_ans, &qInput_Positive, &qInput_Negative,
                              sa_range_allowed_1, sa_range_allowed_2, skip_round_2,
                              numCases, readLength, thisQuery, SRAMismatchModel, SRAMismatchModel_neg, charMap,
                              badReadIndices, badAnswers, badStartOffset, badCountOffset,
                              badReadPos, only_get_round1_ans, isMoreThanSA1 );
        bool isUnique = ( ( rOutput->TotalOccurrences == 1 ) && ( !isMoreThanSA1 ) ); // whether the hit is unique
        bool moreThanOne = ( ( rOutput->TotalOccurrences > 1 ) || isMoreThanSA1 ); // whether there is more than one hit
        bool nohit = ( curr_sa_list->curr_size == 0 && curr_occ_list->curr_size == 0 && ( !isMoreThanSA1 ) ); // whether there is no hit

        // special handling for MAPQ for SINGLE READ
        if ( rOutput->TotalOccurrences > 0 && readType == SINGLE_READ && needOutputMAPQ && currNumMismatch < 4 )
        {
            //Store mismatch statistics and total number of occurrences for later process
            currentTotalOccurrences = rOutput->TotalOccurrences;
            memcpy ( currentWithError, rOutput->WithError, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );

            currentMinNumMismatch = 999;

            for ( int mismatchi = 0; mismatchi <= currNumMismatch; mismatchi++ )
            {
                if ( currentWithError[mismatchi] > 0 )
                {
                    currentMinNumMismatch = mismatchi;
                    break;
                }
            }

            bool need_special_handle = ( currentMinNumMismatch == currNumMismatch );

            if ( need_special_handle )
            {
                if ( isTerminalCase != 1 )
                {
                    unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j;
                    continue;
                }
                else
                {
                    // add the read ID to the array for further process by DP
                    // addReadIDToBothUnalignedPairs(bothUnalignedPairs, batchFirstReadId+batchReadId);
                    // continue;
                    // perform (mismatch+1) alignment on the read
                    // fprintf(stderr, "perform (mismatch+1) alignment on the read\n");
                    SAListReset ( curr_sa_list );
                    OCCListReset ( curr_occ_list );
                    ProcessOneMoreMismatchAllCases ( charMap, &qInput_Positive, &qInput_Negative, SRAMismatchModel2[readLength], SRAMismatchModel2_neg[readLength], curr_sa_list, curr_occ_list, numCases2 );

                    //Store mismatch statistics and total number of occurrences for later process
                    currentTotalOccurrences = rOutput->TotalOccurrences;
                    memcpy ( currentWithError, rOutput->WithError, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );

                    // fprintf(stderr, "size of curr_sa_list = %u; size of curr_occ_list = %u\n", curr_sa_list->curr_size, curr_occ_list->curr_size);
                    isUnique = ( rOutput->TotalOccurrences == 1 ); // whether the hit is unique
                    moreThanOne = ( rOutput->TotalOccurrences > 1 ); // whether there is more than one hit
                    nohit = ( curr_sa_list->curr_size == 0 && curr_occ_list->curr_size == 0 ); //
                }
            }
        }

        // =========================================================== //
        // For LONG READ (read length > 120),                          //
        // for each alignment hit, extend the alignment and            //
        // check the total # of mismatches (by using popcount)         //
        // if the mismatch # <=  MISMATCH_RATIO_FOR_LONG * readlength  //
        // then the alignment is valid.                                //
        // Only keep the valid alignment results                       //
        // =========================================================== //

        if ( longReadMode && ( curr_sa_list->curr_size > 0 || curr_occ_list->curr_size > 0 ) )
        {
            // transfer all sa ranges to occ list.
            transferAllSAToOcc ( curr_sa_list, curr_occ_list, bwt );
            bool only_keep_best_ans = ( readType == SINGLE_READ && reportType != OUTPUT_ALL_VALID && ( !needOutputMAPQ ) );
            int max_mismatch_allowed = ( int ) ceil ( MISMATCH_RATIO_FOR_LONG * readLengths[batchFirstReadId + batchReadId] );

            if ( needOutputMAPQ ) { max_mismatch_allowed = max_mismatch_allowed * 2; }

            int min_seed_mismatch_allowed = ( readType == SINGLE_READ && reportType != OUTPUT_ALL_VALID && ( !needOutputMAPQ ) ) ? currNumMismatch : 0;
            // -----------------------------------------------------------------
            // by adding the following would be more precise, but it will need longer time
            // -----------------------------------------------------------------
            /*
            if (readType==SINGLE_READ && reportType!=OUTPUT_ALL_VALID && (!isTerminalCase)) {
                  max_mismatch_allowed = currNumMismatch;
            }
             */
            // -----------------------------------------------------------------
            validateAlignments ( curr_occ_list, thisQuery, seedLengths[batchFirstReadId + batchReadId], readLengths[batchFirstReadId + batchReadId], hsp, only_keep_best_ans, min_seed_mismatch_allowed, max_mismatch_allowed, origHitLimit, threadID );
            isUnique = ( curr_occ_list->curr_size == 1 );
            moreThanOne = ( curr_occ_list->curr_size > 1 );
            nohit = ( curr_occ_list->curr_size == 0 );

            // restrict on the number of alignment hits
            if ( curr_occ_list->curr_size > origHitLimit )
            { curr_occ_list->curr_size = origHitLimit; }
        }

        if ( readType == SINGLE_READ )
        {
            // set back the original length of the read
            if ( longReadMode )
            { qInfo_Positive.ReadLength = readLengths[batchFirstReadId + batchReadId]; }

            transferAllSAToOcc ( curr_sa_list, curr_occ_list, bwt );

            // report the read with no hit
            if ( nohit )
            {
                unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j;

                if ( needProceedDP )
                {
                    // add the read ID to the array for further process by DP
                    addReadIDToBothUnalignedPairs ( bothUnalignedPairs, batchFirstReadId + batchReadId );
                }
            }

            if ( reportType == OUTPUT_UNIQUE_BEST )
            {
                if ( isUnique )
                {
                    if ( outputFileName != NULL )
                    {
                        if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API )
                        {
                            if ( 1 || occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( &qInput_Positive );} // condition made always true by cx

                            occ->occPositionCache[occ->occPositionCacheCount].tp = ( curr_occ_list->occ[0] ).ambPosition;
                            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = ( curr_occ_list->occ[0] ).strand;
                            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = ( curr_occ_list->occ[0] ).mismatchCount;
                            occ->occPositionCache[occ->occPositionCacheCount].len = readLength; // cx
                            occ->occPositionCacheCount++;
                        }
                        else
                        {
                            SingleAnsOutputSAMAPI ( &qInput_Positive,
                                                    ( curr_occ_list->occ[0] ).strand, ( curr_occ_list->occ[0] ).ambPosition, ( curr_occ_list->occ[0] ).mismatchCount, 1 );
                        }
                    }
                    else
                    {
                        addOCCToArray ( algnResultArrays->algnArrays[threadID], qInfo_Positive.ReadId, ( curr_occ_list->occ[0] ).ambPosition, ( curr_occ_list->occ[0] ).strand, 1, ( curr_occ_list->occ[0] ).mismatchCount );
                    }

                    numOfAnswer++;
                }
                else if ( moreThanOne && outputFormat == SRA_OUTPUT_FORMAT_SAM_API && outputFileName != NULL )
                {
                    // output no hit for SAM format
                    noAnsOutputSAMAPI ( &qInput_Positive );
                    continue;
                }
            }
            else if ( reportType == OUTPUT_RANDOM_BEST )
            {
                if ( !nohit )
                {
                    if ( outputFileName != NULL )
                    {
                        if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API )
                        {
                            if ( 1 || occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( &qInput_Positive );} // condition made always true by cx

                            occ->occPositionCache[occ->occPositionCacheCount].tp = ( curr_occ_list->occ[0] ).ambPosition;
                            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = ( curr_occ_list->occ[0] ).strand;
                            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = ( curr_occ_list->occ[0] ).mismatchCount;
                            occ->occPositionCache[occ->occPositionCacheCount].len = readLength; // cx
                            occ->occPositionCacheCount++;
                        }
                        else
                        {
                            SingleAnsOutputSAMAPI ( &qInput_Positive,
                                                    ( curr_occ_list->occ[0] ).strand, ( curr_occ_list->occ[0] ).ambPosition, ( curr_occ_list->occ[0] ).mismatchCount, 1 );
                        }
                    }
                    else
                    {
                        addOCCToArray ( algnResultArrays->algnArrays[threadID], qInfo_Positive.ReadId, ( curr_occ_list->occ[0] ).ambPosition, ( curr_occ_list->occ[0] ).strand, 1, ( curr_occ_list->occ[0] ).mismatchCount );
                    }

                    numOfAnswer++;
                }
            }
            else
            {
                if ( outputFileName != NULL )
                {
                    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
                    {
                        numOfAnswer += curr_occ_list->curr_size;

                        if ( curr_occ_list->curr_size > 0 )
                        {
                            OCCOutputSAMAPI ( &qInput_Positive, curr_occ_list, charArray, readLength, reportType );
                        }
                    }
                    else
                    {
                        // long read mode and succint
                        for ( i = 0; i < curr_occ_list->curr_size; i++ )
                        {
                            if ( 1 || occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( &qInput_Positive );} // condition made always true by cx

                            occ->occPositionCache[occ->occPositionCacheCount].tp = ( curr_occ_list->occ[i] ).ambPosition;
                            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = ( curr_occ_list->occ[i] ).strand;
                            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = ( curr_occ_list->occ[i] ).mismatchCount;
                            occ->occPositionCache[occ->occPositionCacheCount].len = readLength; // cx
                            occ->occPositionCacheCount++;
                        }

                        numOfAnswer += curr_occ_list->curr_size;
                    }
                }
                else
                {
                    numOfAnswer += curr_occ_list->curr_size;

                    for ( int i = 0; i < curr_occ_list->curr_size; i++ )
                    {
                        addOCCToArray ( algnResultArrays->algnArrays[threadID], qInfo_Positive.ReadId, ( curr_occ_list->occ[i] ).ambPosition, ( curr_occ_list->occ[i] ).strand, 1, ( char ) ( curr_occ_list->occ[i] ).mismatchCount );
                    }
                }
            }

            if ( numOfAnswer > pre_numOfAnswer )
            {
                numOfAlignedRead++;
            }
            else if ( needOutputUnalignedReads && ( outputFileName != NULL ) && ( !needProceedDP ) )
            {
                // output no hit for SAM format
                noAnsOutputSAMAPI ( &qInput_Positive );
            }
        }
        else if ( j % 2 == 0 )
        {
            // FOR FIRST READ OF THE PAIR-END READ
            if ( rOutput->TotalOccurrences == 0 &&
                    ( !needOutputUnalignedReads ) && ( !needProceedDP ) )
            {
                // no hits for the first read
                // and no need to output unaligned reads
                // and no need for unaligned reads to proceed DP
                // then the searching for the second read can be skipped
                j++;
                unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j - 1;
                unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j;
                continue;
            }

            //Store mismatch statistics and total number of occurrences for later process
            previousTotalOccurrences = rOutput->TotalOccurrences;
            memcpy ( previousWithError, rOutput->WithError, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );

            // This block is going to be removed for further performance enhancement.
            if ( needOutputMAPQ && ( !needOutputUnalignedReads ) && ( !needProceedDP ) )
            {
                //int totalNumAns1 = 0;
                currentMinNumMismatch = 999;

                for ( int mismatchi = 0; mismatchi <= currNumMismatch; mismatchi++ )
                {
                    if ( previousWithError[mismatchi] > 0 )
                    {
                        currentMinNumMismatch = mismatchi;
                        break;
                    }
                }

                bool need_handle_first = ( ( currentMinNumMismatch == currNumMismatch ) &&
                                           ( previousTotalOccurrences < qSetting->MaxOutputPerRead ) );

                if ( need_handle_first )
                {
                    j++;
                    unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j - 1;
                    unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j;
                    continue;
                }
            }
        }
        else if ( j % 2 == 1 )
        {
            //Store mismatch statistics and total number of occurrences for later process
            currentTotalOccurrences = rOutput->TotalOccurrences;
            memcpy ( currentWithError, rOutput->WithError, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );

            // FOR SECOND READ OF THE PAIR-END READ
            int first_X0 = 0;
            int first_X1 = 0;
            int second_X0 = 0;
            int second_X1 = 0;

            // get the statistics
            // FOR DEBUG (OCC)
            /*
            if (needProceedDP && previousTotalOccurrences>0 ) {
                  printf("%u %u %u %u %u\n", readId-2, sa_list1->curr_size, totalNumOcc(sa_list1), occ_list1->curr_size, totalNumOcc(sa_list1)+occ_list1->curr_size);
            }
            if (needProceedDP && currentTotalOccurrences>0) {
                  printf("%u %u %u %u %u\n", readId-1, sa_list2->curr_size, totalNumOcc(sa_list2), occ_list2->curr_size, totalNumOcc(sa_list2)+occ_list2->curr_size);
            }
            */
            // if need to proceed the semi-global DP, then select those pairs of reads
            // which one end has alignments but another has not.

            if ( needOutputMAPQ )
            {
                bool need_handle_first = false;
                bool need_handle_second = false;

                if ( previousTotalOccurrences > 0 )
                {
                    previousMinNumMismatch = 999;

                    for ( int mismatchi = 0; mismatchi <= currNumMismatch; mismatchi++ )
                    {
                        if ( previousWithError[mismatchi] > 0 )
                        {
                            previousMinNumMismatch = mismatchi;
                            break;
                        }
                    }

                    need_handle_first = previousTotalOccurrences > 0 &&
                                        previousMinNumMismatch == currNumMismatch &&
                                        previousTotalOccurrences < qSetting->MaxOutputPerRead &&
                                        currNumMismatch < 4;

                    if ( need_handle_first )
                    {
                        // need to perform (mismatch+1) alignment on the first read
                        SAListReset ( sa_list1 );
                        OCCListReset ( occ_list1 );
                        qInfo_Positive.ReadCode = preQuery;
                        qInfo_Positive.ReadLength = preReadLength;
                        ProcessOneMoreMismatchAllCases ( charMap, &qInput_Positive, &qInput_Negative, SRAMismatchModel2[preReadLength], SRAMismatchModel2_neg[preReadLength], sa_list1, occ_list1, numCases2 );
                        // back to the normal values
                        qInfo_Positive.ReadCode = thisQuery;
                        qInfo_Positive.ReadLength = readLength;

                        //Store mismatch statistics and total number of occurrences for later process
                        previousTotalOccurrences = rOutput->TotalOccurrences;
                        memcpy ( previousWithError, rOutput->WithError, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );
                    }

                    first_X0 = previousWithError[previousMinNumMismatch];
                    first_X1 = previousWithError[previousMinNumMismatch + 1];

                    if ( needProceedDP )
                    {
                        hspaux->x0_array[batchFirstReadId + batchReadId - 1] = first_X0;
                        hspaux->x1_array[batchFirstReadId + batchReadId - 1] = first_X1;
                        hspaux->mismatch_array[batchFirstReadId + batchReadId - 1] = previousMinNumMismatch;
                    }
                }

                if ( currentTotalOccurrences > 0 )
                {
                    currentMinNumMismatch = 999;

                    for ( int mismatchi = 0; mismatchi <= currNumMismatch; mismatchi++ )
                    {
                        if ( currentWithError[mismatchi] > 0 )
                        {
                            currentMinNumMismatch = mismatchi;
                            break;
                        }
                    }

                    need_handle_second = currentMinNumMismatch == currNumMismatch &&
                                         currentTotalOccurrences < qSetting->MaxOutputPerRead &&
                                         currNumMismatch < 4;

                    if ( need_handle_second )
                    {
                        // need to perform (mismatch+1) alignment on the second read
                        SAListReset ( sa_list2 );
                        OCCListReset ( occ_list2 );
                        qInfo_Positive.ReadCode = thisQuery;
                        qInfo_Positive.ReadLength = readLength;
                        ProcessOneMoreMismatchAllCases ( charMap, &qInput_Positive, &qInput_Negative, SRAMismatchModel2[readLength], SRAMismatchModel2_neg[readLength], sa_list2, occ_list2, numCases2 );

                        //Store mismatch statistics and total number of occurrences for later process
                        currentTotalOccurrences = rOutput->TotalOccurrences;
                        memcpy ( currentWithError, rOutput->WithError, sizeof ( unsigned int ) * ( MAX_NUM_OF_ERROR + 1 ) );
                    }

                    second_X0 = currentWithError[currentMinNumMismatch];
                    second_X1 = currentWithError[currentMinNumMismatch + 1];

                    if ( needProceedDP )
                    {
                        hspaux->x0_array[batchFirstReadId + batchReadId] = second_X0;
                        hspaux->x1_array[batchFirstReadId + batchReadId] = second_X1;
                        hspaux->mismatch_array[batchFirstReadId + batchReadId] = currentMinNumMismatch;
                    }
                }

                // obtain the value of X0 and X1 for both ends
                // X0: the number of best hits
                // X1: the number of second best hits
                // and record down the x0, x1 and minmismatch
            }

            if ( needProceedDP )
            {
#ifdef SKIP_DEFAULT_DP

                if ( previousTotalOccurrences == 0 || currentTotalOccurrences == 0 )
                {
                    // either first read or second read do not have any hit
                    addReadIDToBothUnalignedPairs ( bothUnalignedPairs, batchFirstReadId + batchReadId - 1 );
                    continue;
                }

#else

                if ( previousTotalOccurrences > 0 && currentTotalOccurrences == 0 )
                {
                    // first read has hit but second read has not
                    if ( occ_list1->curr_size <= maxHitNumForDP && ( tooManyNumOcc ( sa_list1, maxHitNumForDP - occ_list1->curr_size ) == 0 ) )
                    {
                        // the number of alignments <= threshold
                        keepSoap3Results ( batchFirstReadId + batchReadId - 1, sa_list1, occ_list1, dpInput, hspaux );
                        addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId - 1, sa_list1->sa, sa_list1->curr_size,
                                              occ_list1->occ, occ_list1->curr_size );
                        continue;
                    }
                    else
                    {
                        unsigned int hitNum;

                        if ( needOutputMAPQ )
                        {
                            // consider the best and the second best hits
                            hitNum = retainAllBestAndSecBest ( sa_list1, occ_list1 );
                            /*
                            if (hitNum > maxHitNumForDP) {
                                  // consider those best hits
                                  hitNum = retainAllBest(sa_list1, occ_list1);
                            }*/
                        }
                        else
                        {
                            // consider those best hits
                            hitNum = retainAllBest ( sa_list1, occ_list1 );
                        }

                        if ( hitNum <= maxHitNumForDP )
                        {
                            keepSoap3Results ( batchFirstReadId + batchReadId - 1, sa_list1, occ_list1, dpInput, hspaux );
                            addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId - 1, sa_list1->sa, sa_list1->curr_size,
                                                  occ_list1->occ, occ_list1->curr_size );
                        }
                        else
                        {
                            keepSoap3Results ( batchFirstReadId + batchReadId - 1, sa_list1, occ_list1, dpInputForNewDefault, hspaux );
                            addToReadInputForDP ( dpInputForNewDefault, batchFirstReadId + batchReadId - 1, sa_list1->sa, sa_list1->curr_size,
                                                  occ_list1->occ, occ_list1->curr_size );
                        }

                        continue;
                    }
                }
                else if ( currentTotalOccurrences > 0 && previousTotalOccurrences == 0 )
                {
                    // second read has hit but first read has not
                    if ( occ_list2->curr_size <= maxHitNumForDP2 && ( tooManyNumOcc ( sa_list2, maxHitNumForDP2 - occ_list2->curr_size ) == 0 ) )
                    {
                        keepSoap3Results ( batchFirstReadId + batchReadId, sa_list2, occ_list2, dpInput, hspaux );
                        addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId, sa_list2->sa, sa_list2->curr_size,
                                              occ_list2->occ, occ_list2->curr_size );
                        continue;
                    }
                    else
                    {
                        unsigned int hitNum;

                        if ( needOutputMAPQ )
                        {
                            // consider the best and the second best hits
                            hitNum = retainAllBestAndSecBest ( sa_list2, occ_list2 );
                            /*
                            if (hitNum > maxHitNumForDP2) {
                                  // consider those best hits
                                  hitNum = retainAllBest(sa_list2, occ_list2);
                            }*/
                        }
                        else
                        {
                            // consider those best hits
                            hitNum = retainAllBest ( sa_list2, occ_list2 );
                        }

                        if ( hitNum <= maxHitNumForDP2 )
                        {
                            keepSoap3Results ( batchFirstReadId + batchReadId, sa_list2, occ_list2, dpInput, hspaux );
                            addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId, sa_list2->sa, sa_list2->curr_size,
                                                  occ_list2->occ, occ_list2->curr_size );
                        }
                        else
                        {
                            keepSoap3Results ( batchFirstReadId + batchReadId, sa_list2, occ_list2, dpInputForNewDefault, hspaux );
                            addToReadInputForDP ( dpInputForNewDefault, batchFirstReadId + batchReadId, sa_list2->sa, sa_list2->curr_size,
                                                  occ_list2->occ, occ_list2->curr_size );
                        }

                        continue;
                    }
                }
                else if ( previousTotalOccurrences == 0 && currentTotalOccurrences == 0 )
                {
                    // both first read and second read do not have any hit
                    addReadIDToBothUnalignedPairs ( bothUnalignedPairs, batchFirstReadId + batchReadId - 1 );
                    continue;
                }

#endif
            }

            // transfer all the SA ranges in sa_list to occ_list
            transferAllSAToOcc ( sa_list1, occ_list1, bwt );
            transferAllSAToOcc ( sa_list2, occ_list2, bwt );
            char pair_alignment_exist = 0;
            PEPairList * pairList;
            unsigned int numPEAlgnmt = 0;
            // pair-end alignment result exists
            unsigned int peAlgnmtMismatchStats[MAX_NUM_OF_ERROR * 2];
            memset ( peAlgnmtMismatchStats, 0, sizeof ( unsigned int ) *MAX_NUM_OF_ERROR * 2 );

            // fprintf(stderr, "size of occ_list1: %i size of occ_list2: %i\n", occ_list1->curr_size, occ_list2->curr_size);

            if ( previousTotalOccurrences > 0 && currentTotalOccurrences > 0 )
            {
                // FIND THE VALID ALIGNMENT PAIRS
                pe_in->patternLength = readLength;
                PEMappingOccurrences ( pe_in, pe_out,
                                       occ_list1->occ, occ_list1->curr_size,
                                       occ_list2->occ, occ_list2->curr_size );
                unsigned int num_minMismatch = 0; // the number of optimal pairs with the same minimum sum of mismatches
                unsigned int num_soMinMismatch = 0; // the number of sub-optimal pairs with the same minimum sum of mismatches

                if ( pe_out->root->pairsCount > 0 )
                {
                    // pair-end alignment result exists
                    PEPairs * optimal;
                    PEPairs * suboptimal;
                    char min_totalMismatchCount = 127;
                    char secMin_totalMismatchCount = 127;

                    if ( reportType == OUTPUT_RANDOM_BEST || reportType == OUTPUT_ALL_BEST
                            || reportType == OUTPUT_UNIQUE_BEST || outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
                    {
                        // find the minimum total number of mismatches for the resulting paired-end alignments (the best pair)
                        char first_isBestHit = 0;
                        char second_isBestHit = 0;

                        numPEAlgnmt = PEStatsPEOutput ( pe_out, &optimal, &suboptimal, peAlgnmtMismatchStats );

                        if ( optimal != NULL )
                        {
                            min_totalMismatchCount = optimal->totalMismatchCount;
                            num_minMismatch = peAlgnmtMismatchStats[min_totalMismatchCount];
                        }

                        if ( suboptimal != NULL )
                        {
                            secMin_totalMismatchCount = suboptimal->totalMismatchCount;
                            num_soMinMismatch = peAlgnmtMismatchStats[secMin_totalMismatchCount];
                        }

                        // check whether the each end of the best paired-end alignment is
                        // one of the best hits for each end.
                        first_isBestHit = ( previousMinNumMismatch == optimal->mismatch_1 ) ? 1 : 0;
                        second_isBestHit = ( currentMinNumMismatch == optimal->mismatch_2 ) ? 1 : 0;

                        // OUTPUT THE RESULT
                        if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
                        {
                            if ( reportType == OUTPUT_RANDOM_BEST && num_minMismatch > 0 )
                            {
                                pairOutputSAMAPI ( &qInput_Positive, pe_out, optimal,
                                                   preQuery, thisQuery, preQualities, thisQualities,
                                                   preReadLength, readLength, preQueryName, thisQueryName,
                                                   min_totalMismatchCount, secMin_totalMismatchCount, 0, charArray, peMaxOutputPerPair,
                                                   -1, -1, -1, -1, num_minMismatch,
                                                   first_isBestHit, second_isBestHit, numPEAlgnmt );
                                numOfAnswer++;
                                numOfAlignedRead += 2;
                                pair_alignment_exist = 1;
                            }
                            else if ( reportType == OUTPUT_UNIQUE_BEST && num_minMismatch == 1 )
                            {
                                pairOutputSAMAPI ( &qInput_Positive, pe_out, optimal,
                                                   preQuery, thisQuery, preQualities, thisQualities,
                                                   preReadLength, readLength, preQueryName, thisQueryName,
                                                   min_totalMismatchCount, secMin_totalMismatchCount, 0, charArray, peMaxOutputPerPair,
                                                   1, 1, -1, -1, num_minMismatch,
                                                   first_isBestHit, second_isBestHit, numPEAlgnmt );
                                numOfAnswer++;
                                numOfAlignedRead += 2;
                                pair_alignment_exist = 1;
                            }
                            else if ( reportType == OUTPUT_ALL_VALID || reportType == OUTPUT_ALL_BEST )
                            {
                                pairOutputSAMAPI ( &qInput_Positive, pe_out, optimal,
                                                   preQuery, thisQuery, preQualities, thisQualities,
                                                   preReadLength, readLength, preQueryName, thisQueryName,
                                                   min_totalMismatchCount, secMin_totalMismatchCount,
                                                   1, charArray, peMaxOutputPerPair,
                                                   first_X0, second_X0, first_X1, second_X1, num_minMismatch,
                                                   first_isBestHit, second_isBestHit, numPEAlgnmt );

                                if ( numPEAlgnmt <= peMaxOutputPerPair )
                                { numOfAnswer += numPEAlgnmt; }
                                else
                                { numOfAnswer += peMaxOutputPerPair; }

                                numOfAlignedRead += 2;
                                pair_alignment_exist = 1;
                            }
                            else
                            {
                                // reportType == OUTPUT_UNIQUE_BEST but num_minMismatch > 1
                                // here do not report any alignment
                                pairOutputSAMAPI ( &qInput_Positive, pe_out, NULL,
                                                   preQuery, thisQuery, preQualities, thisQualities,
                                                   preReadLength, readLength, preQueryName, thisQueryName,
                                                   min_totalMismatchCount, secMin_totalMismatchCount, 0, charArray, peMaxOutputPerPair,
                                                   -1, -1, -1, -1, num_minMismatch,
                                                   first_isBestHit, second_isBestHit, numPEAlgnmt );
                                pair_alignment_exist = 1;
                            }
                        }
                        else if ( ( reportType == OUTPUT_RANDOM_BEST && num_minMismatch > 0 ) ||
                                  ( reportType == OUTPUT_UNIQUE_BEST && num_minMismatch == 1 ) )
                        {
#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
                            optimal->strand_2 = 2;
#endif
                            OCCReportDelimitor ( &qInput_Positive );

                            // output the alignment of the first read
                            // qInfo_Positive.ReadStrand = strand;
                            if ( 1 || occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( &qInput_Positive );} // condition made always true by cx

                            occ->occPositionCache[occ->occPositionCacheCount].tp = optimal->algnmt_1;
                            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = optimal->strand_1;
                            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = optimal->mismatch_1;
                            occ->occPositionCache[occ->occPositionCacheCount].len = preReadLength; // cx
                            occ->occPositionCacheCount++;

                            // output the alignment of the second read
                            // OCCReportDelimitor(occ,hsp,outputFile,(unsigned long long)readId+accumReadNum);
                            if ( 1 || occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( &qInput_Positive );} // condition made always true by cx

                            occ->occPositionCache[occ->occPositionCacheCount].tp = optimal->algnmt_2;
                            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = optimal->strand_2;
                            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = optimal->mismatch_2;
                            occ->occPositionCache[occ->occPositionCacheCount].len = readLength; // cx
                            occ->occPositionCacheCount++;
                            numOfAnswer++;
                            numOfAlignedRead += 2;
                            pair_alignment_exist = 1;
                        }
                        else
                        {
                            // reportType == OUTPUT_UNIQUE_BEST but num_minMismatch > 1
                            // not proceeding to next step
                            pair_alignment_exist = 1;
                        }
                    }

                    if ( ( reportType == OUTPUT_ALL_BEST || reportType == OUTPUT_ALL_VALID ) &&
                            ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API ) )
                    {
                        unsigned long long dumpOcc = PEDumpAllOccurrence ( pe_out, &qInput_Positive, reportType, min_totalMismatchCount, peMaxOutputPerPair,
                                                     preReadLength, readLength );

                        if ( dumpOcc > 0 )
                        {
                            numOfAnswer += dumpOcc;
                            numOfPairEndAlignment += dumpOcc;
                            numOfAlignedRead += 2;
                            pair_alignment_exist = 1;
                            numPEAlgnmt = numOfPairEndAlignment;
                        }
                    }
                }
            }

            if ( numPEAlgnmt == 0 )
            {
                if ( !needProceedDP )
                {
                    unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j - 1;
                    unAlignedIDs[numOfUnAligned++] = batchFirstReadId + skipFirst + j;
                }
                else
                {
                    // only case remain: both reads can be aligned
                    // but the insert size is out of range.
                    if ( previousTotalOccurrences > 0 && currentTotalOccurrences > 0 )
                    {
                        // first read has hits AND second read has hits
                        if ( previousTotalOccurrences > maxHitNumForDP )
                        {
                            if ( needOutputMAPQ )
                            {
                                // consider the best and the second best hits
                                previousTotalOccurrences = ( int ) retainAllBestAndSecBest ( sa_list1, occ_list1 );
                                /*
                                if (previousTotalOccurrences > maxHitNumForDP) {
                                      // consider those best hits
                                      previousTotalOccurrences = (int) retainAllBest(sa_list1, occ_list1);
                                }*/
                            }
                            else
                            {
                                // consider those best hits
                                previousTotalOccurrences = ( int ) retainAllBest ( sa_list1, occ_list1 );
                            }
                        }

                        if ( currentTotalOccurrences > maxHitNumForDP2 )
                        {
                            if ( needOutputMAPQ )
                            {
                                // consider the best and the second best hits
                                currentTotalOccurrences = ( int ) retainAllBestAndSecBest ( sa_list2, occ_list2 );
                                /*
                                if (currentTotalOccurrences > maxHitNumForDP2) {
                                      // consider those best hits
                                      currentTotalOccurrences = (int) retainAllBest(sa_list2, occ_list2);
                                }*/
                            }
                            else
                            {
                                // consider those best hits
                                currentTotalOccurrences = ( int ) retainAllBest ( sa_list2, occ_list2 );
                            }
                        }

                        if ( previousTotalOccurrences <= maxHitNumForDP && currentTotalOccurrences <= maxHitNumForDP2 )
                        {
                            keepSoap3Results ( batchFirstReadId + batchReadId - 1, sa_list1, occ_list1, dpInput, hspaux );
                            addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId - 1, sa_list1->sa, sa_list1->curr_size,
                                                  occ_list1->occ, occ_list1->curr_size );
                            keepSoap3Results ( batchFirstReadId + batchReadId, sa_list2, occ_list2, dpInput, hspaux );
                            addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId, sa_list2->sa, sa_list2->curr_size,
                                                  occ_list2->occ, occ_list2->curr_size );
                        }
                        else if ( previousTotalOccurrences < currentTotalOccurrences )
                        {
                            if ( previousTotalOccurrences <= maxHitNumForDP )
                            {
                                keepSoap3Results ( batchFirstReadId + batchReadId - 1, sa_list1, occ_list1, dpInput, hspaux );
                                addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId - 1, sa_list1->sa, sa_list1->curr_size,
                                                      occ_list1->occ, occ_list1->curr_size );
                            }
                            else
                            {
                                keepSoap3Results ( batchFirstReadId + batchReadId - 1, sa_list1, occ_list1, dpInputForNewDefault, hspaux );
                                addToReadInputForDP ( dpInputForNewDefault, batchFirstReadId + batchReadId - 1, sa_list1->sa, sa_list1->curr_size,
                                                      occ_list1->occ, occ_list1->curr_size );
                            }

                            // for keeping the results of the second end
                            keepSoap3Results ( batchFirstReadId + batchReadId, sa_list2, occ_list2, otherSoap3Result, hspaux );
                            addToReadInputForDP ( otherSoap3Result, batchFirstReadId + batchReadId, sa_list2->sa, sa_list2->curr_size,
                                                  occ_list2->occ, occ_list2->curr_size );
                        }
                        else
                        {
                            if ( currentTotalOccurrences <= maxHitNumForDP2 )
                            {
                                keepSoap3Results ( batchFirstReadId + batchReadId, sa_list2, occ_list2, dpInput, hspaux );
                                addToReadInputForDP ( dpInput, batchFirstReadId + batchReadId, sa_list2->sa, sa_list2->curr_size,
                                                      occ_list2->occ, occ_list2->curr_size );
                            }
                            else
                            {
                                keepSoap3Results ( batchFirstReadId + batchReadId, sa_list2, occ_list2, dpInputForNewDefault, hspaux );
                                addToReadInputForDP ( dpInputForNewDefault, batchFirstReadId + batchReadId, sa_list2->sa, sa_list2->curr_size,
                                                      occ_list2->occ, occ_list2->curr_size );
                            }

                            // for keeping the results of the first end
                            keepSoap3Results ( batchFirstReadId + batchReadId - 1, sa_list1, occ_list1, otherSoap3Result, hspaux );
                            addToReadInputForDP ( otherSoap3Result, batchFirstReadId + batchReadId - 1, sa_list1->sa, sa_list1->curr_size,
                                                  occ_list1->occ, occ_list1->curr_size );
                        }
                    }
                }

                if ( needOutputUnalignedReads && ( !needProceedDP ) )
                {
                    // for SAM format
                    // output the pairs of reads which is not properly aligned
                    // either half-aligned or both-unaligned
                    unproperlypairOutputSAMAPI ( &qInput_Positive, occ_list1, occ_list2,
                                                 preQuery, thisQuery,
                                                 preQualities, thisQualities,
                                                 preReadLength, readLength, preQueryName, thisQueryName,
                                                 charArray, qSetting->MaxOutputPerRead, reportType );
                }
            }
        }
    }

    // printf("[hostKernel] # of bad reads handled by CPU : %i\n", numReadsHandledByCPU1);
    // printf("[hostKernel] # of super bad reads handled by CPU : %i\n", numReadsHandledByCPU2);
    free ( badReadPos );
#ifdef BGS_OCC_RESULT_BREAKDOWN
    hostRoundAnswer = numOfAnswer - firstRoundAnswer - secondRoundAnswer;
    printf ( "[hostKernel] GPU-Round-1 reported %u occurrences.\n", firstRoundAnswer );
    printf ( "[hostKernel] GPU-Round-2 reported %u occurrences.\n", secondRoundAnswer );
    printf ( "[hostKernel] CPU         reported %u occurrences.\n", hostRoundAnswer );
#endif
    // printf("Child reported %u occurrences.\n",numOfAnswer);

    if ( outputFileName != NULL )
    {
        switch ( outputFormat )
        {
            case SRA_OUTPUT_FORMAT_SAM_API:
                break;

            case SRA_OUTPUT_FORMAT_DEFAULT:
            case SRA_OUTPUT_FORMAT_PLAIN:
                OCCFlushCache ( &qInput_Positive );

            default:
                fclose ( outputFile );
        }
    }

    // printf("numOfAnswer = %u; numOfAlignedRead = %u \n", numOfAnswer, numOfAlignedRead);
    ( *alignedOcc ) = numOfAnswer;
    ( *alignedReads ) = numOfAlignedRead;
    ( *unAlignedOcc ) = numOfUnAligned;
    SAListFree ( sa_list1 );
    SAListFree ( sa_list2 );
    OCCListFree ( occ_list1 );
    OCCListFree ( occ_list2 );
    PEInputFree ( pe_in );
    PEOutputFree ( pe_out );

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    { DynamicUint8ArrayFree ( charArray ); }

#ifdef BGS_ROUND_BREAKDOWN_TIME
    printf ( "[threadID: %i] total time elapsed: %9.4f seconds\n", threadID, getElapsedTime ( start_time ) );
#endif
    return numOfAnswer;
}






//Function hostKernelSingle : core function of output report and host alignment.
// This function is designed for single-end alignment
//
//batchFirstReadId: is the readId of the first read in the batch
//                  w.r.t to the read input file. (same among all threads)
//skipFirst       : (Multi-threading) define the number of reads in the batch to skip for 1st batch-run.
//answers         : is the answer arrays used by GPU to store alignment result of the 1st batch-run.
//badAnswers      : is the answer arrays used by GPU to store alignment result of the 2st batch-run.
//badReadIndices  : is the array of the read id w.r.t to the batch.

inline unsigned long long hostKernelSingle ( char * upkdQueryNames, unsigned int batchFirstReadId,
        unsigned int skipFirst, unsigned int numQueries,
        unsigned int * queries, unsigned int word_per_query,
        SRASetting * qSetting, SRAIndex * aIndex, SRAQueryResultCount * rOutput,
        SRAModel ** SRAMismatchModel, SRAModel ** SRAMismatchModel_neg,
        unsigned int ** answers, unsigned long long * alignedOcc, unsigned int * alignedReads,
        unsigned int ** badReadIndices, unsigned int ** badAnswers,
        unsigned int * badStartOffset, unsigned int * badCountOffset,
        unsigned int maxReadLength, unsigned int * readLengths, unsigned int numCases,
        unsigned int sa_range_allowed_1, unsigned int sa_range_allowed_2,
        char skip_round_2, unsigned int word_per_ans,
        SingleAlgnResult * alignResult, unsigned int maxHitNum,
        unsigned int * readIDs )
{
    ullint numOfAnswer = 0;
    ullint preNumOfAnswer = 0;
    uint numOfAlignedRead = 0;
    unsigned char pStrandQuery[2][MAX_READ_LENGTH];
#ifndef BGS_DISABLE_NEGATIVE_STRAND
    unsigned char oStrandQuery[MAX_READ_LENGTH];
#endif
    uint * ans;
    unsigned int l, r;
    int strand;
    unsigned char * thisQuery = NULL;
    char * thisQueryName = NULL;
    int i, j;
    ullint readId = 0;
    ullint batchReadId;
    int readLength = MAX_READ_LENGTH;
    ullint preSaNum = 0;
    ullint preOccNum = 0;
    int preFlag = 0;
    BWT * bwt = aIndex->bwt;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
#ifdef BGS_ROUND_BREAKDOWN_TIME
    double start_time = setStartTime ();
#endif
    int word_per_ans_2 = sa_range_allowed_2 * 2;
    uint * badReadPos = ( uint * ) malloc ( numCases * sizeof ( uint ) );

    char dummyQuality[MAX_READ_LENGTH];
    memset ( dummyQuality, 1, sizeof ( char ) *MAX_READ_LENGTH );
    SRAQueryInput qInput_Positive;
    SRAQueryInfo qInfo_Positive;
    qInfo_Positive.ReadStrand = QUERY_POS_STRAND;
    qInfo_Positive.ReadQuality = dummyQuality;
    qInput_Positive.QueryInfo = &qInfo_Positive;
    qInput_Positive.QuerySetting = qSetting;
    qInput_Positive.AlgnmtIndex = aIndex;
    qInput_Positive.QueryOutput = rOutput;

    SRAQueryInput qInput_Negative;
    SRAQueryInfo qInfo_Negative;
    qInfo_Negative.ReadStrand = QUERY_NEG_STRAND;
    qInfo_Negative.ReadQuality = dummyQuality;
    qInput_Negative.QueryInfo = &qInfo_Negative;
    qInput_Negative.QuerySetting = qSetting;
    qInput_Negative.AlgnmtIndex = aIndex;
    qInput_Negative.QueryOutput = rOutput;

    /////////////////////////////////
    // Allocate and Fill CharMap
    unsigned char charMap[256];

    INDEXFillCharMap ( charMap );

    for ( j = 0; j < numCases; j++ )
    { badReadPos[j] = 0; }

    for ( j = 0; j < numQueries; ++j )
    {
        rOutput->TotalOccurrences = 0;
        // OBTAIN THE ALIGNMENT RESULT
        batchReadId = skipFirst + j;

        if ( readIDs == NULL )
        { readId = batchFirstReadId + batchReadId; }
        else
        { readId = readIDs[batchFirstReadId + batchReadId] - 1; }

        thisQuery = pStrandQuery[preFlag];
        retrieve2BWTQueryFormat ( queries, batchFirstReadId + skipFirst + j, word_per_query, thisQuery );
        preFlag++;
        preFlag = preFlag % 2;
        thisQueryName = upkdQueryNames + hspaux->maxLenReadName * readId;
        readLength = readLengths[batchFirstReadId + batchReadId];
        qInfo_Positive.ReadName = thisQueryName;
        qInfo_Negative.ReadName = thisQueryName;
        qInfo_Positive.ReadId = readId + 1;
        qInfo_Negative.ReadId = readId + 1;
        qInfo_Positive.ReadCode = thisQuery;
        qInfo_Positive.ReadLength = readLength;
        qInfo_Negative.ReadLength = readLength;
        preSaNum = alignResult->saTotalNum;
        preOccNum = alignResult->occTotalNum;

        for ( int whichCase = 0; whichCase < numCases; whichCase++ )
        {
            ans = answers[whichCase] + ( ( batchReadId ) / 32 * 32 * word_per_ans + ( batchReadId ) % 32 );

            if ( * ( ans ) < 0xFFFFFFFD )
            {
                //Number of SA ranges reported within MAX_SA_RANGES_ALLOWED
                for ( i = 0; i < sa_range_allowed_1; i++ )
                {
                    if ( ( * ( ans + ( i * 2 ) * 32 ) ) >= 0xFFFFFFFD || ( * ( ans + ( i * 2 + 1 ) * 32 ) ) >= 0xFFFFFFFD )
                    {
                        break;
                    }

                    l = * ( ans + ( i * 2 ) * 32 );
                    r = l + ( * ( ans + ( i * 2 + 1 ) * 32 ) &BGS_GPU_ANSWER_OFFSET_MASK );
                    strand = ( ( ( * ( ans + ( i * 2 + 1 ) * 32 ) >> ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) & 1 ) ) + 1;

                    // num_mis = (((*(ans+(i*2+1)*32)>>(BGS_GPU_ANSWER_OFFSET_LENGTH)) & 7));
                    // printf("[inside] ReadID: %u; SA range: [%u - %u]; \n", readId, l, r);
                    if ( l <= r && r <= bwt->textLength ) // assume unique hits
                    {
                        addSAToAlgnResult ( alignResult, l, r, strand );
                        // printf("%u %u %u %i\n", readId, l, r, (int) strand);
                        numOfAnswer += ( r - l + 1 );
                    }
                    else
                    {
                        break;
                    }
                }
            }
            else if ( ( skip_round_2 == 1 ) && ( * ( ans ) > 0xFFFFFFFD ) )
            {
                // for the case when (1) round 2 is skipped and (2) it is a bad read
                // this read will be handled by CPU
#ifndef BGS_DISABLE_NEGATIVE_STRAND
                for ( i = 0; i < readLength; i++ )
                {
                    oStrandQuery[i] = soap3DnaComplement[thisQuery[readLength - i - 1]];
                }

                qInfo_Negative.ReadCode = oStrandQuery;
                //Reversing the read from positive strand to negative strand.
                numOfAnswer += ProcessReadDoubleStrand3 ( &qInput_Positive, &qInput_Negative, SRAMismatchModel[readLength], SRAMismatchModel_neg[readLength], whichCase, alignResult );
#endif
#ifdef BGS_DISABLE_NEGATIVE_STRAND
                numOfAnswer += ProcessReadSingleStrand3 ( &qInput_Positive, SRAMismatchModel[readLength], whichCase, alignResult );
#endif
                // numReadsHandledByCPU1++;
            }
            else if ( ( skip_round_2 == 0 ) && ( * ( ans ) > 0xFFFFFFFD ) )
            {
                // check the second round answer
                bool isSuperBadRead = true;

                while ( badReadPos[whichCase] < badCountOffset[whichCase] && badReadIndices[whichCase][badReadPos[whichCase] + badStartOffset[whichCase]] <= batchReadId )
                {
                    if ( badReadIndices[whichCase][badReadPos[whichCase] + badStartOffset[whichCase]] == batchReadId )
                    {
                        ullint pos = badReadPos[whichCase];
                        ans = badAnswers[whichCase] + ( ( pos + badStartOffset[whichCase] ) / 32 * 32 * word_per_ans_2 + ( pos + badStartOffset[whichCase] ) % 32 );

                        if ( * ( ans ) <= 0xFFFFFFFD )
                        {
                            // either no hit or number of SA ranges reported within sa_range_allowed_2
                            isSuperBadRead = false;

                            if ( * ( ans ) < 0xFFFFFFFD )
                            {
                                //Number of SA ranges reported within sa_range_allowed_2
                                for ( i = 0; i < sa_range_allowed_2; i++ )
                                {
                                    if ( ( * ( ans + ( i * 2 ) * 32 ) ) >= 0xFFFFFFFD || ( * ( ans + ( i * 2 + 1 ) * 32 ) ) >= 0xFFFFFFFD )
                                    {
                                        break;
                                    }

                                    l = * ( ans + ( i * 2 ) * 32 );
                                    r = l + ( * ( ans + ( i * 2 + 1 ) * 32 ) &BGS_GPU_ANSWER_OFFSET_MASK );
                                    strand = ( ( * ( ans + ( i * 2 + 1 ) * 32 ) >> ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) & 1 ) + 1;

                                    // num_mis = ((*(ans+(i*2+1)*32)>>(BGS_GPU_ANSWER_OFFSET_LENGTH)) & 7);
                                    if ( l <= r && r <= bwt->textLength ) // assume unique hits
                                    {
                                        // [attention]: needs to update the corresponding strand and number of mismatches
                                        addSAToAlgnResult ( alignResult, l, r, strand );
                                        // printf("%u %u %u %i\n", readId, l, r, (int) strand);
                                        numOfAnswer += ( r - l + 1 );
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    badReadPos[whichCase]++;
                }

                if ( isSuperBadRead )
                {
                    // super bad read
#ifndef BGS_DISABLE_NEGATIVE_STRAND
                    //Reversing the read from positive strand to negative strand.
                    for ( i = 0; i < readLength; i++ )
                    {
                        oStrandQuery[i] = soap3DnaComplement[thisQuery[readLength - i - 1]];
                    }

                    qInfo_Negative.ReadCode = oStrandQuery;
                    numOfAnswer += ProcessReadDoubleStrand3 ( &qInput_Positive, &qInput_Negative, SRAMismatchModel[readLength], SRAMismatchModel_neg[readLength], whichCase, alignResult );
#endif
#ifdef BGS_DISABLE_NEGATIVE_STRAND
                    numOfAnswer += ProcessReadSingleStrand3 ( &qInput_Positive, SRAMismatchModel[readLength], whichCase, alignResult );
#endif
                }
            }
        }

        if ( hspaux->ProceedDPForTooManyHits == 0 )
        {
            // if too many hits, then skip the seed
            if ( numOfAnswer == preNumOfAnswer )
            {
                // no hit
                addReadInfoToAlgnResult ( alignResult, readId, 0 );
            }
            else if ( numOfAnswer - preNumOfAnswer <= maxHitNum )
            {
                // add read info to alignResult
                unsigned char status = 0;

                if ( preSaNum < alignResult->saTotalNum )
                { status = ( status | 1 ); }

                if ( preOccNum < alignResult->occTotalNum )
                { status = ( status | 2 ); }

                preSaNum = alignResult->saTotalNum;
                preOccNum = alignResult->occTotalNum;
                addReadInfoToAlgnResult ( alignResult, readId, status );
                numOfAlignedRead++;
                preNumOfAnswer = numOfAnswer;
            }
            else
            {
                // too many hits
                // remove all the answers
                alignResult->saTotalNum = preSaNum;
                alignResult->occTotalNum = preOccNum;
                numOfAnswer = preNumOfAnswer;
                addReadInfoToAlgnResult ( alignResult, readId, 4 );
            }
        }
        else
        {
            // if too many hits
            // then only keep "maxHitHum" of hits
            if ( numOfAnswer == preNumOfAnswer )
            {
                // no hit
                addReadInfoToAlgnResult ( alignResult, readId, 0 );
            }
            else
            {
                // have hit
                if ( numOfAnswer - preNumOfAnswer > maxHitNum )
                {
                    // too many hits
                    // then only keep "maxHitHum" of hits
                    if ( alignResult->occTotalNum > preOccNum + maxHitNum )
                    {
                        // only keep the first "maxHitHum" of occurrences
                        alignResult->occTotalNum = preOccNum + maxHitNum;
                        // remove all the SA answers
                        alignResult->saTotalNum = preSaNum;
                    }
                    else
                    {
                        uint currNumOcc = alignResult->occTotalNum - preOccNum;

                        for ( i = preSaNum; i < alignResult->saTotalNum && currNumOcc < maxHitNum; i++ )
                        {
                            unsigned int k, l, r;
                            int strand;
                            l = alignResult->sa_list[i].saLeft;
                            r = alignResult->sa_list[i].saRight;
                            strand = alignResult->sa_list[i].strand;

                            for ( k = l; k <= r; k++ )
                            {
                                addOccToAlgnResult ( alignResult, ( *bwt->_bwtSaValue ) ( bwt, k ), strand );
                            }

                            currNumOcc += ( r - l + 1 );
                        }

                        // only keep the first "maxHitHum" of occurrences
                        if ( alignResult->occTotalNum > preOccNum + maxHitNum )
                        {
                            alignResult->occTotalNum = preOccNum + maxHitNum;
                        }

                        // remove all the SA answers
                        alignResult->saTotalNum = preSaNum;
                    }
                }

                // add read info to alignResult
                unsigned char status = 0;

                if ( preSaNum < alignResult->saTotalNum )
                { status = ( status | 1 ); }

                if ( preOccNum < alignResult->occTotalNum )
                { status = ( status | 2 ); }

                preSaNum = alignResult->saTotalNum;
                preOccNum = alignResult->occTotalNum;
                addReadInfoToAlgnResult ( alignResult, readId, status );
                numOfAlignedRead++;
                preNumOfAnswer = numOfAnswer;
            }
        }

        /*

          if (numOfAnswer==preNumOfAnswer) {
                // no hit
                addReadInfoToAlgnResult(alignResult, readId, 0);
          } else {
              // have hit

              if (numOfAnswer-preNumOfAnswer > maxHitNum) {
                  // too many hits
                  // then only keep "maxHitHum" of hits
                  if (alignResult->occTotalNum > preOccNum + maxHitNum) {
                      // only keep the first "maxHitHum" of occurrences
                      alignResult->occTotalNum = preOccNum + maxHitNum;
                      // remove all the SA answers
                      alignResult->saTotalNum = preSaNum;
                  } else {
                      uint currNumOcc = alignResult->occTotalNum - preOccNum;
                      for (i=preSaNum; i<alignResult->saTotalNum && currNumOcc<maxHitNum; i++) {
                          unsigned int k,l,r;
                          int strand;
                          l = alignResult->sa_list[i].saLeft;
                          r = alignResult->sa_list[i].saRight;
                          strand = alignResult->sa_list[i].strand;

                          for (k=l; k<=r; k++) {
                              addOccToAlgnResult(alignResult, (*bwt->_bwtSaValue)(bwt,k), strand);
                          }
                          currNumOcc+=(r-l+1);
                      }
                      // only keep the first "maxHitHum" of occurrences
                      if (alignResult->occTotalNum > preOccNum + maxHitNum) {
                          alignResult->occTotalNum = preOccNum + maxHitNum;
                      }
                      // remove all the SA answers
                      alignResult->saTotalNum = preSaNum;
                  }
              }

              // add read info to alignResult
              unsigned char status = 0;
              if (preSaNum < alignResult->saTotalNum)
                  status = (status | 1);
              if (preOccNum < alignResult->occTotalNum)
                  status = (status | 2);

              preSaNum = alignResult->saTotalNum;
              preOccNum = alignResult->occTotalNum;
              addReadInfoToAlgnResult(alignResult, readId, status);

              numOfAlignedRead++;
              preNumOfAnswer = numOfAnswer;
          }
         */
    }

    free ( badReadPos );
    // OCCFlushCache(&qInput_Positive);
    // printf("numOfAnswer = %u; numOfAlignedRead = %u \n", numOfAnswer, numOfAlignedRead);
    ( *alignedOcc ) = numOfAnswer;
    ( *alignedReads ) = numOfAlignedRead;
    return numOfAnswer;
}


void bwase_initialize ( int * g_log_n )
{
    int i;

    for ( i = 1; i != 256; ++i ) { g_log_n[i] = ( int ) ( 4.343 * log ( i ) + 0.5 ); }
}

// A Thread Wrapper to perform the CPU task
void * hostKernelThreadWrapper ( void * arg )
{
    HostKernelArguements * myArgs = ( HostKernelArguements * ) arg;
#ifdef BGS_CPU_CASE_BREAKDOWN_TIME
    double startTime = setStartTime ();
#endif
    hostKernel ( myArgs->upkdQualities, myArgs->upkdQueryNames, myArgs->batchFirstReadId, myArgs->skipFirst, myArgs->numQueries,
                 myArgs->queries, myArgs->word_per_query,
                 & ( myArgs->sraQuerySettings ),
                 & ( myArgs->sraIndex ),
                 & ( myArgs->sraAccumulatedResultCount ),
                 myArgs->SRAMismatchModel, myArgs->SRAMismatchModel_neg,
                 myArgs->outputFileName,
                 myArgs->answers,
                 & ( myArgs->alignedOcc ), & ( myArgs->alignedReads ),
                 myArgs->badReadIndices, myArgs->badAnswers, myArgs->badStartOffset, myArgs->badCountOffset,
                 myArgs->maxReadLength, myArgs->readLengths, myArgs->seedLengths, myArgs->numCases,
                 myArgs->sa_range_allowed_1, myArgs->sa_range_allowed_2,
                 myArgs->readIDs, myArgs->skip_round_2, myArgs->accumReadNum,
                 myArgs->word_per_ans, myArgs->readType,
                 myArgs->insert_low, myArgs->insert_high,
                 myArgs->peStrandLeftLeg, myArgs->peStrandRightLeg,
                 myArgs->threadId, myArgs->unAlignedIDs, & ( myArgs->unAlignedOcc ),
                 myArgs->reportType, myArgs->peMaxOutputPerPair, myArgs->isTerminalCase,
                 myArgs->readInput, myArgs->readInputForNewDefault,
                 myArgs->otherSoap3Result,
                 myArgs->bothUnalignedPairs,
                 myArgs->maxHitNum, myArgs->maxHitNum2, myArgs->enableDP,
                 myArgs->SRAMismatchModel2, myArgs->SRAMismatchModel2_neg, myArgs->numCases2 );
#ifdef BGS_CPU_CASE_BREAKDOWN_TIME
    double alignmentTime = getElapsedTime ( startTime );
    printf ( "[Main:Thread_%d] Elapsed time on host : %9.4f seconds\n", myArgs->threadId, alignmentTime );
#endif
    pthread_exit ( 0 );
    return 0;
}

// A Thread Wrapper to perform the CPU task
void * hostKernelThreadWrapperSingle ( void * arg )
{
    HostKernelArguements * myArgs = ( HostKernelArguements * ) arg;
#ifdef BGS_CPU_CASE_BREAKDOWN_TIME
    double startTime = setStartTime ();
#endif
    hostKernelSingle ( myArgs->upkdQueryNames, myArgs->batchFirstReadId, myArgs->skipFirst, myArgs->numQueries,
                       myArgs->queries, myArgs->word_per_query,
                       & ( myArgs->sraQuerySettings ),
                       & ( myArgs->sraIndex ),
                       & ( myArgs->sraAccumulatedResultCount ),
                       myArgs->SRAMismatchModel, myArgs->SRAMismatchModel_neg,
                       myArgs->answers,
                       & ( myArgs->alignedOcc ), & ( myArgs->alignedReads ),
                       myArgs->badReadIndices, myArgs->badAnswers, myArgs->badStartOffset, myArgs->badCountOffset,
                       myArgs->maxReadLength, myArgs->readLengths, myArgs->numCases,
                       myArgs->sa_range_allowed_1, myArgs->sa_range_allowed_2,
                       myArgs->skip_round_2, myArgs->word_per_ans,
                       myArgs->alignResult, myArgs->maxHitNum, myArgs->readIDs );
#ifdef BGS_CPU_CASE_BREAKDOWN_TIME
    double alignmentTime = getElapsedTime ( startTime );
    printf ( "[Main:Thread_%d] Elapsed time on host : %9.4f seconds\n", myArgs->threadId, alignmentTime );
#endif
    pthread_exit ( 0 );
    return 0;
}


// set the multi-threading arguments
// for 1-mismatch single-end alignment
void setHostKernelArguments2 ( HostKernelArguements * hostKernelArguments, pthread_t * threads, Soap3Index * index,
                               uint maxReadLength, uint word_per_query, uint word_per_ans,
                               uint cpuNumThreads )
{

    // MULTI-THREADING arguments
    for ( uint threadId = 0; threadId < cpuNumThreads; ++threadId )
    {
        threads[threadId] = 0;
        hostKernelArguments[threadId].occ = OCCConstruct ();

        for ( uint i = 0 ; i <= MAX_NUM_OF_MISMATCH ; i++ )
        {
            hostKernelArguments[threadId].sraAccumulatedResultCount.WithError[i] = 0;
        }

        // Copy SRAIndex
        memcpy ( & ( hostKernelArguments[threadId].sraIndex ), index->sraIndex, sizeof ( SRAIndex ) );

        hostKernelArguments[threadId].sraQuerySettings.occ = hostKernelArguments[threadId].occ;
        hostKernelArguments[threadId].sraQuerySettings.OutFilePtr = NULL;
        hostKernelArguments[threadId].sraQuerySettings.OutFileName = NULL;
        hostKernelArguments[threadId].sraQuerySettings.OutFilePart = 0;
        hostKernelArguments[threadId].sraQuerySettings.MaxOutputPerRead = 0xFFFFFFFF;
        hostKernelArguments[threadId].sraQuerySettings.writtenCacheCount = 0;
        hostKernelArguments[threadId].sraQuerySettings.ReadStrand = QUERY_BOTH_STRAND;
        hostKernelArguments[threadId].sraQuerySettings.ErrorType = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        hostKernelArguments[threadId].sraQuerySettings.OutputType = SRA_REPORT_ALL;
        hostKernelArguments[threadId].sraQuerySettings.SAMOutFilePtr = NULL;
        hostKernelArguments[threadId].sraAccumulatedResultCount.TotalOccurrences = 0;
        hostKernelArguments[threadId].sraAccumulatedResultCount.RetrievedBySa = 0;
        hostKernelArguments[threadId].sraAccumulatedResultCount.RetrievedByCe = 0;
        hostKernelArguments[threadId].sraAccumulatedResultCount.RetrievedByHocc = 0;
        hostKernelArguments[threadId].batchFirstReadId = 0;
        hostKernelArguments[threadId].skipFirst = 0;
        hostKernelArguments[threadId].numQueries = 0;
        hostKernelArguments[threadId].queries = NULL;
        hostKernelArguments[threadId].word_per_query = word_per_query;
        hostKernelArguments[threadId].answers = NULL;
        hostKernelArguments[threadId].badReadIndices = NULL;
        hostKernelArguments[threadId].badAnswers = NULL;
        hostKernelArguments[threadId].badStartOffset = NULL;
        hostKernelArguments[threadId].badCountOffset = NULL;
        hostKernelArguments[threadId].alignedOcc = 0;
        hostKernelArguments[threadId].alignedReads = 0;
        hostKernelArguments[threadId].threadId = threadId;
        hostKernelArguments[threadId].maxReadLength = maxReadLength;
        hostKernelArguments[threadId].word_per_ans = word_per_ans;
    }
}

