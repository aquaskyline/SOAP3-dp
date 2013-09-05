/*
 *
 *    OutputDPResult.c
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

#include "OutputDPResult.h"


void getReadInfo ( SRAQueryInput * qInput, unsigned int readID, unsigned int displayID,
                   unsigned char * upkdQueries, char * upkdQueryNames, unsigned int * upkdReadLengths,
                   unsigned int maxReadLength, int maxReadNameLength )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    qInfo->ReadName = upkdQueryNames + ( unsigned long long ) readID * maxReadNameLength;
    qInfo->ReadCode = upkdQueries + ( unsigned long long ) readID * maxReadLength;
    qInfo->ReadLength = upkdReadLengths[readID];
    qInfo->ReadId = displayID;
}

void getReadInfoPackedSeq ( SRAQueryInput * qInput, unsigned int readID, unsigned int displayID,
                            unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char * upkdQueryNames,
                            unsigned int word_per_query, unsigned char * unpackedQueries, char * upkdQualities,
                            unsigned int maxReadLength, int maxReadNameLength )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    qInfo->ReadName = upkdQueryNames + ( unsigned long long ) ( upkdReadIDs[readID] - 1 ) * maxReadNameLength;
    retrieve2BWTQueryFormat ( queries, readID, word_per_query, unpackedQueries );
    qInfo->ReadCode = unpackedQueries;
    qInfo->ReadLength = readLengths[readID];
    qInfo->ReadId = displayID;
    qInfo->ReadQuality = upkdQualities + ( unsigned long long ) ( upkdReadIDs[readID] - 1 ) * maxReadLength;
}

void getAlgnInfo ( PEPairs & pePair, AlgnmtDPResult & dpResult, int peStrandLeftLeg, int peStrandRightLeg )
{
    pePair.algnmt_1 = dpResult.algnmt_1;
    pePair.strand_1 = dpResult.strand_1;
    pePair.mismatch_1 = dpResult.score_1;
    pePair.algnmt_2 = dpResult.algnmt_2;
    pePair.strand_2 = dpResult.strand_2;
    pePair.mismatch_2 = dpResult.score_2;
}

void getAlgnInfo2 ( PEPairs & pePair, DeepDPAlignResult & dpResult, int peStrandLeftLeg, int peStrandRightLeg )
{
    pePair.algnmt_1 = dpResult.algnmt_1;
    pePair.strand_1 = dpResult.strand_1;
    pePair.mismatch_1 = dpResult.score_1;
    pePair.algnmt_2 = dpResult.algnmt_2;
    pePair.strand_2 = dpResult.strand_2;
    pePair.mismatch_2 = dpResult.score_2;
}

// output a set of default DP alignment results
void outputRead ( AlgnmtDPResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char * upkdQueryNames,
                  unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                  FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg )
{
    // output a set of alignment results
    // the alignments of the same pair of reads have to be together
    // algnNum: the number of results in "algnResult"
    if ( algnNum == 0 )
    { return; }

    HSPAux * hspaux = index->sraIndex->hspaux;

    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy ( & ( aIndex ), index->sraIndex, sizeof ( SRAIndex ) );

    SRASetting qSetting;
    qSetting.OutFilePtr = outputDPFile;
    qSetting.occ = OCCConstruct ();
    qSetting.OutFileFormat = outputFormat;
    qSetting.SAMOutFilePtr = samOutputDPFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.AlgnmtIndex = &aIndex;
    qInput.QueryInfo = &qInfo;
    OCC * occ = qInput.QuerySetting->occ;
    unsigned int currReadID;
    unsigned int preReadID = 0xFFFFFFFF;
    unsigned int displayID;
    PEPairs pePair;

    for ( unsigned int i = 0; i < algnNum; i++ )
    {
        // check the readID
        /*
        if (algnResult[i].readID%2 > 0) {
              fprintf(stderr, "ERROR! The read ID is not even! i=%u readID=%u\n", i, algnResult[i].readID);
        }*/
        int isFirstUnaligned = ( algnResult[i].algnmt_1 == 0xFFFFFFFF );
        int isSecondUnaligned = ( algnResult[i].algnmt_2 == 0xFFFFFFFF );
        int properlyPaired = ( ( !isFirstUnaligned ) && ( !isSecondUnaligned ) );
        unsigned int readLen1 = upkdReadLengths[algnResult[i].readID]; // length of the first read
        unsigned int readLen2 = upkdReadLengths[algnResult[i].readID + 1]; // length of the second read

        if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT ) && ( algnResult[i].readID != preReadID )
                && properlyPaired )
        {
            // output the information of the (second) read
            currReadID = algnResult[i].readID + 1;
            displayID = ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN ) ? ( currReadID + accumReadNum ) / 2 + 1 : currReadID + accumReadNum + 1;
            getReadInfo ( &qInput, currReadID, displayID, upkdQueries, upkdQueryNames, upkdReadLengths, maxReadLength, hspaux->maxLenReadName );
            OCCReportDelimitorDP ( &qInput );
            preReadID = algnResult[i].readID;
        }

        if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT )
                && properlyPaired )
        {
#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
            algnResult[i].strand_2 = 2;
#endif

            // for plain and binary, we only output the answer for properly aligned pairs
            // output the alignment of the first read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_1;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_1;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_1;

            if ( algnResult[i].whichFromDP == 0 )
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
                occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString;
            }
            else
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 0;
            }

            occ->occPositionCache[occ->occPositionCacheCount].len = readLen1;
            occ->occPositionCacheCount++;

            // output the alignment of the second read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_2;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_2;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_2;

            if ( algnResult[i].whichFromDP == 1 )
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
                occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString;
            }
            else
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 0;
            }

            occ->occPositionCache[occ->occPositionCacheCount].len = readLen2;
            occ->occPositionCacheCount++;
        }
        else if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            // SAM format
            getAlgnInfo ( pePair, algnResult[i], peStrandLeftLeg, peStrandRightLeg );
            // output the alignment of the first read
            currReadID = algnResult[i].readID;
            getReadInfo ( &qInput, currReadID, currReadID, upkdQueries, upkdQueryNames, upkdReadLengths, maxReadLength, hspaux->maxLenReadName );

            if ( algnResult[i].whichFromDP == 0 )
            {
                // this should be properly paired
                OCCDirectWritePairOccSAMAPIwCIGAR ( &qInput, &pePair, isFirstUnaligned, isSecondUnaligned, readLen2, algnResult[i].cigarString );
            }
            else
            {
                OCCDirectWritePairOccSAMAPI ( &qInput, &pePair, isFirstUnaligned, isSecondUnaligned, readLen2, properlyPaired );
            }

            // output the alignment of the second read
            currReadID = algnResult[i].readID + 1;
            getReadInfo ( &qInput, currReadID, currReadID, upkdQueries, upkdQueryNames, upkdReadLengths, maxReadLength, hspaux->maxLenReadName );

            if ( algnResult[i].whichFromDP == 1 )
            {
                // this should be properly paired
                OCCDirectWritePairOccSAMAPI2wCIGAR ( &qInput, &pePair, isFirstUnaligned, isSecondUnaligned, readLen1, algnResult[i].cigarString );
            }
            else
            {
                OCCDirectWritePairOccSAMAPI2 ( &qInput, &pePair, isFirstUnaligned, isSecondUnaligned, readLen1, properlyPaired );
            }
        }
    }

    OCCFlushCacheDP ( &qInput );
    OCCFree ( qSetting.occ );
}

void outputRead2 ( AlgnmtDPResult * algnResult, unsigned int algnNum,
                   unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char * upkdQueryNames, char * upkdQualities,
                   unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                   FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg )
{
    if ( algnNum == 0 )
    {
        return;
    }

    HSPAux * hspaux = index->sraIndex->hspaux;

    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy ( & ( aIndex ), index->sraIndex, sizeof ( SRAIndex ) );

    SRASetting qSetting;
    qSetting.OutFilePtr = outputDPFile;
    qSetting.occ = OCCConstruct ();
    qSetting.OutFileFormat = outputFormat;
    qSetting.SAMOutFilePtr = samOutputDPFilePtr;

    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.AlgnmtIndex = &aIndex;
    qInput.QueryInfo = &qInfo;
    OCC * occ = qSetting.occ;
    unsigned int preReadID = 0xFFFFFFFF;
    unsigned int upkdReadID;
    unsigned int displayID;
    AlgnmtDPResult * bestAlgn = NULL;
    AlgnmtDPResult * currAlgn = NULL;
    char minNumMismatch, currNumMismatch;
    int maxScore, currScore;
    unsigned char upkdPosQuery[MAX_READ_LENGTH];
    unsigned char upkdPosQuery2[MAX_READ_LENGTH];
    char * qualities1, *qualities2;
    unsigned int readLen1, readLen2;
    char * queryName1, *queryName2;
    unsigned int word_per_query = getWordPerQuery ( maxReadLength );
    unsigned int startIndex = 0;
    DynamicUint8Array * charArray = NULL;

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    { charArray = DynamicUint8ArrayConstruct (); }

    for ( unsigned int i = 0; i < algnNum; i++ )
    {
        int isFirstUnaligned = ( algnResult[i].algnmt_1 == 0xFFFFFFFF );
        int isSecondUnaligned = ( algnResult[i].algnmt_2 == 0xFFFFFFFF );
        int properlyPaired = ( ( !isFirstUnaligned ) && ( !isSecondUnaligned ) );

        if ( algnResult[i].readID != preReadID )
        {
            if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
            {
                if ( preReadID != 0xFFFFFFFF )
                {
                    pairDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                                         startIndex, i - startIndex,
                                         upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                         readLen1, readLen2, queryName1, queryName2,
                                         charArray );
                }

                // reset the bestAlgn
                bestAlgn = & ( algnResult[i] );

                if ( isFirstUnaligned && !isSecondUnaligned )
                {
                    minNumMismatch = bestAlgn->score_2;
                    maxScore = -127;
                }
                else if ( isSecondUnaligned && !isFirstUnaligned )
                {
                    minNumMismatch = bestAlgn->score_1;
                    maxScore = -127;
                }
                else if ( properlyPaired )
                {
                    minNumMismatch = ( bestAlgn->whichFromDP == 1 ) ? bestAlgn->score_1 : bestAlgn->score_2;
                    maxScore = ( bestAlgn->whichFromDP == 1 ) ? bestAlgn->score_2 : bestAlgn->score_1;
                }
            }

            // get the information of both reads
            readLen1 = readLengths[algnResult[i].readID]; // length of the first read
            readLen2 = readLengths[algnResult[i].readID + 1]; // length of the second read
            // name of the first read
            queryName1 = upkdQueryNames + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID] - 1 ) * hspaux->maxLenReadName;
            // name of the second read
            queryName2 = upkdQueryNames + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID + 1] - 1 ) * hspaux->maxLenReadName;
            // sequence the first read
            retrieve2BWTQueryFormat ( queries, algnResult[i].readID, word_per_query, upkdPosQuery );
            // sequence the second read
            retrieve2BWTQueryFormat ( queries, algnResult[i].readID + 1, word_per_query, upkdPosQuery2 );
            // base qualities of the first read
            qualities1 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID] - 1 ) * maxReadLength;
            // base qualities of the second read
            qualities2 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID + 1] - 1 ) * maxReadLength;
            // set qInput for the second read
            qInfo.ReadName = queryName2;
            qInfo.ReadCode = upkdPosQuery2;
            qInfo.ReadLength = readLen2;
            qInfo.ReadQuality = qualities2;
            upkdReadID = upkdReadIDs[algnResult[i].readID + 1];
            displayID = ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN ) ? ( upkdReadID + accumReadNum ) / 2 : upkdReadID + accumReadNum;
            qInfo.ReadId = displayID;
            startIndex = i;

            if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT ) && properlyPaired )
            {
                OCCReportDelimitorDP ( &qInput );
            }

            preReadID = algnResult[i].readID;
        }
        else if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            // SAM format
            // update the bestAlgn
            if ( isFirstUnaligned && !isSecondUnaligned )
            {
                if ( algnResult[i].score_2 < minNumMismatch )
                {
                    bestAlgn = & ( algnResult[i] );
                    minNumMismatch = bestAlgn->score_2;
                    maxScore = -127;
                }
            }
            else if ( isSecondUnaligned && !isFirstUnaligned )
            {
                if ( algnResult[i].score_1 < minNumMismatch )
                {
                    bestAlgn = & ( algnResult[i] );
                    minNumMismatch = bestAlgn->score_1;
                    maxScore = -127;
                }
            }
            else if ( properlyPaired )
            {
                currAlgn = & ( algnResult[i] );
                currNumMismatch = ( currAlgn->whichFromDP == 1 ) ? currAlgn->score_1 : currAlgn->score_2;
                currScore = ( currAlgn->whichFromDP == 1 ) ? currAlgn->score_2 : currAlgn->score_1;

                if ( ( currNumMismatch == minNumMismatch && currScore > maxScore ) ||
                        ( currNumMismatch < minNumMismatch ) )
                {
                    bestAlgn = currAlgn;
                    minNumMismatch = currNumMismatch;
                    maxScore = currScore;
                }
            }
        }

        if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT )
                && properlyPaired )
        {
#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
            algnResult[i].strand_2 = 2;
#endif

            // for plain and binary, we only output the answer for properly aligned pairs
            // output the alignment of the first read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_1;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_1;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_1;

            if ( algnResult[i].whichFromDP == 0 )
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
                convertToCigarStr2 ( algnResult[i].cigarString );
                occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString;
            }
            else
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 0;
            }

            occ->occPositionCache[occ->occPositionCacheCount].len = readLen1;
            occ->occPositionCacheCount++;

            // output the alignment of the second read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_2;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_2;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_2;

            if ( algnResult[i].whichFromDP == 1 )
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
                convertToCigarStr2 ( algnResult[i].cigarString );
                occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString;
            }
            else
            {
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 0;
            }

            occ->occPositionCache[occ->occPositionCacheCount].len = readLen2;
            occ->occPositionCacheCount++;
        }
    }

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && preReadID != 0xFFFFFFFF )
    {
        pairDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                             startIndex, algnNum - startIndex,
                             upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                             readLen1, readLen2, queryName1, queryName2,
                             charArray );
    }

    if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API )
    { OCCFlushCacheDP ( &qInput ); }

    OCCFree ( qSetting.occ );

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    { DynamicUint8ArrayFree ( charArray ); }
}


// output a set of deep DP alignment results
void outputDeepDPResult ( DeepDPAlignResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char * upkdQueryNames,
                          unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                          FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg )
{
    // output a set of alignment results
    // the alignments of the same pair of reads have to be together
    // algnNum: the number of results in "algnResult"
    if ( algnNum == 0 )
    { return; }

    HSPAux * hspaux = index->sraIndex->hspaux;

    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy ( & ( aIndex ), index->sraIndex, sizeof ( SRAIndex ) );

    SRASetting qSetting;
    qSetting.OutFilePtr = outputDPFile;
    qSetting.occ = OCCConstruct ();
    qSetting.OutFileFormat = outputFormat;
    qSetting.SAMOutFilePtr = samOutputDPFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.QueryInfo = &qInfo;
    qInput.AlgnmtIndex = &aIndex;
    OCC * occ = qSetting.occ;
    unsigned int currReadID;
    unsigned int preReadID = 0xFFFFFFFF;
    unsigned int displayID;
    PEPairs pePair;

    for ( unsigned int i = 0; i < algnNum; i++ )
    {
        // check the readID
        /*
         if (algnResult[i].readID%2 > 0) {
         fprintf(stderr, "ERROR! The read ID is not even! i=%u readID=%u\n", i, algnResult[i].readID);
         }*/
        int isFirstUnaligned = ( algnResult[i].algnmt_1 == 0xFFFFFFFF );
        int isSecondUnaligned = ( algnResult[i].algnmt_2 == 0xFFFFFFFF );
        int properlyPaired = ( ( !isFirstUnaligned ) & ( !isSecondUnaligned ) );
        unsigned int readLen1 = upkdReadLengths[algnResult[i].readID]; // length of the first read
        unsigned int readLen2 = upkdReadLengths[algnResult[i].readID + 1]; // length of the second read

        if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT ) && ( algnResult[i].readID != preReadID )
                && properlyPaired )
        {
            // output the information of the (second) read
            currReadID = algnResult[i].readID + 1;
            displayID = ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN ) ? ( currReadID + accumReadNum ) / 2 + 1 : currReadID + accumReadNum + 1;
            getReadInfo ( &qInput, currReadID, displayID, upkdQueries, upkdQueryNames, upkdReadLengths, maxReadLength, hspaux->maxLenReadName );
            OCCReportDelimitorDP ( &qInput );
            preReadID = algnResult[i].readID;
        }

        if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT )
                && properlyPaired )
        {
            // output the alignment of the first read
#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
            algnResult[i].strand_2 = 2;
#endif

            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_1;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_1;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_1;
            occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
            convertToCigarStr2 ( algnResult[i].cigarString_1 );
            occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString_1;
            occ->occPositionCache[occ->occPositionCacheCount].len = readLen1;
            occ->occPositionCacheCount++;

            // output the alignment of the second read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_2;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_2;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_2;
            occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
            convertToCigarStr2 ( algnResult[i].cigarString_2 );
            occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString_2;
            occ->occPositionCache[occ->occPositionCacheCount].len = readLen2;
            occ->occPositionCacheCount++;
        }
        else if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            // SAM format
            getAlgnInfo2 ( pePair, algnResult[i], peStrandLeftLeg, peStrandRightLeg );
            // output the alignment of the first read
            currReadID = algnResult[i].readID;
            getReadInfo ( &qInput, currReadID, currReadID, upkdQueries, upkdQueryNames, upkdReadLengths, maxReadLength, hspaux->maxLenReadName );

            if ( properlyPaired )
            { OCCDirectWritePairOccSAMAPIwCIGAR ( &qInput, &pePair, 0, 0, readLen2, algnResult[i].cigarString_1 ); }
            else
            { OCCDirectWritePairOccSAMAPI ( &qInput, &pePair, isFirstUnaligned, isSecondUnaligned, readLen2, properlyPaired ); }

            // output the alignment of the second read
            currReadID = algnResult[i].readID + 1;
            getReadInfo ( &qInput, currReadID, currReadID, upkdQueries, upkdQueryNames, upkdReadLengths, maxReadLength, hspaux->maxLenReadName );

            if ( properlyPaired )
            { OCCDirectWritePairOccSAMAPI2wCIGAR ( &qInput, &pePair, 0, 0, readLen1, algnResult[i].cigarString_2 ); }
            else
            { OCCDirectWritePairOccSAMAPI2 ( &qInput, &pePair, isFirstUnaligned, isSecondUnaligned, readLen1, properlyPaired ); }
        }
    }

    OCCFlushCacheDP ( &qInput );
    OCCFree ( qSetting.occ );
}

void outputDeepDPResult2 ( DeepDPAlignResult * algnResult, unsigned int algnNum,
                           unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char * upkdQueryNames, char * upkdQualities,
                           unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                           FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg )
{
    if ( algnNum == 0 )
    {
        return;
    }

    BWT * bwt = index->sraIndex->bwt;
    HSPAux * hspaux = index->sraIndex->hspaux;

    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy ( & ( aIndex ), index->sraIndex, sizeof ( SRAIndex ) );

    SRASetting qSetting;
    qSetting.OutFilePtr = outputDPFile;
    qSetting.occ = OCCConstruct ();
    qSetting.OutFileFormat = outputFormat;
    qSetting.SAMOutFilePtr = samOutputDPFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.QueryInfo = &qInfo;
    qInput.AlgnmtIndex = &aIndex;
    OCC * occ = qSetting.occ;
    unsigned int preReadID = 0xFFFFFFFF;
    unsigned int upkdReadID;
    unsigned int displayID;
    DeepDPAlignResult * bestAlgn = NULL;
    int maxScore;
    unsigned char upkdPosQuery[MAX_READ_LENGTH];
    unsigned char upkdPosQuery2[MAX_READ_LENGTH];
    char * qualities1, *qualities2;
    unsigned int readLen1, readLen2;
    char * queryName1, *queryName2;
    unsigned int word_per_query = getWordPerQuery ( maxReadLength );
    unsigned int startIndex = 0;
    DynamicUint8Array * charArray = NULL;

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    { charArray = DynamicUint8ArrayConstruct (); }

    for ( unsigned int i = 0; i < algnNum; i++ )
    {
        int isFirstUnaligned = ( algnResult[i].algnmt_1 == 0xFFFFFFFF );
        int isSecondUnaligned = ( algnResult[i].algnmt_2 == 0xFFFFFFFF );
        int properlyPaired = ( ( !isFirstUnaligned ) && ( !isSecondUnaligned ) );

        if ( algnResult[i].readID != preReadID )
        {
            if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
            {
                if ( preReadID != 0xFFFFFFFF )
                {
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS

                    // ==========================================================
                    // for those not properly paired, needs to perform deep-dp for unaligned reads,
                    // ==========================================================
                    if ( bestAlgn->algnmt_1 == 0xFFFFFFFF && bestAlgn->algnmt_2 == 0xFFFFFFFF && ( ( AllHits * ) hspaux->allHits )->readNum < MAX_UNALIGN_READS_NUM_FOR_DEEP_DP )
                    {
                        int readID = preReadID;
                        int first_read_no_hit  = ( ( hspaux->sa_num[readID] == 0 ) && ( hspaux->occ_num[readID] == 0 ) );
                        int second_read_no_hit = ( ( hspaux->sa_num[readID + 1] == 0 ) && ( hspaux->occ_num[readID + 1] == 0 ) );

                        if ( first_read_no_hit )
                        {
                            addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID );
                        }
                        else
                        {
                            inputSoap3AnsToArray ( ( AllHits * ) hspaux->allHits, readID, hspaux, bwt, readLen1 );
                        }

                        if ( second_read_no_hit || ( !first_read_no_hit ) )
                        {
                            addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID + 1 );
                        }
                        else
                        {
                            inputSoap3AnsToArray ( ( AllHits * ) hspaux->allHits, readID + 1, hspaux, bwt, readLen2 );
                        }
                    }
                    else
                    {
                        pairDeepDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                                                 startIndex, i - startIndex,
                                                 upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                                 readLen1, readLen2, queryName1, queryName2,
                                                 charArray );
                    }

                    // ==========================================================
#else
                    pairDeepDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                                             startIndex, i - startIndex,
                                             upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                             readLen1, readLen2, queryName1, queryName2,
                                             charArray );
#endif
                }

                // reset the bestAlgn
                bestAlgn = & ( algnResult[i] );

                if ( isFirstUnaligned && ( !isSecondUnaligned ) )
                {
                    maxScore = bestAlgn->score_2;
                }
                else if ( isSecondUnaligned && ( !isFirstUnaligned ) )
                {
                    maxScore = bestAlgn->score_1;
                }
                else if ( properlyPaired )
                {
                    maxScore = bestAlgn->score_1 + bestAlgn->score_2;
                }
                else
                {
                    maxScore = -127;
                }
            }

            // get the information of both reads
            readLen1 = readLengths[algnResult[i].readID]; // length of the first read
            readLen2 = readLengths[algnResult[i].readID + 1]; // length of the second read
            // name of the first read
            queryName1 = upkdQueryNames + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID] - 1 ) * hspaux->maxLenReadName;
            // name of the second read
            queryName2 = upkdQueryNames + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID + 1] - 1 ) * hspaux->maxLenReadName;
            // sequence the first read
            retrieve2BWTQueryFormat ( queries, algnResult[i].readID, word_per_query, upkdPosQuery );
            // sequence the second read
            retrieve2BWTQueryFormat ( queries, algnResult[i].readID + 1, word_per_query, upkdPosQuery2 );
            // base qualities of the first read
            qualities1 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID] - 1 ) * maxReadLength;
            // base qualities of the second read
            qualities2 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[algnResult[i].readID + 1] - 1 ) * maxReadLength;
            // set qInput for the second read
            qInfo.ReadName = queryName2;
            qInfo.ReadCode = upkdPosQuery2;
            qInfo.ReadLength = readLen2;
            qInfo.ReadQuality = qualities2;
            upkdReadID = upkdReadIDs[algnResult[i].readID + 1];
            displayID = ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN ) ? ( upkdReadID + accumReadNum ) / 2 : upkdReadID + accumReadNum;
            qInfo.ReadId = displayID;
            startIndex = i;

            if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT ) && properlyPaired )
            {
                OCCReportDelimitorDP ( &qInput );
            }

            preReadID = algnResult[i].readID;
        }
        else if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
        {
            // SAM format
            // update the bestAlgn
            if ( isFirstUnaligned && !isSecondUnaligned )
            {
                if ( algnResult[i].score_2 > maxScore )
                {
                    bestAlgn = & ( algnResult[i] );
                    maxScore = algnResult[i].score_2;
                }
            }
            else if ( isSecondUnaligned && !isFirstUnaligned )
            {
                if ( algnResult[i].score_1 > maxScore )
                {
                    bestAlgn = & ( algnResult[i] );
                    maxScore = algnResult[i].score_1;
                }
            }
            else if ( properlyPaired )
            {
                if ( algnResult[i].score_1 + algnResult[i].score_2 > maxScore )
                {
                    bestAlgn = & ( algnResult[i] );
                    maxScore = algnResult[i].score_1 + algnResult[i].score_2;
                }
            }
        }

        if ( ( outputFormat == SRA_OUTPUT_FORMAT_PLAIN || outputFormat == SRA_OUTPUT_FORMAT_DEFAULT )
                && properlyPaired )
        {
#ifdef BGS_FWD_FOR_FIRST_REV_FOR_SECOND
            algnResult[i].strand_2 = 2;
#endif

            // for plain and binary, we only output the answer for properly aligned pairs
            // output the alignment of the first read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_1;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_1;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_1;
            occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
            convertToCigarStr2 ( algnResult[i].cigarString_1 );
            occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString_1;
            occ->occPositionCache[occ->occPositionCacheCount].len = readLen1;
            occ->occPositionCacheCount++;

            // output the alignment of the second read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt_2;
            occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand_2;
            occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
            occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score_2;
            occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
            convertToCigarStr2 ( algnResult[i].cigarString_2 );
            occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString_2;
            occ->occPositionCache[occ->occPositionCacheCount].len = readLen2;
            occ->occPositionCacheCount++;
        }
    }

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && preReadID != 0xFFFFFFFF )
    {
#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS

        // ==========================================================
        // for those not properly paired, needs to perform deep-dp for unaligned reads,
        // ==========================================================
        if ( bestAlgn->algnmt_1 == 0xFFFFFFFF && bestAlgn->algnmt_2 == 0xFFFFFFFF && ( ( AllHits * ) hspaux->allHits )->readNum < MAX_UNALIGN_READS_NUM_FOR_DEEP_DP )
        {
            int readID = preReadID;
            int first_read_no_hit  = ( ( hspaux->sa_num[readID] == 0 ) && ( hspaux->occ_num[readID] == 0 ) );
            int second_read_no_hit = ( ( hspaux->sa_num[readID + 1] == 0 ) && ( hspaux->occ_num[readID + 1] == 0 ) );

            if ( first_read_no_hit )
            {
                addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID );
            }
            else
            {
                inputSoap3AnsToArray ( ( AllHits * ) hspaux->allHits, readID, hspaux, bwt, readLen1 );
            }

            if ( second_read_no_hit || ( !first_read_no_hit ) )
            {
                addReadIDToBothUnalignedPairs ( ( ( UnalignedSinglesArrays * ) ( hspaux->readsIDForSingleDP ) )->array[0], readID + 1 );
            }
            else
            {
                inputSoap3AnsToArray ( ( AllHits * ) hspaux->allHits, readID + 1, hspaux, bwt, readLen2 );
            }
        }
        else
        {
            pairDeepDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                                     startIndex, algnNum - startIndex,
                                     upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                     readLen1, readLen2, queryName1, queryName2,
                                     charArray );
        }

        // ==========================================================
#else
        pairDeepDPOutputSAMAPI ( &qInput, algnResult, bestAlgn,
                                 startIndex, algnNum - startIndex,
                                 upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                 readLen1, readLen2, queryName1, queryName2,
                                 charArray );
#endif
    }

    if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API )
    { OCCFlushCacheDP ( &qInput ); }

    OCCFree ( qSetting.occ );

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
    { DynamicUint8ArrayFree ( charArray ); }
}


// output DP single-end alignment results
void outputDPSingleResult ( SingleAlgnmtResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char * upkdQueryNames,
                            unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                            FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index )
{
    // output a set of alignment results
    // the alignments of the same pair of reads have to be together
    // algnNum: the number of results in "algnResult"
    if ( algnNum == 0 )
    { return; }

    HSP * hsp = index->sraIndex->hsp;
    HSPAux * hspaux = index->sraIndex->hspaux;

    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy ( & ( aIndex ), index->sraIndex, sizeof ( SRAIndex ) );

    SRASetting qSetting;
    qSetting.OutFilePtr = outputDPFile;
    qSetting.occ = OCCConstruct ();
    qSetting.OutFileFormat = outputFormat;
    qSetting.SAMOutFilePtr = samOutputDPFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.QueryInfo = &qInfo;
    qInput.AlgnmtIndex = &aIndex;
    OCC * occ = qSetting.occ;
    unsigned int currReadID;
    unsigned int preReadID = 0xFFFFFFFF;
    unsigned int displayID;

    for ( unsigned int i = 0; i < algnNum; i++ )
    {
        // check the readID
        int isUnaligned = ( algnResult[i].algnmt == 0xFFFFFFFF );
        unsigned int readLen = upkdReadLengths[algnResult[i].readID]; // length of the read
        int needOutput = ( ( !isUnaligned ) || ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API ) );

        if ( ( algnResult[i].readID != preReadID ) && needOutput && outputDPFile != NULL )
        {
            // output the information of the read
            currReadID = algnResult[i].readID;
            displayID = currReadID + accumReadNum + 1;
            getReadInfo ( &qInput, currReadID, displayID, upkdQueries, upkdQueryNames, upkdReadLengths, maxReadLength, hspaux->maxLenReadName );
            OCCReportDelimitorDP ( &qInput );
            preReadID = algnResult[i].readID;
        }

        if ( needOutput && outputDPFile != NULL )
        {
            // output the alignment of the read
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            if ( !isUnaligned )
            {
                occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt;
                occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand;
                occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score;
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
                occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString;
                occ->occPositionCache[occ->occPositionCacheCount].len = readLen;
                occ->occPositionCacheCount++;
            }
            else
            {
                OCCReportNoAlignment ( &qInput );
            }

            if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API )
            { OCCFlushCacheDP ( &qInput ); }
        }

        if ( ( !isUnaligned ) && outputDPFile == NULL )
        {
            // output the alignment to the array
            addOCCToArray ( hspaux->algnResultArrays->algnArrays[hspaux->algnResultArrays->numArrays - 1], currReadID, algnResult[i].algnmt, algnResult[i].strand,
                            2, algnResult[i].score );
        }
    }

    if ( outputDPFile != NULL )
    { OCCFlushCacheDP ( &qInput ); }

    OCCFree ( qSetting.occ );
}

void outputDPSingleResult2 ( SingleAlgnmtResult * algnResult, unsigned int algnNum,
                             unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char * upkdQueryNames, char * upkdQualities,
                             unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                             FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index )
{
    if ( algnNum == 0 )
    { return; }

    HSPAux * hspaux = index->sraIndex->hspaux;

#ifdef PERFORM_DEEP_DP_FOR_UNALIGN_READS

    if ( hspaux->readType == PAIR_END_READ )
    {
        // preforming deep dp for the unaligned pair-end reads
        // save all the alignment results to array
        inputAlgnmtsToArray ( ( AllHits * ) hspaux->allHits, algnResult, algnNum );
        return;
    }

#endif
    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy ( & ( aIndex ), index->sraIndex, sizeof ( SRAIndex ) );

    SRASetting qSetting;
    qSetting.OutFilePtr = outputDPFile;
    qSetting.occ = OCCConstruct ();
    qSetting.OutFileFormat = outputFormat;
    qSetting.SAMOutFilePtr = samOutputDPFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.QueryInfo = &qInfo;
    qInput.AlgnmtIndex = &aIndex;
    // check whether the answers are outputted to file
    // or stored into the global array
    int isAnsToFile = ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API ) || ( outputDPFile != NULL );
    OCC * occ = qSetting.occ;
    unsigned int currReadID;
    unsigned int preReadID = 0xFFFFFFFF;
    unsigned int displayID;
    unsigned int startIndex = 0;
    unsigned char upkdPosQuery[MAX_READ_LENGTH];
    unsigned int word_per_query = getWordPerQuery ( maxReadLength );
    DynamicUint8Array * charArray = NULL;

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && isAnsToFile )
    { charArray = DynamicUint8ArrayConstruct (); }

    for ( unsigned int i = 0; i < algnNum; i++ )
    {
        // check the readID
        int isUnaligned = ( algnResult[i].algnmt == 0xFFFFFFFF );
        unsigned int readLen = readLengths[algnResult[i].readID]; // length of the read
        int needOutput = ( ( !isUnaligned ) || ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API ) );

        if ( ( algnResult[i].readID != preReadID ) && needOutput && isAnsToFile )
        {
            // For SAM output format
            if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && preReadID != 0xFFFFFFFF )
            {
                SingleDPOutputSAMAPI ( &qInput, algnResult,
                                       startIndex, i - startIndex, charArray );
            }

            // get the information of the read
            currReadID = algnResult[i].readID;
            displayID = upkdReadIDs[currReadID] + accumReadNum;
            startIndex = i;
            getReadInfoPackedSeq ( &qInput, currReadID, displayID,
                                   queries, readLengths, upkdReadIDs, upkdQueryNames,
                                   word_per_query, upkdPosQuery, upkdQualities, maxReadLength, hspaux->maxLenReadName );

            if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API )
            { OCCReportDelimitorDP ( &qInput ); }

            preReadID = algnResult[i].readID;
        }

        if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API && isAnsToFile )
        {
            // output the alignment of the read for NON-SAM format
            if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( &qInput );}

            if ( !isUnaligned )
            {
                occ->occPositionCache[occ->occPositionCacheCount].tp = algnResult[i].algnmt;
                occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = algnResult[i].strand;
                occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
                occ->occPositionCache[occ->occPositionCacheCount].occMismatch = algnResult[i].score;
                occ->occPositionCache[occ->occPositionCacheCount].resultSource = 1;
                convertToCigarStr2 ( algnResult[i].cigarString );
                occ->occPositionCache[occ->occPositionCacheCount].cigarString = algnResult[i].cigarString;
                occ->occPositionCache[occ->occPositionCacheCount].len = readLen;
                occ->occPositionCacheCount++;
            }
        }

        if ( ( !isUnaligned ) && ( !isAnsToFile ) )
        {
            // output the alignment to the array
            addOCCToArray ( hspaux->algnResultArrays->algnArrays[hspaux->algnResultArrays->numArrays - 1], upkdReadIDs[algnResult[i].readID], algnResult[i].algnmt, algnResult[i].strand,
                            2, algnResult[i].score );
        }
    }

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && preReadID != 0xFFFFFFFF && isAnsToFile )
    {
        SingleDPOutputSAMAPI ( &qInput, algnResult,
                               startIndex, algnNum - startIndex, charArray );
    }

    if ( outputFormat != SRA_OUTPUT_FORMAT_SAM_API && isAnsToFile )
    { OCCFlushCacheDP ( &qInput ); }

    OCCFree ( qSetting.occ );

    if ( outputFormat == SRA_OUTPUT_FORMAT_SAM_API && isAnsToFile )
    { DynamicUint8ArrayFree ( charArray ); }
}


// output the single-end alignment results for the pair-end reads
// only for SAM format
void outputSingleResultForPairEnds ( AllHits * allHits,
                                     unsigned int * queries, unsigned int * readLengths, unsigned int * upkdReadIDs, char * upkdQueryNames, char * upkdQualities,
                                     unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                                     samfile_t * samOutputUnpairFilePtr, Soap3Index * index )
{
    if ( allHits->readNum == 0 )
    {
        return;
    }

    HSPAux * hspaux = index->sraIndex->hspaux;

    SRAIndex aIndex;
    // Copy SRAIndex
    memcpy ( & ( aIndex ), index->sraIndex, sizeof ( SRAIndex ) );

    SRASetting qSetting;
    qSetting.occ = OCCConstruct ();
    qSetting.SAMOutFilePtr = samOutputUnpairFilePtr;
    SRAQueryInput qInput;
    SRAQueryInfo qInfo;
    qInput.QuerySetting = &qSetting;
    qInput.QueryInfo = &qInfo;
    qInput.AlgnmtIndex = &aIndex;
    OCC * occ = qSetting.occ;
    unsigned int firstReadID;
    unsigned int secondReadID;
    unsigned int upkdReadID;
    unsigned int displayID;
    unsigned char upkdPosQuery[MAX_READ_LENGTH];
    unsigned char upkdPosQuery2[MAX_READ_LENGTH];
    char * qualities1, *qualities2;
    unsigned int readLen1, readLen2;
    char * queryName1, *queryName2;
    unsigned int word_per_query = getWordPerQuery ( maxReadLength );
    DynamicUint8Array * charArray = DynamicUint8ArrayConstruct ();
    Algnmt * algn_list1;
    Algnmt * algn_list2;
    unsigned int i, j;

    for ( i = 0; i < allHits->readNum - 1; i += 2 )
    {
        firstReadID = allHits->readPtrArray[i].readID;
        secondReadID = allHits->readPtrArray[i + 1].readID;

        if ( ( firstReadID + 1 ) != secondReadID )
        {
            // Error! It should not happen
            fprintf ( stderr, "[Inside the function: outputSingleResultForPairEnds] Error appears. Reads are not paired-end.\n" );
        }

        algn_list1 = & ( allHits->hitArray[ ( allHits->readPtrArray[i].startIndex )] );
        algn_list2 = & ( allHits->hitArray[ ( allHits->readPtrArray[i + 1].startIndex )] );
        // get the information of both reads
        readLen1 = readLengths[firstReadID]; // length of the first read
        readLen2 = readLengths[secondReadID]; // length of the second read
        // name of the first read
        queryName1 = upkdQueryNames + ( unsigned long long ) ( upkdReadIDs[firstReadID] - 1 ) * hspaux->maxLenReadName;
        // name of the second read
        queryName2 = upkdQueryNames + ( unsigned long long ) ( upkdReadIDs[secondReadID] - 1 ) * hspaux->maxLenReadName;
        // sequence the first read
        retrieve2BWTQueryFormat ( queries, firstReadID, word_per_query, upkdPosQuery );
        // sequence the second read
        retrieve2BWTQueryFormat ( queries, secondReadID, word_per_query, upkdPosQuery2 );
        // base qualities of the first read
        qualities1 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[firstReadID] - 1 ) * maxReadLength;
        // base qualities of the second read
        qualities2 = upkdQualities + ( unsigned long long ) ( upkdReadIDs[secondReadID] - 1 ) * maxReadLength;
        // set qInput for the second read
        qInfo.ReadName = queryName2;
        qInfo.ReadCode = upkdPosQuery2;
        qInfo.ReadLength = readLen2;
        qInfo.ReadQuality = qualities2;
        upkdReadID = upkdReadIDs[secondReadID];
        displayID = upkdReadID + accumReadNum;
        qInfo.ReadId = displayID;
        unproperlypairDPOutputSAMAPI ( &qInput, algn_list1, algn_list2,
                                       allHits->readPtrArray[i].numAlgnmt, allHits->readPtrArray[i + 1].numAlgnmt,
                                       upkdPosQuery, upkdPosQuery2, qualities1, qualities2,
                                       readLen1, readLen2, queryName1, queryName2,
                                       charArray );
    }

    OCCFree ( qSetting.occ );
    DynamicUint8ArrayFree ( charArray );
}
