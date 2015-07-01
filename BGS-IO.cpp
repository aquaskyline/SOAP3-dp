/*
 *
 *    BGS-IO.c
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

#include "BGS-IO.h"


// mapping score
// dimension 1: # of mismatches (i.e. 0, 1, 2, 3, 4) OR
//              difference % from the full DP score (i.e. 0, 1-5%, 6-10%, 11-15%, 15-20%, >20%)
// dimension 2: average mismatch quality value (i.e. 1-10, 11-20, 21-30, 31-40)
//static int mapping_score[6][4] = {{40,40,40,40},{38,37,35,34},{35,33,30,28},{30,28,25,22},{25,22,19,16},{19,16,13,10}};
// static int mapping_score_orig[6][4] = {{MAPQ_MAX,MAPQ_MAX,MAPQ_MAX,MAPQ_MAX},{MAPQ_MAX*0.875,MAPQ_MAX*0.875,MAPQ_MAX*0.85,MAPQ_MAX*0.85},{MAPQ_MAX*0.75,MAPQ_MAX*0.75,MAPQ_MAX*0.7,MAPQ_MAX*0.7},{MAPQ_MAX*0.625,MAPQ_MAX*0.625,MAPQ_MAX*0.55,MAPQ_MAX*0.55},{MAPQ_MAX*0.475,MAPQ_MAX*0.475,MAPQ_MAX*0.4,MAPQ_MAX*0.4},{MAPQ_MAX*0.325,MAPQ_MAX*0.325,MAPQ_MAX*0.25,MAPQ_MAX*0.25}};
static double mapping_score[6][2] = {{1.0, 1.0}, {0.875, 0.85}, {0.75, 0.7}, {0.625, 0.55}, {0.475, 0.4}, {0.325, 0.25}};


//===========================================================//
// The following is for the updated MAPQ scoring function    //
// for single-end reads for DP module (Date: Oct 19, 2012)   //
//===========================================================//

// The penality score when consideration of the average quality value in mismatch positions
static float penalty_score_avg_mis_qual[41] = {3, 2.85, 2.71, 2.57, 2.43, 2.3, 2.17, 2.04, 1.92, 1.8, 1.69, 1.58, 1.47, 1.37, 1.27, 1.17, 1.08, 0.99, 0.91, 0.83, 0.75, 0.68, 0.61, 0.54, 0.48, 0.42, 0.37, 0.32, 0.27, 0.23, 0.19, 0.15, 0.12, 0.09, 0.07, 0.05, 0.03, 0.02, 0.01, 0, 0};

// The penality ratio when consideration of X1
static float penalty_ratio_x1[101] = {1, 0.5, 0.33, 0.25, 0.2, 0.17, 0.14, 0.13, 0.11, 0.1, 0.09, 0.08, 0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

//


//Define the below parameter to output the alignment result(text position)
// on screen instead of writing into output file
//#define DEBUG_2BWT_OUTPUT_TO_SCREEN

//Define the below parameter to output the alignment result(text position)
// with the old 2BWT-Aligner format
//#define DEBUG_2BWT_OUTPUT_32_BIT

//Define the below parameter to stop cache reported SA range and text position.
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_OUTPUT

//Define the below parameter to skip writing the alignment result(text position)
// into disk. This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_WRITE_FILE

//Define the below parameter to skip translating the alignment result(text position).
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_TRANSLATION

OCC * OCCConstruct ()
{
    OCC * occ = ( OCC * ) malloc ( sizeof ( OCC ) );
    SAMOccurrenceConstruct ( occ );
    occ->occPositionCacheCount = 0;
    return occ;
}


void OCCFree ( OCC * occ )
{
    SAMOccurrenceDestruct ( occ );
    free ( occ );
}

unsigned int OCCWriteOutputHeader ( HSP * hsp, FILE * outFilePtr,
                                    unsigned int maxReadLength,
                                    unsigned int numOfReads,
                                    int outputFormat )
{
    unsigned int * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    unsigned int tp, approxIndex, approxValue;
    unsigned int outputVersion = OCC_OUTPUT_FORMAT;
    unsigned int writtenByteCount = 0;
    int i;

    switch ( outputFormat )
    {
        case SRA_OUTPUT_FORMAT_SAM_API:
            //Ad-hoc writing header is not supported by API
            break;

        case SRA_OUTPUT_FORMAT_SAM:
            fprintf ( outFilePtr, "@HD\tVN:1.4\tSO:unsorted\n" );

            for ( i = 0; i < hsp->numOfSeq; i++ )
            {
                int j;

                for ( j = 0; j < 255; j++ )
                {
                    if ( hsp->annotation[i].text[j] == '\0' ||
                            hsp->annotation[i].text[j] == ' ' ||
                            hsp->annotation[i].text[j] == '\t' ||
                            hsp->annotation[i].text[j] == '\r' ||
                            hsp->annotation[i].text[j] == '\n' )
                    {
                        break;
                    }
                }

                hsp->annotation[i].text[j] = '\0';
                tp = hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos + maxReadLength;
                fprintf ( outFilePtr, "@SQ\tSN:%s\tLN:%u\n", hsp->annotation[i].text, tp );
            }

            fprintf ( outFilePtr, "@PG\tID:%s\tVN:v%d.%d.%d (%s)\n", PROJECT_NAME, PROJECT_MAJOR, PROJECT_MINOR, PROJECT_REV, PROJECT_SPECIAL );
            break;

        case SRA_OUTPUT_FORMAT_PLAIN:
#ifndef SKIP_PLAIN_HEADER
            fprintf ( outFilePtr, "SOAP3 Plain Output File\n" );
            fprintf ( outFilePtr, "========================\n" );
            fprintf ( outFilePtr, "Read#  ChromId#  Offset  Strand  #Mismatch\n" );
#endif
            break;

            //Implicit case SRA_OUTPUT_FORMAT_DEFAULT:
        default:
            fwrite ( & ( outputVersion ), sizeof ( int ), 1, outFilePtr );
            writtenByteCount += sizeof ( int );
            fwrite ( & ( maxReadLength ), sizeof ( int ), 1, outFilePtr );
            writtenByteCount += sizeof ( int );
            fwrite ( & ( numOfReads ), sizeof ( int ), 1, outFilePtr );
            writtenByteCount += sizeof ( int );
            fwrite ( & ( hsp->numOfSeq ), sizeof ( unsigned int ), 1, outFilePtr );
            writtenByteCount += sizeof ( unsigned int );

            for ( i = 0; i < hsp->numOfSeq; i++ )
            {
                fwrite ( &hsp->annotation[i].gi, sizeof ( int ), 1, outFilePtr );
                writtenByteCount += sizeof ( int );
                fwrite ( &hsp->annotation[i].text, sizeof ( char ) * ( MAX_SEQ_NAME_LENGTH + 1 ), 1, outFilePtr );
                writtenByteCount += sizeof ( char ) * ( MAX_SEQ_NAME_LENGTH + 1 );
                tp = hsp->seqOffset[i].endPos;
                approxIndex = tp >> GRID_SAMPLING_FACTOR_2_POWER;
                approxValue = ambiguityMap[approxIndex];

                while ( translate[approxValue].startPos > tp )
                {
                    approxValue--;
                }

                tp -= translate[approxValue].correction;
                fwrite ( &tp, sizeof ( unsigned int ), 1, outFilePtr );
                writtenByteCount += sizeof ( unsigned int );
            }

            // cx Header? I don't know if we need to change this function
            // have HSP

    }

    return writtenByteCount;
}

void AssignCigarStrToSAMIU ( bam1_t * samAlgnmt, int * curSize,
                             char * cigar_string )
{
    char operCode[256];
    operCode['M'] = 0; // match or mismatch
    operCode['I'] = 1; // insertion to the reference
    operCode['D'] = 2; // deletion from the reference
    operCode['S'] = 4; // soft clipping
    samAlgnmt->core.n_cigar = 0;
    int num = 0;

    for ( int i = 0; i < strlen ( cigar_string ); i++ )
    {
        if ( cigar_string[i] >= '0' && cigar_string[i] <= '9' )
        {
            num = num * 10 + cigar_string[i] - '0';
        }
        else if ( num > 0 )
        {
            SAMIUint8ConcatUint32 ( samAlgnmt->data, curSize, num << BAM_CIGAR_SHIFT | operCode[cigar_string[i]] ); //CIGAR
            samAlgnmt->core.n_cigar++;                                //SAM: number of CIGAR operations
            num = 0;
        }
    }
}

void OCCFlushCachePlain ( OCC * occ, HSP * hsp, FILE * outFilePtr, SRAQueryInput * qInput )
{
    char strandStr[4] = "?+-";
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    OCCPositionCache lastDelimEntry;

    if ( occ->occPositionCacheCount > 0 )
    {
        // ALERT: THIS PART IS NOT 64-BIT COMPATIBLE.
        unsigned long long j;
        unsigned int tp;
        unsigned int approxIndex, approxValue;
        unsigned long long k = 0;
        j = 0;

        for ( k = 0; k < occ->occPositionCacheCount; k++ )
        {
            if ( occ->occPositionCache[j].ChromId == OCC_CONST_ALIGNMENT_HEADER )
            {
                lastDelimEntry.ChromId = occ->occPositionCache[j].ChromId;
                lastDelimEntry.occMismatch = occ->occPositionCache[j].occMismatch;
                lastDelimEntry.ReadStrand = occ->occPositionCache[j].ReadStrand;
                lastDelimEntry.tp = occ->occPositionCache[j].tp;
            }
            else if ( occ->occPositionCache[j].ChromId == OCC_CONST_NO_ALIGNMENT )
            {
            }
            else if ( occ->occPositionCache[j].ReadStrand != 0 )
            {
                tp = ( unsigned int ) occ->occPositionCache[j].tp;
#ifndef DEBUG_2BWT_NO_TRANSLATION
                approxIndex = tp >> GRID_SAMPLING_FACTOR_2_POWER;
                approxValue = occAmbiguityMap[approxIndex];

                while ( occTranslate[approxValue].startPos > tp )
                {
                    approxValue--;
                }

                tp -= occTranslate[approxValue].correction;
                // cx non-DP translation check tp
                // have HSP

                // ADD_CX
                unsigned int pacPos = ( unsigned int ) occ->occPositionCache[j].tp;
                unsigned int chrEndPos = hsp->seqOffset[occTranslate[approxValue].chrID - 1].endPos;
                int readLength = occ->occPositionCache[j].len;
                char * buffer;
                int correctedChrID;
                unsigned int correctedChrPos;

                if ( 0 /*BoundaryCheck(pacPos, occTranslate[approxValue].chrID, chrEndPos, readLength, hsp, correctedChrID, buffer)*/ )
                {
                    fprintf ( outFilePtr, "AAA+ %llu %u ???? %c ? %s\n", lastDelimEntry.tp, correctedChrID, /*( unsigned long long ) tp,*/ strandStr[occ->occPositionCache[j].ReadStrand], buffer );
                    free ( buffer );
                }
                else
                {

#endif
#ifdef DEBUG_2BWT_NO_TRANSLATION
                    occ->occPositionCacheToDisk[j].cell[1] = 0;
#endif
#ifdef DEBUG_2BWT_OUTPUT_TO_SCREEN
                    printf ( "[OCCFlushCache] %d %d \n", occTranslate[approxValue].chrID, tp );
#endif
                    fprintf ( outFilePtr, "AAA %llu %u %llu %c %d\n", lastDelimEntry.tp, occTranslate[approxValue].chrID, ( unsigned long long ) tp, strandStr[occ->occPositionCache[j].ReadStrand], occ->occPositionCache[j].occMismatch );
                }
            }

            j++;
        }
    }

    occ->occPositionCache[0].ChromId =       lastDelimEntry.ChromId;
    occ->occPositionCache[0].occMismatch =   lastDelimEntry.occMismatch;
    occ->occPositionCache[0].ReadStrand =    lastDelimEntry.ReadStrand;
    occ->occPositionCache[0].tp =            lastDelimEntry.tp;
    occ->occPositionCache[0].len = lastDelimEntry.len;
    occ->occPositionCacheCount = 1;
}

void OCCFlushCachePlainDP ( OCC * occ, HSP * hsp, FILE * outFilePtr )
{
    char strandStr[4] = "?+-";
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    OCCPositionCache lastDelimEntry;

    if ( occ->occPositionCacheCount > 0 )
    {
        // ALERT: THIS PART IS NOT 64-BIT COMPATIBLE.
        unsigned long long j;
        unsigned int tp;
        unsigned int approxIndex, approxValue;
        unsigned long long k = 0;
        j = 0;

        for ( k = 0; k < occ->occPositionCacheCount; k++ )
        {
            if ( occ->occPositionCache[j].ChromId == OCC_CONST_ALIGNMENT_HEADER )
            {
                lastDelimEntry.ChromId = occ->occPositionCache[j].ChromId;
                lastDelimEntry.occMismatch = occ->occPositionCache[j].occMismatch;
                lastDelimEntry.ReadStrand = occ->occPositionCache[j].ReadStrand;
                lastDelimEntry.tp = occ->occPositionCache[j].tp;
                lastDelimEntry.resultSource = occ->occPositionCache[j].resultSource;
                lastDelimEntry.cigarString = occ->occPositionCache[j].cigarString;
                lastDelimEntry.len = occ->occPositionCache[j].len;
            }
            else if ( occ->occPositionCache[j].ChromId == OCC_CONST_NO_ALIGNMENT )
            {
            }
            else if ( occ->occPositionCache[j].ReadStrand != 0 )
            {
                tp = ( unsigned int ) occ->occPositionCache[j].tp;
#ifndef DEBUG_2BWT_NO_TRANSLATION
                approxIndex = tp >> GRID_SAMPLING_FACTOR_2_POWER;
                approxValue = occAmbiguityMap[approxIndex];

                while ( occTranslate[approxValue].startPos > tp )
                {
                    approxValue--;
                }

                tp -= occTranslate[approxValue].correction;
                // cx both non-DP and DP results (occ->occPositionCache[j].resultSource = 1). check tp
                // have HSP
#endif

                // ADD_CX  .... handle non-DP alignment only
                if ( occ->occPositionCache[j].resultSource == 0 )
                {



                }
                else     // handle DP alignment here
                {


                }

#ifdef DEBUG_2BWT_NO_TRANSLATION
                occ->occPositionCacheToDisk[j].cell[1] = 0;
#endif
#ifdef DEBUG_2BWT_OUTPUT_TO_SCREEN
                printf ( "[OCCFlushCache] %d %d \n", occTranslate[approxValue].chrID, tp );
#endif

                if ( occ->occPositionCache[j].resultSource == 0 )
                {
                    // the result is from SOAP3
                    unsigned int pacPos = ( unsigned int ) occ->occPositionCache[j].tp;
                    unsigned int chrEndPos = hsp->seqOffset[occTranslate[approxValue].chrID - 1].endPos;
                    int readLength = occ->occPositionCache[j].len;
                    char * buffer;
                    int correctedChrID;
                    unsigned int correctedChrPos;

                    if ( 0 /*BoundaryCheck(pacPos, occTranslate[approxValue].chrID, chrEndPos, readLength, hsp, correctedChrID, buffer) */ )
                    {
                        fprintf ( outFilePtr, "XXX+ %llu %u ???? %c ? %d %s\n", lastDelimEntry.tp, correctedChrID, /*( unsigned long long ) tp,*/ strandStr[occ->occPositionCache[j].ReadStrand], occ->occPositionCache[j].resultSource, buffer );
                        free ( buffer );
                    }
                    else
                    {
                        fprintf ( outFilePtr, "XXX %llu %u %llu %c %d %d %dM\n", lastDelimEntry.tp, occTranslate[approxValue].chrID, ( unsigned long long ) tp, strandStr[occ->occPositionCache[j].ReadStrand], occ->occPositionCache[j].occMismatch, occ->occPositionCache[j].resultSource, occ->occPositionCache[j].len );
                    }
                }
                else
                {
                    // the result is from DP
                    unsigned int pacPos = ( unsigned int ) occ->occPositionCache[j].tp;
                    unsigned int chrEndPos = hsp->seqOffset[occTranslate[approxValue].chrID - 1].endPos;
                    int readLength = occ->occPositionCache[j].len;
                    int correctedChrID;
                    unsigned int correctedChrPos;
                    char * buffer;

                    if ( 0 /*BoundaryCheckDP(pacPos, occTranslate[approxValue].chrID, chrEndPos, readLength, occ->occPositionCache[j].cigarString,
                        hsp, correctedChrID, buffer) */ )
                    {
                        fprintf ( outFilePtr, "YYY+ %llu %u ???? %c ? %d %s\n", lastDelimEntry.tp, correctedChrID, strandStr[occ->occPositionCache[j].ReadStrand], occ->occPositionCache[j].resultSource, buffer );
                        free ( buffer );
                    }
                    else
                    {
                        fprintf ( outFilePtr, "YYY %llu %u %llu %c %d %d %s\n", lastDelimEntry.tp, occTranslate[approxValue].chrID, ( unsigned long long ) tp, strandStr[occ->occPositionCache[j].ReadStrand], occ->occPositionCache[j].occMismatch, occ->occPositionCache[j].resultSource, occ->occPositionCache[j].cigarString );
                    }
                }
            }

            j++;
        }
    }

    occ->occPositionCache[0].ChromId =       lastDelimEntry.ChromId;
    occ->occPositionCache[0].occMismatch =   lastDelimEntry.occMismatch;
    occ->occPositionCache[0].ReadStrand =    lastDelimEntry.ReadStrand;
    occ->occPositionCache[0].tp =            lastDelimEntry.tp;
    occ->occPositionCache[0].resultSource =  lastDelimEntry.resultSource;
    occ->occPositionCache[0].cigarString =   lastDelimEntry.cigarString;
    occ->occPositionCache[0].len =           lastDelimEntry.len;
    occ->occPositionCacheCount = 1;
}

// binary output. ignore.
void OCCFlushCacheDefault ( OCC * occ, HSP * hsp, FILE * outFilePtr )
{
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;

    if ( occ->occPositionCacheCount > 0 )
    {
        // ALERT: THIS PART IS NOT 64-BIT COMPATIBLE.
        unsigned long long j;
        unsigned int tp;
        unsigned int approxIndex, approxValue;
        unsigned long long k = 0;
        j = 0;

        for ( k = 0; k < occ->occPositionCacheCount; k++ )
        {
            if ( occ->occPositionCache[j].ChromId == OCC_CONST_ALIGNMENT_HEADER )
            {
#ifdef DEBUG_2BWT_OUTPUT_32_BIT
                ( * ( unsigned int * ) & ( occ->occPositionCacheToDisk[j].cell[1] ) ) = ( unsigned int ) occ->occPositionCache[j].tp;
                occ->occPositionCacheToDisk[j].cell[0] = 0;
#endif
#ifndef DEBUG_2BWT_OUTPUT_32_BIT
                ( * ( unsigned long long * ) & ( occ->occPositionCacheToDisk[j].cell[3] ) ) = ( unsigned long long ) occ->occPositionCache[j].tp;
                occ->occPositionCacheToDisk[j].cell[1] = 0;
                occ->occPositionCacheToDisk[j].cell[0] = 0;
#endif
            }
            else if ( occ->occPositionCache[j].ChromId == OCC_CONST_NO_ALIGNMENT )
            {
            }
            else
            {
                tp = ( unsigned int ) occ->occPositionCache[j].tp;
#ifndef DEBUG_2BWT_NO_TRANSLATION
                approxIndex = tp >> GRID_SAMPLING_FACTOR_2_POWER;
                approxValue = occAmbiguityMap[approxIndex];

                while ( occTranslate[approxValue].startPos > tp )
                {
                    approxValue--;
                }

                tp -= occTranslate[approxValue].correction;
                // cx non-DP??? check tp
                // have HSP
#endif
#ifdef DEBUG_2BWT_NO_TRANSLATION
#ifdef DEBUG_2BWT_OUTPUT_32_BIT
                occ->occPositionCacheToDisk[j].cell[0] = 0;
#endif
#ifndef DEBUG_2BWT_OUTPUT_32_BIT
                occ->occPositionCacheToDisk[j].cell[1] = 0;
#endif
#endif
#ifdef DEBUG_2BWT_OUTPUT_TO_SCREEN
                printf ( "[OCCFlushCache] %d %d \n", occTranslate[approxValue].chrID, tp );
#endif
#ifdef DEBUG_2BWT_OUTPUT_32_BIT
                ( * ( unsigned int * ) & ( occ->occPositionCacheToDisk[j].cell[1] ) ) = ( unsigned int ) tp;
                occ->occPositionCacheToDisk[j].cell[0] = occTranslate[approxValue].chrID;
#endif
#ifndef DEBUG_2BWT_OUTPUT_32_BIT
                ( * ( unsigned long long * ) & ( occ->occPositionCacheToDisk[j].cell[3] ) ) = ( unsigned long long ) tp;
                occ->occPositionCacheToDisk[j].cell[0] = occ->occPositionCache[j].ReadStrand;
                occ->occPositionCacheToDisk[j].cell[1] = occTranslate[approxValue].chrID;
#endif
            }

            j++;
        }

#ifndef DEBUG_2BWT_NO_WRITE_FILE
        fwrite ( occ->occPositionCacheToDisk, sizeof ( OCCPositionCacheToDisk ) *occ->occPositionCacheCount, 1, outFilePtr );
#endif
    }

    occ->occPositionCacheCount = 0;
}

/*
void OCCDirectWritePairUnmapSAM ( SRAQueryInput * qInput, SRAOccurrence * sraOcc, int isFirst )
{
   SRASetting * qSetting = qInput->QuerySetting;
   // int OutFileFormat = qSetting->OutFileFormat;
   FILE * outFilePtr = qSetting->OutFilePtr;
   SRAQueryInfo * qInfo = qInput->QueryInfo;
   SRAIndex * aIndex = qInput->AlgnmtIndex;
   HSP * hsp = aIndex->hsp;
   Translate * occTranslate = hsp->translate;
   unsigned int * occAmbiguityMap = hsp->ambiguityMap;
   // unsigned long long j=0, k=0;
   unsigned long long ambPosition;
   unsigned long long tp;
   unsigned int correctPosition;
   unsigned int approxIndex, approxValue;
   unsigned long long occCount = 0;
   // ATTENTION: THIS PART IS NOT 64-BIT COMPATIBLE.
   {
       ambPosition = sraOcc->ambPosition;
       tp = ambPosition;
       correctPosition = ambPosition;
       approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
       approxValue = occAmbiguityMap[approxIndex];

       while ( occTranslate[approxValue].startPos > ambPosition )
       {
           approxValue--;
       }

       correctPosition -= occTranslate[approxValue].correction;
       tp = correctPosition;
       //QNAME FLAG RNAME POS
       //MAPQ CIGAR RNEXT PNEXT TLEN
       //SEQ QUAL
       unsigned int samFlag = 0;
       samFlag |= SAM_FLAG_MATE_UNMAPPED;
       samFlag |= SAM_FLAG_FIRST_IN_PAIR * isFirst;
       samFlag |= SAM_FLAG_SECOND_IN_PAIR * ( 1 - isFirst );

       if ( sraOcc->strand == QUERY_NEG_STRAND )
       {
           samFlag |= SAM_FLAG_READ_ALGNMT_STRAND * isFirst;
       }

       samFlag |= isFirst *
                  fprintf ( outFilePtr, "%s\t%u\t%s\t%llu\t%d\t%lluM\t%s\t%llu\t%d\t",
                            qInfo->ReadName,         //QNAME
                            samFlag,                  //FLAG
                            hsp->annotation[occTranslate[approxValue].chrID - 1].text,  //RNAME
                            tp,                      //POS
                            SAM_MAPQ_UNAVAILABLE,    //MAPQ
                            qInfo->ReadLength,          //CIGAR
                            "*",                      //MRNM
                            0LL,                      //MPOS
                            0                        //ISIZE
                          );
       int i;

       for ( i = 0; i < qInfo->ReadLength; i++ )
       {
           fprintf ( outFilePtr, "%c", dnaChar[qInfo->ReadCode[i]] ); //SEQ
       }

       fprintf ( outFilePtr, "\t" );

       for ( i = 0; i < qInfo->ReadLength; i++ )
       {
           fprintf ( outFilePtr, "%c", 33 );                           //QUAL
       }

       fprintf ( outFilePtr, "\n" );
       occCount++;
   }
}

void OCCDirectWritePairOccSAM ( SRAQueryInput * qInput, PEPairs * pePair )
{
   SRASetting * qSetting = qInput->QuerySetting;
   SRAIndex * aIndex = qInput->AlgnmtIndex;
   // int OutFileFormat = qSetting->OutFileFormat;
   SRAQueryInfo * qInfo = qInput->QueryInfo;
   FILE * outFilePtr = qSetting->OutFilePtr;
   HSP * hsp = aIndex->hsp;
   Translate * occTranslate = hsp->translate;
   unsigned int * occAmbiguityMap = hsp->ambiguityMap;
   // unsigned long long j=0, k=0;
   unsigned long long ambPosition;
   unsigned long long tp_1, tp_2;
   char chr_1, chr_2;
   unsigned int correctPosition;
   unsigned int approxIndex, approxValue;
   // unsigned long long occCount = 0;
   char * mateChar;
   char equalChr[3] = "=";
   // ATTENTION: THIS PART IS NOT 64-BIT COMPATIBLE.
   {
       ambPosition = pePair->algnmt_1;
       tp_1 = ambPosition;
       correctPosition = ambPosition;
       approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
       approxValue = occAmbiguityMap[approxIndex];

       while ( occTranslate[approxValue].startPos > ambPosition )
       {
           approxValue--;
       }

       correctPosition -= occTranslate[approxValue].correction;
       tp_1 = correctPosition;
       chr_1 = occTranslate[approxValue].chrID - 1;
   }
   {
       ambPosition = pePair->algnmt_2;
       tp_2 = ambPosition;
       correctPosition = ambPosition;
       approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
       approxValue = occAmbiguityMap[approxIndex];

       while ( occTranslate[approxValue].startPos > ambPosition )
       {
           approxValue--;
       }

       correctPosition -= occTranslate[approxValue].correction;
       tp_2 = correctPosition;
       chr_2 = occTranslate[approxValue].chrID - 1;
   }

   if ( chr_1 == chr_2 ) { mateChar = equalChr; }
   else { mateChar = hsp->annotation[chr_2].text; }

   unsigned int samFlag = 0;
   samFlag |= SAM_FLAG_PROPER_PAIR;
   samFlag |= SAM_FLAG_FIRST_IN_PAIR;

   if ( pePair->strand_1 == QUERY_NEG_STRAND )
   {
       samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
   }

   if ( pePair->strand_2 == QUERY_NEG_STRAND )
   {
       samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
   }

   //QNAME FLAG RNAME POS
   //MAPQ CIGAR RNEXT PNEXT TLEN
   //SEQ QUAL
   fprintf ( outFilePtr, "%s\t%u\t%s\t%llu\t%d\t%lluM\t%s\t%llu\t%d\t",
             qInfo->ReadName,         //QNAME
             samFlag,                  //FLAG
             hsp->annotation[chr_1].text,    //RNAME
             tp_1,                       //POS
             SAM_MAPQ_UNAVAILABLE,     //MAPQ  (unavailable)
             qInfo->ReadLength,       //CIGAR
             mateChar,                 //MRNM
             tp_2,                     //MPOS
             pePair->insertion         //ISIZE (unavailable in single-end)
           );
   int i;

   for ( i = 0; i < qInfo->ReadLength; i++ )
   {
       fprintf ( outFilePtr, "%c", dnaChar[qInfo->ReadCode[i]] ); //SEQ
   }

   fprintf ( outFilePtr, "\t" );

   for ( i = 0; i < qInfo->ReadLength; i++ )
   {
       fprintf ( outFilePtr, "%c", 33 );                           //QUAL
   }

   fprintf ( outFilePtr, "\n" );
}
*/

void OCCFlushCacheSAM ( SRAQueryInput * qInput )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    // int OutFileFormat = qSetting->OutFileFormat;
    FILE * outFilePtr = qSetting->OutFilePtr;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    unsigned long long j = 0, k = 0;
    unsigned long long ambPosition;
    unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;
    unsigned long long occCount = 0;

    for ( k = 0; k < occ->occPositionCacheCount; k++ )
    {
        if ( occ->occPositionCache[j].ChromId == OCC_CONST_ALIGNMENT_HEADER )
        {
        }
        else if ( occ->occPositionCache[j].ChromId == OCC_CONST_NO_ALIGNMENT )
        {
        }
        else
        {
            // ATTENTION: THIS PART IS NOT 64-BIT COMPATIBLE.
            ambPosition = occ->occPositionCache[j].tp;
            tp = ambPosition;
            correctPosition = ambPosition;
            approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
            approxValue = occAmbiguityMap[approxIndex];

            while ( occTranslate[approxValue].startPos > ambPosition )
            {
                approxValue--;
            }

            correctPosition -= occTranslate[approxValue].correction;
            tp = correctPosition;
            // cx non-DP, check tp
            // have HSP as local variable
            unsigned int samFlag = 0;

            if ( occ->occPositionCache[j].ReadStrand == QUERY_NEG_STRAND )
            {
                samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
            }

            //QNAME FLAG RNAME POS
            //MAPQ CIGAR RNEXT PNEXT TLEN
            //SEQ QUAL
            fprintf ( outFilePtr, "%s\t%u\t%s\t%llu\t%d\t%lluM\t%s\t%llu\t%d\t",
                      qInfo->ReadName,         //QNAME
                      samFlag,                  //FLAG
                      hsp->annotation[occTranslate[approxValue].chrID - 1].text,  //RNAME
                      tp,                       //POS
                      SAM_MAPQ_UNAVAILABLE,     //MAPQ  (unavailable)
                      qInfo->ReadLength,       //CIGAR
                      "*",                      //MRNM  (unavailable in single-end)
                      0LL,                      //MPOS  (unavailable in single-end)
                      0                         //ISIZE (unavailable in single-end)
                    );
            int i;

            for ( i = 0; i < qInfo->ReadLength; i++ )
            {
                fprintf ( outFilePtr, "%c", dnaChar[qInfo->ReadCode[i]] ); //SEQ
            }

            fprintf ( outFilePtr, "\t" );

            for ( i = 0; i < qInfo->ReadLength; i++ )
            {
                fprintf ( outFilePtr, "%c", 33 );                           //QUAL
            }

            fprintf ( outFilePtr, "\n" );
            occCount++;
        }

        j++;
    }

    if ( occCount == 0 )
    {
        unsigned int samFlag = 0;
        samFlag |= SAM_FLAG_READ_UNMAPPED;
        fprintf ( outFilePtr, "%s\t%u\t%s\t%llu\t%d\t%s\t%s\t%llu\t%d\t%s\t%s\n",
                  qInfo->ReadName,         //QNAME
                  samFlag,                  //FLAG
                  "*",                      //RNAME
                  0LL,                      //POS
                  0,                        //MAPQ
                  "*",                      //CIGAR
                  "*",                      //MRNM
                  0LL,                      //MPOS
                  0,                        //ISIZE
                  "*",                      //SEQ
                  "*"                       //QUAL
                );
    }

    qSetting->writtenCacheCount++;
    occ->occPositionCacheCount = 0;
}

void OCCFlushCacheSAMAPI ( SRAQueryInput * qInput )
{
    //  fprintf(stderr,"OCCFlushCacheSAMAPI\n");
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    // int OutFileFormat = qSetting->OutFileFormat;
    // FILE * outFilePtr = qSetting->OutFilePtr;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long i = 0, j = 0, k = 0;
    unsigned long long ambPosition;
    unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;
    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.n_cigar = 1;                                //SAM: number of CIGAR operations
    samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
    samAlgnmt->core.mtid = -1;
    samAlgnmt->core.mpos = -1;
    samAlgnmt->core.isize = 0;
    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
    SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadLength << BAM_CIGAR_SHIFT ); //CIGAR

    for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
    {
        unsigned char biChar = 0;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    if ( qInfo->ReadLength % 2 == 1 )
    {
        uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
    {
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), 0xff );
    }

    unsigned int auxStart = samAlgnmt->data_len;
    //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
    //----------------------------------------------------------------------------------------------------------------------------
    unsigned long long occCount = 0;

    for ( k = 0; k < occ->occPositionCacheCount; k++ )
    {
        if ( occ->occPositionCache[j].ChromId == OCC_CONST_ALIGNMENT_HEADER )
        {
        }
        else if ( occ->occPositionCache[j].ChromId == OCC_CONST_NO_ALIGNMENT )
        {
            unsigned int samFlag = 0;
            samFlag |= SAM_FLAG_READ_UNMAPPED;
            // ---------------------->  |
            samAlgnmt->core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
            samAlgnmt->core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
            samAlgnmt->core.qual = 0;                                   //SAM: mapping quality
            samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
            samwrite ( samFilePtr, samAlgnmt );
            // ---------------------->  |
        }
        else
        {
            // ATTENTION: THIS PART IS NOT 64-BIT COMPATIBLE.
            ambPosition = occ->occPositionCache[j].tp;
            tp = ambPosition;
            correctPosition = ambPosition;
            approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
            approxValue = occAmbiguityMap[approxIndex];

            while ( occTranslate[approxValue].startPos > ambPosition )
            {
                approxValue--;
            }

            correctPosition -= occTranslate[approxValue].correction;
            tp = correctPosition;
            // cx non-DP, check tp
            // have HSP
            // modify CIGAR string seems hard, probably need call SAMIUintxxxxx functions, or AssignCIGARStrToSAMIU
            unsigned int samFlag = 0;

            if ( occ->occPositionCache[j].ReadStrand == QUERY_NEG_STRAND )
            {
                samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
            }

            // ---------------------->  |
            samAlgnmt->core.tid = occTranslate[approxValue].chrID - 1;  //SAM: chromosome ID, defined by bam_header_t
            samAlgnmt->core.pos = tp - 1;                               //SAM: 0-based leftmost coordinate
            samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
            samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
            samwrite ( samFilePtr, samAlgnmt );
            // ---------------------->  |
            occCount++;
        }

        j++;
    }

    qSetting->writtenCacheCount++;
    occ->occPositionCacheCount = 0;
}

void OCCFlushCacheSAMAPIDP ( SRAQueryInput * qInput )
{
    //  fprintf(stderr,"OCCFlushCacheSAMAPIDP\n");
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    // int OutFileFormat = qSetting->OutFileFormat;
    // FILE * outFilePtr = qSetting->OutFilePtr;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long i = 0, j = 0, k = 0;
    unsigned long long ambPosition;
    unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;
    unsigned long long occCount = 0;

    for ( k = 0; k < occ->occPositionCacheCount; k++ )
    {
        if ( occ->occPositionCache[j].ChromId == OCC_CONST_ALIGNMENT_HEADER )
        {
        }
        else if ( occ->occPositionCache[j].ChromId == OCC_CONST_NO_ALIGNMENT )
        {
            //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
            // Each occurrences reported by SAM_API will share the following aux data.
            //----------------------------------------------------------------------------------------------------------------------------
            samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
            samAlgnmt->core.n_cigar = 1;                                //SAM: number of CIGAR operations
            samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
            samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
            samAlgnmt->core.mtid = -1;
            samAlgnmt->core.mpos = -1;
            samAlgnmt->core.isize = 0;
            samAlgnmt->l_aux = 0;
            samAlgnmt->data_len = 0;
            samAlgnmt->m_data = SAM_MDATA_SIZE;
            SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
            SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadLength << BAM_CIGAR_SHIFT ); //CIGAR

            for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
                biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            if ( qInfo->ReadLength % 2 == 1 )
            {
                uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
            {
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), 0xff );
            }

            unsigned int auxStart = samAlgnmt->data_len;
            //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
            samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
            //----------------------------------------------------------------------------------------------------------------------------
            unsigned int samFlag = 0;
            samFlag |= SAM_FLAG_READ_UNMAPPED;
            // ---------------------->  |
            samAlgnmt->core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
            samAlgnmt->core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
            samAlgnmt->core.qual = 0;                                   //SAM: mapping quality
            samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
            samwrite ( samFilePtr, samAlgnmt );
            // ---------------------->  |
        }
        else
        {
            //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
            // Each occurrences reported by SAM_API will share the following aux data.
            //----------------------------------------------------------------------------------------------------------------------------
            samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
            samAlgnmt->core.n_cigar = 1;                                //SAM: number of CIGAR operations
            samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
            samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
            samAlgnmt->core.mtid = -1;
            samAlgnmt->core.mpos = -1;
            samAlgnmt->core.isize = 0;
            samAlgnmt->l_aux = 0;
            samAlgnmt->data_len = 0;
            samAlgnmt->m_data = SAM_MDATA_SIZE;
            SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
            AssignCigarStrToSAMIU ( samAlgnmt, & ( samAlgnmt->data_len ), occ->occPositionCache[j].cigarString );

            for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
                biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            if ( qInfo->ReadLength % 2 == 1 )
            {
                uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
            {
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), 0xff );
            }

            unsigned int auxStart = samAlgnmt->data_len;
            //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
            samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
            //----------------------------------------------------------------------------------------------------------------------------
            // ATTENTION: THIS PART IS NOT 64-BIT COMPATIBLE.
            ambPosition = occ->occPositionCache[j].tp;
            tp = ambPosition;
            correctPosition = ambPosition;
            approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
            approxValue = occAmbiguityMap[approxIndex];

            while ( occTranslate[approxValue].startPos > ambPosition )
            {
                approxValue--;
            }

            correctPosition -= occTranslate[approxValue].correction;
            tp = correctPosition;
            // cx DP. have HSP
            // need use SAMAPI to modify CIGAR string
            unsigned int samFlag = 0;

            if ( occ->occPositionCache[j].ReadStrand == QUERY_NEG_STRAND )
            {
                samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
            }

            // ---------------------->  |
            samAlgnmt->core.tid = occTranslate[approxValue].chrID - 1;  //SAM: chromosome ID, defined by bam_header_t
            samAlgnmt->core.pos = tp - 1;                               //SAM: 0-based leftmost coordinate
            samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
            samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
            samwrite ( samFilePtr, samAlgnmt );
            // ---------------------->  |
            occCount++;
        }

        j++;
    }

    qSetting->writtenCacheCount++;
    occ->occPositionCacheCount = 0;
}


void OCCDirectWritePairOccSAMAPI ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, int properlyPaired )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // output the pair based on the first read
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    // int OutFileFormat = qSetting->OutFileFormat;
    // FILE * outFilePtr = qSetting->OutFilePtr;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long i = 0; //, j=0, k=0;
    unsigned long long ambPosition;
    unsigned long long tp_1, tp_2;
    char chr_1, chr_2;
    // unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;

    if ( isFirstUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_1;
        tp_1 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_1 = correctPosition;
        chr_1 = occTranslate[approxValue].chrID - 1;
        // cx non-DP. check tp_1. PE
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    if ( isSecondUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_2;
        tp_2 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_2 = correctPosition;
        chr_2 = occTranslate[approxValue].chrID - 1;
        // cx non-DP. check tp_2. PE
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    unsigned int samFlag = 1;

    if ( properlyPaired == 1 )
    {
        samFlag |= SAM_FLAG_PROPER_PAIR;
    }

    if ( isFirstUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( isSecondUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_FIRST_IN_PAIR;

    if ( ( isFirstUnaligned == 0 ) && ( pePair->strand_1 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    if ( ( isSecondUnaligned == 0 ) && ( pePair->strand_2 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.n_cigar = 1;                                //SAM: number of CIGAR operations
    samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
    SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadLength << BAM_CIGAR_SHIFT ); //CIGAR

    for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
    {
        unsigned char biChar = 0;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    if ( qInfo->ReadLength % 2 == 1 )
    {
        uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
    {
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadQuality[i] );
    }

    unsigned int auxStart = samAlgnmt->data_len;
    //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
    //----------------------------------------------------------------------------------------------------------------------------

    //Outputing 1 alignment
    // ----------------------
    if ( isFirstUnaligned == 0 )
    {
        samAlgnmt->core.tid = chr_1;                          //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
    }
    else
    {
        samAlgnmt->core.tid = -1;                             //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                             //SAM: 0-based leftmost coordinate
    }

    samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
    samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag

    if ( isSecondUnaligned == 0 )
    {
        samAlgnmt->core.mtid = chr_2;
        samAlgnmt->core.mpos = tp_2 - 1;
    }
    else
    {
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
    }

    if ( ( isFirstUnaligned == 0 ) && ( isSecondUnaligned == 0 ) && ( chr_1 == chr_2 ) )
    {
        if ( tp_2 > tp_1 )
        { samAlgnmt->core.isize = tp_2 + mateLength - tp_1; }
        else
        { samAlgnmt->core.isize = - ( tp_1 + qInfo->ReadLength - tp_2 ); }
    }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );
    // ----------------------
}

void OCCDirectWritePairOccSAMAPI2 ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, int properlyPaired )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // output the pair based on the second read
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    // int OutFileFormat = qSetting->OutFileFormat;
    // FILE * outFilePtr = qSetting->OutFilePtr;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long i = 0; //, j=0, k=0;
    unsigned long long ambPosition;
    unsigned long long tp_1, tp_2;
    char chr_1, chr_2;
    // unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;

    if ( isSecondUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_2;
        tp_1 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_1 = correctPosition;
        chr_1 = occTranslate[approxValue].chrID - 1;
        // cx non-DP. check tp_1. PE
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    if ( isFirstUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_1;
        tp_2 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_2 = correctPosition;
        chr_2 = occTranslate[approxValue].chrID - 1;
        // cx non-DP. check tp_2. PE
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    unsigned int samFlag = 1;

    if ( properlyPaired == 1 )
    {
        samFlag |= SAM_FLAG_PROPER_PAIR;
    }

    if ( isSecondUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( isFirstUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_SECOND_IN_PAIR;

    if ( ( isFirstUnaligned == 0 ) && ( pePair->strand_1 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    if ( ( isSecondUnaligned == 0 ) && ( pePair->strand_2 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.n_cigar = 1;                                //SAM: number of CIGAR operations
    samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
    SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadLength << BAM_CIGAR_SHIFT ); //CIGAR

    for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
    {
        unsigned char biChar = 0;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    if ( qInfo->ReadLength % 2 == 1 )
    {
        uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
    {
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadQuality[i] );
    }

    unsigned int auxStart = samAlgnmt->data_len;
    //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
    //----------------------------------------------------------------------------------------------------------------------------

    //Outputing 1 alignment
    // ----------------------
    if ( isSecondUnaligned == 0 )
    {
        samAlgnmt->core.tid = chr_1;                          //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
    }
    else
    {
        samAlgnmt->core.tid = -1;                             //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                             //SAM: 0-based leftmost coordinate
    }

    samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
    samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag

    if ( isFirstUnaligned == 0 )
    {
        samAlgnmt->core.mtid = chr_2;
        samAlgnmt->core.mpos = tp_2 - 1;
    }
    else
    {
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
    }

    if ( ( isFirstUnaligned == 0 ) && ( isSecondUnaligned == 0 ) && ( chr_1 == chr_2 ) )
    {
        if ( tp_2 > tp_1 )
        { samAlgnmt->core.isize = tp_2 + mateLength - tp_1; }
        else
        { samAlgnmt->core.isize = - ( tp_1 + qInfo->ReadLength - tp_2 ); }
    }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );
    // ----------------------
}

void OCCDirectWritePairOccSAMAPIwCIGAR ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, char * cigar_str )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // output the pair based on the first read with CIGAR STRING
    SRASetting * qSetting = qInput->QuerySetting;
    // int OutFileFormat = qSetting->OutFileFormat;
    // FILE * outFilePtr = qSetting->OutFilePtr;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long i = 0; //, j=0, k=0;
    unsigned long long ambPosition;
    unsigned long long tp_1, tp_2;
    char chr_1, chr_2;
    // unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;

    if ( isFirstUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_1;
        tp_1 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_1 = correctPosition;
        chr_1 = occTranslate[approxValue].chrID - 1;
        // cx DP?   check tp_1. PE
        // note the cigar string passed in as paramater
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    if ( isSecondUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_2;
        tp_2 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_2 = correctPosition;
        chr_2 = occTranslate[approxValue].chrID - 1;
        // cx DP?   check tp_2. PE
        // note the cigar string passed in as paramater
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    unsigned int samFlag = 1;

    if ( ( isFirstUnaligned == 0 ) && ( isSecondUnaligned == 0 ) )
    {
        samFlag |= SAM_FLAG_PROPER_PAIR;
    }

    if ( isFirstUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( isSecondUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_FIRST_IN_PAIR;

    if ( ( isFirstUnaligned == 0 ) && ( pePair->strand_1 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    if ( ( isSecondUnaligned == 0 ) && ( pePair->strand_2 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
    AssignCigarStrToSAMIU ( samAlgnmt, & ( samAlgnmt->data_len ), cigar_str );

    for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
    {
        unsigned char biChar = 0;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    if ( qInfo->ReadLength % 2 == 1 )
    {
        uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
    {
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadQuality[i] );
    }

    unsigned int auxStart = samAlgnmt->data_len;
    //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
    //----------------------------------------------------------------------------------------------------------------------------

    //Outputing 1 alignment
    // ----------------------
    if ( isFirstUnaligned == 0 )
    {
        samAlgnmt->core.tid = chr_1;                          //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
    }
    else
    {
        samAlgnmt->core.tid = -1;                             //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                             //SAM: 0-based leftmost coordinate
    }

    samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
    samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag

    if ( isSecondUnaligned == 0 )
    {
        samAlgnmt->core.mtid = chr_2;
        samAlgnmt->core.mpos = tp_2 - 1;
    }
    else
    {
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
    }

    if ( ( isFirstUnaligned == 0 ) && ( isSecondUnaligned == 0 ) && ( chr_1 == chr_2 ) )
    {
        if ( tp_2 > tp_1 )
        { samAlgnmt->core.isize = tp_2 + mateLength - tp_1; }
        else
        { samAlgnmt->core.isize = - ( tp_1 + qInfo->ReadLength - tp_2 ); }
    }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );
    // ----------------------
}


void OCCDirectWritePairOccSAMAPI2wCIGAR ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, char * cigar_str )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // output the pair based on the second read with CIGAR STRING
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    // int OutFileFormat = qSetting->OutFileFormat;
    // FILE * outFilePtr = qSetting->OutFilePtr;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long i = 0; //, j=0; //, k=0;
    unsigned long long ambPosition;
    unsigned long long tp_1, tp_2;
    char chr_1, chr_2;
    // unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;

    if ( isSecondUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_2;
        tp_1 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_1 = correctPosition;
        chr_1 = occTranslate[approxValue].chrID - 1;
        // cx DP?   check tp_1. PE
        // note the cigar string passed in as paramater
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    if ( isFirstUnaligned == 0 )
    {
        ambPosition = pePair->algnmt_1;
        tp_2 = ambPosition;
        correctPosition = ambPosition;
        approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
        approxValue = occAmbiguityMap[approxIndex];

        while ( occTranslate[approxValue].startPos > ambPosition )
        {
            approxValue--;
        }

        correctPosition -= occTranslate[approxValue].correction;
        tp_2 = correctPosition;
        chr_2 = occTranslate[approxValue].chrID - 1;
        // cx DP?   check tp_2. PE
        // note the cigar string passed in as paramater
        // have HSP
        // need to use SAMAPI to modify CIGAR string
    }

    unsigned int samFlag = 1;

    if ( ( isFirstUnaligned == 0 ) && ( isSecondUnaligned == 0 ) )
    {
        samFlag |= SAM_FLAG_PROPER_PAIR;
    }

    if ( isSecondUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( isFirstUnaligned == 1 )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_SECOND_IN_PAIR;

    if ( ( isFirstUnaligned == 0 ) && ( pePair->strand_1 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    if ( ( isSecondUnaligned == 0 ) && ( pePair->strand_2 == QUERY_NEG_STRAND ) )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
    AssignCigarStrToSAMIU ( samAlgnmt, & ( samAlgnmt->data_len ), cigar_str );

    for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
    {
        unsigned char biChar = 0;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    if ( qInfo->ReadLength % 2 == 1 )
    {
        uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
    {
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadQuality[i] );
    }

    unsigned int auxStart = samAlgnmt->data_len;
    //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
    //----------------------------------------------------------------------------------------------------------------------------

    //Outputing 1 alignment
    // ----------------------
    if ( isSecondUnaligned == 0 )
    {
        samAlgnmt->core.tid = chr_1;                          //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
    }
    else
    {
        samAlgnmt->core.tid = -1;                             //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                             //SAM: 0-based leftmost coordinate
    }

    samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
    samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag

    if ( isFirstUnaligned == 0 )
    {
        samAlgnmt->core.mtid = chr_2;
        samAlgnmt->core.mpos = tp_2 - 1;
    }
    else
    {
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
    }

    if ( ( isFirstUnaligned == 0 ) && ( isSecondUnaligned == 0 ) && ( chr_1 == chr_2 ) )
    {
        if ( tp_2 > tp_1 )
        { samAlgnmt->core.isize = tp_2 + mateLength - tp_1; }
        else
        { samAlgnmt->core.isize = - ( tp_1 + qInfo->ReadLength - tp_2 ); }
    }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );
    // ----------------------
}


unsigned int getChrAndPos ( SRAQueryInput * qInput, unsigned long long ambPos,
                            unsigned long long * tp, unsigned short * chr_id )
{
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    // get the chromosome and the position for the position on the packed sequence
    HSP * hsp = aIndex->hsp;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    Translate * occTranslate = hsp->translate;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;
    correctPosition = ambPos;
    approxIndex = ambPos >> GRID_SAMPLING_FACTOR_2_POWER;
    approxValue = occAmbiguityMap[approxIndex];
    // printf ( "delay no more!\n");
    // printf ( "%u %u\n", approxIndex, approxValue );
    //for ( int i=0;i<approxIndex; ++i ) 
    //   {  printf ( "%u\n", occAmbiguityMap[i] ); }
    while ( occTranslate[approxValue].startPos > ambPos )
    {
        approxValue--;
    }

    correctPosition -= occTranslate[approxValue].correction;
    ( *tp ) = correctPosition;
    ( *chr_id ) = occTranslate[approxValue].chrID;

    return approxValue < hsp->numOfRemovedSegment - 1 ?
      occTranslate[approxValue + 1].startPos - 1 : hsp->dnaLength;

    //    printf("[getChrAndPos] %u %u %u\n", ambPos, occTranslate[approxValue].startPos, occTranslate[approxValue+1].startPos );
}


int BoundaryCheck ( unsigned int pacPos, int chrID, unsigned int chrEndPos, int readLength, unsigned int segmentEndPos,
                    HSP * hsp, unsigned int & correctedPac, char ** newCIGAR )
{
    segmentEndPos = chrEndPos < segmentEndPos ? chrEndPos : segmentEndPos;

    if ( pacPos + readLength <= segmentEndPos + 1 ) { return 0; }

    int actualAlignedLength = segmentEndPos - pacPos + 1;

    if ( newCIGAR ) { *newCIGAR = ( char * ) malloc ( 15 * sizeof ( char ) ); }
    // printf ( "%d %u %u %u\n", actualAlignedLength, segmentEndPos, pacPos, chrEndPos );
    if ( actualAlignedLength >= ( readLength + 1 ) / 2 )
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%dM%dS", actualAlignedLength, readLength - actualAlignedLength ); }

        correctedPac = pacPos;
        return - ( readLength - actualAlignedLength ); // trim right
    }
    else
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%dS%dM", actualAlignedLength, readLength - actualAlignedLength ); }

        correctedPac = segmentEndPos + 1;
        return actualAlignedLength; // trim left
    }
}

// returns trimmed amount (bp). +ve indicates trimming on left side, -ve indicates right side.
int BoundaryCheckDP ( unsigned int pacPos, int chrID, unsigned int chrEndPos, int readLength, char * cigar, unsigned int segmentEndPos, HSP * hsp,
                      unsigned int & correctedPac, char ** newCIGAR )
{
    if ( pacPos + readLength * 2 <= chrEndPos + 1 && pacPos + readLength * 2 <= segmentEndPos + 1 ) { return 0; } // is this readLength * 2 bound good enough?

    segmentEndPos = chrEndPos < segmentEndPos ? chrEndPos : segmentEndPos;

    // CIGAR string parsing...
    char * buffer = ( char * ) malloc ( ( strlen ( cigar ) + 10 ) * sizeof ( char ) );
    char * leftBuf = ( char * ) malloc ( ( strlen ( cigar ) + 10 ) * sizeof ( char ) );
    char * rightBuf = ( char * ) malloc ( ( strlen ( cigar ) + 10 ) * sizeof ( char ) );
    int leftLen, rightLen;
    int leftOff = 0, rightOff = 0;
    int leftS, rightS;
    char * cigarP = cigar;
    {
        unsigned int refPos;
        refPos = pacPos;
        leftLen = rightLen = 0;
        leftBuf[0] = rightBuf[0] = 0;
        leftS = rightS = 0;

        while ( *cigarP )
        {
            int num = 0;
            char op;

            while ( *cigarP <= '9' )
            {
                num = num * 10 + ( *cigarP - '0' );
                cigarP++;
            }

            op = *cigarP;
            cigarP++;

            if ( op == 'S' )
            {
                // refPos += num;
                sprintf ( buffer, "%dS", num );

                if ( refPos <= segmentEndPos )
                {
                    strcat ( leftBuf, buffer );
                    leftS += num;
                }
                else
                {
                    strcat ( rightBuf, buffer );
                    rightS += num;
                }
            }
            else if ( op == 'M' || op == 'm' )
            {
                if ( refPos > segmentEndPos )
                {
                    rightLen += num;
                    sprintf ( buffer, "%d%c", num, op );
                    strcat ( rightBuf, buffer );
                }
                else if ( refPos + num <= segmentEndPos + 1 )
                {
                    leftLen += num;
                    sprintf ( buffer, "%d%c", num, op );
                    strcat ( leftBuf, buffer );
                }
                else     // crossing boundary
                {
                    leftLen += segmentEndPos - refPos + 1;
                    sprintf ( buffer, "%d%c", segmentEndPos - refPos + 1, op );
                    strcat ( leftBuf, buffer );
                    rightLen += num - segmentEndPos + refPos - 1;
                    sprintf ( buffer, "%d%c", num - segmentEndPos + refPos - 1, op );
                    strcat ( rightBuf, buffer );
                }

                refPos += num;
            }
            else if ( op == 'D' )
            {
                sprintf ( buffer, "%dD", num );

                if ( refPos > segmentEndPos )
                {
                    if ( strlen ( rightBuf ) == 0 )
                    {
                        rightOff = num;
                    }
                    else
                    {
                        strcat ( rightBuf, buffer );
                    }
                }
                else if ( refPos + num <= segmentEndPos + 1 )
                {
                    strcat ( leftBuf, buffer );
                }
                else     // crossing boundary
                {
                    sprintf ( buffer, "%dD", segmentEndPos - refPos + 1 );
                    strcat ( leftBuf, buffer );
                    sprintf ( buffer, "%dD", num - segmentEndPos + refPos - 1 );

                    if ( strlen ( rightBuf ) == 0 )
                    {
                        rightOff = num - segmentEndPos + refPos - 1;
                    }
                    else
                    {
                        strcat ( rightBuf, buffer );
                    }
                }

                refPos += num;
            }
            else if ( op == 'I' )
            {
                if ( refPos == chrEndPos + 1 ) // exactly at boundary
                {
                    // ignore both sides
                }
                else if ( refPos <= segmentEndPos )
                {
                    leftLen += num;
                    sprintf ( buffer, "%dI", num );
                    strcat ( leftBuf, buffer );
                }
                else
                {
                    rightLen += num;
                    sprintf ( buffer, "%dI", num );
                    strcat ( rightBuf, buffer );
                }
            }
        }
    }

    if ( newCIGAR ) { *newCIGAR = buffer; }

    if ( !leftLen || !rightLen )
    {
        free ( leftBuf );
        free ( rightBuf );

        if ( newCIGAR )
        {
            free ( *newCIGAR );
            *newCIGAR = NULL;
        }

        return 0; // no trimming
    }
    else if ( leftLen >= rightLen )
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%s%dS", leftBuf, readLength - leftLen - leftS ); }

        correctedPac = pacPos;
        free ( leftBuf );
        free ( rightBuf );
        return - ( readLength - leftLen - leftS ); /// trim right
    }
    else
    {
        if ( newCIGAR ) { sprintf ( *newCIGAR, "%dS%s", readLength - rightLen - rightS, rightBuf ); }

        correctedPac = segmentEndPos + 1 + rightOff;
        free ( leftBuf );
        free ( rightBuf );
        return readLength - rightLen - rightS; // trim left
    }

}

// returns trimmed amount (bp). +ve indicates trimming on left side, -ve indicates right side.
int getChrAndPosWithBoundaryCheck ( SRAQueryInput * qInput, unsigned int readLength, unsigned long long ambPos,
                                    unsigned long long * tp, unsigned short * chr_id, char ** buffer )
{
    unsigned int segmentEndPos = getChrAndPos ( qInput, ambPos, tp, chr_id );

    HSP * hsp = qInput->AlgnmtIndex->hsp;
    unsigned int chrEndPos = hsp->seqOffset[ *chr_id - 1 ].endPos;
    unsigned int correctedPac;
    int ret;

    if ( ret = BoundaryCheck ( ambPos, *chr_id, chrEndPos, readLength, segmentEndPos, hsp, correctedPac, buffer ) )
    {
        if ( correctedPac > ambPos ) // trimmed left side
        {
            getChrAndPos ( qInput, correctedPac, tp, chr_id );
        }
    }

    /*  if (ret) {
      printf("BAD ALIGNMENT %u\n", ambPos);
      for (int i=0;i<readLength;++i)
        printf("%c", dnaChar[qInput->QueryInfo->ReadCode[i]]);
      printf("\n");
      }*/
    return ret;
}

// returns trimmed amount (bp). +ve indicates trimming on left side, -ve indicates right side.
int getChrAndPosWithBoundaryCheckDP ( SRAQueryInput * qInput, unsigned int readLength, unsigned long long ambPos, char * cigar,
                                      unsigned long long * tp, unsigned short * chr_id, char ** buffer )
{
    unsigned int segmentEndPos = getChrAndPos ( qInput, ambPos, tp, chr_id );
    HSP * hsp = qInput->AlgnmtIndex->hsp;
    unsigned int chrEndPos = hsp->seqOffset[ *chr_id - 1 ].endPos;
    unsigned int correctedPac;
    int ret;

    if ( ret = BoundaryCheckDP ( ambPos, *chr_id, chrEndPos, readLength, cigar, segmentEndPos, hsp, correctedPac, buffer ) )
    {
        if ( correctedPac > ambPos ) // trimmed left side
        {
            getChrAndPos ( qInput, correctedPac, tp, chr_id );
        }
    }

    /*  if (ret) {
      printf("BAD ALIGNMENT DP %u\n", ambPos);
      for (int i=0;i<readLength;++i)
        printf("%c", dnaChar[qInput->QueryInfo->ReadCode[i]]);
      printf("\n");
      }*/

    return ret;
}

void initializeSAMAlgnmt ( bam1_t * samAlgnmt, int readlen, char * queryName, unsigned char * querySeq,
                           char * qualities, char strand, uint8_t * xazStr, int xazlen, char * cigar_str,
                           int isUnmapped, char * readGroup )
{
    int i;
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );      //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.l_qseq = readlen;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( queryName ) + 1; //SAM: length of the query name

    if ( isUnmapped )
    { samAlgnmt->core.qual = 0; }     //SAM: mapping quality
    else
    { samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE; }     //SAM: mapping quality

    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), queryName, strlen ( queryName ) + 1 ); //Name

    if ( isUnmapped )
    {
        samAlgnmt->core.n_cigar = 0;                     //SAM: number of CIGAR operations
    }
    else
    {
        if ( cigar_str != NULL )
        {
            AssignCigarStrToSAMIU ( samAlgnmt, & ( samAlgnmt->data_len ), cigar_str );
        }
        else
        {
            samAlgnmt->core.n_cigar = 1;                     //SAM: number of CIGAR operations
            SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), readlen << BAM_CIGAR_SHIFT ); //CIGAR
        }
    }

    if ( strand == QUERY_NEG_STRAND )
    {
        // for negative strand
        if ( readlen % 2 == 1 )
        {
            for ( i = ( readlen - 1 ) / 2; i > 0 ; i-- )                                                    //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 - 1]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            uint8_t biChar = bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[0]]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }
        else
        {
            for ( i = readlen / 2 - 1; i >= 0 ; i-- )                                                      //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 + 1]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }
        }

        for ( i = readlen - 1; i >= 0 ; i-- )                                                          //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] );
        }
    }
    else
    {
        // for positive strand
        for ( i = 0; i < readlen / 2; i++ )                                                        //Read
        {
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2]]] << 4;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2 + 1]]];
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        if ( readlen % 2 == 1 )
        {
            uint8_t biChar = bam_nt16_table[dnaChar[querySeq[readlen - 1]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        for ( i = 0; i < readlen; i++ )                                                            //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] );
        }
    }

    unsigned int auxStart = samAlgnmt->data_len;
    bam_aux_append ( samAlgnmt, "RG", 'Z', strlen ( readGroup ) + 1, ( uint8_t * ) readGroup );

    if ( xazlen > 0 )
    { bam_aux_append ( samAlgnmt, "XA", 'Z', xazlen + 1, xazStr ); }

    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
}

void initializeSAMAlgnmt2 ( bam1_t * samAlgnmt, int readlen, char * queryName, unsigned char * querySeq,
                            char * qualities, char strand, uint8_t * xazStr, int xazlen, char * cigar_str, int isUnmapped,
                            int mismatchNum, int editDist, int bestHitNum, int secBestHitNum, int gapOpenNum, int gapExtendNum,
                            char * mdStr, int mdlen, int map_qual_score, char * readGroup, bool isPrintMDNM )
{
    // fprintf(stderr, "bestHitNum : %d; secBestHitNum : %d \n", bestHitNum, secBestHitNum);
    // also show the mismatch information

    int i;
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );      //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.l_qseq = readlen;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( queryName ) + 1; //SAM: length of the query name

    if ( isUnmapped )
    { samAlgnmt->core.qual = 0; }     //SAM: mapping quality
    else
    { samAlgnmt->core.qual = map_qual_score; }     //SAM: mapping quality

    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), queryName, strlen ( queryName ) + 1 ); //Name

    if ( isUnmapped )
    {
        samAlgnmt->core.n_cigar = 0;                     //SAM: number of CIGAR operations
    }
    else
    {
        if ( cigar_str != NULL )
        {
            AssignCigarStrToSAMIU ( samAlgnmt, & ( samAlgnmt->data_len ), cigar_str );
        }
        else
        {
            samAlgnmt->core.n_cigar = 1;                     //SAM: number of CIGAR operations
            SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), readlen << BAM_CIGAR_SHIFT ); //CIGAR
        }
    }

    if ( strand == QUERY_NEG_STRAND )
    {
        // for negative strand
        if ( readlen % 2 == 1 )
        {
            for ( i = ( readlen - 1 ) / 2; i > 0 ; i-- )                                                    //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 - 1]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }

            uint8_t biChar = bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[0]]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }
        else
        {
            for ( i = readlen / 2 - 1; i >= 0 ; i-- )                                                      //Read
            {
                unsigned char biChar = 0;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2 + 1]]]] << 4;
                biChar |= bam_nt16_table[dnaChar[soap3DnaComplement[querySeq[i * 2]]]];
                SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
            }
        }

        for ( i = readlen - 1; i >= 0 ; i-- )                                                          //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] );
        }
    }
    else
    {
        // for positive strand
        for ( i = 0; i < readlen / 2; i++ )                                                        //Read
        {
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2]]] << 4;
            biChar |= bam_nt16_table[dnaChar[querySeq[i * 2 + 1]]];
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        if ( readlen % 2 == 1 )
        {
            uint8_t biChar = bam_nt16_table[dnaChar[querySeq[readlen - 1]]] << 4;
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
        }

        for ( i = 0; i < readlen; i++ )                                                            //Quality
        {
            SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qualities[i] );
        }
    }

    unsigned int auxStart = samAlgnmt->data_len;
    bam_aux_append ( samAlgnmt, "RG", 'Z', strlen ( readGroup ) + 1, ( uint8_t * ) readGroup );

    if ( !isUnmapped )
    {
        if ( isPrintMDNM && editDist >= 0 )
        {
            bam_aux_append ( samAlgnmt, "NM", 'i', 4, ( uint8_t * ) &editDist );
        }

        if ( bestHitNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "X0", 'i', 4, ( uint8_t * ) &bestHitNum );
        }

        if ( secBestHitNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "X1", 'i', 4, ( uint8_t * ) &secBestHitNum );
        }

        if ( mismatchNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "XM", 'i', 4, ( uint8_t * ) &mismatchNum );
        }

        if ( gapOpenNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "XO", 'i', 4, ( uint8_t * ) &gapOpenNum );
        }

        if ( gapExtendNum >= 0 )
        {
            bam_aux_append ( samAlgnmt, "XG", 'i', 4, ( uint8_t * ) &gapExtendNum );
        }

        if ( isPrintMDNM && mdlen > 0 )
        {
            bam_aux_append ( samAlgnmt, "MD", 'Z', mdlen + 1, ( uint8_t * ) mdStr );
        }

        if ( xazlen > 0 )
        {
            bam_aux_append ( samAlgnmt, "XA", 'Z', xazlen + 1, xazStr );
        }
    }

    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
}

int getMapQualScore ( int n, int mismatchNum, int avgMismatchQual, int maxMAPQ, int minMAPQ )
{
    // not considering the second best
    if ( n == 1 )
    {
        // n is 1
        int mismatches_index = mismatchNum;

        if ( mismatches_index > 5 )
        { mismatches_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[mismatches_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}

int bwaLikeSingleQualScore ( int x0, int x1, int * g_log_n )
{
    // fprintf(stderr, "x0 : %i; x1 : %i\n", x0, x1);
    int score;

    if ( x0 > 1 )
    { score = 0; }
    else if ( x1 == 0 )
    { score = 37; }
    else
    {
        x1 = ( x1 > 255 ) ? 255 : x1;
        int n = g_log_n[x1];
        score = ( 23 < n ) ? 0 : 23 - n;
    }

    // fprintf(stderr, "score : %i\n", score);
    return score;
}

int getMapQualScoreSingle ( int mismatchNum, int avgMismatchQual, int x0, int x1, int maxMAPQ, int minMAPQ, int isBWALike, int * g_log_n )
{
    if ( isBWALike )
    {
        return bwaLikeSingleQualScore ( x0, x1, g_log_n );
    }

    if ( x0 == 1 )
    {
        if ( x1 > 0 )
        {
            return minMAPQ;
        }

        int mismatches_index = mismatchNum;

        if ( mismatches_index > 5 )
        { mismatches_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[mismatches_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}

int getMapQualScoreForSingleDP ( int maxDPScore, int avgMismatchQual, int x0, int x1_t1, int x1_t2, int bestDPScore, int secondBestDPScore, int maxMAPQ, int minMAPQ, int dpThres, int isBWALike, int * g_log_n )
{
    if ( isBWALike )
    {
        return bwaLikeSingleQualScore ( x0, x1_t1 + x1_t2, g_log_n );
    }

    // considersation of uniqueness
    if ( x0 > 1 || x1_t1 > 0 )
    { return minMAPQ; }

    float R1, R2, R3, P;

    // consideration of suboptimal score (i.e. secondBestDPScore)
    if ( x1_t2 > 0 )
    {
        R1 = 1.0 - ( ( float ) ( secondBestDPScore - dpThres ) ) / ( 0.7 * bestDPScore - dpThres );
    }
    else
    {
        R1 = 1.0;
    };

    // consideration of X1 (i.e. x1_t1 + x1_t2)
    int x1 = x1_t1 + x1_t2;

    R2 = ( x1 > 100 ) ? penalty_ratio_x1[100] : penalty_ratio_x1[x1];

    // consideration of best DP score
    R3 = ( ( float ) ( bestDPScore - dpThres ) ) / ( maxDPScore - dpThres );

    // consideration of average of quality value in mismatch positions
    if ( avgMismatchQual < 0 ) { avgMismatchQual = 0; }
    else if ( avgMismatchQual > 40 ) { avgMismatchQual = 40; }

    P = penalty_score_avg_mis_qual[avgMismatchQual];
    // calculation of the MAPQ score
    int mapqScore = ( int ) ( maxMAPQ * R1 * R2 * R3 - P );

    if ( mapqScore < minMAPQ ) { mapqScore = minMAPQ; }

    return mapqScore;
}


void bwaLikePairQualScore ( int x0_0, int x1_0, int x0_1, int x1_1, int * g_log_n, int op_score, int op_num, int subop_score, int subop_num, int readlen_0, int readlen_1, int * map_score0, int * map_score1 )
{
    // fprintf(stderr, "x0_0:%i x1_0:%i x0_1:%i x1_1:%i op_score:%i op_num:%i subop_score:%i subop_num:%i readlen_0:%i readlen_1:%i \n", x0_0, x1_0, x0_1, x1_1, op_score, op_num, subop_score, subop_num, readlen_0, readlen_1);
    int mapq0 = bwaLikeSingleQualScore ( x0_0, x1_0, g_log_n );
    int mapq1 = bwaLikeSingleQualScore ( x0_1, x1_1, g_log_n );
    // fprintf(stderr, "mapq0 : %i ; mapq1 : %i\n", mapq0, mapq1);
    op_score = op_score * 10;
    subop_score = subop_score * 10;
    int mapq_p = 0;

    if ( mapq0 > 0 && mapq1 > 0 )
    {
        mapq_p = mapq0 + mapq1;

        if ( mapq_p > 60 ) { mapq_p = 60; }

        mapq0 = mapq_p;
        mapq1 = mapq_p;
    }
    else
    {
        if ( op_num == 1 )
        {
            if ( subop_num == 0 )
            { mapq_p = 29; }
            else if ( op_score - subop_score > ( 0.3 * ( ( readlen_0 + readlen_1 ) / 2 ) ) ) { mapq_p = 23; }
            else
            {
                subop_num = ( subop_num > 255 ) ? 255 : subop_num;
                mapq_p = ( op_score - subop_score ) / 2 - g_log_n[subop_num];

                if ( mapq_p < 0 ) { mapq_p = 0; }
            }
        }

        if ( mapq0 == 0 )
        {
            mapq0 = ( mapq_p + 7 < mapq1 ) ? mapq_p + 7 : mapq1;
        }

        if ( mapq1 == 0 )
        {
            mapq1 = ( mapq_p + 7 < mapq0 ) ? mapq_p + 7 : mapq0;
        }
    }

    ( *map_score0 ) = mapq0;
    ( *map_score1 ) = mapq1;
}

int getMapQualScore2 ( int mismatchNum, int avgMismatchQual, int x0, int x1, char isBestHit, unsigned int totalNumValidPairs, int maxMAPQ, int minMAPQ )
{
    if ( x0 == 1 && totalNumValidPairs == 1 )
    {
        if ( isBestHit == 0 && ( x1 > 1 ) )
        {
            return minMAPQ;
        }

        int mismatches_index = mismatchNum;

        if ( mismatches_index > 5 )
        { mismatches_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[mismatches_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}


int getMapQualScoreForDP ( int n, int dpScore, int maxDPScore, int avgMismatchQual, int maxMAPQ, int minMAPQ )
{
    // does not consider second best
    if ( n == 1 )
    {
        // n = 1
        int difference_index = 0;

        if ( dpScore < maxDPScore )
        { difference_index = ( int ) ( ( 1.0 - ( double ) dpScore / maxDPScore ) * 100.0 - 1.0 ) / 5 + 1; }

        if ( difference_index > 5 )
        { difference_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[difference_index][avg_mismatch_qual_index] );

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}

int getMapQualScoreForDP2 ( int dpScore, int maxDPScore, int avgMismatchQual, int x0, int x1, int bestDPScore, int secondBestDPScore, char isBestHit, int totalNumValidPairs, int maxMAPQ, int minMAPQ )
{
    if ( x0 == 1 && totalNumValidPairs == 1 )
    {
        if ( isBestHit == 0 && ( x1 > 1 ) )
        {
            return minMAPQ;
        }

        int difference_index = 0;

        if ( dpScore < maxDPScore )
        {
            difference_index = ( int ) ( ( 1.0 - ( double ) dpScore / maxDPScore ) * 100.0 - 1.0 ) / 4 + 1;
        }

        if ( difference_index > 5 )
        { difference_index = 5; }

        int avg_mismatch_qual_index = ( avgMismatchQual - 1 ) / 20;

        if ( avg_mismatch_qual_index > 1 )
        { avg_mismatch_qual_index = 1; }
        else if ( avg_mismatch_qual_index < 0 )
        { avg_mismatch_qual_index = 0; }

        int mapQualScore = ( int ) ( maxMAPQ * mapping_score[difference_index][avg_mismatch_qual_index] );

        if ( bestDPScore > secondBestDPScore && ( ( double ) bestDPScore - secondBestDPScore ) / maxDPScore < 0.2 )
        { mapQualScore = minMAPQ; }

        if ( mapQualScore < minMAPQ )
        { mapQualScore = minMAPQ; }

        return mapQualScore;
    }
    else
    {
        return minMAPQ;
    }
}


int getMapQualScoreForPair ( int score1, int score2 )
{
    return ( score1 > score2 ) ? ( int ) ( score1 * 0.2 + score2 * 0.8 ) : ( int ) ( score1 * 0.8 + score2 * 0.2 );
}

void unproperlypairOutputSAMAPI ( SRAQueryInput * qInput, OCCList * occ_list1, OCCList * occ_list2,
                                  unsigned char * query1, unsigned char * query2,
                                  char * qualities1, char * qualities2,
                                  int readlen1, int readlen2,
                                  char * queryName1, char * queryName2,
                                  DynamicUint8Array * xazArray,
                                  unsigned int peMaxOutputPerRead, int reportType )
{
    // cx this function... X_x

    //  fprintf(stderr,"unproperlypairOutputSAMAPI\n");

    // output the alignments which are not properly paired
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long tp_1 = 0;
    unsigned long long tp_2 = 0;
    unsigned long long curr_tp = 0;
    unsigned short chr_1 = 0;
    unsigned short chr_2 = 0;
    unsigned short curr_chr = 0;
    unsigned int samFlag;
    unsigned int i;
    char curr_strand;
    int curr_len = 0;
    unsigned int curr_numMisMatch;
    char curr_occStr[500];
    unsigned int total_results;
    char mdStr1[MAX_READ_LENGTH * 2 + 1];
    char mdStr2[MAX_READ_LENGTH * 2 + 1];
    int mdLen1 = 0;
    int mdLen2 = 0;
    int avg_mismatch_qual1;
    int avg_mismatch_qual2;
    int map_qual_score1 = 0;
    int map_qual_score2 = 0;
    // obtain the best occurence
    SRAOccurrence * bestOcc1 = NULL;
    int bestScore1;
    int bestScoreNum1 = 0;
    int secBestScoreNum1 = 0;

    if ( occ_list1->curr_size > 0 )
    {
        bestOcc1 = & ( occ_list1->occ[0] );
        bestScore1 = bestOcc1->mismatchCount;
        bestScoreNum1 = 1;

        for ( i = 1; i < occ_list1->curr_size; i++ )
        {
            if ( occ_list1->occ[i].mismatchCount < bestScore1 )
            {
                if ( bestScore1 == occ_list1->occ[i].mismatchCount + 1 )
                {
                    secBestScoreNum1 = bestScoreNum1;
                }
                else
                {
                    secBestScoreNum1 = 0;
                }

                bestOcc1 = & ( occ_list1->occ[i] );
                bestScore1 = bestOcc1->mismatchCount;
                bestScoreNum1 = 1;
            }
            else if ( occ_list1->occ[i].mismatchCount == bestScore1 )
            {
                bestScoreNum1++;
            }
            else if ( occ_list1->occ[i].mismatchCount == bestScore1 + 1 )
            {
                secBestScoreNum1++;
            }
        }
    }

    SRAOccurrence * bestOcc2 = NULL;
    int bestScore2;
    int bestScoreNum2 = 0;
    int secBestScoreNum2 = 0;

    if ( occ_list2->curr_size > 0 )
    {
        bestOcc2 = & ( occ_list2->occ[0] );
        bestScore2 = bestOcc2->mismatchCount;
        bestScoreNum2 = 1;
    }

    for ( i = 1; i < occ_list2->curr_size; i++ )
    {
        if ( occ_list2->occ[i].mismatchCount < bestScore2 )
        {
            if ( bestScore2 == occ_list2->occ[i].mismatchCount + 1 )
            {
                secBestScoreNum2 = bestScoreNum2;
            }
            else
            {
                secBestScoreNum2 = 0;
            }

            bestOcc2 = & ( occ_list2->occ[i] );
            bestScore2 = bestOcc2->mismatchCount;
            bestScoreNum2 = 1;
        }
        else if ( occ_list2->occ[i].mismatchCount == bestScore2 )
        {
            bestScoreNum2++;
        }
        else if ( occ_list2->occ[i].mismatchCount == bestScore2 + 1 )
        {
            secBestScoreNum2++;
        }
    }

    if ( ( bestOcc1 != NULL ) && ( reportType != OUTPUT_UNIQUE_BEST || bestScoreNum1 == 1 ) )
    {
        getChrAndPos ( qInput, bestOcc1->ambPosition,
                       &tp_1, &chr_1 );
        mdLen1 = getMdStr ( hsp, query1, qualities1, readlen1, bestOcc1->ambPosition, bestOcc1->strand, bestOcc1->mismatchCount, mdStr1, &avg_mismatch_qual1 );

        if ( reportType == OUTPUT_RANDOM_BEST || reportType == OUTPUT_UNIQUE_BEST )
        {
            map_qual_score1 = SAM_MAPQ_UNAVAILABLE;
        }
        else
        {

            map_qual_score1 = getMapQualScoreSingle ( bestOcc1->mismatchCount, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, bestScoreNum1, secBestScoreNum1, hspaux->maxMAPQ, hspaux->minMAPQ, hspaux->bwaLikeScore, hspaux->g_log_n );
            map_qual_score1 = map_qual_score1 >> 1;

            if ( map_qual_score1 < hspaux->minMAPQ )
            {
                map_qual_score1 = hspaux->minMAPQ;
            }
        }
    }
    else
    {
        bestOcc1 = NULL;
    }

    if ( ( bestOcc2 != NULL ) && ( reportType != OUTPUT_UNIQUE_BEST || bestScoreNum2 == 1 ) )
    {
        getChrAndPos ( qInput, bestOcc2->ambPosition,
                       &tp_2, &chr_2 );
        mdLen2 = getMdStr ( hsp, query2, qualities2, readlen2, bestOcc2->ambPosition, bestOcc2->strand, bestOcc2->mismatchCount, mdStr2, &avg_mismatch_qual2 );

        if ( reportType == OUTPUT_RANDOM_BEST || reportType == OUTPUT_UNIQUE_BEST )
        {
            map_qual_score2 = SAM_MAPQ_UNAVAILABLE;
        }
        else
        {
            map_qual_score2 = getMapQualScoreSingle ( bestOcc2->mismatchCount, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, bestScoreNum2, secBestScoreNum2, hspaux->maxMAPQ, hspaux->minMAPQ, hspaux->bwaLikeScore, hspaux->g_log_n );
            map_qual_score2 = map_qual_score2 >> 1;

            if ( map_qual_score2 < hspaux->minMAPQ )
            {
                map_qual_score2 = hspaux->minMAPQ;
            }
        }
    }
    else
    {
        bestOcc2 = NULL;
    }

    //-------------------------------------------//
    // report the first alignment                //
    //-------------------------------------------//
    if ( reportType == OUTPUT_ALL_BEST || reportType == OUTPUT_ALL_VALID )
    {
        // Build the XA:Z tag
        DynamicUint8ArrayReset ( xazArray );
        total_results = 1; // including the best result

        for ( i = 0; i < occ_list1->curr_size && total_results < peMaxOutputPerRead; i++ )
        {
            SRAOccurrence * currOcc = & ( occ_list1->occ[i] );

            if ( bestOcc1 == currOcc )
            { continue; }

            curr_numMisMatch = currOcc->mismatchCount;

            if ( reportType == OUTPUT_ALL_BEST && curr_numMisMatch > bestScore1 )
            { continue; }

            curr_strand = currOcc->strand;
            getChrAndPos ( qInput, currOcc->ambPosition, &curr_tp, &curr_chr );
            curr_len = sprintf ( curr_occStr, "%s,%c%llu,%iM,%i;",
                                 samFilePtr->header->target_name[curr_chr - 1], // chromosome name
                                 ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+', // strand
                                 curr_tp,           // chromosome position
                                 readlen1,
                                 curr_numMisMatch ); // number of mismatches
            appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            total_results++;
        }
    }

    if ( bestOcc1 != NULL )
        initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1,
                               qualities1, bestOcc1->strand, xazArray->charStr, xazArray->length, NULL, 0,
                               bestOcc1->mismatchCount, bestOcc1->mismatchCount,
                               ( reportType == OUTPUT_RANDOM_BEST ) ? -1 : bestScoreNum1,
                               ( reportType == OUTPUT_ALL_VALID || reportType == OUTPUT_ALL_BEST ) ? secBestScoreNum1 : -1, 0, 0, mdStr1, mdLen1, map_qual_score1,
                               hspaux->readGroup, hspaux->isPrintMDNM );
    else
        initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1,
                              qualities1, QUERY_POS_STRAND, xazArray->charStr,
                              xazArray->length, NULL, 1, hspaux->readGroup );

    // compute the value of the flag
    samFlag = 1;

    if ( bestOcc1 == NULL )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( bestOcc2 == NULL )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_FIRST_IN_PAIR;

    if ( bestOcc1 != NULL && bestOcc1->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    if ( bestOcc2 != NULL && bestOcc2->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
    samAlgnmt->core.tid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;    //SAM: chromosome ID, defined by bam_header_t
    samAlgnmt->core.pos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;     //SAM: 0-based leftmost coordinate
    samAlgnmt->core.mtid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;
    samAlgnmt->core.mpos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;

    if ( chr_1 > 0 && chr_1 == chr_2 )
        if ( tp_2 > tp_1 )
        { samAlgnmt->core.isize = tp_2 + readlen2 - tp_1; }
        else
        { samAlgnmt->core.isize = - ( tp_1 + readlen1 - tp_2 ); }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );

    //--------------------------------------------//
    // report the second alignment                //
    //--------------------------------------------//
    if ( reportType == OUTPUT_ALL_BEST || reportType == OUTPUT_ALL_VALID )
    {
        // Build the XA:Z tag
        DynamicUint8ArrayReset ( xazArray );
        total_results = 1; // including the best result

        for ( i = 0; i < occ_list2->curr_size && total_results < peMaxOutputPerRead; i++ )
        {
            SRAOccurrence * currOcc = & ( occ_list2->occ[i] );

            if ( bestOcc2 == currOcc )
            { continue; }

            curr_numMisMatch = currOcc->mismatchCount;

            if ( reportType == OUTPUT_ALL_BEST && curr_numMisMatch > bestScore2 )
            { continue; }

            curr_strand = currOcc->strand;
            getChrAndPos ( qInput, currOcc->ambPosition, &curr_tp, &curr_chr );
            curr_len = sprintf ( curr_occStr, "%s,%c%llu,%iM,%i;",
                                 samFilePtr->header->target_name[curr_chr - 1], // chromosome name
                                 ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+', // strand
                                 curr_tp,           // chromosome position
                                 readlen2,
                                 curr_numMisMatch ); // number of mismatches
            appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            total_results++;
        }
    }

    if ( bestOcc2 != NULL )
        initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2,
                               qualities2, bestOcc2->strand, xazArray->charStr, xazArray->length, NULL, 0,
                               bestOcc2->mismatchCount, bestOcc2->mismatchCount,
                               ( reportType == OUTPUT_RANDOM_BEST ) ? -1 : bestScoreNum2,
                               ( reportType == OUTPUT_ALL_VALID || reportType == OUTPUT_ALL_BEST ) ? secBestScoreNum2 : -1, 0, 0, mdStr2, mdLen2,
                               ( reportType == OUTPUT_RANDOM_BEST ) ? SAM_MAPQ_UNAVAILABLE : map_qual_score2,
                               hspaux->readGroup, hspaux->isPrintMDNM );
    else
        initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2,
                              qualities2, QUERY_POS_STRAND, xazArray->charStr,
                              xazArray->length, NULL, 1, hspaux->readGroup );

    // compute the value of the flag
    samFlag = 1;

    if ( bestOcc2 == NULL )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( bestOcc1 == NULL )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_SECOND_IN_PAIR;

    if ( bestOcc1 != NULL && bestOcc1->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    if ( bestOcc2 != NULL && bestOcc2->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
    samAlgnmt->core.tid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;    //SAM: chromosome ID, defined by bam_header_t
    samAlgnmt->core.pos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;     //SAM: 0-based leftmost coordinate
    samAlgnmt->core.mtid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;
    samAlgnmt->core.mpos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;

    if ( chr_1 > 0 && chr_1 == chr_2 )
        if ( tp_1 > tp_2 )
        { samAlgnmt->core.isize = tp_1 + readlen1 - tp_2; }
        else
        { samAlgnmt->core.isize = - ( tp_2 + readlen2 - tp_1 ); }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );
}


void unproperlypairDPOutputSAMAPI ( SRAQueryInput * qInput, Algnmt * algn_list1, Algnmt * algn_list2,
                                    int hitNum1, int hitNum2,
                                    unsigned char * query1, unsigned char * query2,
                                    char * qualities1, char * qualities2,
                                    int readlen1, int readlen2,
                                    char * queryName1, char * queryName2,
                                    DynamicUint8Array * xazArray )
{
    // cx this function... X_x
    //  fprintf(stderr,"unproperlypairDPOutputSAMAPI\n");

    // output the DP alignments which are not properly paired
    // since the alignments are not properly paired, it follows the single-end scoring scheme
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    OCC * occ = qSetting->occ;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long tp_1 = 0;
    unsigned long long tp_2 = 0;
    unsigned short chr_1 = 0;
    unsigned short chr_2 = 0;
    unsigned int samFlag;
    unsigned int i;
    unsigned int currAmbPos;
    unsigned long long currTP = 0;
    unsigned short currChr = 0;
    char currStrand;
    int currStrLen = 0;
    unsigned int currScore;
    char currOccStr[500];
    char * currSpCigar = NULL;
    int currEditDist = 0;
    char mdStr1[MAX_READ_LENGTH * 2 + 1];
    char mdStr2[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen1 = 0;
    int mdStrLen2 = 0;
    char cigarStr1[MAX_READ_LENGTH + 1];
    char cigarStr2[MAX_READ_LENGTH + 1];
    int cigarStrLen1 = 0;
    int cigarStrLen2 = 0;
    int map_qual_score1 = 0;
    int map_qual_score2 = 0;
    int bestMismatchNum1 = 0;
    int bestGapOpen1 = 0;
    int bestGapExtend1 = 0;
    int avg_mismatch_qual1;
    int bestMismatchNum2 = 0;
    int bestGapOpen2 = 0;
    int bestGapExtend2 = 0;
    int avg_mismatch_qual2;
    // obtain the best hit
    Algnmt * bestAlgn1 = NULL;
    int bestScore1;
    int bestScoreNum1 = 0;
    int boundTrim1 = 0, boundTrim2 = 0;
    char * newCigar1 = NULL, *newCigar2 = NULL;
    int deletedEnd1 = 0, deletedEnd2 = 0;

    if ( hitNum1 > 0 )
    {
        bestAlgn1 = & ( algn_list1[0] );
        bestScore1 = bestAlgn1->score;
        bestScoreNum1 = 1;

        for ( i = 1; i < hitNum1; i++ )
        {
            if ( algn_list1[i].score > bestScore1 )
            {
                bestAlgn1 = & ( algn_list1[i] );
                bestScore1 = bestAlgn1->score;
                bestScoreNum1 = 1;
            }
            else if ( algn_list1[i].score == bestScore1 )
            {
                bestScoreNum1++;
            }
        }
    }

    Algnmt * bestAlgn2 = NULL;
    int bestScore2;
    int bestScoreNum2 = 0;

    if ( hitNum2 > 0 )
    {
        bestAlgn2 = & ( algn_list2[0] );
        bestScore2 = bestAlgn2->score;
        bestScoreNum2 = 1;

        for ( i = 1; i < hitNum2; i++ )
        {
            if ( algn_list2[i].score > bestScore2 )
            {
                bestAlgn2 = & ( algn_list2[i] );
                bestScore2 = bestAlgn2->score;
                bestScoreNum2 = 1;
            }
            else if ( algn_list2[i].score == bestScore2 )
            {
                bestScoreNum2++;
            }
        }
    }

    // obtain x1_t1 and x1_t2, and secondBestScore
    int x1_t1_1 = 0;
    int x1_t2_1 = 0;
    int secondBestScore1 = -9999;
    int subopt_class_thres1 = ( int ) ( 0.7 * bestScore1 );

    for ( i = 0; i < hitNum1; i++ )
    {
        if ( algn_list1[i].score < bestScore1 )
        {
            if ( algn_list1[i].score > secondBestScore1 )
            {
                secondBestScore1 = algn_list1[i].score;
            }

            if ( algn_list1[i].score >= subopt_class_thres1 )
            {
                x1_t1_1++;
            }
            else
            {
                x1_t2_1++;
            }
        }
    }

    int secondBestNum1;

    if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
    { secondBestNum1 = -1; }
    else
    { secondBestNum1 = x1_t1_1 + x1_t2_1; }

    int x1_t1_2 = 0;
    int x1_t2_2 = 0;
    int secondBestScore2 = -9999;
    int subopt_class_thres2 = ( int ) ( 0.7 * bestScore2 );

    for ( i = 0; i < hitNum2; i++ )
    {
        if ( algn_list2[i].score < bestScore2 )
        {
            if ( algn_list2[i].score > secondBestScore2 )
            {
                secondBestScore2 = algn_list2[i].score;
            }

            if ( algn_list2[i].score >= subopt_class_thres2 )
            {
                x1_t1_2++;
            }
            else
            {
                x1_t2_2++;
            }
        }
    }

    int secondBestNum2;

    if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
    { secondBestNum2 = -1; }
    else
    { secondBestNum2 = x1_t1_2 + x1_t2_2; }

    if ( ( bestAlgn1 != NULL ) && ( hspaux->alignmentType != OUTPUT_UNIQUE_BEST || bestScoreNum1 == 1 ) )
    {
        //        getChrAndPos ( qInput, bestAlgn1->algnmt, &tp_1, &chr_1 );

        if ( bestAlgn1->isFromDP == 1 )
        {
            boundTrim1 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen1, bestAlgn1->algnmt, bestAlgn1->cigarString, &tp_1, &chr_1, &newCigar1 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen1 = boundTrim1 ? convertToCigarStr ( newCigar1, cigarStr1 ) : convertToCigarStr ( bestAlgn1->cigarString, cigarStr1, &deletedEnd1 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen1 = getMisInfoForDP ( hsp, query1, qualities1, readlen1, bestAlgn1->algnmt, bestAlgn1->strand,
                                          boundTrim1 ? newCigar1 : bestAlgn1->cigarString, mdStr1, &bestMismatchNum1, &bestGapOpen1, &bestGapExtend1,
                                          &avg_mismatch_qual1, boundTrim1 );
        }
        else
        {
            boundTrim1 = getChrAndPosWithBoundaryCheck ( qInput, readlen1, bestAlgn1->algnmt, &tp_1, &chr_1, &newCigar1 );
            // to get the md str
            mdStrLen1 = getMdStr ( hsp, query1, qualities1, readlen1, bestAlgn1->algnmt, bestAlgn1->strand, bestAlgn1->editdist, mdStr1, &avg_mismatch_qual1, boundTrim1 );
            bestMismatchNum1 = 0;

            for ( int ii = 0; ii < mdStrLen1; ++ii )
            { bestMismatchNum1 += ( mdStr1[ii] > '9' ); } // A C G or T
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
        {
            map_qual_score1 = SAM_MAPQ_UNAVAILABLE;
        }
        else
        {
            map_qual_score1 = getMapQualScoreForSingleDP ( readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, bestScoreNum1, x1_t1_1, x1_t2_1, bestScore1, secondBestScore1, hspaux->maxMAPQ, hspaux->minMAPQ, hspaux->singleDPcutoffThreshold, hspaux->bwaLikeScore, hspaux->g_log_n );

            if ( !hspaux->bwaLikeScore )
            { map_qual_score1 = map_qual_score1 >> 1; }

            if ( map_qual_score1 < hspaux->minMAPQ )
            {
                map_qual_score1 = hspaux->minMAPQ;
            }
        }

        if ( boundTrim1 ) { map_qual_score1 = 0; }
    }
    else
    {
        bestAlgn1 = NULL;
    }

    if ( ( bestAlgn2 != NULL ) && ( hspaux->alignmentType != OUTPUT_UNIQUE_BEST || bestScoreNum2 == 1 ) )
    {
        // getChrAndPos ( qInput, bestAlgn2->algnmt, &tp_2, &chr_2 );

        if ( bestAlgn2->isFromDP == 1 )
        {
            boundTrim2 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen2, bestAlgn2->algnmt, bestAlgn2->cigarString, &tp_2, &chr_2, &newCigar2 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen2 = boundTrim2 ? convertToCigarStr ( newCigar2, cigarStr2 ) : convertToCigarStr ( bestAlgn2->cigarString, cigarStr2, &deletedEnd2 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen2 = getMisInfoForDP ( hsp, query2, qualities2, readlen2, bestAlgn2->algnmt, bestAlgn2->strand,
                                          boundTrim2 ? newCigar2 : bestAlgn2->cigarString, mdStr2, &bestMismatchNum2, &bestGapOpen2, &bestGapExtend2,
                                          &avg_mismatch_qual2, boundTrim2 );

        }
        else
        {
            boundTrim2 = getChrAndPosWithBoundaryCheck ( qInput, readlen2, bestAlgn2->algnmt, &tp_2, &chr_2, &newCigar2 );
            // to get the md str
            mdStrLen2 = getMdStr ( hsp, query2, qualities2, readlen2, bestAlgn2->algnmt, bestAlgn2->strand, bestAlgn2->editdist, mdStr2, &avg_mismatch_qual2, boundTrim2 );
            bestMismatchNum2 = 0;

            for ( int ii = 0; ii < mdStrLen2; ++ii )
            { bestMismatchNum2 += ( mdStr2[ii] > '9' ); } // A C G or T
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST || hspaux->alignmentType == OUTPUT_UNIQUE_BEST )
        {
            map_qual_score2 = SAM_MAPQ_UNAVAILABLE;
        }
        else
        {
            map_qual_score2 = getMapQualScoreForSingleDP ( readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, bestScoreNum2, x1_t1_2, x1_t2_2, bestScore2, secondBestScore2, hspaux->maxMAPQ, hspaux->minMAPQ, hspaux->singleDPcutoffThreshold, hspaux->bwaLikeScore, hspaux->g_log_n );

            if ( !hspaux->bwaLikeScore )
            { map_qual_score2 = map_qual_score2 >> 1; }

            if ( map_qual_score2 < hspaux->minMAPQ )
            {
                map_qual_score2 = hspaux->minMAPQ;
            }

            if ( boundTrim2 ) { map_qual_score2 = 0; }
        }
    }
    else
    {
        bestAlgn2 = NULL;
    }

    //-------------------------------------------//
    // report the first alignment                //
    //-------------------------------------------//
    if ( hspaux->alignmentType == OUTPUT_ALL_BEST || hspaux->alignmentType == OUTPUT_ALL_VALID )
    {
        // Build the XA:Z tag
        DynamicUint8ArrayReset ( xazArray );

        for ( i = 0; i < hitNum1; i++ )
        {
            if ( bestAlgn1 == & ( algn_list1[i] ) )
            { continue; }

            if ( hspaux->alignmentType == OUTPUT_ALL_BEST && algn_list1[i].score < bestScore1 )
            { continue; }

            currScore = algn_list1[i].score;
            currStrand = algn_list1[i].strand;
            currSpCigar = algn_list1[i].cigarString;
            currEditDist = algn_list1[i].editdist;
            currAmbPos = algn_list1[i].algnmt;
            getChrAndPos ( qInput, algn_list1[i].algnmt,
                           &currTP, &currChr );
            char * chr_name = samFilePtr->header->target_name[currChr - 1]; // chromosome name
            memcpy ( currOccStr, chr_name, strlen ( chr_name ) );
            int pos = strlen ( chr_name );
            currOccStr[pos++] = ',';
            currOccStr[pos++] = ( currStrand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
            pos += writeULLToStr ( currTP, & ( currOccStr[pos] ) ); // chromosome position
            currOccStr[pos++] = ',';
            pos += convertToCigarStr ( currSpCigar, & ( currOccStr[pos] ) ); // cigar string
            currOccStr[pos++] = ',';
            pos += writeNumToStr ( currEditDist, & ( currOccStr[pos] ) ); // edit distance
            currOccStr[pos++] = ';';
            currOccStr[pos] = '\0';
            currStrLen = pos;
            appendStringToUint8Array ( xazArray, currOccStr, currStrLen );
        }
    }

    if ( bestAlgn1 != NULL )
        initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1,
                               qualities1, bestAlgn1->strand, xazArray->charStr, xazArray->length, bestAlgn1->isFromDP ? cigarStr1 : newCigar1, 0,
                               bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestScoreNum1, secondBestNum1, bestGapOpen1, bestGapExtend1, mdStr1, mdStrLen1, map_qual_score1, hspaux->readGroup,
                               hspaux->isPrintMDNM );
    else
        initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1,
                              qualities1, QUERY_POS_STRAND, xazArray->charStr,
                              xazArray->length, NULL, 1, hspaux->readGroup );

    // compute the value of the flag
    samFlag = 1;

    if ( bestAlgn1 == NULL )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( bestAlgn2 == NULL )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_FIRST_IN_PAIR;

    if ( bestAlgn1 != NULL && bestAlgn1->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    if ( bestAlgn2 != NULL && bestAlgn2->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
    samAlgnmt->core.tid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;    //SAM: chromosome ID, defined by bam_header_t
    samAlgnmt->core.pos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;     //SAM: 0-based leftmost coordinate
    samAlgnmt->core.mtid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;
    samAlgnmt->core.mpos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;

    if ( chr_1 > 0 && chr_1 == chr_2 )
        if ( tp_2 > tp_1 )
        { samAlgnmt->core.isize = tp_2 - deletedEnd2 + readlen2 - tp_1; }
        else
        { samAlgnmt->core.isize = - ( tp_1 - deletedEnd1 + readlen1 - tp_2 ); }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );

    //--------------------------------------------//
    // report the second alignment                //
    //--------------------------------------------//
    if ( hspaux->alignmentType == OUTPUT_ALL_BEST || hspaux->alignmentType == OUTPUT_ALL_VALID )
    {
        // Build the XA:Z tag
        DynamicUint8ArrayReset ( xazArray );

        for ( i = 0; i < hitNum2; i++ )
        {
            if ( bestAlgn2 == & ( algn_list2[i] ) )
            { continue; }

            if ( hspaux->alignmentType == OUTPUT_ALL_BEST && algn_list2[i].score < bestScore2 )
            { continue; }

            currScore = algn_list2[i].score;
            currStrand = algn_list2[i].strand;
            currSpCigar = algn_list2[i].cigarString;
            currEditDist = algn_list2[i].editdist;
            currAmbPos = algn_list2[i].algnmt;
            getChrAndPos ( qInput, algn_list2[i].algnmt, &currTP, &currChr );
            char * chr_name = samFilePtr->header->target_name[currChr - 1]; // chromosome name
            memcpy ( currOccStr, chr_name, strlen ( chr_name ) );
            int pos = strlen ( chr_name );
            currOccStr[pos++] = ',';
            currOccStr[pos++] = ( currStrand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
            pos += writeULLToStr ( currTP, & ( currOccStr[pos] ) ); // chromosome position
            currOccStr[pos++] = ',';
            pos += convertToCigarStr ( currSpCigar, & ( currOccStr[pos] ) ); // cigar string
            currOccStr[pos++] = ',';
            pos += writeNumToStr ( currEditDist, & ( currOccStr[pos] ) ); // edit distance
            currOccStr[pos++] = ';';
            currOccStr[pos] = '\0';
            currStrLen = pos;
            appendStringToUint8Array ( xazArray, currOccStr, currStrLen );
        }
    }

    if ( bestAlgn2 != NULL )
        initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2,
                               qualities2, bestAlgn2->strand, xazArray->charStr, xazArray->length, bestAlgn2->isFromDP ? cigarStr2 : newCigar2, 0,
                               bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestScoreNum2, secondBestNum2, bestGapOpen2, bestGapExtend2, mdStr2, mdStrLen2, map_qual_score2, hspaux->readGroup,
                               hspaux->isPrintMDNM );
    else
        initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2,
                              qualities2, QUERY_POS_STRAND, xazArray->charStr,
                              xazArray->length, NULL, 1, hspaux->readGroup );

    // compute the value of the flag
    samFlag = 1;

    if ( bestAlgn2 == NULL )
    {
        samFlag |= SAM_FLAG_READ_UNMAPPED;
    }

    if ( bestAlgn1 == NULL )
    {
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
    }

    samFlag |= SAM_FLAG_SECOND_IN_PAIR;

    if ( bestAlgn1 != NULL && bestAlgn1->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
    }

    if ( bestAlgn2 != NULL && bestAlgn2->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
    samAlgnmt->core.tid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;    //SAM: chromosome ID, defined by bam_header_t
    samAlgnmt->core.pos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;     //SAM: 0-based leftmost coordinate
    samAlgnmt->core.mtid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;
    samAlgnmt->core.mpos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;

    if ( chr_1 > 0 && chr_1 == chr_2 )
        if ( tp_1 > tp_2 )
        { samAlgnmt->core.isize = tp_1 - deletedEnd1 + readlen1 - tp_2; }
        else
        { samAlgnmt->core.isize = - ( tp_2 - deletedEnd2 + readlen2 - tp_1 ); }
    else
    { samAlgnmt->core.isize = 0; }

    samwrite ( samFilePtr, samAlgnmt );

    if ( newCigar1 ) { free ( newCigar1 ); }

    if ( newCigar2 ) { free ( newCigar2 ); }

}



OCCList * makeOccList ( unsigned int readID, HSPAux * hspaux, BWT * bwt )
{
    // fprintf(stderr, "[%u] enter the function makeOccList\n", readID);
    OCCList * occ_list = OCCListConstruct ();
    unsigned int sa_num = hspaux->sa_num[readID];
    unsigned int occ_num = hspaux->occ_num[readID];
    // fprintf(stderr, "[%u] sa_num = %u; occ_num = %u\n", readID, sa_num, occ_num);
    ReadInputForDP * ans = ( ( ReadInputForDP ** ) hspaux->soap3AnsArray ) [readID];
    unsigned int i, k, l, r;
    unsigned char strand, num_mis;

    if ( ans != NULL && sa_num > 0 )
    {
        // transfer all the SA ranges to occ_list
        unsigned int sa_start = hspaux->sa_start[readID];

        for ( i = sa_start; i < sa_start + sa_num; i++ )
        {
            if ( ans->sa_list[i].readID != readID )
            {
                fprintf ( stderr, "[sa] Error! readID does not match!" );
                exit ( 1 );
            }

            l = ans->sa_list[i].saIndexLeft;
            r = ans->sa_list[i].saIndexRight;
            strand = ans->sa_list[i].strand;
            num_mis = ans->sa_list[i].mismatchCount;

            for ( k = l; k <= r; k++ )
            {
                addToOCCList ( occ_list, BWTSaValue ( bwt, k ), strand, num_mis );
            }
        }
    }

    if ( ans != NULL && occ_num > 0 )
    {
        // insert all occ to occ_list
        unsigned int occ_start = hspaux->occ_start[readID];

        for ( i = occ_start; i < occ_start + occ_num; i++ )
        {
            if ( ans->occ_list[i].readID != readID )
            {
                fprintf ( stderr, "[occ] Error! readID does not match!" );
                exit ( 1 );
            }

            addToOCCList ( occ_list, ans->occ_list[i].ambPosition, ans->occ_list[i].strand, ans->occ_list[i].mismatchCount );
        }
    }

    // fprintf(stderr, "[%u] leave the function makeOccList\n", readID);
    return occ_list;
}

void unproperlypairOutputSAMAPI2 ( SRAQueryInput * qInput, unsigned int readID,
                                   unsigned char * query1, unsigned char * query2,
                                   char * qualities1, char * qualities2,
                                   int readlen1, int readlen2,
                                   char * queryName1, char * queryName2,
                                   DynamicUint8Array * xazArray )
{
    //    fprintf(stderr,"unproperlypairOutputSAMAPI2\n");
    // output the alignments which are not properly paired
    // the readID refers the first read of the pair
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    HSPAux * hspaux = aIndex->hspaux;
    BWT * bwt = aIndex->bwt;
    // make the arrays of occ_list1 and occ_list2
    OCCList * occ_list1 = makeOccList ( readID, hspaux, bwt );
    OCCList * occ_list2 = makeOccList ( readID + 1, hspaux, bwt );
    unproperlypairOutputSAMAPI ( qInput, occ_list1, occ_list2,
                                 query1, query2,
                                 qualities1, qualities2,
                                 readlen1, readlen2,
                                 queryName1, queryName2,
                                 xazArray,
                                 hspaux->peMaxOutputPerRead, hspaux->alignmentType );
    // free the arrays
    OCCListFree ( occ_list1 );
    OCCListFree ( occ_list2 );
}


void pairOutputSAMAPI ( SRAQueryInput * qInput, PEOutput * pe_out, PEPairs * bestPair,
                        unsigned char * query1, unsigned char * query2,
                        char * qualities1, char * qualities2,
                        int readlen1, int readlen2,
                        char * queryName1, char * queryName2,
                        char min_totalMismatchCount, char secMin_totalMismatchCount,
                        int outputXAZTag, DynamicUint8Array * xazArray,
                        unsigned int peMaxOutputPerPair,
                        int X0_first, int X0_second, int X1_first, int X1_second,
                        int num_minMismatchPair, char isBestHit1, char isBestHit2,
                        unsigned int totalNumValidPairs )
{
    // cx this function... X_x
    //    fprintf(stderr,"pairOutputSAMAPI\n");

    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRASetting * qSetting = qInput->QuerySetting;
    OCC * occ = qSetting->occ;
    // BWT * bwt = aIndex->bwt;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long tp_1, tp_2, curr_tp;
    unsigned short chr_1, chr_2, curr_chr;
    unsigned int samFlag;
    unsigned int i;
    char curr_strand;
    int curr_len = 0;
    unsigned int curr_numMisMatch;
    char curr_occStr[500];
    unsigned int total_results;
    char mdStr1[MAX_READ_LENGTH * 2 + 1];
    char mdStr2[MAX_READ_LENGTH * 2 + 1];
    int avg_mismatch_qual1;
    int avg_mismatch_qual2;
    //int map_qual_score = 0;
    int mapq1, mapq2;
    int boundTrim1 = 0, boundTrim2 = 0;
    char * newCigar1 = NULL, *newCigar2 = NULL;
    int bestMismatch1, bestMismatch2;

    if ( bestPair != NULL )
    {

        boundTrim1 = getChrAndPosWithBoundaryCheck ( qInput, readlen1, bestPair->algnmt_1, &tp_1, &chr_1, &newCigar1 );
        boundTrim2 = getChrAndPosWithBoundaryCheck ( qInput, readlen2, bestPair->algnmt_2, &tp_2, &chr_2, &newCigar2 );

        /*      if ((bestPair->strand_1 == QUERY_POS_STRAND &&
             (bestPair->algnmt_1 > bestPair->algnmt_2 || bestPair->algnmt_1 + readlen1 > bestPair->algnmt_2 + readlen2)) ||
            (bestPair->strand_1 == QUERY_NEG_STRAND &&
             (bestPair->algnmt_2 > bestPair->algnmt_1 || bestPair->algnmt_2 + readlen2 > bestPair->algnmt_1 + readlen1))) {
          // TODO read-through?
          }*/

        // get the MD string
        int mdLen1 = getMdStr ( hsp, query1, qualities1, readlen1, bestPair->algnmt_1, bestPair->strand_1, bestPair->mismatch_1, mdStr1, &avg_mismatch_qual1, boundTrim1 );
        int mdLen2 = getMdStr ( hsp, query2, qualities2, readlen2, bestPair->algnmt_2, bestPair->strand_2, bestPair->mismatch_2, mdStr2, &avg_mismatch_qual2, boundTrim2 );

        bestMismatch1 = bestPair->mismatch_1;
        bestMismatch2 = bestPair->mismatch_2;

        if ( boundTrim1 && bestMismatch1 )
        {
            bestMismatch1 = 0;

            for ( int ii = 0; ii < mdLen1; ++ii )
            { bestMismatch1 += ( mdStr1[ii] > '9' ); } // A C G or T
        }

        if ( boundTrim2 && bestMismatch2 )
        {
            bestMismatch2 = 0;

            for ( int ii = 0; ii < mdLen2; ++ii )
            { bestMismatch2 += ( mdStr2[ii] > '9' ); } // A C G or T
        }

        if ( hspaux->alignmentType == OUTPUT_ALL_VALID || hspaux->alignmentType == OUTPUT_ALL_BEST )
        {
            if ( hspaux->bwaLikeScore )
            {
                int op_score = ( readlen1 + readlen2 - min_totalMismatchCount ) * hspaux->dpMatchScore + min_totalMismatchCount * hspaux->dpMisMatchScore;
                int subop_score = ( readlen1 + readlen2 - secMin_totalMismatchCount ) * hspaux->dpMatchScore + secMin_totalMismatchCount * hspaux->dpMisMatchScore;
                int subop_num = ( totalNumValidPairs - num_minMismatchPair ) > 0 ? totalNumValidPairs - num_minMismatchPair : 0;
                bwaLikePairQualScore ( X0_first, X1_first, X0_second, X1_second, hspaux->g_log_n, op_score, num_minMismatchPair, subop_score, subop_num, readlen1, readlen2, &mapq1, &mapq2 );
            }
            else
            {
                int map_qual_score1 = getMapQualScore2 ( bestPair->mismatch_1, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, X0_first, X1_first, isBestHit1, totalNumValidPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                int map_qual_score2 = getMapQualScore2 ( bestPair->mismatch_2, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, X0_second, X1_second, isBestHit2, totalNumValidPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                mapq1 = mapq2 = getMapQualScoreForPair ( map_qual_score1, map_qual_score2 );
            }

            if ( boundTrim1 ) { mapq1 = 0; }

            if ( boundTrim2 ) { mapq2 = 0; }
        }
        else
        {
            mapq1 = mapq2 = SAM_MAPQ_UNAVAILABLE;
        }


        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( outputXAZTag == 1 )
        {
            // Build the XA:Z tag
            PEPairList * pairList = pe_out->root;
            total_results = 1; // including the best result

            // TODO cx suboptimal alignments?

            while ( pairList != NULL && pairList->pairsCount > 0 && total_results < peMaxOutputPerPair )
            {
                for ( i = 0; i < pairList->pairsCount && total_results < peMaxOutputPerPair; i++ )
                {
                    PEPairs * pePair = & ( pairList->pairs[i] );

                    if ( pePair == bestPair || pePair->totalMismatchCount > min_totalMismatchCount )
                    { continue; }

                    curr_numMisMatch = pePair->mismatch_1;
                    curr_strand = pePair->strand_1;
                    getChrAndPos ( qInput, pePair->algnmt_1, &curr_tp, &curr_chr );
                    char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                    memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                    int pos = strlen ( chr_name );
                    curr_occStr[pos++] = ',';
                    curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                    pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( readlen1, & ( curr_occStr[pos] ) ); // cigar string
                    curr_occStr[pos++] = 'M';
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( curr_numMisMatch, & ( curr_occStr[pos] ) ); // number of mismatches
                    curr_occStr[pos++] = ';';
                    curr_occStr[pos] = '\0';
                    curr_len = pos;

                    appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
                    total_results++;
                }

                pairList = pairList->next;
            }
        }

        initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1,
                               qualities1, bestPair->strand_1, xazArray->charStr, xazArray->length,
                               newCigar1, 0, bestMismatch1, bestMismatch1,
                               X0_first, X1_first, 0, 0, mdStr1, mdLen1, mapq1,
                               hspaux->readGroup, hspaux->isPrintMDNM );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_PROPER_PAIR;
        samFlag |= SAM_FLAG_FIRST_IN_PAIR;

        if ( bestPair->strand_1 == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( bestPair->strand_2 == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = chr_1 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = chr_2 - 1;
        samAlgnmt->core.mpos = tp_2 - 1;

        if ( chr_1 == chr_2 )
            if ( tp_2 > tp_1 )
            { samAlgnmt->core.isize = tp_2 + readlen2 - tp_1; }
            else
            { samAlgnmt->core.isize = - ( tp_1 + readlen1 - tp_2 ); }
        else
        { samAlgnmt->core.isize = 0; }

        samwrite ( samFilePtr, samAlgnmt );
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( outputXAZTag == 1 )
        {
            // Build the XA:Z tag
            PEPairList * pairList = pe_out->root;
            total_results = 1; // including the best result

            // TODO cx suboptimal alignments?

            while ( pairList != NULL && pairList->pairsCount > 0 && total_results < peMaxOutputPerPair )
            {
                for ( i = 0; i < pairList->pairsCount && total_results < peMaxOutputPerPair; i++ )
                {
                    PEPairs * pePair = & ( pairList->pairs[i] );

                    if ( pePair == bestPair || pePair->totalMismatchCount > min_totalMismatchCount )
                    { continue; }

                    curr_numMisMatch = pePair->mismatch_2;
                    curr_strand = pePair->strand_2;
                    getChrAndPos ( qInput, pePair->algnmt_2, &curr_tp, &curr_chr );
                    char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                    memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                    int pos = strlen ( chr_name );
                    curr_occStr[pos++] = ',';
                    curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                    pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( readlen2, & ( curr_occStr[pos] ) ); // cigar string
                    curr_occStr[pos++] = 'M';
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( curr_numMisMatch, & ( curr_occStr[pos] ) ); // number of mismatches
                    curr_occStr[pos++] = ';';
                    curr_occStr[pos] = '\0';
                    curr_len = pos;

                    appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
                    total_results++;
                }

                pairList = pairList->next;
            }
        }

        initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2,
                               qualities2, bestPair->strand_2, xazArray->charStr, xazArray->length,
                               newCigar2, 0, bestMismatch2, bestMismatch2, X0_second, X1_second,
                               0, 0, mdStr2, mdLen2, mapq2, hspaux->readGroup, hspaux->isPrintMDNM );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_PROPER_PAIR;
        samFlag |= SAM_FLAG_SECOND_IN_PAIR;

        if ( bestPair->strand_1 == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        if ( bestPair->strand_2 == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = chr_2 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = tp_2 - 1;                       //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = chr_1 - 1;
        samAlgnmt->core.mpos = tp_1 - 1;

        if ( chr_1 == chr_2 )
        {
            if ( tp_1 > tp_2 )
            { samAlgnmt->core.isize = tp_1 + readlen1 - tp_2; }
            else
            { samAlgnmt->core.isize = - ( tp_2 + readlen2 - tp_1 ); }
        }
        else
        {
            samAlgnmt->core.isize = 0;
        }

        samwrite ( samFilePtr, samAlgnmt );
    }
    else
    {
        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1,
                              qualities1, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_FIRST_IN_PAIR;
        samFlag |= SAM_FLAG_READ_UNMAPPED;
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2,
                              qualities2, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_SECOND_IN_PAIR;
        samFlag |= SAM_FLAG_READ_UNMAPPED;
        samFlag |= SAM_FLAG_MATE_UNMAPPED;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
    }

    if ( newCigar1 ) { free ( newCigar1 ); }

    if ( newCigar2 ) { free ( newCigar2 ); }
}

int readLengthWithCigar ( char * cigar )
{
    char * p = cigar;
    char op = ' ';
    int x = 0;
    int len = 0;

    while ( *p )
    {
        x = 0;

        while ( *p <= '9' )
        {
            x = x * 10 + ( *p - '0' );
            p++;
        }

        op = *p;
        p++;

        if ( op == 'M' || op == 'm' || op == 'D' ) { len += x; }
    }

    if ( op == 'D' ) // last delete
    { len -= x; } // ignore => undo

    return len;
}

void pairDeepDPOutputSAMAPI ( SRAQueryInput * qInput, DeepDPAlignResult * algnResult,
                              DeepDPAlignResult * bestResult,
                              unsigned int start, unsigned int num,
                              unsigned char * query1, unsigned char * query2,
                              char * qualities1, char * qualities2,
                              int readlen1, int readlen2,
                              char * queryName1, char * queryName2,
                              DynamicUint8Array * xazArray )
{
    // cx this function... X_x
    //    fprintf(stderr,"pairDeepDPOutputSAMAPI\n");

    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    OCC * occ = qSetting->occ;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long tp_1 = 0;
    unsigned long long tp_2 = 0;
    unsigned long long curr_tp = 0;
    unsigned short chr_1 = 0;
    unsigned short chr_2 = 0;
    unsigned short curr_chr = 0;
    unsigned int samFlag;
    unsigned int i;
    char curr_strand;
    int curr_len = 0;
    char curr_occStr[500];
    char best_strand1 = 1;
    char best_strand2 = 1;
    char * best_cigar1 = NULL;
    char * best_cigar2 = NULL;
    int best_insert = 0;
    // int mapping_qual_score = 0; // for a valid paired-end alignment
    int mapq1, mapq2;
    char mdStr1[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen1 = 0;
    char cigarStr1[MAX_READ_LENGTH + 1];
    int cigarStrLen1 = 0;
    int bestMismatchNum1 = 0;
    int bestGapOpen1 = 0;
    int bestGapExtend1 = 0;
    int avg_mismatch_qual1;
    char mdStr2[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen2 = 0;
    char cigarStr2[MAX_READ_LENGTH + 1];
    int cigarStrLen2 = 0;
    int bestMismatchNum2 = 0;
    int bestGapOpen2 = 0;
    int bestGapExtend2 = 0;
    int avg_mismatch_qual2;

    int boundTrim1 = 0, boundTrim2 = 0;
    char * newCigar1 = NULL, *newCigar2 = NULL;
    int r1 = readlen1, r2 = readlen2;

    if ( bestResult != NULL )
    {
        // for the reads which do not have any hit for both ends
        if ( bestResult->algnmt_1 == 0xFFFFFFFF && bestResult->algnmt_2 == 0xFFFFFFFF )
        {
            unproperlypairOutputSAMAPI2 ( qInput, bestResult->readID,
                                          query1, query2,
                                          qualities1, qualities2,
                                          readlen1, readlen2,
                                          queryName1, queryName2,
                                          xazArray );
            return;
        }

        // collect the information of the first alignment
        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            //            getChrAndPos ( qInput, bestResult->algnmt_1,
            //                           &tp_1, &chr_1 );
            best_strand1 = bestResult->strand_1;
            best_cigar1 = bestResult->cigarString_1;
        }

        // collect the information of the second alignment
        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            //            getChrAndPos ( qInput, bestResult->algnmt_2,
            //                           &tp_2, &chr_2 );
            best_strand2 = bestResult->strand_2;
            best_cigar2 = bestResult->cigarString_2;
        }

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            best_insert = bestResult->insertSize;
        }

        // obtain the corresponding normal cigar string, # of gap open, # of gap extension, # of mismatch,
        // and the average quality values in the mismatch positions
        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            boundTrim1 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen1, bestResult->algnmt_1, bestResult->cigarString_1, &tp_1, &chr_1, &newCigar1 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen1 = boundTrim1 ? convertToCigarStr ( newCigar1, cigarStr1 ) : convertToCigarStr ( bestResult->cigarString_1, cigarStr1 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen1 = getMisInfoForDP ( hsp, query1, qualities1, readlen1, bestResult->algnmt_1, best_strand1,
                                          boundTrim1 ? newCigar1 : bestResult->cigarString_1, mdStr1, &bestMismatchNum1, &bestGapOpen1, &bestGapExtend1, &avg_mismatch_qual1, boundTrim1 );
        }

        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            boundTrim2 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen2, bestResult->algnmt_2, bestResult->cigarString_2, &tp_2, &chr_2, &newCigar2 );
            // to convert the the best special_cigar into normal cigar string
            cigarStrLen2 = boundTrim2 ? convertToCigarStr ( newCigar2, cigarStr2 ) : convertToCigarStr ( bestResult->cigarString_2, cigarStr2 );
            // to collect the mismatch information including MD string and # of mismatches
            mdStrLen2 = getMisInfoForDP ( hsp, query2, qualities2, readlen2, bestResult->algnmt_2, best_strand2,
                                          boundTrim2 ? newCigar2 : bestResult->cigarString_2, mdStr2, &bestMismatchNum2, &bestGapOpen2, &bestGapExtend2, &avg_mismatch_qual2, boundTrim2 );
        }

        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            r1 = readLengthWithCigar ( bestResult->cigarString_1 );
        }

        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            r2 = readLengthWithCigar ( bestResult->cigarString_2 );
        }

        // check for read-through reads
        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            unsigned long long adjust1 = bestResult->algnmt_1 + ( boundTrim1 > 0 ? boundTrim1 : 0 );
            unsigned long long adjust2 = bestResult->algnmt_2 + ( boundTrim2 > 0 ? boundTrim2 : 0 );

            if ( ( bestResult->strand_1 == QUERY_POS_STRAND &&
                    ( adjust1 > adjust2 || adjust1 + r1 > adjust2 + r2 ) ) ||
                    ( bestResult->strand_1 == QUERY_NEG_STRAND &&
                      ( adjust2 > adjust1 || adjust2 + r2 > adjust1 + r1 ) ) )
            {
                if ( bestMismatchNum1 <= bestMismatchNum2 )
                {
                    bestResult->algnmt_2 = 0xFFFFFFFF;
                    tp_2 = chr_2 = 0;
                }
                else
                {
                    bestResult->algnmt_1 = 0xFFFFFFFF;
                    tp_1 = chr_1 = 0;
                }

                best_insert = 0;
            }
        }

        int bestPairNum = 0;
        int bestPairScore = 0;
        int secBestPairScore = 0;

        // obtain the number of pairs with maximum sum of scores
        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            bestPairNum = 1;
            bestPairScore = bestResult->score_1 + bestResult->score_2;

            if ( num > 1 )
            {
                for ( i = start; i < start + num; i++ )
                {
                    if ( & ( algnResult[i] ) == bestResult )
                    { continue; }

                    if ( algnResult[i].score_1 + algnResult[i].score_2 == bestPairScore )
                    { bestPairNum++; }
                    else if ( algnResult[i].score_1 + algnResult[i].score_2 > secBestPairScore )
                    { secBestPairScore = algnResult[i].score_1 + algnResult[i].score_2; }
                }
            }
        }

        // obtain the number of pairs with similar score to the best score (i.e. difference <= 1 mismatch score)
        int numSimilarBestPairs = 0;

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            numSimilarBestPairs = 1; // includes the best result
            int bestPairScore1 = bestResult->score_1;
            int bestPairScore2 = bestResult->score_2;

            if ( num > 1 )
            {
                for ( i = start; i < start + num; i++ )
                {
                    if ( & ( algnResult[i] ) == bestResult )
                    { continue; }

                    if ( ( algnResult[i].score_1 >= ( bestPairScore1 + hspaux->dpMisMatchScore ) ) &&
                            ( algnResult[i].score_2 >= ( bestPairScore2 + hspaux->dpMisMatchScore ) ) )
                    {
                        numSimilarBestPairs++;
                    }
                }
            }
        }

        // obtain the bestHitNum, bestScore, secBestHitNum and secBestScore for the first read
        int bestHitNum1 = 0;
        int secBestHitNum1 = 0;
        int bestScore1 = 0;
        int secBestScore1 = 0;
        char isBestHit1 = 1;
        unsigned int bestPos1 = 0xFFFFFFFF;
        unsigned int secBestPos1 = 0xFFFFFFFF;

        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            bestScore1 = bestResult->score_1;
            bestPos1 = bestResult->algnmt_1; // alignment position
            bestHitNum1 = 1;
            // bestHitNum1 = bestResult->num_sameScore_1;
        }

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && num > 1 )
        {
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].score_1 >= bestScore1 )
                {
                    if ( algnResult[i].score_1 == bestScore1 )
                    {
                        if ( algnResult[i].algnmt_1 != bestPos1 )
                            // bestHitNum1++;
                        { bestHitNum1 += algnResult[i].num_sameScore_1; }
                    }
                    else
                    {
                        secBestScore1 = bestScore1;
                        secBestHitNum1 = bestHitNum1;
                        secBestPos1 = bestPos1;
                        bestScore1 = algnResult[i].score_1;
                        // bestHitNum1 = 1;
                        bestHitNum1 = algnResult[i].num_sameScore_1;
                        bestPos1 = algnResult[i].algnmt_1;
                        isBestHit1 = 0;
                    }
                }
                else if ( algnResult[i].score_1 >= secBestScore1 )
                {
                    if ( algnResult[i].score_1 == secBestScore1 )
                    {
                        if ( algnResult[i].algnmt_1 != secBestPos1 )
                            // secBestHitNum1++;
                        { secBestHitNum1 += algnResult[i].num_sameScore_1; }
                    }
                    else
                    {
                        secBestScore1 = algnResult[i].score_1;
                        secBestPos1 = algnResult[i].algnmt_1;
                        // secBestHitNum1 = 1;
                        secBestHitNum1 = algnResult[i].num_sameScore_1;
                    }
                }
            }
        }

        if ( hspaux->x0_array[bestResult->readID] > 0 )
        {
            int x0_score = hspaux->mismatch_array[bestResult->readID] * hspaux->dpMisMatchScore
                           + ( readlen1 - hspaux->mismatch_array[bestResult->readID] ) * hspaux->dpMatchScore;

            if ( x0_score >= bestScore1 )
            {
                bestHitNum1 = ( bestHitNum1 > hspaux->x0_array[bestResult->readID] ) ? bestHitNum1 : hspaux->x0_array[bestResult->readID];
                secBestHitNum1 = ( secBestHitNum1 > hspaux->x1_array[bestResult->readID] ) ? secBestHitNum1 : hspaux->x1_array[bestResult->readID];

                if ( x0_score > bestScore1 )
                { isBestHit1 = 0; }
            }
        }

        // obtain the bestHitNum, bestScore, secBestHitNum and secBestScore for the second read
        int bestHitNum2 = 0;
        int secBestHitNum2 = 0;
        int bestScore2 = 0;
        int secBestScore2 = 0;
        char isBestHit2 = 1;
        unsigned int bestPos2 = 0xFFFFFFFF;
        unsigned int secBestPos2 = 0xFFFFFFFF;

        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            bestScore2 = bestResult->score_2;
            bestPos2 = bestResult->algnmt_2; // alignment position
            bestHitNum2 = 1;
            // bestHitNum2 = bestResult->num_sameScore_2;
        }

        if ( bestResult->algnmt_2 != 0xFFFFFFFF && num > 1 )
        {
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].score_2 >= bestScore2 )
                {
                    if ( algnResult[i].score_2 == bestScore2 )
                    {
                        if ( algnResult[i].algnmt_2 != bestPos2 )
                            // bestHitNum2++;
                        { bestHitNum2 += algnResult[i].num_sameScore_2; }
                    }
                    else
                    {
                        secBestScore2 = bestScore2;
                        secBestHitNum2 = bestHitNum2;
                        secBestPos2 = bestPos2;
                        bestScore2 = algnResult[i].score_2;
                        // bestHitNum2 = 1;
                        bestHitNum2 = algnResult[i].num_sameScore_2;
                        bestPos2 = algnResult[i].algnmt_2;
                        isBestHit2 = 0;
                    }
                }
                else if ( algnResult[i].score_2 >= secBestScore2 )
                {
                    if ( algnResult[i].score_2 == secBestScore2 )
                    {
                        if ( algnResult[i].algnmt_2 != secBestPos2 )
                            // secBestHitNum2++;
                        { secBestHitNum2 += algnResult[i].num_sameScore_2; }
                    }
                    else
                    {
                        secBestScore2 = algnResult[i].score_2;
                        secBestPos2 = algnResult[i].algnmt_2;
                        // secBestHitNum2 = 1;
                        secBestHitNum2 = algnResult[i].num_sameScore_2;
                    }
                }
            }
        }

        if ( hspaux->x0_array[bestResult->readID + 1] > 1 )
        {
            int x0_score = hspaux->mismatch_array[bestResult->readID + 1] * hspaux->dpMisMatchScore
                           + ( readlen2 - hspaux->mismatch_array[bestResult->readID + 1] ) * hspaux->dpMatchScore;

            if ( x0_score >= bestScore2 )
            {
                bestHitNum2 = ( bestHitNum2 > hspaux->x0_array[bestResult->readID + 1] ) ? bestHitNum2 : hspaux->x0_array[bestResult->readID + 1];
                secBestHitNum2 = ( secBestHitNum2 > hspaux->x1_array[bestResult->readID + 1] ) ? secBestHitNum2 : hspaux->x1_array[bestResult->readID + 1];

                if ( x0_score > bestScore2 )
                { isBestHit2 = 0; }
            }
        }

        if ( hspaux->alignmentType == OUTPUT_ALL_VALID || hspaux->alignmentType == OUTPUT_ALL_BEST )
        {
            if ( hspaux->bwaLikeScore )
            {
                bwaLikePairQualScore ( bestHitNum1, secBestHitNum1, bestHitNum2, secBestHitNum2, hspaux->g_log_n, bestPairScore, bestPairNum, secBestPairScore, num - bestPairNum, readlen1, readlen2, &mapq1, &mapq2 );
            }
            else
            {
                int mapping_qual_score1 = getMapQualScoreForDP2 ( bestResult->score_1, readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, bestHitNum1, secBestHitNum1, bestScore1, secBestScore1, isBestHit1, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                int mapping_qual_score2 = getMapQualScoreForDP2 ( bestResult->score_2, readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, bestHitNum2, secBestHitNum2, bestScore2, secBestScore2, isBestHit2, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                mapq1 = mapq2 = getMapQualScoreForPair ( mapping_qual_score1, mapping_qual_score2 );
            }

            if ( boundTrim1 ) { mapq1 = 0; }

            if ( boundTrim2 ) { mapq2 = 0; }

        }
        else
        {
            mapq1 = mapq2 = SAM_MAPQ_UNAVAILABLE;
        }

        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && num > 1 )
        {
            // Build the XA:Z tag
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) &&
                        ( algnResult[i].score_1 + algnResult[i].score_2 < bestPairScore ) )
                {
                    continue;
                }

                getChrAndPos ( qInput, algnResult[i].algnmt_1, &curr_tp, &curr_chr );
                curr_strand = algnResult[i].strand_1;
                char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                int pos = strlen ( chr_name );
                curr_occStr[pos++] = ',';
                curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                curr_occStr[pos++] = ',';
                pos += convertToCigarStr ( algnResult[i].cigarString_1, & ( curr_occStr[pos] ) ); // cigar string
                curr_occStr[pos++] = ',';
                pos += writeNumToStr ( algnResult[i].editdist_1, & ( curr_occStr[pos] ) ); // edit distance
                curr_occStr[pos++] = ';';
                curr_occStr[pos] = '\0';
                curr_len = pos;
                appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            }
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST )
        {
            bestHitNum1 = -1;
            secBestHitNum1 = -1;
        }
        else if ( hspaux->alignmentType != OUTPUT_ALL_VALID && hspaux->alignmentType != OUTPUT_ALL_BEST )
        {
            secBestHitNum1 = -1;
        }

        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            if ( bestResult->algnmt_2 != 0xFFFFFFFF )
            {
                initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1,
                                       qualities1, best_strand1, xazArray->charStr, xazArray->length, cigarStr1, 0,
                                       bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestHitNum1, secBestHitNum1, bestGapOpen1,
                                       bestGapExtend1, mdStr1, mdStrLen1, mapq1, hspaux->readGroup, hspaux->isPrintMDNM );
            }
            else
            {
                mapq1 = getMapQualScoreForDP ( bestHitNum1, bestResult->score_1, readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, hspaux->maxMAPQ, hspaux->minMAPQ );

                if ( boundTrim1 ) { mapq1 = 0; }

                initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1,
                                       qualities1, best_strand1, xazArray->charStr, xazArray->length, cigarStr1, 0,
                                       bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestHitNum1, secBestHitNum1, bestGapOpen1,
                                       bestGapExtend1, mdStr1, mdStrLen1, mapq1, hspaux->readGroup, hspaux->isPrintMDNM );
            }
        }
        else
        {
            initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1,
                                  qualities1, best_strand1, NULL, 0, NULL, 1, hspaux->readGroup );
        }

        // compute the value of the flag
        samFlag = 1;

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->algnmt_2 != 0xFFFFFFFF ) )
        { samFlag |= SAM_FLAG_PROPER_PAIR; }

        samFlag |= SAM_FLAG_FIRST_IN_PAIR;

        if ( bestResult->algnmt_1 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestResult->algnmt_2 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->strand_1 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( ( bestResult->algnmt_2 != 0xFFFFFFFF ) && ( bestResult->strand_2 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;
        samAlgnmt->core.mpos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;

        //      samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        //      samAlgnmt->core.tid = chr_1 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        //      samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
        //      samAlgnmt->core.mtid = chr_2 - 1;
        //      samAlgnmt->core.mpos = tp_2 - 1;

        if ( best_insert > 0 )
        {
            if ( tp_1 > tp_2 )
            { samAlgnmt->core.isize = - ( tp_1 + r1 - tp_2 ); }
            else
            { samAlgnmt->core.isize = tp_2 + r2 - tp_1; }
        }
        else
        {
            samAlgnmt->core.isize = 0;
        }

        samwrite ( samFilePtr, samAlgnmt );
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( bestResult->algnmt_2 != 0xFFFFFFFF && num > 1 )
        {
            // Build the XA:Z tag
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) &&
                        ( algnResult[i].score_1 + algnResult[i].score_2 < bestPairScore ) )
                {
                    continue;
                }

                getChrAndPos ( qInput, algnResult[i].algnmt_2, &curr_tp, &curr_chr );
                curr_strand = algnResult[i].strand_2;
                char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                int pos = strlen ( chr_name );
                curr_occStr[pos++] = ',';
                curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                curr_occStr[pos++] = ',';
                pos += convertToCigarStr ( algnResult[i].cigarString_2, & ( curr_occStr[pos] ) ); // cigar string
                curr_occStr[pos++] = ',';
                pos += writeNumToStr ( algnResult[i].editdist_2, & ( curr_occStr[pos] ) ); // edit distance
                curr_occStr[pos++] = ';';
                curr_occStr[pos] = '\0';
                curr_len = pos;
                appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            }
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST )
        {
            bestHitNum2 = -1;
            secBestHitNum2 = -1;
        }
        else if ( hspaux->alignmentType != OUTPUT_ALL_VALID && hspaux->alignmentType != OUTPUT_ALL_BEST )
        {
            secBestHitNum2 = -1;
        }

        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            if ( bestResult->algnmt_1 != 0xFFFFFFFF )
            {
                initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2,
                                       qualities2, best_strand2, xazArray->charStr, xazArray->length, cigarStr2, 0,
                                       bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestHitNum2, secBestHitNum2, bestGapOpen2,
                                       bestGapExtend2, mdStr2, mdStrLen2, mapq2, hspaux->readGroup, hspaux->isPrintMDNM );
            }
            else
            {
                mapq2 = getMapQualScoreForDP ( bestHitNum2, bestResult->score_2, readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, hspaux->maxMAPQ, hspaux->minMAPQ );

                if ( boundTrim2 ) { mapq2 = 0; }

                initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2,
                                       qualities2, best_strand2, xazArray->charStr, xazArray->length, cigarStr2, 0,
                                       bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestHitNum2, secBestHitNum2, bestGapOpen2,
                                       bestGapExtend2, mdStr2, mdStrLen2, mapq2, hspaux->readGroup, hspaux->isPrintMDNM );
            }
        }
        else
        {
            initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2,
                                  qualities2, best_strand2, NULL, 0, NULL, 1, hspaux->readGroup );
        }

        // compute the value of the flag
        samFlag = 1;

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->algnmt_2 != 0xFFFFFFFF ) )
        { samFlag |= SAM_FLAG_PROPER_PAIR; }

        samFlag |= SAM_FLAG_SECOND_IN_PAIR;

        if ( bestResult->algnmt_2 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestResult->algnmt_1 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        if ( ( bestResult->algnmt_2 != 0xFFFFFFFF ) && ( bestResult->strand_2 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->strand_1 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;
        samAlgnmt->core.mpos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;

        //      samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        //      samAlgnmt->core.tid = chr_2 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        //      samAlgnmt->core.pos = tp_2 - 1;                       //SAM: 0-based leftmost coordinate
        //      samAlgnmt->core.mtid = chr_1 - 1;
        //      samAlgnmt->core.mpos = tp_1 - 1;

        if ( best_insert > 0 )
        {
            if ( tp_2 > tp_1 )
            { samAlgnmt->core.isize = - ( tp_2 + r2 - tp_1 ); }
            else
            { samAlgnmt->core.isize = tp_1 + r1 - tp_2; }
        }
        else
        {
            samAlgnmt->core.isize = 0;
        }

        samwrite ( samFilePtr, samAlgnmt );
    }
    else
    {
        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1,
                              qualities1, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_FIRST_IN_PAIR;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2,
                              qualities2, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_SECOND_IN_PAIR;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
    }

    if ( newCigar1 ) { free ( newCigar1 ); }

    if ( newCigar2 ) { free ( newCigar2 ); }

}


void pairDPOutputSAMAPI ( SRAQueryInput * qInput, AlgnmtDPResult * algnResult,
                          AlgnmtDPResult * bestResult,
                          unsigned int start, unsigned int num,
                          unsigned char * query1, unsigned char * query2,
                          char * qualities1, char * qualities2,
                          int readlen1, int readlen2,
                          char * queryName1, char * queryName2,
                          DynamicUint8Array * xazArray )
{
    // cx this function... X_x
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    OCC * occ = qSetting->occ;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long tp_1 = 0;
    unsigned long long tp_2 = 0;
    unsigned long long curr_tp = 0;
    unsigned short chr_1 = 0;
    unsigned short chr_2 = 0;
    unsigned short curr_chr = 0;
    unsigned int samFlag;
    unsigned int i;
    char curr_strand;
    int curr_len = 0;
    char curr_occStr[500];
    char best_strand1 = 1;
    char best_strand2 = 1;
    char * best_cigar1 = NULL;
    char * best_cigar2 = NULL;
    int best_insert = 0;
    // int mapping_qual_score=0; // for a valid paired-end alignment
    int mapq1, mapq2;
    char mdStr1[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen1 = 0;
    char cigarStr1[MAX_READ_LENGTH + 1];
    int cigarStrLen1 = 0;
    int bestMismatchNum1 = 0;
    int bestGapOpen1 = 0;
    int bestGapExtend1 = 0;
    int avg_mismatch_qual1;
    int bestEditDist1 = 0;
    char mdStr2[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen2 = 0;
    char cigarStr2[MAX_READ_LENGTH + 1];
    int cigarStrLen2 = 0;
    int bestMismatchNum2 = 0;
    int bestGapOpen2 = 0;
    int bestGapExtend2 = 0;
    int avg_mismatch_qual2;
    int bestEditDist2 = 0;

    int boundTrim1 = 0, boundTrim2 = 0;
    char * newCigar1 = NULL, *newCigar2 = NULL;
    int r1 = readlen1, r2 = readlen2;

    // int bestHitNum = 0;
    // int secBestHitNum = 0;

    if ( bestResult != NULL )
    {
        // collect the information of the first alignment
        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            best_strand1 = bestResult->strand_1;

            if ( bestResult->whichFromDP == 0 )
            {
                best_cigar1 = bestResult->cigarString;
            }
        }

        // collect the information of the second alignment
        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            best_strand2 = bestResult->strand_2;

            if ( bestResult->whichFromDP == 1 )
            {
                best_cigar2 = bestResult->cigarString;
            }
        }

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            best_insert = bestResult->insertSize;
        }

        // obtain the corresponding normal cigar string, # of gap open, # of gap extension, # of mismatch,
        // and the average quality values in the mismatch positions
        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            if ( bestResult->whichFromDP == 0 )
            {
                // to convert the the best special_cigar into normal cigar string
                char * tempCigar = NULL;
                boundTrim1 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen1, bestResult->algnmt_1, best_cigar1, &tp_1, &chr_1, &tempCigar );
                cigarStrLen1 = boundTrim1 ? convertToCigarStr ( tempCigar, cigarStr1 ) : convertToCigarStr ( best_cigar1, cigarStr1 );
                r1 = readLengthWithCigar ( best_cigar1 );
                // to collect the mismatch information including MD string and # of mismatches
                mdStrLen1 = getMisInfoForDP ( hsp, query1, qualities1, readlen1, bestResult->algnmt_1, best_strand1, boundTrim1 ? tempCigar : best_cigar1, mdStr1, &bestMismatchNum1, &bestGapOpen1, &bestGapExtend1, &avg_mismatch_qual1, boundTrim1 );
                bestEditDist1 = bestResult->editdist;

                if ( boundTrim1 )
                {
                    free ( tempCigar );
                    bestEditDist1 = bestMismatchNum1 + bestGapExtend1;
                }
            }
            else
            {
                boundTrim1 = getChrAndPosWithBoundaryCheck ( qInput, readlen1, bestResult->algnmt_1, &tp_1, &chr_1, &newCigar1 );
                // get the MD string
                mdStrLen1 = getMdStr ( hsp, query1, qualities1, readlen1, bestResult->algnmt_1, best_strand1,
                                       bestResult->score_1, mdStr1, &avg_mismatch_qual1, boundTrim1 );
                bestMismatchNum1 = bestResult->score_1;

                if ( boundTrim1 && bestMismatchNum1 )
                {
                    bestMismatchNum1 = 0;

                    for ( int ii = 0; ii < mdStrLen1; ++ii )
                    { bestMismatchNum1 += ( mdStr1[ii] > '9' ); } // A C G or T
                }

                bestEditDist1 = bestMismatchNum1;
            }
        }

        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            if ( bestResult->whichFromDP == 1 )
            {
                // to convert the the best special_cigar into normal cigar string
                char * tempCigar = NULL;
                boundTrim2 = getChrAndPosWithBoundaryCheckDP ( qInput, readlen2, bestResult->algnmt_2, best_cigar2, &tp_2, &chr_2, &tempCigar );
                cigarStrLen2 = boundTrim2 ? convertToCigarStr ( tempCigar, cigarStr2 ) : convertToCigarStr ( best_cigar2, cigarStr2 );
                r2 = readLengthWithCigar ( best_cigar2 );
                // to collect the mismatch information including MD string and # of mismatches
                mdStrLen2 = getMisInfoForDP ( hsp, query2, qualities2, readlen2, bestResult->algnmt_2, best_strand2, boundTrim2 ? tempCigar : best_cigar2, mdStr2, &bestMismatchNum2, &bestGapOpen2, &bestGapExtend2, &avg_mismatch_qual2, boundTrim2 );
                bestEditDist2 = bestResult->editdist;

                if ( boundTrim2 )
                {
                    free ( tempCigar );
                    bestEditDist2 = bestMismatchNum2 + bestGapExtend2;
                }
            }
            else
            {
                boundTrim2 = getChrAndPosWithBoundaryCheck ( qInput, readlen2, bestResult->algnmt_2, &tp_2, &chr_2, &newCigar2 );
                // get the MD string
                mdStrLen2 = getMdStr ( hsp, query2, qualities2, readlen2, bestResult->algnmt_2, best_strand2,
                                       bestResult->score_2, mdStr2, &avg_mismatch_qual2, boundTrim2 );
                bestMismatchNum2 = bestResult->score_2;

                if ( boundTrim2 && bestMismatchNum2 )
                {
                    bestMismatchNum2 = 0;

                    for ( int ii = 0; ii < mdStrLen2; ++ii )
                    { bestMismatchNum2 += ( mdStr2[ii] > '9' ); } // A C G or T
                }

                bestEditDist2 = bestMismatchNum2;
            }
        }

        // check for read-through reads
        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            unsigned long long adjust1 = bestResult->algnmt_1 + ( boundTrim1 > 0 ? boundTrim1 : 0 );
            unsigned long long adjust2 = bestResult->algnmt_2 + ( boundTrim2 > 0 ? boundTrim2 : 0 );

            if ( ( bestResult->strand_1 == QUERY_POS_STRAND &&
                    ( adjust1 > adjust2 || adjust1 + r1 > adjust2 + r2 ) ) ||
                    ( bestResult->strand_1 == QUERY_NEG_STRAND &&
                      ( adjust2 > adjust1 || adjust2 + r2 > adjust1 + r1 ) ) )
            {
                if ( bestMismatchNum1 <= bestMismatchNum2 )
                {
                    bestResult->algnmt_2 = 0xFFFFFFFF;
                    tp_2 = chr_2 = 0;
                }
                else
                {
                    bestResult->algnmt_1 = 0xFFFFFFFF;
                    tp_1 = chr_1 = 0;
                }

                best_insert = 0;
            }
        }

        int bestPairNum = 0;
        int bestPairScore1 = 0;
        int bestPairScore2 = 0;
        int secBestPairNum = 0;
        int secBestPairScore1;
        int secBestPairScore2;

        // obtain the number of pairs with same mismatch # and dp score
        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            bestPairNum = 1;
            bestPairScore1 = bestResult->score_1;
            bestPairScore2 = bestResult->score_2;

            if ( num > 1 )
            {
                for ( i = start; i < start + num; i++ )
                {
                    if ( & ( algnResult[i] ) == bestResult )
                    { continue; }

                    if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                    { continue; }

                    if ( algnResult[i].score_1 == bestPairScore1 && algnResult[i].score_2 == bestPairScore2 )
                    {
                        bestPairNum++;
                    }
                    else
                    {
                        if ( secBestPairNum == 0 )
                        {
                            secBestPairScore1 = algnResult[i].score_1;
                            secBestPairScore2 = algnResult[i].score_2;
                            secBestPairNum = 1;
                        }
                        else if ( bestResult->whichFromDP == 0 )
                        {
                            // first alignment is from DP
                            if ( algnResult[i].score_2 < secBestPairScore2 )
                            {
                                secBestPairScore1 = algnResult[i].score_1;
                                secBestPairScore2 = algnResult[i].score_2;
                                secBestPairNum = 1;
                            }
                            else if ( ( algnResult[i].score_2 == secBestPairScore2 ) && ( algnResult[i].score_1 > secBestPairScore1 ) )
                            {
                                secBestPairScore1 = algnResult[i].score_1;
                                secBestPairNum = 1;
                            }
                            else if ( ( algnResult[i].score_2 == secBestPairScore2 ) && ( algnResult[i].score_1 == secBestPairScore1 ) )
                            {
                                secBestPairNum++;
                            }
                        }
                        else if ( bestResult->whichFromDP == 1 )
                        {
                            // second alignment is from DP
                            if ( algnResult[i].score_1 < secBestPairScore1 )
                            {
                                secBestPairScore1 = algnResult[i].score_1;
                                secBestPairScore2 = algnResult[i].score_2;
                                secBestPairNum = 1;
                            }
                            else if ( ( algnResult[i].score_1 == secBestPairScore1 ) && ( algnResult[i].score_2 > secBestPairScore2 ) )
                            {
                                secBestPairScore2 = algnResult[i].score_2;
                                secBestPairNum = 1;
                            }
                            else if ( ( algnResult[i].score_1 == secBestPairScore1 ) && ( algnResult[i].score_2 == secBestPairScore2 ) )
                            {
                                secBestPairNum++;
                            }
                        }
                    }
                }
            }
        }

        // obtain the number of pairs with similar score to the best score (i.e. difference <= 1 mismatch score)
        int numSimilarBestPairs = 0;

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            numSimilarBestPairs++;

            if ( num > 1 )
            {
                if ( bestResult->whichFromDP == 0 )
                {
                    // first alignment is from DP
                    for ( i = start; i < start + num; i++ )
                    {
                        if ( & ( algnResult[i] ) == bestResult )
                        { continue; }

                        if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                        { continue; }

                        if ( algnResult[i].score_1 >= bestPairScore1 + hspaux->dpMisMatchScore && algnResult[i].score_2 <= bestPairScore2 + 1 )
                        {
                            numSimilarBestPairs++;
                        }
                    }
                }
                else if ( bestResult->whichFromDP == 1 )
                {
                    // second alignment is from DP
                    for ( i = start; i < start + num; i++ )
                    {
                        if ( & ( algnResult[i] ) == bestResult )
                        { continue; }

                        if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                        { continue; }

                        if ( algnResult[i].score_1 <= bestPairScore1 + 1 && algnResult[i].score_2 >= bestPairScore2 + hspaux->dpMisMatchScore )
                        {
                            numSimilarBestPairs++;
                        }
                    }
                }
                else
                {
                    // both alignment are not from DP
                    // should not happend
                    for ( i = start; i < start + num; i++ )
                    {
                        if ( & ( algnResult[i] ) == bestResult )
                        { continue; }

                        if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                        { continue; }

                        if ( algnResult[i].score_1 <= bestPairScore1 + 1 && algnResult[i].score_2 <= bestPairScore2 + 1 )
                        {
                            numSimilarBestPairs++;
                        }
                    }
                }
            }
        }

        // TODO cx the logic above is so WTF...

        // obtain the bestHitNum, bestScore, secBestHitNum and secBestScore for the first read
        int bestHitNum1 = 0;
        int secBestHitNum1 = 0;
        int bestScore1 = 0;
        int secBestScore1 = 0;
        int isBestHit1 = 1;
        unsigned int bestPos1 = 0xFFFFFFFF;
        unsigned int secBestPos1 = 0xFFFFFFFF;

        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            bestScore1 = bestResult->score_1; // number of mismatches or score from DP
            bestPos1 = bestResult->algnmt_1; // alignment position
            bestHitNum1 = 1;

            // bestHitNum1 = bestResult->num_sameScore;
            if ( bestResult->whichFromDP != 0 )
            {
                secBestScore1 = 999; // mismatch #
            }
            else
            {
                secBestScore1 = 0; // dp score
            }
        }

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && num > 1 )
        {
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                { continue; }

                if ( algnResult[i].whichFromDP != 0 )
                {
                    // mismatch, thus smaller is better
                    if ( algnResult[i].score_1 <= bestScore1 )
                    {
                        if ( algnResult[i].score_1 == bestScore1 )
                        {
                            if ( algnResult[i].algnmt_1 != bestPos1 )
                            { bestHitNum1++; }
                        }
                        else
                        {
                            secBestScore1 = bestScore1;
                            secBestHitNum1 = bestHitNum1;
                            secBestPos1 = bestPos1;
                            bestScore1 = algnResult[i].score_1;
                            bestHitNum1 = 1;
                            bestPos1 = algnResult[i].algnmt_1;
                            isBestHit1 = 0;
                        }
                    }
                    else if ( algnResult[i].score_1 <= secBestScore1 )
                    {
                        if ( algnResult[i].score_1 == secBestScore1 )
                        {
                            if ( algnResult[i].algnmt_1 != secBestPos1 )
                            { secBestHitNum1++; }
                        }
                        else
                        {
                            secBestScore1 = algnResult[i].score_1;
                            secBestPos1 = algnResult[i].algnmt_1;
                            secBestHitNum1 = 1;
                        }
                    }
                }
                else
                {
                    // DP score, thus bigger is better
                    if ( algnResult[i].score_1 >= bestScore1 )
                    {
                        if ( algnResult[i].score_1 == bestScore1 )
                        {
                            if ( algnResult[i].algnmt_1 != bestPos1 )
                                // bestHitNum1++;
                            { bestHitNum1 += algnResult[i].num_sameScore; }
                        }
                        else
                        {
                            secBestScore1 = bestScore1;
                            secBestHitNum1 = bestHitNum1;
                            secBestPos1 = bestPos1;
                            bestScore1 = algnResult[i].score_1;
                            // bestHitNum1 = 1;
                            bestPos1 = algnResult[i].algnmt_1;
                            bestHitNum1 = algnResult[i].num_sameScore;
                            isBestHit1 = 0;
                        }
                    }
                    else if ( algnResult[i].score_1 >= secBestScore1 )
                    {
                        if ( algnResult[i].score_1 == secBestScore1 )
                        {
                            if ( algnResult[i].algnmt_1 != secBestPos1 )
                                // secBestHitNum1++;
                            { secBestHitNum1 += algnResult[i].num_sameScore; }
                        }
                        else
                        {
                            secBestScore1 = algnResult[i].score_1;
                            secBestPos1 = algnResult[i].algnmt_1;
                            secBestHitNum1 = algnResult[i].num_sameScore;
                            // secBestHitNum1 = 1;
                        }
                    }
                }
            }
        }

        if ( hspaux->x0_array[bestResult->readID] > 0 )
        {
            if ( bestResult->whichFromDP != 0 )
            {
                // this result is from SOAP3
                // bestHitNum1 = (bestHitNum1<hspaux->x0_array[bestResult->readID])?bestHitNum1:hspaux->x0_array[bestResult->readID];
                // secBestHitNum1 = (secBestHitNum1<hspaux->x1_array[bestResult->readID])?secBestHitNum1:hspaux->x1_array[bestResult->readID];
                bestHitNum1 = hspaux->x0_array[bestResult->readID];
                secBestHitNum1 = hspaux->x1_array[bestResult->readID];

                if ( hspaux->mismatch_array[bestResult->readID] < bestResult->score_1 )
                { isBestHit1 = 0; }
            }
            else
            {
                // this result is not from SOAP3
                int x0_score = hspaux->mismatch_array[bestResult->readID] * hspaux->dpMisMatchScore
                               + ( readlen1 - hspaux->mismatch_array[bestResult->readID] ) * hspaux->dpMatchScore;

                if ( x0_score >= bestScore1 )
                {
                    bestHitNum1 = ( bestHitNum1 > hspaux->x0_array[bestResult->readID] ) ? bestHitNum1 : hspaux->x0_array[bestResult->readID];
                    secBestHitNum1 = ( secBestHitNum1 > hspaux->x1_array[bestResult->readID] ) ? secBestHitNum1 : hspaux->x1_array[bestResult->readID];

                    if ( x0_score > bestScore1 )
                    { isBestHit1 = 0; }
                }
            }
        }

        // TODO cx the logic above is so WTF

        // obtain the bestHitNum, bestScore, secBestHitNum and secBestScore for the second read
        int bestHitNum2 = 0;
        int secBestHitNum2 = 0;
        int bestScore2 = 0;
        int secBestScore2 = 0;
        int isBestHit2 = 1;
        unsigned int bestPos2 = 0xFFFFFFFF;
        unsigned int secBestPos2 = 0xFFFFFFFF;

        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            bestScore2 = bestResult->score_2; // number of mismatches or score from DP
            bestPos2 = bestResult->algnmt_2; // alignment position
            bestHitNum2 = 1;

            if ( bestResult->whichFromDP != 1 )
            {
                secBestScore2 = 999; // mismatch #
            }
            else
            {
                secBestScore2 = 0; // dp score
            }
        }

        if ( bestResult->algnmt_2 != 0xFFFFFFFF && num > 1 )
        {
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                { continue; }

                if ( algnResult[i].whichFromDP != 1 )
                {
                    // mismatch, thus smaller is better
                    if ( algnResult[i].score_2 <= bestScore2 )
                    {
                        if ( algnResult[i].score_2 == bestScore2 )
                        {
                            if ( algnResult[i].algnmt_2 != bestPos2 )
                            { bestHitNum2++; }
                        }
                        else
                        {
                            secBestScore2 = bestScore2;
                            secBestHitNum2 = bestHitNum2;
                            secBestPos2 = bestPos2;
                            bestScore2 = algnResult[i].score_2;
                            bestHitNum2 = 1;
                            bestPos2 = algnResult[i].algnmt_2;
                            isBestHit2 = 0;
                        }
                    }
                    else if ( algnResult[i].score_2 <= secBestScore2 )
                    {
                        if ( algnResult[i].score_2 == secBestScore2 )
                        {
                            if ( algnResult[i].algnmt_2 != secBestPos2 )
                            { secBestHitNum2++; }
                        }
                        else
                        {
                            secBestScore2 = algnResult[i].score_2;
                            secBestPos2 = algnResult[i].algnmt_2;
                            secBestHitNum2 = 1;
                        }
                    }
                }
                else
                {
                    // DP score, thus bigger is better
                    if ( algnResult[i].score_2 >= bestScore2 )
                    {
                        if ( algnResult[i].score_2 == bestScore2 )
                        {
                            if ( algnResult[i].algnmt_2 != bestPos2 )
                                // bestHitNum2++;
                            { bestHitNum2 += algnResult[i].num_sameScore; }
                        }
                        else
                        {
                            secBestScore2 = bestScore2;
                            secBestHitNum2 = bestHitNum2;
                            secBestPos2 = bestPos2;
                            bestScore2 = algnResult[i].score_2;
                            // bestHitNum2 = 1;
                            bestHitNum2 = algnResult[i].num_sameScore;
                            isBestHit2 = 0;
                        }
                    }
                    else if ( algnResult[i].score_2 >= secBestScore2 )
                    {
                        if ( algnResult[i].score_2 == secBestScore2 )
                        {
                            if ( algnResult[i].algnmt_2 != secBestPos2 )
                                // secBestHitNum2++;
                            { secBestHitNum2 += algnResult[i].num_sameScore; }
                        }
                        else
                        {
                            secBestScore2 = algnResult[i].score_2;
                            secBestPos2 = algnResult[i].algnmt_2;
                            // secBestHitNum2 = 1;
                            secBestHitNum2 = algnResult[i].num_sameScore;
                        }
                    }
                }
            }
        }

        if ( hspaux->x0_array[bestResult->readID + 1] > 1 )
        {
            if ( bestResult->whichFromDP != 1 )
            {
                // this result is from SOAP3
                // bestHitNum2 = (bestHitNum2<hspaux->x0_array[bestResult->readID+1])?bestHitNum2:hspaux->x0_array[bestResult->readID+1];
                // secBestHitNum2 = (secBestHitNum2<hspaux->x1_array[bestResult->readID+1])?secBestHitNum2:hspaux->x1_array[bestResult->readID+1];
                bestHitNum2 = hspaux->x0_array[bestResult->readID + 1];
                secBestHitNum2 = hspaux->x1_array[bestResult->readID + 1];

                if ( hspaux->mismatch_array[bestResult->readID + 1] < bestResult->score_2 )
                { isBestHit2 = 0; }
            }
            else
            {
                // this result is not from SOAP3
                int x0_score = hspaux->mismatch_array[bestResult->readID + 1] * hspaux->dpMisMatchScore
                               + ( readlen2 - hspaux->mismatch_array[bestResult->readID + 1] ) * hspaux->dpMatchScore;

                if ( x0_score >= bestScore2 )
                {
                    bestHitNum2 = ( bestHitNum2 > hspaux->x0_array[bestResult->readID + 1] ) ? bestHitNum2 : hspaux->x0_array[bestResult->readID + 1];
                    secBestHitNum2 = ( secBestHitNum2 > hspaux->x1_array[bestResult->readID + 1] ) ? secBestHitNum2 : hspaux->x1_array[bestResult->readID + 1];

                    if ( x0_score > bestScore2 )
                    { isBestHit2 = 0; }
                }
            }
        }

        // compute the mapping quality score
        if ( bestResult->algnmt_1 != 0xFFFFFFFF && bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            if ( hspaux->alignmentType != OUTPUT_ALL_VALID && hspaux->alignmentType != OUTPUT_ALL_BEST )
            {
                mapq1 = mapq2 = SAM_MAPQ_UNAVAILABLE;
            }
            else
            {
                if ( hspaux->bwaLikeScore )
                {
                    int op_score = 0;
                    int subop_score = 0;

                    if ( bestResult->whichFromDP == 0 )
                    {
                        op_score = bestPairScore2 * hspaux->dpMisMatchScore + ( readlen2 - bestPairScore2 ) * hspaux->dpMatchScore + bestPairScore1;

                        if ( secBestPairNum > 0 )
                        { subop_score = secBestPairScore2 * hspaux->dpMisMatchScore + ( readlen2 - secBestPairScore2 ) * hspaux->dpMatchScore + secBestPairScore1; }
                    }
                    else
                    {
                        op_score = bestPairScore1 * hspaux->dpMisMatchScore + ( readlen1 - bestPairScore1 ) * hspaux->dpMatchScore + bestPairScore2;

                        if ( secBestPairNum > 0 )
                        { subop_score = secBestPairScore1 * hspaux->dpMisMatchScore + ( readlen1 - secBestPairScore1 ) * hspaux->dpMatchScore + secBestPairScore2; }
                    }

                    bwaLikePairQualScore ( bestHitNum1, secBestHitNum1, bestHitNum2, secBestHitNum2, hspaux->g_log_n, op_score, bestPairNum, subop_score, secBestPairNum, readlen1, readlen2, &mapq1, &mapq2 );
                }
                else
                {
                    int mapping_qual_score1 = 0;
                    int mapping_qual_score2 = 0;

                    if ( bestResult->whichFromDP == 0 )
                    {
                        //mapping_qual_score1 = getMapQualScoreForDP2(bestResult->score_1, readlen1*hspaux->dpMatchScore, (hspaux->isFastq==1)?avg_mismatch_qual1:20, bestHitNum1, secBestHitNum1, bestScore1, secBestScore1, isBestHit1, bestResult->num_similarScore, numSimilarBestPairs);
                        mapping_qual_score1 = getMapQualScoreForDP2 ( bestResult->score_1, readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, bestHitNum1, secBestHitNum1, bestScore1, secBestScore1, isBestHit1, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                        mapping_qual_score2 = getMapQualScore2 ( bestResult->score_2, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, bestHitNum2, secBestHitNum2, isBestHit2, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                    }
                    else
                    {
                        mapping_qual_score1 = getMapQualScore2 ( bestResult->score_1, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, bestHitNum1, secBestHitNum1, isBestHit1, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                        //mapping_qual_score2 = getMapQualScoreForDP2(bestResult->score_2, readlen2*hspaux->dpMatchScore, (hspaux->isFastq==1)?avg_mismatch_qual2:20, bestHitNum2, secBestHitNum2, bestScore2, secBestScore2, isBestHit2, bestResult->num_similarScore, numSimilarBestPairs);
                        mapping_qual_score2 = getMapQualScoreForDP2 ( bestResult->score_2, readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, bestHitNum2, secBestHitNum2, bestScore2, secBestScore2, isBestHit2, numSimilarBestPairs, hspaux->maxMAPQ, hspaux->minMAPQ );
                    }

                    mapq1 = mapq2 = getMapQualScoreForPair ( mapping_qual_score1, mapping_qual_score2 );
                    // mapping_qual_score = (mapping_qual_score1 + mapping_qual_score2) / 2;
                }

                if ( boundTrim1 ) { mapq1 = 0; }

                if ( boundTrim2 ) { mapq2 = 0; }
            }
        }

        // TODO cx above T^T sup9 logic finished

        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( bestResult->algnmt_1 != 0xFFFFFFFF && num > 1 )
        {
            // Build the XA:Z tag
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                { continue; }

                if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) &&
                        ( algnResult[i].score_1 != bestPairScore1 ||
                          algnResult[i].score_2 != bestPairScore2 ) )
                { continue; }

                /// TODO cx others
                getChrAndPos ( qInput, algnResult[i].algnmt_1, &curr_tp, &curr_chr );
                curr_strand = algnResult[i].strand_1;

                if ( algnResult[i].whichFromDP != 0 )
                {
                    char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                    memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                    int pos = strlen ( chr_name );
                    curr_occStr[pos++] = ',';
                    curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                    pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( readlen1, & ( curr_occStr[pos] ) ); // cigar string
                    curr_occStr[pos++] = 'M';
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( algnResult[i].score_1, & ( curr_occStr[pos] ) ); // number of mismatches
                    curr_occStr[pos++] = ';';
                    curr_occStr[pos] = '\0';
                    curr_len = pos;
                }
                else
                {
                    char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                    memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                    int pos = strlen ( chr_name );
                    curr_occStr[pos++] = ',';
                    curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                    pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                    curr_occStr[pos++] = ',';
                    pos += convertToCigarStr ( algnResult[i].cigarString, & ( curr_occStr[pos] ) ); // cigar string
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( algnResult[i].editdist, & ( curr_occStr[pos] ) ); // edit distance
                    curr_occStr[pos++] = ';';
                    curr_occStr[pos] = '\0';
                    curr_len = pos;
                }

                appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            }
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST )
        {
            bestHitNum1 = -1;
            secBestHitNum1 = -1;
        }
        else if ( hspaux->alignmentType != OUTPUT_ALL_VALID && hspaux->alignmentType != OUTPUT_ALL_BEST )
        {
            secBestHitNum1 = -1;
        }

        if ( bestResult->algnmt_1 != 0xFFFFFFFF )
        {
            if ( bestResult->algnmt_2 != 0xFFFFFFFF )
            {
                initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1,
                                       qualities1, best_strand1, xazArray->charStr, xazArray->length,
                                       newCigar1 ? newCigar1 : ( bestResult->whichFromDP == 0 ) ? cigarStr1 : NULL, 0,
                                       bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestHitNum1, secBestHitNum1, bestGapOpen1, bestGapExtend1,
                                       mdStr1, mdStrLen1, mapq1, hspaux->readGroup, hspaux->isPrintMDNM );
            }
            else
            {
                // should not come to here (after deep dp is implemented)
                if ( bestResult->whichFromDP == 0 )
                { mapq1 = getMapQualScoreForDP ( bestHitNum1, bestResult->score_1, readlen1 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, hspaux->maxMAPQ, hspaux->minMAPQ ); }
                else
                { mapq1 = getMapQualScore ( bestHitNum1, bestResult->score_1, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual1 : 20, hspaux->maxMAPQ, hspaux->minMAPQ ); }

                if ( boundTrim1 ) { mapq1 = 0; }

                initializeSAMAlgnmt2 ( samAlgnmt, readlen1, queryName1, query1,
                                       qualities1, best_strand1, xazArray->charStr, xazArray->length,
                                       newCigar1 ? newCigar1 : ( bestResult->whichFromDP == 0 ) ? cigarStr1 : NULL, 0,
                                       bestMismatchNum1, bestMismatchNum1 + bestGapExtend1, bestHitNum1, secBestHitNum1, bestGapOpen1, bestGapExtend1,
                                       mdStr1, mdStrLen1, mapq1, hspaux->readGroup, hspaux->isPrintMDNM );
            }
        }
        else
        {
            initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1,
                                  qualities1, best_strand1, xazArray->charStr,
                                  xazArray->length, NULL, 1, hspaux->readGroup );
        }

        // compute the value of the flag
        samFlag = 1;

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->algnmt_2 != 0xFFFFFFFF ) )
        { samFlag |= SAM_FLAG_PROPER_PAIR; }

        samFlag |= SAM_FLAG_FIRST_IN_PAIR;

        if ( bestResult->algnmt_1 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestResult->algnmt_2 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->strand_1 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( ( bestResult->algnmt_2 != 0xFFFFFFFF ) && ( bestResult->strand_2 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;
        samAlgnmt->core.mpos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;

        //      samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        //      samAlgnmt->core.tid = chr_1 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        //      samAlgnmt->core.pos = tp_1 - 1;                       //SAM: 0-based leftmost coordinate
        //      samAlgnmt->core.mtid = chr_2 - 1;
        //      samAlgnmt->core.mpos = tp_2 - 1;

        if ( best_insert > 0 )
        {
            if ( tp_1 > tp_2 )
            { samAlgnmt->core.isize = - ( tp_1 + r1 - tp_2 ); }
            else
            { samAlgnmt->core.isize = ( tp_2 + r2 - tp_1 ); }
        }
        else
        {
            samAlgnmt->core.isize = 0;
        }

        samwrite ( samFilePtr, samAlgnmt );
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        DynamicUint8ArrayReset ( xazArray );

        if ( bestResult->algnmt_2 != 0xFFFFFFFF && num > 1 )
        {
            // Build the XA:Z tag
            for ( i = start; i < start + num; i++ )
            {
                if ( & ( algnResult[i] ) == bestResult )
                { continue; }

                if ( algnResult[i].whichFromDP != bestResult->whichFromDP )
                { continue; }

                if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) &&
                        ( algnResult[i].score_1 != bestPairScore1 ||
                          algnResult[i].score_2 != bestPairScore2 ) )
                { continue; }

                getChrAndPos ( qInput, algnResult[i].algnmt_2, &curr_tp, &curr_chr );
                curr_strand = algnResult[i].strand_2;

                if ( algnResult[i].whichFromDP != 1 )
                {
                    char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                    memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                    int pos = strlen ( chr_name );
                    curr_occStr[pos++] = ',';
                    curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                    pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( readlen2, & ( curr_occStr[pos] ) ); // cigar string
                    curr_occStr[pos++] = 'M';
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( algnResult[i].score_2, & ( curr_occStr[pos] ) ); // number of mismatches
                    curr_occStr[pos++] = ';';
                    curr_occStr[pos] = '\0';
                    curr_len = pos;
                }
                else
                {
                    char * chr_name = samFilePtr->header->target_name[curr_chr - 1]; // chromosome name
                    memcpy ( curr_occStr, chr_name, strlen ( chr_name ) );
                    int pos = strlen ( chr_name );
                    curr_occStr[pos++] = ',';
                    curr_occStr[pos++] = ( curr_strand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                    pos += writeULLToStr ( curr_tp, & ( curr_occStr[pos] ) ); // chromosome position
                    curr_occStr[pos++] = ',';
                    pos += convertToCigarStr ( algnResult[i].cigarString, & ( curr_occStr[pos] ) ); // cigar string
                    curr_occStr[pos++] = ',';
                    pos += writeNumToStr ( algnResult[i].editdist, & ( curr_occStr[pos] ) ); // edit distance
                    curr_occStr[pos++] = ';';
                    curr_occStr[pos] = '\0';
                    curr_len = pos;
                }

                appendStringToUint8Array ( xazArray, curr_occStr, curr_len );
            }
        }

        if ( hspaux->alignmentType == OUTPUT_RANDOM_BEST )
        {
            bestHitNum2 = -1;
            secBestHitNum2 = -1;
        }
        else if ( hspaux->alignmentType != OUTPUT_ALL_VALID && hspaux->alignmentType != OUTPUT_ALL_BEST )
        {
            secBestHitNum2 = -1;
        }

        if ( bestResult->algnmt_2 != 0xFFFFFFFF )
        {
            if ( bestResult->algnmt_1 != 0xFFFFFFFF )
            {
                initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2,
                                       qualities2, best_strand2, xazArray->charStr, xazArray->length,
                                       newCigar2 ? newCigar2 : ( bestResult->whichFromDP == 1 ) ? cigarStr2 : NULL, 0,
                                       bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestHitNum2, secBestHitNum2, bestGapOpen2, bestGapExtend2,
                                       mdStr2, mdStrLen2, mapq2, hspaux->readGroup, hspaux->isPrintMDNM );
            }
            else
            {
                // should not come to here (after deep dp is implemented)
                // int mapping_qual_score2 = 0;
                if ( bestResult->whichFromDP == 1 )
                { mapq2 = getMapQualScoreForDP ( bestHitNum2, bestResult->score_2, readlen2 * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, hspaux->maxMAPQ, hspaux->minMAPQ ); }
                else
                { mapq2 = getMapQualScore ( bestHitNum2, bestResult->score_2, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual2 : 20, hspaux->maxMAPQ, hspaux->minMAPQ ); }

                if ( boundTrim2 ) { mapq2 = 0; }

                initializeSAMAlgnmt2 ( samAlgnmt, readlen2, queryName2, query2,
                                       qualities2, best_strand2, xazArray->charStr, xazArray->length,
                                       newCigar2 ? newCigar2 : ( bestResult->whichFromDP == 1 ) ? cigarStr2 : NULL, 0,
                                       bestMismatchNum2, bestMismatchNum2 + bestGapExtend2, bestHitNum2, secBestHitNum2, bestGapOpen2, bestGapExtend2,
                                       mdStr2, mdStrLen2, mapq2, hspaux->readGroup, hspaux->isPrintMDNM );
            }
        }
        else
        {
            initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2,
                                  qualities2, best_strand2, xazArray->charStr,
                                  xazArray->length, NULL, 1, hspaux->readGroup );
        }

        // compute the value of the flag
        samFlag = 1;

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->algnmt_2 != 0xFFFFFFFF ) )
        { samFlag |= SAM_FLAG_PROPER_PAIR; }

        samFlag |= SAM_FLAG_SECOND_IN_PAIR;

        if ( bestResult->algnmt_2 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_READ_UNMAPPED;
        }

        if ( bestResult->algnmt_1 == 0xFFFFFFFF )
        {
            samFlag |= SAM_FLAG_MATE_UNMAPPED;
        }

        if ( ( bestResult->algnmt_2 != 0xFFFFFFFF ) && ( bestResult->strand_2 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        if ( ( bestResult->algnmt_1 != 0xFFFFFFFF ) && ( bestResult->strand_1 == QUERY_NEG_STRAND ) )
        {
            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
        }

        samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        samAlgnmt->core.tid = ( chr_2 == 0 ) ? ( chr_1 == 0 ? -1 : chr_1 - 1 ) : chr_2 - 1;    //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = ( tp_2 == 0 ) ? ( tp_1 == 0 ? -1 : tp_1 - 1 ) : tp_2 - 1;     //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = ( chr_1 == 0 ) ? ( chr_2 == 0 ? -1 : chr_2 - 1 ) : chr_1 - 1;
        samAlgnmt->core.mpos = ( tp_1 == 0 ) ? ( tp_2 == 0 ? -1 : tp_2 - 1 ) : tp_1 - 1;

        //      samAlgnmt->core.flag = samFlag;                       //SAM: bitwise flag
        //      samAlgnmt->core.tid = chr_2 - 1;                      //SAM: chromosome ID, defined by bam_header_t
        //      samAlgnmt->core.pos = tp_2 - 1;                       //SAM: 0-based leftmost coordinate
        //      samAlgnmt->core.mtid = chr_1 - 1;
        //      samAlgnmt->core.mpos = tp_1 - 1;

        if ( best_insert > 0 )
        {
            if ( tp_2 > tp_1 )
            { samAlgnmt->core.isize = - ( tp_2 + r2 - tp_1 ); }
            else
            { samAlgnmt->core.isize = tp_1 + r1 - tp_2; }
        }
        else
        {
            samAlgnmt->core.isize = 0;
        }

        samwrite ( samFilePtr, samAlgnmt );
    }
    else
    {
        //-------------------------------------------//
        // report the first alignment                //
        //-------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen1, queryName1, query1,
                              qualities1, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_FIRST_IN_PAIR;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
        //--------------------------------------------//
        // report the second alignment                //
        //--------------------------------------------//
        initializeSAMAlgnmt ( samAlgnmt, readlen2, queryName2, query2,
                              qualities2, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
        // compute the value of the flag
        samFlag = 1;
        samFlag |= SAM_FLAG_SECOND_IN_PAIR;
        samAlgnmt->core.flag = samFlag;                  //SAM: bitwise flag
        samAlgnmt->core.tid = -1;                        //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                        //SAM: 0-based leftmost coordinate
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
    }

    if ( newCigar1 ) { free ( newCigar1 ); }

    if ( newCigar2 ) { free ( newCigar2 ); }
}

void OCCOutputSAMAPI ( SRAQueryInput * qInput, OCCList * occ_list,
                       DynamicUint8Array * xazArray, int readlen, int reportType )
{
    // cx this function... X_x
    //  fprintf(stderr,"OCCOutputSAMAPI\n");

    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // For outputting single-read alignment in SAM format
    // Output the best alignment, and for the rest, append the results in the tag XA:Z
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( qSetting->occ->SAMOutBuffer );
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    unsigned long long k = 0;
    unsigned int bestNumMisMatch, currNumMisMatch;
    unsigned long long bestTP, currTP;
    unsigned short bestChr, currChr;
    char bestOccStr[500], currOccStr[500];
    char bestStrand = QUERY_POS_STRAND;
    char currStrand;
    int bestStrLen, currStrLen;
    int bestHitNum = 0;
    int secBestHitNum = 0;
    char mdStr[MAX_READ_LENGTH * 2 + 1];
    int avg_mismatch_qual;
    int map_qual_score = 0;
    int bestOccIndex = -1;

    // obtain the bestHitNum, bestNumMisMatch and secBestHitNum for the read
    int boundTrim = 0;
    char * newCIGAR = NULL;

    if ( occ_list->curr_size > 0 )
    {
        bestOccIndex = 0;
        bestNumMisMatch = occ_list->occ[0].mismatchCount;
        bestHitNum = 1;

        for ( k = 1; k < occ_list->curr_size; k++ )
        {
            if ( occ_list->occ[k].mismatchCount < bestNumMisMatch )
            {

                bestOccIndex = k;
                bestNumMisMatch = occ_list->occ[k].mismatchCount;
                bestHitNum = 1;
            }
            else if ( occ_list->occ[k].mismatchCount == bestNumMisMatch )
            {
                bestHitNum++;

            }
        }

        bestStrand = occ_list->occ[bestOccIndex].strand;

        boundTrim = getChrAndPosWithBoundaryCheck ( qInput, readlen, occ_list->occ[bestOccIndex].ambPosition, &bestTP, &bestChr, NULL ); // check good or not

        if ( boundTrim ) // FOR CROSS-CHROMO BUG. the best alignment sucks... we try to find the best "good" alignment instead
        {
            int storedBestOccIndex = bestOccIndex;
            bestOccIndex = -1;
            bestNumMisMatch = 9999; // hardcode a big number
            bestHitNum = 0;

            for ( k = 0; k < occ_list->curr_size; k++ )
            {
                // check whether it's good
                boundTrim = getChrAndPosWithBoundaryCheck ( qInput, readlen, occ_list->occ[k].ambPosition, &bestTP, &bestChr, NULL );

                if ( ! boundTrim )   // good alignment
                {
                    if ( occ_list->occ[k].mismatchCount < bestNumMisMatch )
                    {
                        bestOccIndex = k;
                        bestNumMisMatch = occ_list->occ[k].mismatchCount;
                        bestHitNum = 1;
                    }
                    else if ( occ_list->occ[k].mismatchCount == bestNumMisMatch )
                    {
                        bestHitNum++;
                    }
                }
            }

            if ( bestOccIndex < 0 ) // damn.. there's no good alignment. so just pick one
            {
                bestOccIndex = storedBestOccIndex;
                bestHitNum = 1;
            }

            boundTrim = getChrAndPosWithBoundaryCheck ( qInput, readlen, occ_list->occ[bestOccIndex].ambPosition, &bestTP, &bestChr, &newCIGAR );
            printf ( "NewCIGAR %x, ambPosition: %u", occ_list->occ[bestOccIndex].ambPosition, newCIGAR);
        } // end CROSS-CHROMO handling

    }

    // updated
    //  secBestHitNum = occ_list->curr_size - bestHitNum;
    secBestHitNum = 0; // counter

    // prepare the XA:Z tag
    if ( occ_list->curr_size > 1 && ! boundTrim )   // if best is cross-chromo, we dont need other bad alignments
    {
        // Build the XA:Z tag
        DynamicUint8ArrayReset ( xazArray );

        for ( k = 0; k < occ_list->curr_size; k++ )
        {
            if ( k == bestOccIndex )
            { continue; }

            if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) && ( occ_list->occ[k].mismatchCount > bestNumMisMatch ) )
            { continue; }

            currNumMisMatch = occ_list->occ[k].mismatchCount;
            currStrand = occ_list->occ[k].strand;
            int boundTrimTemp = getChrAndPosWithBoundaryCheck ( qInput, readlen, occ_list->occ[k].ambPosition, &currTP, &currChr, NULL );

            if ( boundTrimTemp ) { continue; } // don't allow cross-chromo alignments

            char * chr_name = samFilePtr->header->target_name[currChr - 1]; // chromosome name
            memcpy ( currOccStr, chr_name, strlen ( chr_name ) );
            int pos = strlen ( chr_name );
            currOccStr[pos++] = ',';
            currOccStr[pos++] = ( currStrand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
            pos += writeULLToStr ( currTP, & ( currOccStr[pos] ) ); // chromosome position
            currOccStr[pos++] = ',';
            pos += writeNumToStr ( readlen, & ( currOccStr[pos] ) ); // cigar string
            currOccStr[pos++] = 'M';
            currOccStr[pos++] = ',';
            pos += writeNumToStr ( currNumMisMatch, & ( currOccStr[pos] ) ); // number of mismatches
            currOccStr[pos++] = ';';
            currOccStr[pos] = '\0';
            currStrLen = pos;
            appendStringToUint8Array ( xazArray, currOccStr, currStrLen );

            if ( currNumMisMatch > bestNumMisMatch ) { secBestHitNum++; }
        }
    }

    int mdLen = 0;

    if ( bestOccIndex >= 0 )
    {
        mdLen = getMdStr ( hsp, qInfo->ReadCode, qInfo->ReadQuality, qInfo->ReadLength, occ_list->occ[bestOccIndex].ambPosition,
                           bestStrand, bestNumMisMatch, mdStr, &avg_mismatch_qual, boundTrim );

        // map_qual_score = getMapQualScore(bestHitNum, bestNumMisMatch, (hspaux->isFastq==1)?avg_mismatch_qual:20, hspaux->maxMAPQ, hspaux->minMAPQ);
        // recalculate num mismatches
        if ( boundTrim && bestNumMisMatch )
        {
            bestNumMisMatch = 0;

            for ( int ii = 0; ii < mdLen; ++ii )
            { bestNumMisMatch += ( mdStr[ii] > '9' ); } // A C G or T
        }

        map_qual_score = getMapQualScoreSingle ( bestNumMisMatch, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual : 20, bestHitNum, secBestHitNum, hspaux->maxMAPQ, hspaux->minMAPQ, hspaux->bwaLikeScore, hspaux->g_log_n );

    }

    if ( occ_list->curr_size > 1 )
    {
        initializeSAMAlgnmt2 ( samAlgnmt, qInfo->ReadLength, qInfo->ReadName, qInfo->ReadCode,
                               qInfo->ReadQuality, bestStrand, xazArray->charStr, xazArray->length, newCIGAR, 0,
                               bestNumMisMatch, bestNumMisMatch, bestHitNum, secBestHitNum, 0, 0, mdStr, mdLen, map_qual_score, hspaux->readGroup, hspaux->isPrintMDNM );
    }
    else
    {
        initializeSAMAlgnmt2 ( samAlgnmt, qInfo->ReadLength, qInfo->ReadName, qInfo->ReadCode,
                               qInfo->ReadQuality, bestStrand, NULL, 0, newCIGAR, occ_list->curr_size == 0,
                               bestNumMisMatch, bestNumMisMatch, bestHitNum, secBestHitNum, 0, 0, mdStr, mdLen, map_qual_score, hspaux->readGroup, hspaux->isPrintMDNM );
    }

    // compute the flag value and assign the chromosome, position and mapping quality
    if ( bestOccIndex >= 0 )
    {
        unsigned int samFlag = 0;

        if ( bestStrand == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        // ---------------------->  |
        samAlgnmt->core.tid = bestChr - 1;  //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = bestTP - 1;                               //SAM: 0-based leftmost coordinate
        samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
        // ---------------------->  |
    }
    else
    {
        unsigned int samFlag = 0;
        samFlag |= SAM_FLAG_READ_UNMAPPED;
        // ---------------------->  |
        samAlgnmt->core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
        samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
        // ---------------------->  |
    }

    if ( newCIGAR ) { free ( newCIGAR ); printf ( "Free: %x\n", newCIGAR ); }
}

void SingleAnsOutputSAMAPI ( SRAQueryInput * qInput,
                             char strand, unsigned int ambPosition, int bestHitNumOfMismatch, int bestHitNum )
{
    // cx this function... X_x
    //    fprintf(stderr,"SingleAnsOutputSAMAPI\n");

    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // For outputting one answer of single-read alignment in SAM format
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( qSetting->occ->SAMOutBuffer );
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    char mdStr[MAX_READ_LENGTH * 2 + 1];
    unsigned long long TP;
    unsigned short chr;
    int avg_mismatch_qual;
    int map_qual_score = 0;
    getChrAndPos ( qInput, ambPosition,
                   &TP, &chr );
    int mdLen = getMdStr ( hsp, qInfo->ReadCode, qInfo->ReadQuality, qInfo->ReadLength, ambPosition,
                           strand, bestHitNumOfMismatch, mdStr, &avg_mismatch_qual );

    if ( bestHitNum > 0 )
    {
        map_qual_score = getMapQualScore ( bestHitNum, bestHitNumOfMismatch, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual : 20, hspaux->maxMAPQ, hspaux->minMAPQ );
    }
    else
    {
        map_qual_score = SAM_MAPQ_UNAVAILABLE;
    }

    initializeSAMAlgnmt2 ( samAlgnmt, qInfo->ReadLength, qInfo->ReadName, qInfo->ReadCode,
                           qInfo->ReadQuality, strand, NULL, 0, NULL, 0, bestHitNumOfMismatch, bestHitNumOfMismatch,
                           bestHitNum, -1, 0, 0, mdStr, mdLen, map_qual_score, hspaux->readGroup, hspaux->isPrintMDNM );
    // compute the flag value and assign the chromosome, position and mapping quality
    unsigned int samFlag = 0;

    if ( strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
    }

    // ---------------------->  |
    samAlgnmt->core.tid = chr - 1;  //SAM: chromosome ID, defined by bam_header_t
    samAlgnmt->core.pos = TP - 1;                               //SAM: 0-based leftmost coordinate
    samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
    samAlgnmt->core.mtid = -1;
    samAlgnmt->core.mpos = -1;
    samAlgnmt->core.isize = 0;
    samwrite ( samFilePtr, samAlgnmt );
    // ---------------------->  |
}

void noAnsOutputSAMAPI ( SRAQueryInput * qInput )
{
    //    fprintf(stderr,"noAnsOutputSAMAPI\n");
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // For outputting no answer of single-read alignment in SAM format
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( qSetting->occ->SAMOutBuffer );
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    initializeSAMAlgnmt ( samAlgnmt, qInfo->ReadLength, qInfo->ReadName, qInfo->ReadCode,
                          qInfo->ReadQuality, QUERY_POS_STRAND, NULL, 0, NULL, 1, hspaux->readGroup );
    unsigned int samFlag = 0;
    samFlag |= SAM_FLAG_READ_UNMAPPED;
    // ---------------------->  |
    samAlgnmt->core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
    samAlgnmt->core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
    samAlgnmt->core.qual = 0;                                   //SAM: mapping quality
    samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
    samAlgnmt->core.mtid = -1;
    samAlgnmt->core.mpos = -1;
    samAlgnmt->core.isize = 0;
    samwrite ( samFilePtr, samAlgnmt );
    // ---------------------->  |
}


void SingleDPOutputSAMAPI ( SRAQueryInput * qInput, SingleAlgnmtResult * algnResult,
                            unsigned int startIndex, unsigned int numResult,
                            DynamicUint8Array * xazArray )
{
    // cx this function... X_x

    //fprintf(stderr, "SingleDPOutputSAMAPI\n");

    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // For outputting single-read alignment in SAM format
    // Output the best record, and for the rest, append the results in the tag XA:Z
    // numResult should be greater than 0
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( qSetting->occ->SAMOutBuffer );
    HSP * hsp = aIndex->hsp;
    HSPAux * hspaux = aIndex->hspaux;
    long long k = 0;
    unsigned int bestAmbPos, currAmbPos;
    unsigned int bestScore, currScore, secBestScore; // higher score is better
    unsigned long long bestTP, currTP;
    unsigned short bestChr, currChr;
    char bestOccStr[500], currOccStr[500];
    char bestStrand = QUERY_POS_STRAND;
    char currStrand;
    int bestStrLen, currStrLen;
    int isUnaligned = ( algnResult[startIndex].algnmt == 0xFFFFFFFF );
    char * bestSpCigar = NULL;
    char * currSpCigar = NULL;
    int bestEditDist = 0;
    int currEditDist = 0;
    int bestHitNum = 1;
    int secBestHitNum = 0;
    char mdStr[MAX_READ_LENGTH * 2 + 1];
    int mdStrLen = 0;
    char cigarStr[MAX_READ_LENGTH + 1];
    int cigarStrLen = 0;
    int bestMismatchNum = 0;
    int bestGapOpen = 0;
    int bestGapExtend = 0;
    int avg_mismatch_qual;
    int map_qual_score = 0;
    int x1_t1 = 0; // # of suboptimals whose score >= 0.7 * bestDPScore
    int x1_t2 = 0; // # of suboptimals whose score < 0.7 * bestDPScore

    int readlen = qInfo->ReadLength;
    int boundTrim = 0;

    if ( !isUnaligned )
    {

        SingleAlgnmtResult * bestResult = & ( algnResult[startIndex] );
        bestScore = algnResult[startIndex].score;
        secBestScore = 0;

        // obtain the bestHitNum, bestScore, secBestHitNum and secBestScore for the read
        for ( k = 1; k < numResult; k++ )
        {
            currScore = algnResult[startIndex + k].score;
            currStrand = algnResult[startIndex + k].strand;
            currSpCigar = algnResult[startIndex + k].cigarString;
            currEditDist = algnResult[startIndex + k].editdist;
            currAmbPos = algnResult[startIndex + k].algnmt;

            if ( currScore >= bestScore )
            {
                if ( currScore == bestScore )
                {
                    bestHitNum++;
                }
                else
                {
                    secBestScore = bestScore;
                    // secBestHitNum = bestHitNum;
                    bestScore = currScore;
                    bestHitNum = 1;
                    bestResult = & ( algnResult[startIndex + k] );
                }
            }
            else if ( currScore >= secBestScore )
            {
                if ( currScore == secBestScore )
                {
                    // secBestHitNum++;
                }
                else
                {
                    secBestScore = currScore;
                    // secBestHitNum = 1;
                }
            }
        }

        // HANDLE CROSS-CHROMO BUG.
        boundTrim = getChrAndPosWithBoundaryCheckDP ( qInput, readlen, bestResult->algnmt, bestResult->cigarString, &bestTP, &bestChr, NULL ); // check good or not

        if ( boundTrim )
        {
            SingleAlgnmtResult * storedBestResult = bestResult;
            // re-find best alignment
            bestScore = -9998; // hardcode small number
            secBestScore = -9999; // hardcode small number

            for ( k = 0; k < numResult; k++ )
            {
                currAmbPos = algnResult[startIndex + k].algnmt;
                boundTrim = getChrAndPosWithBoundaryCheckDP ( qInput, readlen, currAmbPos, algnResult[startIndex + k].cigarString, &bestTP, &bestChr, NULL );

                if ( !boundTrim ) // good alignment
                {
                    currScore = algnResult[startIndex + k].score;

                    if ( currScore >= bestScore )
                    {
                        if ( currScore == bestScore ) { bestHitNum++; }
                        else { secBestScore = bestScore; bestScore = currScore; bestHitNum = 1; bestResult = & ( algnResult[startIndex + k] ); }
                    }
                    else if ( currScore > secBestScore ) { secBestScore = currScore; }
                }
            }

            if ( ! storedBestResult )
            {
                bestResult = storedBestResult;
                bestScore = bestResult->score;
                secBestScore = 0;
            }

            boundTrim = getChrAndPosWithBoundaryCheckDP ( qInput, readlen, bestResult->algnmt, bestResult->cigarString, &bestTP, &bestChr, NULL ); // check good or not
        }

        // obtain x1_t1 and x1_t2
        int subopt_class_thres = ( int ) ( 0.7 * bestScore );

        if ( numResult > 1 && ( hspaux->alignmentType == OUTPUT_ALL_VALID || hspaux->alignmentType == OUTPUT_ALL_BEST ) && ! boundTrim )
        {
            // Build the XA:Z tag
            DynamicUint8ArrayReset ( xazArray );


            for ( k = 0; k < numResult; k++ )
            {
                if ( & ( algnResult[startIndex + k] ) == bestResult )
                { continue; }

                if ( ( hspaux->alignmentType == OUTPUT_ALL_BEST ) &&
                        ( algnResult[startIndex + k].score < bestScore ) )
                {
                    continue;
                }

                currScore = algnResult[startIndex + k].score;
                currStrand = algnResult[startIndex + k].strand;
                currSpCigar = algnResult[startIndex + k].cigarString;
                currEditDist = algnResult[startIndex + k].editdist;
                currAmbPos = algnResult[startIndex + k].algnmt;

                if ( getChrAndPosWithBoundaryCheckDP ( qInput, readlen, algnResult[startIndex + k].algnmt, currSpCigar, &currTP, &currChr, NULL ) ) { continue; }

                char * chr_name = samFilePtr->header->target_name[currChr - 1]; // chromosome name
                memcpy ( currOccStr, chr_name, strlen ( chr_name ) );
                int pos = strlen ( chr_name );
                currOccStr[pos++] = ',';
                currOccStr[pos++] = ( currStrand == QUERY_NEG_STRAND ) ? '-' : '+'; // strand
                pos += writeULLToStr ( currTP, & ( currOccStr[pos] ) ); // chromosome position
                currOccStr[pos++] = ',';
                pos += convertToCigarStr ( currSpCigar, & ( currOccStr[pos] ) ); // cigar string
                currOccStr[pos++] = ',';
                pos += writeNumToStr ( currEditDist, & ( currOccStr[pos] ) ); // edit distance
                currOccStr[pos++] = ';';
                currOccStr[pos] = '\0';
                currStrLen = pos;
                appendStringToUint8Array ( xazArray, currOccStr, currStrLen );

                if ( currScore < bestScore )
                {
                    if ( currScore >= subopt_class_thres ) { x1_t1++; }
                    else { x1_t2++; }
                }
            }
        }

        // in updated version: secBestHitNum = x1_t1 + x1_t2
        secBestHitNum = x1_t1 + x1_t2;
        bestStrand = bestResult->strand;
        bestSpCigar = bestResult->cigarString;
        bestEditDist = bestResult->editdist;
        bestAmbPos = bestResult->algnmt;

        // to convert the the best special_cigar into normal cigar string
        char * newSpCigar = NULL;
        boundTrim = getChrAndPosWithBoundaryCheckDP ( qInput, readlen, bestAmbPos, bestSpCigar, &bestTP, &bestChr, &newSpCigar );
        cigarStrLen = convertToCigarStr ( boundTrim ? newSpCigar : bestSpCigar, cigarStr );
        // to collect the mismatch information including MD string and # of mismatches
        mdStrLen = getMisInfoForDP ( hsp, qInfo->ReadCode, qInfo->ReadQuality, qInfo->ReadLength, bestAmbPos, bestStrand,
                                     boundTrim ? newSpCigar : bestSpCigar, mdStr, &bestMismatchNum, &bestGapOpen, &bestGapExtend, &avg_mismatch_qual, boundTrim );
        bestEditDist = bestGapExtend + bestMismatchNum;

        if ( newSpCigar ) { free ( newSpCigar ); }

        // getMapQualScoreForSingleDP(int maxDPScore, int avgMismatchQual, int x0, int x1_t1, int x1_t2, int bestDPScore, int secondBestDPScore, int maxMAPQ, int minMAPQ, int dpThres)
        map_qual_score = getMapQualScoreForSingleDP ( qInfo->ReadLength * hspaux->dpMatchScore, ( hspaux->isFastq == 1 ) ? avg_mismatch_qual : 20, bestHitNum, x1_t1, x1_t2, bestScore, secBestScore, hspaux->maxMAPQ, hspaux->minMAPQ, hspaux->singleDPcutoffThreshold, hspaux->bwaLikeScore, hspaux->g_log_n );
        // old version
        // map_qual_score = getMapQualScoreForSingleDP(bestResult->score, qInfo->ReadLength*hspaux->dpMatchScore, (hspaux->isFastq==1)?avg_mismatch_qual:20, bestHitNum, secBestHitNum, bestScore, secBestScore, hspaux->maxMAPQ, hspaux->minMAPQ);
    }

    if ( numResult > 1 && ( hspaux->alignmentType == OUTPUT_ALL_VALID || hspaux->alignmentType == OUTPUT_ALL_BEST ) )
    {
        initializeSAMAlgnmt2 ( samAlgnmt, qInfo->ReadLength, qInfo->ReadName, qInfo->ReadCode,
                               qInfo->ReadQuality, bestStrand, xazArray->charStr, xazArray->length, cigarStr, isUnaligned,
                               bestMismatchNum, bestEditDist, bestHitNum, secBestHitNum, bestGapOpen, bestGapExtend, mdStr, mdStrLen, map_qual_score, hspaux->readGroup, hspaux->isPrintMDNM );
    }
    else
    {
        initializeSAMAlgnmt2 ( samAlgnmt, qInfo->ReadLength, qInfo->ReadName, qInfo->ReadCode,
                               qInfo->ReadQuality, bestStrand, NULL, 0, cigarStr, isUnaligned,
                               bestMismatchNum, bestEditDist, bestHitNum, secBestHitNum, bestGapOpen, bestGapExtend, mdStr, mdStrLen, map_qual_score, hspaux->readGroup, hspaux->isPrintMDNM );
    }

    //----------------------------------------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------------------------------------

    if ( !isUnaligned )
    {
        // ATTENTION: THIS PART IS NOT 64-BIT COMPATIBLE.
        unsigned int samFlag = 0;

        if ( bestStrand == QUERY_NEG_STRAND )
        {
            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
        }

        // ---------------------->  |
        samAlgnmt->core.tid = bestChr - 1;                          //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = bestTP - 1;                           //SAM: 0-based leftmost coordinate
        samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
        // ---------------------->  |
    }
    else
    {
        unsigned int samFlag = 0;
        samFlag |= SAM_FLAG_READ_UNMAPPED;
        // ---------------------->  |
        samAlgnmt->core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
        samAlgnmt->core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
        samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
        samAlgnmt->core.mtid = -1;
        samAlgnmt->core.mpos = -1;
        samAlgnmt->core.isize = 0;
        samwrite ( samFilePtr, samAlgnmt );
        // ---------------------->  |
    }

}



void OCCDirectWritePairUnmapSAMAPI ( SRAQueryInput * qInput, SRAOccurrence * sraOcc, int isFirst )
{
    // cx this function... X_x
    //    fprintf(stderr,"OCCDirectWritePairUnmapSAMAPI\n");
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRASetting * qSetting = qInput->QuerySetting;
    // int OutFileFormat = qSetting->OutFileFormat;
    // FILE * outFilePtr = qSetting->OutFilePtr;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned int * occAmbiguityMap = hsp->ambiguityMap;
    samfile_t * samFilePtr = qSetting->SAMOutFilePtr;
    bam1_t * samAlgnmt = & ( occ->SAMOutBuffer );
    unsigned long long i = 0; // j=0;//, k=0;
    unsigned long long ambPosition;
    unsigned long long tp;
    unsigned int correctPosition;
    unsigned int approxIndex, approxValue;
    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInput
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
    samAlgnmt->core.bin = bam_reg2bin ( 0, 0 );                 //SAM: bin calculated by bam_reg2bin()
    samAlgnmt->core.n_cigar = 1;                                //SAM: number of CIGAR operations
    samAlgnmt->core.l_qseq = qInfo->ReadLength;                //SAM: length of the query sequence (read)
    samAlgnmt->core.l_qname = strlen ( qInfo->ReadName ) + 1;    //SAM: length of the query name
    samAlgnmt->core.mtid = -1;
    samAlgnmt->core.mpos = -1;
    samAlgnmt->core.isize = 0;
    samAlgnmt->l_aux = 0;
    samAlgnmt->data_len = 0;
    samAlgnmt->m_data = SAM_MDATA_SIZE;
    SAMIUint8ConcatString ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadName, strlen ( qInfo->ReadName ) + 1 ); //Name
    SAMIUint8ConcatUint32 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadLength << BAM_CIGAR_SHIFT ); //CIGAR

    for ( i = 0; i < qInfo->ReadLength / 2; i++ )                                                              //Read
    {
        unsigned char biChar = 0;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2]]] << 4;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReadCode[i * 2 + 1]]];
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    if ( qInfo->ReadLength % 2 == 1 )
    {
        uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReadCode[qInfo->ReadLength - 1]]] << 4;
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), biChar );
    }

    for ( i = 0; i < qInfo->ReadLength; i++ )                                                                  //Quality
    {
        SAMIUint8ConcatUint8 ( samAlgnmt->data, & ( samAlgnmt->data_len ), qInfo->ReadQuality[i] );
    }

    unsigned int auxStart = samAlgnmt->data_len;
    //SAMIUint8ConcatString(samAlgnmt->data,&(samAlgnmt->data_len),"XAZTESTING",11);
    samAlgnmt->l_aux += samAlgnmt->data_len - auxStart;
    //----------------------------------------------------------------------------------------------------------------------------
    // ATTENTION: THIS PART IS NOT 64-BIT COMPATIBLE.
    ambPosition = sraOcc->ambPosition;
    tp = ambPosition;
    correctPosition = ambPosition;
    approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
    approxValue = occAmbiguityMap[approxIndex];

    while ( occTranslate[approxValue].startPos > ambPosition )
    {
        approxValue--;
    }

    correctPosition -= occTranslate[approxValue].correction;
    tp = correctPosition;
    unsigned int samFlag = 1; // it is a paired-end read
    samFlag |= SAM_FLAG_MATE_UNMAPPED;
    samFlag |= SAM_FLAG_FIRST_IN_PAIR * isFirst;
    samFlag |= SAM_FLAG_SECOND_IN_PAIR * ( 1 - isFirst );

    if ( sraOcc->strand == QUERY_NEG_STRAND )
    {
        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND * isFirst;
    }

    //Outputing no alignment
    // ----------------------
    samAlgnmt->core.tid = occTranslate[approxValue].chrID - 1;  //SAM: chromosome ID, defined by bam_header_t
    samAlgnmt->core.pos = tp - 1;                               //SAM: 0-based leftmost coordinate4
    samAlgnmt->core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
    samAlgnmt->core.flag = samFlag;                             //SAM: bitwise flag
    samwrite ( samFilePtr, samAlgnmt );
    // ----------------------
}

void OCCFlushCache ( SRAQueryInput * qInput )
{
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    int outputFormat = qSetting->OutFileFormat;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    FILE * outFilePtr = qSetting->OutFilePtr;

    switch ( outputFormat )
    {
        case SRA_OUTPUT_FORMAT_DEFAULT:
            OCCFlushCacheDefault ( occ, hsp, outFilePtr );
            break;

        case SRA_OUTPUT_FORMAT_PLAIN:
            OCCFlushCachePlain ( occ, hsp, outFilePtr, qInput );
            break;

        case SRA_OUTPUT_FORMAT_SAM:
            OCCFlushCacheSAM ( qInput );
            break;

        case SRA_OUTPUT_FORMAT_SAM_API:
            OCCFlushCacheSAMAPI ( qInput );
            break;

        default:
            qSetting->writtenCacheCount++;
            occ->occPositionCacheCount = 0;
            break;
    }
}

void OCCFlushCacheDP ( SRAQueryInput * qInput )
{
    // for outputting the DP result
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    int outputFormat = qSetting->OutFileFormat;
    OCC * occ = qSetting->occ;
    HSP * hsp = aIndex->hsp;
    FILE * outFilePtr = qSetting->OutFilePtr;

    switch ( outputFormat )
    {
        case SRA_OUTPUT_FORMAT_DEFAULT:
            OCCFlushCacheDefault ( occ, hsp, outFilePtr );
            break;

        case SRA_OUTPUT_FORMAT_PLAIN:
            OCCFlushCachePlainDP ( occ, hsp, outFilePtr );
            break;

        case SRA_OUTPUT_FORMAT_SAM:
            OCCFlushCacheSAM ( qInput );
            break;

        case SRA_OUTPUT_FORMAT_SAM_API:
            OCCFlushCacheSAMAPIDP ( qInput );
            break;

        default:
            qSetting->writtenCacheCount++;
            occ->occPositionCacheCount = 0;
            break;
    }
}

void OCCReportDelimitor ( SRAQueryInput * qInput )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    OCC * occ = qSetting->occ;

    if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );}

    occ->occPositionCache[occ->occPositionCacheCount].tp = qInfo->ReadId;
    occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = 0;
    occ->occPositionCache[occ->occPositionCacheCount].ChromId = OCC_CONST_ALIGNMENT_HEADER;
    occ->occPositionCacheCount++;
}

void OCCReportDelimitorDP ( SRAQueryInput * qInput )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    // for outputting the DP result
    SRASetting * qSetting = qInput->QuerySetting;
    OCC * occ = qSetting->occ;

    if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCacheDP ( qInput );}

    occ->occPositionCache[occ->occPositionCacheCount].tp = qInfo->ReadId;
    occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = 0;
    occ->occPositionCache[occ->occPositionCacheCount].ChromId = OCC_CONST_ALIGNMENT_HEADER;
    occ->occPositionCacheCount++;
}


void OCCReportNoAlignment ( SRAQueryInput * qInput )
{
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    OCC * occ = qSetting->occ;

    if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );}

    occ->occPositionCache[occ->occPositionCacheCount].tp = qInfo->ReadId;
    occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = 0;
    occ->occPositionCache[occ->occPositionCacheCount].ChromId = OCC_CONST_NO_ALIGNMENT;
    occ->occPositionCacheCount++;
}

//OCCReportSARange : asumming l<=r
/*
  unsigned long long OCCReportSARange ( SRAQueryInput * qInput,
  unsigned long long l, unsigned long long r,
  int occMismatch )
  {
  SRAQueryInfo * qInfo = qInput->QueryInfo;
  SRASetting * qSetting = qInput->QuerySetting;
  SRAQueryResultCount * rOutput  = qInput->QueryOutput;
  SRAIndex * aIndex = qInput->AlgnmtIndex;
  BWT * bwt = aIndex->bwt;
  // HSP * hsp = aIndex->hsp;
  HOCC * highOcc = aIndex->highOcc;
  OCC * occ = qSetting->occ;
  // FILE * outFilePtr = aIndex->OutFilePtr;
  //Disable this check as the output file pointer does not
  // necessary exist anymore.
  //if (outFilePtr==NULL) return 0;
  #ifdef DEBUG_2BWT_NO_OUTPUT

  if ( l > r ) { return 0; }

  return r - l + 1;
  #endif
  #ifdef DEBUG_2BWT_OUTPUT_TO_SCREEN
  printf ( "[BGS-IO-Read%u] Found SA Range = %u %u\n", qInfo->ReadId, l, r );
  #endif
  unsigned long long j, k;
  unsigned long long saCount = r - l + 1;

  //Retrieve SA values from either SA or HOCC. Depending on the size of the SA Range.
  if ( highOcc != NULL && r - l >= 4 - 1 )
  {
  //By HOCC data structure.
  unsigned long long hOccCachedL, hOccCachedR, tp;
  unsigned long long approxL = HOCCApproximateRank ( highOcc, l + 1 );
  HOCCGetCachedSARange ( highOcc, approxL, &hOccCachedL, &hOccCachedR );

  while ( hOccCachedL > l )
  {
  approxL--;
  HOCCGetCachedSARange ( highOcc, approxL, &hOccCachedL, &hOccCachedR );
  }

  //Case 1: When the entire SA range falls within a cached range of HOOC
  if ( hOccCachedL <= l && hOccCachedR >= r )
  {
  //Get the first index on HOCC
  unsigned long long hOccAnsIndex = highOcc->table[approxL * 2 + 1];
  hOccAnsIndex += ( l - hOccCachedL );
  j = 0;
  k = hOccAnsIndex;

  while ( j < saCount )
  {
  occ->occPositionCache[occ->occPositionCacheCount].tp = highOcc->occ[k];
  occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = qInfo->ReadStrand;
  occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
  occ->occPositionCache[occ->occPositionCacheCount].occMismatch = occMismatch;
  occ->occPositionCacheCount++;
  rOutput->TotalOccurrences++;

  if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );}

  k++;
  j++;
  }

  //Case 2: When the reported SA range enclosed some cached ranges of HOOC
  }
  else if ( hOccCachedR < r )
  {
  unsigned long long hOccAnsIndex = highOcc->table[approxL * 2 + 1];
  j = l;
  k = hOccAnsIndex;

  while ( j <= r )
  {
  if ( j > hOccCachedR && approxL < highOcc->ttlPattern )
  {
  approxL++;
  HOCCGetCachedSARange ( highOcc, approxL, &hOccCachedL, &hOccCachedR );
  hOccAnsIndex = highOcc->table[approxL * 2 + 1];
  k = hOccAnsIndex;
  }

  if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );}

  if ( hOccCachedL <= j && j <= hOccCachedR )
  {
  tp = highOcc->occ[k];
  occ->occPositionCache[occ->occPositionCacheCount].tp = tp;
  occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = qInfo->ReadStrand;
  occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
  occ->occPositionCache[occ->occPositionCacheCount].occMismatch = occMismatch;
  occ->occPositionCacheCount++;
  rOutput->TotalOccurrences++;
  k++;
  }
  else
  {
  tp = BWTSaValue ( bwt, j );
  occ->occPositionCache[occ->occPositionCacheCount].tp = tp;
  occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = qInfo->ReadStrand;
  occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
  occ->occPositionCache[occ->occPositionCacheCount].occMismatch = occMismatch;
  occ->occPositionCacheCount++;
  rOutput->TotalOccurrences++;
  }

  j++;
  }

  //Case 3: Cannot find the records in HOCC (Should not happen)
  //In this case, we compute the occ by Suffix array
  }
  else
  {
  for ( j = l; j <= r; j++ )
  {
  if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );}

  occ->occPositionCache[occ->occPositionCacheCount].tp = BWTSaValue ( bwt, j );
  occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = qInfo->ReadStrand;
  occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
  occ->occPositionCache[occ->occPositionCacheCount].occMismatch = occMismatch;
  occ->occPositionCacheCount++;
  rOutput->TotalOccurrences++;
  }
  }

  }
  else
  {
  //By Suffix array.
  for ( j = l; j <= r; j++ )
  {
  if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );}

  occ->occPositionCache[occ->occPositionCacheCount].tp = BWTSaValue ( bwt, j );
  occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = qInfo->ReadStrand;
  occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
  occ->occPositionCache[occ->occPositionCacheCount].occMismatch = occMismatch;
  occ->occPositionCacheCount++;
  rOutput->TotalOccurrences++;
  }
  }

  return saCount;
  }


  unsigned long long OCCReportTextPositions ( SRAQueryInput * qInput, int posCount )
  {
  //DEBUG : OUTPUT TO SCREEN
  #ifdef DEBUG_2BWT_NO_OUTPUT
  int i;

  for ( i = 0; i < posCount; i++ )
  {
  qInput->QueryOutput->WithError[qInput->QueryOutput->Shared_AccMismatch[i]]++;
  }

  return posCount;
  #endif

  //END OF DEBUG

  if ( posCount <= 0 ) { return 0; }

  SRAQueryInfo * qInfo = qInput->QueryInfo;
  SRASetting * qSetting = qInput->QuerySetting;
  SRAIndex * aIndex = qInput->AlgnmtIndex;
  SRAQueryResultCount * rOutput  = qInput->QueryOutput;
  unsigned long long * saPositions = rOutput->Shared_ReadPositions;
  int * occMismatches = rOutput->Shared_AccMismatch;
  // char * occQualities = rOutput->Shared_AccQualities;
  // int OutputType = qSetting->OutputType;
  // int outputFormat = qSetting->OutFileFormat;
  // FILE * outFilePtr = qSetting->OutFilePtr;
  // HSP * hsp = aIndex->hsp;
  OCC * occ = qSetting->occ;
  int k;

  //
  for ( k = 0; k < posCount; k++ )
  {
  if ( occ->occPositionCacheCount >= OCC_CACHE_SIZE ) {OCCFlushCache ( qInput );}

  occ->occPositionCache[occ->occPositionCacheCount].tp = saPositions[k];
  occ->occPositionCache[occ->occPositionCacheCount].ReadStrand = qInfo->ReadStrand;
  occ->occPositionCache[occ->occPositionCacheCount].ChromId = 0;
  occ->occPositionCache[occ->occPositionCacheCount].occMismatch = occMismatches[k];
  occ->occPositionCacheCount++;
  rOutput->WithError[occMismatches[k]]++;
  rOutput->TotalOccurrences++;
  }

  rOutput->RetrievedByCe += posCount;
  return posCount;
  }

*/
