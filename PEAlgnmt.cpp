/*
 *
 *    PEAlgnmt.c
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

#include "PEAlgnmt.h"

// --------------------------------------------------
// PE_DEBUG #1
// Uncomment the below to turn on logging message for
// PE matching of two occurrences.
//#define PE_DEBUG_PRINT_MATCHING_PROGRESS
// --------------------------------------------------
// --------------------------------------------------
// PE_DEBUG #2
// Uncomment the below to turn on logging message for
// PE Radix Sort of two occ lists.
//#define PE_DEBUG_PRINT_RADIX_PROGRESS
// --------------------------------------------------
// --------------------------------------------------
// PE_DEBUG #3
// Uncomment the below to turn on visual logging for
// PE matching of two occurrences.
//#define PE_DEBUG_PRINT_VISUAL_MATCHING
// --------------------------------------------------

void PERadixSort ( SRAOccurrence * list, unsigned int count,
                   SRAOccurrence * buffer1, SRAOccurrence * buffer2,
                   SRAOccurrence * sortedList );

void PEReportPairResult ( PEOutput * peOutput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int insertion );

inline int PEIsPairEndMatch ( PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int * insertion );
inline int PEIsPairOutOfRange ( PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2 );

void PEPairListInitialise ( PEPairList * pePairList );
void PEOutputInitialise ( PEOutput * peOutput );
PEPairList * PEPairListConstruct ();
void PEPairListFree ( PEPairList * pePairList );

/////////////////////////////////////////////////////
/*
    Function SRAEnrichSARanges
*/
/////////////////////////////////////////////////////
/*
unsigned int SRAEnrichSARanges(BWT * bwt, HSP * hsp,
                                unsigned int saIndexLeft, unsigned int saIndexRight,
                                char * strand, int * mismatchCount) {

}
*/
/////////////////////////////////////////////////////
/*
    Function PERetrievePositionFromSAIndex
    This function is a place-holder for possibly extension of
    PE to use High-Occ data-structure. However unlikely.
*/
/////////////////////////////////////////////////////
inline unsigned int PERetrievePositionFromSAIndex ( BWT * bwt, unsigned int saIndex )
{
    unsigned int ambPosition = ( *bwt->_bwtSaValue ) ( bwt, saIndex );
    return ambPosition;
}

/////////////////////////////////////////////////////
/*
    Function PERetrieveChromoPositioning
*/
/////////////////////////////////////////////////////
void PERetrieveChromoPositioning ( BWT * bwt,
                                   HSP * hsp,
                                   unsigned int ambPosition,
                                   int * sequenceId, unsigned int * offset )
{
    unsigned int * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    unsigned int approxIndex = ambPosition >> GRID_SAMPLING_FACTOR_2_POWER;
    unsigned int approxValue = ambiguityMap[approxIndex];

    while ( translate[approxValue].startPos > ambPosition )
    {
        approxValue--;
    }

    ambPosition -= translate[approxValue].correction;
    ( *sequenceId ) = translate[approxValue].chrID;
    ( *offset ) = ambPosition;
}

/////////////////////////////////////////////////////
/*
    Function PERadixSort
*/
/////////////////////////////////////////////////////
void PERadixSort ( SRAOccurrence * list, unsigned int count,
                   SRAOccurrence * buffer1, SRAOccurrence * buffer2,
                   SRAOccurrence * sortedList )
{
    //How many bits we sort each pass
#define PE_RADIX_BIT_PER_PASS 4
    //How many pass in total we run the sorting
#define PE_RADIX_PASS         8
    //The multiplying result of the two parameters above should be
    //greater than the bit-length of the maximum index on the ref.seq.
    //e.g. 4 and 8 such that the sorting will cover 32-bit
    //     which can represent all index in the ref. seq.
    unsigned int i, j;
    unsigned int base;
    int idxIndexSize = ( 1 << PE_RADIX_BIT_PER_PASS );
    int radixMask = idxIndexSize - 1;
    unsigned int idxRadix[idxIndexSize + 1];
    SRAOccurrence * swap;
    SRAOccurrence * activeBuffer = buffer1;
    SRAOccurrence * oldBuffer = list;
#ifdef PE_DEBUG_PRINT_RADIX_PROGRESS
    printf ( "  idxIndexSize = %d\n", idxIndexSize );
    printf ( "  radixMask = %d\n", radixMask );
    printf ( "  PE_RADIX_BIT_PER_PASS = %d\n", PE_RADIX_BIT_PER_PASS );
    printf ( "  PE_RADIX_PASS = %d\n", PE_RADIX_PASS );
#endif

    for ( base = 0; base < PE_RADIX_PASS; base++ )
    {
        //Initialise idxRadix
        for ( i = 0; i <= idxIndexSize; i++ )
        {
            idxRadix[i] = 0;
        }

        //Count oldBuffer into activeBuffer
        for ( i = 0; i < count; i++ )
        {
            unsigned int key = oldBuffer[i].ambPosition;
            key = key >> base * PE_RADIX_BIT_PER_PASS;
            key &= radixMask;
            idxRadix[key + 1]++;
        }

        //Post-process idxRadix
        for ( i = 1; i <= idxIndexSize; i++ )
        {
            idxRadix[i] += idxRadix[i - 1];
        }

        //Move oldBuffer into activeBuffer with ordering
        for ( i = 0; i < count; i++ )
        {
            unsigned int key = oldBuffer[i].ambPosition;
            key = key >> base * PE_RADIX_BIT_PER_PASS;
            key &= radixMask;
            activeBuffer[idxRadix[key]].ambPosition = oldBuffer[i].ambPosition;
            activeBuffer[idxRadix[key]].mismatchCount = oldBuffer[i].mismatchCount;
            activeBuffer[idxRadix[key]].strand = oldBuffer[i].strand;
            idxRadix[key]++;
        }

        //Swap double buffer
        swap = activeBuffer;

        if ( base == 0 )
        {
            activeBuffer = buffer2;
        }
        else
        {
            activeBuffer = oldBuffer;
        }

        oldBuffer = swap;
    }

    //Count oldBuffer into sortedlist
    for ( i = 0; i < count; i++ )
    {
        sortedList[i].ambPosition = oldBuffer[i].ambPosition;
        sortedList[i].mismatchCount = oldBuffer[i].mismatchCount;
        sortedList[i].strand = oldBuffer[i].strand;
    }
}


inline void PEMappingCore ( PEInput * peInput, PEOutput * peOutput,
                            SRAOccurrence * occList_1, unsigned int occCount_1,
                            SRAOccurrence * occList_2, unsigned int occCount_2 )
{
    // ATTENTION : BRUTE FORCE ENGINE
    //---------------------------------------------------------
    /*unsigned int occIndex_First = 0;
    unsigned int occIndex_Last = 0;
    for (occIndex_1=0;occIndex_1<occCount_1;occIndex_1++) {
        int flag = 0;

        for (occIndex_2=occIndex_First;occIndex_2<occCount_2;occIndex_2++) {
            if (PEIsPairEndMatch(peInput, &(occList_1[occIndex_1]), &(occList_2[occIndex_2]))) {
                if (flag==0) {
                    occIndex_First = occIndex_2;
                    flag=1;
                }
                //Output answer
                PEReportPairResult(peOutput,&(occList_1[occIndex_1]),&(occList_2[occIndex_2]));
                occIndex_Last = occIndex_2;
            }
            if (PEIsPairOutOfRange(peInput, &(occList_1[occIndex_1]), &(occList_2[occIndex_2]))) {
                break;
            }
        }
    }*/
    // ATTENTION : BRUTE FORCE ENGINE - MARK 2
    //---------------------------------------------------------
    unsigned long i;
    unsigned long occIndex_1 = 0;
    unsigned long occIndex_2 = 0;
    int OutputType = peInput->OutputType;
    int strandLeftLeg = peInput->strandLeftLeg;
    int strandRightLeg = peInput->strandRightLeg;
    unsigned int insertion;

    while ( occIndex_1 < occCount_1 && occIndex_2 < occCount_2 )
    {
        if ( occList_1[occIndex_1].ambPosition <= occList_2[occIndex_2].ambPosition )
        {
            if ( occList_1[occIndex_1].strand == strandLeftLeg )
            {
                for ( i = occIndex_2; i < occCount_2; i++ )
                {
                    if ( PEIsPairEndMatch ( peInput, & ( occList_1[occIndex_1] ), & ( occList_2[i] ), &insertion ) )
                    {
                        PEReportPairResult ( peOutput, & ( occList_1[occIndex_1] ), & ( occList_2[i] ), insertion );

                        if ( OutputType == PE_REPORT_ONE ) { break; }
                    }

                    if ( PEIsPairOutOfRange ( peInput, & ( occList_1[occIndex_1] ), & ( occList_2[i] ) ) )
                    {
                        break;
                    }
                }
            }

            occIndex_1++;
        }
        else
        {
#ifndef BGS_FWD_FOR_FIRST_REV_FOR_SECOND

            if ( occList_2[occIndex_2].strand == strandLeftLeg )
            {
                for ( i = occIndex_1; i < occCount_1; i++ )
                {
                    if ( PEIsPairEndMatch ( peInput, & ( occList_2[occIndex_2] ), & ( occList_1[i] ), &insertion ) )
                    {
                        PEReportPairResult ( peOutput, & ( occList_1[i] ), & ( occList_2[occIndex_2] ), insertion );

                        if ( OutputType == PE_REPORT_ONE ) { break; }
                    }

                    if ( PEIsPairOutOfRange ( peInput, & ( occList_2[occIndex_2] ), & ( occList_1[i] ) ) )
                    {
                        break;
                    }
                }
            }

#endif
            occIndex_2++;
        }
    }

    //---------------------------------------------------------
}


/////////////////////////////////////////////////////
/*
    Function PEValidateAndPreparePEInput
*/
/////////////////////////////////////////////////////
inline int PEValidateAndPreparePEInput ( PEInput * peInput )
{
    int isValid = 1;
    int isWarning = 0;

    if ( peInput->bwt == NULL || peInput->hsp == NULL )
    {
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf ( "[PEValidateAndPreparePEInput] BWT/HSP index missing.\n" );
#endif
        isValid = 0;
    }

    if ( peInput->OutputType != PE_REPORT_ALL &&
            peInput->OutputType != PE_REPORT_ONE )
    {
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf ( "[PEValidateAndPreparePEInput] Defaulting Output Type.\n" );
#endif
        peInput->OutputType = PE_REPORT_ALL;
        isWarning = 1;
    }

    if ( peInput->strandLeftLeg == PE_NON_INPUT_VALUE ||
            peInput->strandRightLeg == PE_NON_INPUT_VALUE )
    {
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf ( "[PEValidateAndPreparePEInput] Defaulting strandLeftLeg / strandRightLeg.\n" );
#endif
        peInput->strandLeftLeg = QUERY_POS_STRAND;
        peInput->strandRightLeg = QUERY_NEG_STRAND;
        isWarning = 1;
    }

    if ( peInput->patternLength == PE_NON_INPUT_VALUE )
    {
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf ( "[PEValidateAndPreparePEInput] patternLength missing.\n" );
#endif
        isValid = 0;
    }

    if ( ( peInput->insertSizeMean == PE_NON_INPUT_VALUE ||
            peInput->insertSizeStdDev == PE_NON_INPUT_VALUE ) &&
            ( peInput->insertLbound == PE_NON_INPUT_VALUE ||
              peInput->insertUbound == PE_NON_INPUT_VALUE ) )
    {
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf ( "[PEValidateAndPreparePEInput] insertion size parameter missing.\n" );
#endif
        isValid = 0;
    }

    if ( peInput->insertLbound == PE_NON_INPUT_VALUE ||
            peInput->insertUbound == PE_NON_INPUT_VALUE )
    {
        peInput->insertLbound = peInput->insertSizeMean - 3 * peInput->insertSizeStdDev;
        peInput->insertUbound = peInput->insertSizeMean + 3 * peInput->insertSizeStdDev;
    }

    return isValid;
}


/////////////////////////////////////////////////////
/*
    Function PEMapping
*/
/////////////////////////////////////////////////////
void PEMapping ( PEInput * peInput, PEOutput * peOutput,
                 PESRAAlignmentResult * resultListA, unsigned int resultCountA,
                 PESRAAlignmentResult * resultListB, unsigned int resultCountB )
{
    //2BWT Index parameter
    BWT * bwt = peInput->bwt;
    HSP * hsp = peInput->hsp;
    //Declare variables
    unsigned int i, j;
    int k;
    unsigned int allocMemory = 0;
    //These sort buffer will be used as double buffer for the sorting
    SRAOccurrence * sortBuffer_1;
    SRAOccurrence * sortBuffer_2;
    SRAOccurrence * occList_1;
    SRAOccurrence * occList_2;
    unsigned int occIndex_1, occIndex_2;
    //Initialise the peOutput
    PEOutputInitialise ( peOutput );

    //enrich the peInput
    // [m - 3 * sd, m + 3 * sd], and the
    if ( !PEValidateAndPreparePEInput ( peInput ) )
    {
        peOutput->flag = PE_ALIGNMENT_INPUT_ERROR;
        return;
    }

    //More information about the occurrence list
    unsigned int occCount_1 = 0, occCount_2 = 0;
    unsigned int occMaxCount = 0;

    for ( i = 0; i < resultCountA; i++ )
    {
        occCount_1 += resultListA[i].saIndexRight - resultListA[i].saIndexLeft + 1;
    }

    for ( i = 0; i < resultCountB; i++ )
    {
        occCount_2 += resultListB[i].saIndexRight - resultListB[i].saIndexLeft + 1;
    }

    occMaxCount = occCount_1;

    if ( occMaxCount < occCount_2 )
    {
        occMaxCount = occCount_2;
    }

    //Allocate enough memory for all occurrences
    //We need at most (4 x max(occCount_1,occCount_2) x 10 Bytes (SRAOccurrence) in total for PE mapping.
    //Assuming we have at least 256-Megabyte left in main memory
#define PE_ALLOCATED_BYTE 268435456
    allocMemory = 4 * occMaxCount * sizeof ( SRAOccurrence );

    if ( allocMemory > PE_ALLOCATED_BYTE )
    {
        fprintf ( stderr, "Insufficient memory to proceed. %u required exceeded %u.\n", allocMemory, PE_ALLOCATED_BYTE );
        exit ( 1 );
    }

    sortBuffer_1 = ( SRAOccurrence * ) malloc ( occMaxCount * sizeof ( SRAOccurrence ) ) ;
    sortBuffer_2 = ( SRAOccurrence * ) malloc ( occMaxCount * sizeof ( SRAOccurrence ) ) ;
    occList_1 = ( SRAOccurrence * ) malloc ( occCount_1 * sizeof ( SRAOccurrence ) ) ;
    occList_2 = ( SRAOccurrence * ) malloc ( occCount_2 * sizeof ( SRAOccurrence ) ) ;
    //Retrieve the occurrences
    occIndex_1 = 0;

    for ( i = 0; i < resultCountA; i++ )
    {
        for ( j = resultListA[i].saIndexLeft; j <= resultListA[i].saIndexRight; j++ )
        {
            occList_1[occIndex_1].ambPosition = PERetrievePositionFromSAIndex ( bwt, j );
            occList_1[occIndex_1].strand = resultListA[i].strand;
            occList_1[occIndex_1].mismatchCount = resultListA[i].mismatchCount;
            occIndex_1++;
        }
    }

    occIndex_2 = 0;

    for ( i = 0; i < resultCountB; i++ )
    {
        for ( j = resultListB[i].saIndexLeft; j <= resultListB[i].saIndexRight; j++ )
        {
            occList_2[occIndex_2].ambPosition = PERetrievePositionFromSAIndex ( bwt, j );
            occList_2[occIndex_2].strand = resultListB[i].strand;
            occList_2[occIndex_2].mismatchCount = resultListB[i].mismatchCount;
            occIndex_2++;
        }
    }

#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
    printf ( "%u(%u) positions are retrieved from %u SA ranges.\n", occIndex_1, occCount_1, resultCountA );
    printf ( "%u(%u) positions are retrieved from %u SA ranges.\n", occIndex_2, occCount_2, resultCountB );
#endif
    //Radix sort the lists
    PERadixSort ( occList_1, occCount_1, sortBuffer_1, sortBuffer_2, occList_1 );
    PERadixSort ( occList_2, occCount_2, sortBuffer_1, sortBuffer_2, occList_2 );
    //Calling the core
    PEMappingCore ( peInput, peOutput, occList_1, occCount_1, occList_2, occCount_2 );
    peOutput->flag = PE_ALIGNMENT_COMPLETED;
    free ( sortBuffer_1 );
    free ( sortBuffer_2 );
    free ( occList_1 );
    free ( occList_2 );
}


/////////////////////////////////////////////////////
/*
    Function PEMappingOccurrences
*/
/////////////////////////////////////////////////////
void PEMappingOccurrences ( PEInput * peInput, PEOutput * peOutput,
                            SRAOccurrence * occList_1, unsigned int occCount_1,
                            SRAOccurrence * occList_2, unsigned int occCount_2 )
{
    //2BWT Index parameter
    BWT * bwt = peInput->bwt;
    HSP * hsp = peInput->hsp;
    //Declare variables
    unsigned int i, j;
    int k;
    unsigned int allocMemory = 0;
    //These sort buffer will be used as double buffer for the sorting
    SRAOccurrence * sortBuffer_1;
    SRAOccurrence * sortBuffer_2;
    unsigned int occIndex_1, occIndex_2;
    //Initialise the peOutput
    PEOutputInitialise ( peOutput );

    //enrich the peInput
    // [m - 3 * sd, m + 3 * sd], and the
    if ( !PEValidateAndPreparePEInput ( peInput ) )
    {
        peOutput->flag = PE_ALIGNMENT_INPUT_ERROR;
        return;
    }

    //More information about the occurrence list
    unsigned int occMaxCount = 0;
    occMaxCount = occCount_1;

    if ( occMaxCount < occCount_2 )
    {
        occMaxCount = occCount_2;
    }

    //Allocate enough memory for all occurrences
    //We need at most (4 x max(occCount_1,occCount_2) x 10 Bytes (SRAOccurrence) in total for PE mapping.
    //Assuming we have at least 256-Megabyte left in main memory
#define PE_ALLOCATED_BYTE 268435456
    allocMemory = 4 * occMaxCount * sizeof ( SRAOccurrence );

    if ( allocMemory > PE_ALLOCATED_BYTE )
    {
        fprintf ( stderr, "Insufficient memory to proceed. %u required exceeded %u.\n", allocMemory, PE_ALLOCATED_BYTE );
        exit ( 1 );
    }

    sortBuffer_1 = ( SRAOccurrence * ) malloc ( occMaxCount * sizeof ( SRAOccurrence ) ) ;
    sortBuffer_2 = ( SRAOccurrence * ) malloc ( occMaxCount * sizeof ( SRAOccurrence ) ) ;
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
    printf ( "%u positions in list 1.\n", occCount_1 );
    printf ( "%u positions in list 2.\n", occCount_2 );
#endif

    //Radix sort the lists
    if ( occCount_1 > 1 )
    { PERadixSort ( occList_1, occCount_1, sortBuffer_1, sortBuffer_2, occList_1 ); }

    if ( occCount_2 > 1 )
    { PERadixSort ( occList_2, occCount_2, sortBuffer_1, sortBuffer_2, occList_2 ); }

    //Calling the core
    PEMappingCore ( peInput, peOutput, occList_1, occCount_1, occList_2, occCount_2 );
    peOutput->flag = PE_ALIGNMENT_COMPLETED;
    free ( sortBuffer_1 );
    free ( sortBuffer_2 );
}

/////////////////////////////////////////////////////
/*
    Function PEIsPairEndMatch
    This function takes in 2 occurrence and determind
    if the two occurrences form a pair end mapping.
    Assume occ_1 <= occ_2.

    Input:
        peInput - The input parameter for PE
        occ_1 - The occurrence on left leg
        occ_2 - The occurrence on right leg

    Output:
        1 if the two occurrences form a pair end mapping
        0 otherwise

*/
/////////////////////////////////////////////////////
inline int PEIsPairEndMatch ( PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int * insertion )
{
    //For a valid paired-end reads,
    // let the mean and the standard deviation of the insert size be "m" and "sd",
    // the distance between the leftmost of the left-read alignment and the rightmost
    // of the right-read alignment should fall between [m - 3 * sd, m + 3 * sd], and the
    // strand of the left-read alignment is NOT the same as the strand of the right-read alignment.
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
    printf ( "PEIsPairEndMatch is invoked.\n" );
    printf ( "    pos = %u \t strand = %d \n", occ_1->ambPosition, occ_1->strand );
    printf ( "    pos = %u \t strand = %d \n", occ_2->ambPosition, occ_2->strand );
#endif
    int strandLeftLeg = peInput->strandLeftLeg;
    int strandRightLeg = peInput->strandRightLeg;
    int strandCheck = ( occ_1->strand == strandLeftLeg &&
                        occ_2->strand == strandRightLeg );
    //positions are left aligned, therefore it requires correction
    unsigned int ambPos_2 = occ_2->ambPosition + peInput->patternLength - 1;
    //occ_2->ambPosition += peInput->patternLength - 1;
    unsigned int gap = ambPos_2 - occ_1->ambPosition + 1;
    ( *insertion ) = gap;
#ifdef PE_DEBUG_PRINT_VISUAL_MATCHING
    char DEBUG_sorted[2];
    DEBUG_sorted[0] = 'X';
    DEBUG_sorted[1] = ' ';
    printf ( " %10llu <----- %10llu -----> %10llu    [%d:%d] %c\n", occ_1->ambPosition, gap, ambPos_2, peInput->insertLbound, peInput->insertUbound, DEBUG_sorted[ ( occ_1->ambPosition <= ambPos_2 ) & 1] );
#endif
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
    printf ( "    lbound = %u \t gap = %u \t ubound = %u\n", peInput->insertLbound, gap, peInput->insertUbound );
#endif
    return ( ( peInput->insertLbound <= gap ) && ( gap <= peInput->insertUbound ) && strandCheck );
}

inline int PEIsPairOutOfRange ( PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2 )
{
    if ( occ_1->strand == occ_2->strand ) { return 0; }

    unsigned int occ_1_pos = occ_1->ambPosition + peInput->insertUbound;
    unsigned int occ_2_pos = occ_2->ambPosition + peInput->patternLength - 1;
    return occ_1_pos < occ_2_pos;
}

void PEReportPairResult ( PEOutput * peOutput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int insertion )
{
#ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
    printf ( "    Pair-end alignment is found.\n" );
#endif
    PEPairList * pePairList = peOutput->tail;

    if ( pePairList->pairsCount == PE_MAX_BUCKET_SIZE )
    {
        if ( pePairList->next == NULL )
        {
            pePairList->next = PEPairListConstruct ();
        }

        pePairList = pePairList->next;
        peOutput->tail = pePairList;
    }

    PEPairs * pePair = & ( pePairList->pairs[pePairList->pairsCount] );
    pePair->algnmt_1 = occ_1->ambPosition;
    pePair->strand_1 = occ_1->strand;
    pePair->mismatch_1 = occ_1->mismatchCount;
    pePair->algnmt_2 = occ_2->ambPosition;
    pePair->strand_2 = occ_2->strand;
    pePair->mismatch_2 = occ_2->mismatchCount;
    pePair->insertion = insertion;
    pePair->totalMismatchCount = occ_1->mismatchCount + occ_2->mismatchCount;
    pePairList->pairsCount++;
}

/*
unsigned int PEMappingToDisk(PEInput * peInput, PEPairList * pePairList) {
}
*/
/////////////////////////////////////////////////////
// Constructor and Destructor
/////////////////////////////////////////////////////
PEInput * PEInputConstruct ( BWT * bwt, HSP * hsp )
{
    PEInput * peInput = ( PEInput * ) malloc ( sizeof ( PEInput ) );
    peInput->bwt = bwt;
    peInput->hsp = hsp;
    peInput->patternLength = PE_NON_INPUT_VALUE;
    peInput->strandLeftLeg = PE_NON_INPUT_VALUE;
    peInput->strandRightLeg = PE_NON_INPUT_VALUE;
    peInput->insertSizeMean = PE_NON_INPUT_VALUE;
    peInput->insertSizeStdDev = PE_NON_INPUT_VALUE;
    peInput->insertLbound = PE_NON_INPUT_VALUE;
    peInput->insertUbound = PE_NON_INPUT_VALUE;
    peInput->OutputType = PE_REPORT_ALL;
    return peInput;
}

void PEInputFree ( PEInput * peInput )
{
    free ( peInput );
}

void PEPairListInitialise ( PEPairList * pePairList )
{
    if ( pePairList == NULL ) { return; }

    PEPairListInitialise ( pePairList->next );
    pePairList->pairsCount = 0;
}

void PEOutputInitialise ( PEOutput * peOutput )
{
    PEPairListInitialise ( peOutput->root );
    peOutput->tail = peOutput->root;
    peOutput->flag = PE_ALIGNMENT_INITIALISED;
}

PEPairList * PEPairListConstruct ()
{
    PEPairList * pePairList = ( PEPairList * ) malloc ( sizeof ( PEPairList ) );
    pePairList->pairsCount = 0;
    pePairList->next = NULL;
    return pePairList;
}

PEOutput * PEOutputConstruct ()
{
    PEOutput * peOutput = ( PEOutput * ) malloc ( sizeof ( PEOutput ) );
    peOutput->root = PEPairListConstruct ();
    peOutput->tail = peOutput->root;
    return peOutput;
}

void PEPairListFree ( PEPairList * pePairList )
{
    if ( pePairList == NULL ) { return; }

    PEPairListFree ( pePairList->next );
    free ( pePairList );
}

void PEOutputFree ( PEOutput * peOutput )
{
    PEPairListFree ( peOutput->root );
    free ( peOutput );
}

/////////////////////////////////////////////////////
// Utilities
/////////////////////////////////////////////////////
void PEPrintOccurrenceList ( const char * listName, SRAOccurrence * list, unsigned int length )
{
    unsigned int i;
    printf ( "List of occurrences %s =\n", listName );

    for ( i = 0; i < length; i++ )
    {
        printf ( "%u ", list[i].ambPosition );

        if ( ( i + 1 ) % 8 == 0 ) {printf ( "\n" );}
    }

    printf ( "\n" );
}
void PEPrintPEPairList ( PEPairList * pairList, int maxOutput )
{
    int i;
    int j = 0;
    char strand[3];
    strand[1] = '+';
    strand[2] = '-';

    while ( pairList != NULL && pairList->pairsCount > 0 )
    {
        for ( i = 0; i < pairList->pairsCount; i++ )
        {
            PEPairs * pePair = & ( pairList->pairs[i] );
            printf ( "algnmt_1 = %u(%c)      algnmt_2 = %u(%c)\n",
                     pePair->algnmt_1, strand[pePair->strand_1],
                     pePair->algnmt_2, strand[pePair->strand_2] );
            j++;

            if ( j >= maxOutput && maxOutput >= 0 ) { return; }
        }

        pairList = pairList->next;
    }
}
void PEPrintPEOutput ( PEOutput * peOutput )
{
    PEPrintPEPairList ( peOutput->root, -1 );
}
void PEPrintPEOutputWithLimit ( PEOutput * peOutput, int maxOutput )
{
    PEPrintPEPairList ( peOutput->root, maxOutput );
}
unsigned int PECountPEPairList ( PEPairList * pairList )
{
    unsigned int count = 0;

    while ( pairList != NULL && pairList->pairsCount > 0 )
    {
        count += pairList->pairsCount;
        pairList = pairList->next;
    }

    return count;
}
unsigned int PECountPEOutput ( PEOutput * peOutput )
{
    return PECountPEPairList ( peOutput->root );
}

unsigned int PEStatsPEPairList ( PEPairList * pairList, PEPairs ** optimal, PEPairs ** suboptimal, unsigned int * mismatchStats )
{
    unsigned int count = 0;
    unsigned int i;

    PEPairs * iOptimal = NULL;
    PEPairs * iSubOptimal = NULL;
    uint8_t iOptimal_MismatchCount = 255;
    uint8_t iOptimal_MismatchDiff = 255;
    uint8_t iSubOptimal_MismatchCount = 255;
    uint8_t iSubOptimal_MismatchDiff = 255;

    while ( pairList != NULL && pairList->pairsCount > 0 )
    {
        count += pairList->pairsCount;

        for ( i = 0; i < pairList->pairsCount; i++ )
        {
            PEPairs * pePair = & ( pairList->pairs[i] );
            int numMismatchOnPE = pePair->totalMismatchCount;
            mismatchStats[numMismatchOnPE]++;

            int diffMismatch = pePair->mismatch_1 - pePair->mismatch_2;

            if ( pePair->mismatch_2 > pePair->mismatch_1 )
            {
                diffMismatch = -diffMismatch;
            }

            if ( numMismatchOnPE < iOptimal_MismatchCount )
            {
                iSubOptimal_MismatchCount = iOptimal_MismatchCount;
                iSubOptimal_MismatchDiff = iOptimal_MismatchDiff;
                iSubOptimal = iOptimal;
                iOptimal = pePair;
                iOptimal_MismatchCount = numMismatchOnPE;
                iOptimal_MismatchDiff = diffMismatch;
            }
            else if ( numMismatchOnPE == iOptimal_MismatchCount &&
                      diffMismatch < iOptimal_MismatchDiff )
            {
                iOptimal = pePair;
                iOptimal_MismatchCount = numMismatchOnPE;
                iOptimal_MismatchDiff = diffMismatch;
            }

        }

        pairList = pairList->next;
    }

    ( *optimal ) = iOptimal;
    ( *suboptimal ) = iSubOptimal;

    return count;
}

unsigned int PEStatsPEOutput ( PEOutput * peOutput, PEPairs ** optimal, PEPairs ** suboptimal, unsigned int * mismatchStats )
{
    return PEStatsPEPairList ( peOutput->root, optimal, suboptimal, mismatchStats );
}

// Construct the structure ReadInputForDP for each CPU thread
ReadInputForDP * constructReadInputForDP ( int cpuNum )
{
    ReadInputForDP * readInput;
    readInput = ( ReadInputForDP * ) malloc ( sizeof ( ReadInputForDP ) );
    readInput->readNum = 0;
    unsigned int initial_size_sa_list = INITIAL_SIZE_SA_LIST_FOR_DP / cpuNum;
    unsigned int initial_size_occ_list = INITIAL_SIZE_OCC_LIST_FOR_DP / cpuNum;
    readInput->sa_list = ( PESRAAlignmentResult * ) malloc ( initial_size_sa_list * sizeof ( PESRAAlignmentResult ) );
    readInput->occ_list = ( SRAOccurrence * ) malloc ( initial_size_occ_list * sizeof ( SRAOccurrence ) );
    readInput->saRangeTotalNum = 0;
    readInput->occTotalNum = 0;
    readInput->saListSize = initial_size_sa_list;
    readInput->occListSize = initial_size_occ_list;
    readInput->lastReadID = 0xFFFFFFFF;
    return readInput;
}

// Reset the structure ReadInputForDP
void resetReadInputForDP ( ReadInputForDP * readInput )
{
    readInput->readNum = 0;
    readInput->saRangeTotalNum = 0;
    readInput->occTotalNum = 0;
}

// Release the memory for the structure ReadInputForDP
void freeReadInputForDP ( ReadInputForDP * readInput )
{
    free ( readInput->sa_list );
    free ( readInput->occ_list );
    free ( readInput );
}


// To add the alignment results of a read to ReadInputForDP
void addToReadInputForDP ( ReadInputForDP * readInput, unsigned int readid, PESRAAlignmentResult * saList, unsigned int saNum,
                           SRAOccurrence * occList, unsigned int occNum )
{
    if ( saNum > 0 )
    {
        if ( readInput->saRangeTotalNum + saNum > readInput->saListSize )
        {
            // enlarge the arrays "sa_list" by at least double
            unsigned int new_size = readInput->saListSize * 2;

            while ( readInput->saRangeTotalNum + saNum > new_size )
            {
                new_size *= 2;
            }

            PESRAAlignmentResult * new_sa_list = ( PESRAAlignmentResult * ) malloc ( sizeof ( PESRAAlignmentResult ) * new_size );
            memcpy ( new_sa_list, readInput->sa_list, sizeof ( PESRAAlignmentResult ) *readInput->saRangeTotalNum );
            free ( readInput->sa_list );
            readInput->sa_list = new_sa_list;
            readInput->saListSize = new_size;
        }

        // PESRAAlignmentResult* saListPtr = readInput->sa_list+readInput->saRangeTotalNum;
        // memcpy(saListPtr,saList,sizeof(PESRAAlignmentResult)*saNum);

        for ( unsigned int i = 0; i < saNum; i++ )
        {
            readInput->sa_list[i + readInput->saRangeTotalNum] = saList[i];
            readInput->sa_list[i + readInput->saRangeTotalNum].readID = readid;
        }

        readInput->saRangeTotalNum += saNum;
    }

    if ( occNum > 0 )
    {
        if ( readInput->occTotalNum + occNum > readInput->occListSize )
        {
            // enlarge the arrays "occ_list" by at least double
            unsigned int new_size = readInput->occListSize * 2;

            while ( readInput->occTotalNum + occNum > new_size )
            {
                new_size *= 2;
            }

            SRAOccurrence * new_occ_list = ( SRAOccurrence * ) malloc ( sizeof ( SRAOccurrence ) * new_size );
            memcpy ( new_occ_list, readInput->occ_list, sizeof ( SRAOccurrence ) *readInput->occTotalNum );
            free ( readInput->occ_list );
            readInput->occ_list = new_occ_list;
            readInput->occListSize = new_size;
        }

        // SRAOccurrence* occListPtr = readInput->occ_list+readInput->occTotalNum;
        // memcpy(occListPtr,occList,sizeof(SRAOccurrence)*occNum);

        for ( unsigned int i = 0; i < occNum; i++ )
        {
            readInput->occ_list[i + readInput->occTotalNum] = occList[i];
            readInput->occ_list[i + readInput->occTotalNum].readID = readid;
        }

        readInput->occTotalNum += occNum;
    }

    // if (saNum > 0 || occNum > 0)
    unsigned int evenid = readid / 2 * 2;

    if ( evenid != readInput->lastReadID )
    {
        readInput->readNum++;
        readInput->lastReadID = evenid;
    }
}

DynamicUint8Array * DynamicUint8ArrayConstruct ()
{
    DynamicUint8Array * newuint8Array = ( DynamicUint8Array * ) malloc ( sizeof ( DynamicUint8Array ) );
    newuint8Array->charStr = ( uint8_t * ) malloc ( sizeof ( uint8_t ) * INITIAL_SIZE_CHAR_ARRAY );
    newuint8Array->charStr[0] = '\0';
    newuint8Array->size = INITIAL_SIZE_CHAR_ARRAY;
    newuint8Array->length = 0;
    return newuint8Array;
}

void DynamicUint8ArrayFree ( DynamicUint8Array * uint8Array )
{
    free ( uint8Array->charStr );
    free ( uint8Array );
}

void DynamicUint8ArrayReset ( DynamicUint8Array * uint8Array )
{
    uint8Array->length = 0;
    uint8Array->charStr[0] = '\0';
}

void appendStringToUint8Array ( DynamicUint8Array * uint8Array, char * charString, int len )
{
    while ( uint8Array->length + len >= uint8Array->size )
    {
        // double the size
        uint8_t * newArray = ( uint8_t * ) malloc ( sizeof ( uint8_t ) * uint8Array->size * 2 );
        memcpy ( newArray, uint8Array->charStr, uint8Array->length );
        free ( uint8Array->charStr );
        uint8Array->charStr = newArray;
        uint8Array->size = uint8Array->size * 2;
    }

    memcpy ( uint8Array->charStr + uint8Array->length, charString, len );
    uint8Array->length += len;
    uint8Array->charStr[uint8Array->length] = '\0';
}

// ================================================
// To hold a set of single-end alignment results
// ================================================

AllHits * constructAllHits ()
{
    AllHits * newAllHits = ( AllHits * ) malloc ( sizeof ( AllHits ) );
    newAllHits->hitArray = ( Algnmt * ) malloc ( sizeof ( Algnmt ) * SIZE_1_M );
    newAllHits->hitNum = 0;
    newAllHits->hitArrayAvailSize = SIZE_1_M;
    newAllHits->readPtrArray = ( ReadPtr * ) malloc ( sizeof ( ReadPtr ) * SIZE_1_M );
    newAllHits->readNum = 0;
    newAllHits->readArrayAvailSize = SIZE_1_M;
    return newAllHits;
}

void increaseHitArraySize ( AllHits * allHits, int new_size )
{
    // increase the size of allHits->hitArray
    Algnmt * newHitArray = ( Algnmt * ) malloc ( sizeof ( Algnmt ) * new_size );

    if ( allHits->hitNum > 0 )
    {
        memcpy ( newHitArray, allHits->hitArray, allHits->hitNum * sizeof ( Algnmt ) );
    }

    free ( allHits->hitArray );
    allHits->hitArray = newHitArray;
    allHits->hitArrayAvailSize = new_size;
}

void increaseReadPtrArraySize ( AllHits * allHits, int new_size )
{
    // increase the size of allHits->readPtrArray
    ReadPtr * newReadPtrArray = ( ReadPtr * ) malloc ( sizeof ( ReadPtr ) * new_size );

    if ( allHits->readNum > 0 )
    { memcpy ( newReadPtrArray, allHits->readPtrArray, allHits->readNum * sizeof ( ReadPtr ) ); }

    free ( allHits->readPtrArray );
    allHits->readPtrArray = newReadPtrArray;
    allHits->readArrayAvailSize = allHits->readArrayAvailSize * 2;
}

void inputAlgnmtsToArray ( AllHits * allHits, SingleAlgnmtResult * algnResults, int algnNum )
{
    int i;
    int new_size;

    // check whether hitArray's size is large enough
    if ( allHits->hitNum + algnNum > allHits->hitArrayAvailSize )
    {
        // too small, double the size of allHits->hitArray
        new_size = ( allHits->hitNum + algnNum ) * 2;
        increaseHitArraySize ( allHits, new_size );
    }

    unsigned int preReadID = 0xFFFFFFFF;

    // insert the alignment into the array
    for ( i = 0; i < algnNum; i++ )
    {
        if ( algnResults[i].readID != preReadID )
        {
            // new read ID
            if ( allHits->readNum >= allHits->readArrayAvailSize )
            {
                new_size = 2 * allHits->readArrayAvailSize;
                increaseReadPtrArraySize ( allHits, new_size );
            }

            // create the new readPtr
            allHits->readPtrArray[allHits->readNum].numAlgnmt = 0;
            allHits->readPtrArray[allHits->readNum].readID = algnResults[i].readID;
            allHits->readPtrArray[allHits->readNum].startIndex = allHits->hitNum;
            allHits->readNum++;
            preReadID = algnResults[i].readID;
        }

        if ( algnResults[i].algnmt != 0xFFFFFFFF )
        {
            allHits->hitArray[allHits->hitNum].cigarString = ( char * ) malloc ( strlen ( algnResults[i].cigarString ) + 1 );
            strcpy ( allHits->hitArray[allHits->hitNum].cigarString, algnResults[i].cigarString );
            allHits->hitArray[allHits->hitNum].algnmt = algnResults[i].algnmt;
            allHits->hitArray[allHits->hitNum].strand = algnResults[i].strand;
            allHits->hitArray[allHits->hitNum].score = algnResults[i].score;
            allHits->hitArray[allHits->hitNum].editdist = algnResults[i].editdist;
            allHits->hitArray[allHits->hitNum].num_sameScore = algnResults[i].num_sameScore;
            allHits->hitArray[allHits->hitNum].isFromDP = 1;
            allHits->hitNum++;
            allHits->readPtrArray[allHits->readNum - 1].numAlgnmt++;
        }
    }
}

void sortReadPtrs ( AllHits * allHits )
{
    // sort allHits->readPtrArray
    // define the radix sort
#define NUM_BUCKETS 65536
#define SORT_DIGIT(array, var, shift_bits, num_elements, buckets, buffer) \
    do { memset(buckets, 0, NUM_BUCKETS * sizeof(unsigned int));          \
        for (unsigned int i=0; i<num_elements; ++i) buckets[(array[i].var>>shift_bits)&0xFFFF]++; \
        for (unsigned int b=0, acc=0; b<NUM_BUCKETS; ++b) {unsigned int t=acc; acc+=buckets[b]; buckets[b]=t;} \
        for (unsigned int i=0; i<num_elements; ++i) buffer[buckets[(array[i].var>>shift_bits)&0xFFFF]++] = array[i]; \
        { ReadPtr *t = buffer; buffer = array; array = t; } } while (0)
    // use radix sort to sort the numbers
    unsigned int * buckets = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * NUM_BUCKETS );
    ReadPtr * buffer = ( ReadPtr * ) malloc ( sizeof ( ReadPtr ) * allHits->readArrayAvailSize );
    // first pass: pos (lower 16 bits)
    SORT_DIGIT ( allHits->readPtrArray, readID, 0, allHits->readNum, buckets, buffer );
    // second pass: pos (upper 16 bits)
    SORT_DIGIT ( allHits->readPtrArray, readID, 16, allHits->readNum, buckets, buffer );
    // release the array
    free ( buckets );
    free ( buffer );
}

void resetAllHits ( AllHits * allHits )
{
    int i;

    for ( i = 0; i < allHits->hitNum; i++ )
    {
        free ( allHits->hitArray[i].cigarString );
    }

    allHits->hitNum = 0;
    allHits->readNum = 0;
}

void releaseAllHits ( AllHits * allHits )
{
    int i;

    for ( i = 0; i < allHits->hitNum; i++ )
    {
        free ( allHits->hitArray[i].cigarString );
    }

    free ( allHits->hitArray );
    free ( allHits->readPtrArray );
}

void printOutHits ( AllHits * allHits ) // for debugging
{
    int i, j;

    for ( i = 0; i < allHits->readNum; i++ )
    {
        fprintf ( stderr, "%u - ", allHits->readPtrArray[i].readID );

        for ( j = 0; j < allHits->readPtrArray[i].numAlgnmt; j++ )
        {
            Algnmt curr_algn = allHits->hitArray[allHits->readPtrArray[i].startIndex + j];
            fprintf ( stderr, "algnmt: %u; strand: %i; editdist: %i;", curr_algn.algnmt, ( int ) curr_algn.strand, curr_algn.editdist );
        }

        fprintf ( stderr, "\n" );
    }
}

int writeNumToStr2 ( int num, char * str )
{
    // write number to string
    // return size of the text of the number
    char * currpos = str;

    if ( num < 0 )
    {
        num = -num;
        ( *currpos ) = '-';
        currpos++;
    }

    int digit = 1;

    if ( num > 0 )
    {
        digit = ( ( int ) log10 ( num ) ) + 1;
    }

    for ( int i = digit - 1; i >= 0; i-- )
    {
        currpos[i] = num % 10 + '0';
        num = num / 10;
    }

    return digit + ( currpos - str );
}

void inputSoap3AnsToArray ( AllHits * allHits, unsigned int readID, HSPAux * hspaux, BWT * bwt, int readLength )
{
    // input soap3 answer to array
    // create the new readPtr
    int new_size;

    if ( allHits->readNum >= allHits->readArrayAvailSize )
    {
        new_size = 2 * allHits->readArrayAvailSize;
        increaseReadPtrArraySize ( allHits, new_size );
    }

    allHits->readPtrArray[allHits->readNum].numAlgnmt = 0;
    allHits->readPtrArray[allHits->readNum].readID = readID;
    allHits->readPtrArray[allHits->readNum].startIndex = allHits->hitNum;
    allHits->readNum++;
    unsigned int sa_num = hspaux->sa_num[readID];
    unsigned int occ_num = hspaux->occ_num[readID];
    // fprintf(stderr, "[%u] sa_num = %u; occ_num = %u\n", readID, sa_num, occ_num);
    ReadInputForDP * ans = ( ( ReadInputForDP ** ) hspaux->soap3AnsArray ) [readID];
    unsigned int i, k, l, r, n;
    unsigned char strand, num_mis;
    char * cigar = NULL;

    if ( ans != NULL && ( sa_num > 0 || occ_num > 0 ) )
    {
        // create the cigar
        cigar = ( char * ) malloc ( ( ( int ) log10 ( readLength ) ) + 3 );
        int pos = writeNumToStr2 ( readLength, cigar );
        cigar[pos++] = 'M';
        cigar[pos++] = '\0';
    }

    if ( ans != NULL && sa_num > 0 )
    {
        // insert all the SA ranges
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

            if ( r >= l )
            {
                n = r - l + 1;

                if ( allHits->hitNum + n > allHits->hitArrayAvailSize )
                {
                    // too small, increase the size of allHits->hitArray
                    new_size = ( allHits->hitNum + n ) * 2;
                    increaseHitArraySize ( allHits, new_size );
                }

                for ( k = l; k <= r; k++ )
                {
                    allHits->hitArray[allHits->hitNum].cigarString = ( char * ) malloc ( strlen ( cigar ) + 1 );
                    strcpy ( allHits->hitArray[allHits->hitNum].cigarString, cigar );
                    allHits->hitArray[allHits->hitNum].algnmt = ( *bwt->_bwtSaValue ) ( bwt, k );
                    allHits->hitArray[allHits->hitNum].strand = strand;
                    allHits->hitArray[allHits->hitNum].score = readLength * hspaux->dpMatchScore + num_mis * hspaux->dpMisMatchScore;
                    allHits->hitArray[allHits->hitNum].editdist = num_mis;
                    allHits->hitArray[allHits->hitNum].num_sameScore = 0;
                    allHits->hitArray[allHits->hitNum].isFromDP = 0;
                    allHits->hitNum++;
                    allHits->readPtrArray[allHits->readNum - 1].numAlgnmt++;
                }
            }
        }
    }

    if ( ans != NULL && occ_num > 0 )
    {
        // insert all occ
        if ( allHits->hitNum + occ_num > allHits->hitArrayAvailSize )
        {
            // too small, increase the size of allHits->hitArray
            new_size = ( allHits->hitNum + occ_num ) * 2;
            increaseHitArraySize ( allHits, new_size );
        }

        unsigned int occ_start = hspaux->occ_start[readID];

        for ( i = occ_start; i < occ_start + occ_num; i++ )
        {
            if ( ans->occ_list[i].readID != readID )
            {
                fprintf ( stderr, "[occ] Error! readID does not match!" );
                exit ( 1 );
            }

            allHits->hitArray[allHits->hitNum].cigarString = ( char * ) malloc ( strlen ( cigar ) + 1 );
            strcpy ( allHits->hitArray[allHits->hitNum].cigarString, cigar );
            allHits->hitArray[allHits->hitNum].algnmt = ans->occ_list[i].ambPosition;
            allHits->hitArray[allHits->hitNum].strand = ans->occ_list[i].strand;
            allHits->hitArray[allHits->hitNum].score = readLength * hspaux->dpMatchScore + ans->occ_list[i].mismatchCount * hspaux->dpMisMatchScore;
            allHits->hitArray[allHits->hitNum].editdist = ans->occ_list[i].mismatchCount;
            allHits->hitArray[allHits->hitNum].num_sameScore = 0;
            allHits->hitArray[allHits->hitNum].isFromDP = 0;
            allHits->hitNum++;
            allHits->readPtrArray[allHits->readNum - 1].numAlgnmt++;
        }
    }

    if ( ans != NULL && ( sa_num > 0 || occ_num > 0 ) )
    {
        free ( cigar );
    }
}
