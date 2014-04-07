/*
 *
 *    AlgnResult.h
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

#include "AlgnResult.h"

///////////////////////////////////////////////////////////
// Structures for storing single-end alignment result    //
///////////////////////////////////////////////////////////

// Construct the structure SingleAlgnResultArray
SingleAlgnResultArray * constructSingleAlgnResultArray ( unsigned int numCPUThreads )
{
    SingleAlgnResultArray * newArray;
    newArray = ( SingleAlgnResultArray * ) malloc ( sizeof ( SingleAlgnResultArray ) );
    newArray->arrayNum = numCPUThreads;
    newArray->array = ( SingleAlgnResult ** ) malloc ( sizeof ( SingleAlgnResult * ) * numCPUThreads );

    for ( unsigned int i = 0; i < numCPUThreads; i++ )
    {
        newArray->array[i] = constructSingleAlgnResult ();
    }

    return newArray;
}

// Release the memory for the structure SingleAlgnResultArray
void freeSingleAlgnResultArray ( SingleAlgnResultArray * algnResultArray )
{
    for ( unsigned int i = 0; i < algnResultArray->arrayNum; i++ )
    {
        freeSingleAlgnResult ( algnResultArray->array[i] );
    }

    free ( algnResultArray );
}


// Construct the structure SingleAlgnResult for each CPU thread
SingleAlgnResult * constructSingleAlgnResult ()
{
    SingleAlgnResult * algnResult;
    algnResult = ( SingleAlgnResult * ) malloc ( sizeof ( SingleAlgnResult ) );
    algnResult->readNum = 0;
    // SA List and OCC List
    algnResult->sa_list = ( SARecord * ) malloc ( INITIAL_SIZE_SA_LIST_FOR_SINGLE * sizeof ( SARecord ) );
    algnResult->occ_list = ( OccRecord * ) malloc ( INITIAL_SIZE_OCC_LIST_FOR_SINGLE * sizeof ( OccRecord ) );
    algnResult->saTotalNum = 0;
    algnResult->occTotalNum = 0;
    algnResult->saSize = INITIAL_SIZE_SA_LIST_FOR_SINGLE;
    algnResult->occSize = INITIAL_SIZE_OCC_LIST_FOR_SINGLE;
    // List for read information
    algnResult->readIDs = ( unsigned int * ) malloc ( INITIAL_SIZE_READ_FOR_SINGLE * sizeof ( unsigned int ) );
    algnResult->saEnds = ( unsigned int * ) malloc ( INITIAL_SIZE_READ_FOR_SINGLE * sizeof ( unsigned int ) );
    algnResult->occEnds = ( unsigned int * ) malloc ( INITIAL_SIZE_READ_FOR_SINGLE * sizeof ( unsigned int ) );
    algnResult->isTooManyHit = ( unsigned char * ) malloc ( INITIAL_SIZE_READ_FOR_SINGLE * sizeof ( unsigned int ) );
    algnResult->readSize = INITIAL_SIZE_READ_FOR_SINGLE;
    return algnResult;
}

// Release the memory for the structure SingleAlgnResultP
void freeSingleAlgnResult ( SingleAlgnResult * algnResult )
{
    free ( algnResult->sa_list );
    free ( algnResult->occ_list );
    free ( algnResult->readIDs );
    free ( algnResult->saEnds );
    free ( algnResult->occEnds );
    free ( algnResult->isTooManyHit );
    free ( algnResult );
}

// Reset the structure SingleAlgnResult
void resetSingleAlgnResult ( SingleAlgnResult * algnResult )
{
    algnResult->readNum = 0;
    algnResult->saTotalNum = 0;
    algnResult->occTotalNum = 0;
}

// Add a SA Range to the structure SingleAlgnResult
void addSAToAlgnResult ( SingleAlgnResult * algnResult, unsigned int l, unsigned int r, unsigned char strand )
{
    if ( algnResult->saTotalNum == algnResult->saSize )
    {
        // enlarge the arrays "sa_list" by at least double
        unsigned int new_size = algnResult->saSize * 2;
        SARecord * new_sa_list = ( SARecord * ) malloc ( sizeof ( SARecord ) * new_size );
        memcpy ( new_sa_list, algnResult->sa_list, sizeof ( SARecord ) *algnResult->saTotalNum );
        free ( algnResult->sa_list );
        algnResult->sa_list = new_sa_list;
        algnResult->saSize = new_size;
    }

    algnResult->sa_list[algnResult->saTotalNum].saLeft = l;
    algnResult->sa_list[algnResult->saTotalNum].saRight = r;
    algnResult->sa_list[algnResult->saTotalNum].strand = strand;
    algnResult->saTotalNum++;
}

// Add an occurrence to the structure SingleAlgnResult
void addOccToAlgnResult ( SingleAlgnResult * algnResult, unsigned int pos, unsigned char strand )
{
    if ( algnResult->occTotalNum == algnResult->occSize )
    {
        // enlarge the arrays "occ_list" by at least double
        unsigned int new_size = algnResult->occSize * 2;
        OccRecord * new_occ_list = ( OccRecord * ) malloc ( sizeof ( OccRecord ) * new_size );
        memcpy ( new_occ_list, algnResult->occ_list, sizeof ( OccRecord ) *algnResult->occTotalNum );
        free ( algnResult->occ_list );
        algnResult->occ_list = new_occ_list;
        algnResult->occSize = new_size;
    }

    algnResult->occ_list[algnResult->occTotalNum].pos = pos;
    algnResult->occ_list[algnResult->occTotalNum].strand = strand;
    algnResult->occTotalNum++;
}

// Add read info to the structure SingleAlgnResult
void addReadInfoToAlgnResult ( SingleAlgnResult * algnResult, unsigned int readID, unsigned char status )
{
    // status: 0 - no hit
    // status: 1 - there exists records only in sa list
    // status: 2 - there exists records only in occ list
    // status: 3 - there exists reconds in both sa lista nd occ list
    // status: 4 - too many hits, thus no record stored in sa list or occ list
    if ( algnResult->readNum == algnResult->readSize )
    {
        // enlarge the arrays "readIDs", "saEnds" and "occEnds" by at least double
        unsigned int new_size = algnResult->readSize * 2;
        unsigned int * new_readIDs = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * new_size );
        unsigned int * new_saEnds = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * new_size );
        unsigned int * new_occEnds = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * new_size );
        memcpy ( new_readIDs, algnResult->readIDs, sizeof ( unsigned int ) *algnResult->readNum );
        memcpy ( new_saEnds, algnResult->saEnds, sizeof ( unsigned int ) *algnResult->readNum );
        memcpy ( new_occEnds, algnResult->occEnds, sizeof ( unsigned int ) *algnResult->readNum );
        free ( algnResult->readIDs );
        free ( algnResult->saEnds );
        free ( algnResult->occEnds );
        algnResult->readIDs = new_readIDs;
        algnResult->saEnds = new_saEnds;
        algnResult->occEnds = new_occEnds;
        algnResult->readSize = new_size;
    }

    algnResult->readIDs[algnResult->readNum] = readID;
    algnResult->saEnds[algnResult->readNum] = ( status & 1 ) ? algnResult->saTotalNum - 1 : 0xFFFFFFFF;
    algnResult->occEnds[algnResult->readNum] = ( ( status & 2 ) >> 1 ) ? algnResult->occTotalNum - 1 : 0xFFFFFFFF;
    algnResult->isTooManyHit[algnResult->readNum] = ( status == 4 ) ? 1 : 0;
    algnResult->readNum++;
}

// print the SA information inside the structure SingleAlgnResult
void printSA ( SingleAlgnResult * algnResult )
{
    unsigned int currRead = 0;

    for ( unsigned int i = 0; i < algnResult->saTotalNum; i++ )
    {
        while ( ( algnResult->saEnds[currRead] == 0xFFFFFFFF ) || ( algnResult->saEnds[currRead] < i ) )
        {
            currRead++;
        }

        printf ( "%u %u %u %i\n", algnResult->readIDs[currRead], algnResult->sa_list[i].saLeft,
                 algnResult->sa_list[i].saRight, ( int ) algnResult->sa_list[i].strand );
    }
}

// print all SA information inside the structure SingleAlgnResultArray
void printAllSA ( SingleAlgnResultArray * algnResultArray )
{
    for ( unsigned int i = 0; i < algnResultArray->arrayNum; i++ )
    {
        printSA ( algnResultArray->array[i] );
    }
}


// Construct the structure BothUnalignedPairs
BothUnalignedPairs * constructBothUnalignedPairs ()
{
    BothUnalignedPairs * bothUnalignedPairs = ( BothUnalignedPairs * ) malloc ( sizeof ( BothUnalignedPairs ) );
    bothUnalignedPairs->readIDs = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * INITIAL_SIZE_READ_ID_FOR_BOTH_UNALIGN );
    bothUnalignedPairs->totalNum = 0;
    bothUnalignedPairs->size = INITIAL_SIZE_READ_ID_FOR_BOTH_UNALIGN;
    return bothUnalignedPairs;
}

// Release the memory for the structure BothUnalignedPairs
void freeBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs )
{
    free ( bothUnalignedPairs->readIDs );
    free ( bothUnalignedPairs );
}

// Reset the structure BothUnalignedPairsArrays
void resetBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArrays )
{
    for ( unsigned int i = 0; i < bothUnalignedPairsArrays->arrayNum; i++ )
    {
        bothUnalignedPairsArrays->array[i]->totalNum = 0;
    }
}

// Add read id to BothUnalignedPairs
void addReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int readID )
{
    if ( bothUnalignedPairs->totalNum == bothUnalignedPairs->size )
    {
        // enlarge the arrays "readIDs" by at least double
        unsigned int new_size = bothUnalignedPairs->size * 2;
        unsigned int * new_readIDs = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * new_size );
        memcpy ( new_readIDs, bothUnalignedPairs->readIDs, sizeof ( unsigned int ) *bothUnalignedPairs->totalNum );
        free ( bothUnalignedPairs->readIDs );
        bothUnalignedPairs->readIDs = new_readIDs;
        bothUnalignedPairs->size = new_size;
    }

    bothUnalignedPairs->readIDs[bothUnalignedPairs->totalNum] = readID;
    bothUnalignedPairs->totalNum++;
}

// Add all first reads of the pairs to BothUnalignedPairs (i.e. 0, 2, 4, ..., totalReadNum-2)
// assume totalReadNum is an even number
void addAllFirstReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int totalReadNum )
{
    unsigned int newReadNum = totalReadNum / 2;

    if ( bothUnalignedPairs->totalNum + newReadNum > bothUnalignedPairs->size )
    {
        // enlarge the arrays such that it is enough for storing all the new read IDs
        unsigned int new_size = ( bothUnalignedPairs->totalNum + newReadNum ) * 1.2;
        unsigned int * new_readIDs = ( unsigned int * ) malloc ( sizeof ( unsigned int ) * new_size );
        memcpy ( new_readIDs, bothUnalignedPairs->readIDs, sizeof ( unsigned int ) *bothUnalignedPairs->totalNum );
        free ( bothUnalignedPairs->readIDs );
        bothUnalignedPairs->readIDs = new_readIDs;
    }

    unsigned int i;

    for ( i = 0; i <= totalReadNum - 2; i += 2 )
    {
        bothUnalignedPairs->readIDs[bothUnalignedPairs->totalNum] = i;
        bothUnalignedPairs->totalNum++;
    }
}

// Construct the structure BothUnalignedPairsArrays
BothUnalignedPairsArrays * constructBothUnalignedPairsArrays ( unsigned int numCPUThreads )
{
    BothUnalignedPairsArrays * newArray;
    newArray = ( BothUnalignedPairsArrays * ) malloc ( sizeof ( BothUnalignedPairsArrays ) );
    newArray->arrayNum = numCPUThreads;
    newArray->array = ( BothUnalignedPairs ** ) malloc ( sizeof ( BothUnalignedPairs * ) * numCPUThreads );

    for ( unsigned int i = 0; i < numCPUThreads; i++ )
    {
        newArray->array[i] = constructBothUnalignedPairs ();
    }

    return newArray;
}

// Release the memory for the structure SingleAlgnResultArray
void freeBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArray )
{
    for ( unsigned int i = 0; i < bothUnalignedPairsArray->arrayNum; i++ )
    {
        freeBothUnalignedPairs ( bothUnalignedPairsArray->array[i] );
    }

    free ( bothUnalignedPairsArray );
}


// print out the read id inside BothUnalignedPairs
void printReadIDs ( BothUnalignedPairs * bothUnalignedPairs, unsigned long long accumReadNum, unsigned int * readIDs )
{
    for ( unsigned int i = 0; i < bothUnalignedPairs->totalNum; i++ )
    {
        fprintf ( stderr, "%llu\n", ( unsigned long long ) readIDs[bothUnalignedPairs->readIDs[i]] + accumReadNum );
    }
}


// print out the read id inside BothUnalignedPairsArrays
void printAllReadIDs ( BothUnalignedPairsArrays * bothUnalignedPairsArrays, unsigned long long accumReadNum, unsigned int * readIDs )
{
    for ( unsigned int i = 0; i < bothUnalignedPairsArrays->arrayNum; i++ )
    {
        printReadIDs ( bothUnalignedPairsArrays->array[i], accumReadNum, readIDs );
    }
}
