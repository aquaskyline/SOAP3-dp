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

///////////////////////////////////////////////////////////
// Structures for storing single-end alignment result    //
///////////////////////////////////////////////////////////

#ifndef __ALIGN_RESULT_H__
#define __ALIGN_RESULT_H__


#define INITIAL_SIZE_READ_ID_FOR_BOTH_UNALIGN 10485760 // 10M

#define INITIAL_SIZE_SA_LIST_FOR_SINGLE 10485760 // 10M
#define INITIAL_SIZE_OCC_LIST_FOR_SINGLE 10485760 // 10M
#define INITIAL_SIZE_READ_FOR_SINGLE 1048576 // 1M

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"


// To hold one occurrence.
typedef struct OccRecord
{
    unsigned int pos;
    char strand;
} OccRecord;

// To hold one SA range.
typedef struct SARecord
{
    unsigned int saLeft;
    unsigned int saRight;
    unsigned char strand;
} SARecord;

// To hold single-end SOAP3 alignment results of the reads
typedef struct SingleAlgnResult
{
    unsigned int readNum; // number of reads
    SARecord * sa_list; // the list of SA ranges for all reads
    OccRecord * occ_list; // the list of occurrences for all reads
    unsigned int saTotalNum; // the number of SA ranges in sa_list
    unsigned int occTotalNum; // the number of occurrences in occ_list
    unsigned int saSize; // available size of the array sa_list
    unsigned int occSize; // available size of the array occ_list

    // number of items in the following arrays: "readNum"
    unsigned int * readIDs; // ID of each read
    unsigned int * saEnds; // (if not exist, then 0xFFFFFFFF)
    // the index of the LAST SA range on sa_list for each read (note: the index starts from zero)
    unsigned int * occEnds; // (if not exist, then 0xFFFFFFFF)
    // the index of the LAST occurrence on occ_list for each read (note: the index starts from zero)
    unsigned char * isTooManyHit; // 1: too many hits; 0: otherwise

    // available size of the above three arrays
    unsigned int readSize;
} SingleAlgnResult;

// Array of SingleAlgnResult
typedef struct SingleAlgnResultArray
{
    SingleAlgnResult ** array;
    unsigned int arrayNum; // number of arrays
} SingleAlgnResultArray;


// Construct the structure SingleAlgnResultArray
SingleAlgnResultArray * constructSingleAlgnResultArray ( unsigned int numCPUThreads );

// Release the memory for the structure SingleAlgnResultArray
void freeSingleAlgnResultArray ( SingleAlgnResultArray * algnResultArray );

// Construct the structure SingleAlgnResult for each CPU thread
SingleAlgnResult * constructSingleAlgnResult ();

// Release the memory for the structure ReadInputForDP
void freeSingleAlgnResult ( SingleAlgnResult * algnResult );

// Reset the structure ReadInputForDP
void resetSingleAlgnResult ( SingleAlgnResult * algnResult );

// Add a SA Range to the structure SingleAlgnResult
void addSAToAlgnResult ( SingleAlgnResult * algnResult, unsigned int l, unsigned int r, unsigned char strand );

// Add an occurrence to the structure SingleAlgnResult
void addOccToAlgnResult ( SingleAlgnResult * algnResult, unsigned int pos, unsigned char strand );

// Add read info to the structure SingleAlgnResult
void addReadInfoToAlgnResult ( SingleAlgnResult * algnResult, unsigned int readID, unsigned char status );

// print the SA information inside the structure SingleAlgnResult
void printSA ( SingleAlgnResult * algnResult );

// print all SA information inside the structure SingleAlgnResultArray
void printAllSA ( SingleAlgnResultArray * algnResultArray );


// To contain the first (even) read id of the both-unaligned pairs
// (i.e. both ends cannot be aligned)
typedef struct BothUnalignedPairs
{
    unsigned int * readIDs; // contain the first read id of the unaligned pairs
    unsigned int totalNum; // the number of read id inside the array "readIDs".
    unsigned int size; // available size of the array
} BothUnalignedPairs;

// An array of BothUnalignedPairs
typedef struct BothUnalignedPairsArrays
{
    BothUnalignedPairs ** array;
    unsigned int arrayNum; // number of arrays
} BothUnalignedPairsArrays;

typedef BothUnalignedPairs UnalignedSingles;
typedef BothUnalignedPairsArrays UnalignedSinglesArrays;

// Construct the structure BothUnalignedPairs
BothUnalignedPairs * constructBothUnalignedPairs ();

// Release the memory for the structure BothUnalignedPairs
void freeBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs );

// Reset the structure BothUnalignedPairsArrays
void resetBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArrays );

// Add read id to BothUnalignedPairs
void addReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int readID );

// Add all first reads of the pairs to BothUnalignedPairs (i.e. 0, 2, 4, ..., totalReadNum-2)
void addAllFirstReadIDToBothUnalignedPairs ( BothUnalignedPairs * bothUnalignedPairs, unsigned int totalReadNum );

// Construct the structure BothUnalignedPairsArrays
BothUnalignedPairsArrays * constructBothUnalignedPairsArrays ( unsigned int numCPUThreads );

// Release the memory for the structure SingleAlgnResultArray
void freeBothUnalignedPairsArrays ( BothUnalignedPairsArrays * bothUnalignedPairsArray );

// print out the read id inside BothUnalignedPairs
void printReadIDs ( BothUnalignedPairs * bothUnalignedPairs, unsigned long long accumReadNum, unsigned int * readIDs );

// print out the read id inside BothUnalignedPairsArrays
void printAllReadIDs ( BothUnalignedPairsArrays * bothUnalignedPairsArrays, unsigned long long accumReadNum, unsigned int * readIDs );


#endif

