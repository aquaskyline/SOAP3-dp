/*
 *
 *    global_arrays.h
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

#ifndef __GLOBAL_ARRAYS_H__
#define __GLOBAL_ARRAYS_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INITIAL_SIZE_GLOBAL_ARRAY 10485760 // 10M

//---------------------------------------------------//
// Global Arrays                                     //
// They are used when storing all results in memory, //
// instead of outputting to the files.               //
//---------------------------------------------------//

// a record of occurrence
typedef struct occRec
{
    unsigned int readID; // read ID
    unsigned int ambPosition; // position on packed sequence
    unsigned char strand; // (i.e. 1: +ve; 2: -ve)
    unsigned char source; // (i.e. 1: from SOAP3; 2: from DP)
    char score; // (i.e. mismatch # if from SOAP3; DP score if from DP
    // char* cigar_str; // cigar string
} occRec;

// To hold SOAP3 alignment results of the reads
typedef struct AlgnResult
{
    // unsigned int readNum; // number of reads
    occRec * occ_list; // the list of occurrences for all reads
    unsigned int occTotalNum; // the number of occurrences in occ_list
    unsigned int availableSize; // the available size in occ_list
} AlgnResult;

// To hold a set of arrays of AlgnResult
typedef struct AlgnResultArrays
{
    AlgnResult ** algnArrays;
    unsigned int numArrays; // number of arrays
} AlgnResultArrays;

// functions

// construct the arrays
AlgnResultArrays * resultArraysConstruct ( int num );

// free the arrays
void resultArraysFree ( AlgnResultArrays * algnResultArrays );

// reset the arrays
void resultArraysReset ( AlgnResultArrays * algnResultArrays );

// add occ to the array
void addOCCToArray ( AlgnResult * algnResult, unsigned int readID, unsigned int ambPosition, unsigned char strand,
                     unsigned char source, char score );


#endif

