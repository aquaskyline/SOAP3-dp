/*
 *
 *    DV-SemiDP.h
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

#ifndef _SEMIDP_H_
#define _SEMIDP_H_

#include "PEAlgnmt.h"
#include "definitions.h"
#include "IndexHandler.h"

///////////////////////////////////////////////////////////
// The functions are for semi-global DP                  //
///////////////////////////////////////////////////////////


// To perform semi-global DP
void semiGlobalDP2 ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                     unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                     unsigned int maxReadLength, int peStrandLeftLeg, int peStrandRightLeg,
                     Soap3Index * index,
                     unsigned int * _bwt, unsigned int * _revBwt,
                     unsigned int * _occ, unsigned int * _revOcc,
                     int alignmentType, DPParameters * dpParameters,
                     unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                     unsigned int accumReadNum, int outputFormat,
                     FILE * outputFile, samfile_t * samOutputDPFilePtr,
                     BothUnalignedPairs * unalignedReads );

// To perform new semi-global DP
void newSemiGlobalDP ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                       unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                       unsigned int maxReadLength, int peStrandLeftLeg, int peStrandRightLeg,
                       Soap3Index * index,
                       unsigned int * _bwt, unsigned int * _revBwt,
                       unsigned int * _occ, unsigned int * _revOcc,
                       /* end */
                       int alignmentType, DPParameters * dpParameters,
                       unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                       unsigned int accumReadNum, int outputFormat,
                       FILE * outputFile, samfile_t * samOutputDPFilePtr,
                       BothUnalignedPairs * unalignedReads );


// readInputArray: a set of ReadInputForDP array contains SOAP3 alignment results of the reads for DP.
// insert_high : maximum value of insert size
// insert_low: minimum value of insert size

// queries: the sequence of all reads. Each character is represented by 2 bits.
// like: [read 1] [read 2] ... [read 32] [read 1] ....... [read 32] [read 33] [read 34] ...
//        16...1   16...1       16....1   32..17           32...17   16....1   16....1

// readLengths: the read lengths of all reads. The read length of read with id i: readLengths[i]
// origReadIDs: the ORIGINAL read IDs of all reads
// upkdQueryNames: the description of all reads (based on original read IDs)
// upkdQualities: the base quality values of all reads (based on original read IDs)

// peStrandLeftLeg: the required strand of the left alignment (i.e. 1: +ve; 2: -ve)
// peStrandRightLeg: the required strand of the right alignment (i.e. 1: +ve; 2: -ve)
// BWT: bwt structure with SA table inside
// HSP: hsp structure with packed sequence inside
// alignmentType: 1: All valid alignment; 4: random best alignment
// dpParameters: the parameters for DP
// Output:
// numDPAlignedRead: the total number of reads aligned by DP
// numDPAlignment: the total number of resulting alignments by DP
// accumReadNum: the accumulated number of reads being processed
// outputFormat: the format of output
// outputFile: the file pointer for outputing the results



// To perform semi-global DP
void semiGlobalDP ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                    unsigned char * upkdQueries, char * upkdQueryNames, unsigned int * upkdReadLengths,
                    unsigned int maxReadLength, int peStrandLeftLeg, int peStrandRightLeg,
                    Soap3Index * index,
                    int alignmentType, DPParameters * dpParameters,
                    unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                    unsigned int accumReadNum, int outputFormat,
                    FILE * outputFile, samfile_t * samOutputDPFilePtr );


// readInputArray: a set of ReadInputForDP array contains SOAP3 alignment results of the reads for DP.
// insert_high : maximum value of insert size
// insert_low: minimum value of insert size
// upkdQueries: the sequence of all reads. The first character of the read with ID "i"
//              is on the position "i*maxReadLength" of the array
// upkdQueryNames: the description of all reads.
// upkdReadLengths: the read lengths of all reads
// peStrandLeftLeg: the required strand of the left alignment (i.e. 1: +ve; 2: -ve)
// peStrandRightLeg: the required strand of the right alignment (i.e. 1: +ve; 2: -ve)
// BWT: bwt structure with SA table inside
// HSP: hsp structure with packed sequence inside
// alignmentType: 1: All valid alignment; 4: random best alignment
// dpParameters: the parameters for DP
// Output:
// numDPAlignedRead: the total number of reads aligned by DP
// numDPAlignment: the total number of resulting alignments by DP
// accumReadNum: the accumulated number of reads being processed
// outputFormat: the format of output
// outputFile: the file pointer for outputing the results



//*********************************************//
// The following functions should be obselete  //
//*********************************************//


// To perform semi-global DP
AlgnmtDPResult * semiGlobalDP ( ReadInputForDP * readInput, int insert_high, int insert_low,
                                unsigned char * upkdQueries, unsigned int * readLengths, unsigned int maxReadLength,
                                int peStrandLeftLeg, int peStrandRightLeg,
                                Soap3Index * index,
                                int alignmentType, unsigned int & numDPResult );

// To release the memory occupied by "AlgnmtDPResult"
void algnmtDPResultFree ( AlgnmtDPResult * algnDPResult );

// To show the result
// This function is just for testing.
// It will not be used in the future
void showResult ( AlgnmtDPResult * algnDPResult, unsigned int numDPResult,
                  unsigned int & numDPAlignedRead );

// To print an alignment result
void printDPResult ( AlgnmtDPResult & result,
                     unsigned char * upkdQueries, unsigned int * readLengths, int maxReadLength,
                     unsigned int * packedDNA, int peStrandLeftLeg, int peStrandRightLeg,
                     FILE * outputFile = stdout );

#endif

