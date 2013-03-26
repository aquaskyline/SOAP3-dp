/*
 *
 *    OutputDPResult.h
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

#ifndef _OUTPUTDPRESULT_H_
#define _OUTPUTDPRESULT_H_

#include "PEAlgnmt.h"
#include "CPUfunctions.h"
#include "BGS-IO.h"
#include "definitions.h"

#define MC_CeilDivide16(x) ((x+15)>>4)


// output a set of alignment results
// the alignments of the same pair of reads have to be together
// algnNum: the number of results in "algnResult"
void outputRead ( AlgnmtDPResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char * upkdQueryNames,
                  unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                  FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg );


void outputRead2 ( AlgnmtDPResult * algnResult, unsigned int algnNum,
                   unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                   unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                   FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg );


// output a set of deep DP alignment results
// the alignments of the same pair of reads have to be together
void outputDeepDPResult ( DeepDPAlignResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char * upkdQueryNames,
                          unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                          FILE * outputDPFile, samfile_t * samOutputDPFilePtr, HSP * hsp, int peStrandLeftLeg, int peStrandRightLeg );

void outputDeepDPResult2 ( DeepDPAlignResult * algnResult, unsigned int algnNum,
                           unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                           unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                           FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index, int peStrandLeftLeg, int peStrandRightLeg );


// output DP single-end alignment results
// the alignments of the same pair of reads have to be together
// algnNum: the number of results in "algnResult"
void outputDPSingleResult ( SingleAlgnmtResult * algnResult, unsigned int algnNum, unsigned char * upkdQueries, char * upkdQueryNames,
                            unsigned int * upkdReadLengths, unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                            FILE * outputDPFile, samfile_t * samOutputDPFilePtr, HSP * hsp );

void outputDPSingleResult2 ( SingleAlgnmtResult * algnResult, unsigned int algnNum,
                             unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                             unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                             FILE * outputDPFile, samfile_t * samOutputDPFilePtr, Soap3Index * index );


// output the single-end alignment results for the pair-end reads
void outputSingleResultForPairEnds ( AllHits * allHits,
                                     unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                                     unsigned int maxReadLength, unsigned int accumReadNum, int outputFormat,
                                     samfile_t * samOutputDPFilePtr, Soap3Index * index );

#endif
