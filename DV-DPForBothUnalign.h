/*
 *
 *    DV-DPForBothUnalign.h
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

#ifndef _DP_FOR_BOTH_UNALIGN_H_
#define _DP_FOR_BOTH_UNALIGN_H_

#include "AlgnResult.h"
#include "definitions.h"
#include "IndexHandler.h"


//////////////////////////////////////////////////////////////////////////////////////
// The functions are for the DP designed for pairs of which both ends are unaligned //
//////////////////////////////////////////////////////////////////////////////////////

// To perform semi-global DP for both-unaligned pairs (i.e. both ends cannot be aligned)
void DPForUnalignPairs2 ( BothUnalignedPairsArrays * unalignedReads, int insert_high, int insert_low,
                          unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                          unsigned int maxReadLength,
                          int peStrandLeftLeg, int peStrandRightLeg,
                          Soap3Index * index,
                          unsigned int * _bwt, unsigned int * _revBwt,
                          unsigned int * _occ, unsigned int * _revOcc,
                          int alignmentType, DPParameters * dpParameters,
                          unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                          unsigned int accumReadNum, int outputFormat,
                          FILE * outputFile, samfile_t * samOutputDPFilePtr );



// To perform semi-global DP for both-unaligned pairs (i.e. both ends cannot be aligned)
void DPForUnalignPairs ( BothUnalignedPairsArrays * unalignedReads, int insert_high, int insert_low,
                         unsigned char * upkdQueries, char * upkdQueryNames,
                         unsigned int * upkdReadLengths, unsigned int maxReadLength,
                         int peStrandLeftLeg, int peStrandRightLeg,
                         Soap3Index * index,
                         unsigned int * _bwt, unsigned int * _revBwt,
                         unsigned int * _occ, unsigned int * _revOcc,
                         int alignmentType, DPParameters * dpParameters,
                         unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                         unsigned int accumReadNum, int outputFormat,
                         FILE * outputFile, samfile_t * samOutputDPFilePtr );


#endif
