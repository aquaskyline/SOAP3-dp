/*
 *
 *    PE.h
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

#ifndef _PE_H_
#define _PE_H_

#define DEFAULT_QUAL_VALUE 20

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "2bwt-flex/SRACore.h"

#include "definitions.h"
#include "ssse3_popcount.h"

// make the packed DNA for the query
void createQueryPackedDNA ( unsigned char * query, unsigned int query_length, unsigned int * outPackedDNA );

// make the packed DNA for the reverse strand of the query
void createRevQueryPackedDNA ( unsigned char * query, unsigned int query_length, unsigned int * outPackedDNA );

// make the packed DNA for the target
void createTargetPackedDNA ( unsigned int * target, unsigned int pos, unsigned int target_length,
                             unsigned int * outPackedDNA );

// check the number of mismatches (at most max_mismatch)
char numMismatch ( unsigned int * query, unsigned int * target, unsigned int query_length, unsigned int max_mismatch );
int numMismatchNew ( unsigned int * query, unsigned int * target, unsigned int query_length, int threadID );
// SRAEnrichSARanges takes a SA range as input and return the mismatch count and the strand of the SA range.
// void SRAEnrichSARanges(HSP* hsp, BWT* bwt, unsigned int* query, unsigned int* rev_query, unsigned int query_length, unsigned int l, unsigned int r, char max_mismatch, char* num_mis, char* strand);


// It takes a SA range [l,r] as well as the strand as input and return the number of mismatches of the SA range.
// char SRAEnrichSARangeswithStrand(HSP* hsp, BWT* bwt, unsigned int* query, unsigned int* rev_query, unsigned int query_length, unsigned int l, unsigned int r, char strand);
// strand = 1: forward; strand = 2: reverse

// MismatchInPos takes a position and the strand as input and return the mismatch count
// char MismatchInPos(HSP* hsp, BWT* bwt, unsigned int* query, unsigned int* rev_query, unsigned int query_length, unsigned int pos, char strand);

int writeULLToStr ( unsigned long long num, char * str );
// write unsigned long long to string
// return size of the text of the number

int writeNumToStr ( int num, char * str );
// write number to string
// return size of the text of the number

int getMdStr ( HSP * hsp, unsigned char * query, char * qualities, unsigned int query_length, unsigned int pos, char strand, char mismatchNum, char * md_str, int * avg_mismatch_qual, int trim = 0 );
// compute the MD string
// return the size of MD string

int getMisInfoForDP ( HSP * hsp, unsigned char * query, char * qualities, unsigned int query_length,
                      unsigned int pos, char strand, char * special_cigar, char * md_str, int * numMismatch,
                      int * gapOpen, int * gapExt, int * avg_mismatch_qual, int trim = 0 );
// to compute the MD string, the number of mismatches, the number of gap open and the number of gap extension
// return the size of MD string
// this function is designed for DP module
// the md string is a char array with length > query's length * 2
// special_cigar includes the "m" character for mismatches

int convertToCigarStr ( char * special_cigar, char * cigar, int * deletedEnd = NULL );
// to convert the special_cigar into cigar string

int convertToCigarStr2 ( char * cigar );
// to convert the special_cigar into cigar string
// and save back to the cigar variable

#endif
