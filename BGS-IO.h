/*
 *
 *    BGS-IO.h
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

#ifndef __BGS_IO__
#define __BGS_IO__

#include <stdio.h>
#include <stdlib.h>
#include "PEAlgnmt.h"
#include "2bwt-flex/SRAArguments.h"
#include "SAM.h"
#include "SAList.h"
#include "PE.h"
#include "Release.h"

//Define the below parameter to output the alignment result(text position)
// on screen instead of writing into output file
//#define DEBUG_2BWT_OUTPUT_TO_SCREEN

//Define the below parameter to stop cache reported SA range and text position.
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_OUTPUT

//Define the below parameter to skip writing the alignment result(text position)
// into disk. This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_WRITE_FILE

//Define the below parameter to skip translating the alignment result(text position).
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_TRANSLATION

#define OCC_CONST_ALIGNMENT_HEADER 65535
#define OCC_CONST_NO_ALIGNMENT     65534

//Define below parameter to skip outputing the header for plain output
#define SKIP_PLAIN_HEADER

OCC * OCCConstruct ();
void OCCFree ( OCC * occ );



unsigned int OCCWriteOutputHeader ( HSP * hsp, FILE * outFilePtr,
                                    unsigned int maxReadLength,
                                    unsigned int numOfReads,
                                    int outputFormat );


// succinct format
void OCCFlushCache ( SRAQueryInput * qInput );
void OCCFlushCacheDP ( SRAQueryInput * qInput ); // for outputting the DP result

// binary
void OCCFlushCacheDefault ( OCC * occ, HSP * hsp, FILE * outFilePtr );

// obseleted (NO, DONT BELIEVE)
// -- start --
void OCCFlushCacheSAM ( SRAQueryInput * qInput );
/*
void OCCDirectWritePairOccSAM ( SRAQueryInput * qInput, PEPairs * pePair );
void OCCDirectWritePairUnmapSAM ( SRAQueryInput * qInput, SRAOccurrence * sraOcc, int isFirst );
*/
void OCCDirectWritePairOccSAMAPI ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, int properlyPaired );
void OCCDirectWritePairOccSAMAPI2 ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, int properlyPaired );
void OCCDirectWritePairUnmapSAMAPI ( SRAQueryInput * qInput, SRAOccurrence * sraOcc, int isFirst );
void OCCDirectWritePairOccSAMAPIwCIGAR ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, char * cigar_str );
void OCCDirectWritePairOccSAMAPI2wCIGAR ( SRAQueryInput * qInput, PEPairs * pePair, int isFirstUnaligned, int isSecondUnaligned, int mateLength, char * cigar_str );
// -- end --
void OCCReportDelimitor ( SRAQueryInput * qInput );
void OCCReportDelimitorDP ( SRAQueryInput * qInput ); // for outputting the DP result
/*
//////////////////////////
//    FOR SAM FORMAT    //
//////////////////////////


// single-end alignment
// -- start --
*/
void OCCOutputSAMAPI ( SRAQueryInput * qInput, OCCList * occ_list,
                       DynamicUint8Array * xazArray, int readlen, int reportType );
/*
// For outputting single-read alignment in SAM format
// Output the best record, and for the rest, append the results in the tag XA:Z
*/
void SingleAnsOutputSAMAPI ( SRAQueryInput * qInput,
                             char strand, unsigned int ambPosition, int bestHitNumOfMismatch, int bestHitNum );
/*
// For outputting one answer of single-read alignment in SAM format
// for unique-best and random-best single-read alignment
*/
void noAnsOutputSAMAPI ( SRAQueryInput * qInput );

// For outputting no answer of single-read alignment in SAM format

void SingleDPOutputSAMAPI ( SRAQueryInput * qInput, SingleAlgnmtResult * algnResult,
                            unsigned int startIndex, unsigned int numResult,
                            DynamicUint8Array * xazArray );
// For outputting single-read alignment in SAM format
// Output the best record, and for the rest, append the results in the tag XA:Z
// -- end --


void pairOutputSAMAPI ( SRAQueryInput * qInput, PEOutput * pe_out, PEPairs * bestPair,
                        unsigned char * query1, unsigned char * query2,
                        char * qualities1, char * qualities2,
                        int readlen1, int readlen2,
                        char * queryName1, char * queryName2,
                        char min_totalMismatchCount, char secMin_totalMismatchCount,
                        int outputXAZTag, DynamicUint8Array * xazArray,
                        unsigned int peMaxOutputPerPair,
                        int X0_first, int X0_second, int X1_first, int X1_second, int num_minMismatchPair,
                        char isBestHit1, char isBestHit2, unsigned int totalNumValidPairs );
/*
// For outputting pair-end alignment in SAM format
// Output the best record, and for the rest, append the results in the tag XA:Z

*/
void pairDPOutputSAMAPI ( SRAQueryInput * qInput, AlgnmtDPResult * algnResult,
                          AlgnmtDPResult * bestResult,
                          unsigned int start, unsigned int num,
                          unsigned char * query1, unsigned char * query2,
                          char * qualities1, char * qualities2,
                          int readlen1, int readlen2,
                          char * queryName1, char * queryName2,
                          DynamicUint8Array * xazArray );
// For outputting DP pair-end alignment in SAM format

void pairDeepDPOutputSAMAPI ( SRAQueryInput * qInput, DeepDPAlignResult * algnResult,
                              DeepDPAlignResult * bestResult,
                              unsigned int start, unsigned int num,
                              unsigned char * query1, unsigned char * query2,
                              char * qualities1, char * qualities2,
                              int readlen1, int readlen2,
                              char * queryName1, char * queryName2,
                              DynamicUint8Array * xazArray );
// For outputting Deep DP pair-end alignment in SAM format

void unproperlypairOutputSAMAPI ( SRAQueryInput * qInput, OCCList * occ_list1, OCCList * occ_list2,
                                  unsigned char * query1, unsigned char * query2,
                                  char * qualities1, char * qualities2,
                                  int readlen1, int readlen2,
                                  char * queryName1, char * queryName2,
                                  DynamicUint8Array * xazArray,
                                  unsigned int peMaxOutputPerRead, int reportType );
// output the alignments which are not properly paired

void unproperlypairDPOutputSAMAPI ( SRAQueryInput * qInput, Algnmt * algn_list1, Algnmt * algn_list2,
                                    int hitNum1, int hitNum2,
                                    unsigned char * query1, unsigned char * query2,
                                    char * qualities1, char * qualities2,
                                    int readlen1, int readlen2,
                                    char * queryName1, char * queryName2,
                                    DynamicUint8Array * xazArray );
// output the DP alignments which are not properly paired



void OCCReportNoAlignment ( SRAQueryInput * qInput );
/*
//OCCReportSARange : asumming l<=r

unsigned long long OCCReportSARange ( SRAQueryInput * qInput,
                                      unsigned long long l, unsigned long long r,
                                      int occMismatch );

unsigned long long OCCReportTextPositions ( SRAQueryInput * qInput, int posCount );
*/
#endif
