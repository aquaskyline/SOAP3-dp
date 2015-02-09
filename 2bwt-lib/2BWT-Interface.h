/*

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Enhancing the 2BWT library with interface functions
            for basic BWT search.

*/

#ifndef __2BWT_INTERFACE_H__
#define __2BWT_INTERFACE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BWT.h"
#include "HSP.h"
#include "MemManager.h"

#define MAX_INDEX_FILENAME_LENGTH 1024
#define POOLSIZE                  2097152

typedef struct Idx2BWT {
    MMPool * mmPool;
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    unsigned char charMap[256];
    unsigned char complementMap[256];
} Idx2BWT;



//=============================================================================
//  BWTLoad2BWT
//  This function loads up the 2BWT index from memory and return the placeholder
//  Input :     indexFilePrefix - the file name prefix of the 2BWT index files
//                                e.g. "ncbi.genome.fa.index"
//              saFileNameExtension - the file name extension of suffix array in
//                                    the 2BWT index. This string should correspond
//                                    to the last part of the SaValueFileName value 
//                                    in 2bwt-builder.ini when the index was built.
//                                    e.g. ".sa"
//  Output :    function returns the placeholder of 2BWT index
//=============================================================================
Idx2BWT * BWTLoad2BWT(const char * indexFilePrefix, const char * saFileNameExtension);

//=============================================================================
//  BWTFree2BWT
//  This function frees the 2BWT index from memory and destory the placeholder
//  Input :     idx2BWT - placeholder of 2BWT index
//  Output :    none
//=============================================================================
void BWTFree2BWT(Idx2BWT * idx2BWT);


//=============================================================================
//  BWTConvertPattern
//  This function converts a human readable nucleotide sequence (e.g. ACCAACAT) 
//  into a coding scheme recognised by 2BWT index (e.g. 01100103)
//  Input :     idx2BWT - placeholder of 2BWT index
//              patternSource - the human readable nucleotide sequence
//              patternLength - the length of the pattern, obviously
//  Output :    patternDestination - the converted pattern
//=============================================================================
void BWTConvertPattern(Idx2BWT * idx2BWT, const char * patternSource, int patternLength, unsigned char * patternDestination);

//=============================================================================
//  BWTSARangeInitial
//  This function performs the 2BWT search of the first character 
//  of the pattern and returns the resulting SA range to the input parameter 
//  resultSaIndexLeft and resultSaIndexRight.
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//  Output :    resultSaIndexLeft  - lower bound of the resulting SA range
//              resultSaIndexRight - upper bound of the resulting SA range
//=============================================================================
void BWTSARangeInitial(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);
     
//=============================================================================
//  BWTSARangeBackward
//  This function performs the 2BWT backward search for 1 character
//  and returns the resulting SA range to the input parameter 
//  resultSaIndexLeft and resultSaIndexRight.
//  Input :     idx2BWT - placeholder of 2BWT index
//              c   - the character to search for
//  Output :    saIndexLeft  - lower bound of the resulting SA range
//              saIndexRight - upper bound of the resulting SA range
//=============================================================================
void BWTSARangeBackward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight);

//=============================================================================
//  BWTSARangeBackward_Bidirection
//  This function performs the 2BWT backward search for 1 character
//  and returns the resulting SA ranges to the input parameter 
//  saIndexLeft and saIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    saIndexLeft  - lower bound of the resulting SA range
//              saIndexRight - upper bound of the resulting SA range
//              rev_saIndexLeft  - lower bound of the resulting SA range
//              rev_saIndexRight - upper bound of the resulting SA range
//=============================================================================
void BWTSARangeBackward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight);


//=============================================================================
//  BWTSARangeForward_Bidirection
//  This function performs the 2BWT forward search for 1 character
//  and returns the resulting SA ranges to the input parameter 
//  saIndexLeft and saIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    saIndexLeft  - lower bound of the resulting SA range
//              saIndexRight - upper bound of the resulting SA range
//              rev_saIndexLeft  - lower bound of the resulting SA range
//              rev_saIndexRight - upper bound of the resulting SA range
//=============================================================================

void BWTSARangeForward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight);
                        
//=============================================================================
//  BWTAllSARangesBackward
//  This function performs the 2BWT backward search for all characters in the 
//  alphabet set defined and supported, and returns the resulting SA ranges 
//  to the input parameter resultSaIndexLeft and resultSaIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//
//  Output :    resultSaIndexLefts  - lower bound of the resulting SA ranges.
//                                    Array in size of ALPHABET_SIZE.
//              resultSaIndexRights - upper bound of the resulting SA ranges.
//                                    Array in size of ALPHABET_SIZE.
//=============================================================================

void BWTAllSARangesBackward(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        unsigned int *resultSaIndexesLeft, unsigned int *resultSaIndexesRight);
//=============================================================================
//  BWTAllSARangesBackward_Bidirection
//  This function performs the 2BWT backward search for all characters in the 
//  alphabet set defined and supported and returns the resulting SA ranges 
//  to the input parameter resultSaIndexLeft and resultSaIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    resultSaIndexLeft  - lower bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              resultSaIndexRight - upper bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexLeft  - lower bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexRight - upper bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//=============================================================================
void BWTAllSARangesBackward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        const unsigned int rev_saIndexLeft, const unsigned int rev_saIndexRight,
                        unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,
                        unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight);

//=============================================================================
//  BWTAllSARangesForward_Bidirection
//  This function performs the 2BWT forward search for all characters in the 
//  alphabet set defined and supported and returns the resulting SA ranges 
//  to the input parameter resultSaIndexLeft and resultSaIndexRight.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              c       - the character to search for
//
//              saIndexLeft  - lower bound of the intermediate SA ranges
//              saIndexRight - upper bound of the intermediate SA ranges
//              rev_saIndexLeft  - lower bound of the intermediate Rev SA ranges
//              rev_saIndexRight - upper bound of the intermediate Rev SA ranges
//
//  Output :    resultSaIndexLeft  - lower bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              resultSaIndexRight - upper bound of the resulting SA range
//                                   Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexLeft  - lower bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//              rev_resultSaIndexRight - upper bound of the resulting SA range
//                                       Array in size of ALPHABET_SIZE.
//=============================================================================
void BWTAllSARangesForward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        const unsigned int rev_saIndexLeft, const unsigned int rev_saIndexRight,
                        unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,
                        unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight);
                        

//=============================================================================
//  BWTRetrievePositionFromSAIndex
//  This function translates an SA index (an index within a reported SA range)
//  into a text position in the original reference sequence in the format of
//  1-based sequenceId and offset.
//
//  Input :     idx2BWT - placeholder of 2BWT index
//              saIndex - the SA index within an SA range
//  Output :    sequenceId  - 1-based sequence id in the original reference
//                            sequence. Value can be 1,2,3....
//              offset      - 1-based offset in the above sequence. Value 
//                            can be 1,2,3...
//=============================================================================
void BWTRetrievePositionFromSAIndex(Idx2BWT * idx2BWT, 
                                    unsigned int saIndex, 
                                    int * sequenceId, unsigned int * offset);
                                    
#endif

