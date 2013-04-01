//
//    SRA2BWTCheckAndExtend.h
//
//    SOAP2 / 2BWT
//
//    Copyright (C) 2012, HKU
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 2
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/*

   Modification History 
   
   Date   : 8th May 2011
   Author : Edward MK Wu
   Change : New file.

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Thus, changing all references to 2bwt lib to subdirectory.

*/
/////////////////////////////////////////////////////

#ifndef __SRA_2BWT_CHECK_AND_EXTEND__
#define __SRA_2BWT_CHECK_AND_EXTEND__

#include "SRACore.h"

//Sensitive to BIT_PER_CHAR
#define SRA_CE_BIT_MASK        0x5555555555555555ull
#define SRA_CE_BIT_PER_64        64
#define SRA_CE_BIT_PER_WORD    32

void CEDebugPrintPackedSequence(unsigned long long * packedKey,int len);

//GetPackedPatternLong : packs the convertedKey into packedKey, each character taking only BIT_PER_CHAR.
//The function returns the number of cell the packedKey took.
//  The packedKey(the uneven key is right aligned on the last cell) is also returned; 
//  thus it should be an unsigned long long array with at least (SRA_MAX_READ_LENGTH/CHAR_PER_64) elements.
//The preliminary design : A cell = 64bit unsigned long long.
int CEPackPattern(const unsigned char * convertedKey,
                         int start, int len,
                         unsigned long long * packedKey);

//CEPackedMismatchMatching gives the number of mismatch between packedKey[keyStart..keyStart+len-1] and
//  hsp->packedDNA[seqStart...seqStart+len-1]. packedKeyLength is necessary for packedKey to be read correctly.
int CEPackedMismatchMatching(unsigned long long * packedKey, int packedKeyLength,
                             int keyStart, HSP * hsp, unsigned long long seqStart, int len);


//CEPackedMismatchMatching gives the number of mismatch between packedKey[keyStart..keyStart+len-1] and
//  hsp->packedDNA[seqStart...seqStart+len-1]. packedKeyLength is necessary for packedKey to be read correctly.
// ** keyStart and packedKeyLength is not used.
int CEPackedMismatchMatchingWithQuality(unsigned long long * packedKey, int packedKeyLength,
                             int keyStart, HSP * hsp, unsigned long long seqStart, int len, 
                             SRAError errors[]);

















//====================PRIMITIVE (SHOULD NOT BE USED)======================
//GetPackedPatternLong : packs the convertedKey into packedKey, each character taking only BIT_PER_CHAR.
//The function returns the number of cell the packedKey took.
//  The packedKey(the uneven key is right aligned on the last cell) is also returned; 
//  thus it should be an unsigned long long array with at least (SRA_MAX_READ_LENGTH/CHAR_PER_64) elements.
//The preliminary design : A cell = 64bit unsigned long long.
int GetPackedPatternLong(const unsigned char * convertedKey,
                         int start, int len,
                         unsigned long long * packedKey);

unsigned int PackedDifference64(unsigned long long seq,HSP *hsp,unsigned int tp, int offset, unsigned int keyLength);

unsigned int PackedDifferenceLong(unsigned long long * exSeq,HSP *hsp,unsigned int start, unsigned int len);

#endif
