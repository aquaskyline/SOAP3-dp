//
//    SRA2BWTCheckAndExtend.c
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

*/
/////////////////////////////////////////////////////

#include "SRA2BWTCheckAndExtend.h"

void CEDebugPrintPackedSequence(unsigned long long * packedKey,int len) {
    int bufferCount = len / CHAR_PER_64 + (len % CHAR_PER_64>0);
    int i,k,j;

    j=0;
    for (i=0;i<bufferCount;i++) {
        unsigned long long packed64 = packedKey[i];
        for (k=CHAR_PER_64-1;k>=0;k--) {
            printf("%c",dnaChar[(packed64>>k*2) & 3]);
            j++;
            if (j % 4 == 0 ) {printf (" ");}
        }
        printf(" ");
    }
    printf("\n");
}

//GetPackedPatternLong : packs the convertedKey into packedKey, each character taking only BIT_PER_CHAR.
//The function returns the number of cell the packedKey took.
//  The packedKey(the uneven key is right aligned on the last cell) is also returned; 
//  thus it should be an unsigned long long array with at least (SRA_MAX_READ_LENGTH/CHAR_PER_64) elements.
//The preliminary design : A cell = 64bit unsigned long long.
int CEPackPattern(const unsigned char * convertedKey,
                         int start, int len,
                         unsigned long long * packedKey) {

    if (len>SRA_MAX_READ_LENGTH) {
        fprintf(stderr,"[GetPackedPatternLong] Read length longer than the maximum short read length defined.\n");
        return 0;
    }
    if (len==0) return 0;

    int pos = start;
    int i = 0;
    int bufferCount=0;
    
    packedKey[bufferCount]=convertedKey[pos++];

    for (i=1;i<len;i++) {
        if (i % CHAR_PER_64==0) {
            bufferCount++;
            packedKey[bufferCount]=convertedKey[pos];
        } else {
            packedKey[bufferCount]<<=BIT_PER_CHAR;
            packedKey[bufferCount]|=convertedKey[pos];
        }
        pos++;
    }
    return bufferCount+1;
//The pattern bit patterns are aligned to the right if the buffer is unfilled.
}

//CEPackedMismatchMatching gives the number of mismatch between packedKey[keyStart..keyStart+len-1] and
//  hsp->packedDNA[seqStart...seqStart+len-1]. packedKeyLength is necessary for packedKey to be read correctly.
int CEPackedMismatchMatching(unsigned long long * packedKey, int packedKeyLength,
                             int keyStart, HSP * hsp, unsigned long long seqStart, int len) {
    //int bufferCount = packedKeyLength / CHAR_PER_64 + (packedKeyLength % CHAR_PER_64>0);
    int matchBufferCount = len / CHAR_PER_64 + (len % CHAR_PER_64>0);

    int i = 0;
    int offset = 0;
    unsigned long long seqPos = seqStart;
    unsigned long long seqIndex;
    int seqShift = seqPos % CHAR_PER_WORD;
    unsigned int * packedDNA = hsp->packedDNA;
    int mismatchInserted = 0;

    //Bit masks
    unsigned long long packedSeqExtract;
    unsigned long long diffBitVector = 0;


    for (i=0;i<matchBufferCount-1;i++) {
        seqIndex = seqPos/CHAR_PER_WORD;
        packedSeqExtract=packedDNA[seqIndex];
        packedSeqExtract=(packedSeqExtract<<SRA_CE_BIT_PER_WORD) | packedDNA[seqIndex+1];
        if (seqShift>0) {
            packedSeqExtract <<= seqShift * BIT_PER_CHAR;
            packedSeqExtract |= packedDNA[seqIndex+2] >> ((CHAR_PER_WORD-seqShift)*BIT_PER_CHAR);
        }
        diffBitVector = packedSeqExtract ^ packedKey[i];
        diffBitVector = (diffBitVector | (diffBitVector >> 1)) & SRA_CE_BIT_MASK;
        mismatchInserted += __builtin_popcountll(diffBitVector);
        seqPos+=CHAR_PER_64;
    }

    seqIndex = seqPos/CHAR_PER_WORD;
    packedSeqExtract=packedDNA[seqIndex];
    packedSeqExtract=(packedSeqExtract<<SRA_CE_BIT_PER_WORD) | packedDNA[seqIndex+1];
    if (seqShift>0) {
        packedSeqExtract <<= seqShift * BIT_PER_CHAR;
        packedSeqExtract |= packedDNA[seqIndex+2] >> ((CHAR_PER_WORD-seqShift)*BIT_PER_CHAR);
    }
    packedSeqExtract >>= ( SRA_CE_BIT_PER_64 - ((len % CHAR_PER_64) * BIT_PER_CHAR) ) % SRA_CE_BIT_PER_64;
    diffBitVector = packedSeqExtract ^ packedKey[i];
    diffBitVector = (diffBitVector | (diffBitVector >> 1)) & SRA_CE_BIT_MASK;
    mismatchInserted += __builtin_popcountll(diffBitVector);

    return mismatchInserted;
}


//CEPackedMismatchMatching gives the number of mismatch between packedKey[keyStart..keyStart+len-1] and
//  hsp->packedDNA[seqStart...seqStart+len-1]. packedKeyLength is necessary for packedKey to be read correctly.
// ** keyStart and packedKeyLength is not used.
int CEPackedMismatchMatchingWithQuality(unsigned long long * packedKey, int packedKeyLength,
                             int keyStart, HSP * hsp, unsigned long long seqStart, int len, 
                             SRAError errors[]) {
    //int bufferCount = packedKeyLength / CHAR_PER_64 + (packedKeyLength % CHAR_PER_64>0);
    int matchBufferCount = len / CHAR_PER_64 + (len % CHAR_PER_64>0);

    int i = 0;
    int k = 0;
    int offset = 0;
    unsigned long long seqPos = seqStart;
    unsigned long long seqIndex;
    int seqShift = seqPos % CHAR_PER_WORD;
    unsigned int * packedDNA = hsp->packedDNA;
    int mismatchInserted = 0;
    int mismatch = 0;
    int leadingZero = 0;
    int shiftedAway = 0;
    int shiftedAwayLocal = 0;
    int mismatchOcc = 0;

    //Bit masks
    unsigned long long packedSeqExtract;
    unsigned long long countMask=0;
    unsigned long long diffBitVector = 0;

    //printf("CEPackedMismatchMatchingWithQuality %llu %d\n",seqStart,len);

    for (i=0;i<matchBufferCount-1;i++) {
        seqIndex = seqPos/CHAR_PER_WORD;
        packedSeqExtract=packedDNA[seqIndex];
        packedSeqExtract=(packedSeqExtract<<SRA_CE_BIT_PER_WORD) | packedDNA[seqIndex+1];
        if (seqShift>0) {
            packedSeqExtract <<= seqShift * BIT_PER_CHAR;
            packedSeqExtract |= packedDNA[seqIndex+2] >> ((CHAR_PER_WORD-seqShift)*BIT_PER_CHAR);
        }
        //CEDebugPrintPackedSequence(&packedSeqExtract,CHAR_PER_64);
        //CEDebugPrintPackedSequence(&(packedKey[i]),CHAR_PER_64);
        diffBitVector = packedSeqExtract ^ packedKey[i];
        diffBitVector = (diffBitVector | (diffBitVector >> 1)) & SRA_CE_BIT_MASK;
        mismatch = __builtin_popcountll(diffBitVector);

        countMask = diffBitVector;
        shiftedAwayLocal = shiftedAway;
        for (k=0;k<mismatch;k++) {
            leadingZero = __builtin_clzll(countMask);
            mismatchOcc = (leadingZero-BIT_PER_CHAR+1)/BIT_PER_CHAR;
            //printf("Mismatch at position = %d on read.\n",shiftedAwayLocal+mismatchOcc);
            errors[mismatchInserted+k].type = SRA_CHAR_ERROR_TYPE_MISMATCH;
            errors[mismatchInserted+k].position = shiftedAwayLocal+mismatchOcc;
            countMask <<= leadingZero+1;
            shiftedAwayLocal += mismatchOcc+1;
        }

        seqPos+=CHAR_PER_64;
        mismatchInserted+=mismatch;
        shiftedAway+=CHAR_PER_64;
    }

    seqIndex = seqPos/CHAR_PER_WORD;
    packedSeqExtract=packedDNA[seqIndex];
    packedSeqExtract=(packedSeqExtract<<SRA_CE_BIT_PER_WORD) | packedDNA[seqIndex+1];
    if (seqShift>0) {
        packedSeqExtract <<= seqShift * BIT_PER_CHAR;
        packedSeqExtract |= packedDNA[seqIndex+2] >> ((CHAR_PER_WORD-seqShift)*BIT_PER_CHAR);
    }
    packedSeqExtract >>= ( SRA_CE_BIT_PER_64 - ((len % CHAR_PER_64) * BIT_PER_CHAR) ) % SRA_CE_BIT_PER_64;

    //CEDebugPrintPackedSequence(&packedSeqExtract,CHAR_PER_64);
    //CEDebugPrintPackedSequence(&(packedKey[i]),CHAR_PER_64);

    diffBitVector = packedSeqExtract ^ packedKey[i];
    diffBitVector = (diffBitVector | (diffBitVector >> 1)) & SRA_CE_BIT_MASK;
    mismatch = __builtin_popcountll(diffBitVector);

    countMask = diffBitVector << (SRA_CE_BIT_PER_64 - ((len % CHAR_PER_64) * BIT_PER_CHAR));
    shiftedAwayLocal = shiftedAway;
    for (k=0;k<mismatch;k++) {
        leadingZero = __builtin_clzll(countMask);
        mismatchOcc = (leadingZero-BIT_PER_CHAR+1)/BIT_PER_CHAR;
        //printf("Mismatch at position = %d on read.\n",shiftedAwayLocal+mismatchOcc);
        errors[mismatchInserted+k].type = SRA_CHAR_ERROR_TYPE_MISMATCH;
        errors[mismatchInserted+k].position = shiftedAwayLocal+mismatchOcc;
        countMask <<= leadingZero+1;
        shiftedAwayLocal += mismatchOcc+1;
    }

    mismatchInserted+=mismatch;
    return mismatchInserted;
}

















//====================PRIMITIVE (SHOULD NOT BE USED)======================
//GetPackedPatternLong : packs the convertedKey into packedKey, each character taking only BIT_PER_CHAR.
//The function returns the number of cell the packedKey took.
//  The packedKey(the uneven key is right aligned on the last cell) is also returned; 
//  thus it should be an unsigned long long array with at least (SRA_MAX_READ_LENGTH/CHAR_PER_64) elements.
//The preliminary design : A cell = 64bit unsigned long long.
int GetPackedPatternLong(const unsigned char * convertedKey,
                         int start, int len,
                         unsigned long long * packedKey) {

    if (len>SRA_MAX_READ_LENGTH) {
        fprintf(stderr,"[GetPackedPatternLong] Read length longer than the maximum short read length defined.\n");
        return 0;
    }
    if (len==0) return 0;

    int pos = start;
    int i = 0;
    int bufferCount=0;
    
    packedKey[bufferCount]=convertedKey[pos++];

    for (i=1;i<len;i++) {
        if (i % CHAR_PER_64==0) {
            bufferCount++;
            packedKey[bufferCount]=convertedKey[pos];
        } else {
            packedKey[bufferCount]<<=BIT_PER_CHAR;
            packedKey[bufferCount]|=convertedKey[pos];
        }
        pos++;
    }
    return bufferCount;
//The pattern bit patterns are aligned to the right if the buffer is unfilled.
}

unsigned int PackedDifference64(unsigned long long seq,HSP *hsp,unsigned int tp, int offset, unsigned int keyLength) {
    
    //printf("TP:%u OFFSET:%d\n",tp,offset);
    //unsigned int textPosition = BWTSaValue(bwt,saIndex);
    unsigned int pos = (tp + offset);
    //printf("POS:%u\n",pos);
    unsigned int * packedDNA = hsp->packedDNA;
    if (pos+keyLength > hsp->dnaLength || pos<0) {
        printf("Invalid pos: %u\n", pos);
    }
    unsigned long long verify = packedDNA[pos / 16];
    verify = (verify << 32) | packedDNA[pos / 16 + 1];
    
    if (pos%16>0) {
        verify <<= (pos % 16) * 2;
        verify |= ((packedDNA[pos / 16 + 2])>>(16-(pos % 16))*2);
    }
    verify >>= 64 - (keyLength * 2);
    unsigned long long diff = verify ^ seq;
    diff = (diff | (diff >> 1)) & 0x5555555555555555ull;
    return __builtin_popcountll(diff);
}

unsigned int PackedDifferenceLong(unsigned long long * exSeq,HSP *hsp,unsigned int start, unsigned int len) {
    unsigned int bufferCount;
    unsigned int offset=len % CHAR_PER_64;
    
    if (offset ==0) {
        bufferCount = len/CHAR_PER_64;
    } else {
        bufferCount = len/CHAR_PER_64 + 1;
    }
    
    unsigned int flag=0;
    int k=0;
    while (k<bufferCount-1) {
        flag+=PackedDifference64(exSeq[k],hsp,start,CHAR_PER_64*k,CHAR_PER_64);
        k++;
    }
    flag+=PackedDifference64(exSeq[k],hsp,start,CHAR_PER_64*k,offset);
    
    return flag;
}

