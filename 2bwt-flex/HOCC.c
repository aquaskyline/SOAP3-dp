//
//    HOCC.c
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

#include "HOCC.h"

HOCC * HOCCLoad(const char * patternFileName,  const char * occFileName,  const char * bitFileName,  const char * auxValueFileName) {
    FILE *inFile;
    HOCC* occTable;
    unsigned int i;
    unsigned int highOccPatternCount;
    unsigned int bpWordCount;
    unsigned int accOccCount;
    unsigned int majorCount;
    unsigned int minorCount;
    
    occTable=(HOCC*) malloc(sizeof(HOCC));
    occTable->table=NULL;
    occTable->occ=NULL;
    occTable->packed=NULL;
    occTable->majorValue=NULL;
    occTable->minorValue=NULL;
    occTable->ttlPattern = 0;
    occTable->ttlOcc = 0;
    
    
    if(!(inFile = fopen(patternFileName, "r"))) return NULL;
    fread((unsigned int *)&highOccPatternCount,sizeof(unsigned int),1,inFile);
    occTable->table = (unsigned int*) malloc(sizeof(unsigned int) * highOccPatternCount *2);
    for (i=0;i<highOccPatternCount *2;i++) {
        fread((unsigned int *)&(occTable->table[i]),sizeof(unsigned int),1,inFile);
    }
    fclose(inFile);
    
    if(!(inFile = fopen(occFileName, "r"))) return NULL;
    fread((unsigned int *)&accOccCount,sizeof(unsigned int),1,inFile);
    occTable->occ = (unsigned int*) malloc(sizeof(unsigned int) * accOccCount);
    for (i=0;i<accOccCount;i++) {
        fread((unsigned int *)&(occTable->occ[i]),sizeof(unsigned int),1,inFile);
    }
    fclose(inFile);
    
    if(!(inFile = fopen(bitFileName, "r"))) return NULL;
    fread((unsigned int *)&bpWordCount,sizeof(unsigned int),1,inFile);
    occTable->packed = (unsigned int*) malloc(sizeof(unsigned int) * bpWordCount);
    for (i=0;i<bpWordCount;i++) {
        fread((unsigned int *)&(occTable->packed[i]),sizeof(unsigned int),1,inFile);
    }
    fclose(inFile);
    
    if(!(inFile = fopen(auxValueFileName, "r"))) return NULL;
    fread((unsigned int *)&majorCount,sizeof(unsigned int),1,inFile);
    fread((unsigned int *)&minorCount,sizeof(unsigned int),1,inFile);
    occTable->majorValue = (unsigned int*) malloc(sizeof(unsigned int) * majorCount);
    occTable->minorValue = (unsigned short*) malloc(sizeof(unsigned short) * minorCount);
    for (i=0;i<majorCount;i++) {
        fread((unsigned int *)&(occTable->majorValue[i]),sizeof(unsigned int),1,inFile);
    }
    for (i=0;i<minorCount;i++) {
        fread((unsigned int *)&(occTable->minorValue[i]),sizeof(unsigned short),1,inFile);
    }
    fclose(inFile);
    
    occTable->ttlPattern = highOccPatternCount;
    occTable->ttlOcc = accOccCount;
    return occTable;
}

void HOCCFree(HOCC * occTable) {
    unsigned int i;
    if (occTable!=NULL) {
        if (occTable->packed!=NULL) free(occTable->packed);
        if (occTable->majorValue!=NULL) free(occTable->majorValue);
        if (occTable->minorValue!=NULL) free(occTable->minorValue);
        if (occTable->occ!=NULL) free(occTable->occ);
        if (occTable->table!=NULL) free(occTable->table);
        free(occTable);
    }
} 

void HOCCGetCachedSARange(HOCC* highOcc, unsigned long long index,unsigned long long * l,unsigned long long * r) {
    unsigned int hOccCachedL = highOcc->table[index*2];
    unsigned int hOccCachedR;
    if (index<(highOcc->ttlPattern)-1) {
        hOccCachedR = hOccCachedL+ ((highOcc->table[index*2+3])-(highOcc->table[index*2+1]))-1;
    } else { 
        hOccCachedR = hOccCachedL+ ((highOcc->ttlOcc)-(highOcc->table[index*2+1]))-1;
    }
    (*l)=hOccCachedL;
    (*r)=hOccCachedR;
}


unsigned int HOCCApproximateRank(HOCC* highOcc, unsigned long long l) {
    unsigned int * bitPattern=highOcc->packed;
    unsigned int * occMajorValue=highOcc->majorValue;
    unsigned short * occMinorValue=highOcc->minorValue;
    
    int BIT_PER_WORD = 32;
    int OCC_MAJOR_INTERVAL = 65536;
    int OCC_MINOR_INTERVAL = 256;
    unsigned int SAMPLE_FACTOR=4;
    
    unsigned int sampledIndex = l/SAMPLE_FACTOR;
    
    unsigned int majorIndex = sampledIndex / OCC_MAJOR_INTERVAL;
    unsigned int minorIndex = sampledIndex / OCC_MINOR_INTERVAL;
    
    int offset = sampledIndex % OCC_MINOR_INTERVAL;    
    unsigned int i = minorIndex * OCC_MINOR_INTERVAL;
    int wordToCount = offset / BIT_PER_WORD;
    wordToCount += (offset % BIT_PER_WORD > 0);
    
    unsigned int popCount = 0;
    int j;
    for (j=0;j<wordToCount-1;j++) {
        popCount+=__builtin_popcount(bitPattern[i/BIT_PER_WORD+j]);
        //printBinary(bitPattern[i/BIT_PER_WORD+j],BIT_PER_WORD);
    }
    unsigned int offsetBitPattern = bitPattern[i/BIT_PER_WORD+j] >> (BIT_PER_WORD - (offset % BIT_PER_WORD));
    //printBinary(bitPattern[i/BIT_PER_WORD+j],BIT_PER_WORD);
    popCount+=__builtin_popcount(offsetBitPattern);
    
    unsigned int adsCount = 0;
    if (majorIndex>0) adsCount+=occMajorValue[majorIndex-1];
    if (minorIndex>0 && minorIndex%(OCC_MAJOR_INTERVAL/OCC_MINOR_INTERVAL)>0) adsCount+=occMinorValue[minorIndex-1];
    //printf("ADS count = %u\n",adsCount);
    
    
    return adsCount+popCount;
}