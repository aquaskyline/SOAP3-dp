//
//    SRACore.c
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
   
   Date   : 30th December 2011
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#include "SRACore.h"

/////////////////////////////////////////////////////
/*
    Function SRAOccurrenceCopy
*/
/////////////////////////////////////////////////////
void SRAOccurrenceCopy(SRAOccurrence * source, SRAOccurrence * destination) {
    int i;
    destination->type                       = source->type;
    destination->ambPosition                = source->ambPosition;
    destination->strand                     = source->strand;
    destination->mismatchCount              = source->mismatchCount;
    destination->matchLen                   = source->matchLen;
    
    for (i=0;i<source->mismatchCount&&i<MAX_NUM_OF_ERROR;i++) {
        destination->errors[i].type         = source->errors[i].type;
        destination->errors[i].position     = source->errors[i].position;
    }
}


/////////////////////////////////////////////////////
/*
    Function SRARadixSort
*/
/////////////////////////////////////////////////////
void SRARadixSort(SRAOccurrence * list, unsigned long long * _count,
                 SRAOccurrence * buffer1, SRAOccurrence * buffer2,
                 SRAOccurrence * sortedList) {
                 
    
    //How many bits we sort each pass
    #define PE_RADIX_BIT_PER_PASS 4
    //How many pass in total we run the sorting
    #define PE_RADIX_PASS         8
    //The multiplying result of the two parameters above should be 
    //greater than the bit-length of the maximum index on the ref.seq.
    //e.g. 4 and 8 such that the sorting will cover 32-bit
    //     which can represent all index in the ref. seq.
    
    unsigned long long i,j;
    unsigned int base;
    int idxIndexSize = (1<<PE_RADIX_BIT_PER_PASS);
    int radixMask = idxIndexSize-1;
    unsigned int idxRadix[idxIndexSize+1];
    
    unsigned long long count = (*_count);
    
    if (count==0) return;
    
    SRAOccurrence * swap;
    SRAOccurrence * activeBuffer = buffer1;
    SRAOccurrence * oldBuffer = list;

#ifdef PE_DEBUG_PRINT_RADIX_PROGRESS
    printf("  idxIndexSize = %d\n",idxIndexSize);
    printf("  radixMask = %d\n",radixMask);
    printf("  PE_RADIX_BIT_PER_PASS = %d\n",PE_RADIX_BIT_PER_PASS);
    printf("  PE_RADIX_PASS = %d\n",PE_RADIX_PASS);
#endif
    
    for (base=0;base<PE_RADIX_PASS;base++) {
        
        //Initialise idxRadix
        for (i=0;i<=idxIndexSize;i++) {
            idxRadix[i] = 0;
        }
        
        //Count oldBuffer into activeBuffer
        for (i=0;i<count;i++) {
            unsigned int key = oldBuffer[i].ambPosition;
            key = key >> base*PE_RADIX_BIT_PER_PASS;
            key &= radixMask;
            idxRadix[key+1]++;
        }
        
        //Post-process idxRadix
        for (i=1;i<=idxIndexSize;i++) {
            idxRadix[i] += idxRadix[i-1];
        }
        
        //Move oldBuffer into activeBuffer with ordering
        for (i=0;i<count;i++) {
            unsigned int key = oldBuffer[i].ambPosition;
            key = key >> base*PE_RADIX_BIT_PER_PASS;
            key &= radixMask;
            
            SRAOccurrenceCopy(&(oldBuffer[i]),&(activeBuffer[idxRadix[key]]));
            //activeBuffer[idxRadix[key]].ambPosition = oldBuffer[i].ambPosition;
            //activeBuffer[idxRadix[key]].mismatchCount = oldBuffer[i].mismatchCount;
            //activeBuffer[idxRadix[key]].strand = oldBuffer[i].strand;
            
            idxRadix[key]++;
        }
                
        //Swap double buffer
        swap = activeBuffer;
        if (base==0) {
            activeBuffer = buffer2;
        } else {
            activeBuffer = oldBuffer;
        }
        oldBuffer = swap;
    }
    
    //Count oldBuffer into sortedlist
    unsigned long long lastOccIdx = 0;
    unsigned long long sortedListIdx = 0;
    SRAOccurrenceCopy(&(oldBuffer[0]),&(sortedList[sortedListIdx++]));
    
    for (i=1;i<count;i++) {
        if (oldBuffer[i].ambPosition == oldBuffer[lastOccIdx].ambPosition &&
            oldBuffer[i].strand == oldBuffer[lastOccIdx].strand) {
            
        } else {
            SRAOccurrenceCopy(&(oldBuffer[i]),&(sortedList[sortedListIdx++]));
            //sortedList[i].ambPosition = oldBuffer[i].ambPosition;
            //sortedList[i].mismatchCount = oldBuffer[i].mismatchCount;
            //sortedList[i].strand = oldBuffer[i].strand;
            lastOccIdx = i;
        }
    }
    
    (*_count) = sortedListIdx;
}


/////////////////////////////////////////////////////
/*
    Function SRAOccurrencesSort
*/
/////////////////////////////////////////////////////
void SRAOccurrencesSort(SRAOccurrence * occList_1, unsigned long long * _occCount_1) {
                        
                        
    //Declare variables
    unsigned long long allocMemory = 0;
    unsigned long long occCount_1 = (*_occCount_1);
    
    //These sort buffer will be used as double buffer for the sorting
    SRAOccurrence * sortBuffer_1;
    SRAOccurrence * sortBuffer_2;
    
    //Allocate enough memory for all occurrences
    #define PE_ALLOCATED_BYTE 268435456
    allocMemory = 4 * occCount_1 * sizeof(SRAOccurrence);
    if (allocMemory>PE_ALLOCATED_BYTE) {
        fprintf(stderr, "Insufficient memory to proceed. %llu required exceeded %u.\n",allocMemory,PE_ALLOCATED_BYTE);
        exit (1);
    }
    
    sortBuffer_1 = (SRAOccurrence*) malloc( occCount_1 * sizeof(SRAOccurrence) ) ;
    sortBuffer_2 = (SRAOccurrence*) malloc( occCount_1 * sizeof(SRAOccurrence) ) ;
    
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
    printf("%llu positions in list 1.\n",occCount_1);
    #endif
    
    //Radix sort the lists
    SRARadixSort(occList_1,_occCount_1,sortBuffer_1,sortBuffer_2,occList_1);
    
    free(sortBuffer_1);
    free(sortBuffer_2);
}







/////////////////////////////////////////////////////
/*
    DEBUG-Functions Below
*/
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
/*
    Function SRAOccurrencesPrint
*/
/////////////////////////////////////////////////////
void SRAOccurrencesPrint(SRAOccurrence * occList_1, unsigned long long occCount_1) {

    if (occCount_1==0) return;
    
    unsigned int i;
    char strandChar[5] = ".+-?";
    
    for (i=0;i<occCount_1;i++) {
        printf("%d\t%c(%d)\t%llu %u\n",occList_1[i].mismatchCount,strandChar[occList_1[i].strand],occList_1[i].strand,occList_1[i].ambPosition,occList_1[i].matchLen);
    }
}
