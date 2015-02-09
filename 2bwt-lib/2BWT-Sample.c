/*

   2BWT-Sample.c         SAMPLE CODE FOR 2BWT-LIB
   
   Copyright (C) 2011, Edward Wu.
   Website: http://www.cs.hku.hk/~mkewu/2bwt
   
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
   
   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Sample file to use 2BWT index.
   
            This program performs the following with 2BWT library
            1. BWT backward search on a given pattern
            2. BWT forward search on a given pattern
            3. BWT bi-directional search on a given pattern
            3. BWT 1-mismatch search on a given pattern
            4. Report the first occurrences of (1).
   
            For complete guide on what the 2BWT index library provided
            please check out the file 2BWT-Interface.h for the description
            of all library calls.

*/

#include <stdio.h>
#include <stdlib.h>
#include "2BWT-Interface.h"

int main() {
    int i,j,k,c;

    //Variables for backward and forward search
    unsigned int l,r,rev_l,rev_r;
    
    //Variables for search all sa ranges functions
    unsigned int result_l[ALPHABET_SIZE];
    unsigned int result_r[ALPHABET_SIZE];
    unsigned int result_rev_l[ALPHABET_SIZE];
    unsigned int result_rev_r[ALPHABET_SIZE];
    
    //Variables for result
    unsigned int offset;
    int sequenceId;
    unsigned int saCount;
    
    //Variables for pattern
    char pattern[1024];
    strcpy(pattern,"AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA");
    int patternLength = strlen(pattern);

    // Load up the index with the below statement
    printf("Loading index ... "); 
    fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT("ncbi.genome.fa.index",".sa");
    printf("DONE\n\n"); 
    
    // Convert the pattern into 2BWT recognised coding scheme
    unsigned char packedPattern[1024];
    BWTConvertPattern(idx2BWT,pattern,patternLength,packedPattern);
    
    
    
    
    
    
// The following performs a backward search of the pattern
// ===================================================================================
// |
    printf("Performing backward search of the pattern..\n");
    BWTSARangeInitial(idx2BWT,packedPattern[patternLength-1],&l,&r);
    for (i=patternLength-2;i>=0;i--) {
        BWTSARangeBackward(idx2BWT,packedPattern[i],&l,&r);
    }
    printf("SA Range being = %u %u (%u)\n\n",l,r,r-l+1);
// |
// ===================================================================================








    
// The following performs a forward search of the pattern
// ===================================================================================
// |
    printf("Performing forward search of the pattern..\n");
    BWTSARangeInitial(idx2BWT,packedPattern[0],&l,&r);
    BWTSARangeInitial(idx2BWT,packedPattern[0],&rev_l,&rev_r);
    for (i=1;i<patternLength;i++) {
        BWTSARangeForward_Bidirection(idx2BWT,packedPattern[i],&l,&r,&rev_l,&rev_r);
    }
    printf("SA Range being = %u %u %u %u (%u)\n\n",l,r,rev_l,rev_r,r-l+1);
// |
// ===================================================================================
    
    
// The following performs a bi-directional search of the pattern
// Starting from the middle of the pattern, first move right, then move left.
// ===================================================================================
// |
    printf("Performing bi-directional search of the pattern..\n");
    j = patternLength / 2;
    BWTSARangeInitial(idx2BWT,packedPattern[j],&l,&r);
    BWTSARangeInitial(idx2BWT,packedPattern[j],&rev_l,&rev_r);
    for (i=j+1;i<patternLength;i++) {
        BWTSARangeForward_Bidirection(idx2BWT,packedPattern[i],&l,&r,&rev_l,&rev_r);
    }
    for (i=j-1;i>=0;i--) {
        BWTSARangeBackward_Bidirection(idx2BWT,packedPattern[i],&l,&r,&rev_l,&rev_r);
    }
    printf("SA Range being = %u %u %u %u (%u)\n\n",l,r,rev_l,rev_r,r-l+1);
// |
// ===================================================================================
    
    
// The following performs a 1-mismatch search of the pattern
// ===================================================================================
// |
// |
    printf("Performing 1-mismatch search of the pattern..\n");
    saCount = 0;
    j = patternLength / 2;
    BWTSARangeInitial(idx2BWT,packedPattern[patternLength-1],&l,&r);
    for (i=patternLength-2;i>j-1;i--) { BWTSARangeBackward(idx2BWT,packedPattern[i],&l,&r); }
    
    for (i=j-1;i>=0;i--) {
        BWTAllSARangesBackward(idx2BWT,l,r,result_l,result_r);
        for (c=0;c<ALPHABET_SIZE;c++) {
            if (c==packedPattern[i]) continue;
            unsigned int err_l=result_l[c];
            unsigned int err_r=result_r[c];
            for (k=i-1;k>=0;k--) {
                if (err_l>err_r) break;
                BWTSARangeBackward(idx2BWT,packedPattern[k],&err_l,&err_r);
            }
            if (err_l<=err_r && k<0) {
                //An SA range of occurrence is found (err_l,err_r)
                saCount+=err_r-err_l+1;
            }
        }
        l=result_l[packedPattern[i]];
        r=result_r[packedPattern[i]];
    }
    
    BWTSARangeInitial(idx2BWT,packedPattern[0],&l,&r);
    BWTSARangeInitial(idx2BWT,packedPattern[0],&rev_l,&rev_r);
    for (i=1;i<j;i++) { BWTSARangeForward_Bidirection(idx2BWT,packedPattern[i],&l,&r,&rev_l,&rev_r); }
    for (i=j;i<patternLength;i++) {
        BWTAllSARangesForward_Bidirection(idx2BWT,l,r,rev_l,rev_r,result_l,result_r,result_rev_l,result_rev_r);
        for (c=0;c<ALPHABET_SIZE;c++) {
            if (c==packedPattern[i]) continue;
            unsigned int err_l=result_l[c];
            unsigned int err_r=result_r[c];
            unsigned int rev_err_l=result_rev_l[c];
            unsigned int rev_err_r=result_rev_r[c];
            for (k=i+1;k<patternLength;k++) {
                if (err_l>err_r) break;
                BWTSARangeForward_Bidirection(idx2BWT,packedPattern[k],&err_l,&err_r,&rev_err_l,&rev_err_r);
            }
            if (err_l<=err_r && k>=patternLength) {
                //An SA range of occurrence is found (err_l,err_r)
                saCount+=err_r-err_l+1;
            }
        }
        l=result_l[packedPattern[i]];
        r=result_r[packedPattern[i]];
        rev_l=result_rev_l[packedPattern[i]];
        rev_r=result_rev_r[packedPattern[i]];
    }
    printf("%u SA-indexes/occurrences were found.\n\n",saCount);
// |
// |
// ===================================================================================
    
    
    
    
    
    
// The following output the first 5 position of the pattern
// ===================================================================================
// |
// |
    j=(r-l+1<5)?r-l+1:5;
    printf("Reporting %d arbitrary occurrences..\n",j);
    for (i=0;i<j;i++) {
        BWTRetrievePositionFromSAIndex(idx2BWT,l+i,&sequenceId,&offset);
        printf("Occurrence found in sequence #%d with offset %u\n",sequenceId,offset);
    }
// |
// |
// ===================================================================================
    
    
    
    
    
    
    // Free up the 2BWT index
    printf("\nFree index ... "); 
    fflush(stdout);
    BWTFree2BWT(idx2BWT);
    printf("DONE\n"); 
    
    return 0;
}