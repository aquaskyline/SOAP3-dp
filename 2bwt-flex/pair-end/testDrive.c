
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include "PEAlgnmt.h"
#include "../2bwt-lib/2BWT-Interface.h"
#include "../2bwt-lib/Timing.h"

// AATTCTGTCCTGCCTTTTTCTTCCTGTATTTAAGTCTCCGGGGGCTGGGGGAACCAGGGTTTCCCA
// Positive = 267221340 267221341 (2)
// Negative = 2493137843 2493137845 (3)

// CTGACAGTCTCAGTTGCACACACGAGCCAGCAGAGGGGTTTTGTGCCACTTCTGGATGCTAGGGTT
// Positive = 1318044704 1318044704 (1)
// Negative = 132201094 132201093 (0)

int main() {

    #define DEFINE_INSERTION_SIZE 2000

    PESRAAlignmentResult * resultA = NULL;
    PESRAAlignmentResult * resultB = NULL;
    SRAOccurrence * occList_1 = NULL;
    SRAOccurrence * occList_2 = NULL;
    SRAOccurrence * sortBuffer_1 = NULL;
    SRAOccurrence * sortBuffer_2 = NULL;

    printf("Loading index ... "); 
    fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT("../../../../db/ncbi.genome.fa.index",".sa4");
    printf("DONE\n\n"); 
    
    PEInput * peInput = PEInputConstruct(idx2BWT->bwt,idx2BWT->hsp);
    
    //char d;
    //scanf("%c",&d);
    /*{
        resultA = malloc(sizeof(PESRAAlignmentResult) * 10240);
        resultB = malloc(sizeof(PESRAAlignmentResult) * 10240);
        
        int i=0;
        
        {
            resultA[i].saIndexLeft = 267221340;
            resultA[i].saIndexRight = 267221340;
            resultA[i].strand = QUERY_NEG_STRAND;
            resultA[i].mismatchCount = 0;
            
            i++;
            resultA[i].saIndexLeft = 2493137843;
            resultA[i].saIndexRight = 2493137845;
            resultA[i].strand = QUERY_NEG_STRAND;
            resultA[i].mismatchCount = 0;

            i=0;
            resultB[i].saIndexLeft = 1318044704;
            resultB[i].saIndexRight = 1318044704;
            resultB[i].strand = QUERY_POS_STRAND;
            resultB[i].mismatchCount = 0;
        }
        
        
        occList_1 = malloc(sizeof(SRAOccurrence) * 1024);
        sortBuffer_1 = malloc(sizeof(SRAOccurrence) * 1024);
        sortBuffer_2 = malloc(sizeof(SRAOccurrence) * 1024);
        
        peInput->bwt = idx2BWT->bwt;
        peInput->hsp = idx2BWT->hsp;
        peInput->insertSizeMean = DEFINE_INSERTION_SIZE;
        peInput->insertSizeStdDev = 5;
        
        PEPairList * resultList = PEMapping(peInput,resultA,10240,resultB,10240);
        
        PEPrintPEPairList(resultList);
    }*/
    
    {
        #define LIST_SIZE 102400
        int i;
        
        occList_1 = malloc(sizeof(SRAOccurrence) * LIST_SIZE);
        occList_2 = malloc(sizeof(SRAOccurrence) * LIST_SIZE);
        sortBuffer_1 = malloc(sizeof(SRAOccurrence) * LIST_SIZE);
        sortBuffer_2 = malloc(sizeof(SRAOccurrence) * LIST_SIZE);
        PEOutput * resultList = PEOutputConstruct();
        
        {
            for (i=0;i<LIST_SIZE;i++) {
                occList_1[i].ambPosition = i+DEFINE_INSERTION_SIZE;
                occList_1[i].strand = QUERY_NEG_STRAND;
                occList_1[i].mismatchCount = 0;
            }
            for (i=0;i<LIST_SIZE;i++) {
                occList_2[i].ambPosition = i;
                occList_2[i].strand = QUERY_POS_STRAND;
                occList_2[i].mismatchCount = 0;
            }
        }
        
        
        peInput->bwt = idx2BWT->bwt;
        peInput->hsp = idx2BWT->hsp;
        peInput->patternLength = 1;
        //peInput->insertSizeMean = DEFINE_INSERTION_SIZE;
        //peInput->insertSizeStdDev = 5;
        peInput->insertLbound = DEFINE_INSERTION_SIZE - 3 * 5;
        peInput->insertUbound = DEFINE_INSERTION_SIZE + 3 * 5;
        peInput->strandLeftLeg  = QUERY_POS_STRAND;
        peInput->strandRightLeg = QUERY_NEG_STRAND;
            
        //Radix sort the lists
        //PERadixSort(occList_1,LIST_SIZE,sortBuffer_1,sortBuffer_2,occList_1);
        //PERadixSort(occList_2,LIST_SIZE,sortBuffer_1,sortBuffer_2,occList_2);
        
        //Merge
        unsigned int occIndex_1 = 0;
        unsigned int occIndex_2 = 0;
        
        unsigned int occCount_1 = LIST_SIZE;
        unsigned int occCount_2 = LIST_SIZE;
        
        
        PEOutput * peOutput = PEOutputConstruct();
        
        {
            peInput->OutputType = PE_REPORT_ONE;
            PEOutputInitialise(peOutput);
            unsigned int occIndex_First = 0;
            unsigned int occIndex_Last = 0;
            double startTime;
            double elapsedTime = 0;
            startTime = setStartTime();
            PEMappingOccurrences(peInput,peOutput,occList_1,LIST_SIZE,occList_2,LIST_SIZE);
            elapsedTime = getElapsedTime(startTime);
            printf("Pair-end Alignment Time = %9.4f s\n", elapsedTime);
            unsigned int count = PECountPEOutput(peOutput);
            printf("Pair-end Alignment Found = %u\n", count);
        }
        
        printf("PE_REPORT_ONE\n");
        PEPrintPEOutputWithLimit(peOutput,10);
        
        {
            peInput->OutputType = PE_REPORT_ONE;
            PEOutputInitialise(peOutput);
            unsigned int occIndex_First = 0;
            unsigned int occIndex_Last = 0;
            double startTime;
            double elapsedTime = 0;
            startTime = setStartTime();
            PEMappingOccurrences(peInput,peOutput,occList_1,LIST_SIZE,occList_2,LIST_SIZE);
            elapsedTime = getElapsedTime(startTime);
            printf("Pair-end Alignment Time = %9.4f s\n", elapsedTime);
            unsigned int count = PECountPEOutput(peOutput);
            printf("Pair-end Alignment Found = %u\n", count);
        }
        
        printf("PE_REPORT_ONE\n");
        PEPrintPEOutputWithLimit(peOutput,10);
        
        {
            peInput->OutputType = PE_REPORT_ALL;
            PEOutputInitialise(peOutput);
            unsigned int occIndex_First = 0;
            unsigned int occIndex_Last = 0;
            double startTime;
            double elapsedTime = 0;
            startTime = setStartTime();
            PEMappingOccurrences(peInput,peOutput,occList_1,LIST_SIZE,occList_2,LIST_SIZE);
            elapsedTime = getElapsedTime(startTime);
            printf("Pair-end Alignment Time = %9.4f s\n", elapsedTime);
            unsigned int count = PECountPEOutput(peOutput);
            printf("Pair-end Alignment Found = %u\n", count);
        }
        
        printf("PE_REPORT_ALL\n");
        PEPrintPEOutputWithLimit(peOutput,10);
        
        /*{
            peInput->insertSizeMean = DEFINE_INSERTION_SIZE;
            peInput->insertSizeStdDev = 5;
            peInput->insertLbound = peInput->insertSizeMean - 3 * peInput->insertSizeStdDev;
            peInput->insertUbound = peInput->insertSizeMean + 3 * peInput->insertSizeStdDev;
                
            //Radix sort the lists
            PERadixSort(occList_1,LIST_SIZE,sortBuffer_1,sortBuffer_2,occList_1);
            PERadixSort(occList_2,LIST_SIZE,sortBuffer_1,sortBuffer_2,occList_2);
            
            PEOutputInitialise(peOutput);
            unsigned int occIndex_First = 0;
            unsigned int occIndex_Last = 0;
            double startTime;
            double elapsedTime = 0;
            startTime = setStartTime();
            
            // ATTENTION : BRUTE FORCE ENGINE
            //---------------------------------------------------------
            for (occIndex_1=0;occIndex_1<occCount_1;occIndex_1++) {
                int flag = 0;
                
                for (occIndex_2=occIndex_First;occIndex_2<occCount_2;occIndex_2++) {
                    if (PEIsPairEndMatch(peInput, &(occList_1[occIndex_1]), &(occList_2[occIndex_2]))) {
                        if (flag==0) {
                            occIndex_First = occIndex_2;
                            flag=1; 
                        }
                        //Output answer
                        PEReportPairResult(peOutput,&(occList_1[occIndex_1]),&(occList_2[occIndex_2]));
                        
                        occIndex_Last = occIndex_2;
                    }
                        
                    if (PEIsPairOutOfRange(peInput, &(occList_1[occIndex_1]), &(occList_2[occIndex_2]))) {
                        break;
                    }
                }
            }
            //---------------------------------------------------------
            
            elapsedTime = getElapsedTime(startTime);
            printf("Pair-end Alignment Time = %9.4f s\n", elapsedTime);
            unsigned int count = PECountPEOutput(peOutput);
            printf("Pair-end Alignment Found = %u\n", count);
        }*/
        
        PEOutputFree(peOutput);
    }
    
    /*{
        occList_1 = malloc(sizeof(SRAOccurrence) * 1024);
        sortBuffer_1 = malloc(sizeof(SRAOccurrence) * 1024);
        sortBuffer_2 = malloc(sizeof(SRAOccurrence) * 1024);
        
        int i = 0;
        int count=1024;
        srand(time(NULL));
        for (i=0;i<count;i++) {
            //Random
            occList_1[i].ambPosition = rand();
            
            //Reverse
            //occList_1[i].ambPosition = count-i;
        }
            
        printf("List of occurrences A =\n");
        for (i=0;i<count;i++) {
            printf("%u ", occList_1[i].ambPosition);
            if ((i+1)%8==0) {printf("\n");}
        } printf("\n");
        
        PERadixSort(occList_1,1024,sortBuffer_1,sortBuffer_2,occList_1);
        
        printf("List of occurrences A =\n");
        int flag=0;
        for (i=0;i<count;i++) {
            printf("%u ", occList_1[i].ambPosition);
            if ((i+1)%8==0) {printf("\n");}
            if (i>0 && occList_1[i-1].ambPosition>occList_1[i].ambPosition) {
                flag=1;
            }
        } printf("\n");
        if (flag==1) {printf("ERROR\n");}
    }//*/
                        
    
    if (resultA!=NULL) free(resultA);
    if (resultB!=NULL) free(resultB);
    if (occList_1!=NULL) free(occList_1);
    if (occList_2!=NULL) free(occList_2);
    if (sortBuffer_1!=NULL) free(sortBuffer_1);
    if (sortBuffer_2!=NULL) free(sortBuffer_2);
    
    // Free up the 2BWT index
    printf("\nFree index ... "); 
    fflush(stdout);
    BWTFree2BWT(idx2BWT);
    printf("DONE\n"); 
    
    return 0;

}
