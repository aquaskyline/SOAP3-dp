//
//    PEAlgnmt.c
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
/*

    Modification History 

    Date  : 22nd June 2011
    Author: Edward MK Wu
    Change: New file.

    Date  : 9th September 2011
    Author: Edward MK Wu
    Change: Free allocated memory after PEMapping.
            Invented PEMapping Mark 2 with double performance.

    Date  : 17th July 2011
    Author: Edward MK Wu
    Change: Add checking for query parameters verification.
            Add support for POS/NEG, POS/POS, NEG/NEG or NEG/POS PE Alignment.

    NOTE  : Before this file is released.
            Developer must remove all ATTENTION tokens in the content.
            These ATTENTION tokens are in place for purpose.
            Rectify them before release.

*/
/////////////////////////////////////////////////////

#include "PEAlgnmt.h"

// --------------------------------------------------
// PE_DEBUG #1
// Uncomment the below to turn on logging message for
// PE matching of two occurrences.
//#define PE_DEBUG_PRINT_MATCHING_PROGRESS
// --------------------------------------------------
// --------------------------------------------------
// PE_DEBUG #2
// Uncomment the below to turn on logging message for
// PE Radix Sort of two occ lists.
//#define PE_DEBUG_PRINT_RADIX_PROGRESS
// --------------------------------------------------
// --------------------------------------------------
// PE_DEBUG #3
// Uncomment the below to turn on visual logging for
// PE matching of two occurrences.
//#define PE_DEBUG_PRINT_VISUAL_MATCHING
// --------------------------------------------------

void PEReportPairResult(PEOutput * peOutput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int insertion);

inline int PEIsPairEndMatch(PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int * insertion);
inline int PEIsPairOutOfRange(PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2);

void PEPairListInitialise(PEPairList * pePairList);
void PEOutputInitialise(PEOutput * peOutput);
PEPairList * PEPairListConstruct();
void PEPairListFree(PEPairList * pePairList);

/////////////////////////////////////////////////////
/*
    Function PERetrievePositionFromSAIndex
    This function is a place-holder for possibly extension of
    PE to use High-Occ data-structure. However unlikely.
*/
/////////////////////////////////////////////////////
inline unsigned long long PERetrievePositionFromSAIndex(BWT * bwt, unsigned long long saIndex) {
    unsigned long long ambPosition = BWTSaValue(bwt,saIndex);
    return ambPosition;
}

/////////////////////////////////////////////////////
/*
    Function PERetrieveChromoPositioning
*/
/////////////////////////////////////////////////////
void PERetrieveChromoPositioning(BWT * bwt, 
                                HSP * hsp,
                                unsigned long long ambPosition, 
                                int * sequenceId, unsigned long long * offset) {
                                
    unsigned short * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;

    unsigned int approxIndex = ambPosition>>GRID_SAMPLING_FACTOR_2_POWER;
    unsigned int approxValue = ambiguityMap[approxIndex];
    while (translate[approxValue].startPos>ambPosition) {
        approxValue--;
    }
    ambPosition-=translate[approxValue].correction;
    
    (*sequenceId) = translate[approxValue].chrID;
    (*offset) = ambPosition;
}

inline void PEMappingCore(PEInput * peInput, PEOutput * peOutput,
                        SRAOccurrence * occList_1, unsigned long long occCount_1,
                        SRAOccurrence * occList_2, unsigned long long occCount_2) {
                        
    // ATTENTION : BRUTE FORCE ENGINE
    //---------------------------------------------------------
    /*unsigned long long occIndex_First = 0;
    unsigned long long occIndex_Last = 0;
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
    }*/
    
    // ATTENTION : BRUTE FORCE ENGINE - MARK 2
    //---------------------------------------------------------
    unsigned long i;
    unsigned long occIndex_1 = 0;
    unsigned long occIndex_2 = 0;
    uint8_t OutputType = peInput->OutputType;
    uint8_t strandLeftLeg = peInput->strandLeftLeg;
    uint8_t strandRightLeg = peInput->strandRightLeg;
    unsigned int insertion;
    
    int bestCollectorStatus = PE_BEST_NOT_AVAILABLE;
    unsigned long long bestPairOccIdx1;
    unsigned long long bestPairOccIdx2;
    unsigned int bestPairInsertion;
    
    while (occIndex_1<occCount_1 && occIndex_2<occCount_2) {
        if (occList_1[occIndex_1].ambPosition < occList_2[occIndex_2].ambPosition) {
            if (occList_1[occIndex_1].strand == strandLeftLeg) {
                for (i=occIndex_2;i<occCount_2;i++) {
                    if (PEIsPairEndMatch(peInput, &(occList_1[occIndex_1]), &(occList_2[i]),&insertion)) {
                        if (OutputType == PE_REPORT_RANDOM_BEST) {
                            if (bestCollectorStatus==PE_BEST_NOT_AVAILABLE ||
                                occList_1[bestPairOccIdx1].mismatchCount + occList_2[bestPairOccIdx2].mismatchCount >
                                occList_1[occIndex_1].mismatchCount + occList_2[i].mismatchCount) {
                                
                                bestCollectorStatus = PE_BEST_AVAILABLE;
                                bestPairOccIdx1 = occIndex_1;
                                bestPairOccIdx2 = i;
                                bestPairInsertion = insertion;
                            }
                            
                        } else if (OutputType == PE_REPORT_ALL) {
                                   
                            PEReportPairResult(peOutput,&(occList_1[occIndex_1]),&(occList_2[i]),insertion);
                            
                        } else if (OutputType == PE_REPORT_ONE) {
                        
                            PEReportPairResult(peOutput,&(occList_1[occIndex_1]),&(occList_2[i]),insertion);
                            return;
                        }
                    }
                    if (PEIsPairOutOfRange(peInput, &(occList_1[occIndex_1]), &(occList_2[i]))) {
                        break;
                    }
                }
            }
            occIndex_1++;
        } else {
            if (occList_2[occIndex_2].strand == strandLeftLeg) {
                for (i=occIndex_1;i<occCount_1;i++) {
                    if (PEIsPairEndMatch(peInput, &(occList_2[occIndex_2]), &(occList_1[i]),&insertion)) {
                        if (OutputType == PE_REPORT_RANDOM_BEST) {
                            if (bestCollectorStatus==PE_BEST_NOT_AVAILABLE ||
                                occList_1[bestPairOccIdx1].mismatchCount + occList_2[bestPairOccIdx2].mismatchCount >
                                occList_1[i].mismatchCount + occList_2[occIndex_2].mismatchCount) {
                                
                                bestCollectorStatus = PE_BEST_AVAILABLE;
                                bestPairOccIdx1 = i;
                                bestPairOccIdx2 = occIndex_2;
                                bestPairInsertion = insertion;                   
                            }
                            
                        } else if (OutputType == PE_REPORT_ALL) {
                                   
                            PEReportPairResult(peOutput, &(occList_1[i]), &(occList_2[occIndex_2]),insertion);
                        } else if (OutputType == PE_REPORT_ONE) {
                        
                            PEReportPairResult(peOutput, &(occList_1[i]), &(occList_2[occIndex_2]),insertion);
                            return;
                        }
                    }
                    if (PEIsPairOutOfRange(peInput, &(occList_2[occIndex_2]), &(occList_1[i]))) {
                        break;
                    }
                }
            }
            occIndex_2++;
        }
    }
    
    if (OutputType == PE_REPORT_RANDOM_BEST && bestCollectorStatus==PE_BEST_AVAILABLE) {
        PEReportPairResult(peOutput, &(occList_1[bestPairOccIdx1]), &(occList_2[bestPairOccIdx2]),bestPairInsertion);
    }
    //---------------------------------------------------------
}


/////////////////////////////////////////////////////
/*
    Function PEValidateAndPreparePEInput
*/
/////////////////////////////////////////////////////
inline int PEValidateAndPreparePEInput(PEInput * peInput) {
    int isValid = 1;
    int isWarning = 0;
    
    if (peInput->bwt==NULL || peInput->hsp==NULL) {
        #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("[PEValidateAndPreparePEInput] BWT/HSP index missing.\n");
        #endif
        isValid = 0;
    }
    
    if (peInput->OutputType!=PE_REPORT_ALL &&
        peInput->OutputType!=PE_REPORT_RANDOM_BEST &&
        peInput->OutputType!=PE_REPORT_ONE) {
        #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("[PEValidateAndPreparePEInput] Defaulting Output Type.\n");
        #endif
        peInput->OutputType = PE_REPORT_ALL;
        isWarning = 1;
    }
    
    if (peInput->strandLeftLeg==PE_NON_INPUT_VALUE || 
        peInput->strandRightLeg==PE_NON_INPUT_VALUE) {
        #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("[PEValidateAndPreparePEInput] Defaulting strandLeftLeg / strandRightLeg.\n");
        #endif
        peInput->strandLeftLeg = QUERY_POS_STRAND;
        peInput->strandRightLeg = QUERY_NEG_STRAND;
        isWarning = 1;
    }
    
    if ((peInput->insertSizeMean==PE_NON_INPUT_VALUE ||
        peInput->insertSizeStdDev==PE_NON_INPUT_VALUE) &&
        (peInput->insertLbound==PE_NON_INPUT_VALUE ||
        peInput->insertUbound==PE_NON_INPUT_VALUE )) {
        #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("[PEValidateAndPreparePEInput] insertion size parameter missing.\n");
        #endif
        isValid = 0;
    }
    
    if (peInput->insertLbound==PE_NON_INPUT_VALUE || 
        peInput->insertUbound==PE_NON_INPUT_VALUE) {
        
        peInput->insertLbound = peInput->insertSizeMean - 3 * peInput->insertSizeStdDev;
        peInput->insertUbound = peInput->insertSizeMean + 3 * peInput->insertSizeStdDev;
    }
    
    return isValid;
}


/////////////////////////////////////////////////////
/*
    Function PEMapping
*/
/////////////////////////////////////////////////////
void PEMapping(PEInput * peInput, PEOutput * peOutput,
                        PESRAAlignmentResult * resultListA, unsigned long long resultCountA,
                        PESRAAlignmentResult * resultListB, unsigned long long resultCountB,
                        int resultListsQualifier) {
                        
    //2BWT Index parameter
    BWT * bwt = peInput->bwt;
    HSP * hsp = peInput->hsp;
                        
    //Declare variables
    unsigned long long i,j;
    int k;
    unsigned long long allocMemory = 0;
    
    //These sort buffer will be used as double buffer for the sorting
    SRAOccurrence * sortBuffer_1;
    SRAOccurrence * sortBuffer_2;
    SRAOccurrence * occList_1;
    SRAOccurrence * occList_2;
    unsigned long long occIndex_1, occIndex_2;
    
    //Initialise the peOutput
    PEOutputInitialise(peOutput);
        
    //enrich the peInput
    // [m - 3 * sd, m + 3 * sd], and the 
    if (!PEValidateAndPreparePEInput(peInput)) {
        peOutput->flag = PE_ALIGNMENT_INPUT_ERROR;
        return;
    }
    
    //More information about the occurrence list
    unsigned long long occCount_1 = 0, occCount_2 = 0;
    unsigned long long occMaxCount = 0;
    for (i=0;i<resultCountA;i++) {
        occCount_1 += resultListA[i].saIndexRight - resultListA[i].saIndexLeft + 1;
    }
    for (i=0;i<resultCountB;i++) {
        occCount_2 += resultListB[i].saIndexRight - resultListB[i].saIndexLeft + 1;
    }
    occMaxCount = occCount_1;
    if (occMaxCount < occCount_2) {
        occMaxCount = occCount_2;
    }

    //Allocate enough memory for all occurrences
    //We need at most (4 x max(occCount_1,occCount_2) x 10 Bytes (SRAOccurrence) in total for PE mapping.
    //Assuming we have at least 256-Megabyte left in main memory
    #define PE_ALLOCATED_BYTE 268435456
    allocMemory = 4 * occMaxCount * sizeof(SRAOccurrence);
    if (allocMemory>PE_ALLOCATED_BYTE) {
        fprintf(stderr, "Insufficient memory to proceed. %llu required exceeded %u.\n",allocMemory,PE_ALLOCATED_BYTE);
        exit (1);
    }
    
    sortBuffer_1 = (SRAOccurrence*) malloc( occMaxCount * sizeof(SRAOccurrence) ) ;
    sortBuffer_2 = (SRAOccurrence*) malloc( occMaxCount * sizeof(SRAOccurrence) ) ;
    occList_1 = (SRAOccurrence*) malloc( occCount_1 * sizeof(SRAOccurrence) ) ;
    occList_2 = (SRAOccurrence*) malloc( occCount_2 * sizeof(SRAOccurrence) ) ;    
    
    //Retrieve the occurrences
    occIndex_1 = 0;
    for (i=0;i<resultCountA;i++) {
        for (j=resultListA[i].saIndexLeft;j<=resultListA[i].saIndexRight;j++) {
            occList_1[occIndex_1].ambPosition = PERetrievePositionFromSAIndex(bwt,j);
            occList_1[occIndex_1].strand = resultListA[i].strand;
            occList_1[occIndex_1].mismatchCount = resultListA[i].mismatchCount;
            
            occIndex_1++;
        }
    }
    
    occIndex_2 = 0;
    for (i=0;i<resultCountB;i++) {
        for (j=resultListB[i].saIndexLeft;j<=resultListB[i].saIndexRight;j++) {
            occList_2[occIndex_2].ambPosition = PERetrievePositionFromSAIndex(bwt,j);
            occList_2[occIndex_2].strand = resultListB[i].strand;
            occList_2[occIndex_2].mismatchCount = resultListB[i].mismatchCount;
            
            occIndex_2++;
        }
    }

    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
    printf("%llu(%llu) positions are retrieved from %llu SA ranges.\n",occIndex_1,occCount_1,resultCountA);
    printf("%llu(%llu) positions are retrieved from %llu SA ranges.\n",occIndex_2,occCount_2,resultCountB);
    #endif
    
    if (resultListsQualifier==PE_INPUTOCC_UNSORTED) {
        //Radix sort the lists
        SRARadixSort(occList_1,&occCount_1,sortBuffer_1,sortBuffer_2,occList_1);
        SRARadixSort(occList_2,&occCount_2,sortBuffer_1,sortBuffer_2,occList_2);
    }
    
    //Calling the core
    PEMappingCore(peInput,peOutput,occList_1,occCount_1,occList_2,occCount_2);
    
    peOutput->flag = PE_ALIGNMENT_COMPLETED;
    
    free(sortBuffer_1);
    free(sortBuffer_2);
    free(occList_1);
    free(occList_2);
}


/////////////////////////////////////////////////////
/*
    Function PEMappingOccurrences
*/
/////////////////////////////////////////////////////
void PEMappingOccurrences(PEInput * peInput, PEOutput * peOutput,
                        SRAOccurrence * occList_1, unsigned long long occCount_1,
                        SRAOccurrence * occList_2, unsigned long long occCount_2,
                        int occListsQualifier) {
                        
    //2BWT Index parameter
    BWT * bwt = peInput->bwt;
    HSP * hsp = peInput->hsp;
                        
    //Declare variables
    unsigned long long i,j;
    int k;
    unsigned long long allocMemory = 0;
    
    //These sort buffer will be used as double buffer for the sorting
    SRAOccurrence * sortBuffer_1;
    SRAOccurrence * sortBuffer_2;
    unsigned long long occIndex_1, occIndex_2;
    
    //Initialise the peOutput
    PEOutputInitialise(peOutput);
    
    //enrich the peInput
    // [m - 3 * sd, m + 3 * sd], and the 
    if (!PEValidateAndPreparePEInput(peInput)) {
        peOutput->flag = PE_ALIGNMENT_INPUT_ERROR;
        return;
    }
    
    
    if (occListsQualifier==PE_INPUTOCC_UNSORTED) {
        //More information about the occurrence list
        unsigned long long occMaxCount = 0;
        occMaxCount = occCount_1;
        if (occMaxCount < occCount_2) {
            occMaxCount = occCount_2;
        }
        
        //Allocate enough memory for all occurrences
        //We need at most (4 x max(occCount_1,occCount_2) x 10 Bytes (SRAOccurrence) in total for PE mapping.
        //Assuming we have at least 256-Megabyte left in main memory
        #define PE_ALLOCATED_BYTE 268435456
        allocMemory = 4 * occMaxCount * sizeof(SRAOccurrence);
        if (allocMemory>PE_ALLOCATED_BYTE) {
            fprintf(stderr, "Insufficient memory to proceed. %llu required exceeded %u.\n",allocMemory,PE_ALLOCATED_BYTE);
            exit (1);
        }
        
        sortBuffer_1 = (SRAOccurrence*) malloc( occMaxCount * sizeof(SRAOccurrence) ) ;
        sortBuffer_2 = (SRAOccurrence*) malloc( occMaxCount * sizeof(SRAOccurrence) ) ;
        
        #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("%llu positions in list 1.\n",occCount_1);
        printf("%llu positions in list 2.\n",occCount_2);
        #endif
        
        //Radix sort the lists
        SRARadixSort(occList_1,&occCount_1,sortBuffer_1,sortBuffer_2,occList_1);
        SRARadixSort(occList_2,&occCount_2,sortBuffer_1,sortBuffer_2,occList_2);
        
        free(sortBuffer_1);
        free(sortBuffer_2);
    }
    
    //Calling the core
    PEMappingCore(peInput,peOutput,occList_1,occCount_1,occList_2,occCount_2);
    
    peOutput->flag = PE_ALIGNMENT_COMPLETED;
    
}

/////////////////////////////////////////////////////
/*
    Function PEIsPairEndMatch
    This function takes in 2 occurrence and determind
    if the two occurrences form a pair end mapping.
    Assume occ_1 <= occ_2.
    
    Input:
        peInput - The input parameter for PE
        occ_1 - The occurrence on left leg
        occ_2 - The occurrence on right leg
        
    Output:
        1 if the two occurrences form a pair end mapping
        0 otherwise
        
*/
/////////////////////////////////////////////////////
int PEIsPairEndMatch(PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int * insertion) {
    //For a valid paired-end reads,
    // let the mean and the standard deviation of the insert size be "m" and "sd",
    // the distance between the leftmost of the left-read alignment and the rightmost 
    // of the right-read alignment should fall between [m - 3 * sd, m + 3 * sd], and the 
    // strand of the left-read alignment is NOT the same as the strand of the right-read alignment.
    
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("PEIsPairEndMatch is invoked.\n");
        printf("    pos = %llu \t strand = %d \n",occ_1->ambPosition,occ_1->strand);
        printf("    pos = %llu \t strand = %d \n",occ_2->ambPosition,occ_2->strand);
    #endif
    
    int strandLeftLeg = peInput->strandLeftLeg;
    int strandRightLeg = peInput->strandRightLeg;
    
    int strandCheck = (occ_1->strand == strandLeftLeg &&
                        occ_2->strand == strandRightLeg);
    
    //positions are left aligned, therefore it requires correction
    unsigned long long ambPos_2 = occ_2->ambPosition + occ_2->matchLen - 1;

    unsigned long long gap = ambPos_2 - occ_1->ambPosition + 1;
    (*insertion) = gap;
    
    #ifdef PE_DEBUG_PRINT_VISUAL_MATCHING
    char DEBUG_sorted[2];
    DEBUG_sorted[0]='X';DEBUG_sorted[1]=' ';
    printf(" %10llu <----- %10llu -----> %10llu    [%d:%d] %c\n", occ_1->ambPosition, gap, ambPos_2, peInput->insertLbound,peInput->insertUbound,DEBUG_sorted[(occ_1->ambPosition<=ambPos_2)&1]);
    #endif
    
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("    lbound = %u \t gap = %llu \t ubound = %u\n",peInput->insertLbound,gap,peInput->insertUbound);
    #endif
    
    return ((peInput->insertLbound <= gap) && (gap <= peInput->insertUbound) && strandCheck);
}

int PEDPIsPairEndMatch(PEInput * peInput, DPOccurrence * occ_1, DPOccurrence * occ_2, unsigned int * insertion) {
    //For a valid paired-end reads,
    // let the mean and the standard deviation of the insert size be "m" and "sd",
    // the distance between the leftmost of the left-read alignment and the rightmost 
    // of the right-read alignment should fall between [m - 3 * sd, m + 3 * sd], and the 
    // strand of the left-read alignment is NOT the same as the strand of the right-read alignment.
    
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("PEDPIsPairEndMatch is invoked.\n");
        printf("    pos = %llu \t strand = %d \n",occ_1->ambPosition,occ_1->strand);
        printf("    pos = %llu \t strand = %d \n",occ_2->ambPosition,occ_2->strand);
    #endif
    
    int strandLeftLeg = peInput->strandLeftLeg;
    int strandRightLeg = peInput->strandRightLeg;
    
    int strandCheck = (occ_1->strand == strandLeftLeg &&
                        occ_2->strand == strandRightLeg);
    
    //positions are left aligned, therefore it requires correction
    unsigned long long ambPos_2 = occ_2->ambPosition + occ_2->matchLen - 1;

    unsigned long long gap = ambPos_2 - occ_1->ambPosition + 1;
    (*insertion) = gap;
    
    #ifdef PE_DEBUG_PRINT_VISUAL_MATCHING
    char DEBUG_sorted[2];
    DEBUG_sorted[0]='X';DEBUG_sorted[1]=' ';
    printf(" %10llu <----- %10llu -----> %10llu    [%d:%d] %c\n", occ_1->ambPosition, gap, ambPos_2, peInput->insertLbound,peInput->insertUbound,DEBUG_sorted[(occ_1->ambPosition<=ambPos_2)&1]);
    #endif
    
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("    lbound = %u \t gap = %llu \t ubound = %u\n",peInput->insertLbound,gap,peInput->insertUbound);
        printf("    lstrand = %d rstrand= %d\n",strandLeftLeg,strandRightLeg);
    #endif
    
    return ((peInput->insertLbound <= gap) && (gap <= peInput->insertUbound) && strandCheck);
}

inline int PEIsPairOutOfRange(PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2) {
    if (occ_1->strand == occ_2->strand) { return 0; }
    
    unsigned long long occ_1_pos = occ_1->ambPosition + peInput->insertUbound;
    unsigned long long occ_2_pos = occ_2->ambPosition + occ_2->matchLen - 1;
    
    return occ_1_pos<occ_2_pos;
}

void PEReportPairResult(PEOutput * peOutput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int insertion) {
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("    Pair-end alignment is found.\n");
    #endif
    
    PEPairList * pePairList = peOutput->tail;

    if (pePairList->pairsCount == PE_MAX_BUCKET_SIZE) {
        if (pePairList->next == NULL) {
            pePairList->next = PEPairListConstruct();
        }
        pePairList = pePairList->next;
        peOutput->tail = pePairList;
    }
    
    PEPairs * pePair = &(pePairList->pairs[pePairList->pairsCount]);
    pePair->occ_1 = occ_1;
    pePair->occ_2 = occ_2;
    pePair->insertion = insertion;
    
    pePairList->pairsCount++;
}

/////////////////////////////////////////////////////
// Constructor and Destructor
/////////////////////////////////////////////////////
PEInput * PEInputConstruct(BWT * bwt, HSP * hsp) {
    PEInput * peInput = (PEInput*) malloc(sizeof(PEInput));
    
    peInput->bwt = bwt;
    peInput->hsp = hsp;
    
    peInput->strandLeftLeg = (uint8_t) PE_NON_INPUT_VALUE;
    peInput->strandRightLeg = (uint8_t) PE_NON_INPUT_VALUE;
    
    peInput->insertSizeMean = PE_NON_INPUT_VALUE;
    peInput->insertSizeStdDev = PE_NON_INPUT_VALUE;
    peInput->insertLbound = PE_NON_INPUT_VALUE;
    peInput->insertUbound = PE_NON_INPUT_VALUE;
    
    peInput->OutputType = PE_REPORT_ALL;
    
    return peInput;
}

void PEInputFree(PEInput * peInput) {
    free(peInput);
}

void PEPairListInitialise(PEPairList * pePairList) {
    if (pePairList==NULL) return;
    PEPairListInitialise(pePairList->next);
    pePairList->pairsCount = 0;
}

void PEOutputInitialise(PEOutput * peOutput) {
    PEPairListInitialise(peOutput->root);
    peOutput->tail = peOutput->root;
    peOutput->flag = PE_ALIGNMENT_INITIALISED;
}

PEPairList * PEPairListConstruct() {
    PEPairList * pePairList = (PEPairList*) malloc(sizeof(PEPairList));
    pePairList->pairsCount = 0;
    pePairList->next = NULL;
    return pePairList;
}

PEOutput * PEOutputConstruct() {
    PEOutput * peOutput = (PEOutput*) malloc(sizeof(PEOutput));
    peOutput->root = PEPairListConstruct();
    peOutput->tail = peOutput->root;
    peOutput->flag = PE_ALIGNMENT_INITIALISED;
    
    return peOutput;
}

void PEPairListFree(PEPairList * pePairList) {
    if (pePairList==NULL) return;
    PEPairListFree(pePairList->next);
    free(pePairList);
}

void PEOutputFree(PEOutput * peOutput) {
    PEPairListFree(peOutput->root);
    free(peOutput);
}

int PEAligned(PEOutput * peOutput) {
    if (peOutput->root==NULL) return 0;
    if (peOutput->root->pairsCount>0) return 1;
    return 0;
}

/////////////////////////////////////////////////////
// Utilities
/////////////////////////////////////////////////////
void PEPrintOccurrenceList(const char * listName, SRAOccurrence * list, unsigned long long length) {
    unsigned long long i;
    printf("List of occurrences %s =\n",listName);
    for (i=0;i<length;i++) {
        printf("%llu ", list[i].ambPosition);
        if ((i+1)%8==0) {printf("\n");}
    } printf("\n");
}
void PEPrintPEPairList(PEPairList * pairList, int maxOutput) {
    int i;
    int j = 0;
    char strand[3];
    strand[1]='+';
    strand[2]='-';
    while (pairList!=NULL && pairList->pairsCount>0) {
        for (i=0;i<pairList->pairsCount;i++) {
            PEPairs * pePair = &(pairList->pairs[i]);
            SRAOccurrence * occ_1 = pePair->occ_1;
            SRAOccurrence * occ_2 = pePair->occ_2;
            printf("algnmt_1 = %llu(%c)      algnmt_2 = %llu(%c)\n",
                    occ_1->ambPosition,strand[occ_1->strand],
                    occ_2->ambPosition,strand[occ_2->strand]);
            j++;
            if (j>=maxOutput && maxOutput>=0) return;
        }
        pairList = pairList->next;
    }
}
void PEPrintPEOutput(PEOutput * peOutput) {
    PEPrintPEPairList(peOutput->root, -1);
}
void PEPrintPEOutputWithLimit(PEOutput * peOutput, int maxOutput) {
    PEPrintPEPairList(peOutput->root, maxOutput);
}
unsigned long long PECountPEPairList(PEPairList * pairList) {
    unsigned long long count = 0;
    while (pairList!=NULL && pairList->pairsCount>0) {
        count += pairList->pairsCount;
        pairList = pairList->next;
    }
    return count;
}
unsigned long long PECountPEOutput(PEOutput * peOutput) {
    return PECountPEPairList(peOutput->root);
}


/*

WRONG BECAUSE OF STRAND
inline int PEIsPairEndMatch(PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2) {
    //For a valid paired-end reads,
    // let the mean and the standard deviation of the insert size be "m" and "sd",
    // the distance between the leftmost of the left-read alignment and the rightmost 
    // of the right-read alignment should fall between [m - 3 * sd, m + 3 * sd], and the 
    // strand of the left-read alignment is NOT the same as the strand of the right-read alignment.
    
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("PEIsPairEndMatch is invoked.\n");
        printf("    pos = %u \t strand = %d \n",occ_1->ambPosition,occ_1->strand);
        printf("    pos = %u \t strand = %d \n",occ_2->ambPosition,occ_2->strand);
    #endif
        
    if (occ_1->strand == occ_2->strand) { return 0; }
    
    unsigned int occ_1_pos = occ_1->ambPosition;
    unsigned int occ_2_pos = occ_2->ambPosition;
        
    unsigned int r = (occ_1_pos - occ_2_pos) & -(occ_1_pos < occ_2_pos);
    unsigned int minIndex = occ_2_pos + r;
    unsigned int maxIndex = occ_1_pos - r;
    
    unsigned int gap = maxIndex - minIndex + 1;
    #ifdef PE_DEBUG_PRINT_MATCHING_PROGRESS
        printf("    lbound = %u \t gap = %u \t ubound = %u\n",peInput->insertLbound,gap,peInput->insertUbound);
    #endif
    
    return ((peInput->insertLbound <= gap) && (gap <= peInput->insertUbound));
}

*/
