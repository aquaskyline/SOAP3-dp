//
//    highOccBuilder.c
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

#include "HOCCConstruct.h"

HOCCBuilder * HOCCInitialHighOccWorkGround() {
    
    HOCCBuilder * hoccBuilder = (HOCCBuilder*)  malloc (sizeof(HOCCBuilder));
    
    hoccBuilder->ho_Size = 104857600;
    hoccBuilder->ho_root=(struct HOCCCellSpace*) malloc(sizeof(struct HOCCCellSpace));
    hoccBuilder->ho_root->saL = (unsigned int*) malloc(sizeof(unsigned int) * hoccBuilder->ho_Size );
    hoccBuilder->ho_root->saR = (unsigned int*) malloc(sizeof(unsigned int) * hoccBuilder->ho_Size );
    hoccBuilder->ho_root->count=0;
    hoccBuilder->ho_root->next=NULL;
    
    hoccBuilder->ho_currentRoot=hoccBuilder->ho_root;
    
    hoccBuilder->ho_ttlItem=0;
    hoccBuilder->ho_ttlOccurrence=0;
    
    return hoccBuilder;
}

void AppendHighOccPattern(HOCCBuilder * hoccBuilder, unsigned int saL, unsigned int saR) {
    if (hoccBuilder->ho_currentRoot->count >= hoccBuilder->ho_Size) {
        hoccBuilder->ho_currentRoot->next = (struct HOCCCellSpace*) malloc(sizeof(struct HOCCCellSpace));
        hoccBuilder->ho_currentRoot =hoccBuilder->ho_currentRoot->next;
        
        hoccBuilder->ho_currentRoot->saL=(unsigned int*) malloc(sizeof(unsigned int) *hoccBuilder->ho_Size);
        hoccBuilder->ho_currentRoot->saR=(unsigned int*) malloc(sizeof(unsigned int) *hoccBuilder->ho_Size);
        hoccBuilder->ho_currentRoot->count=0;
        hoccBuilder->ho_currentRoot->next=NULL;
    }
    hoccBuilder->ho_currentRoot->saL[hoccBuilder->ho_currentRoot->count]=saL;
    hoccBuilder->ho_currentRoot->saR[hoccBuilder->ho_currentRoot->count]=saR;
    hoccBuilder->ho_currentRoot->count++;
    hoccBuilder->ho_ttlItem++;
    hoccBuilder->ho_ttlOccurrence+=saR-saL+1;
}

void ho_recurFree(ho_cellSpacePtr node) {
    if (node!=NULL) {
        ho_recurFree(node->next);
        free(node);
    }
}

unsigned int CountHighOccPattern(HOCCBuilder * hoccBuilder) {
    return hoccBuilder->ho_ttlItem;
}
void HOCCFreeHighOccWorkGround(HOCCBuilder * hoccBuilder) {
    ho_recurFree(hoccBuilder->ho_root);
    free(hoccBuilder);
}

void HOCCExtractionHighOccPattern(HOCCBuilder * hoccBuilder, const BWT *bwt,const BWT * rev_bwt, 
                                int index, unsigned int HashOccThreshold,
                                unsigned int l, unsigned int r,unsigned int rev_l, unsigned int rev_r) {
    if (l<=r && r-l+1>=HashOccThreshold) {
        if (index>0) {
            //If the packing is not finished
            unsigned int new_rev_l,new_rev_r;
            unsigned int new_l,new_r;
            int numofsagroups;
            
            unsigned int occCount_pstart[ALPHABET_SIZE];
            unsigned int occCount_pend[ALPHABET_SIZE];
            unsigned int occCount[ALPHABET_SIZE];
            
            BWTAllOccValue(rev_bwt,rev_l,occCount_pstart);
            BWTAllOccValue(rev_bwt,rev_r + 1,occCount_pend);
            occCount[3]=0;
            int k;
            for (k=2;k>=0;k--) {
                occCount[k]=occCount[k+1]+occCount_pend[k+1]-occCount_pstart[k+1];
            }
            
            unsigned char ec;
            for (ec=0;ec<ALPHABET_SIZE;ec++) {
                new_rev_l = rev_bwt->cumulativeFreq[ec] + occCount_pstart[ec] + 1;
                new_rev_r = rev_bwt->cumulativeFreq[ec] + occCount_pend[ec];
                
                new_r = r - occCount[ec];
                new_l = new_r - (new_rev_r-new_rev_l);

                HOCCExtractionHighOccPattern(hoccBuilder,bwt,rev_bwt,index-1,HashOccThreshold,new_l,new_r,new_rev_l,new_rev_r);
            }
        } else if (index==0) {
            //If the packing pattern reaches the maximum length 
            //l,r
            AppendHighOccPattern(hoccBuilder,l,r);
        }
    }
}

void HOCCPopulateHighOccPatternToBitPattern(HOCCBuilder * hoccBuilder, BWT * bwt, HSP* hsp,
                                        const char * HighOccPatternFileName,const char * HighOccPackedFileName, 
                                        const char * HighOccAuxValueFileName, const char * HighOccFileName) {
    int BIT_PER_WORD = 32;
    int OCC_MAJOR_INTERVAL = 65536;
    int OCC_MINOR_INTERVAL = 256;
    unsigned int SAMPLE_FACTOR=4;
    
    
    unsigned int entryIndex=0;
    unsigned int accumulateOcc=0;
    
    //Build up the lookup table (bitTemplate).
    unsigned int i = 0;
    int j;
    unsigned int bitTemplate[32];
    bitTemplate[31]=1;
    for (j=30;j>=0;j--) {
        bitTemplate[j]=bitTemplate[j+1]<<1;
    }
    
    //Build up the bit pattern (bitPattern).
    unsigned int bitPatternLength = bwt->textLength/SAMPLE_FACTOR;
    if (bwt->textLength % SAMPLE_FACTOR==0) bitPatternLength++;
    
    unsigned int wordCount = (bitPatternLength/BIT_PER_WORD);
    if (bitPatternLength % BIT_PER_WORD != 0) wordCount++;
    
    unsigned int * bitPattern = (unsigned int*) malloc(sizeof(unsigned int) * wordCount);
    for (i=0;i<wordCount;i++) {bitPattern[i]=0;}
    
    //printf("word count = %u\n",wordCount);
    
    unsigned int occMajorCount = bitPatternLength/OCC_MAJOR_INTERVAL;
    if (bitPatternLength % OCC_MAJOR_INTERVAL !=0 ) occMajorCount++;
    unsigned int * occMajorValue = (unsigned int*) malloc(sizeof(unsigned int) * occMajorCount);
    for (i=0;i<occMajorCount;i++) {occMajorValue[i]=0;}
    
    unsigned int occMinorCount = bitPatternLength/OCC_MINOR_INTERVAL;
    if (bitPatternLength % OCC_MINOR_INTERVAL !=0 ) occMinorCount++;
    unsigned short * occMinorValue = (unsigned short*) malloc(sizeof(unsigned short) * occMinorCount);
    for (i=0;i<occMinorCount;i++) {occMinorValue[i]=0;}
    
    //Initialized the data table (dataTable)
    unsigned int highOccPatternCount = CountHighOccPattern(hoccBuilder);
    unsigned int * accTable = (unsigned int*) malloc(sizeof(unsigned int) * highOccPatternCount);
    
    //Populate the data structure, as well as in the same time..
    //Build up the auxiliary data structure for rank query
    ho_cellSpacePtr node =  hoccBuilder->ho_root;
    while (node!=NULL) {
        for (i=0;i<node->count;i++) {
            //For every node, we have saL and saR
            unsigned sampledPos = (node->saL[i])/SAMPLE_FACTOR;
            
            unsigned int offset = sampledPos % BIT_PER_WORD;
            unsigned int pos = (sampledPos-offset) / BIT_PER_WORD;
            bitPattern[pos]|=bitTemplate[offset];
            
            /*dataTable[entryIndex++]=node->saL[i];
            dataTable[entryIndex++]=node->saR[i];
            dataTable[entryIndex++]=accumulateOcc;*/
            accTable[entryIndex++]=accumulateOcc;
            accumulateOcc+=node->saR[i]-node->saL[i]+1;
            
            occMajorValue[sampledPos/OCC_MAJOR_INTERVAL]++;
            occMinorValue[sampledPos/OCC_MINOR_INTERVAL]++;
//printf("Here we have the sa range [%u,%u]\n",node->saL[i],node->saR[i]);
//printf("L falls on the %u-th unsigned int, offset = %u\n",pos,offset);
//printf("bit pattern %u is overwriting onto it.\n",bitTemplate[offset]);
//scanf("%u");
        }
        node=node->next;
    }
    
    //Second phrase to set up the auxiliary data structure for rank query
    for (i=1;i<occMajorCount;i++) {occMajorValue[i]+=occMajorValue[i-1];}
    for (i=0;i<occMinorCount;i++) {
        if (i*OCC_MINOR_INTERVAL % OCC_MAJOR_INTERVAL != 0) {
            occMinorValue[i]+=occMinorValue[i-1];
        }
    }
    
    //Third phrase to save down all the data structures
    //1. Accumulate Occ Table
    FILE * outFile;
    if(!(outFile = fopen(HighOccPatternFileName, "w"))) return;
    fwrite(&highOccPatternCount,sizeof(unsigned int),1,outFile);
    entryIndex=0;
    node =  hoccBuilder->ho_root;
    while (node!=NULL) {
        for (i=0;i<node->count;i++) {
            fwrite(&(node->saL[i]),sizeof(unsigned int),1,outFile);
            fwrite(&(accTable[entryIndex]),sizeof(unsigned int),1,outFile);
            entryIndex++;
        }
        node=node->next;
    }
    fclose(outFile);
    
    //2. Bit Pattern
    if(!(outFile = fopen(HighOccPackedFileName, "w"))) return;
    fwrite(&wordCount,sizeof(unsigned int),1,outFile);
    for (i=0;i<wordCount;i++) {
        fwrite(&(bitPattern[i]),sizeof(unsigned int),1,outFile);
    }
    fclose(outFile);
    
    
    //3. Auxiliary Data Structure
    if(!(outFile = fopen(HighOccAuxValueFileName, "w"))) return;
    fwrite(&occMajorCount,sizeof(unsigned int),1,outFile);
    fwrite(&occMinorCount,sizeof(unsigned int),1,outFile);
    for (i=0;i<occMajorCount;i++) {
        fwrite(&(occMajorValue[i]),sizeof(unsigned int),1,outFile);
    }
    for (i=0;i<occMinorCount;i++) {
        fwrite(&(occMinorValue[i]),sizeof(unsigned short),1,outFile);
    }
    fclose(outFile);
    
    
    unsigned int chromosomeSampleFactor=-1;
    unsigned int chromosomeSampleTable[150];
    unsigned int acc = 0;
    for (i=0;i<150;i++) chromosomeSampleTable[i]=0;
    for (i=0;i<hsp->numOfSeq;i++) {
        unsigned int len = hsp->seqOffset[i].endPos-hsp->seqOffset[i].startPos+1;
        if (len<chromosomeSampleFactor) {chromosomeSampleFactor=len;}
    }
    for (i=0;i<hsp->numOfSeq;i++) {
        unsigned int len = hsp->seqOffset[i].endPos-hsp->seqOffset[i].startPos+1;
        acc+=len;
        chromosomeSampleTable[acc/chromosomeSampleFactor]+=1;
    }
    for (i=1;i<150;i++) {
        chromosomeSampleTable[i]+=chromosomeSampleTable[i-1];
    }
    //4. Genuine Occ Values (Inevitable)
    if(!(outFile = fopen(HighOccFileName, "w"))) return;
    fwrite(&accumulateOcc,sizeof(unsigned int),1,outFile);
    node =  hoccBuilder->ho_root;
    while (node!=NULL) {
        for (i=0;i<node->count;i++) {
            unsigned int k=0;
            for (k=node->saL[i];k<=node->saR[i];k++) {
                unsigned int tp = BWTSaValue(bwt,k);
                //unsigned int chrNo = DetermineChromosomeRank(hsp,chromosomeSampleTable,chromosomeSampleFactor,tp);
                //unsigned int utp = TranslateUnambiguousTextPosition(hsp,tp,chrNo);
                fwrite(&tp,sizeof(unsigned int),1,outFile);
            }
        }
        node=node->next;
    }
    fclose(outFile);
    
    free(occMajorValue);
    free(occMinorValue);
    free(bitPattern);
    free(accTable);
}
