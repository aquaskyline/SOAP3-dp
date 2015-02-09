/*

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Enhancing the 2BWT library with interface functions
            for basic BWT search.

*/

#include "2BWT-Interface.h"

Idx2BWT * BWTLoad2BWT(const char * indexFilePrefix, const char * saFileNameExtension) {

    Idx2BWT * idx2BWT = (Idx2BWT*) malloc(sizeof(Idx2BWT));
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    
    char bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char bwtOccFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char saFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtOccFilename[MAX_INDEX_FILENAME_LENGTH];
    char packedDnaFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char annotationFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char ambiguityFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char translateFilename[MAX_INDEX_FILENAME_LENGTH]; 
    
    strcpy(bwtFilename,indexFilePrefix);
    strcpy(bwtOccFilename,indexFilePrefix);
    strcpy(saFilename,indexFilePrefix);
    strcpy(rev_bwtFilename,indexFilePrefix);
    strcpy(rev_bwtOccFilename,indexFilePrefix);
    strcpy(packedDnaFilename,indexFilePrefix);
    strcpy(annotationFilename,indexFilePrefix);
    strcpy(ambiguityFilename,indexFilePrefix);
    strcpy(translateFilename,indexFilePrefix);
    
    strcat(bwtFilename,".bwt");
    strcat(bwtOccFilename,".fmv");
    strcat(saFilename,saFileNameExtension);
    strcat(rev_bwtFilename,".rev.bwt");
    strcat(rev_bwtOccFilename,".rev.fmv");
    strcat(packedDnaFilename,".pac");
    strcat(annotationFilename,".ann");
    strcat(ambiguityFilename,".amb");
    strcat(translateFilename,".tra");
    
    MMMasterInitialize(3, 0, FALSE, NULL);
    MMPool * mmPool = MMPoolCreate(2097152);
    
    bwt = BWTLoad(mmPool, bwtFilename, bwtOccFilename, saFilename, NULL, NULL, NULL, 0);
    rev_bwt = BWTLoad(mmPool, rev_bwtFilename, rev_bwtOccFilename, NULL, NULL, NULL, NULL, 0);
    hsp = HSPLoad(mmPool, packedDnaFilename, annotationFilename, ambiguityFilename,translateFilename, 1, 0);
    
    HSPFillCharMap(idx2BWT->charMap);
    HSPFillComplementMap(idx2BWT->complementMap);
     
    idx2BWT->bwt = bwt;
    idx2BWT->rev_bwt = rev_bwt;
    idx2BWT->hsp = hsp;
    idx2BWT->mmPool = mmPool;
    
    return idx2BWT;
}

void BWTFree2BWT(Idx2BWT * idx2BWT) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    HSP * hsp = idx2BWT->hsp;
    MMPool * mmPool = idx2BWT->mmPool;
    
    HSPFree(mmPool, hsp, 1, 0);
    BWTFree(mmPool, bwt, 0);
    BWTFree(mmPool, rev_bwt, 0);
    MMPoolFree(mmPool);
    
    free(idx2BWT);
}

void BWTConvertPattern(Idx2BWT * idx2BWT, const char * patternSource, int patternLength, unsigned char * patternDestination) {

    int i;
    for (i=0;i<patternLength;i++) {
        patternDestination[i] = idx2BWT->charMap[patternSource[i]];
    }
    patternDestination[i]='\0';
}


void BWTSARangeInitial(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    (*saIndexLeft) = bwt->cumulativeFreq[c]+1;
    (*saIndexRight) = bwt->cumulativeFreq[c+1];

}


void BWTSARangeBackward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int l = (*saIndexLeft);
    unsigned int r = (*saIndexRight);
    (*saIndexLeft)  = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
    (*saIndexRight) = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r + 1, c);

}

void BWTSARangeBackward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = (*saIndexLeft);
    unsigned int r = (*saIndexRight);
    unsigned int rev_l = (*rev_saIndexLeft);
    unsigned int rev_r = (*rev_saIndexRight);
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }

    l = bwt->cumulativeFreq[c] + oL[c] + 1;
    r = bwt->cumulativeFreq[c] + oR[c];
    rev_r = rev_r - oCount[c];
    rev_l = rev_r - (r-l);

    (*saIndexLeft) = l;
    (*saIndexRight) = r;
    (*rev_saIndexLeft) = rev_l;
    (*rev_saIndexRight) = rev_r;
    
}
void BWTSARangeForward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = (*saIndexLeft);
    unsigned int r = (*saIndexRight);
    unsigned int rev_l = (*rev_saIndexLeft);
    unsigned int rev_r = (*rev_saIndexRight);
    int k;
    
    BWTAllOccValue(rev_bwt,rev_l,oL);
    BWTAllOccValue(rev_bwt,rev_r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
    rev_r = bwt->cumulativeFreq[c] + oR[c];
    r = r - oCount[c];
    l = r - (rev_r-rev_l);
    
    (*saIndexLeft) = l;
    (*saIndexRight) = r;
    (*rev_saIndexLeft) = rev_l;
    (*rev_saIndexRight) = rev_r;
}

void BWTAllSARangesBackward(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        unsigned int *resultSaIndexesLeft, unsigned int *resultSaIndexesRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = saIndexLeft;
    unsigned int r = saIndexRight;
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    
    for (k=0;k<ALPHABET_SIZE;k++) {
        resultSaIndexesLeft[k]  = bwt->cumulativeFreq[k] + oL[k] + 1;
        resultSaIndexesRight[k] = bwt->cumulativeFreq[k] + oR[k];
    }

}


void BWTAllSARangesBackward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        const unsigned int rev_saIndexLeft, const unsigned int rev_saIndexRight,
                        unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,
                        unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = saIndexLeft;
    unsigned int r = saIndexRight;
    unsigned int rev_l = rev_saIndexLeft;
    unsigned int rev_r = rev_saIndexRight;
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    for (k=0;k<ALPHABET_SIZE;k++) {

        resultSaIndexLeft[k] = bwt->cumulativeFreq[k] + oL[k] + 1;
        resultSaIndexRight[k] = bwt->cumulativeFreq[k] + oR[k];
        rev_resultSaIndexRight[k] = rev_r - oCount[k];
        rev_resultSaIndexLeft[k] = rev_resultSaIndexRight[k] - 
                                    (resultSaIndexRight[k]-resultSaIndexLeft[k]);

    }

}

void BWTAllSARangesForward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        const unsigned int rev_saIndexLeft, const unsigned int rev_saIndexRight,
                        unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,
                        unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = saIndexLeft;
    unsigned int r = saIndexRight;
    unsigned int rev_l = rev_saIndexLeft;
    unsigned int rev_r = rev_saIndexRight;
    int k;
    
    BWTAllOccValue(rev_bwt,rev_l,oL);
    BWTAllOccValue(rev_bwt,rev_r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    for (k=0;k<ALPHABET_SIZE;k++) {

        rev_resultSaIndexLeft[k] = bwt->cumulativeFreq[k] + oL[k] + 1;
        rev_resultSaIndexRight[k] = bwt->cumulativeFreq[k] + oR[k];
        resultSaIndexRight[k] = r - oCount[k];
        resultSaIndexLeft[k] = resultSaIndexRight[k] - 
                                    (rev_resultSaIndexRight[k]-rev_resultSaIndexLeft[k]);

    }
}

void BWTRetrievePositionFromSAIndex(Idx2BWT * idx2BWT, 
                                    unsigned int saIndex, 
                                    int * sequenceId, unsigned int * offset) {
    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    HSP * hsp = idx2BWT->hsp;
    unsigned int * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    
    unsigned int ambPosition = BWTSaValue(bwt,saIndex);
    unsigned int approxIndex = ambPosition>>GRID_SAMPLING_FACTOR_2_POWER;
    unsigned int approxValue = ambiguityMap[approxIndex];
    while (translate[approxValue].startPos>ambPosition) {
        approxValue--;
    }
    ambPosition-=translate[approxValue].correction;
    
    (*sequenceId) = translate[approxValue].chrID;
    (*offset) = ambPosition;
}

