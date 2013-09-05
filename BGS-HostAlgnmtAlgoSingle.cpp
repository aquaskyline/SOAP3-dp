/*
 *
 *    BGS-HostAlgnmtAlgoSingle.c
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#include "BGS-HostAlgnmtAlgoSingle.h"




//====================MODELIZED BELOW======================
// BWTxxxModelxxx functions are fully generalised BWT search algorithm that searchs for reads contain any number of edit/mismatch
// However, calling these functions requires user to define themselves a 'searching model'.
// The searching model requires each BWT step to be defined. The search algorithm will then follow the defined model.

// BWTMismatchModelAnyDirection_CE matches steps with check and extend mechanism.
// It allows starting off CE in the middle of a step and recursive CE until SRACase completes.
unsigned long long BWTMismatchModelAnyDirection_CE3 ( SRAQueryInput * qInput, int i, int mismatchInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occMismatch, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    // This function should not be invoked because the constant MAX_CE_THRESHOLD is set to 0
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    int OutputType = qSetting->OutputType;
    HSP * hsp = aIndex->hsp;
    BWT * ceBwt = aIndex->bwt;
    int * leftMostAligned = alignmentCase->leftMostAligned;
    int * occMismatches = rOutput->Shared_AccMismatch;
    unsigned long long * saPositions = rOutput->Shared_ReadPositions;
    char * occQualities = rOutput->Shared_AccQualities;
    unsigned long long extensionSeq[MAX_CE_BUFFER_SIZE];
    unsigned long long extensionSeq2[MAX_CE_BUFFER_SIZE];
    SRAError errorsTmp[SRA_MAX_READ_LENGTH];
    int occNBMismatch[MAX_CE_THRESHOLD];
    unsigned int readPos;
    int mismatch;
    uint8_t mismatchCe;
    uint8_t minMismatch;
    uint8_t maxMismatch;
    int step, len, start, end, pos, newPosCount, swap, left, right, posCount;
    int bufferCount;
    int k = 0;
    unsigned long long saCount = 0;
    step = alignmentCase->steps[stepInCase].step;
    start = alignmentCase->steps[stepInCase].start;
    pos = start + step * i;
    int leftAlign = leftMostAligned[pos];
    unsigned long long j;
    int param_i = i;
    int param_stepInCase = stepInCase;
    int flag;
    int errorsIdx;

    for ( j = l; j <= r; j++ )
    {
        occMismatches[k] = occMismatch;
        occNBMismatch[k] = nbMismatch;
        occQualities[k] = occQuality;
        saPositions[k] = ( *ceBwt->_bwtSaValue ) ( ceBwt, j ) - leftAlign;
        k++;
    }

    posCount = k;
    newPosCount = 0;

    while ( alignmentCase->steps[stepInCase].type != SRA_STEP_TYPE_COMPLETE )
    {
        minMismatch = alignmentCase->steps[stepInCase].MinError;
        maxMismatch = alignmentCase->steps[stepInCase].MaxError;
        step = alignmentCase->steps[stepInCase].step;
        len = alignmentCase->steps[stepInCase].len;
        start = alignmentCase->steps[stepInCase].start;
        end = alignmentCase->steps[stepInCase].end;
        pos = start + step * i;
        swap = ( pos - end ) & - ( pos < end );
        left = end + swap;
        right = pos - swap;
        len = right - left + 1;
        //bufferCount = GetPackedPatternLong(convertedKey,left,len,extensionSeq);
        bufferCount = CEPackPattern ( convertedKey, left, len, extensionSeq );

        //printf("> ");
        //for (j=0;j<len;j++) {
        //  printf("%c",dnaChar[convertedKey[left+j]]);
        //  if ((j+1) % 4==0) {printf(" ");}
        //}printf("\n");

        for ( j = 0; j < posCount; j++ )
        {
            readPos = saPositions[j];

            //mismatch = PackedDifferenceLong(extensionSeq,hsp,readPos+left,len);
            //mismatch = CEPackedMismatchMatching(extensionSeq,len,0,hsp,readPos+left,len);
            if ( readPos + left > ceBwt->textLength )
            {
                mismatchCe = MAX_NUM_OF_ERROR + 1;
            }
            else
            {
                mismatchCe = CEPackedMismatchMatchingWithQuality ( extensionSeq, len, 0, hsp, readPos + left, len, errorsTmp );
            }

            //int mismatch2 = CEPackedMismatchMatching(extensionSeq,len,0,hsp,readPos+left,len);
            //if (mismatch != mismatch2) {
            //  printf("Return number of mismatch is different : %d %d\n",mismatch,mismatch2);
            //}
            mismatch = mismatchCe + mismatchInserted;

            if ( mismatch >= minMismatch && mismatch <= ( maxMismatch + qSetting->MaxNBMismatch - occNBMismatch[j] ) )
            {
                occNBMismatch[newPosCount] = occNBMismatch[j] + ( ( mismatch - maxMismatch ) * ( mismatch > maxMismatch ) );

                occQualities[newPosCount] = occQualities[j];

                for ( errorsIdx = 0; errorsIdx < mismatchCe; errorsIdx++ )
                {
                    occQualities[newPosCount] += keyQuality[errorsTmp[errorsIdx].position];
                }

                occMismatches[newPosCount] = occMismatches[j] + mismatch;
                saPositions[newPosCount] = readPos;
                newPosCount++;
            }
        }

        posCount = newPosCount;
        newPosCount = 0;
        stepInCase++;
        i = 0;
        mismatchInserted = 0;
    }

    // collect the answers to occ_list
    for ( j = 0; j < posCount; j++ )
    {
        addOccToAlgnResult ( algnResult, saPositions[j], qInfo->ReadStrand );
    }

    return posCount;
    // return OCCReportTextPositions(qInput, posCount);
}
/*unsigned long long BWTMismatchModelAnyDirection_CE(SRAQueryInput * qInput, int i, int mismatchInserted,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges,
                                    int occMismatch, uint8_t nbMismatch, char occQuality) {

    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    unsigned char * convertedKey = qInfo->ReadCode;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    int OutputType = qSetting->OutputType;
    HSP * hsp = aIndex->hsp;
    BWT * ceBwt = aIndex->bwt;

    int * leftMostAligned = alignmentCase->leftMostAligned;
    int * occMismatches = rOutput->Shared_AccMismatch;
    unsigned long long * saPositions = rOutput->Shared_ReadPositions;
    char * occQualities = rOutput->Shared_AccQualities;
    unsigned long long extensionSeq[4];
    unsigned int readPos;
    int mismatch;
    uint8_t minMismatch;
    uint8_t maxMismatch;
    int step,len,start,end,pos,newPosCount,swap,left,right,posCount;
    int bufferCount;

    int k=0;
    unsigned long long saCount = 0;
    step = alignmentCase->steps[stepInCase].step;
    start = alignmentCase->steps[stepInCase].start;
    pos = start + step * i;
    int leftAlign = leftMostAligned[pos];
    unsigned long long j;
    int param_i = i;
    int param_stepInCase = stepInCase;
    int flag;

    for (j=l;j<=r;j++) {
        occMismatches[0]=occMismatch;
        occQualities[0]=occQuality;
        saPositions[0]=(*ceBwt->_bwtSaValue)(ceBwt,j)-leftAlign;

        i = param_i;
        stepInCase = param_stepInCase;
        flag = 1;
        mismatch = mismatchInserted;

        while (alignmentCase->steps[stepInCase].type != SRA_STEP_TYPE_COMPLETE) {

            minMismatch = alignmentCase->steps[stepInCase].MinError;
            maxMismatch = alignmentCase->steps[stepInCase].MaxError;
            step = alignmentCase->steps[stepInCase].step;
            len = alignmentCase->steps[stepInCase].len;
            start = alignmentCase->steps[stepInCase].start;
            end = alignmentCase->steps[stepInCase].end;

            pos = start + step * i;

            swap = (pos - end) & -(pos < end);
            left = end + swap;
            right = pos - swap;
            len = right - left + 1;

            bufferCount = GetPackedPatternLong(convertedKey,left,len,extensionSeq);
            mismatch += PackedDifferenceLong(extensionSeq,hsp,saPositions[0]+left,len);
            if (mismatch>=minMismatch && mismatch<=maxMismatch) {
                flag=1;
                occMismatches[0]+=mismatch;
            } else {
                flag=0;
                break;
            }

            stepInCase++;mismatch=0;
            i=0;
        }

        if (flag==1) {
            saCount+=OCCReportTextPositions(qInput, 1);

            if (rOutput->IsClosed==1) {
                return saCount;
            }
        }
    }
    return saCount;
}//*/



// BWTExactModelForward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
unsigned long long BWTExactModelForward_Lookup3 ( SRAQueryInput * qInput,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges, SingleAlgnResult * algnResult )
{
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;
    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long rev_l = 0;
    unsigned long long rev_r = 0;
    unsigned long long j;
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    LT * lookupTable = aIndex->lookupTable;
    LT * rev_lookupTable = aIndex->rev_lookupTable;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int step = alignmentCase->steps[stepInCase].step;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int pos = start;
    int k;
    int branchCount = 0;
    char branchChar = 0;
    uint8_t nbMismatch = 0;

    unsigned int lookupLength = ( lookupTable->tableSize > len ) ? len : lookupTable->tableSize;
    int i;

    for ( i = 0; i < lookupLength ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= ( convertedKey[pos++]  & LOOKUP_CHAR_MASK );
    }

    r_packedPattern = l_packedPattern;

    for ( i = lookupLength; i < lookupTable->tableSize ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK ;
    }

    l = l_packedPattern ? lookupTable->table[l_packedPattern - 1] + 1 : 1;
    r = lookupTable->table[r_packedPattern];
    l_packedPattern = 0;
    r_packedPattern = 0;

    for ( i = 0; i < lookupLength ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= ( convertedKey[--pos] & LOOKUP_CHAR_MASK );
    }

    r_packedPattern = l_packedPattern;

    for ( i = lookupLength; i < rev_lookupTable->tableSize ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK ;
    }

    rev_l = l_packedPattern ? rev_lookupTable->table[l_packedPattern - 1] + 1 : 1;
    rev_r = rev_lookupTable->table[r_packedPattern];

    if ( r - l < rev_r - rev_l )
    {
        rev_l = rev_r - ( r - l );
    }
    else if ( r - l > rev_r - rev_l )
    {
        l = r - ( rev_r - rev_l );
    }

    i = lookupLength;
    pos = start + lookupLength;

    while ( i < len && l <= r )
    {
        if ( r - l < ceThreshold && pos >= ceStart && pos <= ceEnd )
        {
            saRanges[0] = l;
            saRanges[1] = r;
            return BWTMismatchModelAnyDirection_CE3 ( qInput, i, 0,
                    alignmentCase, stepInCase,
                    saRanges, 0, nbMismatch, 0, algnResult );
        }

        unsigned char c = convertedKey[pos];
        BWTAllOccValue ( bwt, rev_l, oL );
        BWTAllOccValue ( bwt, rev_r + 1, oR );
        oCount[ALPHABET_SIZE - 1] = 0;
        branchCount = ( oR[ALPHABET_SIZE - 1] > oL[ALPHABET_SIZE - 1] );
        branchChar = ( ALPHABET_SIZE - 1 ) *  branchCount;

        for ( k = ALPHABET_SIZE - 2; k >= 0; k-- )
        {
            oCount[k] = oCount[k + 1] + oR[k + 1] - oL[k + 1];

            branchCount += ( oR[k] > oL[k] );
            branchChar |= k * ( oR[k] > oL[k] );
        }

        // Non-branching Mismatch Handler
        if ( oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < qSetting->MaxNBMismatch )
        {
            rev_l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            rev_r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            r = r - oCount[branchChar];
            l = r - ( rev_r - rev_l );
            nbMismatch++;
        }
        else
        {
            rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
            rev_r = bwt->cumulativeFreq[c] + oR[c];
            r = r - oCount[c];
            l = r - ( rev_r - rev_l );
        }

        i++;
        pos++;
    }

    if ( l <= r )
    {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saRanges[2] = rev_l;
        saRanges[3] = rev_r;

        if ( alignmentCase->steps[stepInCase + 1].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
        {
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchBackward\n");
            return BWTModelSwitchBackward3 ( qInput, 0, 0,
                                             alignmentCase, stepInCase + 1,
                                             saRanges, 0, nbMismatch, 0, algnResult );
        }
        else
        {
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchAnyDirection\n");
            return BWTModelSwitchAnyDirection3 ( qInput, 0, 0,
                                                 alignmentCase, stepInCase + 1,
                                                 saRanges, 0, nbMismatch, 0, algnResult );
        }
    }

    return 0;
}


// BWTExactModelBackward_Lookup lookup your pattern in lookup table, single direction and assuming backward
unsigned long long BWTExactModelBackward_Lookup3 ( SRAQueryInput * qInput,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges, SingleAlgnResult * algnResult )
{
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;
    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long j = 0;
    int k = 0;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    LT * lookupTable = aIndex->lookupTable;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int step = alignmentCase->steps[stepInCase].step;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---

    unsigned int lookupLength = ( lookupTable->tableSize > len ) ? len : lookupTable->tableSize;
    int branchCount = 0;
    unsigned char branchChar = 0;
    uint8_t nbMismatch = 0;

    int i;
    int pos = start - lookupLength + 1;

    for ( i = 0; i < lookupLength ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= ( convertedKey[pos++]  & LOOKUP_CHAR_MASK );
    }

    r_packedPattern = l_packedPattern;

    for ( i = lookupLength; i < lookupTable->tableSize ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK ;
    }

    l = l_packedPattern ? lookupTable->table[l_packedPattern - 1] + 1 : 1;
    r = lookupTable->table[r_packedPattern];
    i = lookupLength;
    pos = start - lookupLength;

    while ( i < len && l <= r )
    {
        if ( r - l < ceThreshold && pos >= ceStart && pos <= ceEnd )
        {
            saRanges[0] = l;
            saRanges[1] = r;
            return BWTMismatchModelAnyDirection_CE3 ( qInput, i, 0,
                    alignmentCase, stepInCase,
                    saRanges, 0, nbMismatch, 0, algnResult );
        }

        unsigned char c = convertedKey[pos];
        BWTAllOccValue ( bwt, l, oL );
        BWTAllOccValue ( bwt, r + 1, oR );

        branchCount = 0;
        branchChar = 0;

        for ( k = ALPHABET_SIZE - 1; k >= 0; k-- )
        {
            branchCount += ( oR[k] > oL[k] );
            branchChar |= k * ( oR[k] > oL[k] );
        }

        // Non-branching Mismatch Handler
        if ( oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < qSetting->MaxNBMismatch )
        {
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            nbMismatch++;
        }
        else
        {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
        }

        i++;
        pos--;
    }

    if ( l <= r )
    {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        return BWTModelSwitchBackward3 ( qInput, 0, 0,
                                         alignmentCase, stepInCase + 1,
                                         saRanges, 0, nbMismatch, 0, algnResult );
    }

    return 0;
}

// BWTExactModelBackwardAnyDirection_Lookup lookup your pattern in lookup table, bi-directional and assuming backward
unsigned long long BWTExactModelBackwardAnyDirection_Lookup3 ( SRAQueryInput * qInput,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges, SingleAlgnResult * algnResult )
{
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;

    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long rev_l = 0;
    unsigned long long rev_r = 0;
    unsigned long long j;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---

    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    LT * lookupTable = aIndex->lookupTable;
    LT * rev_lookupTable = aIndex->rev_lookupTable;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int step = alignmentCase->steps[stepInCase].step;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int k;
    int branchCount = 0;
    char branchChar = 0;
    uint8_t nbMismatch = 0;

    unsigned int lookupLength = ( lookupTable->tableSize > len ) ? len : lookupTable->tableSize;


    int i;
    int pos = start - lookupLength + 1;

    for ( i = 0; i < lookupLength ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= ( convertedKey[pos++]  & LOOKUP_CHAR_MASK );
    }

    r_packedPattern = l_packedPattern;

    for ( i = lookupLength; i < lookupTable->tableSize ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK ;
    }

    l = l_packedPattern ? lookupTable->table[l_packedPattern - 1] + 1 : 1;
    r = lookupTable->table[r_packedPattern];
    l_packedPattern = 0;
    r_packedPattern = 0;

    for ( i = 0; i < lookupLength ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= ( convertedKey[--pos] & LOOKUP_CHAR_MASK );
    }

    r_packedPattern = l_packedPattern;

    for ( i = lookupLength; i < rev_lookupTable->tableSize ; i++ )
    {
        l_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern <<= LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK ;
    }

    rev_l = l_packedPattern ? rev_lookupTable->table[l_packedPattern - 1] + 1 : 1;
    rev_r = rev_lookupTable->table[r_packedPattern];

    if ( r - l < rev_r - rev_l )
    {
        rev_l = rev_r - ( r - l );
    }
    else if ( r - l > rev_r - rev_l )
    {
        l = r - ( rev_r - rev_l );
    }

    i = lookupLength;
    pos = start - lookupLength;

    while ( i < len && l <= r )
    {
        if ( r - l < ceThreshold && pos >= ceStart && pos <= ceEnd )
        {
            saRanges[0] = l;
            saRanges[1] = r;
            return BWTMismatchModelAnyDirection_CE3 ( qInput, i, 0,
                    alignmentCase, stepInCase,
                    saRanges, 0, nbMismatch, 0, algnResult );
        }

        unsigned char c = convertedKey[pos];
        BWTAllOccValue ( bwt, l, oL );
        BWTAllOccValue ( bwt, r + 1, oR );
        oCount[ALPHABET_SIZE - 1] = 0;
        branchCount = ( oR[ALPHABET_SIZE - 1] > oL[ALPHABET_SIZE - 1] );
        branchChar = ( ALPHABET_SIZE - 1 ) *  branchCount;

        for ( k = ALPHABET_SIZE - 2; k >= 0; k-- )
        {
            oCount[k] = oCount[k + 1] + oR[k + 1] - oL[k + 1];

            branchCount += ( oR[k] > oL[k] );
            branchChar |= k * ( oR[k] > oL[k] );
        }

        // Non-branching Mismatch Handler
        if ( oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < qSetting->MaxNBMismatch )
        {
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            rev_r = rev_r - oCount[branchChar];
            rev_l = rev_r - ( r - l );
            nbMismatch++;
        }
        else
        {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
            rev_r = rev_r - oCount[c];
            rev_l = rev_r - ( r - l );
        }

        i++;
        pos--;
    }

    if ( l <= r )
    {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saRanges[2] = rev_l;
        saRanges[3] = rev_r;

        if ( alignmentCase->steps[stepInCase + 1].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT )
        {
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchBackward\n");
            return BWTModelSwitchBackward3 ( qInput, 0, 0,
                                             alignmentCase, stepInCase + 1,
                                             saRanges, 0, nbMismatch, 0, algnResult );
        }
        else
        {
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchAnyDirection\n");
            return BWTModelSwitchAnyDirection3 ( qInput, 0, 0,
                                                 alignmentCase, stepInCase + 1,
                                                 saRanges, 0, nbMismatch, 0, algnResult );
        }
    }

    return 0;
}


// BWTExactModelBackward matches pattern on text without using any other aux, e.g. lookup table.
unsigned long long BWTExactModelBackward3 ( SRAQueryInput * qInput, int i, int errorInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occError, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int step = -1;
    int pos, k;
    unsigned long long j;
    pos = start + step * i;
    step = -1;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---

    int branchCount = 0;
    char branchChar = 0;

    //printf("[BWTExactModelBackward] from %d to %d allowing %u<%u<%u mismatches.\n", start,end, minError,errorInserted, maxError);
    //printf("[BWTExactModelBackward] %d/%d\n", i,len);

    while ( i < len && l <= r )
    {
        unsigned char c = convertedKey[pos];

        if ( r - l < ceThreshold  && pos >= ceStart && pos <= ceEnd )
        {
            saRanges[0] = l;
            saRanges[1] = r;
            return BWTMismatchModelAnyDirection_CE3 ( qInput, i, errorInserted,
                    alignmentCase, stepInCase,
                    saRanges, occError, nbMismatch, occQuality, algnResult );
        }

        BWTAllOccValue ( bwt, l, oL );
        BWTAllOccValue ( bwt, r + 1, oR );
        branchCount = 0;
        branchChar = 0;

        for ( k = ALPHABET_SIZE - 1; k >= 0; k-- )
        {
            branchCount += ( oR[k] > oL[k] );
            branchChar |= k * ( oR[k] > oL[k] );
        }

        // Non-branching Mismatch Handler
        if ( oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < qSetting->MaxNBMismatch )
        {
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            nbMismatch++;
        }
        else
        {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
        }

        i++;
        pos += step;
    }

    if ( errorInserted >= minError && l <= r )
    {
        saRanges[0] = l;
        saRanges[1] = r;
        return BWTModelSwitchBackward3 ( qInput, 0, 0,
                                         alignmentCase, stepInCase + 1,
                                         saRanges, occError + errorInserted, nbMismatch, occQuality, algnResult );
    }

    return 0;
}

// BWTMismatchModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTExactModelAnyDirection3 ( SRAQueryInput * qInput, int i, int errorInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occError, uint8_t nbMismatch,  char occQuality, SingleAlgnResult * algnResult )
{
    int step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[ ( step < 0 ) * 2];
    unsigned long long r = saRanges[ ( step < 0 ) * 2 + 1];
    unsigned long long rev_l = saRanges[ ( step > 0 ) * 2];
    unsigned long long rev_r = saRanges[ ( step > 0 ) * 2 + 1];
    unsigned long long saCount = 0;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    int k, pos;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    unsigned long long err_rev_l;
    unsigned long long err_rev_r;
    unsigned long long err_r;
    unsigned long long err_l;
    unsigned long long errSaRange[4];
    unsigned char c;
    unsigned long long j;
    int branchCount = 0;
    char branchChar = 0;

    //printf("Alignment from %d to %d allowing %u-%u mismatches.\n", start,end, minError, maxError);
    //printf("[BWTExactModelAnyDirection] from %d to %d allowing %u<%u<%u mismatches.\n", start,end, minError,errorInserted, maxError);
    //printf("[BWTExactModelAnyDirection] %d/%d\n", i,len);

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---

    pos = start + step * i;
    len -= ( minError - errorInserted - 1 ) * ( minError > errorInserted + 1 );

    while ( i < len && l <= r )
    {
        if ( r - l < ceThreshold && pos >= ceStart && pos <= ceEnd )
        {
            saRanges[ ( step < 0 ) * 2] = l;
            saRanges[ ( step < 0 ) * 2 + 1] = r;
            saRanges[ ( step > 0 ) * 2] = rev_l;
            saRanges[ ( step > 0 ) * 2 + 1] = rev_r;
            k = 0;
            return BWTMismatchModelAnyDirection_CE3 ( qInput, i, errorInserted,
                    alignmentCase, stepInCase,
                    saRanges, occError, nbMismatch, occQuality, algnResult );
        }

        c = convertedKey[pos];
        BWTAllOccValue ( bwt, rev_l, oL );
        BWTAllOccValue ( bwt, rev_r + 1, oR );
        oCount[ALPHABET_SIZE - 1] = 0;
        branchCount = ( oR[ALPHABET_SIZE - 1] > oL[ALPHABET_SIZE - 1] );
        branchChar = ( ALPHABET_SIZE - 1 ) *  branchCount;

        for ( k = ALPHABET_SIZE - 2; k >= 0; k-- )
        {
            oCount[k] = oCount[k + 1] + oR[k + 1] - oL[k + 1];
            branchCount += ( oR[k] > oL[k] );
            branchChar |= k * ( oR[k] > oL[k] );
        }


        // Non-branching Mismatch Handler
        if ( oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < qSetting->MaxNBMismatch )
        {
            rev_l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            rev_r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            r = r - oCount[branchChar];
            l = r - ( rev_r - rev_l );
            nbMismatch++;
        }
        else
        {
            rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
            rev_r = bwt->cumulativeFreq[c] + oR[c];
            r = r - oCount[c];
            l = r - ( rev_r - rev_l );
        }

        i += 1;
        pos += step;
    }

    if ( errorInserted >= minError && l <= r )
    {
        //Next Step
        saRanges[ ( step < 0 ) * 2] = l;
        saRanges[ ( step < 0 ) * 2 + 1] = r;
        saRanges[ ( step > 0 ) * 2] = rev_l;
        saRanges[ ( step > 0 ) * 2 + 1] = rev_r;
        return BWTModelSwitchAnyDirection3 ( qInput, 0, 0,
                                             alignmentCase, stepInCase + 1,
                                             saRanges, occError + errorInserted, nbMismatch, occQuality, algnResult );
    }

    return 0;
}

// BWTMismatchModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTMismatchModelAnyDirection3 ( SRAQueryInput * qInput, int i, int mismatchInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occMismatch, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    int step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[ ( step < 0 ) * 2];
    unsigned long long r = saRanges[ ( step < 0 ) * 2 + 1];
    unsigned long long rev_l = saRanges[ ( step > 0 ) * 2];
    unsigned long long rev_r = saRanges[ ( step > 0 ) * 2 + 1];
    unsigned long long saCount = 0;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    int k, pos;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    uint8_t minMismatch = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxMismatch = alignmentCase->steps[stepInCase].MaxError;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    unsigned long long err_rev_l;
    unsigned long long err_rev_r;
    unsigned long long err_r;
    unsigned long long err_l;
    unsigned long long errSaRange[4];
    unsigned long long j;
    unsigned char c;
    int branchCount = 0;
    char branchChar = 0;
    //printf("Alignment from %d to %d allowing %u-%u mismatches.\n", start,end, minMismatch, maxMismatch);
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    pos = start + step * i;
    len -= ( minMismatch - mismatchInserted - 1 ) * ( minMismatch > mismatchInserted + 1 );

    while ( i < len && l <= r )
    {
        //printf("Current progress: %d(%d) within %d characters. stepping by %d.\n", i, pos, len, step);
        if ( r - l < ceThreshold && pos >= ceStart && pos <= ceEnd )
        {
            saRanges[ ( step < 0 ) * 2] = l;
            saRanges[ ( step < 0 ) * 2 + 1] = r;
            saRanges[ ( step > 0 ) * 2] = rev_l;
            saRanges[ ( step > 0 ) * 2 + 1] = rev_r;
            saCount += BWTMismatchModelAnyDirection_CE3 ( qInput, i, mismatchInserted,
                       alignmentCase, stepInCase,
                       saRanges, occMismatch, nbMismatch, occQuality, algnResult );
            return saCount;
        }

        c = convertedKey[pos];
        BWTAllOccValue ( bwt, rev_l, oL );
        BWTAllOccValue ( bwt, rev_r + 1, oR );
        oCount[ALPHABET_SIZE - 1] = 0;
        branchCount = ( oR[ALPHABET_SIZE - 1] > oL[ALPHABET_SIZE - 1] );
        branchChar = ( ALPHABET_SIZE - 1 ) *  branchCount;

        for ( k = ALPHABET_SIZE - 2; k >= 0; k-- )
        {
            oCount[k] = oCount[k + 1] + oR[k + 1] - oL[k + 1];
            branchCount += ( oR[k] > oL[k] );
            branchChar |= k * ( oR[k] > oL[k] );
        }

        //if (mismatchInserted<maxMismatch) {
        for ( k = 0; k < ALPHABET_SIZE; k++ )
        {
            if ( k != c )
            {
                err_rev_l = bwt->cumulativeFreq[k] + oL[k] + 1;
                err_rev_r = bwt->cumulativeFreq[k] + oR[k];
                err_r = r - oCount[k];
                err_l = err_r - ( err_rev_r - err_rev_l );

                if ( err_l <= err_r )
                {
                    errSaRange[ ( step < 0 ) * 2] = err_l;
                    errSaRange[ ( step < 0 ) * 2 + 1] = err_r;
                    errSaRange[ ( step > 0 ) * 2] = err_rev_l;
                    errSaRange[ ( step > 0 ) * 2 + 1] = err_rev_r;
                    saCount += BWTModelSwitchAnyDirection3 ( qInput, i + 1, mismatchInserted + 1,
                               alignmentCase, stepInCase,
                               errSaRange, occMismatch, nbMismatch,
                               occQuality + keyQuality[pos], algnResult );

                    if ( rOutput->IsClosed == 1 ) { return saCount; }
                }
            }
        }

        //}
        // Non-branching Mismatch Handler
        if ( oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < qSetting->MaxNBMismatch )
        {
            rev_l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            rev_r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            r = r - oCount[branchChar];
            l = r - ( rev_r - rev_l );
            nbMismatch++;
        }
        else
        {
            rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
            rev_r = bwt->cumulativeFreq[c] + oR[c];
            r = r - oCount[c];
            l = r - ( rev_r - rev_l );
        }

        i += 1;
        pos += step;
    }

    if ( mismatchInserted >= minMismatch && l <= r )
    {
        //Next Step
        saRanges[ ( step < 0 ) * 2] = l;
        saRanges[ ( step < 0 ) * 2 + 1] = r;
        saRanges[ ( step > 0 ) * 2] = rev_l;
        saRanges[ ( step > 0 ) * 2 + 1] = rev_r;
        saCount += BWTModelSwitchAnyDirection3 ( qInput, 0, 0,
                   alignmentCase, stepInCase + 1,
                   saRanges, occMismatch + mismatchInserted, nbMismatch, occQuality, algnResult );
    }

    return saCount;
}

// BWTMismatchModelBackward matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTMismatchModelBackward3 ( SRAQueryInput * qInput, int i, int mismatchInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occMismatch, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    unsigned long long saCount = 0;
    int k, pos;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    uint8_t minMismatch = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxMismatch = alignmentCase->steps[stepInCase].MaxError;
    int step = -1; //alignmentCase->steps[stepInCase].step;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    unsigned long long errSaRange[4];
    unsigned long long err_l;
    unsigned long long err_r;
    unsigned char c;
    //printf("Alignment from %d to %d allowing %u-%u mismatches. i = %d\n", start,end, minMismatch, maxMismatch,i);
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    int branchCount = 0;
    char branchChar = 0;

    pos = start + step * i;
    len -= ( minMismatch - mismatchInserted - 1 ) * ( minMismatch > ( mismatchInserted + 1 ) );

    while ( i < len && l <= r )
    {
        if ( r - l < ceThreshold && pos >= ceStart && pos <= ceEnd )
        {
            saRanges[0] = l;
            saRanges[1] = r;
            saCount += BWTMismatchModelAnyDirection_CE3 ( qInput, i, mismatchInserted,
                       alignmentCase, stepInCase,
                       saRanges, occMismatch, nbMismatch, occQuality, algnResult );
            return saCount;
        }

        c = convertedKey[pos];
        BWTAllOccValue ( bwt, l, oL );
        BWTAllOccValue ( bwt, r + 1, oR );

        //if (mismatchInserted<maxMismatch) {
        branchCount = 0;
        branchChar = 0;

        for ( k = 0; k < ALPHABET_SIZE; k++ )
        {
            if ( k != c )
            {
                errSaRange[0] = bwt->cumulativeFreq[k] + oL[k] + 1;
                errSaRange[1] = bwt->cumulativeFreq[k] + oR[k];

                if ( errSaRange[0] <= errSaRange[1] )
                {
                    saCount += BWTModelSwitchBackward3 ( qInput, i + 1, mismatchInserted + 1,
                                                         alignmentCase, stepInCase,
                                                         errSaRange, occMismatch, nbMismatch, occQuality + keyQuality[pos], algnResult );

                    if ( rOutput->IsClosed == 1 ) { return saCount; }
                }
            }

            branchCount += ( oR[k] > oL[k] );
            branchChar |= k * ( oR[k] > oL[k] );
        }

        //}
        // Non-branching Mismatch Handler
        if ( oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < qSetting->MaxNBMismatch )
        {
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            nbMismatch++;
        }
        else
        {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
        }

        i += 1;
        pos += step;
    }

    if ( mismatchInserted >= minMismatch && l <= r )
    {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saCount += BWTModelSwitchBackward3 ( qInput, 0, 0,
                                             alignmentCase, stepInCase + 1,
                                             saRanges, occMismatch + mismatchInserted, nbMismatch, occQuality, algnResult );
    }

    return saCount;
}


// BWTEditModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTEditModelAnyDirection3 ( SRAQueryInput * qInput, int i, int editInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occEdit, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    int step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[ ( step < 0 ) * 2];
    unsigned long long r = saRanges[ ( step < 0 ) * 2 + 1];
    unsigned long long rev_l = saRanges[ ( step > 0 ) * 2];
    unsigned long long rev_r = saRanges[ ( step > 0 ) * 2 + 1];
    unsigned long long saCount = 0;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    int pos, misLen, delLen;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minEdit = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxEdit = alignmentCase->steps[stepInCase].MaxError;
    int residueEdit = maxEdit - editInserted;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    unsigned long long err_rev_l;
    unsigned long long err_rev_r;
    unsigned long long err_r;
    unsigned long long err_l;
    unsigned long long errSaRange[4];
    unsigned long long j;
    unsigned char c, nc;
    int k;
    //printf("[BWTEditModelAnyDirection] Read %llu Alignment from %d to %d, i = %d\n", qInfo->ReadId,start,end,i);
    //printf("[BWTEditModelAnyDirection] ErrorType=%d  allowing %u-%u edits.\n", errType, minEdit, maxEdit);
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    nc = ALPHABET_SIZE + 1;
    pos = start + step * i;
    misLen = len - ( minEdit - editInserted - 1 ) * ( minEdit > ( editInserted + 1 ) );
    delLen = len - ( maxEdit - editInserted - 1 ) * ( maxEdit > ( editInserted + 1 ) );

    while ( i < len && l <= r )
    {
        //printf("Current progress: %d(%d) within %d characters. stepping by %d.\n", i, pos, len, step);
        /*if ( r-l < ceThreshold && pos >=ceStart && pos<=ceEnd ) {
            k=0;
            saRanges[(step<0)*2] = l;
            saRanges[(step<0)*2+1] = r;
            saRanges[(step>0)*2] = rev_l;
            saRanges[(step>0)*2+1] = rev_r;
            for (j=saRanges[0];j<=saRanges[1];j++) {
                occEdits[k]=occEdit;
                occQualities[k]=occQuality;
                saPositions[k++]=(*ceBwt->_bwtSaValue)(ceBwt,j)-leftMostAligned[pos];

            }
            saCount+=BWTMismatchModelAnyDirection_CE(qInput,rOutput,i,editInserted,
                                    hsp,
                                    alignmentCase,stepInCase,k);
            return saCount;
        }*/
        c = convertedKey[pos];

        if ( pos + step >= 0 && pos + step <= qInfo->ReadLength - 1 ) {nc = convertedKey[pos + step];}
        else {nc = ALPHABET_SIZE + 1;}

        BWTAllOccValue ( bwt, rev_l, oL );
        BWTAllOccValue ( bwt, rev_r + 1, oR );
        oCount[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; k-- )
        {
            oCount[k] = oCount[k + 1] + oR[k + 1] - oL[k + 1];
        }

        if ( editInserted < maxEdit )
        {
            //Delete
            if ( c != nc )
            {
                //printf("[BWTEditModelAnyDirection] Delete at position %u (i = %d)\n", pos, i);
                if ( i < len - 1 )
                {
                    err_rev_l = bwt->cumulativeFreq[nc] + oL[nc] + 1;
                    err_rev_r = bwt->cumulativeFreq[nc] + oR[nc];
                    err_r = r - oCount[nc];
                    err_l = err_r - ( err_rev_r - err_rev_l );
                }
                else
                {
                    err_rev_l = rev_l;
                    err_rev_r = rev_r;
                    err_r = r;
                    err_l = l;
                }

                errSaRange[ ( step < 0 ) * 2] = err_l;
                errSaRange[ ( step < 0 ) * 2 + 1] = err_r;
                errSaRange[ ( step > 0 ) * 2] = err_rev_l;
                errSaRange[ ( step > 0 ) * 2 + 1] = err_rev_r;
                unsigned long long result = BWTModelSwitchAnyDirection3 ( qInput, i + 2, editInserted + 1,
                                            alignmentCase, stepInCase,
                                            errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                if ( result )
                { printf ( "[BWTEditModelAnyDirection:Case%d] Delete at position %u (i = %d)\n", alignmentCase->id, pos, i ); }

                saCount += result;

                if ( rOutput->IsClosed == 1 ) { return saCount; }
            }

            //Long delete
            for ( k = residueEdit; k > 1; k-- )
            {
                if ( len - i >= k && c != convertedKey[pos + k * step] )
                {
                    //printf("[BWTEditModelAnyDirection] Long Delete at position %u (i = %d, gap len = %d)\n", pos, i, k);
                    errSaRange[ ( step < 0 ) * 2] = l;
                    errSaRange[ ( step < 0 ) * 2 + 1] = r;
                    errSaRange[ ( step > 0 ) * 2] = rev_l;
                    errSaRange[ ( step > 0 ) * 2 + 1] = rev_r;
                    unsigned long long result = BWTModelSwitchAnyDirection3 ( qInput, i + k, editInserted + k,
                                                alignmentCase, stepInCase,
                                                errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                    if ( result )
                    { printf ( "[BWTEditModelAnyDirection:Case%d] Long Delete at position %u (i = %d, gap len = %d)\n", alignmentCase->id, pos, i, k ); }

                    saCount += result;

                    if ( rOutput->IsClosed == 1 ) { return saCount; }
                }
            }

            for ( k = 0; k < ALPHABET_SIZE; k++ )
            {
                err_rev_l = bwt->cumulativeFreq[k] + oL[k] + 1;
                err_rev_r = bwt->cumulativeFreq[k] + oR[k];
                err_r = r - oCount[k];
                err_l = err_r - ( err_rev_r - err_rev_l );

                if ( err_l <= err_r )
                {
                    if ( k != c && i < misLen )
                    {
                        //Mismatch
                        //printf("[BWTEditModelAnyDirection] Mismatch at position %u (i = %d)\n", pos, i);
                        errSaRange[ ( step < 0 ) * 2] = err_l;
                        errSaRange[ ( step < 0 ) * 2 + 1] = err_r;
                        errSaRange[ ( step > 0 ) * 2] = err_rev_l;
                        errSaRange[ ( step > 0 ) * 2 + 1] = err_rev_r;
                        unsigned long long result = BWTModelSwitchAnyDirection3 ( qInput, i + 1, editInserted + 1,
                                                    alignmentCase, stepInCase,
                                                    errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                        if ( result )
                        { printf ( "[BWTEditModelAnyDirection:Case%d] Mismatch at position %u (i = %d)\n", alignmentCase->id, pos, i ); }

                        saCount += result;

                        if ( rOutput->IsClosed == 1 ) { return saCount; }
                    }

                    //Insert
                    if ( ( k != c && ( i != 0 || errType == SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY ) ) )
                    {
                        //printf("[BWTEditModelAnyDirection] Insert at position %u (i = %d)\n", pos, i);
                        errSaRange[ ( step < 0 ) * 2] = err_l;
                        errSaRange[ ( step < 0 ) * 2 + 1] = err_r;
                        errSaRange[ ( step > 0 ) * 2] = err_rev_l;
                        errSaRange[ ( step > 0 ) * 2 + 1] = err_rev_r;
                        unsigned long long result = BWTModelSwitchAnyDirection3 ( qInput, i, editInserted + 1,
                                                    alignmentCase, stepInCase,
                                                    errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                        if ( result )
                        { printf ( "[BWTEditModelAnyDirection:Case%d] Insert at position %u (i = %d)\n", alignmentCase->id, pos, i ); }

                        saCount += result;

                        if ( rOutput->IsClosed == 1 ) { return saCount; }
                    }
                }
            }
        }

        rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
        rev_r = bwt->cumulativeFreq[c] + oR[c];
        r = r - oCount[c];
        l = r - ( rev_r - rev_l );
        i += 1;
        pos += step;
    }

    if ( editInserted >= minEdit && l <= r )
    {
        //Next Step
        saRanges[ ( step < 0 ) * 2] = l;
        saRanges[ ( step < 0 ) * 2 + 1] = r;
        saRanges[ ( step > 0 ) * 2] = rev_l;
        saRanges[ ( step > 0 ) * 2 + 1] = rev_r;
        saCount += BWTModelSwitchAnyDirection3 ( qInput, 0, 0,
                   alignmentCase, stepInCase + 1,
                   saRanges, occEdit + editInserted, nbMismatch, occQuality, algnResult );
    }

    return saCount;
}

// BWTEditModelBackward matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTEditModelBackward3 ( SRAQueryInput * qInput, int i, int editInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occEdit, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    unsigned long long saCount = 0;
    int pos, delLen, misLen;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = qInfo->ReadQuality;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minEdit = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxEdit = alignmentCase->steps[stepInCase].MaxError;
    int residueEdit = maxEdit - editInserted;
    int step = -1; //alignmentCase->steps[stepInCase].step;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    unsigned long long errSaRange[4];
    unsigned long long err_l;
    unsigned long long err_r;
    unsigned char nc, c, k;
    //printf("[BWTEditModelBackward] Read %llu Alignment from %d to %d, i = %d\n", qInfo->ReadId,start,end,i);
    //printf("[BWTEditModelBackward] ErrorType=%d  allowing %u-%u edits.\n", errType, minEdit, maxEdit);
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    pos = start + step * i;
    misLen = len - ( minEdit - editInserted - 1 ) * ( minEdit > ( editInserted + 1 ) );
    delLen = len - ( maxEdit - editInserted - 1 ) * ( maxEdit > ( editInserted + 1 ) );
    nc = ALPHABET_SIZE + 1;

    while ( i < len && l <= r )
    {
        /*if ( r-l < ceThreshold && pos >=ceStart && pos<=ceEnd) {
            k=0;
            for (err_l=l;err_l<=r;err_l++) {
                occEdits[k]=occEdit;
                occQualities[k]=occQuality;
                saPositions[k++]=(*ceBwt->_bwtSaValue)(ceBwt,err_l)-leftMostAligned[pos];
            }
            saCount += BWTMismatchModelAnyDirection_CE(qInput,rOutput,i,editInserted,
                                    hsp,
                                    alignmentCase,stepInCase,k);
            return saCount;
        }*/
        //printf("Current progress: %d(%d) within %d characters. stepping by %d.\n", i, pos, len, step);
        c = convertedKey[pos];

        if ( pos + step >= 0 && pos + step <= qInfo->ReadLength - 1 ) {nc = convertedKey[pos + step];}
        else {nc = ALPHABET_SIZE + 1;}

        BWTAllOccValue ( bwt, l, oL );
        BWTAllOccValue ( bwt, r + 1, oR );

        if ( editInserted < maxEdit )
        {
            //Delete
            if ( c != nc )
            {
                if ( i < len - 1 )
                {
                    err_l = bwt->cumulativeFreq[nc] + oL[nc] + 1;
                    err_r = bwt->cumulativeFreq[nc] + oR[nc];
                }
                else
                {
                    err_r = r;
                    err_l = l;
                }

                errSaRange[0] = err_l;
                errSaRange[1] = err_r;
                unsigned long long result = BWTModelSwitchBackward3 ( qInput, i + 2, editInserted + 1,
                                            alignmentCase, stepInCase,
                                            errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                if ( result )
                { printf ( "[BWTEditModelBackward:Case%d] Delete at position %u (i = %d)\n", alignmentCase->id, pos, i ); }

                saCount += result;

                if ( rOutput->IsClosed == 1 ) { return saCount; }
            }

            //Long delete
            for ( k = residueEdit; k > 1; k-- )
            {
                if ( len - i >= k && c != convertedKey[pos + k * step] )
                {
                    errSaRange[0] = l;
                    errSaRange[1] = r;
                    unsigned long long result = BWTModelSwitchBackward3 ( qInput, i + k, editInserted + k,
                                                alignmentCase, stepInCase,
                                                errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                    if ( result )
                    { printf ( "[BWTEditModelBackward:Case%d] Long Delete at position %u (i = %d, gap len = %d)\n", alignmentCase->id, pos, i, k ); }

                    saCount += result;

                    if ( rOutput->IsClosed == 1 ) { return saCount; }
                }
            }

            for ( k = 0; k < ALPHABET_SIZE; k++ )
            {
                err_l = bwt->cumulativeFreq[k] + oL[k] + 1;
                err_r = bwt->cumulativeFreq[k] + oR[k];

                if ( err_l <= err_r )
                {
                    if ( k != c && i < misLen )
                    {
                        //Mismatch
                        //printf("[BWTEditModelBackward] Mismatch at position %u (i = %d)\n", pos, i);
                        errSaRange[0] = err_l;
                        errSaRange[1] = err_r;
                        unsigned long long result = BWTModelSwitchBackward3 ( qInput, i + 1, editInserted + 1,
                                                    alignmentCase, stepInCase,
                                                    errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                        if ( result )
                        { printf ( "[BWTEditModelBackward:Case%d] Mismatch at position %u (i = %d)\n", alignmentCase->id, pos, i ); }

                        saCount += result;

                        if ( rOutput->IsClosed == 1 ) { return saCount; }
                    }

                    //Insert
                    if ( ( k != c && ( errType == SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY || i != 0 ) ) )
                    {
                        //printf("[BWTEditModelBackward] Insert at position %u (i = %d)\n", pos, i);
                        errSaRange[0] = err_l;
                        errSaRange[1] = err_r;
                        unsigned long long result = BWTModelSwitchBackward3 ( qInput, i, editInserted + 1,
                                                    alignmentCase, stepInCase,
                                                    errSaRange, occEdit, nbMismatch, occQuality, algnResult );

                        if ( result )
                        { printf ( "[BWTEditModelBackward:Case%d] Insert at position %u (i = %d)\n", alignmentCase->id, pos, i ); }

                        saCount += result;

                        if ( rOutput->IsClosed == 1 ) { return saCount; }
                    }
                }
            }
        }

        l = bwt->cumulativeFreq[c] + oL[c] + 1;
        r = bwt->cumulativeFreq[c] + oR[c];
        i += 1;
        pos += step;
    }

    if ( editInserted >= minEdit && l <= r )
    {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saCount += BWTModelSwitchBackward3 ( qInput, 0, 0,
                                             alignmentCase, stepInCase + 1,
                                             saRanges, occEdit + editInserted, nbMismatch, occQuality, algnResult );
    }

    return saCount;
}


unsigned long long BWTModelSwitchAnyDirection3 ( SRAQueryInput * qInput, int i, int errorInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occError, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    unsigned long long j;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    int OutputType = qSetting->OutputType;
    HSP * hsp = aIndex->hsp;
    BWT * bwt = aIndex->bwt;
    HOCC * highOcc = aIndex->highOcc;
    OCC * occ = qSetting->occ;
    FILE * outFilePtr = qSetting->OutFilePtr;

    if ( alignmentCase->steps[stepInCase].type == SRA_STEP_TYPE_COMPLETE )
    {
        //for (j=l;j<=r;j++) {
        //  printf("ReadId=%u Occ=%llu Error=%d\n",qInfo->ReadId,(*ceBwt->_bwtSaValue)(ceBwt,j),errorInserted+occError);
        //}
        //rOutput->WithError[occError+errorInserted]+=r-l+1;
        //return OCCReportSARange(qInput,l,r,occError+errorInserted,occQuality);
        if ( r >= l )
        {
            addSAToAlgnResult ( algnResult, l, r, qInfo->ReadStrand );
            return r - l + 1;
        }
        else
        {
            return 0;
        }

        // return OCCReportSARange(occ,bwt,hsp,highOcc,outFilePtr,qInfo->ReadId,qInfo->ReadStrand, l,r);
    }

    int errorType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;

    if ( qInput->QueryOutput->IsClosed == 1 )
    {
        return 0;
    }

    if ( errorType == SRA_STEP_ERROR_TYPE_NO_ERROR || errorInserted >= maxError )
    {
        return BWTExactModelAnyDirection3 ( qInput, i, errorInserted, alignmentCase, stepInCase, saRanges, occError,  nbMismatch, occQuality, algnResult );
    }
    else if ( errorType == SRA_STEP_ERROR_TYPE_MISMATCH_ONLY )
    {
        return BWTMismatchModelAnyDirection3 ( qInput, i, errorInserted, alignmentCase, stepInCase, saRanges, occError,  nbMismatch, occQuality, algnResult );
    }
    else if ( errorType == SRA_STEP_ERROR_TYPE_EDIT_DISTANCE || errorType == SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY )
    {
        return BWTEditModelAnyDirection3 ( qInput, i, errorInserted, alignmentCase, stepInCase, saRanges, occError,  nbMismatch, occQuality, algnResult );
    }

    return 0;
}

unsigned long long BWTModelSwitchBackward3 ( SRAQueryInput * qInput,  int i, int errorInserted,
        SRACase * alignmentCase, int stepInCase,
        unsigned long long * saRanges,
        int occError, uint8_t nbMismatch, char occQuality, SingleAlgnResult * algnResult )
{
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    unsigned long long j;
    SRAQueryInfo * qInfo = qInput->QueryInfo;
    SRASetting * qSetting = qInput->QuerySetting;
    SRAIndex * aIndex = qInput->AlgnmtIndex;
    SRAQueryResultCount * rOutput  = qInput->QueryOutput;
    int OutputType = qSetting->OutputType;
    HSP * hsp = aIndex->hsp;
    BWT * bwt = aIndex->bwt;
    HOCC * highOcc = aIndex->highOcc;
    OCC * occ = qSetting->occ;
    FILE * outFilePtr = qSetting->OutFilePtr;

    //printf("BWTModelSwitchBackward: Step %d in case %d. i = %d\n", stepInCase, alignmentCase->id, i);

    if ( alignmentCase->steps[stepInCase].type == SRA_STEP_TYPE_COMPLETE )
    {
        //for (j=l;j<=r;j++) {
        //  printf("ReadId=%u Occ=%llu Error=%d\n",qInfo->ReadId,(*ceBwt->_bwtSaValue)(ceBwt,j),errorInserted+occError);
        //}
        //rOutput->WithError[occError+errorInserted]+=r-l+1;
        //return OCCReportSARange(qInput,l,r,occError+errorInserted,occQuality);
        if ( r >= l )
        {
            addSAToAlgnResult ( algnResult, l, r, qInfo->ReadStrand );
            return r - l + 1;
        }
        else
        {
            return 0;
        }

        // return OCCReportSARange(occ,bwt,hsp,highOcc,outFilePtr,qInfo->ReadId,qInfo->ReadStrand, l,r);
    }

    int errorType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;

    if ( qInput->QueryOutput->IsClosed == 1 )
    {
        return 0;
    }

    if ( errorType == SRA_STEP_ERROR_TYPE_NO_ERROR || errorInserted >= maxError )
    {
        return BWTExactModelBackward3 ( qInput,  i, errorInserted, alignmentCase, stepInCase, saRanges, occError, nbMismatch, occQuality, algnResult );
    }
    else  if ( errorType == SRA_STEP_ERROR_TYPE_MISMATCH_ONLY )
    {
        return BWTMismatchModelBackward3 ( qInput, i, errorInserted, alignmentCase, stepInCase, saRanges, occError, nbMismatch, occQuality, algnResult );
    }
    else if ( errorType == SRA_STEP_ERROR_TYPE_EDIT_DISTANCE || errorType == SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY )
    {
        return BWTEditModelBackward3 ( qInput, i, errorInserted, alignmentCase, stepInCase, saRanges, occError, nbMismatch, occQuality, algnResult );
    }

    return 0;
}
