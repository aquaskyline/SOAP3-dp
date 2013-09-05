/*
 *
 *    DV-SemiDP.cu
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

#include "DV-SemiDP.h"
#include "DV-DPfunctions.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

using namespace DP_Space;
class SemiDPWrapper
{
        Soap3Index * index;
        ReadInputForDPArrays * readInputArrays;
        BothUnalignedPairs * unalignedReads;
        int insert_high, insert_low;
        uint * queries, *upkdReadLengths, *origReadIDs;
        char * upkdQueryNames, *upkdQualities;
        uint maxReadLength;
        int peStrandLeftLeg, peStrandRightLeg;
        uint * _bwt, *_revBwt, *_occ, *_revOcc;
        int alignmentType;
        uint accumReadNum;
        int outputFormat;
        FILE * outputFile;
        samfile_t * samOutputDPFilePtr;

        DPParameters * dpParameters;
        SOAP3Wrapper<void> * soap3Wrapper;

    public:
        uint numDPAlignedRead, numDPAlignment;

        SemiDPWrapper ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                        uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                        char * upkdQueryNames, char * upkdQualities, uint maxReadLength,
                        int peStrandLeftLeg, int peStrandRightLeg,
                        Soap3Index * index,
                        int alignmentType, DPParameters * dpParameters,
                        uint accumReadNum, int outputFormat,
                        FILE * outputFile, samfile_t * samOutputDPFilePtr,
                        BothUnalignedPairs * unalignedReads )
        {
            MC_MemberCopy3 ( this->, , readInputArrays, dpParameters, unalignedReads );
            MC_MemberCopy4 ( this->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
            MC_MemberCopy3 ( this->, , queries, upkdReadLengths, origReadIDs );
            MC_MemberCopy3 ( this->, , upkdQueryNames, upkdQualities, maxReadLength );
            MC_MemberCopy ( this->, , index );
            MC_MemberCopy2 ( this->, , alignmentType, accumReadNum );
            MC_MemberCopy3 ( this->, , outputFormat, outputFile, samOutputDPFilePtr );
            soap3Wrapper = NULL;
        }

        SemiDPWrapper ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                        uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                        char * upkdQueryNames, char * upkdQualities, uint maxReadLength,
                        int peStrandLeftLeg, int peStrandRightLeg,
                        Soap3Index * index,
                        uint * _bwt, uint * _revBwt,
                        uint * _occ, uint * _revOcc,
                        int alignmentType, DPParameters * dpParameters,
                        uint accumReadNum, int outputFormat,
                        FILE * outputFile, samfile_t * samOutputDPFilePtr,
                        BothUnalignedPairs * unalignedReads )
        {
            MC_MemberCopy3 ( this->, , readInputArrays, dpParameters, unalignedReads );
            MC_MemberCopy4 ( this->, , insert_high, insert_low, peStrandLeftLeg, peStrandRightLeg );
            MC_MemberCopy3 ( this->, , queries, upkdReadLengths, origReadIDs );
            MC_MemberCopy3 ( this->, , upkdQueryNames, upkdQualities, maxReadLength );
            MC_MemberCopy ( this->, , index );
            MC_MemberCopy4 ( this->, , _bwt, _revBwt, _occ, _revOcc );
            MC_MemberCopy2 ( this->, , alignmentType, accumReadNum );
            MC_MemberCopy3 ( this->, , outputFormat, outputFile, samOutputDPFilePtr );
            soap3Wrapper =
                new SOAP3Wrapper<void> ( index,
                                         _bwt, _revBwt,
                                         _occ, _revOcc );
        }
        ~SemiDPWrapper ()
        {
            if ( soap3Wrapper )
            { delete soap3Wrapper; }
        }

        void alignment ( HalfEndOccStream * canStream,
                         QueryIDStream * unalignedIDStream )
        {
            uint alignedRead = 0, alignment = 0;
            HalfEndAlignmentEngine::
            performAlignment (
                /* input */
                canStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                upkdQueryNames, origReadIDs, upkdQualities,
                index,
                alignmentType,
                accumReadNum, outputFormat,
                outputFile, samOutputDPFilePtr,
                /* output */
                unalignedIDStream, alignedRead, alignment );
            numDPAlignedRead += alignedRead;
            numDPAlignment += alignment;
        }

        void passUnalignedToDeepDP ( QueryIDStream * unalignedIDStream )
        {
            for ( int i = 0; i < unalignedIDStream->data->size (); i++ )
            {
                addReadIDToBothUnalignedPairs ( unalignedReads, ( * ( unalignedIDStream->data ) ) [i] );
            }
        }

        // for new default dp & deep dp
        void seeding ( HalfEndOccStream * inputStream,
                       DeepDP_Space::CandidateStream * canStream, QueryIDStream * unseededIDStream,
                       int stage )
        {
            DeepDP_Space::PairEndSeedingEngine::
            performSeeding (
                /* input */
                inputStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                stage,
                soap3Wrapper, index,
                /* output */
                canStream, unseededIDStream );
        }

        void seeding ( QueryIDStream * inputStream,
                       DeepDP_Space::CandidateStream * canStream, QueryIDStream * unseededIDStream,
                       int stage )
        {
            DeepDP_Space::PairEndSeedingEngine::
            performSeeding (
                /* input */
                inputStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                stage,
                soap3Wrapper, index,
                /* output */
                canStream, unseededIDStream );
        }

        void seeding_ext ( QueryIDStream * inputStream,
                           DeepDP_Space::CandidateStream * canStream, QueryIDStream * unseededIDStream )
        {
            QueryIDStream * tmpStream = new QueryIDStream;
            dpParameters->paramRead[0].maxHitNum = 100;
            dpParameters->paramRead[1].maxHitNum = 100;
            //      seeding(inputStream, canStream, tmpStream, STAGE_NEW_DEFAULT_DP_ROUND1);
            // printf("unseeded 1: %d\n", (int) tmpStream->data->size());
            dpParameters->paramRead[0].maxHitNum = 500;
            dpParameters->paramRead[1].maxHitNum = 500;
            // seeding(tmpStream, canStream, unseededIDStream, STAGE_NEW_DEFAULT_DP_ROUND2);
            // printf("unseeded 2: %d\n", (int) unseededIDStream->data->size());
            delete tmpStream;
        }

        void alignment ( DeepDP_Space::CandidateStream * canStream,
                         QueryIDStream * unalignedIDStream )
        {
            uint alignedRead = 0, alignment = 0;
            DeepDP_Space::PairEndAlignmentEngine::
            performAlignment (
                /* input */
                canStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                upkdQueryNames, origReadIDs, upkdQualities,
                index,
                alignmentType,
                accumReadNum, outputFormat,
                outputFile, samOutputDPFilePtr,
                /* output */
                unalignedIDStream, alignedRead, alignment );
            numDPAlignedRead += alignedRead;
            numDPAlignment += alignment;
        }

        void outputUnaligned ( QueryIDStream * unalignedIDStream )
        {
            DeepDP_Space::DP2OutputUnalignedReads (
                unalignedIDStream,
                queries, upkdReadLengths, maxReadLength,
                index, peStrandLeftLeg, peStrandRightLeg,
                upkdQueryNames, origReadIDs, upkdQualities,
                accumReadNum, outputFormat,
                outputFile, samOutputDPFilePtr
            );
        }

        void newDefaultDP ( HalfEndOccStream * inputStream,
                            QueryIDStream * unseededIDStream, QueryIDStream * unalignedIDStream )
        {
            DeepDP_Space::CandidateStream * canStream = new DeepDP_Space::CandidateStream;
            soap3Wrapper->copyIndex ();
            seeding ( inputStream, canStream, unseededIDStream, STAGE_NEW_DEFAULT_DP );
            soap3Wrapper->freeIndex ();
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }

        void newDefaultDP2 ( HalfEndOccStream * inputStream,
                             QueryIDStream * unseededIDStream, QueryIDStream * unalignedIDStream )
        {
            DeepDP_Space::CandidateStream * canStream = new DeepDP_Space::CandidateStream;
            QueryIDStream * tmpStream = new QueryIDStream;
            soap3Wrapper->copyIndex ();
            seeding ( inputStream, canStream, tmpStream, STAGE_NEW_DEFAULT_DP );
            //      printf("unseeded 1: %d\n", (int) tmpStream->data->size());
            seeding_ext ( tmpStream, canStream, unseededIDStream );
            delete tmpStream;
            soap3Wrapper->freeIndex ();
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }

        void run ()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            HalfEndOccStream * canStream = new HalfEndOccStream ( readInputArrays, index->sraIndex->bwt );
            QueryIDStream * unalignedIDStream = new QueryIDStream;
            alignment ( canStream, unalignedIDStream );
            // printf("total unaligned: %d\n", (int) unalignedIDStream->data->size());
            passUnalignedToDeepDP ( unalignedIDStream );
            delete canStream;
            delete unalignedIDStream;
        }

        void run2 ()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            HalfEndOccStream * occStream = new HalfEndOccStream ( readInputArrays, index->sraIndex->bwt );
            QueryIDStream * unalignedIDStream = new QueryIDStream;
            newDefaultDP ( occStream, unalignedIDStream, unalignedIDStream );
            // printf("total unaligned: %d\n", (int) unalignedIDStream->data->size());
            passUnalignedToDeepDP ( unalignedIDStream );
            delete occStream;
            delete unalignedIDStream;
        }

        void run3 ()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            HalfEndOccStream * occStream = new HalfEndOccStream ( readInputArrays, index->sraIndex->bwt );
            QueryIDStream * unalignedIDStream = new QueryIDStream;
            newDefaultDP2 ( occStream, unalignedIDStream, unalignedIDStream );
            // printf("total unaligned: %d\n", (int) unalignedIDStream->data->size());
            // passUnalignedToDeepDP(unalignedIDStream);
            outputUnaligned ( unalignedIDStream );
            delete occStream;
            delete unalignedIDStream;
        }

};

// TODO
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MAIN start here ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

// To perform semi-global DP
void semiGlobalDP2 ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                     unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                     unsigned int maxReadLength, int peStrandLeftLeg, int peStrandRightLeg,
                     Soap3Index * index,
                     unsigned int * _bwt, unsigned int * _revBwt,
                     unsigned int * _occ, unsigned int * _revOcc,
                     int alignmentType, DPParameters * dpParameters,
                     unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                     unsigned int accumReadNum, int outputFormat,
                     FILE * outputFile, samfile_t * samOutputDPFilePtr,
                     BothUnalignedPairs * unalignedReads )
{
    // readInput: input reads for DP3
    // insert_high : maximum value of insert size
    // insert_low: minimum value of insert size
    // upkdQueries: the sequence of all reads. The first character of the read with ID "i"
    //              is on the position "i*maxReadLength" of the array
    // readLengths: the read length of all reads
    // peStrandLeftLeg: the required strand of the left alignment (i.e. 1: +ve; 2: -ve)
    // peStrandRightLeg: the required strand of the right alignment (i.e. 1: +ve; 2: -ve)
    // BWT: bwt structure with SA table inside
    // HSP: hsp structure with packed sequence inside
    // alignmentType: 1: All valid alignment; 2: all best alignment
    // numDPResult: the total number of alignment results outputted by DP
    SemiDPWrapper
    semiDPWrapper (
        readInputArrays, insert_high, insert_low,
        queries, readLengths, origReadIDs,
        upkdQueryNames, upkdQualities, maxReadLength,
        peStrandLeftLeg, peStrandRightLeg,
        index,
        _bwt,  _revBwt,
        _occ,  _revOcc,
        alignmentType, dpParameters,
        accumReadNum, outputFormat,
        outputFile, samOutputDPFilePtr,
        unalignedReads );
    semiDPWrapper.run ();
    numDPAlignedRead = semiDPWrapper.numDPAlignedRead;
    numDPAlignment = semiDPWrapper.numDPAlignment;
    return;
}


// To perform new semi-global DP
void newSemiGlobalDP ( ReadInputForDPArrays * readInputArrays, int insert_high, int insert_low,
                       unsigned int * queries, unsigned int * readLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                       unsigned int maxReadLength, int peStrandLeftLeg, int peStrandRightLeg,
                       Soap3Index * index,
                       unsigned int * _bwt, unsigned int * _revBwt,
                       unsigned int * _occ, unsigned int * _revOcc,
                       int alignmentType, DPParameters * dpParameters,
                       unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                       unsigned int accumReadNum, int outputFormat,
                       FILE * outputFile, samfile_t * samOutputDPFilePtr,
                       BothUnalignedPairs * unalignedReads )
{
    // readInput: input reads for DP
    // insert_high : maximum value of insert size
    // insert_low: minimum value of insert size
    // upkdQueries: the sequence of all reads. The first character of the read with ID "i"
    //              is on the position "i*maxReadLength" of the array
    // readLengths: the read length of all reads
    // peStrandLeftLeg: the required strand of the left alignment (i.e. 1: +ve; 2: -ve)
    // peStrandRightLeg: the required strand of the right alignment (i.e. 1: +ve; 2: -ve)
    // BWT: bwt structure with SA table inside
    // HSP: hsp structure with packed sequence inside
    // alignmentType: 1: All valid alignment; 2: all best alignment
    // numDPResult: the total number of alignment results outputted by DP
    SemiDPWrapper
    semiDPWrapper (
        readInputArrays, insert_high, insert_low,
        queries, readLengths, origReadIDs,
        upkdQueryNames, upkdQualities, maxReadLength,
        peStrandLeftLeg, peStrandRightLeg,
        index,
        _bwt,  _revBwt,
        _occ,  _revOcc,
        alignmentType, dpParameters,
        accumReadNum, outputFormat,
        outputFile, samOutputDPFilePtr,
        unalignedReads );
    semiDPWrapper.run2 ();
    numDPAlignedRead = semiDPWrapper.numDPAlignedRead;
    numDPAlignment = semiDPWrapper.numDPAlignment;
    return;
}

