/*
 *
 *    DV-DPForBothUnalign.cu
 *    Soap3(gpu)
 *
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


#include "DV-DPForBothUnalign.h"
#include "DV-DPfunctions.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>
using namespace std;

#define DP2_SEEDING_BATCH_SIZE 128 * 1024
#define DP2_MARGIN(l) ((l>100) ? (l>>2) : 25)

#define MC_Max(x,y) (x > y ? x : y)

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

using namespace DeepDP_Space;
class DeepDPWrapper
{
        BothUnalignedPairsArrays * unalignedReads;
        int insert_high, insert_low;
        uint * queries, *upkdReadLengths, *origReadIDs;
        char * upkdQueryNames, *upkdQualities;
        uint maxReadLength;
        int peStrandLeftLeg, peStrandRightLeg;
        Soap3Index * index;
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

        DeepDPWrapper ( BothUnalignedPairsArrays * unalignedReads, int insert_high, int insert_low,
                        uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                        char * upkdQueryNames, char * upkdQualities, uint maxReadLength,
                        int peStrandLeftLeg, int peStrandRightLeg,
                        Soap3Index * index,
                        uint * _bwt, uint * _revBwt,
                        uint * _occ, uint * _revOcc,
                        int alignmentType, DPParameters * dpParameters,
                        uint accumReadNum, int outputFormat,
                        FILE * outputFile, samfile_t * samOutputDPFilePtr )
        {
            MC_MemberCopy3 ( this->, , unalignedReads, insert_high, insert_low );
            MC_MemberCopy3 ( this->, , queries, upkdReadLengths, origReadIDs );
            MC_MemberCopy3 ( this->, , upkdQueryNames, upkdQualities, maxReadLength );
            MC_MemberCopy3 ( this->, , dpParameters, peStrandLeftLeg, peStrandRightLeg );
            MC_MemberCopy ( this->, , index );
            MC_MemberCopy4 ( this->, , _bwt, _revBwt, _occ, _revOcc );
            MC_MemberCopy2 ( this->, , alignmentType, accumReadNum );
            MC_MemberCopy3 ( this->, , outputFormat, outputFile, samOutputDPFilePtr );
            soap3Wrapper =
                new SOAP3Wrapper<void> ( index,
                                         _bwt, _revBwt,
                                         _occ, _revOcc );
        }
        ~DeepDPWrapper ()
        {
            delete soap3Wrapper;
        }

        void seeding ( QueryIDStream * inputStream,
                       CandidateStream * canStream,
                       QueryIDStream * unseededIDStream,
                       int stage )
        {
            PairEndSeedingEngine::
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
                       CandidateStream * canStream,
                       QueryIDStream * tooManyHitIDStream, QueryIDStream * unseededIDStream,
                       int stage )
        {
            PairEndSeedingEngine::
            performSeeding (
                /* input */
                inputStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                insert_high, insert_low,
                peStrandLeftLeg, peStrandRightLeg,
                stage,
                soap3Wrapper, index,
                /* output */
                canStream,
                tooManyHitIDStream, unseededIDStream );
        }

        void seeding_ext ( QueryIDStream * inputStream,
                           CandidateStream * canStream, QueryIDStream * unseededIDStream )
        {
            QueryIDStream * tooManyHitIDStream = new QueryIDStream;
            seeding ( inputStream, canStream, tooManyHitIDStream, unseededIDStream, STAGE_DEEP_DP_ROUND1 );
            // printf("tooManyHit 1: %d\n", (int) tooManyHitIDStream->data->size());
            /* set deep dp round2 maxHitNum here */
            dpParameters->paramRead[0].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ_2;
            dpParameters->paramRead[1].maxHitNum = MAX_SEED_HITS_DEEP_DP_FOR_NORMAL_READ_2;
            seeding ( tooManyHitIDStream, canStream, unseededIDStream, STAGE_DEEP_DP_ROUND2 );
            // printf("unseeded 2: %d\n", (int) unseededIDStream->data->size());
            delete tooManyHitIDStream;
        }
        /*
            void seeding_ext(QueryIDStream *inputStream,
                             CandidateStream *canStream, QueryIDStream *unseededIDStream) {
                QueryIDStream *tmpStream = new QueryIDStream;
                seeding(inputStream, canStream, tmpStream, STAGE_DEEP_DP_ROUND1);
                printf("unseeded 1: %d\n", (int) tmpStream->data->size());

                seeding(tmpStream, canStream, unseededIDStream, STAGE_DEEP_DP_ROUND2);
                printf("unseeded 2: %d\n", (int) unseededIDStream->data->size());
                delete tmpStream;
            }
        */

        void alignment ( CandidateStream * canStream,
                         QueryIDStream * unalignedIDStream )
        {
            uint alignedRead = 0, alignment = 0;
            PairEndAlignmentEngine::
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

        void deepDP ( QueryIDStream * inputStream,
                      QueryIDStream * unseededIDStream, QueryIDStream * unalignedIDStream )
        {
            CandidateStream * canStream = new CandidateStream;
            soap3Wrapper->copyIndex ();
            seeding ( inputStream, canStream, unseededIDStream, STAGE_DEEP_DP_ROUND2 );
            soap3Wrapper->freeIndex ();
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }

        void deepDP_ext ( QueryIDStream * inputStream,
                          QueryIDStream * unseededIDStream, QueryIDStream * unalignedIDStream )
        {
            CandidateStream * canStream = new CandidateStream;
            soap3Wrapper->copyIndex ();
            seeding_ext ( inputStream, canStream, unseededIDStream );
            soap3Wrapper->freeIndex ();
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }

        void outputUnaligned ( QueryIDStream * unalignedIDStream )
        {
            DP2OutputUnalignedReads (
                unalignedIDStream,
                queries, upkdReadLengths, maxReadLength,
                index, peStrandLeftLeg, peStrandRightLeg,
                upkdQueryNames, origReadIDs, upkdQualities,
                accumReadNum, outputFormat,
                outputFile, samOutputDPFilePtr
            );
        }

        void run ()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            QueryIDStream * input = new QueryIDStream ( unalignedReads );
            QueryIDStream * unaligned_round1 = new QueryIDStream;
            deepDP ( input, unaligned_round1, unaligned_round1 );
            delete input;
            // printf("total unaligned: %d\n", (int) unaligned_round1->data->size());
            outputUnaligned ( unaligned_round1 );
            delete unaligned_round1;
        }

        void run2 ()
        {
            numDPAlignedRead = 0;
            numDPAlignment = 0;
            QueryIDStream * input = new QueryIDStream ( unalignedReads );
            QueryIDStream * unaligned_round1 = new QueryIDStream;
            deepDP_ext ( input, unaligned_round1, unaligned_round1 );
            delete input;
            // printf("total unaligned: %d\n", (int) unaligned_round1->data->size());
            outputUnaligned ( unaligned_round1 );
            delete unaligned_round1;
        }

};

//////////////////////////////////////////////////////////////////////////////////////
// The functions are for the DP designed for pairs of which both ends are unaligned //
//////////////////////////////////////////////////////////////////////////////////////

void DPForUnalignPairs2 ( BothUnalignedPairsArrays * unalignedReads, int insert_high, int insert_low,
                          unsigned int * queries, unsigned int * upkdReadLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                          unsigned int maxReadLength,
                          int peStrandLeftLeg, int peStrandRightLeg,
                          Soap3Index * index,
                          unsigned int * _bwt, unsigned int * _revBwt,
                          unsigned int * _occ, unsigned int * _revOcc,
                          int alignmentType, DPParameters * dpParameters,
                          unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                          unsigned int accumReadNum, int outputFormat,
                          FILE * outputFile, samfile_t * samOutputDPFilePtr )
{
    DeepDPWrapper
    deepDPWrapper (
        unalignedReads, insert_high, insert_low,
        queries, upkdReadLengths, origReadIDs,
        upkdQueryNames, upkdQualities, maxReadLength,
        peStrandLeftLeg, peStrandRightLeg,
        index,
        _bwt,  _revBwt,
        _occ,  _revOcc,
        alignmentType, dpParameters,
        accumReadNum, outputFormat,
        outputFile, samOutputDPFilePtr );
    deepDPWrapper.run2 ();
    numDPAlignedRead = deepDPWrapper.numDPAlignedRead;
    numDPAlignment = deepDPWrapper.numDPAlignment;
}


