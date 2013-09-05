/*
 *
 *    DV-DPForSingleReads.cu
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

#include "DV-DPForSingleReads.h"
#include "OutputDPResult.h"
#include "DV-DPfunctions.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <semaphore.h>
#include <cuda.h>

#include <functional>
#include <vector>
using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;

using namespace SingleDP_Space;
class SingleDPWrapper
{
        UnalignedSinglesArrays * unalignedReads;
        uint * queries, *upkdReadLengths, *origReadIDs;
        char * upkdQueryNames, *upkdQualities;
        uint maxReadLength;
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

        SingleDPWrapper ( UnalignedSinglesArrays * unalignedReads,
                          uint * queries, uint * upkdReadLengths, uint * origReadIDs,
                          char * upkdQueryNames, char * upkdQualities, uint maxReadLength,
                          Soap3Index * index,
                          uint * _bwt, uint * _revBwt,
                          uint * _occ, uint * _revOcc,
                          int alignmentType, DPParameters * dpParameters,
                          uint accumReadNum, int outputFormat,
                          FILE * outputFile, samfile_t * samOutputDPFilePtr )
        {
            MC_MemberCopy2 ( this->, , unalignedReads, dpParameters );
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
        ~SingleDPWrapper ()
        {
            delete soap3Wrapper;
        }

        void seeding ( QueryIDStream * inputStream,
                       CandidateStream * canStream, QueryIDStream * unseededIDStream )
        {
            SingleEndSeedingEngine::
            performSeeding (
                /* input */
                inputStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
                soap3Wrapper, index,
                /* output */
                canStream, unseededIDStream );
        }
        void alignment ( CandidateStream * canStream,
                         QueryIDStream * unalignedIDStream )
        {
            uint alignedRead = 0, alignment = 0;
            SingleEndAlignmentEngine::
            performAlignment (
                /* input */
                canStream, dpParameters,
                queries, upkdReadLengths, maxReadLength,
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
        void singleDPOneRound ( QueryIDStream * inputStream,
                                QueryIDStream * unseededIDStream, QueryIDStream * unalignedIDStream )
        {
            CandidateStream * canStream = new CandidateStream;
            soap3Wrapper->copyIndex ();
            seeding ( inputStream, canStream, unseededIDStream );
            soap3Wrapper->freeIndex ();
            alignment ( canStream, unalignedIDStream );
            delete canStream;
        }
        void outputUnaligned ( QueryIDStream * unalignedIDStream )
        {
            DPSOutputUnalignedReads (
                unalignedIDStream,
                queries, upkdReadLengths, maxReadLength,
                index, upkdQueryNames, origReadIDs, upkdQualities,
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
            singleDPOneRound ( input, unaligned_round1, unaligned_round1 );
            delete input;
            outputUnaligned ( unaligned_round1 );
            delete unaligned_round1;
        }
};


void DPForUnalignSingle2 ( UnalignedSinglesArrays * unalignedReads,
                           unsigned int * queries, unsigned int * upkdReadLengths, unsigned int * origReadIDs, char * upkdQueryNames, char * upkdQualities,
                           unsigned int maxReadLength,
                           Soap3Index * index,
                           unsigned int * _bwt, unsigned int * _revBwt,
                           unsigned int * _occ, unsigned int * _revOcc,
                           int alignmentType, DPParameters * dpParameters,
                           unsigned int & numDPAlignedRead, unsigned int & numDPAlignment,
                           unsigned int accumReadNum, int outputFormat,
                           FILE * outputFile, samfile_t * samOutputDPFilePtr )
{
    using namespace SingleDP_Space;
    SingleDPWrapper
    singleDPWrapper (
        unalignedReads,
        queries, upkdReadLengths, origReadIDs,
        upkdQueryNames, upkdQualities, maxReadLength,
        index, _bwt,  _revBwt,
        _occ,  _revOcc,
        alignmentType, dpParameters,
        accumReadNum, outputFormat,
        outputFile, samOutputDPFilePtr );
    singleDPWrapper.run ();
    numDPAlignedRead = singleDPWrapper.numDPAlignedRead;
    numDPAlignment = singleDPWrapper.numDPAlignment;
}
