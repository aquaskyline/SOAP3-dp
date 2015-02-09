/*

   BWTLoad.h        BWTLoad - Loading BWT index into memory

   Copyright (C) 2006, Wong Chi Kwong.

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

*/

#ifndef __BWT_LOAD_H__
#define __BWT_LOAD_H__

#include "BWT.h"
#include "HSP.h"
#include "Socket.h"

// This is the only message send by the resident process to a client
// If bwt == NULL signals an error condition

typedef struct BWTIndexPointer {
    BWT *bwt;
    SaIndexRange *saIndexRange;
    int saIndexRangeNumOfChar;
    unsigned int *packedDNA;
    Annotation *annotation;
    SeqOffset *seqOffset;
    int numOfSeq;
    Ambiguity *ambiguity;
    int numOfAmbiguity;
    unsigned int minSeqLength;
} BWTIndexPointer;

// This is the only message send by a client to the resident process
// If BWTCodeFileName[0] == '\0', the resident process will stop

typedef struct BWTIndexFileName {
    char BWTCodeFileName[MAX_FILENAME_LEN+1];
    char BWTOccValueFileName[MAX_FILENAME_LEN+1];
    char SaValueFileName[MAX_FILENAME_LEN+1];
    char SaIndexFileName[MAX_FILENAME_LEN+1];
    char PackedDNAFileName[MAX_FILENAME_LEN+1];
    char AnnotationFileName[MAX_FILENAME_LEN+1];
    char AmbiguityFileName[MAX_FILENAME_LEN+1];
} BWTIndexFileName;

#define BWT_INDEX_SOCKET_TYPE    LOCAL_SOCKET
#define BWT_INDEX_SOCKET_NAME    "BWT_INDEX_SOCKET"

#endif
