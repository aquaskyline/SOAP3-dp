//
//    SRACore.h
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
/////////////////////////////////////////////////////
/*

   Modification History 
   
   Date   : 30th October 2011
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __SRA_CORE_H__
#define __SRA_CORE_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../Release.h"
#include "../2bwt-lib/BWT.h"
#include "../2bwt-lib/HSP.h"
#include "OCC.h"
#include "HOCC.h"
#include "LT.h"
#include "../samtools-0.1.18/sam.h"
#include "../HSPAux.h"

static const char soap3DnaComplement[ALPHABET_SIZE]        = { 3, 2, 1, 0 };

#define SRA_MIN_READ_LENGTH                         25
#define SRA_MAX_READ_LENGTH                         1024
#define MAX_NUM_OF_MISMATCH                         5
#define MAX_NUM_OF_INDEL                            5
#define MAX_NUM_OF_NBM_ERROR                        10

// MAX_NUM_OF_ERROR must be greater than both
// MAX(MAX_NUM_OF_MISMATCH and MAX_NUM_OF_INDEL) + MAX_NUM_OF_NBM_ERROR
#define MAX_NUM_OF_ERROR                            15

#define SRA_CHAR_ERROR_TYPE_NONE                    0
#define SRA_CHAR_ERROR_TYPE_MISMATCH                1
#define SRA_CHAR_ERROR_TYPE_INSERT                  2
#define SRA_CHAR_ERROR_TYPE_DELETE                  3

#define SRA_REPORT_UNIQUE_BEST                           0
#define SRA_REPORT_RANDOM_BEST                           1
#define SRA_REPORT_ALL                                   2
#define SRA_REPORT_ALL_BEST                              3
#define SRA_REPORT_BEST_QUALITY                          4
#define SRA_REPORT_DETAIL                                5
#define SRA_REPORT_EXACT_NUM_ERROR                       6

#define SRA_OUTPUT_FORMAT_DEFAULT               0
#define SRA_OUTPUT_FORMAT_PLAIN                 1
#define SRA_OUTPUT_FORMAT_SAM                   4
#define SRA_OUTPUT_FORMAT_BAM                   3
#define SRA_OUTPUT_FORMAT_SAM_API               2

typedef struct SRAError {
    uint16_t type: 4, position: 12;
} SRAError;

typedef struct SRAOccurrence
{
    unsigned int readID; // read ID
    uint8_t type;
    unsigned long long ambPosition;
    uint8_t strand;
    uint8_t mismatchCount;
    SRAError errors[MAX_NUM_OF_ERROR];
    uint16_t matchLen;
} SRAOccurrence;

typedef struct SRAQueryInfo
{
    unsigned long long ReadId;
    char * ReadName;
    unsigned long long ReadLength;
    unsigned char * ReportingReadCode; //This is identical to ReadCode in QueryInputPos; But different in QueryInputNeg
    unsigned char * ReadCode;
    unsigned char ReadStrand;
    char * ReadQuality;
} SRAQueryInfo;

typedef struct SRASetting
{
    int ReadStrand;
    int OutputType;
    int MaxError;
    uint8_t MaxNBMismatch;
    int ErrorType;
    FILE * OutFilePtr;
    char * OutFileName;
    char   OutFileFormat;

    //Parameter for file splitting
    unsigned int OutFilePart;
    unsigned int writtenCacheCount;
    unsigned int MaxOutputPerRead;

    OCC * occ;

    //Parameter for SAM
    samfile_t * SAMOutFilePtr;

} SRASetting;

typedef struct SRAIndex
{
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    HSPAux * hspaux;
    HOCC * highOcc;
    LT * lookupTable;
    LT * rev_lookupTable;
} SRAIndex;

#endif
