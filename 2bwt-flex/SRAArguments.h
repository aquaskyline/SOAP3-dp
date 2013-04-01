/*
 *
 *    SRAArguments.h
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
 
#ifndef __SRA_ARGUMENTS_H__
#define __SRA_ARGUMENTS_H__

#include "SRACore.h"
#include "SRA2BWTMdl.h"

typedef struct SRAQueryResultCount
{
    unsigned long long Shared_ReadPositions[MAX_CE_THRESHOLD];
    int Shared_AccMismatch[MAX_CE_THRESHOLD];
    char Shared_AccQualities[MAX_CE_THRESHOLD];

    //IsOpened      1 : Yes
    //              0 : No
    int IsOpened;

    //IsClosed      1 : Yes
    //              0 : No
    int IsClosed;

    //IsUnique      1 : Yes
    //              0 : No
    int IsUnique;

    //Best alignment position
    //IsBest        0 : Unknown
    //              1 : Exists
    int IsBest;

    unsigned long long BestAlignmentPosition;
    int BestAlignmentStrand;
    int BestAlignmentNumOfError;
    char BestAlignmentQuality;

    //WithError stores the number of occurrences found with the corresponding number of mismatch/edits
    unsigned int WithError[ MAX_NUM_OF_ERROR + 1];
    unsigned int TotalOccurrences;

    //Position retrieved by SA, CE or HOCC
    unsigned int RetrievedBySa;
    unsigned int RetrievedByCe;
    unsigned int RetrievedByHocc;
} SRAQueryResultCount;


typedef struct SRAQueryInput
{
    SRAQueryInfo * QueryInfo;
    SRASetting * QuerySetting;
    SRAQueryResultCount * QueryOutput;
    SRAIndex * AlgnmtIndex;
} SRAQueryInput;

#endif
