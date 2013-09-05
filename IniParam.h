/*
 *
 *    IniParam.c
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

#ifndef __INI_PARAM_H__
#define __INI_PARAM_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "2bwt-lib/dictionary.h"
#include "2bwt-lib/iniparser.h"
#include "definitions.h"
#include "2bwt-flex/SRA2BWTMdl.h"
#include "FileUtilities.h"
#include "UsageInterface.h"

// define a structure for storing the inputs for pair-multi model
typedef struct MultiInputItem
{
    char queryFile1[MAX_FIELD_LEN];
    char queryFile2[MAX_FIELD_LEN];
    int insert_low;
    int insert_high;
    char outputPrefix[MAX_FIELD_LEN];
    char readGrpID[MAX_FIELD_LEN];
    char sampleName[MAX_FIELD_LEN];
    char readGrpOpt[MAX_FIELD_LEN];
} MultiInputItem;

typedef struct InputOptions
{
    char * indexName;
    char * queryFileName;
    char * queryFileName2;
    int maxReadLength;

    int outputFormat;
    int isOutputBinary;

    // Three cases:
    // Succinct: outputFormat = 1, isOutputBinary = 0;
    // SAM: outputFormat = 2, isOutputBinary = 0;
    // BAM: outputFormat = 2, isOutputBinary = 1;

    char * outputPrefix;

    char readType; // 1: single; 2: pair-end
    char isReadList; // 0: each input file is a read file; 1: each input file contains a list of names of read files.
    char isReadBAM; // 1: the input file is in BAM format; 0: otherwise.

    int numMismatch; // the maximum number of mismatch allowed for the alignment of input reads
    int alignmentType; // 1: All valid alignment; 2: all best alignment

    int insert_low;
    int insert_high;

    char enableDP; // 1: enable; 0: disable dynamic programming for the unmapped mate

    int maxHitNum; // for default DP (even read)
    int maxHitNum2; // for default DP (odd read)


    char isIlluminaQual; // 0: Pred+33; 1: Pred+64

    int GPUDeviceID; // GPU device ID

    char * readGroup; // read group

    char * sampleName; // sample name

    char * readGrpOption; // read group option

    bool isPrintMDNM; // print MD string and NM tag

} InputOptions;


typedef struct IniParams
{
    char Ini_SaValueFileExt[MAX_FILEEXT_LEN];
    int  Ini_NumOfCpuThreads;
    char Ini_HostAlignmentModelStr[4];
    int  Ini_HostAlignmentModel;
    int  Ini_GPUMemory;
    int  Ini_PEStrandLeftLeg;
    int  Ini_PEStrandRightLeg;
    unsigned int Ini_MaxOutputPerRead;
    unsigned int Ini_PEMaxOutputPerPair;
    unsigned int Ini_MaxHitsEachEndForPairing;
    int Ini_MatchScore;
    int Ini_MismatchScore;
    int Ini_GapOpenScore;
    int Ini_GapExtendScore;
    int Ini_DPScoreThreshold;
    int Ini_isDefaultThreshold;
    int Ini_Soap3MisMatchAllow;
    int Ini_maxMAPQ;
    int Ini_minMAPQ;
    int Ini_shareIndex;
    int Ini_maxReadNameLen;
    int Ini_maxFrontLenClipped;
    int Ini_maxEndLenClipped;
    int Ini_proceedDPForTooManyHits;
    int Ini_skipSOAP3Alignment;
    int Ini_bwaLikeScore;
} IniParams;


// This function is to load the input file
MultiInputItem * loadMultiInputFile ( char * fileName, int isPair, int isBAM, int & lineNum );


// This function replace the special character \n or \t with a newline or a tab
void updateUserInputText ( char * str );

// This function is to parse the ini file and collect the paramaters of
// 1. The file extension of SA file
// 2. The number of CPU threads
// 3. The alignment model
// 4. The memory size in the GPU card
int ParseIniFile ( char * iniFileName, IniParams & ini_params );


// This function is to check and parse the input arguments
// arguments: <program> single bwtCodeIndex queryFileName numQueries maxReadLength [options]
// OR         <program> pair bwtCodeIndex queryFileName1 queryFileName2 numQueries maxReadLength [options]
// If the argument is not correct, then it returns FALSE. Else it returns TRUE
bool parseInputArgs ( int argc, char ** argv, InputOptions & input_options );


// print out the parameters
void printParameters ( InputOptions input_options, IniParams ini_params );

// print out the parameters
void printDPParameters ( DPParameters dp_params );

// This function is to update the values inside InputOptions
// according to the i-th set of values inside the array of MultiInputItem
void updateInputOption ( InputOptions * input_options, MultiInputItem * multiInputArray, int i );

#endif
