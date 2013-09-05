/*
 *
 *    IndexHandler.h
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

#ifndef __INDEX_HANDLER_H__
#define __INDEX_HANDLER_H__

#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>

#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "2bwt-flex/LT.h"
#include "2bwt-lib/Timing.h"
#include "2bwt-lib/dictionary.h"
#include "2bwt-lib/iniparser.h"
#include "HSPAux.h"

#include "definitions.h"
#include "IniParam.h"

#define INDEX_MMPOOL_SIZE 2097152                // 2M  - fixed; not configurable through ini 

#define GPU_OCC_PAYLOAD_OFFSET ( 1 + ALPHABET_SIZE )

typedef struct Soap3Index
{
    //CPU Index
    MMPool * mmPool;
    unsigned char * charMap;

    //GPU Index
    uint gpu_numOfOccValue;
    uint * gpu_occValue;
    uint * gpu_revOccValue;

    SRAIndex * sraIndex;

} Soap3Index;

typedef struct IndexFileNames
{
    char * iniFileName;
    char * bwtCodeFileName;
    char * occValueFileName;
    char * gpuOccValueFileName;
    char * lookupTableFileName;
    char * revBwtCodeFileName;
    char * revOccValueFileName;
    char * revGpuOccValueFileName;
    char * revLookupTableFileName;
    char * saCodeFileName;
    char * memControlFileName;

    //For HSP
    char * packedDnaFileName;
    char * annotationFileName;
    char * ambiguityFileName;
    char * translateFileName;

    // For mmap use
    char * mmapOccValueFileName;
    char * mmapRevOccValueFileName;
    char * mmapPackedDnaFileName;

} IndexFileNames;

void INDEXFillCharMap ( unsigned char charMap[255] );

// process index file name
void INDEXProcessFilenames ( IndexFileNames * index, char * indexName, IniParams * ini_params );

void INDEXLoad ( MMPool * mmPool, IndexFileNames indexFilenames, BWT ** bwt, BWT ** revBwt, HSP ** hsp, LT ** lkt, LT ** revLkt,
                 uint ** revOccValue, uint ** occValue, uint & numOfOccValue,
                 char isShareIndex );

// To load the index
Soap3Index * INDEXLoad ( IniParams * ini_params, char * indexName, char isShareIndex );

// To free the index
void INDEXFree ( Soap3Index * index, char isShareIndex );

#endif
