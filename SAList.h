/*
 *
 *    SAList.h
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

#ifndef _SALIST_H_
#define _SALIST_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PEAlgnmt.h"
#include "definitions.h"

#define INITIAL_SIZE 10240

typedef struct SAList
{
    PESRAAlignmentResult * sa;
    unsigned int curr_size;
    unsigned int available_size;
} SAList;


typedef struct OCCList
{
    SRAOccurrence * occ;
    unsigned int curr_size;
    unsigned int available_size;
} OCCList;


SAList * SAListConstruct ();
void SAListReset ( SAList * sa_list );
void addSAToSAList ( SAList * sa_list, unsigned int l_value, unsigned int r_value, unsigned char strand, unsigned char mismatchCount );
void SAListFree ( SAList * sa_list );
unsigned int totalNumOcc ( SAList * sa_list );
char tooManyNumOcc ( SAList * sa_list, unsigned int maxOcc );
unsigned int numOccWithCap ( SAList * sa_list, unsigned int maxOcc );
int getMinMatchAndNumAns ( SAList * sa_list, OCCList * occ_list, int * totNumAns );
// obtain the minimum number of mismatches and total number of answers for sa_list and occ_list
unsigned int retainAllBest ( SAList * sa_list, OCCList * occ_list );
unsigned int retainAllBestWithCap ( SAList * sa_list, OCCList * occ_list, int max_num );
// retain all best answers
unsigned int retainAllBestAndSecBest ( SAList * sa_list, OCCList * occ_list );

OCCList * OCCListConstruct ();
void OCCListReset ( OCCList * occ_list );
void addToOCCList ( OCCList * occ_list, unsigned int pos, char strand, char mismatchCount );
void OCCListFree ( OCCList * occ_list );

void transferAllSAToOcc ( SAList * sa_list, OCCList * occ_list, BWT * bwt );

#endif
