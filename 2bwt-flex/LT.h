//
//    LT.h
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
   
   Date   : 8th May 2011
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __LOOKUPTABLE_H__
#define __LOOKUPTABLE_H__

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>

// LOOKUP_BIT_PER_CHAR defines as the actual bit per char
// used by lookup table.
#define LOOKUP_BIT_PER_CHAR   2
#define LOOKUP_CHAR_MASK      3

#define LOOKUP_LOAD_STEP      1048576

#define LOOKUP_SIZE 13

typedef struct LT {
       int tableSize;
       unsigned int ltSizeInWord;
       unsigned int * table;
} LT;


LT *  LTLoad(const char * fileName);

void LTFree(LT * lookupTable);

#endif

