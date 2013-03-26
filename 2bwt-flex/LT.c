//
//    LT.c
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

#include "LT.h"

LT *  LTLoad(const char * fileName)  {
    LT * lookupTable = (LT*) malloc(sizeof(LT));
    FILE * fin = (FILE*)fopen64(fileName, "rb");
    if (fin==NULL) { 
        fprintf(stderr,"Cannot open lookupTableFile %s\n",fileName); 
        exit(1);
    }
    int tableSize;
    fread(&(tableSize),sizeof(int),1,fin);
    unsigned long long lookupWordSize = 1 << (tableSize * LOOKUP_BIT_PER_CHAR);
    lookupTable->table = (unsigned int*) malloc(sizeof(unsigned int) * lookupWordSize);
    lookupTable->ltSizeInWord = lookupWordSize;
    lookupTable->tableSize = tableSize;
    unsigned step = 1048576;
    unsigned long long i;
    for (i = 0; i < lookupWordSize; i += LOOKUP_LOAD_STEP) {
        fread(lookupTable->table+i, LOOKUP_LOAD_STEP * sizeof(unsigned int), 1, fin);
    }
    fclose(fin);
    return lookupTable;
}

void LTFree(LT * lookupTable) {
    free(lookupTable->table);
    free(lookupTable);
}

