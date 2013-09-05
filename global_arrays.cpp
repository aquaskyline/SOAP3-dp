/*
 *
 *    global_arrays.cpp
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


#include "global_arrays.h"

// construct the array
AlgnResult * arrayConstruct ()
{
    AlgnResult * newArray;
    newArray = ( AlgnResult * ) malloc ( sizeof ( AlgnResult ) );
    newArray->occ_list = ( occRec * ) malloc ( sizeof ( occRec ) * INITIAL_SIZE_GLOBAL_ARRAY );
    newArray->occTotalNum = 0;
    newArray->availableSize = INITIAL_SIZE_GLOBAL_ARRAY;
    return newArray;
}

// free the array
void arrayFree ( AlgnResult * algnResult )
{
    free ( algnResult->occ_list );
    free ( algnResult );
}

// reset the array
void arrayReset ( AlgnResult * algnResult )
{
    algnResult->occTotalNum = 0;
}

// construct the array
AlgnResultArrays * resultArraysConstruct ( int num )
{
    AlgnResultArrays * newArray;
    newArray = ( AlgnResultArrays * ) malloc ( sizeof ( AlgnResultArrays ) );
    newArray->algnArrays = ( AlgnResult ** ) malloc ( sizeof ( AlgnResult * ) * num );

    for ( int i = 0; i < num; i++ )
    {
        newArray->algnArrays[i] = arrayConstruct ();
    }

    newArray->numArrays = num;
    return newArray;
}

// free the array
void resultArraysFree ( AlgnResultArrays * algnResultArrays )
{
    for ( int i = 0; i < algnResultArrays->numArrays; i++ )
    {
        arrayFree ( algnResultArrays->algnArrays[i] );
    }

    free ( algnResultArrays->algnArrays );
}

// reset the array
void resultArraysReset ( AlgnResultArrays * algnResultArrays )
{
    for ( int i = 0; i < algnResultArrays->numArrays; i++ )
    {
        arrayReset ( algnResultArrays->algnArrays[i] );
    }
}



// add occ to the array
void addOCCToArray ( AlgnResult * algnResult, unsigned int readID, unsigned int ambPosition, unsigned char strand,
                     unsigned char source, char score )
{
    if ( algnResult->occTotalNum >= algnResult->availableSize )
    {
        // increase the size of the array by double
        unsigned int new_size = algnResult->availableSize * 2;
        occRec * new_occ = ( occRec * ) malloc ( sizeof ( occRec ) * new_size );
        memcpy ( new_occ, algnResult->occ_list, sizeof ( occRec ) *algnResult->occTotalNum );
        free ( algnResult->occ_list );
        algnResult->occ_list = new_occ;
        algnResult->availableSize = new_size;
    }

    // add the values to the array
    algnResult->occ_list[algnResult->occTotalNum].readID = readID;
    algnResult->occ_list[algnResult->occTotalNum].ambPosition = ambPosition;
    algnResult->occ_list[algnResult->occTotalNum].strand = strand;
    algnResult->occ_list[algnResult->occTotalNum].source = source;
    algnResult->occ_list[algnResult->occTotalNum].score = score;
    algnResult->occTotalNum++;
}

