/*
 *
 *    SAList.c
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

#include "SAList.h"

SAList * SAListConstruct ()
{
    SAList * sa_list;
    sa_list = ( SAList * ) malloc ( sizeof ( SAList ) );
    sa_list->sa = ( PESRAAlignmentResult * ) malloc ( sizeof ( PESRAAlignmentResult ) * INITIAL_SIZE );
    sa_list->curr_size = 0;
    sa_list->available_size = INITIAL_SIZE;
    return sa_list;
}

void SAListReset ( SAList * sa_list )
{
    sa_list->curr_size = 0;
}

void SAListFree ( SAList * sa_list )
{
    free ( sa_list->sa );
    free ( sa_list );
}

void addSAToSAList ( SAList * sa_list, unsigned int l_value, unsigned int r_value, unsigned char strand, unsigned char mismatchCount )
{
    if ( l_value <= r_value && r_value != -1 )
    {
        if ( sa_list->curr_size >= sa_list->available_size )
        {
            // increase the size of the array by double
            unsigned int new_size = sa_list->available_size * 2;
            PESRAAlignmentResult * new_sa = ( PESRAAlignmentResult * ) malloc ( sizeof ( PESRAAlignmentResult ) * new_size );
            memcpy ( new_sa, sa_list->sa, sizeof ( PESRAAlignmentResult ) *sa_list->available_size );
            free ( sa_list->sa );
            sa_list->sa = new_sa;
            sa_list->available_size = new_size;
        }

        // add the values to the array
        sa_list->sa[sa_list->curr_size].saIndexLeft = l_value;
        sa_list->sa[sa_list->curr_size].saIndexRight = r_value;
        sa_list->sa[sa_list->curr_size].strand = strand;
        sa_list->sa[sa_list->curr_size].mismatchCount = mismatchCount;
        sa_list->curr_size++;
    }
}

unsigned int totalNumOcc ( SAList * sa_list )
{
    // total number of occurrences
    unsigned int sum = 0;

    for ( unsigned int i = 0; i < sa_list->curr_size; i++ )
    {
        sum += sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;
    }

    return sum;
}

char tooManyNumOcc ( SAList * sa_list, unsigned int maxOcc )
{
    // return true if the number of occurences is more than "maxOcc"
    unsigned int sum = 0;

    for ( unsigned int i = 0; i < sa_list->curr_size && sum <= maxOcc; i++ )
    {
        sum += sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;
    }

    return ( sum > maxOcc ) ? 1 : 0;
}

unsigned int numOccWithCap ( SAList * sa_list, unsigned int maxOcc )
{
    unsigned int sum = 0;

    for ( unsigned int i = 0; i < sa_list->curr_size && sum <= maxOcc; i++ )
    {
        sum += sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;
    }

    return sum;
}

int getMinMatchAndNumAns ( SAList * sa_list, OCCList * occ_list, int * totNumAns )
{
    // obtain the minimum number of mismatches and total number of answers for sa_list and occ_list
    int minMatch = 999;
    ( *totNumAns ) = 0;

    // scan SA List
    for ( int i = 0; i < sa_list->curr_size; i++ )
    {
        if ( sa_list->sa[i].mismatchCount < minMatch )
        {
            minMatch = sa_list->sa[i].mismatchCount;
        }

        ( *totNumAns ) += sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;
    }

    // scan occ List
    for ( int i = 0; i < occ_list->curr_size; i++ )
    {
        if ( occ_list->occ[i].mismatchCount < minMatch )
        {
            minMatch = occ_list->occ[i].mismatchCount;
        }

        ( *totNumAns ) ++;
    }

    return minMatch;
}

unsigned int retainAllBest ( SAList * sa_list, OCCList * occ_list )
{
    // retain all best answers
    int minMatch = 999;
    int newSAListSize = 0;
    int newOccListSize = 0;
    unsigned int num = 0;

    // scan SA List
    for ( int i = 0; i < sa_list->curr_size; i++ )
    {
        if ( sa_list->sa[i].mismatchCount < minMatch )
        {
            minMatch = sa_list->sa[i].mismatchCount;
            newSAListSize = 0;
            num = sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;
            sa_list->sa[newSAListSize].readID = sa_list->sa[i].readID;
            sa_list->sa[newSAListSize].saIndexLeft = sa_list->sa[i].saIndexLeft;
            sa_list->sa[newSAListSize].saIndexRight = num + sa_list->sa[i].saIndexLeft - 1;
            sa_list->sa[newSAListSize].strand = sa_list->sa[i].strand;
            sa_list->sa[newSAListSize].mismatchCount = sa_list->sa[i].mismatchCount;
            newSAListSize++;
        }
        else if ( sa_list->sa[i].mismatchCount == minMatch )
        {
            int curr_num = sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;
            sa_list->sa[newSAListSize].readID = sa_list->sa[i].readID;
            sa_list->sa[newSAListSize].saIndexLeft = sa_list->sa[i].saIndexLeft;
            sa_list->sa[newSAListSize].saIndexRight = curr_num + sa_list->sa[i].saIndexLeft - 1;
            sa_list->sa[newSAListSize].strand = sa_list->sa[i].strand;
            sa_list->sa[newSAListSize].mismatchCount = sa_list->sa[i].mismatchCount;
            num += curr_num;
            newSAListSize++;
        }
    }

    // scan occ List
    for ( int i = 0; i < occ_list->curr_size; i++ )
    {
        if ( occ_list->occ[i].mismatchCount < minMatch )
        {
            minMatch = occ_list->occ[i].mismatchCount;
            newSAListSize = 0;
            newOccListSize = 0;
            occ_list->occ[newOccListSize].readID = occ_list->occ[i].readID;
            occ_list->occ[newOccListSize].ambPosition = occ_list->occ[i].ambPosition;
            occ_list->occ[newOccListSize].strand = occ_list->occ[i].strand;
            occ_list->occ[newOccListSize].mismatchCount = occ_list->occ[i].mismatchCount;
            num = 1;
            newOccListSize++;
        }
        else if ( occ_list->occ[i].mismatchCount == minMatch )
        {
            occ_list->occ[newOccListSize].readID = occ_list->occ[i].readID;
            occ_list->occ[newOccListSize].ambPosition = occ_list->occ[i].ambPosition;
            occ_list->occ[newOccListSize].strand = occ_list->occ[i].strand;
            occ_list->occ[newOccListSize].mismatchCount = occ_list->occ[i].mismatchCount;
            num++;
            newOccListSize++;
        }
    }

    sa_list->curr_size = newSAListSize;
    occ_list->curr_size = newOccListSize;
    // fprintf(stderr, "%u\n", num);
    return num;
}


unsigned int retainAllBestWithCap ( SAList * sa_list, OCCList * occ_list, int max_num )
{
    // retain all best answers
    int minMatch = 999;
    int newSAListSize = 0;
    int newOccListSize = 0;
    unsigned int num = 0;

    // scan SA List
    for ( int i = 0; i < sa_list->curr_size; i++ )
    {
        if ( sa_list->sa[i].mismatchCount < minMatch )
        {
            minMatch = sa_list->sa[i].mismatchCount;
            newSAListSize = 0;
            int curr_num = sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;

            if ( curr_num > max_num )
            {
                curr_num = max_num;
            }

            sa_list->sa[newSAListSize].readID = sa_list->sa[i].readID;
            sa_list->sa[newSAListSize].saIndexLeft = sa_list->sa[i].saIndexLeft;
            sa_list->sa[newSAListSize].saIndexRight = curr_num + sa_list->sa[i].saIndexLeft - 1;
            sa_list->sa[newSAListSize].strand = sa_list->sa[i].strand;
            sa_list->sa[newSAListSize].mismatchCount = sa_list->sa[i].mismatchCount;
            num = curr_num;
            newSAListSize++;
        }
        else if ( sa_list->sa[i].mismatchCount == minMatch && num < max_num )
        {
            int curr_num = sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;

            if ( num + curr_num > max_num )
            {
                curr_num = max_num - num;
            }

            sa_list->sa[newSAListSize].readID = sa_list->sa[i].readID;
            sa_list->sa[newSAListSize].saIndexLeft = sa_list->sa[i].saIndexLeft;
            sa_list->sa[newSAListSize].saIndexRight = curr_num + sa_list->sa[i].saIndexLeft - 1;
            sa_list->sa[newSAListSize].strand = sa_list->sa[i].strand;
            sa_list->sa[newSAListSize].mismatchCount = sa_list->sa[i].mismatchCount;
            num += curr_num;
            newSAListSize++;
        }
    }

    // scan occ List
    for ( int i = 0; i < occ_list->curr_size; i++ )
    {
        if ( occ_list->occ[i].mismatchCount < minMatch )
        {
            minMatch = occ_list->occ[i].mismatchCount;
            newSAListSize = 0;
            newOccListSize = 0;
            occ_list->occ[newOccListSize].readID = occ_list->occ[i].readID;
            occ_list->occ[newOccListSize].ambPosition = occ_list->occ[i].ambPosition;
            occ_list->occ[newOccListSize].strand = occ_list->occ[i].strand;
            occ_list->occ[newOccListSize].mismatchCount = occ_list->occ[i].mismatchCount;
            num = 1;
            newOccListSize++;
        }
        else if ( occ_list->occ[i].mismatchCount == minMatch && num < max_num )
        {
            occ_list->occ[newOccListSize].readID = occ_list->occ[i].readID;
            occ_list->occ[newOccListSize].ambPosition = occ_list->occ[i].ambPosition;
            occ_list->occ[newOccListSize].strand = occ_list->occ[i].strand;
            occ_list->occ[newOccListSize].mismatchCount = occ_list->occ[i].mismatchCount;
            num++;
            newOccListSize++;
        }
    }

    sa_list->curr_size = newSAListSize;
    occ_list->curr_size = newOccListSize;
    // fprintf(stderr, "%u\n", num);
    return num;
}

unsigned int retainAllBestAndSecBest ( SAList * sa_list, OCCList * occ_list )
{
    // retain all best (k-mismatch) and second best (i.e. k+1-mismatch) answers
    int minMatch = 999;
    int newSAListSize = 0;
    int newOccListSize = 0;
    unsigned int num = 0;

    // first scan SA List and OCC List and obtain the minMatch
    for ( int i = 0; i < sa_list->curr_size; i++ )
    {
        if ( sa_list->sa[i].mismatchCount < minMatch )
        {
            minMatch = sa_list->sa[i].mismatchCount;
        }
    }

    for ( int i = 0; i < occ_list->curr_size; i++ )
    {
        if ( occ_list->occ[i].mismatchCount < minMatch )
        {
            minMatch = occ_list->occ[i].mismatchCount;
        }
    }

    // then collect the best and second best hits for SA List
    for ( int i = 0; i < sa_list->curr_size; i++ )
    {
        if ( sa_list->sa[i].mismatchCount <= minMatch + 1 )
        {
            int curr_num = sa_list->sa[i].saIndexRight - sa_list->sa[i].saIndexLeft + 1;
            sa_list->sa[newSAListSize].readID = sa_list->sa[i].readID;
            sa_list->sa[newSAListSize].saIndexLeft = sa_list->sa[i].saIndexLeft;
            sa_list->sa[newSAListSize].saIndexRight = curr_num + sa_list->sa[i].saIndexLeft - 1;
            sa_list->sa[newSAListSize].strand = sa_list->sa[i].strand;
            sa_list->sa[newSAListSize].mismatchCount = sa_list->sa[i].mismatchCount;
            num += curr_num;
            newSAListSize++;
        }
    }

    // scan occ List
    for ( int i = 0; i < occ_list->curr_size; i++ )
    {
        if ( occ_list->occ[i].mismatchCount <= minMatch + 1 )
        {
            occ_list->occ[newOccListSize].readID = occ_list->occ[i].readID;
            occ_list->occ[newOccListSize].ambPosition = occ_list->occ[i].ambPosition;
            occ_list->occ[newOccListSize].strand = occ_list->occ[i].strand;
            occ_list->occ[newOccListSize].mismatchCount = occ_list->occ[i].mismatchCount;
            num++;
            newOccListSize++;
        }
    }

    sa_list->curr_size = newSAListSize;
    occ_list->curr_size = newOccListSize;
    return num;
}

OCCList * OCCListConstruct ()
{
    OCCList * occ_list;
    occ_list = ( OCCList * ) malloc ( sizeof ( OCCList ) );
    occ_list->occ = ( SRAOccurrence * ) malloc ( sizeof ( SRAOccurrence ) * INITIAL_SIZE );
    occ_list->curr_size = 0;
    occ_list->available_size = INITIAL_SIZE;
    return occ_list;
}

void OCCListReset ( OCCList * occ_list )
{
    occ_list->curr_size = 0;
}

void addToOCCList ( OCCList * occ_list, unsigned int pos, char strand, char mismatchCount )
{
    if ( occ_list->curr_size >= occ_list->available_size )
    {
        // increase the size of the array by double
        unsigned int new_size = occ_list->available_size * 2;
        SRAOccurrence * new_occ = ( SRAOccurrence * ) malloc ( sizeof ( SRAOccurrence ) * new_size );
        memcpy ( new_occ, occ_list->occ, sizeof ( SRAOccurrence ) *occ_list->available_size );
        free ( occ_list->occ );
        occ_list->occ = new_occ;
        occ_list->available_size = new_size;
    }

    // add the values to the array
    occ_list->occ[occ_list->curr_size].ambPosition = pos;
    occ_list->occ[occ_list->curr_size].strand = strand;
    occ_list->occ[occ_list->curr_size].mismatchCount = mismatchCount;
    occ_list->curr_size++;
}

void OCCListFree ( OCCList * occ_list )
{
    free ( occ_list->occ );
    free ( occ_list );
}


void transferAllSAToOcc ( SAList * sa_list, OCCList * occ_list, BWT * bwt )
{
    // transfer all the SA ranges in sa_list to occ_list
    // then reset sa_list
    unsigned int i, k, l, r;
    int strand;
    char num_mis;

    if ( sa_list->curr_size > 0 )
    {
        for ( i = 0; i < sa_list->curr_size; i++ )
        {
            l = sa_list->sa[i].saIndexLeft;
            r = sa_list->sa[i].saIndexRight;
            strand = sa_list->sa[i].strand;
            num_mis = sa_list->sa[i].mismatchCount;

            for ( k = l; k <= r; k++ )
            {
                addToOCCList ( occ_list, ( *bwt->_bwtSaValue ) ( bwt, k ), strand, num_mis );
            }
        }
    }

    // reset sa_list
    SAListReset ( sa_list );
}

