/*
 *
 *    SAM.c
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

#include "SAM.h"

void SAMCreateCIGAR ( int * mis_pos, int num_mis, char * cigar )
{
    int i = 0;
    char stack[10];
    int j;
    int k;
    int tmp;
    int l = 0;

    for ( i = 0; i < num_mis; i++ )
    {
        tmp = mis_pos[i];
        j = 0;

        while ( tmp > 0 )
        {
            stack[j++] = '0' + tmp % 10;
        }

        for ( k = j - 1; k >= 0; k-- )
        {
            cigar[l++] = stack[k];
        }
    }

    cigar[l] = '\0';
}

int SAMIUint8ConcatUint8 ( uint8_t * data, int * curSize,
                           uint8_t key )
{
    data[ ( *curSize ) ++] = key;
}
int SAMIUint8ConcatUint32 ( uint8_t * data, int * curSize,
                            uint32_t key )
{
    int i;
    int len = sizeof ( uint32_t ) / sizeof ( uint8_t );
    uint8_t * key_8 = ( uint8_t * ) &key;

    for ( i = 0; i < len; i++ )
    {
        data[ ( *curSize ) ++] = key_8[i];
    }
}
int SAMIUint8ConcatString ( uint8_t * data, int * curSize,
                            char * key, int len )
{
    int i;

    for ( i = 0; i < len; i++ )
    {
        data[ ( *curSize ) ++] = key[i];
    }
}

void SAMOutputHeaderConstruct ( bam_header_t * sheader, HSP * hsp, HSPAux * hspaux, int maxReadLength )
{
    unsigned int * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    unsigned int tp, approxIndex, approxValue;
    int i, j;
    //Number of sequences
    sheader->n_targets = hsp->numOfSeq;
    sheader->target_name = ( char ** ) malloc ( sizeof ( char * ) * hsp->numOfSeq );
    sheader->target_len = ( uint32_t * ) malloc ( sizeof ( uint32_t ) * hsp->numOfSeq );

    char * snTags = ( char * ) malloc ( MAX_CHROMOSOME_NUM * MAX_SEQ_NAME_LENGTH * 2 ); // maximum chromosomes * (number of chars max each*2)
    memset ( snTags, '\0', MAX_CHROMOSOME_NUM * MAX_SEQ_NAME_LENGTH * 2 );
    char * snTagsTmp = ( char * ) malloc ( MAX_SEQ_NAME_LENGTH * 2 );
    memset ( snTagsTmp, '\0', MAX_SEQ_NAME_LENGTH * 2 );

    for ( i = 0; i < hsp->numOfSeq; i++ )
    {
        for ( j = 0; j < 255; j++ )
        {
            if ( hsp->annotation[i].text[j] == '\0' ||
                    hsp->annotation[i].text[j] == ' ' ||
                    hsp->annotation[i].text[j] == '\t' ||
                    hsp->annotation[i].text[j] == '\r' ||
                    hsp->annotation[i].text[j] == '\n' )
            {
                break;
            }
        }

        hsp->annotation[i].text[j] = '\0';
        //Fill in the Ref.Seq names and lengths.
        sheader->target_name[i] = hsp->annotation[i].text;
        sheader->target_len[i] = hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos;
        sprintf ( snTagsTmp, "@SQ\tSN:%s\tLN:%u\n", hsp->annotation[i].text, hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos );
        strcat ( snTags, snTagsTmp );
    }

    char programInfo[1024];
    sprintf ( programInfo, "@PG\tID:%s\tPN:%s\tVN:v%d.%d.%d (%s)\n", PROJECT_NAME, PROJECT_NAME, PROJECT_MAJOR, PROJECT_MINOR, PROJECT_REV, PROJECT_SPECIAL );
    int textLen = strlen ( hspaux->readGroup ) + strlen ( hspaux->sampleName ) + strlen ( hspaux->readGrpOption ) + strlen ( programInfo ) + strlen ( snTags ) + 15 + 23;
    sheader->text = ( char * ) malloc ( textLen * sizeof ( char ) );
    sprintf ( sheader->text, "@HD\tVN:1.3\tSO:unsorted\n@RG\tID:%s\tSM:%s\t%s\n%s%s", hspaux->readGroup, hspaux->sampleName, hspaux->readGrpOption, snTags, programInfo );
    sheader->l_text = strlen ( sheader->text );
    //Given up unknown parameters.
    //If someone ever found out what the hell is this pls update.
    free ( snTags );
    free ( snTagsTmp );
    sheader->hash = NULL;
    sheader->rg2lib = NULL;
}

void SAMOutputHeaderDestruct ( bam_header_t * sheader )
{
    free ( sheader->target_name );
    free ( sheader->target_len );
    free ( sheader->text );
}


void SAMOccurrenceConstruct ( OCC * occ )
{
    occ->SAMOutBuffer.data = ( uint8_t * ) malloc ( sizeof ( uint8_t ) * SAM_MDATA_SIZE );
}

void SAMOccurrenceDestruct ( OCC * occ )
{
    free ( occ->SAMOutBuffer.data );
}
