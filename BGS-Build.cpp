/*
 *
 *    SOAP3-Builder.cpp
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "2bwt-lib/HSPstatistic.h"

#include "Release.h"

int PoolSize = 2097152;                // 2M  - fixed; not configurable through ini

void printUsage()
{

    printf("\n[Main] %s 2BWT Index Enhancer v%d.%d.%d (%s) - Usage guide:\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
    printf("\n");
    printf("Syntax: BGS-Build <2bwtIndex>\n");
    printf("        2bwtIndex : The bwt index built by %s.\n", PROJECT_BUILDER_BINARY );
}

void reAlignFile ( const char * fn, const char * to )
{
    size_t i;
    char * buf = ( char * ) malloc ( 1UL << 10 + 16 );
    FILE * fhi = fopen ( fn, "r" );
    fread ( buf, sizeof ( unsigned int ), 5, fhi );
    fread ( buf, sizeof ( char ), 12, fhi );

    if ( strcmp ( ECOR1, buf ) == 0 )
    {
        fprintf ( stderr, "%s modified already.\n", fn );
    }
    else
    {
        FILE * fho = fopen ( "realignTmp", "w-" );
        fseek ( fhi, 0L, SEEK_END );
        long pos = ftell ( fhi );
        //fprintf(stderr, "Original size of %s: %ld\n", fn, pos);
        fseek ( fhi, 0L, SEEK_SET );
        fread ( buf, sizeof ( unsigned int ), 5, fhi );
        fwrite ( buf, sizeof ( unsigned int ), 5, fho );
        fwrite ( ECOR1, sizeof ( char ), 12, fho );

        while ( !feof ( fhi ) )
        {
            i = fread ( ( void * ) buf, 1, 1UL << 10, fhi );
            fwrite ( buf, 1, i, fho );
        }

        pos = ftell ( fho );
        //fprintf(stderr, "New size of %s: %ld\n", fn, pos);
        fclose ( fho );
        rename ( "realignTmp", to );
    }

    fclose ( fhi );
    free ( buf );
}

int main ( int argc, char ** argv )
{
    if ( argc != 2 )
    {
        fprintf ( stderr, "Invalid number of command-line arguments\n" );
        printUsage();
        return 1;
    }

    char * packedDNAFileName = ( char * ) malloc ( strlen ( argv[1] ) + 5 );
    char * bwtCodeFileName = ( char * ) malloc ( strlen ( argv[1] ) + 5 );
    char * occValueFileName = ( char * ) malloc ( strlen ( argv[1] ) + 5 );
    char * revBwtCodeFileName = ( char * ) malloc ( strlen ( argv[1] ) + 15 );
    char * revOccValueFileName = ( char * ) malloc ( strlen ( argv[1] ) + 15 );
    char * newPackedDNAFileName = ( char * ) malloc ( strlen ( argv[1] ) + 25 );
    char * newOccValueFileName = ( char * ) malloc ( strlen ( argv[1] ) + 25 );
    char * newRevOccValueFileName = ( char * ) malloc ( strlen ( argv[1] ) + 35 );
    char * gpuOccValueFileName = ( char * ) malloc ( strlen ( argv[1] ) + 25 );
    char * revGpuOccValueFileName = ( char * ) malloc ( strlen ( argv[1] ) + 35 );
    strcpy ( packedDNAFileName, argv[1] );
    strcpy ( packedDNAFileName + strlen ( argv[1] ), ".pac" );
    strcpy ( bwtCodeFileName, argv[1] );
    strcpy ( bwtCodeFileName + strlen ( argv[1] ), ".bwt" );
    strcpy ( occValueFileName, argv[1] );
    strcpy ( occValueFileName + strlen ( argv[1] ), ".fmv" );
    strcpy ( revBwtCodeFileName, argv[1] );
    strcpy ( revBwtCodeFileName + strlen ( argv[1] ), ".rev.bwt" );
    strcpy ( revOccValueFileName, argv[1] );
    strcpy ( revOccValueFileName + strlen ( argv[1] ), ".rev.fmv" );
    strcpy ( newPackedDNAFileName, argv[1] );
    strcpy ( newPackedDNAFileName + strlen ( argv[1] ), ".pac.mmap" );
    strcpy ( newOccValueFileName, argv[1] );
    strcpy ( newOccValueFileName + strlen ( argv[1] ), ".fmv.mmap" );
    strcpy ( newRevOccValueFileName, argv[1] );
    strcpy ( newRevOccValueFileName + strlen ( argv[1] ), ".rev.fmv.mmap" );
    strcpy ( gpuOccValueFileName, argv[1] );
    strcpy ( gpuOccValueFileName + strlen ( argv[1] ), ".fmv.gpu" );
    strcpy ( revGpuOccValueFileName, argv[1] );
    strcpy ( revGpuOccValueFileName + strlen ( argv[1] ), ".rev.fmv.gpu" );
    MMPool * mmPool;
    MMMasterInitialize ( 3, 0, FALSE, NULL );
    mmPool = MMPoolCreate ( PoolSize );
    BWT * bwt = BWTLoad ( mmPool, bwtCodeFileName, occValueFileName, NULL, NULL, NULL, NULL, 0 );
    BWT * revBwt = BWTLoad ( mmPool, revBwtCodeFileName, revOccValueFileName, NULL, NULL, NULL, NULL, 0 );
    // Start building the BGS-Index
    printf ( "[BGS-Build] Writing GPU occurrence table to %s..\n", gpuOccValueFileName );
    FILE * fout = ( FILE * ) fopen64 ( gpuOccValueFileName, "wb" );

    if ( fout == NULL )
    {
        fprintf ( stderr, "[BGS-Build] Cannot open output file!\n" );
        exit ( 1 );
    }

    fwrite ( &bwt->inverseSa0, sizeof ( unsigned int ), 1, fout );
    fwrite ( bwt->cumulativeFreq + 1, sizeof ( unsigned int ), ALPHABET_SIZE, fout );
    unsigned int numOfOccValue = ( bwt->textLength + 128 - 1 ) / 128 + 1;
    int i, j;
    unsigned int x;
    int adj = 0;

    for ( x = 0; x < numOfOccValue; x++ )
    {
        unsigned int sumTmp = 0;

        for ( j = 0; j < ALPHABET_SIZE; j++ )
        {
            if ( x * 128 > bwt->inverseSa0 ) {adj = 1;}

            unsigned int tmp = BWTOccValue ( bwt, x * 128 + adj, j );
            sumTmp += tmp;
            tmp += bwt->cumulativeFreq[j];
            fwrite ( &tmp, sizeof ( unsigned int ), 1, fout );
        }

        if ( sumTmp % 128 != 0 ) {fprintf ( stderr, "[BGS-Build] ERROR %u %u\n", x, sumTmp ); return 0;}
    }

    fclose ( fout );
    printf ( "[BGS-Build] Writing reversed GPU occurrence table to %s..\n", revGpuOccValueFileName );
    fout = ( FILE * ) fopen64 ( revGpuOccValueFileName, "wb" );

    if ( fout == NULL )
    {
        fprintf ( stderr, "[BGS-Build] Cannot open output file!\n" );
        exit ( 1 );
    }

    fwrite ( &revBwt->inverseSa0, sizeof ( unsigned int ), 1, fout );
    fwrite ( revBwt->cumulativeFreq + 1, sizeof ( unsigned int ), ALPHABET_SIZE, fout );
    adj = 0;

    for ( x = 0; x < numOfOccValue; x++ )
    {
        unsigned int sumTmp = 0;

        for ( j = 0; j < ALPHABET_SIZE; j++ )
        {
            if ( x * 128 > revBwt->inverseSa0 ) {adj = 1;}

            unsigned int tmp = BWTOccValue ( revBwt, x * 128 + adj, j );
            sumTmp += tmp;
            tmp += revBwt->cumulativeFreq[j];
            fwrite ( &tmp, sizeof ( unsigned int ), 1, fout );
        }

        if ( sumTmp % 128 != 0 ) {fprintf ( stderr, "[BGS-Build] ERROR %u %u\n", x, sumTmp ); return 0;}
    }

    fclose ( fout );
    printf ( "[BGS-Build] Building of BGS-Index has completed.\n" );
    BWTFree ( mmPool, bwt, 0 );
    BWTFree ( mmPool, revBwt, 0 );
    MMPoolFree ( mmPool );
    reAlignFile ( occValueFileName, newOccValueFileName );
    reAlignFile ( revOccValueFileName, newRevOccValueFileName );
    HSP * hsp = ( HSP * ) malloc ( sizeof ( HSP ) );
    hsp->packedDNA = ( unsigned int * ) DNALoadPacked ( packedDNAFileName, &hsp->dnaLength, TRUE, 1 );
    FILE * fho = fopen ( "realignTmp", "w-" );
    unsigned int trailerBufferInWord = 1;
    unsigned int wordToProcess = ( hsp->dnaLength + 64 - 1 ) / 64 * 4 + ( ( ( trailerBufferInWord + 3 ) / 4 * 4 ) * 4 );
    fwrite ( & ( hsp->dnaLength ), sizeof ( unsigned int ), 1, fho );
    fwrite ( ECOR1, sizeof ( char ), DNAPACK_HOLLOW_BETWEEN_METADATA_PAYLOAD, fho );
    fwrite ( hsp->packedDNA, sizeof ( unsigned int ), wordToProcess, fho );
    fclose ( fho );
    rename ( "realignTmp", newPackedDNAFileName );
    DNAFreePacked ( hsp->packedDNA, hsp->dnaLength, trailerBufferInWord, 0 );
    free ( hsp );
    fprintf ( stderr, "[BGS-Build] Building index files for direct memory mapping completed.\n" );
    printf ( "[BGS-Build] GPUBWT terminated.\n" );
    free ( packedDNAFileName );
    free ( bwtCodeFileName );
    free ( occValueFileName );
    free ( revBwtCodeFileName );
    free ( revOccValueFileName );
    free ( newPackedDNAFileName );
    free ( newOccValueFileName );
    free ( newRevOccValueFileName );
    free ( gpuOccValueFileName );
    free ( revGpuOccValueFileName );
    return 0;
}
