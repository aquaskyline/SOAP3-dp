/*
 *
 *    IndexHandler.c
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

#include "IndexHandler.h"

void INDEXFillCharMap ( unsigned char charMap[255] )
{

    int i;

    for ( i = 0; i < 255; i++ )
    {
        charMap[i] = 0;
    }

    //Replacing DNA_CHAR_SIZE with ALPHABET_SIZE as SOAP3 has been developed as
    // a restrictive tool..
    for ( i = 0; i < ALPHABET_SIZE; i++ )
    {
        charMap[dnaChar[i]] = ( unsigned char ) i;
        charMap[dnaChar[i] - 'A' + 'a'] = ( unsigned char ) i;
    }

    charMap['U'] = charMap['T'];
    charMap['N'] = charMap['G']; // N -> G

    charMap['u'] = charMap['t'];
    charMap['n'] = charMap['g']; // N -> G
}

// To load the index
Soap3Index * INDEXLoad ( IniParams * ini_params, char * indexName, char isShareIndex )
{
    MMMasterInitialize ( 3, 0, FALSE, NULL );
    printf ( "[Main] loading index into host...\n" );
    double startTime = setStartTime ();
    Soap3Index * index = ( Soap3Index * ) malloc ( sizeof ( Soap3Index ) );
    SRAIndex * sraIndex = ( SRAIndex * ) malloc ( sizeof ( SRAIndex ) );

    index->sraIndex = sraIndex;
    index->mmPool = MMPoolCreate ( INDEX_MMPOOL_SIZE );
    IndexFileNames indexFilenames;
    INDEXProcessFilenames ( &indexFilenames, indexName, ini_params );
    INDEXLoad ( index->mmPool, indexFilenames, & ( sraIndex->bwt ), & ( sraIndex->rev_bwt ), & ( sraIndex->hsp ), & ( sraIndex->lookupTable ), & ( sraIndex->rev_lookupTable ),
                & ( index->gpu_revOccValue ), & ( index->gpu_occValue ), index->gpu_numOfOccValue, isShareIndex );


    /////////////////////////////////
    // Allocate and Fill CharMap
    index->charMap = ( unsigned char * ) malloc ( sizeof ( unsigned char ) * 256 );
    INDEXFillCharMap ( index->charMap );

    /////////////////////////////////
    // Allocate HSPAuxilliary
    index->sraIndex->hspaux = ( HSPAux * ) malloc ( sizeof ( HSPAux ) );

    double loadIndexTime = getElapsedTime ( startTime );
    printf ( "[Main] Finished loading index into host.\n" );
    printf ( "[Main] Loading time : %9.4f seconds\n\n", loadIndexTime );
    return index;
}


// To free the index
void INDEXFree ( Soap3Index * index, char isShareIndex )
{
    SRAIndex * sraIndex = index->sraIndex;
    uint numOfOccValue = ( sraIndex->bwt->textLength + GPU_OCC_INTERVAL - 1 ) / GPU_OCC_INTERVAL + 1;
    BWTFree ( index->mmPool, sraIndex->bwt, isShareIndex );
    BWTFree ( index->mmPool, sraIndex->rev_bwt, isShareIndex );
    HSPFree ( index->mmPool, sraIndex->hsp, 1, isShareIndex );

    if ( isShareIndex )
    {
        munmap ( sraIndex->lookupTable->table, sizeof ( uint ) * ( sraIndex->lookupTable->ltSizeInWord + 1 ) );
        munmap ( sraIndex->rev_lookupTable->table, sizeof ( uint ) * ( sraIndex->rev_lookupTable->ltSizeInWord + 1 ) );
        free ( sraIndex->lookupTable );
        free ( sraIndex->rev_lookupTable );
    }
    else
    {
        LTFree ( sraIndex->lookupTable );
        LTFree ( sraIndex->rev_lookupTable );
    }

    MMPoolFree ( index->mmPool );
    free ( index->sraIndex->hspaux );

    if ( isShareIndex )
    {
        munmap ( index->gpu_revOccValue, sizeof ( unsigned int ) * ( numOfOccValue * ALPHABET_SIZE + GPU_OCC_PAYLOAD_OFFSET ) );
        munmap ( index->gpu_occValue, sizeof ( unsigned int ) * ( numOfOccValue * ALPHABET_SIZE + GPU_OCC_PAYLOAD_OFFSET ) );
    }
    else
    {
        free ( index->gpu_revOccValue );
        free ( index->gpu_occValue );
    }

    free ( index->charMap );
    free ( sraIndex );
    free ( index );
}

void INDEXLoad ( MMPool * mmPool, IndexFileNames indexFilenames, BWT ** bwt, BWT ** revBwt, HSP ** hsp, LT ** lkt, LT ** revLkt,
                 uint ** revOccValue, uint ** occValue, uint & numOfOccValue,
                 char isShareIndex )
{
    if ( isShareIndex )
    { ( *bwt ) = BWTLoad ( mmPool, indexFilenames.bwtCodeFileName, indexFilenames.mmapOccValueFileName, indexFilenames.saCodeFileName, NULL, NULL, NULL, isShareIndex ); }
    else
    { ( *bwt ) = BWTLoad ( mmPool, indexFilenames.bwtCodeFileName, indexFilenames.occValueFileName, indexFilenames.saCodeFileName, NULL, NULL, NULL, isShareIndex ); }

    free ( indexFilenames.bwtCodeFileName );
    free ( indexFilenames.occValueFileName );
    free ( indexFilenames.saCodeFileName );
    free ( indexFilenames.mmapOccValueFileName );

    if ( isShareIndex )
    { ( *revBwt ) = BWTLoad ( mmPool, indexFilenames.revBwtCodeFileName, indexFilenames.mmapRevOccValueFileName, NULL, NULL, NULL, NULL, isShareIndex ); }
    else
    { ( *revBwt ) = BWTLoad ( mmPool, indexFilenames.revBwtCodeFileName, indexFilenames.revOccValueFileName, NULL, NULL, NULL, NULL, isShareIndex ); }

    free ( indexFilenames.revBwtCodeFileName );
    free ( indexFilenames.revOccValueFileName );
    free ( indexFilenames.mmapRevOccValueFileName );

    if ( isShareIndex )
    { ( *hsp ) = HSPLoad ( mmPool, indexFilenames.mmapPackedDnaFileName, indexFilenames.annotationFileName, indexFilenames.ambiguityFileName, indexFilenames.translateFileName, 1, isShareIndex ); }
    else
    { ( *hsp ) = HSPLoad ( mmPool, indexFilenames.packedDnaFileName, indexFilenames.annotationFileName, indexFilenames.ambiguityFileName, indexFilenames.translateFileName, 1, isShareIndex ); }

    free ( indexFilenames.packedDnaFileName );
    free ( indexFilenames.mmapPackedDnaFileName );
    free ( indexFilenames.annotationFileName );
    free ( indexFilenames.ambiguityFileName );
    free ( indexFilenames.translateFileName );
    numOfOccValue = ( ( *bwt )->textLength + GPU_OCC_INTERVAL - 1 ) / GPU_OCC_INTERVAL + 1; // Value at both end for bi-directional encoding
    unsigned int lookupWordSize = 1 << ( LOOKUP_SIZE * 2 );
    FILE * gpuOccValueFile, *revGpuOccValueFile, *lookupTableFile, *revLookupTableFile;
    // GPU OCC TABLE
    gpuOccValueFile = ( FILE * ) fopen ( indexFilenames.gpuOccValueFileName, "rb" );

    if ( gpuOccValueFile == NULL ) { fprintf ( stderr, "Cannot open gpuOccValueFile %s\n", indexFilenames.gpuOccValueFileName ); exit ( 1 );}

    free ( indexFilenames.gpuOccValueFileName );
    // consume useless parts of OCC file
    unsigned int tmp;
    fread ( &tmp, sizeof ( unsigned int ), 1, gpuOccValueFile );

    // read cumulative frequency
    for ( uint i = 1; i <= ALPHABET_SIZE; i++ )
    { fread ( &tmp, sizeof ( unsigned int ), 1, gpuOccValueFile ); }

    // read OCC from file. copied from BWT.c
    if ( isShareIndex == 1 )
    {
        //fprintf ( stderr, "Using MMAP GPU OCC\n" );
        int gpuOccValueFD = fileno ( gpuOccValueFile );
        ( *occValue ) = ( unsigned int * ) mmap ( NULL, ( sizeof ( unsigned int ) * ( numOfOccValue * ALPHABET_SIZE + GPU_OCC_PAYLOAD_OFFSET ) ), PROT_READ, MAP_SHARED | MAP_NORESERVE, gpuOccValueFD, 0 );
        int rtVal;
        rtVal = mlock ( ( *occValue ), ( sizeof ( unsigned int ) * ( numOfOccValue * ALPHABET_SIZE + GPU_OCC_PAYLOAD_OFFSET ) ) );

        if ( ( void * ) rtVal == MAP_FAILED )
        {
            fprintf ( stderr, "\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:" );
            perror ( "" );
            exit ( EXIT_FAILURE );
        }

        ( *occValue ) += GPU_OCC_PAYLOAD_OFFSET;
    }
    else
    {
        //fprintf ( stderr, "Using malloc GPU OCC\n" );
        ( *occValue ) = ( unsigned int * ) malloc ( numOfOccValue * ALPHABET_SIZE * sizeof ( unsigned int ) );
        fread ( ( *occValue ), sizeof ( unsigned int ), numOfOccValue * ALPHABET_SIZE, gpuOccValueFile );
    }

    fclose ( gpuOccValueFile );

    if ( isShareIndex == 1 )
    {
        // LOOK UP TABLE
        lookupTableFile = ( FILE * ) fopen ( indexFilenames.lookupTableFileName, "rb" );

        if ( lookupTableFile == NULL ) { fprintf ( stderr, "Cannot open lookupTableFile %s\n", indexFilenames.lookupTableFileName ); exit ( 1 );}

        // read lookup table
        int tableSize;
        fread ( &tableSize, sizeof ( unsigned int ), 1, lookupTableFile );

        lookupWordSize = 1 << tableSize * BIT_PER_CHAR;
        //fprintf ( stderr, "Using MMAP LKT\n" );
        int lookupTableFD = fileno ( lookupTableFile );
        ( *lkt ) = ( LT * ) malloc ( sizeof ( LT ) );
        ( *lkt )->table = ( unsigned int * ) mmap ( NULL, ( sizeof ( uint ) * ( lookupWordSize + 1 ) ), PROT_READ, MAP_SHARED | MAP_NORESERVE, lookupTableFD, 0 );
        ( *lkt )->ltSizeInWord = lookupWordSize;
        ( *lkt )->tableSize = tableSize;
        int rtVal;
        rtVal = mlock ( ( *lkt )->table, ( sizeof ( uint ) * ( lookupWordSize + 1 ) ) );

        if ( ( void * ) rtVal == MAP_FAILED )
        {
            fprintf ( stderr, "\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:" );
            perror ( "" );
            exit ( EXIT_FAILURE );
        }

        ( *lkt )->table += 1;
        fclose ( lookupTableFile );
    }
    else
    {
        ( *lkt ) = LTLoad ( indexFilenames.lookupTableFileName );
    }

    if ( ( *lkt )->tableSize != LOOKUP_SIZE )
    {
        fprintf ( stderr, "SOAP3 only supports lookup table size of %u, please re-build index.\n", LOOKUP_SIZE );
        exit ( 1 );
    }

    free ( indexFilenames.lookupTableFileName );
    /////////////////// READ REVERSE ////////////////////////////////////
    // GPU OCC TABLE
    revGpuOccValueFile = ( FILE * ) fopen ( indexFilenames.revGpuOccValueFileName, "rb" );

    if ( revGpuOccValueFile == NULL ) { fprintf ( stderr, "Cannot open revGpuOccValueFile\n" ); exit ( 1 );}

    free ( indexFilenames.revGpuOccValueFileName );
    // consume useless parts of OCC file
    fread ( &tmp, sizeof ( unsigned int ), 1, revGpuOccValueFile );

    // read cumulative frequency
    for ( uint i = 1; i <= ALPHABET_SIZE; i++ )
    { fread ( &tmp, sizeof ( unsigned int ), 1, revGpuOccValueFile ); }

    // read OCC from file. copied from BWT.c
    if ( isShareIndex == 1 )
    {
        //fprintf ( stderr, "Using MMAP rev GPU OCC\n" );
        int revGpuOccValueFD = fileno ( revGpuOccValueFile );
        ( *revOccValue ) = ( unsigned int * ) mmap ( NULL, ( sizeof ( unsigned int ) * ( numOfOccValue * ALPHABET_SIZE + GPU_OCC_PAYLOAD_OFFSET ) ), PROT_READ, MAP_SHARED | MAP_NORESERVE, revGpuOccValueFD, 0 );
        int rtVal;
        rtVal = mlock ( ( *revOccValue ), ( sizeof ( unsigned int ) * ( numOfOccValue * ALPHABET_SIZE + GPU_OCC_PAYLOAD_OFFSET ) ) );

        if ( ( void * ) rtVal == MAP_FAILED )
        {
            fprintf ( stderr, "\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:" );
            perror ( "" );
            exit ( EXIT_FAILURE );
        }

        ( *revOccValue ) += GPU_OCC_PAYLOAD_OFFSET;
    }
    else
    {
        //fprintf ( stderr, "Using malloc rev GPU OCC\n" );
        ( *revOccValue ) = ( uint * ) malloc ( numOfOccValue * ALPHABET_SIZE * sizeof ( uint ) );
        fread ( ( *revOccValue ), sizeof ( uint ), numOfOccValue * ALPHABET_SIZE, revGpuOccValueFile );
    }

    fclose ( revGpuOccValueFile );

    if ( isShareIndex == 1 )
    {
        // LOOK UP TABLE
        revLookupTableFile = ( FILE * ) fopen ( indexFilenames.revLookupTableFileName, "rb" );

        if ( revLookupTableFile == NULL ) { fprintf ( stderr, "Cannot open revLookupTableFile\n" ); exit ( 1 );}

        // read lookup table
        int tableSize;
        fread ( &tableSize, sizeof ( unsigned int ), 1, revLookupTableFile );

        lookupWordSize = 1 << tableSize * BIT_PER_CHAR;
        //fprintf ( stderr, "Using MMAP rev LKT\n" );
        int revLookupTableFD = fileno ( revLookupTableFile );
        ( *revLkt ) = ( LT * ) malloc ( sizeof ( LT ) );
        ( *revLkt )->table = ( unsigned int * ) mmap ( NULL, ( sizeof ( uint ) * ( lookupWordSize + 1 ) ), PROT_READ, MAP_SHARED | MAP_NORESERVE, revLookupTableFD, 0 );
        ( *revLkt )->ltSizeInWord = lookupWordSize;
        ( *revLkt )->tableSize = tableSize;
        int rtVal;
        rtVal = mlock ( ( *revLkt )->table, ( sizeof ( uint ) * ( lookupWordSize + 1 ) ) );

        if ( ( void * ) rtVal == MAP_FAILED )
        {
            fprintf ( stderr, "\nMemory allocation error.\nPlease check:\n1. You've got enough memory left.\n2. Kernel memlock size was set as 'unlimited'.\n\tType 'unlimit -l' to check the value.\n\tType 'unlimit -l unlimited' to solve the problem.\n\tUse root account to run command \"echo '* - memlock unlimited' >> /etc/security/limits.conf\" for permanent solution.\nErrno:" );
            perror ( "" );
            exit ( EXIT_FAILURE );
        }

        ( *revLkt )->table += 1;
        fclose ( revLookupTableFile );
    }
    else
    {
        ( *revLkt ) = LTLoad ( indexFilenames.revLookupTableFileName );
    }

    if ( ( *revLkt )->tableSize != LOOKUP_SIZE )
    {
        fprintf ( stderr, "SOAP3 only supports lookup table size of %u, please re-build index.\n", LOOKUP_SIZE );
        exit ( 1 );
    }

    free ( indexFilenames.revLookupTableFileName );
    // SHOW THE MEMORY USAGE ON HOST
#ifdef BGS_OUTPUT_MEMORY_USAGE
    // for BWT
    printf ( "BWT code size    : %lu bytes (%9.4f M)\n", ( *bwt )->bwtSizeInWord * sizeof ( unsigned int ), ( float ) ( *bwt )->bwtSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    printf ( "Occ value size   : %lu bytes (%9.4f M)\n", ( ( *bwt )->occSizeInWord + ( *bwt )->occMajorSizeInWord ) * sizeof ( unsigned int ), ( float ) ( ( *bwt )->occSizeInWord + ( *bwt )->occMajorSizeInWord ) * sizeof ( unsigned int ) / 1024.0 / 1024.0 );

    if ( ( *bwt )->saValueSizeInWord > 0 )
    {
        printf ( "SA value size    : %lu bytes (%9.4f M)\n", ( *bwt )->saValueSizeInWord * sizeof ( unsigned int ), ( float ) ( *bwt )->saValueSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    if ( ( *bwt )->inverseSaSizeInWord > 0 )
    {
        printf ( "Inverse SA size  : %lu bytes (%9.4f M)\n", ( *bwt )->inverseSaSizeInWord * sizeof ( unsigned int ), ( float ) ( *bwt )->inverseSaSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    if ( ( *bwt )->cachedSaIndex > 0 )
    {
        printf ( "SA index range  : %lu bytes (%9.4f M)\n", ( *bwt )->cachedSaIndexSizeInWord * sizeof ( unsigned int ), ( float ) ( *bwt )->cachedSaIndexSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    // for REVBWT
    printf ( "Rev BWT code size    : %lu bytes (%9.4f M)\n", ( *revBwt )->bwtSizeInWord * sizeof ( unsigned int ), ( float ) ( *revBwt )->bwtSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    printf ( "Rev Occ value size   : %lu bytes (%9.4f M)\n", ( ( *revBwt )->occSizeInWord + ( *revBwt )->occMajorSizeInWord ) * sizeof ( unsigned int ), ( float ) ( ( *revBwt )->occSizeInWord + ( *revBwt )->occMajorSizeInWord ) * sizeof ( unsigned int ) / 1024.0 / 1024.0 );

    if ( ( *revBwt )->saValueSizeInWord > 0 )
    {
        printf ( "Rev SA value size    : %lu bytes (%9.4f M)\n", ( *revBwt )->saValueSizeInWord * sizeof ( unsigned int ), ( float ) ( *revBwt )->saValueSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    if ( ( *revBwt )->inverseSaSizeInWord > 0 )
    {
        printf ( "Rev Inverse SA size  : %lu bytes (%9.4f M)\n", ( *revBwt )->inverseSaSizeInWord * sizeof ( unsigned int ), ( float ) ( *revBwt )->inverseSaSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    if ( ( *revBwt )->cachedSaIndex > 0 )
    {
        printf ( "Rev SA index range  : %lu bytes (%9.4f M)\n", ( *revBwt )->cachedSaIndexSizeInWord * sizeof ( unsigned int ), ( float ) ( *revBwt )->cachedSaIndexSizeInWord * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    }

    // for HSP
    uint trailerBufferInWord = 1;
    uint trailerBufferIn128 = ( trailerBufferInWord + 3 ) / 4 * 4;
    uint textLength = ( *hsp )->dnaLength;
    // hsp->packedDNA
    uint hsp_mem_usage = ( ( textLength + 64 - 1 ) / 64 * 4 + trailerBufferIn128 * 4 ) * sizeof ( unsigned int );
    // hsp->seqOffset
    hsp_mem_usage += ( ( *hsp )->numOfSeq + 1 ) * sizeof ( SeqOffset );
    // hsp->annotation
    hsp_mem_usage += ( ( *hsp )->numOfSeq + 1 ) * sizeof ( Annotation );
    // hsp->ambiguity
    hsp_mem_usage += ( ( *hsp )->numOfAmbiguity + 2 ) * sizeof ( Ambiguity );
    printf ( "[Memory usage in Host] hsp: %u bytes (%9.4f M)\n", hsp_mem_usage, ( float ) hsp_mem_usage / 1024.0 / 1024.0 );
    printf ( "[Memory usage in Host] occValue: %i bytes (%9.4f M)\n", numOfOccValue * ALPHABET_SIZE * sizeof ( unsigned int ), ( float ) numOfOccValue * ALPHABET_SIZE * sizeof ( unsigned int ) / 1024.0 / 1024.0 );
    printf ( "[Memory usage in Host] lkt: %i bytes (%9.4f M)\n", sizeof ( unsigned int ) * lookupWordSize, ( float ) sizeof ( unsigned int ) * lookupWordSize / 1024.0 / 1024.0 );
    printf ( "[Memory usage in Host] revOccValue: %i bytes (%9.4f M)\n", numOfOccValue * ALPHABET_SIZE * sizeof ( uint ), ( float ) numOfOccValue * ALPHABET_SIZE * sizeof ( uint ) / 1024.0 / 1024.0 );
    printf ( "[Memory usage in Host] revLkt: %i bytes (%9.4f M)\n", sizeof ( uint ) * lookupWordSize, ( float ) sizeof ( uint ) * lookupWordSize / 1024.0 / 1024.0 );
    printf ( "[Memory usage in Host] occ for all cpu threads: %i bytes (%9.4f M)\n", sizeof ( OCC ) *ini_params.Ini_NumOfCpuThreads, ( float ) sizeof ( OCC ) *ini_params.Ini_NumOfCpuThreads / 1024.0 / 1024.0 );
#endif
}



void INDEXProcessFilenames ( IndexFileNames * indexFilenames, char * indexName, IniParams * ini_params )
{
    char saExtension[10] = ".sa";

    if ( ini_params != NULL )
    {
        strcpy ( saExtension, ini_params->Ini_SaValueFileExt );
    }

    int indexNameLen = strlen ( indexName );
    indexFilenames->iniFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->bwtCodeFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->occValueFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->gpuOccValueFileName = ( char * ) malloc ( indexNameLen + 25 );
    indexFilenames->lookupTableFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->revBwtCodeFileName = ( char * ) malloc ( indexNameLen + 15 );
    indexFilenames->revOccValueFileName = ( char * ) malloc ( indexNameLen + 15 );
    indexFilenames->revGpuOccValueFileName = ( char * ) malloc ( indexNameLen + 35 );
    indexFilenames->revLookupTableFileName = ( char * ) malloc ( indexNameLen + 15 );
    indexFilenames->saCodeFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->memControlFileName = ( char * ) malloc ( indexNameLen + 5 );
    //For HSP
    indexFilenames->packedDnaFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->annotationFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->ambiguityFileName = ( char * ) malloc ( indexNameLen + 5 );
    indexFilenames->translateFileName = ( char * ) malloc ( indexNameLen + 5 );
    // For mmap
    indexFilenames->mmapOccValueFileName = ( char * ) malloc ( indexNameLen + 25 );
    indexFilenames->mmapRevOccValueFileName = ( char * ) malloc ( indexNameLen + 35 );
    indexFilenames->mmapPackedDnaFileName = ( char * ) malloc ( indexNameLen + 25 );
    strcpy ( indexFilenames->bwtCodeFileName, indexName );
    strcpy ( indexFilenames->bwtCodeFileName + indexNameLen, ".bwt" );
    strcpy ( indexFilenames->occValueFileName, indexName );
    strcpy ( indexFilenames->occValueFileName + indexNameLen, ".fmv" );
    strcpy ( indexFilenames->lookupTableFileName, indexName );
    strcpy ( indexFilenames->lookupTableFileName + indexNameLen, ".lkt" );
    strcpy ( indexFilenames->gpuOccValueFileName, indexName );
    strcpy ( indexFilenames->gpuOccValueFileName + indexNameLen, ".fmv.gpu" );
    strcpy ( indexFilenames->revBwtCodeFileName, indexName );
    strcpy ( indexFilenames->revBwtCodeFileName + indexNameLen, ".rev.bwt" );
    strcpy ( indexFilenames->revOccValueFileName, indexName );
    strcpy ( indexFilenames->revOccValueFileName + indexNameLen, ".rev.fmv" );
    strcpy ( indexFilenames->revLookupTableFileName, indexName );
    strcpy ( indexFilenames->revLookupTableFileName + indexNameLen, ".rev.lkt" );
    strcpy ( indexFilenames->revGpuOccValueFileName, indexName );
    strcpy ( indexFilenames->revGpuOccValueFileName + indexNameLen, ".rev.fmv.gpu" );
    strcpy ( indexFilenames->saCodeFileName, indexName );
    strcpy ( indexFilenames->saCodeFileName + indexNameLen, saExtension );
    strcpy ( indexFilenames->packedDnaFileName, indexName );
    strcpy ( indexFilenames->packedDnaFileName + indexNameLen, ".pac" );
    strcpy ( indexFilenames->annotationFileName, indexName );
    strcpy ( indexFilenames->annotationFileName + indexNameLen, ".ann" );
    strcpy ( indexFilenames->ambiguityFileName, indexName );
    strcpy ( indexFilenames->ambiguityFileName + indexNameLen, ".amb" );
    strcpy ( indexFilenames->translateFileName, indexName );
    strcpy ( indexFilenames->translateFileName + indexNameLen, ".tra" );
    strcpy ( indexFilenames->mmapOccValueFileName, indexName );
    strcpy ( indexFilenames->mmapOccValueFileName + indexNameLen, ".fmv.mmap" );
    strcpy ( indexFilenames->mmapRevOccValueFileName, indexName );
    strcpy ( indexFilenames->mmapRevOccValueFileName + indexNameLen, ".rev.fmv.mmap" );
    strcpy ( indexFilenames->mmapPackedDnaFileName, indexName );
    strcpy ( indexFilenames->mmapPackedDnaFileName + indexNameLen, ".pac.mmap" );
#ifdef BGS_OUTPUT_FILENAMES_TO_SCREEN
    printf ( "[Main] Using BWT index ....                %s\n", bwtCodeFileName );
    printf ( "[Main] Using Occ index ....                %s\n", occValueFileName );
    printf ( "[Main] Using Lookup index ....             %s\n", lookupTableFileName );
    printf ( "[Main] Using Gpu-Occ index ....            %s\n", gpuOccValueFileName );
    printf ( "[Main] Using Rev-BWT index ....            %s\n", revBwtCodeFileName );
    printf ( "[Main] Using Rev-Occ index ....            %s\n", revOccValueFileName );
    printf ( "[Main] Using Rev-Lookup index ....         %s\n", revLookupTableFileName );
    printf ( "[Main] Using Rev-Gpu-Occ index ....        %s\n", revGpuOccValueFileName );
    printf ( "[Main] Using SA index ....                 %s\n", saCodeFileName );
    printf ( "[Main] Using PackedDna index ....          %s\n", packedDnaFileName );
    printf ( "[Main] Using Annotation index ....         %s\n", annotationFileName );
    printf ( "[Main] Using Ambiguity index ....          %s\n", ambiguityFileName );
    printf ( "[Main] Using Translation index ....        %s\n", translateFileName );
#endif
}
