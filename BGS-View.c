/*
 *
 *    BGS-View.c
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
#include "2bwt-lib/MiscUtilities.h"
#include "2bwt-lib/MemManager.h"
#include "2bwt-lib/TextConverter.h"
#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "2bwt-lib/HSPstatistic.h"
#include "BGS-IO.h"

#define BGS_VIEWER_QUERY_NAME_LENGTH            60
#define BGS_VIEWER_INPUT_BUFFER_SIZE            1048576 // 1024 x 1024

//Define the below parameter to output the alignment result(text position)
// with the old 2BWT-Aligner format
//#define DEBUG_2BWT_OUTPUT_32_BIT

//Define the debug flag to turn on extra output.
// #define DEBUG_2BWT

char AlignmentTypeString[2][5] = {"MISM", "EDIT"};

int GetStringLength ( char * str, int max )
{
    int i;
    int len = 0;

    for ( i = 0; i < max; i++ )
    {
        if ( str[i] == '\0' )
        {
            return len;
        }

        len++;
    }

    return max;
}

int main ( int argc, char ** argv )
{
    if ( argc < 3 )
    {
        fprintf ( stderr, "\nSOAP3-GPU Output Viewer/Convertor v1.0 - Usage guide:\n" );
        fprintf ( stderr, "Windows Platform:\n\t BGS-View <query_file> <output_file>\n" );
        fprintf ( stderr, "Linux Platform:\n\t ./BGS-View <query_file> <output_file>\n" );
        fprintf ( stderr, "Note: <output_file> must be the corresponding output file of running soap3_aligner with <query_file>\n\n" );
        return 0;
    }

#ifdef DEBUG_2BWT_OUTPUT_32_BIT
    fprintf ( stderr, "[OCCDebug] Expect output in 32-bit structures.\n" );
#endif
    char * readFileName = argv[1];
    char * outputFileName = argv[2];
    unsigned int i, j, k;
    //Open the necessary files for decoding
    // fprintf(stderr, "Opening %s .. ",readFileName);
    FILE * queryFile = ( FILE * ) fopen ( readFileName, "r" );

    if ( queryFile == NULL || queryFile == 0 )
    {
        fprintf ( stderr, "Opening %s failed!\n", readFileName );
        return 0;
    }

    //Open the necessary files for decoding
    // fprintf(stderr, "Opening %s ..",outputFileName);
    FILE * outputFile = ( FILE * ) fopen64 ( outputFileName, "rb" );

    if ( outputFile == NULL || outputFile == 0 )
    {
        fprintf ( stderr, "Opening %s failed!\n", outputFileName );
        return 0;
    }

    fseek ( outputFile, -1, SEEK_END );
    unsigned long long outputFileLength = ftell ( outputFile );
    outputFileLength++;
    // fprintf(stderr,"File size = %llu\n",outputFileLength);
    fseek ( outputFile, 0, SEEK_SET );
    //Load output version
    unsigned int outputVersion;

    if ( fread ( &outputVersion, 1, sizeof ( unsigned int ), outputFile ) != sizeof ( unsigned int ) )
    {
        fprintf ( stderr, "Bad file content in %s!\n", outputFileName );
        return 0;
    }

#ifdef DEBUG_2BWT
    fprintf ( stderr, "Output file version = %u\n", outputVersion );
#endif

    if ( outputVersion != OCC_OUTPUT_FORMAT )
    {
        fprintf ( stderr, "Unsupported output format %u expects %u\n", outputVersion, OCC_OUTPUT_FORMAT );
    }

    //Load read length from output file.
    int maxReadLength;

    if ( fread ( &maxReadLength, 1, sizeof ( int ), outputFile ) != sizeof ( int ) )
    {
        fprintf ( stderr, "Bad file content in %s!\n", outputFileName );
        return 0;
    }

#ifdef DEBUG_2BWT
    fprintf ( stderr, "Max read length = %u\n", maxReadLength );
#endif
    //Load read count from output file.
    unsigned int numQueries;

    if ( fread ( &numQueries, 1, sizeof ( int ), outputFile ) != sizeof ( int ) )
    {
        fprintf ( stderr, "Bad file content in %s!\n", outputFileName );
        return 0;
    }

#ifdef DEBUG_2BWT
    fprintf ( stderr, "Number of queries in the read file = %u\n", numQueries );
#endif
    //Load Annotation from output file
    unsigned int genomeNumOfSequence;

    if ( fread ( &genomeNumOfSequence, 1, sizeof ( unsigned int ), outputFile ) != sizeof ( unsigned int ) )
    {
        fprintf ( stderr, "Bad file content in %s!\n", outputFileName );
        return 0;
    }

    unsigned long long * genomeLength = ( unsigned long long * ) malloc ( sizeof ( unsigned long long ) * genomeNumOfSequence );
    Annotation * genomeAnnotation = ( Annotation * ) malloc ( sizeof ( Annotation ) * genomeNumOfSequence );

    for ( i = 0; i < genomeNumOfSequence; i++ )
    {
        if ( fread ( & ( genomeAnnotation[i].gi ), 1, sizeof ( int ), outputFile ) != sizeof ( int ) ) { return 0; }

        if ( fread ( & ( genomeAnnotation[i].text ), 1, sizeof ( char ) * ( MAX_SEQ_NAME_LENGTH + 1 ), outputFile ) != sizeof ( char ) * ( MAX_SEQ_NAME_LENGTH + 1 ) ) { return 0; }

        for ( j = 0; j < MAX_SEQ_NAME_LENGTH + 1; j++ )
        {
            if ( genomeAnnotation[i].text[j] == '\n' ) { break; }

            if ( genomeAnnotation[i].text[j] == ' ' ) { break; }

            if ( genomeAnnotation[i].text[j] == '\0' ) { break; }

            genomeAnnotation[i].decoratedText[j] = genomeAnnotation[i].text[j];
        }

        genomeAnnotation[i].decoratedText[j] = '\0';

        if ( fread ( & ( genomeLength[i] ), 1, sizeof ( unsigned int ), outputFile ) != sizeof ( unsigned int ) ) { return 0; }
    }

    unsigned long long fileHeaderSize = ftell ( outputFile );
#ifdef DEBUG_2BWT
    fprintf ( stderr, "Output file header size = %llu\n", ( unsigned long long ) fileHeaderSize );
    fprintf ( stderr, "Output file record size = %llu\n", ( unsigned long long ) sizeof ( OCCPositionCacheToDisk ) );
#endif

    //Check file body size
    if ( ( outputFileLength - fileHeaderSize ) % sizeof ( OCCPositionCacheToDisk ) > 0 )
    {
        fprintf ( stderr, "Bad file format in %s.\n[Header:%llu Content:%llu LeftOver:%llu]!\n",
                  outputFileName,
                  fileHeaderSize,
                  outputFileLength - fileHeaderSize,
                  ( outputFileLength - fileHeaderSize ) % sizeof ( OCCPositionCacheToDisk ) );
        return 0;
    }

    char * upkdQueries = ( char * ) malloc ( ( long long int ) numQueries * ( maxReadLength + 1 ) * sizeof ( char ) ); // a large array to store all queries
    char * upkdQueriesName = ( char * ) malloc ( ( long long int ) numQueries * ( BGS_VIEWER_QUERY_NAME_LENGTH + 1 ) * sizeof ( char ) ); // a large array to store all queries
    long long int queriesRead = 0;
    char * queryFileBuffer = ( char * ) malloc ( sizeof ( char ) * BGS_VIEWER_INPUT_BUFFER_SIZE );
    size_t bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
    unsigned int bufferIndex = 0;
    char queryChar = queryFileBuffer[bufferIndex++];
    int isFastq = 0;

    while ( bufferSize != 0 )
    {
        //Read everything before the entry point of a read, character ">"
        while ( bufferSize != 0 && queryChar != '>' && queryChar != '@' )
        {
            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }
        }

        isFastq = 0;

        if ( queryChar == '@' ) {isFastq = 1;}

        //Read the header of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
            bufferIndex = 0;
        }

        int headerCharId = 0;
        int cropped = 0;

        while ( bufferSize != 0  && queryChar != '\n' )
        {
            if ( queryChar == ' ' || queryChar == '\t' ) { cropped = 1; }

            if ( !cropped && headerCharId < BGS_VIEWER_QUERY_NAME_LENGTH )
            { upkdQueriesName[queriesRead * ( BGS_VIEWER_QUERY_NAME_LENGTH + 1 ) + ( headerCharId++ )] = queryChar; }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }
        }

        upkdQueriesName[queriesRead * ( BGS_VIEWER_QUERY_NAME_LENGTH + 1 ) + ( headerCharId++ )] = '\0';

        //Read the pattern body of a read
        if ( bufferSize == 0 ) { break; }

        queryChar = queryFileBuffer[bufferIndex++];

        if ( bufferIndex >= bufferSize )
        {
            bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
            bufferIndex = 0;
        }

        int nucleoId = 0;

        while ( bufferSize != 0  && queryChar != '>' && queryChar != '@' && queryChar != '+' )
        {
            if ( queryChar != '\n' )
            {
                if ( nucleoId < maxReadLength )
                {
                    upkdQueries[queriesRead * ( maxReadLength + 1 ) + nucleoId] = queryChar;
                }

                nucleoId++;
            }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }
        }

        if ( nucleoId < maxReadLength )
        { upkdQueries[queriesRead * ( maxReadLength + 1 ) + nucleoId] = '\0'; }
        else
        { upkdQueries[queriesRead * ( maxReadLength + 1 ) + maxReadLength] = '\0'; }

        //Parse the quality part of the read
        if ( isFastq )
        {
            //Read everything before the entry point of a read, character "+"
            while ( bufferSize != 0 && queryChar != '+' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                    bufferIndex = 0;
                }
            }

            //Read the header of a read
            if ( bufferSize == 0 ) { break; }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }

            while ( bufferSize != 0 && queryChar != '\n' )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                    bufferIndex = 0;
                }
            }

            //Read the quality of a read
            if ( bufferSize == 0 ) { break; }

            queryChar = queryFileBuffer[bufferIndex++];

            if ( bufferIndex >= bufferSize )
            {
                bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                bufferIndex = 0;
            }

            for ( int i = 0; i < nucleoId; i++ )
            {
                queryChar = queryFileBuffer[bufferIndex++];

                if ( bufferIndex >= bufferSize )
                {
                    bufferSize = fread ( queryFileBuffer, sizeof ( char ), BGS_VIEWER_INPUT_BUFFER_SIZE, queryFile );
                    bufferIndex = 0;
                }

                //In SOAP3-GPU/SOAP3-CPU we do not have handling for read quality.
                //therefore the read quality is CURRENTLY read from the input file
                //and discarded immediately.
            }
        }

        queriesRead++;

        if ( queriesRead >= numQueries ) { break; }
    }

    free ( queryFileBuffer );
    fclose ( queryFile );
    // fprintf(stderr,"[Main] %lld short reads loaded from query file\n", queriesRead); numQueries=queriesRead;
    unsigned int outputNumberOfRecords = ( outputFileLength - fileHeaderSize ) / sizeof ( OCCPositionCacheToDisk );
    //char queryName[MAX_SEQ_NAME_LENGTH];
    //char querySequence[QUERY_LENGTH];
    unsigned querySequenceLength = 0;
    OCCPositionCacheToDisk outputResultBuffer;
    unsigned int numberOfAlignments = 0, currentQueryRecord = 0;
    unsigned char cacheFlag1, cacheFlag2;
    unsigned long long cacheFlag3;
    unsigned char ChromId, ReadStrand;
    unsigned long long ReadId, tp;
    char * alignedRead;
    char * alignedReadHeader;
    i = 0;

    // printf("#Read_Pos\tRead_Name\tRead\tChromosome\tOffset\tRead_Length\tStrand\tType");
    while ( i < outputNumberOfRecords )
    {
        fread ( &outputResultBuffer, 1, sizeof ( OCCPositionCacheToDisk ), outputFile );
        cacheFlag1 = outputResultBuffer.cell[0];
        cacheFlag2 = outputResultBuffer.cell[1];
        cacheFlag3 = ( * ( unsigned int * ) &outputResultBuffer.cell[3] );
#ifdef DEBUG_2BWT
        // fprintf(stderr,"\nRecord is read : %u %u %llu\n",cacheFlag1, cacheFlag2, cacheFlag3);
#endif

        if ( cacheFlag2 == 0 )
        {
            ReadId = cacheFlag3;
            //if (ReadId>currentQueryRecord)
            //    currentQueryRecord = GetQueryFromFile(queryFile,currentQueryRecord,ReadId-currentQueryRecord-1,queryName,querySequence,&querySequenceLength);
            alignedRead = upkdQueries + ( ReadId - 1 ) * ( maxReadLength + 1 );
            alignedReadHeader = upkdQueriesName + ( ReadId - 1 ) * ( BGS_VIEWER_QUERY_NAME_LENGTH + 1 );
            querySequenceLength = GetStringLength ( alignedRead, maxReadLength );
        }
        else
        {
            //get read
            ChromId = cacheFlag2 - 1;
            tp = cacheFlag3;
            ReadStrand = cacheFlag1;
            printf ( "%llu\t%s\t%s\t", ReadId, alignedReadHeader, alignedRead );

            if ( genomeAnnotation[ChromId].gi > 0 ) { printf ( "%u|", genomeAnnotation[ChromId].gi ); }

            printf ( "%s\t%llu\t%u\t", genomeAnnotation[ChromId].decoratedText, tp, querySequenceLength );

            if ( ReadStrand == QUERY_POS_STRAND ) {printf ( "+\t" );}
            else if ( ReadStrand == QUERY_NEG_STRAND ) {printf ( "-\t" );}
            else {printf ( "?\t" );}

            //printf("?\t");
            // printf("%d%s",MaxError,AlignmentTypeString[0]);
            numberOfAlignments++;
            printf ( "\n" );
        }

        i++;
    }

    //printf("\n");
    // fprintf(stderr,"\nThere were %u alignments.\n",numberOfAlignments);
    fclose ( outputFile );
    //fclose(queryFile);
    free ( genomeLength );
    free ( genomeAnnotation );
    return 0;
}
