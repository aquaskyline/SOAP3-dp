/*
 *
 *    FileUtilities.c
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

#include "FileUtilities.h"


// This function is to get the number of lines inside a file
int FUGetNumOfLines ( char * fileName )
{
    FILE * filein;
    char buffer[INPUT_BUFFER_SIZE];
    filein = ( FILE * ) fopen ( fileName, "r" );
    size_t bufferSize = fread ( buffer, sizeof ( char ), INPUT_BUFFER_SIZE, filein );
    int numlines = 0;
    int i;

    while ( bufferSize > 0 )
    {
        // count how many '\n' inside the buffer
        for ( i = 0; i < bufferSize; i++ )
        {
            if ( buffer[i] == '\n' )
            { numlines++; }
        }

        bufferSize = fread ( buffer, sizeof ( char ), INPUT_BUFFER_SIZE, filein );
    }

    fclose ( filein );
    return numlines;
}


// This function is to get the next field
// The fields are sepaated by tab
void FUGetNextField ( FILE * filein, char * buffer, size_t & bufferSize, int & bufferIndex, char * field,
                      int & isEndOfLine )
{
    isEndOfLine = 0;
    int index = 0;
    char c;

    while ( bufferSize > 0 )
    {
        // not end of file
        if ( bufferIndex >= bufferSize )
        {
            size_t bufferSize = fread ( buffer, sizeof ( char ), INPUT_BUFFER_SIZE, filein );
            bufferIndex = 0;

            if ( bufferSize == 0 ) // reach end of file
            { break; }
        }

        // scan until it reaches '\t' or '\n'
        while ( bufferIndex < bufferSize )
        {
            c = buffer[bufferIndex++];

            switch ( c )
            {
                case '\n':
                    isEndOfLine = 1;

                case '\t':
                    field[index] = '\0';
                    return;
                    break;

                default:
                    if ( index < MAX_FIELD_LEN - 1 )
                    {
                        field[index++] = c;
                    }

                    break;
            }
        }
    }

    field[index] = '\0';
}

// This function is to skip the characters until the end of line
void FUSkipToEOL ( FILE * filein, char * buffer, size_t & bufferSize, int & bufferIndex )
{
    while ( bufferSize > 0 )
    {
        // not end of file
        if ( bufferIndex >= bufferSize )
        {
            size_t bufferSize = fread ( buffer, sizeof ( char ), INPUT_BUFFER_SIZE, filein );
            bufferIndex = 0;

            if ( bufferSize == 0 ) // reach end of file
            { break; }
        }

        // scan until it reaches '\n'
        while ( bufferIndex < bufferSize )
        {
            if ( buffer[bufferIndex++] == '\n' )
            {
                return;
            }
        }
    }
}

