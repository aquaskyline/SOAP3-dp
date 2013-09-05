/*
 *
 *    PE.c
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

#include "PE.h"

unsigned int patternidt[MAX_NUM_CPU_THREADS][ ( ( MAX_READ_LENGTH >> 6 ) + 1 ) << 4] __attribute__ ( ( aligned ( 16 ) ) );

void createQueryPackedDNA ( unsigned char * query, unsigned int query_length, unsigned int * outPackedDNA )
{
    unsigned int numOfWord = ( query_length + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD;
    // initialize outPackedDNA
    unsigned int i;

    for ( i = 0; i < numOfWord; i++ )
    {
        outPackedDNA[i] = 0;
    }

    for ( i = 0; i < query_length; i++ )
    {
        outPackedDNA[i / CHAR_PER_WORD] |= ( ( ( unsigned int ) query[i] ) << ( ( CHAR_PER_WORD - 1 - i % CHAR_PER_WORD ) * BIT_PER_CHAR ) );
    }
}

void createRevQueryPackedDNA ( unsigned char * query, unsigned int query_length, unsigned int * outPackedDNA )
{
    unsigned int numOfWord = ( query_length + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD;
    // initialize outPackedDNA
    unsigned int i;

    for ( i = 0; i < numOfWord; i++ )
    {
        outPackedDNA[i] = 0;
    }

    for ( i = 0; i < query_length; i++ )
    {
        outPackedDNA[i / CHAR_PER_WORD] |= ( ( soap3DnaComplement[ ( unsigned int ) query[query_length - 1 - i]] ) << ( ( CHAR_PER_WORD - 1 - i % CHAR_PER_WORD ) * BIT_PER_CHAR ) );
    }
}

int writeULLToStr ( unsigned long long num, char * str )
{
    // write unsigned long long to string
    // return size of the text of the number
    char * currpos = str;
    int digit = 1;

    if ( num > 0 )
    {
        digit = ( ( int ) log10 ( num ) ) + 1;
    }

    for ( int i = digit - 1; i >= 0; i-- )
    {
        currpos[i] = num % 10 + '0';
        num = num / 10;
    }

    return digit + ( currpos - str );
}

int writeNumToStr ( int num, char * str )
{
    // write number to string
    // return size of the text of the number
    char * currpos = str;

    if ( num < 0 )
    {
        num = -num;
        ( *currpos ) = '-';
        currpos++;
    }

    int digit = 1;

    if ( num > 0 )
    {
        digit = ( ( int ) log10 ( num ) ) + 1;
    }

    for ( int i = digit - 1; i >= 0; i-- )
    {
        currpos[i] = num % 10 + '0';
        num = num / 10;
    }

    return digit + ( currpos - str );
}

char numMismatch ( unsigned int * query, unsigned int * target, unsigned int query_length, unsigned int max_mismatch )
{
    // 2 bit per char
    unsigned int i;
    unsigned int numOfWord = ( query_length + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD;
    /*
    unsigned int j;
    fprintf(stderr, "Query:\n");
    for (i = 0; i < numOfWord * 16; ++i) {
          j = (query[i / 16] >> ((15 - i%16)*2)) & 3;
          fprintf(stderr, "%u ", j);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "Target:\n");
    for (i = 0; i < numOfWord * 16; ++i) {
          j = (target[i / 16] >> ((15 - i%16)*2)) & 3;
          fprintf(stderr, "%u ", j);
    }
    fprintf(stderr, "\n");
     */
    unsigned int pattern;
    char mismatches = 0;

    for ( i = 0; i < numOfWord && mismatches <= max_mismatch; i++ )
    {
        pattern = query[i] ^ target[i];
        pattern = ( ( ( pattern & 0xAAAAAAAA ) >> 1 ) | ( pattern & 0x55555555 ) );
        // 0xAAAAAAAA = 101010....10
        // 0x55555555 = 010101....01
        mismatches += __builtin_popcount ( pattern );
    }

    // fprintf(stderr, "Number of mismatches: %u\n", mismatches);
    return mismatches;
}

int numMismatchNew ( unsigned int * query, unsigned int * target, unsigned int query_length, int threadID )
{
    // 2 bit per char
    unsigned int i;
    unsigned int numOfWord = ( query_length + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD;
    /*
     unsigned int j;
     fprintf(stderr, "Query:\n");
     for (i = 0; i < numOfWord * 16; ++i) {
     j = (query[i / 16] >> ((15 - i%16)*2)) & 3;
     fprintf(stderr, "%u ", j);
     }
     fprintf(stderr, "\n");
     fprintf(stderr, "Target:\n");
     for (i = 0; i < numOfWord * 16; ++i) {
     j = (target[i / 16] >> ((15 - i%16)*2)) & 3;
     fprintf(stderr, "%u ", j);
     }
     fprintf(stderr, "\n");
     */
#define USE_OLD_POPCOUNT
#ifdef USE_OLD_POPCOUNT
    unsigned int pattern;
    int mismatches = 0;

    for ( i = 0; i < numOfWord; i++ )
    {
        pattern = query[i] ^ target[i];
        pattern = ( ( ( pattern & 0xAAAAAAAA ) >> 1 ) | ( pattern & 0x55555555 ) );
        // 0xAAAAAAAA = 101010....10
        // 0x55555555 = 010101....01
        mismatches += __builtin_popcount ( pattern );
    }

    return mismatches;
#else

    for ( i = 1; i < numOfWord; i += 2 )
    {
        const int j = i - 1;
        patternidt[threadID][j] = query[j] ^ target[j];
        patternidt[threadID][i] = query[i] ^ target[i];
        patternidt[threadID][j] = ( ( ( patternidt[threadID][j] & 0xAAAAAAAA ) >> 1 ) | ( patternidt[threadID][j] & 0x55555555 ) );
        patternidt[threadID][i] = ( ( ( patternidt[threadID][i] & 0xAAAAAAAA ) >> 1 ) | ( patternidt[threadID][i] & 0x55555555 ) );
    }

    if ( numOfWord & 1 )
    {
        const int j = numOfWord - 1;
        patternidt[threadID][j] = query[j] ^ target[j];
        patternidt[threadID][j] = ( ( ( patternidt[threadID][j] & 0xAAAAAAAA ) >> 1 ) | ( patternidt[threadID][j] & 0x55555555 ) );
    }

    memset ( patternidt[threadID] + numOfWord, 0, ( ( ALPHABET_SIZE - ( numOfWord & CHAR_MASK ) ) & CHAR_MASK ) << BIT_PER_CHAR );
    // calculate how many 16-bytes for this query
    int num_16_bytes = ( ( query_length + 63 ) >> 6 );
    return ssse3_popcount2 ( ( uint8_t * ) patternidt[threadID], num_16_bytes );
#endif
}


int mdStr ( unsigned int * query, char * qualities, unsigned int * target, unsigned int query_length, char mismatchNum, char * md_str, int * avg_mismatch_qual )
{
    // return the MD:Z string
    // 2 bit per char
    // the md string is a char array with length > query's length * 2
    // return the size of MD string
    int md_size = 0;

    if ( mismatchNum == 0 )
    {
        // printf("mismatchNum = %i\n", mismatchNum);
        // printf("MD:Z:%u\n", query_length);
        md_size += writeNumToStr ( query_length, md_str );
        md_str[md_size] = '\0';
        return md_size;
    }

    //char nucl[4] = {'A', 'C', 'G', 'T'};
    double sum_mismatch_qual = 0.0;
    // 2863311530 = 101010....10
    // 1431655765 = 010101....01
    // 2147483648 = 100000....00
    unsigned int i;
    unsigned int j;
    unsigned int numOfWord = ( query_length + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD;
    /*
    printf("Query:\n");
    for (i = 0; i < numOfWord * 16; ++i) {
           j = (query[i / 16] >> ((15 - i%16)*2)) & 3;
           printf("%c", nucl[j]);
    }
    printf("\n");
    printf("Target:\n");
    for (i = 0; i < numOfWord * 16; ++i) {
           j = (target[i / 16] >> ((15 - i%16)*2)) & 3;
           printf("%c", nucl[j]);
    }
    printf("\n");
    printf("mismatchNum = %i\n", mismatchNum);
    */
    // compute the MD string
    //printf("MD:Z:");
    unsigned int pattern;
    char mismatches = 0;
    int preMismatchPos = -1;
    int leadingZero;
    int mismatchPos;

    for ( i = 0; i < numOfWord; i++ )
    {
        pattern = query[i] ^ target[i];
        pattern = ( ( pattern & 2863311530 ) | ( ( pattern & 1431655765 ) << 1 ) );
        mismatches = __builtin_popcount ( pattern );

        for ( j = 0; j < mismatches; j++ )
        {
            leadingZero = __builtin_clz ( pattern ) >> 1;
            mismatchPos = leadingZero + ( i * CHAR_PER_WORD );
            md_size += writeNumToStr ( mismatchPos - preMismatchPos - 1, & ( md_str[md_size] ) );
            // md_str[md_size++] = nucl[(query[mismatchPos / 16] >> ((15 - mismatchPos%16)*2)) & 3];
            md_str[md_size++] = dnaChar[ ( target[mismatchPos / CHAR_PER_WORD] >> ( ( CHAR_PER_WORD - 1 - mismatchPos % CHAR_PER_WORD ) * BIT_PER_CHAR ) ) & CHAR_MASK];
            // printf("%i", mismatchPos-preMismatchPos-1);
            // printf("%c", nucl[(query[mismatchPos / 16] >> ((15 - mismatchPos%16)*2)) & 3]);
            // remove the corresponding "1" inside the "pattern"
            pattern ^= ( 2147483648 >> ( leadingZero << 1 ) );
            preMismatchPos = mismatchPos;
            sum_mismatch_qual += qualities[mismatchPos];
        }
    }

    md_size += writeNumToStr ( query_length - preMismatchPos - 1, & ( md_str[md_size] ) );
    md_str[md_size] = '\0';
    // printf("%i\n", query_length-preMismatchPos-1);
    ( *avg_mismatch_qual ) = ( int ) ( sum_mismatch_qual / mismatchNum );
    return md_size;
}


void createTargetPackedDNA ( unsigned int * target, unsigned int pos, unsigned int target_length,
                             unsigned int * outPackedDNA )
{
    unsigned int i;
    /*
    unsigned int j;
    printf("Target (before packed):\n");
    printf("pos = %u; pos mod 16 = %u\n", pos, pos%16);
    for (i=0; i<target_length; i++) {
          j = (target[(pos+i)/16] >> (15-(pos+i)%16)*2) & 3;
          printf("%u ", j);
    }
    printf("\n");
    */
    unsigned int numOfWord = ( target_length + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD;

    // initialize targetPackedDNA
    for ( i = 0; i < numOfWord; i++ )
    {
        outPackedDNA[i] = 0;
    }

    // get the position from l  (hsp->packedDNA)
    // unsigned int s = (*bwt->_bwtSaValue)(bwt,l);
    for ( i = 0; i < numOfWord; i++ )
    {
        outPackedDNA[i] |= ( target[ ( pos + i * CHAR_PER_WORD ) / CHAR_PER_WORD] << ( ( pos + i * CHAR_PER_WORD ) % CHAR_PER_WORD * BIT_PER_CHAR ) );

        if ( ( pos % CHAR_PER_WORD ) > 0 )
        { outPackedDNA[i] |= ( target[ ( pos + i * CHAR_PER_WORD ) / CHAR_PER_WORD + 1] >> ( ( CHAR_PER_WORD - ( pos + i * CHAR_PER_WORD ) % CHAR_PER_WORD ) * BIT_PER_CHAR ) ); }
    }

    // the length of the last word
    unsigned int lastWordLen = ( target_length - 1 ) % CHAR_PER_WORD + 1; // range ~ [1,16]

    if ( lastWordLen < CHAR_PER_WORD )
        // 4294967295 = 11111...1
    { outPackedDNA[numOfWord - 1] &= ( 4294967295 << ( ( CHAR_PER_WORD - lastWordLen ) * BIT_PER_CHAR ) ); }
}

void SRAEnrichSARanges ( HSP * hsp, BWT * bwt, unsigned int * query, unsigned int * rev_query, unsigned int query_length,
                         unsigned int l, unsigned int r, char max_mismatch, char * num_mis, char * strand )
{
    // SRAEnrichSARanges takes a SA range as input and return the mismatch count and the strand of the SA range.
    // create the packedDNA for the target sequence;
    unsigned int targetPackedDNA[ ( MAX_READ_LENGTH + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD];
    createTargetPackedDNA ( hsp->packedDNA, ( *bwt->_bwtSaValue ) ( bwt, l ), query_length, targetPackedDNA );
    ( *num_mis ) = numMismatch ( query, targetPackedDNA, query_length, max_mismatch );
    ( *strand ) = QUERY_POS_STRAND;

    if ( ( *num_mis ) > max_mismatch )
    {
        ( *num_mis ) = numMismatch ( rev_query, targetPackedDNA, query_length, max_mismatch );
        ( *strand ) = QUERY_NEG_STRAND;
    }
}

/*
char SRAEnrichSARangeswithStrand(HSP* hsp, BWT* bwt, unsigned int* query, unsigned int* rev_query, unsigned int query_length,
                       unsigned int l, unsigned int r, char strand) {
      // It takes a SA range [l,r] as well as the strand as input and return the number of mismatches of the SA range.
      // strand = 1: forward; strand = 2: reverse

      // create the packedDNA for the target sequence;
      unsigned int targetPackedDNA[(MAX_READ_LENGTH+15)/16];
      createTargetPackedDNA(hsp->packedDNA, (*bwt->_bwtSaValue)(bwt,l), query_length, targetPackedDNA);
      if (strand==1)
            return numMismatchAll(query, targetPackedDNA, query_length);
      else
            return numMismatchAll(rev_query, targetPackedDNA, query_length);
}
*/

/*
char MismatchInPos(HSP* hsp, BWT* bwt, unsigned int* query, unsigned int* rev_query, unsigned int query_length, unsigned int pos, char strand) {
      // MismatchInPos takes a position and the strand as input and return the mismatch count

      // create the packedDNA for the target sequence;
      unsigned int targetPackedDNA[(MAX_READ_LENGTH+15)/16];
      createTargetPackedDNA(hsp->packedDNA, pos, query_length, targetPackedDNA);
      if (strand==1)
            return numMismatchAll(query, targetPackedDNA, query_length);
      else
            return numMismatchAll(rev_query, targetPackedDNA, query_length);
}
 */

int getMdStr ( HSP * hsp, unsigned char * query, char * qualities, unsigned int query_length, unsigned int pos, char strand, char mismatchNum, char * md_str, int * avg_mismatch_qual, int trim )
{

    if ( strand == QUERY_POS_STRAND )
    {
        if ( trim > 0 )
        {
            query += trim;
            pos += trim;
            query_length -= trim;
        }
        else if ( trim < 0 )
        {
            query_length -= -trim;
        }
    }
    else
    {
        if ( trim > 0 )
        {
            pos += trim;
            query_length -= trim;
        }
        else if ( trim < 0 )
        {
            query += -trim;
            query_length -= -trim;
        }
    }

    // compute the MD string
    // return the size of MD string
    ( *avg_mismatch_qual ) = DEFAULT_QUAL_VALUE;
    // create query packed DNA
    unsigned int outPackedDNA[ ( MAX_READ_LENGTH + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD];

    if ( strand == QUERY_NEG_STRAND )
    { createRevQueryPackedDNA ( query, query_length, outPackedDNA ); }
    else
    { createQueryPackedDNA ( query, query_length, outPackedDNA ); }

    // create the packedDNA for the target sequence;
    unsigned int targetPackedDNA[ ( MAX_READ_LENGTH + CHAR_PER_WORD - 1 ) / CHAR_PER_WORD];
    createTargetPackedDNA ( hsp->packedDNA, pos, query_length, targetPackedDNA );
    return mdStr ( outPackedDNA, qualities, targetPackedDNA, query_length, mismatchNum, md_str, avg_mismatch_qual );
}

int convertToCigarStr ( char * special_cigar, char * cigar, int * deletedEnd )
{
    // to convert the special_cigar into cigar string
    int currInt = 0;
    int currM = 0;
    int cigar_len = 0;

    if ( deletedEnd ) { *deletedEnd = 0; }

    for ( int i = 0; i < strlen ( special_cigar ); i++ )
    {
        if ( isdigit ( special_cigar[i] ) )
        {
            currInt = currInt * 10 + ( special_cigar[i] - '0' );
        }
        else
        {
            switch ( special_cigar[i] )
            {
                case 'M':
                case 'm':
                    currM += currInt;
                    currInt = 0;
                    break;

                case 'D':
                    if ( ( cigar_len == 0 && currM == 0 ) || ( i == strlen ( special_cigar ) - 1 ) )
                    {
                        // the deletion in front or at the end is ignored
                        if ( i == strlen ( special_cigar ) - 1 )
                        {
                            if ( deletedEnd ) { *deletedEnd = currM; }
                        }

                        break;
                    }

                case 'I':
                case 'S':
                    if ( currM > 0 )
                    {
                        cigar_len += writeNumToStr ( currM, & ( cigar[cigar_len] ) );
                        cigar[cigar_len++] = 'M';
                        currM = 0;
                    }

                    cigar_len += writeNumToStr ( currInt, & ( cigar[cigar_len] ) );
                    cigar[cigar_len++] = special_cigar[i];
                    currInt = 0;
                    break;
            }
        }
    }

    if ( currM > 0 )
    {
        cigar_len += writeNumToStr ( currM, & ( cigar[cigar_len] ) );
        cigar[cigar_len++] = 'M';
        currM = 0;
    }

    cigar[cigar_len] = '\0';
    // printf("%s\n", cigar);
    return cigar_len;
}


int convertToCigarStr2 ( char * special_cigar )
{
    // to convert the special_cigar into cigar string
    // and save back to the cigar variable
    char cigar[MAX_READ_LENGTH + 1];
    int cigar_len = convertToCigarStr ( special_cigar, cigar );
    strcpy ( special_cigar, cigar );
    return cigar_len;
}


int getMisInfoForDP ( HSP * hsp, unsigned char * query, char * qualities, unsigned int query_length,
                      unsigned int pos, char strand, char * special_cigar, char * md_str, int * numMismatch,
                      int * gapOpen, int * gapExt, int * avg_mismatch_qual, int trim )
{

    if ( strand == QUERY_POS_STRAND )
    {
        if ( trim > 0 )
        {
            query += trim;
            pos += trim;
            query_length -= trim;
        }
        else if ( trim < 0 )
        {
            query_length -= -trim;
        }
    }
    else
    {
        if ( trim > 0 )
        {
            pos += trim;
            query_length -= trim;
        }
        else if ( trim < 0 )
        {
            query += -trim;
            query_length -= -trim;
        }
    }

    // to compute the MD string, the number of mismatches, the number of gap open and the number of gap extension
    // return the size of MD string
    // this function is designed for DP module
    // the md string is a char array with length > query's length * 2
    // special_cigar includes the "m" character for mismatches
    ( *avg_mismatch_qual ) = DEFAULT_QUAL_VALUE;
    // create corresponding query DNA
    unsigned char revQuery[MAX_READ_LENGTH];
    unsigned char * currQuery;
    unsigned int i, j, k;
    //char nucl[4] = {'A', 'C', 'G', 'T'};
    double sum_mismatch_qual = 0.0;

    if ( strand == QUERY_NEG_STRAND )
    {
        for ( i = 0; i < query_length; i++ )
        { revQuery[i] = soap3DnaComplement[query[query_length - i - 1]]; }

        currQuery = revQuery;
    }
    else
    {
        currQuery = query;
    }

    /*
    // print the query sequence
    printf("Query:\n");
    for (i=0; i<query_length; i++)
    printf("%c", nucl[currQuery[i]]);
    printf("\n");

    printf("Target:\n");
    for (i=0; i<query_length; i++) {
    j = (hsp->packedDNA[(pos+i)/16] >> (15-(pos+i)%16)*2) & 3;
    printf("%c", nucl[j]);
    }
    printf("\n");

    printf("special CIGAR: %s\n", special_cigar);
    */
    *numMismatch = 0;
    *gapOpen = 0;
    *gapExt = 0;
    int currInt = 0;
    int currMatch = 0;
    int qPos = 0;
    unsigned int tPos = pos;
    // compute the MD string
    int md_len = 0;

    // printf("MD:Z: ");
    int l = strlen ( special_cigar );

    for ( i = 0; i < l; i++ )
    {
        if ( isdigit ( special_cigar[i] ) )
        {
            currInt = currInt * 10 + ( special_cigar[i] - '0' );
        }
        else
        {
            switch ( special_cigar[i] )
            {
                case 'M':
                    currMatch += currInt;
                    qPos += currInt;
                    tPos += currInt;
                    currInt = 0;
                    break;

                case 'm':
                    md_len += writeNumToStr ( currMatch, & ( md_str[md_len] ) );
                    md_str[md_len++] = dnaChar[ ( hsp->packedDNA[ ( tPos ) / CHAR_PER_WORD] >> ( CHAR_PER_WORD - 1 - ( tPos ) % CHAR_PER_WORD ) * BIT_PER_CHAR ) & CHAR_MASK];
                    // md_str[md_len++] = nucl[currQuery[qPos++]];
                    sum_mismatch_qual += qualities[qPos];

                    for ( j = 1; j < currInt; j++ )
                    {
                        md_str[md_len++] = '0';
                        md_str[md_len++] = dnaChar[ ( hsp->packedDNA[ ( tPos + j ) / CHAR_PER_WORD] >> ( CHAR_PER_WORD - 1 - ( tPos + j ) % CHAR_PER_WORD ) * BIT_PER_CHAR ) & CHAR_MASK];
                        // md_str[md_len++] = nucl[currQuery[qPos++]];
                        sum_mismatch_qual += qualities[qPos + j];
                    }

                    qPos += currInt;
                    tPos += currInt;
                    ( *numMismatch ) += currInt;
                    currMatch = 0;
                    currInt = 0;
                    break;

                case 'I':
                    qPos += currInt;
                    ( *gapOpen ) ++;
                    ( *gapExt ) += currInt;
                    currInt = 0;
                    break;

                case 'D':
                    if ( i == l - 1 ) { break; } // last delete, ignore

                    md_len += writeNumToStr ( currMatch, & ( md_str[md_len] ) );
                    md_str[md_len++] = '^';

                    for ( j = 0; j < currInt; j++ )
                    {
                        k = ( hsp->packedDNA[ ( tPos + j ) / CHAR_PER_WORD] >> ( CHAR_PER_WORD - 1 - ( tPos + j ) % CHAR_PER_WORD ) * BIT_PER_CHAR ) & CHAR_MASK;
                        md_str[md_len++] = dnaChar[k];
                    }

                    tPos += currInt;
                    ( *gapOpen ) ++;
                    ( *gapExt ) += currInt;
                    currMatch = 0;
                    currInt = 0;
                    break;

                case 'S':
                    qPos += currInt;
                    currInt = 0;
                    break;
            }
        }
    }

    md_len += writeNumToStr ( currMatch, & ( md_str[md_len] ) );
    md_str[md_len] = '\0';

    // printf("%s\n", md_str);
    //   printf("numMismatch: %i\n", (*numMismatch));
    if ( ( *numMismatch ) > 0 )
    { ( *avg_mismatch_qual ) = ( int ) ( sum_mismatch_qual / ( *numMismatch ) ); }

    return md_len;
}
