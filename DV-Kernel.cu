/*
 *
 *    DV-Kernel.cu
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

#include "DV-Kernel.h"

// copied from DNACount.c
__forceinline__ __device__ uint GPUDNAOccCount ( uint * dna, uint index, char character, uint backward )
{
    uint wordToCount, charToCount;
    uint i;
    uint sum = 0;
    wordToCount = index / 32;
    charToCount = index - wordToCount * 32;
    dna -= backward * 4;
    ulonglong2 dd;
    //ulonglong2 is a 16Byte = 128bit CUDA vector
    //.x is the first 8Byte while .y is the last.
    dd = * ( ( ulonglong2 * ) dna );
    i = 0;

    if ( wordToCount > 0 ) // i = 0
    {
        unsigned long long d = backward ? dd.y : dd.x;
        unsigned long long b = ( d >> 1 );
        b ^= ( ( character & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
        unsigned long long a = d ^ ( ( character & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
        a &= b;
        a &= 0x5555555555555555;
        sum += __popcll ( a );
        i = 1;
    }

    if ( wordToCount == 2 || charToCount > 0 )
    {
        if ( backward ) { i = 3 - i; }

        if ( wordToCount == 2 ) { charToCount = 32; }

        unsigned long long d = ( i & 1 ) ? dd.y : dd.x;
        unsigned long long b = ( d >> 1 );
        b ^= ( ( character & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
        unsigned long long a = d ^ ( ( character & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
        a &= b;
        unsigned long long mask;

        if ( charToCount < 16 )
        {
            mask = ( 0x00000000FFFFFFFF << ( 32 - charToCount * 2 ) ) & 0x00000000FFFFFFFF;
        }
        else
        {
            mask = ( 0xFFFFFFFF00000000 << ( 64 - charToCount * 2 ) ) | 0x00000000FFFFFFFF;
        }

        if ( backward )
        { mask = __brevll ( mask ); }

        a &= ( mask & 0x5555555555555555 );
        sum += __popcll ( a );
    }

    return sum;
}

__forceinline__ __device__ void GPUDNAAllOccCount ( uint * dna, uint index, uint backward, uint occCount[] )
{
    uint wordToCount, charToCount;
    uint i;
    wordToCount = index / 32;
    charToCount = index - wordToCount * 32;
    dna -= backward * 4;
    ulonglong2 dd;
    dd = * ( ( ulonglong2 * ) dna );
    i = 0;

    if ( wordToCount > 0 ) // i = 0
    {
        unsigned long long d = backward ? dd.y : dd.x;

        for ( int j = 0; j < ALPHABET_SIZE; ++j )
        {
            unsigned long long b = ( d >> 1 );
            b ^= ( ( j & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            unsigned long long a = d ^ ( ( j & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            a &= b;
            a &= 0x5555555555555555;

            if ( backward )
            { occCount[j] -= __popcll ( a ); }
            else
            { occCount[j] += __popcll ( a ); }
        }

        i = 1;
    }

    if ( wordToCount == 2 || charToCount > 0 )
    {
        if ( backward ) { i = 3 - i; }

        if ( wordToCount == 2 ) { charToCount = 32; }

        unsigned long long d = ( i & 1 ) ? dd.y : dd.x;

        for ( int j = 0; j < ALPHABET_SIZE; ++j )
        {
            unsigned long long b = ( d >> 1 );
            b ^= ( ( j & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            unsigned long long a = d ^ ( ( j & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            a &= b;
            unsigned long long mask;

            if ( charToCount < 16 )
            {
                mask = ( 0x00000000FFFFFFFF << ( 32 - charToCount * 2 ) ) & 0x00000000FFFFFFFF;
            }
            else
            {
                mask = ( 0xFFFFFFFF00000000 << ( 64 - charToCount * 2 ) ) | 0x00000000FFFFFFFF;
            }

            if ( backward )
            { mask = __brevll ( mask ); }

            a &= ( mask & 0x5555555555555555 );

            if ( backward )
            { occCount[j] -= __popcll ( a ); }
            else
            { occCount[j] += __popcll ( a ); }
        }
    }
}

__forceinline__ __device__ uint GPUDNAOccCountWithCumu ( uint * dna, uint index, char c, uint backward, uint & cumu )
{
    uint wordToCount, charToCount;
    uint i;
    uint sum = 0;
    wordToCount = index / 32;
    charToCount = index - wordToCount * 32;
    dna -= backward * 4;
    ulonglong2 dd;
    dd = * ( ( ulonglong2 * ) dna );
    i = 0;

    if ( wordToCount > 0 ) // i = 0
    {
        unsigned long long d = backward ? dd.y : dd.x;

        for ( int j = 3; j > c; --j ) // TODO hardcode
        {
            unsigned long long b = ( d >> 1 );
            b ^= ( ( j & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            unsigned long long a = d ^ ( ( j & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            a &= b;
            a &= 0x5555555555555555;
            cumu += __popcll ( a );
        }

        {
            // copy&paste for c
            unsigned long long b = ( d >> 1 );
            b ^= ( ( c & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            unsigned long long a = d ^ ( ( c & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            a &= b;
            a &= 0x5555555555555555;
            sum += __popcll ( a );
        }

        i = 1;
    }

    if ( wordToCount == 2 || charToCount > 0 )
    {
        if ( backward ) { i = 3 - i; }

        if ( wordToCount == 2 ) { charToCount = 32; }

        unsigned long long d = ( i & 1 ) ? dd.y : dd.x;

        for ( int j = 3; j > c; --j ) // hardcode
        {
            unsigned long long b = ( d >> 1 );
            b ^= ( ( j & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            unsigned long long a = d ^ ( ( j & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            a &= b;
            unsigned long long mask;

            if ( charToCount < 16 )
            {
                mask = ( 0x00000000FFFFFFFF << ( 32 - charToCount * 2 ) ) & 0x00000000FFFFFFFF;
            }
            else
            {
                mask = ( 0xFFFFFFFF00000000 << ( 64 - charToCount * 2 ) ) | 0x00000000FFFFFFFF;
            }

            if ( backward )
            { mask = __brevll ( mask ); }

            a &= ( mask & 0x5555555555555555 );
            cumu += __popcll ( a );
        }

        {
            // copy&paste for c
            unsigned long long b = ( d >> 1 );
            b ^= ( ( c & 0x2 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            unsigned long long a = d ^ ( ( c & 0x1 ) ? 0 : 0xFFFFFFFFFFFFFFFF );
            a &= b;
            unsigned long long mask;

            if ( charToCount < 16 )
            {
                mask = ( 0x00000000FFFFFFFF << ( 32 - charToCount * 2 ) ) & 0x00000000FFFFFFFF;
            }
            else
            {
                mask = ( 0xFFFFFFFF00000000 << ( 64 - charToCount * 2 ) ) | 0x00000000FFFFFFFF;
            }

            if ( backward )
            { mask = __brevll ( mask ); }

            a &= ( mask & 0x5555555555555555 );
            sum += __popcll ( a );
        }
    }

    return sum;
}


// copied from BWT.c
__forceinline__ __device__ uint GPUBWTOccValue ( uint * bwt, uint * occ, uint index, char c, uint inverseSa0 )
{
    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    index -= ( index > inverseSa0 );
    uint occExplicitIndex = ( index + GPU_OCC_INTERVAL / 2 - 1 ) / GPU_OCC_INTERVAL; // Bidirectional encoding
    uint occIndex = occExplicitIndex * GPU_OCC_INTERVAL;
    uint occValue = occ[occExplicitIndex * ALPHABET_SIZE + c];

    if ( occIndex != index )
    {
        // __usad(x,y,z) = |x-y| + z
        // GPUDNAOccCount(explicitIndex, len, character, isBackwardDirectionCount)
        //   - pointer to position of explicitIndex on BWT: Reference point on sampled occ
        //   - len: number of characters on BWT to count
        //   - character: the character to count
        //   - isBackwardDirectionCount: a boolean flag to identify backward counting
        uint cnt = GPUDNAOccCount ( bwt + occIndex / CHAR_PER_WORD, __usad ( index, occIndex, 0 ), c, occIndex > index );
        return occIndex > index ? occValue - cnt : occValue + cnt;
    }
    else
    {
        return occValue;
    }
}

__forceinline__ __device__ void GPUBWTAllOccValue ( uint * bwt,
        uint * occ,
        uint inverseSa0,
        uint index,
        uint occCount[] )
{
    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    index -= ( index > inverseSa0 );
    uint occExplicitIndex = ( index + GPU_OCC_INTERVAL / 2 - 1 ) / GPU_OCC_INTERVAL; // Bidirectional encoding
    uint occIndex = occExplicitIndex * GPU_OCC_INTERVAL;
    * ( ( uint4 * ) occCount ) = * ( ( uint4 * ) ( occ + occExplicitIndex * ALPHABET_SIZE ) );

    if ( occIndex != index )
    {
        GPUDNAAllOccCount ( bwt + occIndex / CHAR_PER_WORD, __usad ( index, occIndex, 0 ), occIndex > index, occCount );
    }
}

__forceinline__ __device__ uint GPUBWTOccValueWithCumu ( uint * bwt, uint * occ, uint index, char c, uint inverseSa0, uint & cumu )
{
    // $ is supposed to be positioned at inverseSa0 but it is not encoded
    // therefore index is subtracted by 1 for adjustment
    index -= ( index > inverseSa0 );
    uint occExplicitIndex = ( index + GPU_OCC_INTERVAL / 2 - 1 ) / GPU_OCC_INTERVAL; // Bidirectional encoding
    uint occIndex = occExplicitIndex * GPU_OCC_INTERVAL;
    uint occValues[ALPHABET_SIZE];
    * ( ( uint4 * ) occValues ) = * ( ( uint4 * ) ( occ + occExplicitIndex * ALPHABET_SIZE ) );
    uint occValue = occValues[c];
    cumu = 0;

    for ( int j = 3; j > c; --j ) // TODO hardcode
    {
        cumu += occValues[j];
    }

    uint cum = 0;

    if ( occIndex != index )
    {
        uint cnt = GPUDNAOccCountWithCumu ( bwt + occIndex / CHAR_PER_WORD, __usad ( index, occIndex, 0 ), c, occIndex > index, cum );

        if ( occIndex > index )
        {
            cumu -= cum;
            return occValue - cnt;
        }

        cumu += cum;
        return occValue + cnt;
    }

    return occValue;
}

__forceinline__ __device__
void contBackwardSearch ( uint * query, uint start, uint len,
                          uint * bwt, uint * occ, uint inverseSa0,
                          uint saL, uint saR,
                          uint & saCount, uint * output,
                          uint maxSARangesAllowed,
                          uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL = GPUBWTOccValue ( bwt, occ, saL, c, inverseSa0 ) + 1;
        saR = GPUBWTOccValue ( bwt, occ, saR + 1, c, inverseSa0 );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[contBackwardSearch] Reporting of SA ranges %u %u ...", saL, saR );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = saL;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( saR - saL ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}

__forceinline__ __device__
void revContForwardSearch ( uint * query, uint start, uint len,
                            uint * revBwt, uint * revOcc, uint revInverseSa0,
                            uint saL, uint saR,
                            uint revSaL, uint revSaR,
                            uint & saCount, uint * output,
                            uint maxSARangesAllowed,
                            uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint start, end;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, end );
        saR = saR + start - end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[revContForwardSearch] Reporting of SA ranges %u %u ...", saL, saR );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = saL;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( saR - saL ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}


__forceinline__ __device__
void backward1Mismatch ( uint * query, uint start, uint len,
                         uint * bwt, uint * occ, uint inverseSa0,
                         uint pl, uint pr,
                         uint & saCount, uint * output,
                         uint maxSARangesAllowed,
                         uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr  && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];
            contBackwardSearch ( query, start, i - start - 1,
                                 bwt, occ, inverseSa0,
                                 mkL, mkR,
                                 saCount, output,
                                 maxSARangesAllowed,
                                 strand, accumMismatches + 1 );

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }
}

__forceinline__ __device__
void backward1MismatchAndExact ( uint * query, uint start, uint len,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint pl, uint pr,
                                 uint & saCount, uint * output,
                                 uint maxSARangesAllowed,
                                 uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];

            if ( mkL <= mkR )
            {
                contBackwardSearch ( query, start, i - start - 1,
                                     bwt, occ, inverseSa0,
                                     mkL, mkR,
                                     saCount, output,
                                     maxSARangesAllowed, strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }

    if ( pl <= pr && i == start && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[backward1MismatchOrExact] Reporting of SA ranges %u %u ...", pl, pr );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = pl;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( pr - pl ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}

__forceinline__ __device__
void backward2Mismatch ( uint * query, uint start, uint len,
                         uint * bwt, uint * occ, uint inverseSa0,
                         uint pl, uint pr,
                         uint & saCount, uint * output,
                         uint maxSARangesAllowed,
                         uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];
            backward1Mismatch ( query, start, i - start - 1,
                                bwt, occ, inverseSa0,
                                mkL, mkR,
                                saCount, output,
                                maxSARangesAllowed,
                                strand, accumMismatches + 1 );

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }
}

__forceinline__ __device__
void backward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
                                     uint * bwt, uint * occ, uint inverseSa0,
                                     uint pl, uint pr,
                                     uint & saCount, uint * output,
                                     uint maxSARangesAllowed,
                                     uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];

            if ( mkL <= mkR )
            {
                backward1MismatchAndExact ( query, start, i - start - 1,
                                            bwt, occ, inverseSa0,
                                            mkL, mkR,
                                            saCount, output,
                                            maxSARangesAllowed,
                                            strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }
}

__forceinline__ __device__
void backward2MismatchAnd1MismatchAndExact ( uint * query, uint start, uint len,
        uint * bwt, uint * occ, uint inverseSa0,
        uint pl, uint pr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed,
        uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];

            if ( mkL <= mkR )
            {
                backward1MismatchAndExact ( query, start, i - start - 1,
                                            bwt, occ, inverseSa0,
                                            mkL, mkR,
                                            saCount, output,
                                            maxSARangesAllowed,
                                            strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }

    if ( pl <= pr && i == start && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[backward1MismatchOrExact] Reporting of SA ranges %u %u ...", pl, pr );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = pl;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( pr - pl ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}

__forceinline__ __device__
void backward3Mismatch ( uint * query, uint start, uint len,
                         uint * bwt, uint * occ, uint inverseSa0,
                         uint pl, uint pr,
                         uint & saCount, uint * output,
                         uint maxSARangesAllowed,
                         uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];
            backward2Mismatch ( query, start, i - start - 1,
                                bwt, occ, inverseSa0,
                                mkL, mkR,
                                saCount, output,
                                maxSARangesAllowed,
                                strand, accumMismatches + 1 );

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }
}

__forceinline__ __device__
void backward3MismatchAnd2Mismatch ( uint * query, uint start, uint len,
                                     uint * bwt, uint * occ, uint inverseSa0,
                                     uint pl, uint pr,
                                     uint & saCount, uint * output,
                                     uint maxSARangesAllowed,
                                     uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];

            if ( mkL <= mkR )
            {
                backward2MismatchAnd1Mismatch ( query, start, i - start - 1,
                                                bwt, occ, inverseSa0,
                                                mkL, mkR,
                                                saCount, output,
                                                maxSARangesAllowed,
                                                strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }
}

__forceinline__ __device__
void backward3MismatchAnd2MismatchAnd1MismatchAndExact ( uint * query, uint start, uint len,
        uint * bwt, uint * occ, uint inverseSa0,
        uint pl, uint pr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed,
        uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];

            if ( mkL <= mkR )
            {
                backward2MismatchAnd1MismatchAndExact ( query, start, i - start - 1,
                                                        bwt, occ, inverseSa0,
                                                        mkL, mkR,
                                                        saCount, output,
                                                        maxSARangesAllowed,
                                                        strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }

    if ( pl <= pr && i == start && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[backward3MismatchAnd2MismatchAnd1MismatchAndExact] Reporting of SA ranges %u %u ...", pl, pr );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = pl;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( pr - pl ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}


__forceinline__ __device__
void backward4MismatchAnd3MismatchAnd2MismatchAnd1MismatchAndExact ( uint * query, uint start, uint len,
        uint * bwt, uint * occ, uint inverseSa0,
        uint pl, uint pr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed,
        uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];
            backward3MismatchAnd2MismatchAnd1MismatchAndExact ( query, start, i - start - 1,
                    bwt, occ, inverseSa0,
                    mkL, mkR,
                    saCount, output,
                    maxSARangesAllowed,
                    strand, accumMismatches + 1 );

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }

    if ( pl <= pr && i == start && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[backward4MismatchAnd3MismatchAnd2MismatchAnd1MismatchAndExact] Reporting of SA ranges %u %u ...", pl, pr );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = pl;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( pr - pl ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}


__forceinline__ __device__
void revForward1Mismatch ( uint * query, uint start, uint len,
                           uint * revBwt, uint * revOcc, uint revInverseSa0,
                           uint pl, uint pr,
                           uint revPl, uint revPr,
                           uint & saCount, uint * output,
                           uint maxSARangesAllowed,
                           uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearch ( query, i + 1, start + len - i - 1,
                                       // search range = query[i+1, start+len)
                                       revBwt, revOcc, revInverseSa0,
                                       mkL, mkR,
                                       revMkL, revMkR,
                                       saCount, output,
                                       maxSARangesAllowed,
                                       strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void revForward1MismatchAndExact ( uint * query, uint start, uint len,
                                   uint * revBwt, uint * revOcc, uint revInverseSa0,
                                   uint pl, uint pr,
                                   uint revPl, uint revPr,
                                   uint & saCount, uint * output,
                                   uint maxSARangesAllowed,
                                   uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearch ( query, i + 1, start + len - i - 1,
                                       // search range = query[i+1, start+len)
                                       revBwt, revOcc, revInverseSa0,
                                       mkL, mkR,
                                       revMkL, revMkR,
                                       saCount, output,
                                       maxSARangesAllowed,
                                       strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }

    if ( pl <= pr && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[backward3MismatchAndExact] Reporting of SA ranges %u %u ...", pl, pr );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = pl;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( pr - pl ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}

__forceinline__ __device__
void revForward2Mismatch ( uint * query, uint start, uint len,
                           uint * revBwt, uint * revOcc, uint revInverseSa0,
                           uint pl, uint pr,
                           uint revPl, uint revPr,
                           uint & saCount, uint * output,
                           uint maxSARangesAllowed,
                           uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revForward1Mismatch ( query, i + 1, start + len - i - 1,
                                      // search range = query[i+1, start+len)
                                      revBwt, revOcc, revInverseSa0,
                                      mkL, mkR,
                                      revMkL, revMkR,
                                      saCount, output,
                                      maxSARangesAllowed,
                                      strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void revForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
                                       uint * revBwt, uint * revOcc, uint revInverseSa0,
                                       uint pl, uint pr,
                                       uint revPl, uint revPr,
                                       uint & saCount, uint * output,
                                       uint maxSARangesAllowed,
                                       uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revForward1MismatchAndExact ( query, i + 1, start + len - i - 1,
                                              // search range = query[i+1, start+len)
                                              revBwt, revOcc, revInverseSa0,
                                              mkL, mkR,
                                              revMkL, revMkR,
                                              saCount, output,
                                              maxSARangesAllowed,
                                              strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void revForward2MismatchAnd1MismatchAndExact ( uint * query, uint start, uint len,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed,
        uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revForward1MismatchAndExact ( query, i + 1, start + len - i - 1,
                                              // search range = query[i+1, start+len)
                                              revBwt, revOcc, revInverseSa0,
                                              mkL, mkR,
                                              revMkL, revMkR,
                                              saCount, output,
                                              maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }

    if ( pl <= pr && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[backward3MismatchAndExact] Reporting of SA ranges %u %u ...", pl, pr );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = pl;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( pr - pl ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}

__forceinline__ __device__
void revForward3MismatchAnd2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revForward2MismatchAnd1MismatchAndExact ( query, i + 1, start + len - i - 1,
                        // search range = query[i+1, start+len)
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void revForward3MismatchAnd2MismatchAnd1MismatchAndExact ( uint * query, uint start, uint len,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revForward2MismatchAnd1MismatchAndExact ( query, i + 1, start + len - i - 1,
                        // search range = query[i+1, start+len)
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed,
                        strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }

    if ( pl <= pr && saCount <= maxSARangesAllowed )
    {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", output );
        printf ( "[revForward3MismatchAnd2MismatchAnd1MismatchAndExact] Reporting of SA ranges %u %u ...", pl, pr );
#endif

        if ( saCount < maxSARangesAllowed )
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "DONE\n" );
#endif
            * ( output + 32 * 2 * saCount ) = pl;
            * ( output + 32 * ( 2 * saCount + 1 ) ) = ( pr - pl ) + ( strand << ( BGS_GPU_ANSWER_OFFSET_LENGTH + 3 ) ) +
                    ( accumMismatches << BGS_GPU_ANSWER_OFFSET_LENGTH );
        }
        else
        {
#ifdef BGS_OUTPUT_GENERAL_DEBUG_MESSAGE
            printf ( "FAILED\n" );
#endif
        }

        ++saCount;
    }
}

__forceinline__ __device__
void revForward4MismatchAnd3MismatchAnd2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revForward3MismatchAnd2MismatchAnd1MismatchAndExact ( query, i + 1, start + len - i - 1,
                        // search range = query[i+1, start+len)
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed,
                        strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}


__forceinline__ __device__
void revContForwardSearchAndBackward1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        backward1Mismatch ( query, start2, len2,
                            bwt, occ, inverseSa0,
                            saL, saR,
                            saCount, output,
                            maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void revContForwardSearchAndForward1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward1Mismatch ( query, start2, len2,
                              revBwt, revOcc, revInverseSa0,
                              saL, saR,
                              revSaL, revSaR,
                              saCount, output,
                              maxSARangesAllowed,
                              strand, accumMismatches );
    }
}


__forceinline__ __device__
void revContForwardSearchAndBackward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        backward2MismatchAnd1Mismatch ( query, start2, len2,
                                        bwt, occ, inverseSa0,
                                        saL, saR,
                                        saCount, output,
                                        maxSARangesAllowed,
                                        strand, accumMismatches );
    }
}


__forceinline__ __device__
void revContForwardSearchAndBackward3Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        backward3Mismatch ( query, start2, len2,
                            bwt, occ, inverseSa0,
                            saL, saR,
                            saCount, output,
                            maxSARangesAllowed,
                            strand, accumMismatches );
    }
}


__forceinline__ __device__
void revForward1MismatchAndBackward1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndBackward1Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for backward
                        revBwt, revOcc, revInverseSa0,
                        bwt, occ, inverseSa0, // for backward
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed,
                        strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void revForward1MismatchAndForward1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndForward1Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2,
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}


__forceinline__ __device__
void revContForwardSearchAndBackward3MismatchAnd2Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        backward3MismatchAnd2Mismatch ( query, start2, len2,
                                        bwt, occ, inverseSa0,
                                        saL, saR,
                                        saCount, output,
                                        maxSARangesAllowed,
                                        strand, accumMismatches );
    }
}

__forceinline__ __device__
void revForward1MismatchAndBackward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndBackward2MismatchAnd1Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for backward
                        revBwt, revOcc, revInverseSa0,
                        bwt, occ, inverseSa0, // for backward
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed,
                        strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void contBackwardSearchAndBackward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint saL, uint saR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL = GPUBWTOccValue ( bwt, occ, saL, c, inverseSa0 ) + 1;
        saR = GPUBWTOccValue ( bwt, occ, saR + 1, c, inverseSa0 );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        backward2MismatchAnd1Mismatch ( query, start2, len2,
                                        bwt, occ, inverseSa0,
                                        saL, saR,
                                        saCount, output,
                                        maxSARangesAllowed,
                                        strand, accumMismatches );
    }
}


__forceinline__ __device__
void backward1MismatchAndBackward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint pl, uint pr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];

            if ( mkL <= mkR )
            {
                contBackwardSearchAndBackward2MismatchAnd1Mismatch ( query, start, i - start - 1,
                        start2, len2,
                        bwt, occ, inverseSa0,
                        mkL, mkR,
                        saCount, output,
                        maxSARangesAllowed,
                        strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }
}


__forceinline__ __device__
void contBackwardSearchAndForward2Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL =
            GPUBWTOccValueWithCumu ( bwt, occ, saL, c, inverseSa0, cum_start ) + 1;
        saR =
            GPUBWTOccValueWithCumu ( bwt, occ, saR + 1, c, inverseSa0, cum_end );
        revSaR = revSaR + cum_start - cum_end;
        revSaL = revSaR - ( saR - saL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward2Mismatch ( query, start2, len2,
                              revBwt, revOcc, revInverseSa0,
                              saL, saR,
                              revSaL, revSaR,
                              saCount, output,
                              maxSARangesAllowed,
                              strand, accumMismatches );
    }
}

__forceinline__ __device__
void backward1MismatchAndForward2Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            mkR = occCountPEnd[ec];
            revMkR = revPr - occCountP[ec];
            revMkL = revMkR - ( mkR - mkL );

            if ( mkL <= mkR )
            {
                contBackwardSearchAndForward2Mismatch ( query, start, i - start - 1,
                                                        start2, len2,
                                                        bwt, occ, inverseSa0,
                                                        revBwt, revOcc, revInverseSa0,
                                                        mkL, mkR,
                                                        revMkL, revMkR,
                                                        saCount, output,
                                                        maxSARangesAllowed,
                                                        strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
        revPr = revPr - occCountP[c];
    }
}

__forceinline__ __device__
void revContForwardSearchAndForward2Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for second-step-forward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward2Mismatch ( query, start2, len2,
                              revBwt, revOcc, revInverseSa0,
                              saL, saR,
                              revSaL, revSaR,
                              saCount, output,
                              maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void revContForwardSearchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for second-step-forward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward2MismatchAnd1Mismatch ( query, start2, len2,
                                          revBwt, revOcc, revInverseSa0,
                                          saL, saR,
                                          revSaL, revSaR,
                                          saCount, output,
                                          maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void revContForwardSearchAndForward3MismatchAnd2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for second-step-forward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start; i < start + len && saL <= saR; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        revSaL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaL, c, revInverseSa0, cum_start ) + 1;
        revSaR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revSaR + 1, c, revInverseSa0, cum_end );
        saR = saR + cum_start - cum_end;
        saL = saR - ( revSaR - revSaL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward3MismatchAnd2MismatchAnd1Mismatch ( query, start2, len2,
                revBwt, revOcc, revInverseSa0,
                saL, saR,
                revSaL, revSaR,
                saCount, output,
                maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void revForward1MismatchAndBackward3Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndBackward3Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for backward
                        revBwt, revOcc, revInverseSa0,
                        bwt, occ, inverseSa0, // for backward
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}


__forceinline__ __device__
void revForward1MismatchAndForward2Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for second-step-forward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndForward2Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for second-step-forward
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}


__forceinline__ __device__
void revForward1MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for second-step-forward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndForward2MismatchAnd1Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for second-step-forward
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void revForward1MismatchAndForward3MismatchAnd2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for second-step-forward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndForward3MismatchAnd2MismatchAnd1Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for second-step-forward
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void revForward2MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for second-step-forward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revForward1MismatchAndForward2MismatchAnd1Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for second-step-forward
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void contBackwardSearchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL =
            GPUBWTOccValueWithCumu ( bwt, occ, saL, c, inverseSa0, cum_start ) + 1;
        saR =
            GPUBWTOccValueWithCumu ( bwt, occ, saR + 1, c, inverseSa0, cum_end );
        revSaR = revSaR + cum_start - cum_end;
        revSaL = revSaR - ( saR - saL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward2MismatchAnd1Mismatch ( query, start2, len2,
                                          revBwt, revOcc, revInverseSa0,
                                          saL, saR,
                                          revSaL, revSaR,
                                          saCount, output,
                                          maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void backward1MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            mkR = occCountPEnd[ec];
            revMkR = revPr - occCountP[ec];
            revMkL = revMkR - ( mkR - mkL );

            if ( mkL <= mkR )
            {
                contBackwardSearchAndForward2MismatchAnd1Mismatch ( query, start, i - start - 1,
                        start2, len2,
                        bwt, occ, inverseSa0,
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
        revPr = revPr - occCountP[c];
    }
}

__forceinline__ __device__
void contBackwardSearchAndForward1MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint start3, uint len3,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL =
            GPUBWTOccValueWithCumu ( bwt, occ, saL, c, inverseSa0, cum_start ) + 1;
        saR =
            GPUBWTOccValueWithCumu ( bwt, occ, saR + 1, c, inverseSa0, cum_end );
        revSaR = revSaR + cum_start - cum_end;
        revSaL = revSaR - ( saR - saL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward1MismatchAndForward2MismatchAnd1Mismatch ( query, start2, len2, start3, len3,
                revBwt, revOcc, revInverseSa0,
                saL, saR,
                revSaL, revSaR,
                saCount, output,
                maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void contBackwardSearchAndBackward1MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint start3, uint len3,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL =
            GPUBWTOccValueWithCumu ( bwt, occ, saL, c, inverseSa0, cum_start ) + 1;
        saR =
            GPUBWTOccValueWithCumu ( bwt, occ, saR + 1, c, inverseSa0, cum_end );
        revSaR = revSaR + cum_start - cum_end;
        revSaL = revSaR - ( saR - saL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        backward1MismatchAndForward2MismatchAnd1Mismatch ( query, start2, len2, start3, len3,
                bwt, occ, inverseSa0,
                revBwt, revOcc, revInverseSa0,
                saL, saR,
                revSaL, revSaR,
                saCount, output,
                maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void backward1MismatchAndForward1MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint start3, uint len3,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            mkR = occCountPEnd[ec];
            revMkR = revPr - occCountP[ec];
            revMkL = revMkR - ( mkR - mkL );

            if ( mkL <= mkR )
            {
                contBackwardSearchAndForward1MismatchAndForward2MismatchAnd1Mismatch ( query, start, i - start - 1,
                        start2, len2,
                        start3, len3,
                        bwt, occ, inverseSa0,
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
        revPr = revPr - occCountP[c];
    }
}



__forceinline__ __device__
void backward1MismatchAndBackward1MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint start3, uint len3,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            mkR = occCountPEnd[ec];
            revMkR = revPr - occCountP[ec];
            revMkL = revMkR - ( mkR - mkL );

            if ( mkL <= mkR )
            {
                contBackwardSearchAndBackward1MismatchAndForward2MismatchAnd1Mismatch ( query, start, i - start - 1,
                        start2, len2,
                        start3, len3,
                        bwt, occ, inverseSa0,
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
        revPr = revPr - occCountP[c];
    }
}

__forceinline__ __device__
void revForward1MismatchAndBackward3MismatchAnd2Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2, // for backward
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint * bwt, uint * occ, uint inverseSa0, // for backward
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start; i < start + len && pl <= pr && saCount <= maxSARangesAllowed; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPl, occCountPStart );
        GPUBWTAllOccValue ( revBwt, revOcc, revInverseSa0, revPr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Forward manner
        for ( ec = 0; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            if ( c == ec ) { continue; } // TODO optimize as before

            revMkL = occCountPStart[ec] + 1;
            revMkR = occCountPEnd[ec];
            mkR = pr - occCountP[ec];
            mkL = mkR - ( revMkR - revMkL );

            if ( mkL <= mkR )
            {
                revContForwardSearchAndBackward3MismatchAnd2Mismatch ( query, i + 1, start + len - i - 1,
                        start2, len2, // for backward
                        revBwt, revOcc, revInverseSa0,
                        bwt, occ, inverseSa0, // for backward
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }
        }

        revPl = occCountPStart[c] + 1;
        revPr = occCountPEnd[c];
        pr = pr - occCountP[c];
        pl = pr - ( revPr - revPl );
    }
}

__forceinline__ __device__
void contBackwardSearchAndBackward3Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint saL, uint saR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL = GPUBWTOccValue ( bwt, occ, saL, c, inverseSa0 ) + 1;
        saR = GPUBWTOccValue ( bwt, occ, saR + 1, c, inverseSa0 );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        backward3Mismatch ( query, start2, len2,
                            bwt, occ, inverseSa0,
                            saL, saR,
                            saCount, output,
                            maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void backward1MismatchAndBackward3Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint pl, uint pr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            uint mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            uint mkR = occCountPEnd[ec];

            if ( mkL <= mkR )
            {
                contBackwardSearchAndBackward3Mismatch ( query, start, i - start - 1,
                        start2, len2,
                        bwt, occ, inverseSa0,
                        mkL, mkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
    }
}

__forceinline__ __device__
void contBackwardSearchAndForward3MismatchAnd2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint saL, uint saR,
        uint revSaL, uint revSaR,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    unsigned char c;
    uint i;
    uint cum_start, cum_end;

    for ( i = start + len; i > start && saL <= saR; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        saL =
            GPUBWTOccValueWithCumu ( bwt, occ, saL, c, inverseSa0, cum_start ) + 1;
        saR =
            GPUBWTOccValueWithCumu ( bwt, occ, saR + 1, c, inverseSa0, cum_end );
        revSaR = revSaR + cum_start - cum_end;
        revSaL = revSaR - ( saR - saL );
    }

    if ( saL <= saR && saCount <= maxSARangesAllowed )
    {
        revForward3MismatchAnd2MismatchAnd1Mismatch ( query, start2, len2,
                revBwt, revOcc, revInverseSa0,
                saL, saR,
                revSaL, revSaR,
                saCount, output,
                maxSARangesAllowed, strand, accumMismatches );
    }
}

__forceinline__ __device__
void backward2MismatchAndForward2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            mkR = occCountPEnd[ec];
            revMkR = revPr - occCountP[ec];
            revMkL = revMkR - ( mkR - mkL );

            if ( mkL <= mkR )
            {
                backward1MismatchAndForward2MismatchAnd1Mismatch ( query, start, i - start - 1,
                        start2, len2,
                        bwt, occ, inverseSa0,
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
        revPr = revPr - occCountP[c];
    }
}

__forceinline__ __device__
void backward1MismatchAndForward3MismatchAnd2MismatchAnd1Mismatch ( uint * query, uint start, uint len,
        uint start2, uint len2,
        uint * bwt, uint * occ, uint inverseSa0,
        uint * revBwt, uint * revOcc, uint revInverseSa0,
        uint pl, uint pr,
        uint revPl, uint revPr,
        uint & saCount, uint * output,
        uint maxSARangesAllowed, uint strand, uint accumMismatches )
{
    uint mkL, mkR, revMkL, revMkR;
    uint __align__(16) occCountPStart[ALPHABET_SIZE];
    uint __align__(16) occCountPEnd[ALPHABET_SIZE];
    uint __align__(16) occCountP[ALPHABET_SIZE];
    unsigned char c;
    unsigned char ec;
    uint i;

    for ( i = start + len; i > start && pl <= pr && saCount <= maxSARangesAllowed; --i )
    {
        // note that we use i-1 here to prevent counter overflow
        c = ( query[ ( ( i - 1 ) / CHAR_PER_WORD ) * 32] >>
              ( ( i - 1 ) % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pl, occCountPStart );
        GPUBWTAllOccValue ( bwt, occ, inverseSa0, pr + 1, occCountPEnd );
        int k;
        occCountP[ALPHABET_SIZE - 1] = 0;

        for ( k = ALPHABET_SIZE - 2; k >= 0; --k )
        {
            occCountP[k] = occCountP[k + 1] + occCountPEnd[k + 1] - occCountPStart[k + 1];
        }

        // Backward manner
        for ( ec = c ? 0 : 1; ec < ALPHABET_SIZE && saCount <= maxSARangesAllowed; ++ec )
        {
            //      if (c == ec) continue; // TODO optimize as before
            mkL = occCountPStart[ec] + 1; // compute SA range if query[i] was ec
            mkR = occCountPEnd[ec];
            revMkR = revPr - occCountP[ec];
            revMkL = revMkR - ( mkR - mkL );

            if ( mkL <= mkR )
            {
                contBackwardSearchAndForward3MismatchAnd2MismatchAnd1Mismatch ( query, start, i - start - 1,
                        start2, len2,
                        bwt, occ, inverseSa0,
                        revBwt, revOcc, revInverseSa0,
                        mkL, mkR,
                        revMkL, revMkR,
                        saCount, output,
                        maxSARangesAllowed, strand, accumMismatches + 1 );
            }

            if ( c == ec + 1 )
            { ec = c; }
        }

        pl = occCountPStart[c] + 1;
        pr = occCountPEnd[c];
        revPr = revPr - occCountP[c];
    }
}

__forceinline__ __device__
void matchQueryCaseA_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    int i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case A   Backward search
    // 1. exact cell4+5
    // 2. 4/3/2/1/0-mismatch cell1+2+3
    //=======================================
    l = 1;
    r = textLength;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );

    // exact cell4+5
    for ( i = readLength - 1;
            i >= cell1 + cell2 + cell3 && l <= r; --i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        l = GPUBWTOccValue ( bwt, occ, l, c, inverseSa0 ) + 1;
        r = GPUBWTOccValue ( bwt, occ, r + 1, c, inverseSa0 );
    }

    // 4/3/2/1/0-mismatch cell1+2+3
    if ( l <= r )
    {
        // search range = query[0, cell1 + cell2 + cell3 - 1]
        backward4MismatchAnd3MismatchAnd2MismatchAnd1MismatchAndExact ( query, 0, cell1 + cell2 + cell3,
                bwt, occ, inverseSa0,
                l, r,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}


__forceinline__ __device__
void matchQueryCaseB_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    int i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case B   Forward search
    // 1. cell1+2+3 exact
    // 2. 4/3/2/1-mismatch cell4+5
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // cell1+2+3 exact
    for ( i = 0; i < cell1 + cell2 + cell3 && l <= r; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint start, end;
        revL = GPUBWTOccValueWithCumu ( revBwt, revOcc, revL,
                                        c, revInverseSa0, start ) + 1;
        revR = GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1,
                                        c, revInverseSa0, end );
        r = r + start - end;
        l = r - ( revR - revL );
    }

    // 4-mismatch cell4+5
    if ( l <= r )
    {
        // search range = query[cell1 + cell2 + cell3, QUERY_LENGTH-1]
        revForward4MismatchAnd3MismatchAnd2MismatchAnd1Mismatch ( query, cell1 + cell2 + cell3, cell4 + cell5,
                revBwt, revOcc, revInverseSa0,
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}



__forceinline__ __device__
void matchQueryCaseC_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    int i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case C
    // 1. exact cell1 (forward)
    // 2. 1-mismatch cell2+3 (forward)
    // 3. 3/2/1-mismatch cell4+5 (forward)
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // exact cell1 (forward)
    for ( i = 0; i < cell1 && l <= r; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint cum_start, cum_end;
        revL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
        revR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
        r = r + cum_start - cum_end;
        l = r - ( revR - revL );
    }

    // 1-mismatch cell2+3 (forward)
    // 3/2/1-mismatch cell4+5 (forward)
    if ( l <= r )
    {
        revForward1MismatchAndForward3MismatchAnd2MismatchAnd1Mismatch ( query, cell1, cell2 + cell3,
                cell1 + cell2 + cell3, cell4 + cell5, // for second-step-forward
                revBwt, revOcc, revInverseSa0,
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}



__forceinline__ __device__
void matchQueryCaseD_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    int i;
    unsigned char c;
    //=======================================
    //  FOR FOUR MISMATCH
    //=======================================
    // Case D
    // 1. exact cell2+3 (forward)
    // 2. 1-mismatch cell1 (backward)
    // 3. 3/2/1-mismatch cell4+5 (forward)
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // exact cell2+3 (forward)
    for ( i = cell1; i < cell1 + cell2 + cell3 && l <= r; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint cum_start, cum_end;
        revL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
        revR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
        r = r + cum_start - cum_end;
        l = r - ( revR - revL );
    }

    // 1-mismatch cell1 (backward)
    // 3/2/1-mismatch cell4+5 (forward)
    if ( l <= r )
    {
        backward1MismatchAndForward3MismatchAnd2MismatchAnd1Mismatch ( query, 0, cell1,
                cell1 + cell2 + cell3, cell4 + cell5,
                bwt, occ, inverseSa0,
                revBwt, revOcc, revInverseSa0,
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}



__forceinline__ __device__
void matchQueryCaseE_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    int i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case E
    // 1. exact cell1 (forward)
    // 2. 2-mismatch cell2+3 (forward)
    // 3. 1/2-mismatch cell4+5 (forward)
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // exact cell1 (forward)
    for ( i = 0; i < cell1 && l <= r; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint cum_start, cum_end;
        revL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
        revR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
        r = r + cum_start - cum_end;
        l = r - ( revR - revL );
    }

    // 2-mismatch cell2+3 (forward)
    // 1/2-mismatch cell4+5 (forward)
    if ( l <= r )
    {
        revForward2MismatchAndForward2MismatchAnd1Mismatch ( query, cell1, cell2 + cell3,
                cell1 + cell2 + cell3, cell4 + cell5, // for second-step-forward
                revBwt, revOcc, revInverseSa0,
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}


__forceinline__ __device__
void matchQueryCaseF_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    int i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case F
    // 1. exact cell2+3 (forward)
    // 2. 2-mismatch cell1 (backward)
    // 3. 1/2-mismatch cell4+5 (forward)
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // exact cell3 (forward)
    for ( i = cell1; i < cell1 + cell2 + cell3 && l <= r; ++i )
    {
        // TODO compute c here
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint cum_start, cum_end;
        revL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
        revR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
        r = r + cum_start - cum_end;
        l = r - ( revR - revL );
    }

    // 2-mismatch cell1 (backward)
    // 1/2-mismatch cell4+5 (forward)
    if ( l <= r )
    {
        backward2MismatchAndForward2MismatchAnd1Mismatch ( query, 0, cell1,
                cell1 + cell2 + cell3, cell4 + cell5,
                bwt, occ, inverseSa0,
                revBwt, revOcc, revInverseSa0,
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}


__forceinline__ __device__
void matchQueryCaseG_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    uint i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case G
    // 1. exact cell2 (forward)
    // 2. 1-mismatch cell1 (backward)
    // 3. 1-mismatch cell3 (forward)
    // 4. 1/2-mismatch cell4+5 (forward)
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // exact cell2 (forward)
    for ( i = cell1; i < cell1 + cell2 && l <= r; ++i )
    {
        // TODO compute c here
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint cum_start, cum_end;
        revL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
        revR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
        r = r + cum_start - cum_end;
        l = r - ( revR - revL );
    }

    // 1-mismatch cell1 (backward)
    // 1-mismatch cell3 (forward)
    // 1/2-mismatch cell4+5 (forward)
    if ( l <= r )
    {
        backward1MismatchAndForward1MismatchAndForward2MismatchAnd1Mismatch ( query, 0, cell1,
                cell1 + cell2, cell3,
                cell1 + cell2 + cell3, cell4 + cell5,
                bwt, occ, inverseSa0,
                revBwt, revOcc, revInverseSa0,
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}



__forceinline__ __device__
void matchQueryCaseH_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    uint i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case H
    // 1. exact cell3 (forward)
    // 2. 1-mismatch cell2 (backward)
    // 3. 1-mismatch cell1 (backward)
    // 4. 1/2-mismatch cell4+5 (forward)
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // exact cell3 (forward)
    for ( i = cell1 + cell2; i < cell1 + cell2 + cell3 && l <= r; ++i )
    {
        // TODO compute c here
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint cum_start, cum_end;
        revL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
        revR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
        r = r + cum_start - cum_end;
        l = r - ( revR - revL );
    }

    // 1-mismatch cell2 (backward)
    // 1-mismatch cell1 (backward)
    // 1/2-mismatch cell4+5 (forward)
    if ( l <= r )
    {
        backward1MismatchAndBackward1MismatchAndForward2MismatchAnd1Mismatch ( query, cell1, cell2,
                0, cell1,
                cell1 + cell2 + cell3, cell4 + cell5,
                bwt, occ, inverseSa0,
                revBwt, revOcc, revInverseSa0,
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}




__forceinline__ __device__
void matchQueryCaseI_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint revL, revR;
    uint i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case I
    // 1. exact cell4 (forward)
    // 2. 1-mismatch cell5 (forward)
    // 3. 3-mismatch cell1+2+3 (backward)
    //=======================================
    l = 0;
    r = textLength;
    revL = 0;
    revR = r;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    uint cell5 = readLength - cell1 - cell2 - cell3 - cell4;

    // exact cell4 (forward)
    for ( i = cell1 + cell2 + cell3; i < cell1 + cell2 + cell3 + cell4 && l <= r; ++i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        uint cum_start, cum_end;
        revL =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
        revR =
            GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
        r = r + cum_start - cum_end;
        l = r - ( revR - revL );
    }

    // 1-mismatch cell5 (forward)
    // 3-mismatch cell1+2+3 (backward)
    if ( l <= r )
    {
        revForward1MismatchAndBackward3Mismatch ( query, cell1 + cell2 + cell3 + cell4, cell5,
                0, cell1 + cell2 + cell3, // for backward
                revBwt, revOcc, revInverseSa0,
                bwt,  occ,  inverseSa0, // for backward
                l, r,
                revL, revR,
                saCount, output,
                maxSARangesAllowed, strand, 0 );
    }
}

__forceinline__ __device__
void matchQueryCaseJ_4mismatch ( uint * query, uint readLength,
                                 uint * bwt, uint * occ, uint inverseSa0,
                                 uint * revBwt, uint * revOcc, uint revInverseSa0,
                                 uint & saCount, uint * output, unsigned int textLength,
                                 uint maxSARangesAllowed, uint strand )
{
    uint l, r;
    uint i;
    unsigned char c;
    //======================================================
    // FOR FOUR MISMATCH
    //======================================================
    // Case J
    // 1. exact cell5 (backward)
    // 2. 1-mismatch cell4 (backward)
    // 3. 3-mismatch cell1+2+3 (backward)
    //=======================================
    l = 0;
    r = textLength;
    uint cell1 = ( int ) ( readLength * SIZE_A_RATIO );
    uint cell2 = ( int ) ( readLength * SIZE_B_RATIO );
    uint cell3 = ( int ) ( readLength * SIZE_C_RATIO );
    uint cell4 = ( int ) ( readLength * SIZE_D_RATIO );
    l = 1;
    r = textLength;

    // exact cell5 (backward)
    for ( i = readLength - 1;
            i >= cell1 + cell2 + cell3 + cell4 && l <= r; --i )
    {
        c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
        l = GPUBWTOccValue ( bwt, occ, l, c, inverseSa0 ) + 1;
        r = GPUBWTOccValue ( bwt, occ, r + 1, c, inverseSa0 );
    }

    // 1-mismatch cell4 (backward)
    // 3-mismatch cell1+2+3 (backward)
    if ( l <= r )
    {
        backward1MismatchAndBackward3Mismatch ( query, cell1 + cell2 + cell3, cell4,
                                                0, cell1 + cell2 + cell3,
                                                bwt, occ, inverseSa0,
                                                l, r,
                                                saCount, output, maxSARangesAllowed, strand, 0 );
    }
}



__forceinline__ __device__
void matchQueryCaseA ( uint * query, uint readLength,
                       uint * bwt, uint * occ, uint inverseSa0,
                       uint * revBwt, uint * revOcc, uint revInverseSa0,
                       uint & saCount, uint * output, uint textLength,
                       uint maxSARangesAllowed, uint strand, uint numMismatch,
                       bool isExactNumMismatch )
{
    uint l, r;
    int i;
    unsigned char c;

    if ( numMismatch == 0 )
    {
        //======================================================
        // FOR EXACT MATCH
        //======================================================
        l = 1;
        r = textLength;
        // backward search with BWT
        contBackwardSearch ( query, 0, readLength,
                             bwt, occ, inverseSa0,
                             l, r,
                             saCount, output,
                             maxSARangesAllowed,
                             strand, 0 );
    }

    if ( numMismatch == 1 )
    {
        //======================================================
        // FOR ONE MISMATCH
        //======================================================
        // Case A   Backward search
        // 1. exact cellY
        // 2. 1-mismatch cellX
        //=======================================
        l = 1;
        r = textLength;
        uint sizeX = ( int ) ( readLength * 0.5 );

        // backward search with BWT in cellY
        for ( i = readLength - 1; i >= sizeX && l <= r; --i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            l = GPUBWTOccValue ( bwt, occ, l, c, inverseSa0 ) + 1;
            r = GPUBWTOccValue ( bwt, occ, r + 1, c, inverseSa0 );
        }

        // 1/0-mismatch cellX
        if ( l <= r )
        {
            if ( isExactNumMismatch )
                backward1Mismatch ( query, 0, sizeX, bwt, occ, inverseSa0,
                                    l, r, saCount, output, maxSARangesAllowed, strand, 0 );
            else
                backward1MismatchAndExact ( query, 0, sizeX, bwt, occ, inverseSa0,
                                            l, r, saCount, output, maxSARangesAllowed, strand, 0 );
        }
    }
    else if ( numMismatch == 2 )
    {
        //======================================================
        // FOR TWO MISMATCH
        //======================================================
        // Case A   Backward search
        // 1. exact cellZ
        // 2. 2/1/0-mismatch cellX+Y
        //=======================================
        l = 1;
        r = textLength;
        uint sizeX = ( int ) ( readLength * SIZE_X_RATIO );
        uint sizeY = ( int ) ( readLength * SIZE_Y_RATIO );

        // backward search with BWT in cellZ
        for ( i = readLength - 1; i >= sizeX + sizeY && l <= r; --i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            l = GPUBWTOccValue ( bwt, occ, l, c, inverseSa0 ) + 1;
            r = GPUBWTOccValue ( bwt, occ, r + 1, c, inverseSa0 );
        }

        // 2. 2/1/0-mismatch cellX+Y
        if ( l <= r )
        {
            backward2MismatchAnd1MismatchAndExact ( query, 0, sizeX + sizeY,
                                                    bwt, occ, inverseSa0,
                                                    l, r,
                                                    saCount, output,
                                                    maxSARangesAllowed, strand, 0 );
        }
    }
    else if ( numMismatch == 3 )
    {
        //======================================================
        // FOR THREE MISMATCH
        //======================================================
        // Case A   Backward search
        // 1. exact cell3+4
        // 2. 3/2/1/0-mismatch cell1+2
        //=======================================
        l = 1;
        r = textLength;
        uint cell1 = ( int ) ( readLength * SIZE_1_RATIO );
        uint cell2 = ( int ) ( readLength * SIZE_2_RATIO );

        // backward search with BWT in cellZ
        for ( i = readLength - 1;
                i >= cell1 + cell2 && l <= r; --i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            l = GPUBWTOccValue ( bwt, occ, l, c, inverseSa0 ) + 1;
            r = GPUBWTOccValue ( bwt, occ, r + 1, c, inverseSa0 );
        }

        // 3/2/1/0-mismatch cell1+2
        if ( l <= r )
        {
            // search range = query[0, cell1 + cell2 - 1]
            backward3MismatchAnd2MismatchAnd1MismatchAndExact ( query, 0, cell1 + cell2,
                    bwt, occ, inverseSa0,
                    l, r,
                    saCount, output,
                    maxSARangesAllowed, strand, 0 );
        }
    }
}

__forceinline__ __device__
void matchQueryCaseB ( uint * query, uint readLength,
                       uint * bwt, uint * occ, uint inverseSa0,
                       uint * revBwt, uint * revOcc, uint revInverseSa0,
                       uint & saCount, uint * output, uint textLength,
                       uint maxSARangesAllowed, uint strand, uint numMismatch )
{
    uint l, r;
    uint revL, revR;
    int i;
    unsigned char c;

    if ( numMismatch == 1 )
    {
        //======================================================
        // FOR ONE MISMATCH
        //======================================================
        // Case B   Forward search
        // 1. cellX exact
        // 2. 1-mismatch cellY
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint sizeX = ( int ) ( readLength * 0.5 );
        // uint sizeY = (int)(readLength * 0.5);
        uint sizeY = readLength - sizeX;

        // forward search with BWT until the end of forward depth section
        for ( i = 0; i < sizeX && l <= r; ++i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint cum_start, cum_end;
            revL =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
            revR =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
            r = r + cum_start - cum_end;
            l = r - ( revR - revL );
        }

        // 2. 1-mismatch cellY
        if ( l <= r )
        {
            revForward1Mismatch ( query, sizeX, sizeY,
                                  revBwt, revOcc, revInverseSa0,
                                  l, r,
                                  revL, revR,
                                  saCount, output,
                                  maxSARangesAllowed, strand, 0 );
        }
    }
    else if ( numMismatch == 2 )
    {
        //======================================================
        // FOR TWO MISMATCH
        //======================================================
        // Case B   Forward search
        // 1. cellX+Y exact
        // 2. 2/1-mismatch cellZ
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint sizeX = ( int ) ( readLength * SIZE_X_RATIO );
        uint sizeY = ( int ) ( readLength * SIZE_Y_RATIO );
        uint sizeZ = readLength - sizeX - sizeY;

        // forward search with BWT until the end of forward depth section
        for ( i = 0; i < sizeX + sizeY && l <= r; ++i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint cum_start, cum_end;
            revL =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
            revR =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
            r = r + cum_start - cum_end;
            l = r - ( revR - revL );
        }

        // 2 errors in cellZ
        if ( l <= r )
        {
            revForward2MismatchAnd1Mismatch ( query, sizeX + sizeY, sizeZ,
                                              revBwt, revOcc, revInverseSa0,
                                              l, r,
                                              revL, revR,
                                              saCount, output,
                                              maxSARangesAllowed, strand, 0 );
        }
    }
    else if ( numMismatch == 3 )
    {
        //======================================================
        // FOR THREE MISMATCH
        //======================================================
        // Case B   Forward search
        // 1. cell1+2 exact
        // 2. 3/2/1-mismatch cell3+4
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint cell1 = ( int ) ( readLength * SIZE_1_RATIO );
        uint cell2 = ( int ) ( readLength * SIZE_2_RATIO );
        uint cell3 = ( int ) ( readLength * SIZE_3_RATIO );
        uint cell4 = readLength - cell1 - cell2 - cell3;

        // forward search with BWT until the end of forward depth section
        for ( i = 0; i < cell1 + cell2 && l <= r; ++i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint start, end;
            revL = GPUBWTOccValueWithCumu ( revBwt, revOcc, revL,
                                            c, revInverseSa0, start ) + 1;
            revR = GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1,
                                            c, revInverseSa0, end );
            r = r + start - end;
            l = r - ( revR - revL );
        }

        // 3/2/1-mismatch cell3+4
        if ( l <= r )
        {
            // search range = query[cell1 + cell2, QUERY_LENGTH-1]
            revForward3MismatchAnd2MismatchAnd1Mismatch ( query, cell1 + cell2, cell3 + cell4,
                    revBwt, revOcc, revInverseSa0,
                    l, r,
                    revL, revR,
                    saCount, output,
                    maxSARangesAllowed, strand, 0 );
        }
    }
}

__forceinline__ __device__
void matchQueryCaseC ( uint * query, uint readLength,
                       uint * bwt, uint * occ, uint inverseSa0,
                       uint * revBwt, uint * revOcc, uint revInverseSa0,
                       uint & saCount, uint * output, uint textLength,
                       uint maxSARangesAllowed, uint strand, uint numMismatch )
{
    uint l, r;
    uint revL, revR;
    uint i;
    unsigned char c;

    if ( numMismatch == 2 )
    {
        //======================================================
        // FOR TWO MISMATCH
        //======================================================
        // Case C
        // 1. cellX (forward)
        // 2. 1-mismatch cellY (forward)
        // 3. 1-mismatch cellZ (forward)
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint sizeX = ( int ) ( readLength * SIZE_X_RATIO );
        uint sizeY = ( int ) ( readLength * SIZE_Y_RATIO );
        uint sizeZ = readLength - sizeX - sizeY;

        // forward search with BWT until the end of forward depth section
        for ( i = 0; i < sizeX && l <= r; ++i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint cum_start, cum_end;
            revL =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
            revR =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
            r = r + cum_start - cum_end;
            l = r - ( revR - revL );
        }

        // 2. 1-mismatch cellY (forward)
        // 3. 1-mismatch cellZ (forward)
        if ( l <= r )
        {
            revForward1MismatchAndForward1Mismatch ( query, sizeX, sizeY,
                    sizeX + sizeY, sizeZ, // for backward
                    revBwt, revOcc, revInverseSa0,
                    l, r,
                    revL, revR,
                    saCount, output,
                    maxSARangesAllowed, strand, 0 );
        }
    }
    else if ( numMismatch == 3 )
    {
        //======================================================
        // FOR THREE MISMATCH
        //======================================================
        // Case C
        // 1. exact cell1 (forward)
        // 2. 1-mismatch cell2 (forward)
        // 3. 2-mismatch cell3+4 (forward)
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint cell1 = ( int ) ( readLength * SIZE_1_RATIO );
        uint cell2 = ( int ) ( readLength * SIZE_2_RATIO );
        uint cell3 = ( int ) ( readLength * SIZE_3_RATIO );
        uint cell4 = readLength - cell1 - cell2 - cell3;

        // exact cell1 (forward)
        for ( i = 0; i < cell1 && l <= r; ++i )
        {
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint cum_start, cum_end;
            revL =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
            revR =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
            r = r + cum_start - cum_end;
            l = r - ( revR - revL );
        }

        // 1-mismatch cell2 (forward)
        // 2-mismatch cell3+4 (forward)

        if ( l <= r )
        {
            revForward1MismatchAndForward2Mismatch ( query, cell1, cell2,
                    cell1 + cell2, cell3 + cell4,
                    revBwt, revOcc, revInverseSa0,
                    l, r,
                    revL, revR,
                    saCount, output,
                    maxSARangesAllowed, strand, 0 );
        }
    }
}

__forceinline__ __device__
void matchQueryCaseD ( uint * query, uint readLength,
                       uint * bwt, uint * occ, uint inverseSa0,
                       uint * revBwt, uint * revOcc, uint revInverseSa0,
                       uint & saCount, uint * output, uint textLength,
                       uint maxSARangesAllowed, uint strand, uint numMismatch )
{
    uint l, r;
    uint revL, revR;
    uint i;
    unsigned char c;

    if ( numMismatch == 2 )
    {
        //======================================================
        // FOR TWO MISMATCH
        //======================================================
        // Case D
        // 1. cellY (forward)
        // 2. 1-mismatch cellZ (forward)
        // 3. 1-mismatch cellX (backward)
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint sizeX = ( int ) ( readLength * SIZE_X_RATIO );
        uint sizeY = ( int ) ( readLength * SIZE_Y_RATIO );
        uint sizeZ = readLength - sizeX - sizeY;

        // forward search with BWT until the end of forward depth section
        for ( i = sizeX; i < sizeX + sizeY && l <= r; ++i )
        {
            // TODO compute c here
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint cum_start, cum_end;
            revL =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
            revR =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
            r = r + cum_start - cum_end;
            l = r - ( revR - revL );
        }

        // 2. 1-mismatch cellZ (forward)
        // 3. 1-mismatch cellX (backward)
        if ( l <= r )
        {
            revForward1MismatchAndBackward1Mismatch ( query, sizeX + sizeY, sizeZ,
                    0, sizeX, // for backward
                    revBwt, revOcc, revInverseSa0,
                    bwt, occ, inverseSa0, // for backward
                    l, r,
                    revL, revR,
                    saCount, output,
                    maxSARangesAllowed, strand, 0 );
        }
    }
    else if ( numMismatch == 3 )
    {
        //=======================================
        //  FOR THREE MISMATCH
        //=======================================
        // Case D
        // 1. exact cell3 (forward)
        // 2. 1-mismatch cell4 (forward)
        // 3. 2/1-mismatch cell1+2 (backward)
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint cell1 = ( int ) ( readLength * SIZE_1_RATIO );
        uint cell2 = ( int ) ( readLength * SIZE_2_RATIO );
        uint cell3 = ( int ) ( readLength * SIZE_3_RATIO );
        uint cell4 = readLength - cell1 - cell2 - cell3;

        // forward search with BWT until the end of forward depth section
        for ( i = cell1 + cell2; i < cell1 + cell2 + cell3 && l <= r; ++i )
        {
            // TODO compute c here
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint cum_start, cum_end;
            revL =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
            revR =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
            r = r + cum_start - cum_end;
            l = r - ( revR - revL );
        }

        // 1-mismatch cell4 (forward)
        // 2/1-mismatch cell1+2 (backward)
        if ( l <= r )
        {
            revForward1MismatchAndBackward2MismatchAnd1Mismatch ( query, cell1 + cell2 + cell3, cell4,
                    0, cell1 + cell2, // for backward
                    revBwt, revOcc, revInverseSa0,
                    bwt,  occ,  inverseSa0, // for backward
                    l, r,
                    revL, revR,
                    saCount, output,
                    maxSARangesAllowed, strand, 0 );
        }
    }
}

__forceinline__ __device__
void matchQueryCaseE ( uint * query, uint readLength,
                       uint * bwt, uint * occ, uint inverseSa0,
                       uint * revBwt, uint * revOcc, uint revInverseSa0,
                       uint & saCount, uint * output, uint textLength,
                       uint maxSARangesAllowed, uint strand, uint numMismatch )
{
    uint l, r;
    int i;
    unsigned char c;

    if ( numMismatch == 3 )
    {
        //======================================================
        // FOR THREE MISMATCH
        //======================================================
        // Case E
        // 1. exact cell4 (backward)
        // 2. 1-mismatch cell3 (backward)
        // 3. 2/1-mismatch cell1+2 (backward)
        //=======================================
        l = 1;
        r = textLength;
        uint cell1 = ( int ) ( readLength * SIZE_1_RATIO );
        uint cell2 = ( int ) ( readLength * SIZE_2_RATIO );
        uint cell3 = ( int ) ( readLength * SIZE_3_RATIO );

        // backward search with BWT in cell4
        for ( i = readLength - 1;
                i >= cell1 + cell2 + cell3 && l <= r; --i )
        {
            /* this optimization is like 10 ms faster only....
             if (i % 16 == 15) curWord = *(query + i / CHAR_PER_WORD * 32);
             c = (curWord >> (i % CHAR_PER_WORD * 2)) & 3; */
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            l = GPUBWTOccValue ( bwt, occ, l, c, inverseSa0 ) + 1;
            r = GPUBWTOccValue ( bwt, occ, r + 1, c, inverseSa0 );
        }

        // 1-mismatch cell3 (backward)
        // 2/1-mismatch cell1+2 (backward)
        if ( l <= r )
        {
            backward1MismatchAndBackward2MismatchAnd1Mismatch ( query, cell1 + cell2, cell3,
                    0, cell1 + cell2,
                    bwt, occ, inverseSa0,
                    l, r,
                    saCount, output, maxSARangesAllowed, strand, 0 );
        }
    }
}

__forceinline__ __device__
void matchQueryCaseF ( uint * query, uint readLength,
                       uint * bwt, uint * occ, uint inverseSa0,
                       uint * revBwt, uint * revOcc, uint revInverseSa0,
                       uint & saCount, uint * output, uint textLength,
                       uint maxSARangesAllowed, uint strand, uint numMismatch )
{
    uint l, r;
    uint revL, revR;
    uint i;
    unsigned char c;

    if ( numMismatch == 3 )
    {
        //======================================================
        // FOR THREE MISMATCH
        //======================================================
        // Case F
        // 1. exact cell2 (forward)
        // 2. 2-mismatch cell3+4 (forward)
        // 3. 1-mismatch cell1 (backward)
        //=======================================
        l = 0;
        r = textLength;
        revL = 0;
        revR = r;
        uint cell1 = ( int ) ( readLength * SIZE_1_RATIO );
        uint cell2 = ( int ) ( readLength * SIZE_2_RATIO );
        uint cell3 = ( int ) ( readLength * SIZE_3_RATIO );
        uint cell4 = readLength - cell1 - cell2 - cell3;

        // forward search with BWT until the end of forward depth section
        for ( i = cell1; i < cell1 + cell2 && l <= r; ++i )
        {
            // TODO compute c here
            c = ( query[ ( i / CHAR_PER_WORD ) * 32] >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
            uint cum_start, cum_end;
            revL =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revL, c, revInverseSa0, cum_start ) + 1;
            revR =
                GPUBWTOccValueWithCumu ( revBwt, revOcc, revR + 1, c, revInverseSa0, cum_end );
            r = r + cum_start - cum_end;
            l = r - ( revR - revL );
        }

        // 1-mismatch cell1 (backward)
        // 2-mismatch cell3+4 (forward)

        if ( l <= r )
        {
            backward1MismatchAndForward2Mismatch ( query, 0, cell1,
                                                   cell1 + cell2, cell3 + cell4,
                                                   bwt, occ, inverseSa0,
                                                   revBwt, revOcc, revInverseSa0,
                                                   l, r,
                                                   revL, revR,
                                                   saCount, output,
                                                   maxSARangesAllowed, strand, 0 );
        }
    }
}


// entry point of kernel
__global__ void kernel ( uint whichCase, uint * queries, uint * readLengths, uint numQueries,
                         uint wordPerQuery,
                         uint * bwt, uint * occ, uint inverseSa0,
                         uint * revBwt, uint * revOcc, uint revInverseSa0,
                         uint textLength,
                         uint * answers,
                         bool * isBad, uint round, uint numMismatch,
                         uint sa_range_allowed, uint wordPerAnswer,
                         bool isExactNumMismatch )
{
    uint queryId = ( blockIdx.x * THREADS_PER_BLOCK + threadIdx.x ) *
                   QUERIES_PER_THREAD;

    if ( queryId < numQueries )
    {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[KERNEL] Starting..\n" );
#endif
        ///////////////// Process read //////////////////
        ullint idx = queryId % 32;
        ullint queryOffset = queryId / 32 * 32 * wordPerQuery + idx;
        uint * query = queries + queryOffset;
        ullint answerOffset = queryId / 32 * 32 * wordPerAnswer + idx;
        uint * answer = answers + answerOffset;
        uint readLength = readLengths[queryId];

        for ( uint i = 0; i < wordPerAnswer; ++i )
        { answer[i * 32] = 0xFFFFFFFF; }

        uint numSARanges = 0;
#ifndef BGS_DISABLE_NEGATIVE_STRAND
        uint strand = ( round > 0 ) ? 0 : ( whichCase % 2 );
#else
        uint strand = 0;
#endif

        if ( ( round > 0 ) || ( isBad[queryId] == 0 ) )
        {
            // match the positive strand
            if ( whichCase == 0 )
            {
#ifndef BGS_DISABLE_CASE_A
                matchQueryCaseA ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch, isExactNumMismatch );
#endif
            }
            else if ( whichCase == 1 )
            {
#ifndef BGS_DISABLE_CASE_B
                matchQueryCaseB ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 2 )
            {
#ifndef BGS_DISABLE_CASE_C
                matchQueryCaseC ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 3 )
            {
#ifndef BGS_DISABLE_CASE_D
                matchQueryCaseD ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 4 )
            {
#ifndef BGS_DISABLE_CASE_E
                matchQueryCaseE ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 5 )
            {
#ifndef BGS_DISABLE_CASE_F
                matchQueryCaseF ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }

#ifndef BGS_DISABLE_NEGATIVE_STRAND
            // reverse the read
            uint leftWordIndex = 0, rightWordIndex = ( readLength - 1 ) / CHAR_PER_WORD;
            uint leftWord = query[leftWordIndex * 32];
            uint rightWord = query[rightWordIndex * 32];

            for ( uint i = 0, j = readLength - 1; i <= j; ++i, --j )
            {
                // check if need to move to next word
                if ( i / CHAR_PER_WORD != leftWordIndex )
                {
                    // write back leftword
                    query[leftWordIndex * 32] = leftWord;
                    // load next leftword
                    leftWordIndex++;
                    leftWord = query[leftWordIndex * 32];
                }

                if ( j / CHAR_PER_WORD != rightWordIndex )
                {
                    // write back rightword
                    query[rightWordIndex * 32] = rightWord;
                    // load next rightword
                    rightWordIndex--;
                    rightWord = query[rightWordIndex * 32];
                }

                // swap left and right characters
                unsigned char leftChar = ( leftWord >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
                unsigned char rightChar = ( rightWord >> ( j % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
                leftWord ^= ( leftChar ^ ( _soap3DnaComplement[rightChar] ) ) << ( i % CHAR_PER_WORD * BIT_PER_CHAR );
                rightWord ^= ( ( _soap3DnaComplement[leftChar] ) ^ rightChar ) << ( j % CHAR_PER_WORD * BIT_PER_CHAR );
            }

            // write back
            if ( leftWordIndex == rightWordIndex )
            {
                uint numLeftBits = ( ( readLength - 1 ) / 2 % CHAR_PER_WORD + 1 ) * BIT_PER_CHAR;
                uint numRightBits = 32 - numLeftBits;
                query[leftWordIndex * 32] = ( ( leftWord << numRightBits ) >> numRightBits ) |
                                            ( ( rightWord >> numLeftBits ) << numLeftBits );
            }
            else
            {
                query[leftWordIndex * 32] = leftWord;
                query[rightWordIndex * 32] = rightWord;
            }

            strand = 1 - strand;

            // match the negative strand
            if ( whichCase == 0 )
            {
#ifndef BGS_DISABLE_CASE_A
                matchQueryCaseA ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch, isExactNumMismatch );
#endif
            }
            else if ( whichCase == 1 )
            {
#ifndef BGS_DISABLE_CASE_B
                matchQueryCaseB ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 2 )
            {
#ifndef BGS_DISABLE_CASE_C
                matchQueryCaseC ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 3 )
            {
#ifndef BGS_DISABLE_CASE_D
                matchQueryCaseD ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 4 )
            {
#ifndef BGS_DISABLE_CASE_E
                matchQueryCaseE ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }
            else if ( whichCase == 5 )
            {
#ifndef BGS_DISABLE_CASE_F
                matchQueryCaseF ( query, readLength,
                                  bwt, occ, inverseSa0,
                                  revBwt, revOcc, revInverseSa0,
                                  numSARanges, answer, textLength,
                                  sa_range_allowed, strand, numMismatch );
#endif
            }

#endif

            // write error code if there are too many SA ranges
            if ( numSARanges == 0 )
            {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
                printf ( "[KERNEL] This query does not return result\n" );
#endif
                * ( answer ) = 0xFFFFFFFD;
            }
            else if ( numSARanges > sa_range_allowed )
            {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
                printf ( "[KERNEL] This query returns too many result\n" );
#endif
                * ( answer ) = 0xFFFFFFFE;

                if ( round == 0 )
                { isBad[queryId] = 1; }
            }
        }
        else
        {
            numSARanges = sa_range_allowed + 1;
            * ( answer ) = 0xFFFFFFFE;
        }

#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", answer );
        printf ( "[KERNEL] Sending out %u SA ranges..\n", numSARanges );

        for ( uint i = 0; i < wordPerQuery / 2; ++i )
        {
            printf ( "%d: %u %u\n", i, * ( answer + ( i * 2 ) * 32 ), * ( answer + ( i * 2 + 1 ) * 32 ) );
        }

#endif
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[KERNEL] Exiting successfully.\n" );
#endif
    }
}

// entry point of kernel
__global__ void kernel_4mismatch_1 ( uint whichCase, uint * queries, uint * readLengths, uint numQueries,
                                     uint wordPerQuery,
                                     uint * bwt, uint * occ, uint inverseSa0,
                                     uint * revBwt, uint * revOcc, uint revInverseSa0,
                                     uint textLength,
                                     uint * answers,
                                     bool * isBad, uint round, uint sa_range_allowed,
                                     uint wordPerAnswer, bool isExactNumMismatch )
{
    uint queryId = ( blockIdx.x * THREADS_PER_BLOCK + threadIdx.x ) *
                   QUERIES_PER_THREAD;

    if ( queryId < numQueries )
    {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[KERNEL] Starting..\n" );
#endif
        ///////////////// Process read //////////////////
        ullint idx = queryId % 32;
        ullint queryOffset = queryId / 32 * 32 * wordPerQuery + idx;
        uint * query = queries + queryOffset;
        ullint answerOffset = queryId / 32 * 32 * wordPerAnswer + idx;
        uint * answer = answers + answerOffset;
        uint readLength = readLengths[queryId];

        for ( uint i = 0; i < wordPerAnswer; ++i )
        { answer[i * 32] = 0xFFFFFFFF; }

        uint numSARanges = 0;
#ifndef BGS_DISABLE_NEGATIVE_STRAND
        uint strand = ( round > 0 ) ? 0 : ( whichCase % 2 );
#else
        uint strand = 0;
#endif

        if ( ( round > 0 ) || ( isBad[queryId] == 0 ) )
        {
            // match the positive strand
            if ( whichCase == 0 )
            {
#ifndef BGS_DISABLE_CASE_A
                matchQueryCaseA_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 1 )
            {
#ifndef BGS_DISABLE_CASE_B
                matchQueryCaseB_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 2 )
            {
#ifndef BGS_DISABLE_CASE_C
                matchQueryCaseC_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 3 )
            {
#ifndef BGS_DISABLE_CASE_D
                matchQueryCaseD_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 4 )
            {
#ifndef BGS_DISABLE_CASE_E
                matchQueryCaseE_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }

#ifndef BGS_DISABLE_NEGATIVE_STRAND
            // reverse the read
            uint leftWordIndex = 0, rightWordIndex = ( readLength - 1 ) / CHAR_PER_WORD;
            uint leftWord = query[leftWordIndex * 32];
            uint rightWord = query[rightWordIndex * 32];

            for ( uint i = 0, j = readLength - 1; i <= j; ++i, --j )
            {
                // check if need to move to next word
                if ( i / CHAR_PER_WORD != leftWordIndex )
                {
                    // write back leftword
                    query[leftWordIndex * 32] = leftWord;
                    // load next leftword
                    leftWordIndex++;
                    leftWord = query[leftWordIndex * 32];
                }

                if ( j / CHAR_PER_WORD != rightWordIndex )
                {
                    // write back rightword
                    query[rightWordIndex * 32] = rightWord;
                    // load next rightword
                    rightWordIndex--;
                    rightWord = query[rightWordIndex * 32];
                }

                // swap left and right characters
                unsigned char leftChar = ( leftWord >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
                unsigned char rightChar = ( rightWord >> ( j % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
                leftWord ^= ( leftChar ^ ( _soap3DnaComplement[rightChar] ) ) << ( i % CHAR_PER_WORD * BIT_PER_CHAR );
                rightWord ^= ( ( _soap3DnaComplement[leftChar] ) ^ rightChar ) << ( j % CHAR_PER_WORD * BIT_PER_CHAR );
            }

            // write back
            if ( leftWordIndex == rightWordIndex )
            {
                uint numLeftBits = ( ( readLength - 1 ) / 2 % CHAR_PER_WORD + 1 ) * BIT_PER_CHAR;
                uint numRightBits = 32 - numLeftBits;
                query[leftWordIndex * 32] = ( ( leftWord << numRightBits ) >> numRightBits ) |
                                            ( ( rightWord >> numLeftBits ) << numLeftBits );
            }
            else
            {
                query[leftWordIndex * 32] = leftWord;
                query[rightWordIndex * 32] = rightWord;
            }

            strand = 1 - strand;

            // match the negative strand
            if ( whichCase == 0 )
            {
#ifndef BGS_DISABLE_CASE_A
                matchQueryCaseA_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 1 )
            {
#ifndef BGS_DISABLE_CASE_B
                matchQueryCaseB_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 2 )
            {
#ifndef BGS_DISABLE_CASE_C
                matchQueryCaseC_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 3 )
            {
#ifndef BGS_DISABLE_CASE_D
                matchQueryCaseD_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 4 )
            {
#ifndef BGS_DISABLE_CASE_E
                matchQueryCaseE_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }

#endif

            // write error code if there are too many SA ranges
            if ( numSARanges == 0 )
            {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
                printf ( "[KERNEL] This query does not return result\n" );
#endif
                * ( answer ) = 0xFFFFFFFD;
            }
            else if ( numSARanges > sa_range_allowed )
            {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
                printf ( "[KERNEL] This query returns too many result\n" );
#endif
                * ( answer ) = 0xFFFFFFFE;

                if ( round == 0 )
                { isBad[queryId] = 1; }
            }
        }
        else
        {
            numSARanges = sa_range_allowed + 1;
            * ( answer ) = 0xFFFFFFFE;
        }

#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", answer );
        printf ( "[KERNEL] Sending out %u SA ranges..\n", numSARanges );

        for ( uint i = 0; i < wordPerQuery / 2; ++i )
        {
            printf ( "%d: %u %u\n", i, * ( answer + ( i * 2 ) * 32 ), * ( answer + ( i * 2 + 1 ) * 32 ) );
        }

#endif
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[KERNEL] Exiting successfully.\n" );
#endif
    }
}


// entry point of kernel
__global__ void kernel_4mismatch_2 ( uint whichCase, uint * queries, uint * readLengths, uint numQueries,
                                     uint wordPerQuery,
                                     uint * bwt, uint * occ, uint inverseSa0,
                                     uint * revBwt, uint * revOcc, uint revInverseSa0,
                                     uint textLength,
                                     uint * answers,
                                     bool * isBad, uint round, uint sa_range_allowed,
                                     uint wordPerAnswer, bool isExactNumMismatch )
{
    uint queryId = ( blockIdx.x * THREADS_PER_BLOCK + threadIdx.x ) *
                   QUERIES_PER_THREAD;

    if ( queryId < numQueries )
    {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[KERNEL] Starting..\n" );
#endif
        ///////////////// Process read //////////////////
        ullint idx = queryId % 32;
        ullint queryOffset = queryId / 32 * 32 * wordPerQuery + idx;
        uint * query = queries + queryOffset;
        ullint answerOffset = queryId / 32 * 32 * wordPerAnswer + idx;
        uint * answer = answers + answerOffset;
        uint readLength = readLengths[queryId];

        for ( uint i = 0; i < wordPerAnswer; ++i )
        { answer[i * 32] = 0xFFFFFFFF; }

        uint numSARanges = 0;
#ifndef BGS_DISABLE_NEGATIVE_STRAND
        uint strand = ( round > 0 ) ? 0 : ( whichCase % 2 );
#else
        uint strand = 0;
#endif

        if ( ( round > 0 ) || ( isBad[queryId] == 0 ) )
        {
            // match the positive strand
            if ( whichCase == 5 )
            {
#ifndef BGS_DISABLE_CASE_F
                matchQueryCaseF_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 6 )
            {
#ifndef BGS_DISABLE_CASE_G
                matchQueryCaseG_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 7 )
            {
#ifndef BGS_DISABLE_CASE_H
                matchQueryCaseH_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 8 )
            {
#ifndef BGS_DISABLE_CASE_I
                matchQueryCaseI_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 9 )
            {
#ifndef BGS_DISABLE_CASE_J
                matchQueryCaseJ_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }

#ifndef BGS_DISABLE_NEGATIVE_STRAND
            // reverse the read
            uint leftWordIndex = 0, rightWordIndex = ( readLength - 1 ) / CHAR_PER_WORD;
            uint leftWord = query[leftWordIndex * 32];
            uint rightWord = query[rightWordIndex * 32];

            for ( uint i = 0, j = readLength - 1; i <= j; ++i, --j )
            {
                // check if need to move to next word
                if ( i / CHAR_PER_WORD != leftWordIndex )
                {
                    // write back leftword
                    query[leftWordIndex * 32] = leftWord;
                    // load next leftword
                    leftWordIndex++;
                    leftWord = query[leftWordIndex * 32];
                }

                if ( j / CHAR_PER_WORD != rightWordIndex )
                {
                    // write back rightword
                    query[rightWordIndex * 32] = rightWord;
                    // load next rightword
                    rightWordIndex--;
                    rightWord = query[rightWordIndex * 32];
                }

                // swap left and right characters
                unsigned char leftChar = ( leftWord >> ( i % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
                unsigned char rightChar = ( rightWord >> ( j % CHAR_PER_WORD * BIT_PER_CHAR ) ) & CHAR_MASK;
                leftWord ^= ( leftChar ^ ( _soap3DnaComplement[rightChar] ) ) << ( i % CHAR_PER_WORD * BIT_PER_CHAR );
                rightWord ^= ( ( _soap3DnaComplement[leftChar] ) ^ rightChar ) << ( j % CHAR_PER_WORD * BIT_PER_CHAR );
            }

            // write back
            if ( leftWordIndex == rightWordIndex )
            {
                uint numLeftBits = ( ( readLength - 1 ) / 2 % CHAR_PER_WORD + 1 ) * BIT_PER_CHAR;
                uint numRightBits = 32 - numLeftBits;
                query[leftWordIndex * 32] = ( ( leftWord << numRightBits ) >> numRightBits ) |
                                            ( ( rightWord >> numLeftBits ) << numLeftBits );
            }
            else
            {
                query[leftWordIndex * 32] = leftWord;
                query[rightWordIndex * 32] = rightWord;
            }

            strand = 1 - strand;

            // match the negative strand
            if ( whichCase == 5 )
            {
#ifndef BGS_DISABLE_CASE_F
                matchQueryCaseF_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 6 )
            {
#ifndef BGS_DISABLE_CASE_G
                matchQueryCaseG_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 7 )
            {
#ifndef BGS_DISABLE_CASE_H
                matchQueryCaseH_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 8 )
            {
#ifndef BGS_DISABLE_CASE_I
                matchQueryCaseI_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }
            else if ( whichCase == 9 )
            {
#ifndef BGS_DISABLE_CASE_J
                matchQueryCaseJ_4mismatch ( query, readLength,
                                            bwt, occ, inverseSa0,
                                            revBwt, revOcc, revInverseSa0,
                                            numSARanges, answer, textLength,
                                            sa_range_allowed, strand );
#endif
            }

#endif

            // write error code if there are too many SA ranges
            if ( numSARanges == 0 )
            {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
                printf ( "[KERNEL] This query does not return result\n" );
#endif
                * ( answer ) = 0xFFFFFFFD;
            }
            else if ( numSARanges > sa_range_allowed )
            {
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
                printf ( "[KERNEL] This query returns too many result\n" );
#endif
                * ( answer ) = 0xFFFFFFFE;

                if ( round == 0 )
                { isBad[queryId] = 1; }
            }
        }
        else
        {
            numSARanges = sa_range_allowed + 1;
            * ( answer ) = 0xFFFFFFFE;
        }

#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[Global] Address %u\n", answer );
        printf ( "[KERNEL] Sending out %u SA ranges..\n", numSARanges );

        for ( uint i = 0; i < wordPerQuery / 2; ++i )
        {
            printf ( "%d: %u %u\n", i, * ( answer + ( i * 2 ) * 32 ), * ( answer + ( i * 2 + 1 ) * 32 ) );
        }

#endif
#ifdef BGS_OUTPUT_KERNEL_DEBUG_MESSAGE
        printf ( "[KERNEL] Exiting successfully.\n" );
#endif
    }
}
