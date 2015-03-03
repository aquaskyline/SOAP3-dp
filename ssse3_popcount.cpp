/*
 *
 *    ssse3_popcount.h
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
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <time.h>

#ifdef __x86_64
//#define ALIGN_DATA
#ifdef ALIGN_DATA
#   define __aligned__ __attribute__((aligned(16)))
#else
#   define __aligned__
#endif

#define MAX_CHUNKS 2048
// lookup for SSE
uint8_t POPCOUNT_4bit[16] __aligned__ =
{
    /* 0 */ 0,
    /* 1 */ 1,
    /* 2 */ 1,
    /* 3 */ 2,
    /* 4 */ 1,
    /* 5 */ 2,
    /* 6 */ 2,
    /* 7 */ 3,
    /* 8 */ 1,
    /* 9 */ 2,
    /* a */ 2,
    /* b */ 3,
    /* c */ 2,
    /* d */ 3,
    /* e */ 3,
    /* f */ 4
};

// ---- SSSE3 - naive approach --------------------------------------------
uint32_t ssse3_popcount1 ( uint8_t * buffer, int chunks16 )
{
    uint32_t dummy __attribute__ ( ( unused ) );
    static char MASK_4bit[16] = {0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf};
    uint32_t result, tmp;
    int n;
    __asm__ volatile ( "movdqu (%%eax), %%xmm7" : : "a" ( POPCOUNT_4bit ) );
    __asm__ volatile ( "movdqu (%%eax), %%xmm6" : : "a" ( MASK_4bit ) );
    result = 0;

    for ( n = 0; n < chunks16; n++ )
    {
        __asm__ volatile (
#ifdef ALIGN_DATA
            "movdqa	  (%%ebx), %%xmm0	\n"
#else
            "movdqu	  (%%ebx), %%xmm0	\n"
#endif
            "movdqa    %%xmm0, %%xmm1	\n"
            "psrlw         $4, %%xmm1	\n"
            "pand      %%xmm6, %%xmm0	\n"   // xmm0 := lower nibbles
            "pand      %%xmm6, %%xmm1	\n"   // xmm1 := higher nibbles
            "movdqa    %%xmm7, %%xmm2	\n"
            "movdqa    %%xmm7, %%xmm3	\n"   // get popcount
            "pshufb    %%xmm0, %%xmm2	\n"   // for all nibbles
            "pshufb    %%xmm1, %%xmm3	\n"   // using PSHUFB
            "paddb     %%xmm2, %%xmm3	\n"   // popcount for all bytes
            "pxor      %%xmm0, %%xmm0	\n"
            "psadbw    %%xmm3, %%xmm0	\n"   // sum popcounts
            "movhlps   %%xmm0, %%xmm1	\n"
            "paddd     %%xmm0, %%xmm1	\n"
            "movd      %%xmm1, %%eax	\n"
            : "=a" ( tmp )
            : "b" ( &buffer[n*16] )
        );
        result += tmp;
    }

    return result;
}


// ---- SSSE3 - better alorithm, minimized psadbw usage -------------------

uint32_t ssse3_popcount2 ( uint8_t * buffer, int chunks16 )
{
    static char MASK_4bit[16] = {0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf};
    uint32_t result = 0;
    __asm__ volatile ( "movdqu (%%eax), %%xmm7" : : "a" ( POPCOUNT_4bit ) );
    __asm__ volatile ( "movdqu (%%eax), %%xmm6" : : "a" ( MASK_4bit ) );
    __asm__ volatile ( "pxor    %%xmm5, %%xmm5" : : ); // xmm5 -- global accumulator
    int n, i = 0;
    __asm__ volatile ( "pxor %xmm4, %xmm4" ); // xmm4 -- local accumulator

    for ( n = 0; n < chunks16; n++ )
    {
        __asm__ volatile (
#if ALIGN_DATA
            "movdqa   (%%eax), %%xmm0       \n"
#else
            "movdqu   (%%eax), %%xmm0       \n"
#endif
            "movdqa    %%xmm0, %%xmm1       \n"
            "psrlw         $4, %%xmm1       \n"
            "pand      %%xmm6, %%xmm0       \n"     // xmm0 := lower nibbles
            "pand      %%xmm6, %%xmm1       \n"     // xmm1 := higher nibbles
            "movdqa    %%xmm7, %%xmm2       \n"
            "movdqa    %%xmm7, %%xmm3       \n"     // get popcount
            "pshufb    %%xmm0, %%xmm2       \n"     // for all nibbles
            "pshufb    %%xmm1, %%xmm3       \n"     // using PSHUFB
            "paddb     %%xmm2, %%xmm4       \n"     // update local
            "paddb     %%xmm3, %%xmm4       \n"     // accumulator
            :
            : "a" ( &buffer[i] )
        );
        i += 16;
    }

    // update global accumulator (two 32-bits counters)
    __asm__ volatile (
        "pxor   %xmm0, %xmm0            \n"
        "psadbw %xmm0, %xmm4            \n"
        "paddd  %xmm4, %xmm5            \n"
    );
    // finally add together 32-bits counters stored in global accumulator
    __asm__ volatile (
        "movhlps   %%xmm5, %%xmm0       \n"
        "paddd     %%xmm5, %%xmm0       \n"
        "movd      %%xmm0, %%eax        \n"
        : "=a" ( result )
    );
    return result;
}

// ---- SSSE3 - better alorithm, inner loop unrolled ----------------------
uint32_t ssse3_popcount3 ( uint8_t * buffer, int chunks16 )
{
    static char MASK_4bit[16] = {0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf};
    uint32_t result = 0;
    __asm__ volatile ( "movdqu (%%eax), %%xmm7" : : "a" ( POPCOUNT_4bit ) );
    __asm__ volatile ( "movdqu (%%eax), %%xmm6" : : "a" ( MASK_4bit ) );
    __asm__ volatile ( "pxor    %%xmm5, %%xmm5" : : ); // xmm5 -- global accumulator
    int k, n, i;
    i = 0;

    while ( chunks16 > 0 )
    {
        // max(POPCOUNT_8bit) = 8, thus byte-wise addition could be done
        // for floor(255/8) = 31 iterations
#define MAX (7*4)
        if ( chunks16 > MAX )
        {
            k = MAX;
            chunks16 -= MAX;
        }
        else
        {
            k = chunks16;
            chunks16 = 0;
        }

#undef MAX
        __asm__ volatile ( "pxor %xmm4, %xmm4" ); // xmm4 -- local accumulator
        const int c = 4;

        for ( n = 0; n < k; n += c )
        {
#define body(index) \
    __asm__ volatile( \
                      "movdqa	  (%%eax), %%xmm0	\n" \
                      "movdqa    %%xmm0, %%xmm1	\n" \
                      "psrlw         $4, %%xmm1	\n" \
                      "pand      %%xmm6, %%xmm0	\n" \
                      "pand      %%xmm6, %%xmm1	\n" \
                      "movdqa    %%xmm7, %%xmm2	\n" \
                      "movdqa    %%xmm7, %%xmm3	\n" \
                      "pshufb    %%xmm0, %%xmm2	\n" \
                      "pshufb    %%xmm1, %%xmm3	\n" \
                      "paddb     %%xmm2, %%xmm4	\n" \
                      "paddb     %%xmm3, %%xmm4	\n" \
                      : : "a" (&buffer[index]));
            body ( i );
            body ( i + 1 * 16 );
            body ( i + 2 * 16 );
            body ( i + 3 * 16 );
#undef body
            //          i += 4*16;
            i += c * 16;
        }

        // update global accumulator (two 32-bits counters)
        __asm__ volatile (
            "pxor	%xmm0, %xmm0		\n"
            "psadbw	%xmm0, %xmm4		\n"
            "paddd	%xmm4, %xmm5		\n"
        );
    }

    // finally add together 32-bits counters stored in global accumulator
    __asm__ volatile (
        "movhlps   %%xmm5, %%xmm0	\n"
        "paddd     %%xmm5, %%xmm0	\n"
        "movd      %%xmm0, %%eax	\n"
        : "=a" ( result )
    );
    return result;
}
#endif
