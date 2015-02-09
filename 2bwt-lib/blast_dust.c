/*

   blast_dust.c        Adaptation of dust filtering in BLAST

*/

/* $Id: blast_dust.c,v 1.37 2005/11/16 14:27:03 madden Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ==========================================================================
 *
 * Authors: Richa Agarwala (based upon versions variously worked upon by Roma 
 *          Tatusov, John Kuzio, and Ilya Dondoshansky).
 *   
 * ==========================================================================
 */

/** @file blast_dust.c
 * A utility to find low complexity NA regions. This parallels functionality 
 * of dust.c from the C toolkit, but without using the structures generated 
 * from ASN.1 spec.
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: blast_dust.c,v 1.37 2005/11/16 14:27:03 madden Exp $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blast_dust.h"


const int kDustLevel = 20;
const int kDustWindow = 64;
const int kDustLinker = 1;


/* local, file scope, structures and variables */

/** localcurrents 
 * @todo expand documentation
 */
typedef struct DCURLOC { 
    int    cursum, curstart, curend;
    short    curlength;
} DCURLOC;

unsigned char *SUM_THRESHOLD;

/* local functions */

static void wo (int, unsigned char*, int, DCURLOC*, unsigned char*, unsigned char, int, const unsigned char*);
static unsigned char wo1 (int, unsigned char*, int, DCURLOC*);
static int dust_triplet_find (unsigned char*, int, int, unsigned char*, const unsigned char*);

/** Safe free a pointer: belongs to a higher level header. */
#ifndef sfree
#define sfree(x) __sfree((void**)&(x))
#endif


void __sfree(void **x)
{
  free(*x);
  *x=NULL;
  return;
}



/* entry point for dusting */

int blast_dust (unsigned char* __restrict sequence, const int length, int level, int window, int linker)
{
   int    len;
   int    i;
   unsigned char* seq;
   DCURLOC    cloc;

   int to, from;
   int j;
   unsigned char lowercase[256];
   unsigned char blastnaMap[256];

   // initialize blastnaMap
   for (i=0; i<256; i++) {
      blastnaMap[i] = 15;
   }
   blastnaMap['A'] = blastnaMap['a'] = 0;
   blastnaMap['C'] = blastnaMap['c'] = 1;
   blastnaMap['M'] = blastnaMap['m'] = 6;
   blastnaMap['G'] = blastnaMap['g'] = 2;
   blastnaMap['R'] = blastnaMap['r'] = 4;
   blastnaMap['S'] = blastnaMap['s'] = 9;
   blastnaMap['V'] = blastnaMap['v'] = 13;
   blastnaMap['T'] = blastnaMap['t'] = 3;
   blastnaMap['W'] = blastnaMap['w'] = 8;
   blastnaMap['Y'] = blastnaMap['y'] = 5;
   blastnaMap['H'] = blastnaMap['h'] = 12;
   blastnaMap['K'] = blastnaMap['k'] = 7;
   blastnaMap['D'] = blastnaMap['d'] = 11;
   blastnaMap['B'] = blastnaMap['b'] = 10;
   blastnaMap['N'] = blastnaMap['n'] = 14;

   for (i=0; i<256; i++) {
      lowercase[i] = 'l';
   }
   for (i=0; i<26; i++) {
      lowercase['a'+i] = (unsigned char)('a'+i);
      lowercase['A'+i] = (unsigned char)('a'+i);
   }
   
   to = 0 - linker - 1;
   from = 0;

   /* defaults are more-or-less in keeping with original dust */
   if (level < 2 || level > 64) level = kDustLevel;
   if (window < 8 || window > 64) window = kDustWindow;
   if (linker < 1 || linker > 32) linker = kDustLinker;
   
   seq = (unsigned char*) calloc(1, window);            /* triplets */
   if (!seq) {
      return -1;
   }
   SUM_THRESHOLD = (unsigned char *) calloc(window, sizeof(unsigned char));  
   if (!SUM_THRESHOLD) {
      return -1;
   }
   SUM_THRESHOLD[0] = 1;
   for (i=1; i < window; i++)
       SUM_THRESHOLD[i] = (unsigned char)((level * i)/10);

   if (length < window) window = length;

   /* Consider smaller windows in beginning of the sequence */
   for (i = 2; i <= window-1; i++) {
      len = i-1;
      wo (len, sequence, 0, &cloc, seq, 1, level, blastnaMap);
      
      if (cloc.cursum*10 > level*cloc.curlength) {
         if (to + linker >= cloc.curstart &&
             from <= cloc.curend + linker) {
            /* overlap windows nicely if needed */
            if (to < cloc.curend)
                to = cloc.curend;
            if (from > cloc.curstart)
                from = cloc.curstart;
         } else    {
            /* new window or dusted regions do not overlap */
            for (j=from; j<=to; j++) {
               sequence[j] = lowercase[sequence[j]];
            }
            from = cloc.curstart;
            to = cloc.curend;
         }
      }                /* end 'if' high score    */
   }                    /* end for */

   for (i = 1; i < length-2; i++) {
      len = (int) ((length > i+window) ? window : length-i);
      len -= 2;
      if (length >= i+window)
          wo (len, sequence, i, &cloc, seq, 2, level, blastnaMap);
      else /* remaining portion of sequence is less than windowsize */
          wo (len, sequence, i, &cloc, seq, 3, level, blastnaMap);
      
      if (cloc.cursum*10 > level*cloc.curlength) {
         if (to + linker >= cloc.curstart+i &&
             from <= cloc.curend + i + linker) {
            /* overlap windows nicely if needed */
            if (to < cloc.curend + i)
                to = cloc.curend + i;
            if (from > cloc.curstart + i)
                from = cloc.curstart + i;
         } else    {
            /* new window or dusted regions do not overlap */
            for (j=from; j<=to; j++) {
               sequence[j] = lowercase[sequence[j]];
            }
            from = cloc.curstart + i;
            to = cloc.curend + i;
         }
      }                /* end 'if' high score    */
   }                    /* end for */
   for (j=from; j<=to; j++) {
      sequence[j] = lowercase[sequence[j]];
   }
   sfree (seq);
   sfree(SUM_THRESHOLD);
   return 0;
}

static void wo (int len, unsigned char* seq_start, int iseg, DCURLOC* cloc, 
                unsigned char* seq, unsigned char FIND_TRIPLET, int level, const unsigned char* blastnaMap)
{
    int smaller_window_start, mask_window_end;
        unsigned char SINGLE_TRIPLET;

    cloc->cursum = 0;
    cloc->curlength = 1;
    cloc->curstart = 0;
    cloc->curend = 0;

    if (len < 1)
        return;

        /* get the chunk of sequence in triplets */
    if (FIND_TRIPLET==1) /* Append one */
    {
        seq[len-1] = seq[len] = seq[len+1] = 0;
        dust_triplet_find (seq_start, iseg+len-1, 1, seq+len-1, blastnaMap);
    }
    if (FIND_TRIPLET==2) /* Copy suffix as prefix and find one */
    {
        memmove(seq,seq+1,(len-1)*sizeof(unsigned char));
        seq[len-1] = seq[len] = seq[len+1] = 0;
        dust_triplet_find (seq_start, iseg+len-1, 1, seq+len-1, blastnaMap);
    }
    if (FIND_TRIPLET==3) /* Copy suffix */
        memmove(seq,seq+1,len*sizeof(unsigned char));

        /* dust the chunk */
    SINGLE_TRIPLET = wo1 (len, seq, 0, cloc); /* dust at start of window */

        /* consider smaller windows only if anything interesting 
           found for starting position  and smaller windows have potential of
           being at higher level */
        if ((cloc->cursum*10 > level*cloc->curlength) && (!SINGLE_TRIPLET)) {
        mask_window_end = cloc->curend-1;
        smaller_window_start = 1;
                while ((smaller_window_start < mask_window_end) &&
                       (!SINGLE_TRIPLET)) {
            SINGLE_TRIPLET = wo1(mask_window_end-smaller_window_start,
                             seq+smaller_window_start, smaller_window_start, cloc);
                    smaller_window_start++;
            }
    }

    cloc->curend += cloc->curstart;
}

/** returns TRUE if there is single triplet in the sequence considered */
static unsigned char wo1 (int len, unsigned char* seq, int iwo, DCURLOC* cloc)
{
   unsigned int sum;
    int loop;

    short* countsptr;
    short counts[4*4*4];
    unsigned char triplet_count = 0;

    memset (counts, 0, sizeof (counts));
/* zero everything */
    sum = 0;

/* dust loop -- specific for triplets    */
    for (loop = 0; loop < len; loop++)
    {
        countsptr = &counts[*seq++];
        if (*countsptr)
        {
            sum += (unsigned int)(*countsptr);

            if (sum >= SUM_THRESHOLD[loop])
            {
            if ((unsigned int)cloc->cursum*loop < sum*cloc->curlength)
            {
                cloc->cursum = sum;
                cloc->curlength = (short)loop;
                cloc->curstart = iwo;
                cloc->curend = loop + 2; /* triplets */
            }
            }
        }
        else
            triplet_count++;
        (*countsptr)++;
    }

    if (triplet_count > 1)
        return(FALSE);
    return(TRUE);
}

/** Fill an array with 2-bit encoded triplets.
 * @param seq_start Pointer to the start of the sequence in blastna 
 *                  encoding [in]
 * @param icur Offset at which to start extracting triplets [in]
 * @param max Maximal length of the sequence segment to be processed [in]
 * @param s1 Array of triplets [out]
 * @return How far was the sequence processed?
 */
static int 
dust_triplet_find (unsigned char* seq_start, int icur, int max, unsigned char* s1, const unsigned char* blastnaMap)
{
   int n;
   unsigned char* s2,* s3;
   short c;
   unsigned char* seq = &seq_start[icur];
   unsigned char end_byte = 15;    // 15 is not used in DNA sequence

   n = 0;
   
   s2 = s1 + 1;
   s3 = s1 + 2;
   
   /* set up 1 */
   if ((c = blastnaMap[*seq++]) == end_byte)
      return n;
   c &= NCBI2NA_MASK;
   *s1 |= c;
   *s1 <<= 2;
   
   /* set up 2 */
   if ((c = blastnaMap[*seq++]) == end_byte)
      return n;
   c &= NCBI2NA_MASK;
   *s1 |= c;
   *s2 |= c;
   
   /* triplet fill loop */
   while (n < max && (c = blastnaMap[*seq++]) != end_byte) {
      c &= NCBI2NA_MASK;
      *s1 <<= 2;
      *s2 <<= 2;
      *s1 |= c;
      *s2 |= c;
      *s3 |= c;
      s1++;
      s2++;
      s3++;
      n++;
   }
   
   return n;
}

