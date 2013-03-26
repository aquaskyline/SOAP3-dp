//
//    LTConstruct.c
//
//    SOAP2 / 2BWT
//
//    Copyright (C) 2012, HKU
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 2
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
///////////////////////////////////////////////////////////////////////////////////////////////
// vim:set autoindent shiftwidth=4 tabstop=4 noexpandtab:
#include "LTConstruct.h"

int LEN;
unsigned long long NR_TOP;
unsigned long long ALL_MASK;

const unsigned IN_BUF_SIZE = 100u Mibi;
const unsigned step = 1 Mibi;

unsigned char * ibuf;

unsigned long long text_pos; // Text position, 0-based
unsigned long long window; // Last length-LEN characters

int fin_src, fout_int;

unsigned int * otop;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline void process(int v) {
    ++text_pos;
    window = (window << LOOKUP_BIT_PER_CHAR) | v;
    if (text_pos >= -((unsigned long long)LEN)) return; // First LEN-1 characters
    unsigned long long top = (window & ALL_MASK);
    ++otop[top];
}

inline void gen() {
    // init
    int charMask = (1 << LOOKUP_BIT_PER_CHAR) - 1;
    text_pos = -(unsigned long long)LEN;
    window = 0;
    // start
    TRY(lseek64(fin_src, 0, SEEK_SET));
    while (1) {
        unsigned int i;
        unsigned nbuf = read(fin_src, ibuf, sizeof(*ibuf) * IN_BUF_SIZE);
        if (nbuf < IN_BUF_SIZE) nbuf -= 2;
        for (i = 0; i < nbuf; ++i) {
            unsigned short c = ibuf[i];
            int j;
            for (j = 8-BIT_PER_CHAR; j >= 0; j -= BIT_PER_CHAR) process((c >> j) & charMask);
        }
        if (nbuf < IN_BUF_SIZE) {
            int rem = ibuf[nbuf + 1];
            unsigned short c = ibuf[nbuf];
            int j;
            for (j = 8-BIT_PER_CHAR; j >= 8 - (rem * BIT_PER_CHAR); j -= BIT_PER_CHAR) process((c >> j) & charMask);
            break;
        }
    }
}

int BuildLookupTable(char * PackedDNAFileName, char * LookupTableFileName, int lookupTableSize) {
    
    unsigned long long i;

    LEN = lookupTableSize;
    NR_TOP = 1LL << LEN * LOOKUP_BIT_PER_CHAR;
    ALL_MASK = NR_TOP - 1;

    ibuf = (unsigned char*) malloc(sizeof(unsigned char)*IN_BUF_SIZE);
    otop = (unsigned int*) malloc(sizeof(unsigned int)*NR_TOP);

    TRY(fin_src = open64(PackedDNAFileName, O_RDONLY ));
    TRY(fout_int = open64(LookupTableFileName, O_CREAT | O_WRONLY , 0664));
    TRYEQ(write(fout_int, &lookupTableSize, sizeof(int)), sizeof(int));
    gen();
    for (i = 1; i < LEN; ++i) process(0);
    for (i = 1; i < NR_TOP; ++i) otop[i] += otop[i-1];
    for (i = 0; i < NR_TOP; i += step) {
        TRYEQ(write(fout_int, otop + i, step * sizeof(*otop)), step * sizeof(*otop));
    }
    TRY(close(fin_src));
    TRY(close(fout_int));

    free(ibuf);
    free(otop);
}

