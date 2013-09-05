/*
 *
 *    QueryParser.h
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

#ifndef __QUERY_PARSER_H__
#define __QUERY_PARSER_H__

#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <zlib.h>
#include "SAM.h"
#include "definitions.h"
#include "samtools-0.1.18/bam.h"

// An updated function is to load the pair-end reads for at most "maxNumQueries/2" # of pairs of reads
// faster than the previous version
int loadPairReadsGz2 ( gzFile queryFile, char * queryFileBuffer, gzFile queryFile2, char * queryFileBuffer2, unsigned char * charMap,
                       uint * queries, uint * readLengths, uint * readIDs, char * upkdQualities,  char * upkdQueryNames,
                       uint maxReadLength, uint maxNumQueries,
                       size_t & bufferSize, char & queryChar, uint & bufferIndex, size_t & bufferSize2, char & queryChar2, uint & bufferIndex2,
                       uint accumReadNum, uint wordPerQuery, int qualityConstant, char & isFastq, int maxLenReadName );


// This function is to load the single reads for at most "maxNumQueries" # of reads
int loadSingleReadsGz ( gzFile queryFile, char * queryFileBuffer, unsigned char * charMap,
                        uint * queries, uint * readLengths, uint * readIDs,
                        char * upkdQualities, char * upkdQueryNames, uint maxReadLength, uint maxNumQueries,
                        size_t & bufferSize, char & queryChar, uint & bufferIndex, uint accumReadNum, uint wordPerQuery,
                        int qualityConstant, char & isFastq, int maxLenReadName );

// This function is to load reads for BAM format
int loadBAMReads ( bamFile bamQueryFile, bam_header_t * bamHeader, bam1_t * b, unsigned char * charMap,
                   uint * queries, uint * readLengths, uint * readIDs, char * upkdQualities, char * upkdQueryNames, uint maxReadLength, uint maxNumQueries,
                   uint wordPerQuery, int qualityConstant, int maxLenReadName );

// scan the first ten reads
// and get the max read length among them
uint GetReadLength ( uint * readLengths, uint numQueries, int sample );

#endif

