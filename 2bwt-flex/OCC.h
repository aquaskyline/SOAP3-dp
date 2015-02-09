/*
 *
 *    OCC.h
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

#ifndef __OCC_H__
#define __OCC_H__

#include "../samtools-0.1.18/sam.h"

//Define the below parameter to output the alignment result(text position)
// with the old 2BWT-Aligner format
//#define DEBUG_2BWT_OUTPUT_32_BIT

#define OCC_CACHE_SIZE               81920
#define OCC_OUTPUT_FORMAT            20110320

typedef struct OCCPositionCache
{
	unsigned short ChromId;
	unsigned char ReadStrand;
	unsigned long long tp;
	int occMismatch;
	int resultSource;   // 0: SOAP3; 1: Semi-global DP
	char * cigarString; // pattern of the alignment
	unsigned int len;   // length of the read
} OCCPositionCache;

#ifdef DEBUG_2BWT_OUTPUT_32_BIT
typedef struct OCCPositionCacheToDisk
{
	unsigned char cell[5];
} OCCPositionCacheToDisk;
#endif

#ifndef DEBUG_2BWT_OUTPUT_32_BIT
typedef struct OCCPositionCacheToDisk
{
	unsigned char cell[11];
} OCCPositionCacheToDisk;
#endif

typedef struct OCC
{
	unsigned int occPositionCacheCount;
	OCCPositionCache occPositionCache[OCC_CACHE_SIZE];
	OCCPositionCacheToDisk occPositionCacheToDisk[OCC_CACHE_SIZE];

	bam1_t SAMOutBuffer;
} OCC;

#endif
