/*

   IndexCtl.h        Index Control

   This module provides helper function to maintain index header.

   Copyright (C) 2013, Edward Man Kit Wu.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef __INDEX_CONTROL_H__
#define __INDEX_CONTROL_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define INDEXCTL_MAX_HEADER_BLOCK 10
#define INDEXCTL_RESERVE_BYTE     256
#define INDEXCTL_BLOCK_BYTE       4096

// IndexHeaderReserve should be maintained to be smaller than INDEXCTL_RESERVE_BYTE-byte in size
typedef struct IndexHeaderReserve {
    uint32_t totalSize;       // Auto-populated
    char indexName[16];
    uint32_t versionMajor;    //4-byte
    uint32_t versionMinor;    //4-byte
    uint32_t versionRev;      //4-byte
} IndexHeaderReserve;

typedef struct IndexHeader {
    uint8_t buffer[INDEXCTL_MAX_HEADER_BLOCK * INDEXCTL_BLOCK_BYTE];
    uint32_t totalSize;
    uint32_t readSize;
} IndexHeader;

// IDXHEnqueueUint32 enqueue a uint32_t value to the back of the header block.
uint8_t IDXHEnqueueUint32(IndexHeader * idxHeader, uint32_t key32);

// IDXHEnqueueUint32Array enqueue an array of uint32_t value to the back of the header block.
uint8_t IDXHEnqueueUint32Array(IndexHeader * idxHeader, uint32_t * key32, uint32_t count);

// IDXHEnqueueRaw enqueue a byte-array to the back of the header block.
uint8_t IDXHEnqueueRaw(IndexHeader * idxHeader, uint8_t * key8, uint32_t count);

// IDXHDequeueUint32 dequeue a uint32_t value to the back of the header block.
// and stores it to key32.
uint8_t IDXHDequeueUint32(IndexHeader * idxHeader, uint32_t * key32);

// IDXHDequeueUint32Array dequeue an array of uint32_t value to the back of the header block.
// and stores it to key32.
uint8_t IDXHDequeueUint32Array(IndexHeader * idxHeader, uint32_t * key32, uint32_t count);

// IDXHDequeueRaw dequeue a byte-array to the back of the header block.
uint8_t IDXHDequeueRaw(IndexHeader * idxHeader, uint8_t * key8, uint32_t count);

// IDXHCreate creates a header structure and initialise it
IndexHeader * IDXHCreate();

// IDXHFree free the created IndexHeader structure
void  IDXHFree(IndexHeader * indexHeader);

// IDXHWriteReserved writes the reserve header structure onto the IndexHeader
void IDXHWriteReserved (IndexHeader * idxHeader, IndexHeaderReserve * reserv);

// IDXHWrite writes the header block into the file opened as fp.
void IDXHWrite(FILE * fp, IndexHeader * idxHeader, IndexHeaderReserve * reserv);

// IDXHRead reads the header data from dp and populate IndexHeaderReserve
IndexHeader * IDXHRead(FILE * fp, IndexHeaderReserve * reserv);

#endif
