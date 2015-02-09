/*

   IndexCtl.c        Index Control

   This module provides helper function to maintain index header.

   Copyright (C) 2004, Edward Man Kit Wu.

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

#include "IndexCtl.h"

uint8_t IDXHEnqueueUint32(IndexHeader * idxHeader, uint32_t key32) {
    IDXHEnqueueRaw(idxHeader,(uint8_t*)(&key32),sizeof(uint32_t)/sizeof(uint8_t));
}

uint8_t IDXHEnqueueUint32Array(IndexHeader * idxHeader, uint32_t * key32, uint32_t count) {
    IDXHEnqueueRaw(idxHeader,(uint8_t*)key32,count*sizeof(uint32_t)/sizeof(uint8_t));
}

uint8_t IDXHEnqueueRaw(IndexHeader * idxHeader, uint8_t * key8, uint32_t count) {

    uint32_t byteLeftOnBuffer = INDEXCTL_MAX_HEADER_BLOCK * INDEXCTL_BLOCK_BYTE - idxHeader->totalSize;
    uint8_t * bufferPtr = idxHeader->buffer + idxHeader->totalSize;
    
    if ( count > byteLeftOnBuffer ) {
        fprintf(stderr,"Header total size exceeds %u bytes\n", INDEXCTL_MAX_HEADER_BLOCK * INDEXCTL_BLOCK_BYTE);
        exit(1);
    }
    
    memcpy(bufferPtr,key8,count);
    idxHeader->totalSize += count;
}

uint8_t IDXHDequeueUint32(IndexHeader * idxHeader, uint32_t * key32) {
    IDXHDequeueRaw(idxHeader,(uint8_t*)key32,sizeof(uint32_t)/sizeof(uint8_t));
}

uint8_t IDXHDequeueUint32Array(IndexHeader * idxHeader, uint32_t * key32, uint32_t count) {
    IDXHDequeueRaw(idxHeader,(uint8_t*)key32,count*sizeof(uint32_t)/sizeof(uint8_t));
}

uint8_t IDXHDequeueRaw(IndexHeader * idxHeader, uint8_t * key8, uint32_t count) {

    uint32_t byteLeftOnBuffer = idxHeader->totalSize - idxHeader->readSize;
    uint8_t * bufferPtr = idxHeader->buffer + idxHeader->readSize;
    
    if ( count > byteLeftOnBuffer ) {
        fprintf(stderr,"Header total size exceeds %u bytes\n", INDEXCTL_MAX_HEADER_BLOCK * INDEXCTL_BLOCK_BYTE);
        exit(1);
    }
    
    memcpy(key8,bufferPtr,count);
    idxHeader->readSize += count;
}

IndexHeader * IDXHCreate() {
    IndexHeader * idxHeader = (IndexHeader *) malloc(sizeof(IndexHeader));
    memset(idxHeader->buffer,0,INDEXCTL_MAX_HEADER_BLOCK * INDEXCTL_RESERVE_BYTE);
    
    idxHeader->readSize = INDEXCTL_RESERVE_BYTE;
    idxHeader->totalSize = INDEXCTL_RESERVE_BYTE;
    
    return idxHeader;
}

void  IDXHFree(IndexHeader * idxHeader) {
    free(idxHeader);
}

void IDXHWriteReserved (IndexHeader * idxHeader, IndexHeaderReserve * reserv) {
    if (sizeof(IndexHeaderReserve)>INDEXCTL_RESERVE_BYTE) {
        fprintf(stderr,"IndexHeaderReserve structure is larger than INDEXCTL_RESERVE_BYTE\n");
        fprintf(stderr,"Caution: Changing INDEXCTL_RESERVE_BYTE harms backward compatibility of index\n");
        exit(1);
    }
    reserv->totalSize = idxHeader->totalSize;
    memcpy(idxHeader->buffer,reserv,sizeof(IndexHeaderReserve));
    memset(idxHeader->buffer+sizeof(IndexHeaderReserve),0,
        INDEXCTL_RESERVE_BYTE-sizeof(IndexHeaderReserve));
}

void IDXHWrite(FILE * fp, IndexHeader * idxHeader, IndexHeaderReserve * reserv) {
    int i;
    IDXHWriteReserved(idxHeader,reserv);
    uint8_t blockCount = (idxHeader->totalSize + INDEXCTL_BLOCK_BYTE - 1) / INDEXCTL_BLOCK_BYTE;
    uint8_t * bufferPtr = idxHeader->buffer;
    for (i=0;i<blockCount;i++) {
        fwrite(bufferPtr,sizeof(uint8_t),INDEXCTL_BLOCK_BYTE,fp);
        bufferPtr+=INDEXCTL_BLOCK_BYTE;
    }
}

IndexHeader * IDXHRead(FILE * fp, IndexHeaderReserve * reserv) {
    int i;
    IndexHeader * idxHeader = IDXHCreate();
    
    uint8_t * bufferPtr = idxHeader->buffer;
    fread(bufferPtr,sizeof(uint8_t),INDEXCTL_BLOCK_BYTE,fp);
    memcpy(reserv,bufferPtr,sizeof(IndexHeaderReserve));
    bufferPtr+=INDEXCTL_BLOCK_BYTE;
    
    idxHeader->totalSize = reserv->totalSize;
    uint8_t blockCount = (idxHeader->totalSize + INDEXCTL_BLOCK_BYTE - 1) / INDEXCTL_BLOCK_BYTE;
    uint32_t byteLeftToRead = idxHeader->totalSize;
    
    for (i=1;i<blockCount;i++) {
        fread(bufferPtr,sizeof(uint8_t),INDEXCTL_BLOCK_BYTE,fp);
        bufferPtr+=INDEXCTL_BLOCK_BYTE;
    }
    
    return idxHeader;
}

