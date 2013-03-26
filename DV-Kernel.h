/*
 *
 *    DV-Kernel.h
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

#ifndef _GPUFUNCTIONS_H_
#define _GPUFUNCTIONS_H_

#include "definitions.h"

// The offset mask to retrieve the least significant 24-bit from the 32-bit word.
#define  BGS_GPU_ANSWER_OFFSET_LENGTH   24
#define  BGS_GPU_ANSWER_OFFSET_MASK     ((1<<BGS_GPU_ANSWER_OFFSET_LENGTH)-1)

//Constant Memory Declaration
__constant__ __device__ char gpuCharMap[256];
__constant__ __device__ char _soap3DnaComplement[ALPHABET_SIZE]        = { 3, 2, 1, 0 };

__global__ void kernel ( uint whichCase, uint * queries, uint * readLengths, uint numQueries,
                         uint wordPerQuery,
                         uint * bwt, uint * occ, uint inverseSa0,
                         uint * revBwt, uint * revOcc, uint revInverseSa0,
                         uint textLength,
                         uint * answers,
                         bool * isBad, uint round, uint numMismatch,
                         uint sa_range_allowed, uint wordPerAnswer,
                         bool isExactNumMismatch );


__global__ void kernel_4mismatch_1 ( uint whichCase, uint * queries, uint * readLengths, uint numQueries,
                                     uint wordPerQuery,
                                     uint * bwt, uint * occ, uint inverseSa0,
                                     uint * revBwt, uint * revOcc, uint revInverseSa0,
                                     uint textLength,
                                     uint * answers,
                                     bool * isBad, uint round, uint sa_range_allowed,
                                     uint wordPerAnswer, bool isExactNumMismatch );


__global__ void kernel_4mismatch_2 ( uint whichCase, uint * queries, uint * readLengths, uint numQueries,
                                     uint wordPerQuery,
                                     uint * bwt, uint * occ, uint inverseSa0,
                                     uint * revBwt, uint * revOcc, uint revInverseSa0,
                                     uint textLength,
                                     uint * answers,
                                     bool * isBad, uint round, uint sa_range_allowed,
                                     uint wordPerAnswerbool, bool isExactNumMismatch );

#endif
