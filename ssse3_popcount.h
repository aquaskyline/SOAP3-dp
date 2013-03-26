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

///////////////////////////////////////////////////////////
// Structures for storing single-end alignment result    //
///////////////////////////////////////////////////////////

#ifndef __SSSE3_POPCOUNT_H__
#define __SSSE3_POPCOUNT_H__


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <time.h>


uint32_t ssse3_popcount1 ( uint8_t * buffer, int chunks16 );
uint32_t ssse3_popcount2 ( uint8_t * buffer, int chunks16 );
uint32_t ssse3_popcount3 ( uint8_t * buffer, int chunks16 );
#endif
