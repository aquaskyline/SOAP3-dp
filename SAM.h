/*
 *
 *    SAM.h
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

#ifndef __SAM_H__
#define __SAM_H__

#include <stdio.h>
#include <stdlib.h>
#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "2bwt-flex/HOCC.h"
#include "2bwt-flex/OCC.h"
#include "definitions.h"
#include "Release.h"

#include "samtools-0.1.18/sam.h"

#define SAM_FLAG_IS_PAIR_READ           1
#define SAM_FLAG_PROPER_PAIR            2
#define SAM_FLAG_READ_UNMAPPED          4
#define SAM_FLAG_MATE_UNMAPPED          8

#define SAM_FLAG_READ_ALGNMT_STRAND     16
#define SAM_FLAG_MATE_ALGNMT_STRAND     32
#define SAM_FLAG_FIRST_IN_PAIR          64
#define SAM_FLAG_SECOND_IN_PAIR         128

#define SAM_FLAG_SECONDARY_ALGNMT       256
#define SAM_FLAG_QC_CHECK_FAILED        512
#define SAM_FLAG_DUPLICATE              1024

#define SAM_MDATA_SIZE                  2048

// SAM_MAPQ_UNAVAILABLE should be set to 255
// according to SAM v1.4 specification.
// However we are making an exception here
// as some of the downstream programs are not
// capable of handling 255.
#define SAM_MAPQ_UNAVAILABLE            255

int SAMIUint8ConcatUint8 ( uint8_t * data, int * curSize,
                           uint8_t key );
int SAMIUint8ConcatUint32 ( uint8_t * data, int * curSize,
                            uint32_t key );
int SAMIUint8ConcatString ( uint8_t * data, int * curSize,
                            char * key, int len );

void SAMOutputHeaderConstruct ( bam_header_t * sheader, HSP * hsp, HSPAux * hspaux, int maxReadLength );
void SAMOutputHeaderDestruct ( bam_header_t * sheader );
void SAMOccurrenceConstruct ( OCC * occ );
void SAMOccurrenceDestruct ( OCC * occ );

#endif
