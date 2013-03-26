//
//    HOCCConstruct.h
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

#ifndef __HOCC_CONSTRUCT_H__
#define __HOCC_CONSTRUCT_H__

#include <stdio.h>
#include <stdlib.h>

#include "../2bwt-lib/BWT.h"
#include "../2bwt-lib/HSP.h"

typedef struct HOCCCellSpace * ho_cellSpacePtr;
struct HOCCCellSpace {
    unsigned int count;
    unsigned int *saL;
    unsigned int *saR;
    struct HOCCCellSpace * next;
};

typedef struct HOCCBuilder {
    //hashtable essentials
    unsigned int ho_ttlItem;
    unsigned int ho_ttlOccurrence;
    unsigned int ho_Size;
    
    ho_cellSpacePtr ho_root;
    ho_cellSpacePtr ho_currentRoot;
} HOCCBuilder;


void HOCCExtractionHighOccPattern(HOCCBuilder * hoccBuilder, const BWT *bwt,const BWT * rev_bwt, 
                                int index, unsigned int HashOccThreshold,
                                unsigned int l, unsigned int r,unsigned int rev_l, unsigned int rev_r);

void HOCCPopulateHighOccPatternToBitPattern(HOCCBuilder * hoccBuilder, BWT * bwt, HSP* hsp,
                                const char * HighOccPatternFileName,const char * HighOccPackedFileName, 
                                const char * HighOccAuxValueFileName, const char * HighOccFileName);
HOCCBuilder * HOCCInitialHighOccWorkGround();
void HOCCFreeHighOccWorkGround(HOCCBuilder * hoccBuilder);

#endif
