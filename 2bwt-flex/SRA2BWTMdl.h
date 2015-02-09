//
//    SRA2BWTMdl.h
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
/////////////////////////////////////////////////////
/*

   Modification History 
   
   Date   : 8th May 2011
   Author : Edward MK Wu
   Change : New file.

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Thus, changing all references to 2bwt lib to subdirectory.

   Date   : 3rd July 2011
   Author : Edward MK Wu
   Change : Add OutputHandling to SRASetting for pair-end alignment.
   
*/
/////////////////////////////////////////////////////

#ifndef __SRA_2BWT_MODEL_H__
#define __SRA_2BWT_MODEL_H__

#include "SRACore.h"

#define MAX_NUM_OF_SRA_STEPS                         16
#define MAX_NUM_OF_SRA_CASES                         64

#define MAX_CE_THRESHOLD                             32
#define MAX_CE_BUFFER_SIZE                           SRA_MAX_READ_LENGTH / CHAR_PER_64 + 1

#define SRA_MODEL_8G                                 0
#define SRA_MODEL_16G                                1

// SRA Model Build Flag indicates how the SRA Model should be constructed.
// SRA_MODEL_BUILD_COVER_ALL constructs a model that is simple and cover
// all cases to find <= k-mismatch alignment (0<=k<=MaxError). However,
// there is no guarantee that the alignments are reported in the order of
// their number of mismatches.
// SRA_MODEL_BUILD_EXACT_NUM_ERROR constructs a model that only finds exactly
// k-mismatch alignments
// SRA_MODEL_BUILD_INCREASING_NUM_ERROR constructs a model is similar to
// COVER_ALL but the alignments are guaranteed to output in the order of their
// number of mismatches. Performance is affected inevitably.
#define SRA_MODEL_BUILD_COVER_ALL                    0
#define SRA_MODEL_BUILD_EXACT_NUM_ERROR              1
#define SRA_MODEL_BUILD_INCREASING_NUM_ERROR         2

#define SRA_MODEL_TYPE_MISMATCH_ONLY                 0
#define SRA_MODEL_TYPE_EDIT_DISTANCE                 1

static const char SRA_MODEL_TYPE[][128] = {
    "SRA_MODEL_TYPE_MISMATCH_ONLY",
    "SRA_MODEL_TYPE_EDIT_DISTANCE" };
    
#define SRA_CASE_TYPE_NOT_INITALISED                 0
#define SRA_CASE_TYPE_ALIGNMENT                      1
#define SRA_CASE_TYPE_NEXT_STAGE                     2

static const char SRA_CASE_TYPE[][128] = {
    "SRA_CASE_TYPE_NOT_INITALISED",
    "SRA_CASE_TYPE_ALIGNMENT",
    "SRA_CASE_TYPE_NEXT_STAGE" };

#define SRA_STEP_TYPE_COMPLETE                       0
#define SRA_STEP_TYPE_COMPLETE_ERROR                 1
#define SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP           2
#define SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP 3
#define SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP  4
#define SRA_STEP_TYPE_BACKWARD_ONLY_BWT              5
#define SRA_STEP_TYPE_BI_DIRECTIONAL_BWT             6

static const char SRA_STEP_TYPE[][128] = {
    "SRA_STEP_TYPE_COMPLETE",
    "SRA_STEP_TYPE_COMPLETE_ERROR",
    "SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP",
    "SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP",
    "SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP",
    "SRA_STEP_TYPE_BACKWARD_ONLY_BWT",
    "SRA_STEP_TYPE_BI_DIRECTIONAL_BWT" };

#define SRA_STEP_ERROR_TYPE_NO_ERROR                             0
#define SRA_STEP_ERROR_TYPE_MISMATCH_ONLY                        1
#define SRA_STEP_ERROR_TYPE_EDIT_DISTANCE                        2
#define SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY               3
#define SRA_STEP_ERROR_TYPE_NO_OR_NBM_MISMATCH_ONLY              4

static const char SRA_STEP_ERROR_TYPE[][128] = {
    "SRA_STEP_ERROR_TYPE_NO_ERROR",
    "SRA_STEP_ERROR_TYPE_MISMATCH_ONLY",
    "SRA_STEP_ERROR_TYPE_EDIT_DISTANCE",
    "SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY" };

//----------------------------
// DEBUG FLAGS
//----------------------------
//Define the debug flag to turn on extra output.
//#define DEBUG_2BWT
//
//Define the debug flag to disable check and extend.
//#define DEBUG_2BWT_DISABLE_CHECK_AND_EXTEND
//----------------------------


typedef struct SRAStep {
    uint8_t type;
    int start;
    int end;

    int ceThreshold;
    int ceStart;
    int ceEnd;
    
    uint8_t ErrorType;
    uint8_t MinError;
    uint8_t MaxError;

    char step;         //computed
    int len;           //computed
    BWT * bwt;         //computed
} SRAStep;

typedef struct SRACase {
    int id;
    uint8_t type;
    uint8_t ErrorType;
    uint8_t MaxError;
    SRAStep steps[MAX_NUM_OF_SRA_STEPS];
    int leftMostAligned[SRA_MAX_READ_LENGTH];
} SRACase;

typedef struct SRAModel {
    unsigned int ModelReadLength;
    SRACase     cases[MAX_NUM_OF_SRA_CASES];
} SRAModel;

void DebugPrintModel(SRAModel * SRAMismatchModel,FILE * stream);

void SRAModelFree(SRAModel * sraModel);


int SRAEditCasesPopulate(SRAModel * SRAMismatchModel, int nextCase, 
                             unsigned int ReadLength, uint8_t MaxError, uint8_t smallErrorFirst, 
                              BWT * bwt, BWT * rev_bwt);
                              
int SRAMismatchCasesPopulate(SRAModel * SRAMismatchModel, int nextCase, 
                             unsigned int ReadLength, uint8_t MaxError, uint8_t smallErrorFirst, 
                              BWT * bwt, BWT * rev_bwt);

int SRAEditCasesPopulate8G(SRAModel * SRAMismatchModel, int nextCase, 
                             unsigned int ReadLength, uint8_t MaxError, uint8_t smallErrorFirst, 
                              BWT * bwt, BWT * rev_bwt);
int SRAMismatchCasesPopulate8G(SRAModel * SRAMismatchModel, int nextCase, 
                             unsigned int ReadLength, uint8_t MaxError, uint8_t smallErrorFirst, 
                              BWT * bwt, BWT * rev_bwt);
 
SRAModel *  SRAModelConstruct(unsigned int ReadLength, int ReadStrand,
                              SRASetting * qSettings, SRAIndex * aIndex, int modelId);
                              
void SRAModelPopulate(SRAModel * SRAMismatchModel,
                        unsigned int ReadLength, int ReadStrand,
                        SRASetting * qSettings, SRAIndex * aIndex, int modelId);

#endif
