//
//    SRA2BWTMdl.c
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

*/
/////////////////////////////////////////////////////

#include "SRA2BWTMdl.h"

//Define DEBUG_2BWT_DISABLE_LOOKUP_TABLE to
// disable the use of lookup table in BWT search.
//#define DEBUG_2BWT_DISABLE_LOOKUP_TABLE

//Define DEBUG_2BWT_DISABLE_CHECK_AND_EXTEND to
// disable the use of check and extend in BWT search.
//#define DEBUG_2BWT_DISABLE_CHECK_AND_EXTEND

void DebugPrintModel ( SRAModel * SRAMismatchModel, FILE * stream )
{
    int i, j, k;
    fprintf ( stream, "Debug printing alignment model:\n" );

    for ( j = 0; j < MAX_NUM_OF_SRA_CASES; j++ )
    {
        if ( SRAMismatchModel->cases[j].type == SRA_CASE_TYPE_NOT_INITALISED ) { continue; }

        fprintf ( stream, "    Case %u\n", j );
        fprintf ( stream, "    Case Type %s\n", SRA_CASE_TYPE[SRAMismatchModel->cases[j].type] );

        for ( k = 0; k < MAX_NUM_OF_SRA_STEPS; k++ )
        {
            fprintf ( stream, "        Step type : %s\n", SRA_STEP_TYPE[SRAMismatchModel->cases[j].steps[k].type] );
            fprintf ( stream, "        Step ErrorType : %s\n", SRA_STEP_ERROR_TYPE[SRAMismatchModel->cases[j].steps[k].ErrorType] );

            if ( SRAMismatchModel->cases[j].steps[k].type == SRA_STEP_TYPE_COMPLETE )
            { break; }

            fprintf ( stream, "        From %d to %d allowing mismatch %d-%d\n",
                      SRAMismatchModel->cases[j].steps[k].start,
                      SRAMismatchModel->cases[j].steps[k].end,
                      SRAMismatchModel->cases[j].steps[k].MinError,
                      SRAMismatchModel->cases[j].steps[k].MaxError );
            fprintf ( stream, "        CE enabled %d-%d (range %d)\n",
                      SRAMismatchModel->cases[j].steps[k].ceStart,
                      SRAMismatchModel->cases[j].steps[k].ceEnd,
                      SRAMismatchModel->cases[j].steps[k].ceThreshold );
            fprintf ( stream, "        Computed value :\n" );
            fprintf ( stream, "        - step %u\n", SRAMismatchModel->cases[j].steps[k].step );
            fprintf ( stream, "        - len %u\n", SRAMismatchModel->cases[j].steps[k].len );
            fprintf ( stream, "        - bwt (ptr) %llu\n", ( unsigned long long ) SRAMismatchModel->cases[j].steps[k].bwt );
        }

        // Print left aligned for CE
        //for (k=0;k<SRAMismatchModel->ModelReadLength;k++) {
        //    fprintf(stream,"%d ",SRAMismatchModel->cases[j].leftMostAligned[k]);
        //}
        //fprintf(stream,"\n");
    }
}


void SRAModelFree ( SRAModel * sraModel )
{
    if ( sraModel == NULL ) { return; }

    free ( sraModel );
}


int SRAEditCasesPopulate ( SRAModel * SRAMismatchModel, int nextCase,
                           unsigned int ReadLength, unsigned char MaxError, uint8_t buildMode,
                           BWT * bwt, BWT * rev_bwt )
{
    int firstCase = nextCase;
    int j, k;
    unsigned int region_m1a_start = 0;
    unsigned int region_m1a_end = ReadLength - ( ReadLength * 0.5f );
    unsigned int region_m1b_start = region_m1a_end + 1;
    unsigned int region_m1b_end = ReadLength - 1;
    unsigned int region_m2a_start = 0;
    unsigned int region_m2a_end = ( ReadLength * 0.3f );
    unsigned int region_m2b_start = region_m2a_end + 1;
    unsigned int region_m2b_end = region_m2b_start + ( ReadLength * 0.3f ) - 1;
    unsigned int region_m2c_start = region_m2b_end + 1;
    unsigned int region_m2c_end = ReadLength - 1;
    unsigned int region_m3a_start = 0;
    unsigned int region_m3a_end = ( ReadLength * 0.25f );
    unsigned int region_m3b_start = region_m3a_end + 1;
    unsigned int region_m3b_end = region_m3b_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3c_start = region_m3b_end + 1;
    unsigned int region_m3c_end = region_m3c_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3d_start = region_m3c_end + 1;
    unsigned int region_m3d_end = ReadLength - 1;
    unsigned int region_m4a_start = 0;
    unsigned int region_m4a_end = ( ReadLength * 0.2f );
    unsigned int region_m4b_start = region_m4a_end + 1;
    unsigned int region_m4b_end = region_m4b_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4c_start = region_m4b_end + 1;
    unsigned int region_m4c_end = region_m4c_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4d_start = region_m4c_end + 1;
    unsigned int region_m4d_end = region_m4d_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4e_start = region_m4d_end + 1;
    unsigned int region_m4e_end = ReadLength - 1;
    unsigned int region_m5a_start = 0;
    unsigned int region_m5a_end = ( ReadLength * 0.16f );
    unsigned int region_m5b_start = region_m5a_end + 1;
    unsigned int region_m5b_end = region_m5b_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5c_start = region_m5b_end + 1;
    unsigned int region_m5c_end = region_m5c_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5d_start = region_m5c_end + 1;
    unsigned int region_m5d_end = region_m5d_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5e_start = region_m5d_end + 1;
    unsigned int region_m5e_end = region_m5e_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5f_start = region_m5e_end + 1;
    unsigned int region_m5f_end = ReadLength - 1;
    


    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 0 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].MaxError = 0;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = ReadLength - 1;
        SRAMismatchModel->cases[nextCase].steps[0].end   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].ceThreshold   = 1;
        SRAMismatchModel->cases[nextCase].steps[0].ceStart       = 10;
        SRAMismatchModel->cases[nextCase].steps[0].ceEnd         = 34;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 1 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 11;
        SRAMismatchModel->cases[nextCase].MaxError = 1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        //SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m1a_end-2;
        //SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m1a_end;
        //SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 12;
        SRAMismatchModel->cases[nextCase].MaxError = 1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        //SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m1b_start;
        //SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m1b_start+2;
        //SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;//*/
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 2 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 21;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m2b_end - 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 22;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = 0;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m2c_start + 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 23;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = 0;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 24;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 3 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 31;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 32;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 33;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 34;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 35;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 36;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 4 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 41;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 42;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 43;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 44;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 45;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 46;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 47;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 48;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 49;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 410;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 5 ) ) 
    {
        // ___________________________________________
        // |_____________5_____________|______0______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 51;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |_____________0_____________|______5______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______0______|______1______|______4______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______1______|______0______|______4______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______ ______4______ ______|__0___|__1___|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______ ______4______ ______|__1___|__0___|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______0______|______2______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______2______|______0______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |___0__|__1___|______1______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |___1__|__0___|______1______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______0______|______3______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______3______|______0______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |__0___|__1___|______2______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |__1___|__0___|______2______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______2______|__0___|__1___|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______2______|__1___|__0___|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    return nextCase;
}

int SRAMismatchCasesPopulate ( SRAModel * SRAMismatchModel, int nextCase,
                               unsigned int ReadLength, unsigned char MaxError, uint8_t buildMode,
                               BWT * bwt, BWT * rev_bwt )
{
    int firstCase = nextCase;
    int j, k;
    unsigned int region_m1a_start = 0;
    unsigned int region_m1a_end = ( ReadLength * 0.5f ) - 1;
    unsigned int region_m1b_start = region_m1a_end + 1;
    unsigned int region_m1b_end = ReadLength - 1;
    unsigned int region_m2a_start = 0;
    unsigned int region_m2a_end = ( ReadLength * 0.3f ) - 1;
    unsigned int region_m2b_start = region_m2a_end + 1;
    unsigned int region_m2b_end = region_m2b_start + ( ReadLength * 0.3f ) - 1;
    unsigned int region_m2c_start = region_m2b_end + 1;
    unsigned int region_m2c_end = ReadLength - 1;
    unsigned int region_m3a_start = 0;
    unsigned int region_m3a_end = ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3b_start = region_m3a_end + 1;
    unsigned int region_m3b_end = region_m3b_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3c_start = region_m3b_end + 1;
    unsigned int region_m3c_end = region_m3c_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3d_start = region_m3c_end + 1;
    unsigned int region_m3d_end = ReadLength - 1;
    unsigned int region_m4a_start = 0;
    unsigned int region_m4a_end = ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4b_start = region_m4a_end + 1;
    unsigned int region_m4b_end = region_m4b_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4c_start = region_m4b_end + 1;
    unsigned int region_m4c_end = region_m4c_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4d_start = region_m4c_end + 1;
    unsigned int region_m4d_end = region_m4d_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4e_start = region_m4d_end + 1;
    unsigned int region_m4e_end = ReadLength - 1;
    unsigned int region_m5a_start = 0;
    unsigned int region_m5a_end = ( ReadLength * 0.16f );
    unsigned int region_m5b_start = region_m5a_end + 1;
    unsigned int region_m5b_end = region_m5b_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5c_start = region_m5b_end + 1;
    unsigned int region_m5c_end = region_m5c_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5d_start = region_m5c_end + 1;
    unsigned int region_m5d_end = region_m5d_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5e_start = region_m5d_end + 1;
    unsigned int region_m5e_end = region_m5e_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5f_start = region_m5e_end + 1;
    unsigned int region_m5f_end = ReadLength - 1;

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 0 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].MaxError = 0;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = ReadLength - 1;
        SRAMismatchModel->cases[nextCase].steps[0].end   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 1 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 11;
        SRAMismatchModel->cases[nextCase].MaxError = 1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 12;
        SRAMismatchModel->cases[nextCase].MaxError = 1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 2 ) ) 
    {
        // | 0.3 | 0.3 | 0.4 |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 21;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart   = region_m2b_end - 9;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 22;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = 0;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd   = region_m2c_start + 9;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 23;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = 0;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart   = ( region_m2a_end >= 17 ) ? SRAMismatchModel->cases[nextCase].steps[1].ceStart : 17;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].ceEnd   = region_m2c_start + 9;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 24;
        SRAMismatchModel->cases[nextCase].MaxError = 2;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd   = ( region_m2b_end - region_m2b_start + 1 >= 17 ) ? SRAMismatchModel->cases[nextCase].steps[1].ceEnd : region_m2b_end - 18;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].ceEnd   = region_m2c_start + 9;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 3 ) ) 
    {
        // | 0.25 | 0.25 | 0.25 | 0.25 |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 31;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 32;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 33;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 35;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 36;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 34;
        SRAMismatchModel->cases[nextCase].MaxError = 3;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 4 ) ) 
    {
        // | 0.2 | 0.2 | 0.2 | 0.2 | 0.2 |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 41;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 42;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 43;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 44;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 45;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 46;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 47;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 48;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 49;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 410;
        SRAMismatchModel->cases[nextCase].MaxError = 4;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 5 ) ) 
    {
        // | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
        // ___________________________________________
        // |_____________5_____________|______0______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 51;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |_____________0_____________|______5______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______0______|______1______|______4______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______1______|______0______|______4______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______ ______4______ ______|__0___|__1___|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______ ______4______ ______|__1___|__0___|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______0______|______2______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______2______|______0______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |___0__|__1___|______1______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |___1__|__0___|______1______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______0______|______3______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______3______|______0______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |__0___|__1___|______2______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |__1___|__0___|______2______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______2______|__0___|__1___|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______2______|__1___|__0___|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].MaxError = 5;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    return nextCase;
}


int SRAEditCasesPopulate8G ( SRAModel * SRAMismatchModel, int nextCase,
                             unsigned int ReadLength, unsigned char MaxError, uint8_t buildMode,
                             BWT * bwt, BWT * rev_bwt )
{
    int firstCase = nextCase;
    int j, k;
    unsigned int region_m1a_start = 0;
    unsigned int region_m1a_end = ReadLength - ( ReadLength * 0.5f );
    unsigned int region_m1b_start = region_m1a_end + 1;
    unsigned int region_m1b_end = ReadLength - 1;
    unsigned int region_m2a_start = 0;
    unsigned int region_m2a_end = ( ReadLength * 0.3f );
    unsigned int region_m2b_start = region_m2a_end + 1;
    unsigned int region_m2b_end = region_m2b_start + ( ReadLength * 0.3f ) - 1;
    unsigned int region_m2c_start = region_m2b_end + 1;
    unsigned int region_m2c_end = ReadLength - 1;
    unsigned int region_m3a_start = 0;
    unsigned int region_m3a_end = ( ReadLength * 0.25f );
    unsigned int region_m3b_start = region_m3a_end + 1;
    unsigned int region_m3b_end = region_m3b_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3c_start = region_m3b_end + 1;
    unsigned int region_m3c_end = region_m3c_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3d_start = region_m3c_end + 1;
    unsigned int region_m3d_end = ReadLength - 1;
    unsigned int region_m4a_start = 0;
    unsigned int region_m4a_end = ( ReadLength * 0.2f );
    unsigned int region_m4b_start = region_m4a_end + 1;
    unsigned int region_m4b_end = region_m4b_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4c_start = region_m4b_end + 1;
    unsigned int region_m4c_end = region_m4c_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4d_start = region_m4c_end + 1;
    unsigned int region_m4d_end = region_m4d_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4e_start = region_m4d_end + 1;
    unsigned int region_m4e_end = ReadLength - 1;
    unsigned int region_m5a_start = 0;
    unsigned int region_m5a_end = ( ReadLength * 0.16f );
    unsigned int region_m5b_start = region_m5a_end + 1;
    unsigned int region_m5b_end = region_m5b_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5c_start = region_m5b_end + 1;
    unsigned int region_m5c_end = region_m5c_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5d_start = region_m5c_end + 1;
    unsigned int region_m5d_end = region_m5d_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5e_start = region_m5d_end + 1;
    unsigned int region_m5e_end = region_m5e_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5f_start = region_m5e_end + 1;
    unsigned int region_m5f_end = ReadLength - 1;

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 0 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = ReadLength - 1;
        SRAMismatchModel->cases[nextCase].steps[0].end   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].ceThreshold   = 1;
        SRAMismatchModel->cases[nextCase].steps[0].ceStart       = 10;
        SRAMismatchModel->cases[nextCase].steps[0].ceEnd         = 34;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 1 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 11;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        //SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m1a_end-2;
        //SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m1a_end;
        //SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 12;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        //SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m1b_start;
        //SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m1b_start+2;
        //SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;//*/
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 2 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 21;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m2b_end - 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 22;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = 0;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m2c_start + 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 23;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = 0;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 24;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 3 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 31;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 32;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 33;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 34;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 35;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 36;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 4 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 41;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 42;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 43;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 44;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 45;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 46;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 47;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 48;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 49;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 410;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 5 ) ) 
    {
        // ___________________________________________
        // |_____________5_____________|______0______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 51;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |_____________0_____________|______5______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______0______|______1______|______4______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______1______|______0______|______4______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______ ______4______ ______|__0___|__1___|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______ ______4______ ______|__1___|__0___|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______0______|______2______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______2______|______0______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |___0__|__1___|______1______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |___1__|__0___|______1______|______3______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______0______|______3______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______3______|______0______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |__0___|__1___|______2______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |__1___|__0___|______2______|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______2______|__0___|__1___|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        // ___________________________________________
        // |______2______|__1___|__0___|______2______|
        // |      |      |      |      |      |      |
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    return nextCase;
}

int SRAMismatchCasesPopulate8G ( SRAModel * SRAMismatchModel, int nextCase,
                                 unsigned int ReadLength, unsigned char MaxError, uint8_t buildMode,
                                 BWT * bwt, BWT * rev_bwt )
{
    int firstCase = nextCase;
    int j, k;
    unsigned int region_m1a_start = 0;
    unsigned int region_m1a_end = ( ReadLength * 0.5f ) - 1;
    unsigned int region_m1b_start = region_m1a_end + 1;
    unsigned int region_m1b_end = ReadLength - 1;
    unsigned int region_m2a_start = 0;
    unsigned int region_m2a_end = ( ReadLength * 0.3f ) - 1;
    unsigned int region_m2b_start = region_m2a_end + 1;
    unsigned int region_m2b_end = region_m2b_start + ( ReadLength * 0.3f ) - 1;
    unsigned int region_m2c_start = region_m2b_end + 1;
    unsigned int region_m2c_end = ReadLength - 1;
    unsigned int region_m3a_start = 0;
    unsigned int region_m3a_end = ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3b_start = region_m3a_end + 1;
    unsigned int region_m3b_end = region_m3b_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3c_start = region_m3b_end + 1;
    unsigned int region_m3c_end = region_m3c_start + ( ReadLength * 0.25f ) - 1;
    unsigned int region_m3d_start = region_m3c_end + 1;
    unsigned int region_m3d_end = ReadLength - 1;
    unsigned int region_m4a_start = 0;
    unsigned int region_m4a_end = ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4b_start = region_m4a_end + 1;
    unsigned int region_m4b_end = region_m4b_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4c_start = region_m4b_end + 1;
    unsigned int region_m4c_end = region_m4c_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4d_start = region_m4c_end + 1;
    unsigned int region_m4d_end = region_m4d_start + ( ReadLength * 0.2f ) - 1;
    unsigned int region_m4e_start = region_m4d_end + 1;
    unsigned int region_m4e_end = ReadLength - 1;
    unsigned int region_m5a_start = 0;
    unsigned int region_m5a_end = ( ReadLength * 0.16f );
    unsigned int region_m5b_start = region_m5a_end + 1;
    unsigned int region_m5b_end = region_m5b_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5c_start = region_m5b_end + 1;
    unsigned int region_m5c_end = region_m5c_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5d_start = region_m5c_end + 1;
    unsigned int region_m5d_end = region_m5d_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5e_start = region_m5d_end + 1;
    unsigned int region_m5e_end = region_m5e_start + ( ReadLength * 0.16f ) - 1;
    unsigned int region_m5f_start = region_m5e_end + 1;
    unsigned int region_m5f_end = ReadLength - 1;

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 0 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = ReadLength - 1;
        SRAMismatchModel->cases[nextCase].steps[0].end   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].ceStart    = 9;
        SRAMismatchModel->cases[nextCase].steps[0].ceThreshold = 1;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;//*/
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 1 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 11;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m1a_end - 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 12;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart   = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd   = region_m1b_start + 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;//*/
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 2 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 2 ) ) 
    {
        // | 0.3 | 0.3 | 0.4 |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 21;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart   = region_m2b_end - 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd     = region_m2b_end;//-1;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;//*/
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 22;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart       = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd         = region_m2c_start + 4;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;//*/
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 23;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceStart   = ( region_m2a_end >= 17 ) ? SRAMismatchModel->cases[nextCase].steps[1].ceStart : 17;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[2].ceEnd   = region_m2c_start + 9;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;//*/
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 24;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m2b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m2b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m2c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m2c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].ceEnd   = ( region_m2b_end - region_m2b_start + 1 >= 17 ) ? SRAMismatchModel->cases[nextCase].steps[1].ceEnd : region_m2b_end - 18;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m2a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m2a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[2].ceEnd   = region_m2c_start + 9;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;//*/
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 3 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 3 ) ) 
    {
        // | 0.25 | 0.25 | 0.25 | 0.25 |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 31;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 32;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 33;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 35;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 36;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 34;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m3b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m3b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m3a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m3a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m3c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m3d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 4 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 4 ) ) 
    {
        // | 0.2 | 0.2 | 0.2 | 0.2 | 0.2 |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 41;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 42;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 43;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 44;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 45;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 46;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 47;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 48;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4b_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4a_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 49;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 410;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m4e_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m4e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m4d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m4d_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m4c_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m4a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 5 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 5 ) ) 
    {
        // | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 | 0.16 |
        // ___________________________________________
        // |_____________5_____________|______0______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 51;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |_____________0_____________|______5______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______0______|______1______|______4______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______1______|______0______|______4______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______ ______4______ ______|__0___|__1___|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______ ______4______ ______|__1___|__0___|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5f_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5e_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 4;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______0______|______2______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______2______|______0______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |___0__|__1___|______1______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |___1__|__0___|______1______|______3______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______0______|______3______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______3______|______0______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 3;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |__0___|__1___|______2______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |__1___|__0___|______2______|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______2______|__0___|__1___|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        // ___________________________________________
        // |______2______|__1___|__0___|______2______|
        // |      |      |      |      |      |      |
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id = 52;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m5d_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m5d_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m5c_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m5c_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 10;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[2].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[2].start = region_m5b_end;
        SRAMismatchModel->cases[nextCase].steps[2].end   = region_m5a_start;
        SRAMismatchModel->cases[nextCase].steps[2].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[2].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[3].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[3].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[3].start = region_m5e_start;
        SRAMismatchModel->cases[nextCase].steps[3].end   = region_m5f_end;
        SRAMismatchModel->cases[nextCase].steps[3].MinError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].MaxError   = 2;
        SRAMismatchModel->cases[nextCase].steps[3].ceThreshold   = 0;
        SRAMismatchModel->cases[nextCase].steps[4].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        SRAMismatchModel->cases[nextCase].id = 1;
        SRAMismatchModel->cases[nextCase].type = SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    return nextCase;
}



int SRANBMismatchCasesPopulate(SRAModel * SRAMismatchModel, int nextCase, 
                             unsigned int ReadLength, int ReadStrand, unsigned char MaxError, uint8_t buildMode, 
                              BWT * bwt, BWT * rev_bwt) {
    
    int firstCase = nextCase;
    int j,k;

    unsigned int region_m1a_start = 0;
    unsigned int region_m1a_end = ( ReadLength * 0.5f ) - 1;
    unsigned int region_m1b_start = region_m1a_end + 1;
    unsigned int region_m1b_end = ReadLength - 1;
    
    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 0 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 0 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type=SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id=1;
        SRAMismatchModel->cases[nextCase].MaxError=0;
        if (ReadStrand == QUERY_POS_STRAND) {
            SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
            SRAMismatchModel->cases[nextCase].steps[0].start = 0;
            SRAMismatchModel->cases[nextCase].steps[0].end   = ReadLength - 1;
        } else {
            SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP;
            SRAMismatchModel->cases[nextCase].steps[0].start = ReadLength - 1;
            SRAMismatchModel->cases[nextCase].steps[0].end   = 0;
        }
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].ceThreshold   = 5;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        
        SRAMismatchModel->cases[nextCase].id=1;
        SRAMismatchModel->cases[nextCase].type=SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }

    if ( ( buildMode == SRA_MODEL_BUILD_COVER_ALL && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR && MaxError == 1 ) ||
         ( buildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR && MaxError >= 1 ) ) 
    {
        SRAMismatchModel->cases[nextCase].type=SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id=11;
        SRAMismatchModel->cases[nextCase].MaxError=1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;

        SRAMismatchModel->cases[nextCase].type=SRA_CASE_TYPE_ALIGNMENT;
        SRAMismatchModel->cases[nextCase].id=12;
        SRAMismatchModel->cases[nextCase].MaxError=1;
        SRAMismatchModel->cases[nextCase].steps[0].type  = SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP;
        SRAMismatchModel->cases[nextCase].steps[0].start = region_m1b_start;
        SRAMismatchModel->cases[nextCase].steps[0].end   = region_m1b_end;
        SRAMismatchModel->cases[nextCase].steps[0].MinError   = 0;
        SRAMismatchModel->cases[nextCase].steps[0].MaxError   = 0;
        SRAMismatchModel->cases[nextCase].steps[1].type  = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
        SRAMismatchModel->cases[nextCase].steps[1].ErrorType  = SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        SRAMismatchModel->cases[nextCase].steps[1].start = region_m1a_end;
        SRAMismatchModel->cases[nextCase].steps[1].end   = region_m1a_start;
        SRAMismatchModel->cases[nextCase].steps[1].MinError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].MaxError   = 1;
        SRAMismatchModel->cases[nextCase].steps[1].ceThreshold   = 20;
        SRAMismatchModel->cases[nextCase].steps[2].type  = SRA_STEP_TYPE_COMPLETE;
        nextCase++;
        
        SRAMismatchModel->cases[nextCase].id=1;
        SRAMismatchModel->cases[nextCase].type=SRA_CASE_TYPE_NEXT_STAGE;
        nextCase++;
    }
}

SRAModel  * SRAModelConstruct ( unsigned int ReadLength, int ReadStrand,
                                SRASetting * aSettings, SRAIndex * aIndex, int modelId )
{
    SRAModel * SRAMismatchModel = ( SRAModel * ) malloc ( sizeof ( SRAModel ) );
    SRAModelPopulate ( SRAMismatchModel, ReadLength, ReadStrand, aSettings, aIndex, modelId );
    return SRAMismatchModel;
}

void SRAModelPopulate ( SRAModel * SRAMismatchModel,
                        unsigned int ReadLength, int ReadStrand,
                        SRASetting * aSettings, SRAIndex * aIndex, int modelId )
{
#ifdef DEBUG_2BWT
    printf ( "\n[SRAModelConstruct] Size of SRAModel = %u\n", sizeof ( SRAModel ) );
#endif
    BWT * bwt = aIndex->bwt;
    BWT * rev_bwt = aIndex->rev_bwt;
    uint8_t ErrorType = aSettings->ErrorType;
    uint8_t MaxError = aSettings->MaxError;
    uint8_t sraModelBuildMode = SRA_MODEL_BUILD_COVER_ALL;

    // SOAP3-DP COMMENT OUT
    // The following codes are commented out as ALL_BEST
    // handling in SOAP3-DP is different to SOAP2-DP
    if ( aSettings->OutputType == SRA_REPORT_EXACT_NUM_ERROR ) {
        sraModelBuildMode = SRA_MODEL_BUILD_EXACT_NUM_ERROR;
    }
    /*
    if ( aSettings->OutputType == SRA_REPORT_UNIQUE_BEST ||
            aSettings->OutputType == SRA_REPORT_RANDOM_BEST ||
            aSettings->OutputType == SRA_REPORT_ALL_BEST )
    {
        SmallErrorFirst = 1;
    }
    */
    // END OF SOAP3-DP COMMENT OUT

    //DIU
    //SmallErrorFirst=1;
    //EODIU
    //_____________________________________________________________________________________________________
    //Initialize it
    int i, j, k, l;
    SRAMismatchModel->ModelReadLength = ReadLength;

    if ( ReadLength == 0 ) { return; }

    for ( j = 0; j < MAX_NUM_OF_SRA_CASES; j++ )
    {
        SRAMismatchModel->cases[j].id = 0;
        SRAMismatchModel->cases[j].type = SRA_CASE_TYPE_NOT_INITALISED;
        SRAMismatchModel->cases[j].ErrorType = ErrorType;

        for ( k = 0; k < MAX_NUM_OF_SRA_STEPS; k++ )
        {
            SRAMismatchModel->cases[j].steps[k].type = SRA_STEP_TYPE_COMPLETE;
            SRAMismatchModel->cases[j].steps[k].start = 0;
            SRAMismatchModel->cases[j].steps[k].end = 0;
            SRAMismatchModel->cases[j].steps[k].ceThreshold = 0;
            SRAMismatchModel->cases[j].steps[k].ceStart = 0;
            SRAMismatchModel->cases[j].steps[k].ceEnd = ReadLength - 1;
            SRAMismatchModel->cases[j].steps[k].ErrorType = SRA_STEP_ERROR_TYPE_NO_ERROR;
            SRAMismatchModel->cases[j].steps[k].MinError = 0;
            SRAMismatchModel->cases[j].steps[k].MaxError = 0;
        }
    }

    j = 0;

    if ( ErrorType == SRA_STEP_ERROR_TYPE_NO_ERROR || ErrorType == SRA_STEP_ERROR_TYPE_MISMATCH_ONLY )
    {
        if ( modelId == SRA_MODEL_8G )
        {
#ifdef DEBUG_2BWT
            printf ( "[SRAModelConstruct] Constructed SRA_MODEL_8G\n", j );
#endif
            //Direction Heuristic for NBMismatch matching
            //
            if (aSettings->MaxNBMismatch>0 && MaxError<=1) {
                j=SRANBMismatchCasesPopulate(SRAMismatchModel,j,ReadLength,ReadStrand,MaxError,sraModelBuildMode,bwt,rev_bwt);
            } else {
                j=SRAMismatchCasesPopulate8G(SRAMismatchModel,j,ReadLength,MaxError,sraModelBuildMode,bwt,rev_bwt);
            }
        }
        else
        {
#ifdef DEBUG_2BWT
            printf ( "[SRAModelConstruct] Constructed SRA_MODEL_16G\n", j );
#endif 
            if (aSettings->MaxNBMismatch>0 && MaxError<=1) {
                j=SRANBMismatchCasesPopulate(SRAMismatchModel,j,ReadLength,ReadStrand,MaxError,sraModelBuildMode,bwt,rev_bwt);
            } else {
                j=SRAMismatchCasesPopulate(SRAMismatchModel,j,ReadLength,MaxError,sraModelBuildMode,bwt,rev_bwt);
            }
        }
    }
    else if ( ErrorType == SRA_STEP_ERROR_TYPE_EDIT_DISTANCE )
    {
        if ( modelId == SRA_MODEL_8G )
        {
#ifdef DEBUG_2BWT
            printf ( "[SRAModelConstruct] Constructed SRA_MODEL_8G\n", j );
#endif
            j = SRAEditCasesPopulate8G ( SRAMismatchModel, j, ReadLength, MaxError, sraModelBuildMode, bwt, rev_bwt );
        }
        else
        {
#ifdef DEBUG_2BWT
            printf ( "[SRAModelConstruct] Constructed SRA_MODEL_16G\n", j );
#endif
            j = SRAEditCasesPopulate ( SRAMismatchModel, j, ReadLength, MaxError, sraModelBuildMode, bwt, rev_bwt );
        }
    }

#ifdef DEBUG_2BWT
    printf ( "[SRAModelConstruct] Constructed a %d cases model\n", j );
#endif

    //_____________________________________________________________________________________________________
    //Post-compute
    for ( j = 0; j < MAX_NUM_OF_SRA_CASES; j++ )
    {
        if ( SRAMismatchModel->cases[j].type != SRA_CASE_TYPE_NOT_INITALISED &&
            SRAMismatchModel->cases[j].type != SRA_CASE_TYPE_NEXT_STAGE)
        {
            int leftMost = -1;

            for ( k = 0; k < MAX_NUM_OF_SRA_STEPS; k++ )
            {
                if ( SRAMismatchModel->cases[j].steps[k].type != SRA_STEP_TYPE_COMPLETE )
                {
                    int start = SRAMismatchModel->cases[j].steps[k].start;
                    int end = SRAMismatchModel->cases[j].steps[k].end;
                    int step = -2 * ( start > end ) + 1;
                    SRAMismatchModel->cases[j].steps[k].step = step;
                    SRAMismatchModel->cases[j].steps[k].len = ( end - start ) * step + 1;

                    if ( step > 0 )
                    {
                        //Forward alignment
                        if (leftMost == -1) {
                            leftMost = start;
                        }
                        SRAMismatchModel->cases[j].steps[k].bwt = rev_bwt;

                        for ( l = start; l <= end; l += step )
                        {
                            SRAMismatchModel->cases[j].leftMostAligned[l] = leftMost;
                        }
                    }
                    else
                    {
                        //Backward alignment
                        if (leftMost == -1) {
                            leftMost = start+1;
                        }
                        SRAMismatchModel->cases[j].steps[k].bwt = bwt;

                        for ( l = start; l >= end; l += step )
                        {
                            SRAMismatchModel->cases[j].leftMostAligned[l] = leftMost--;
                        }
                    }

                    if ( sraModelBuildMode == SRA_MODEL_BUILD_EXACT_NUM_ERROR || 
                         sraModelBuildMode == SRA_MODEL_BUILD_INCREASING_NUM_ERROR )
                    {
                        SRAMismatchModel->cases[j].steps[k].MinError = SRAMismatchModel->cases[j].steps[k].MaxError;
                    }

                    if ( SRAMismatchModel->cases[j].steps[k].ceThreshold > MAX_CE_THRESHOLD )
                    {
                        SRAMismatchModel->cases[j].steps[k].ceThreshold = MAX_CE_THRESHOLD;
                    }

                    if ( SRAMismatchModel->cases[j].steps[k].ceStart > SRAMismatchModel->cases[j].steps[k].ceEnd )
                    {
                        int swap = SRAMismatchModel->cases[j].steps[k].ceEnd;
                        SRAMismatchModel->cases[j].steps[k].ceEnd = SRAMismatchModel->cases[j].steps[k].ceStart;
                        SRAMismatchModel->cases[j].steps[k].ceStart = swap;
                    }

#ifdef DEBUG_2BWT_DISABLE_CHECK_AND_EXTEND
                    SRAMismatchModel->cases[j].steps[k].ceThreshold = 0;
#endif
#ifdef DEBUG_2BWT_DISABLE_LOOKUP_TABLE

                    if (SRAMismatchModel->cases[j].steps[k].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP ||
                        SRAMismatchModel->cases[j].steps[k].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP )
                    {
                        SRAMismatchModel->cases[j].steps[k].type = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
                    }
                    else if ( SRAMismatchModel->cases[j].steps[k].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP )
                    {
                        SRAMismatchModel->cases[j].steps[k].type = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
                    }

#endif
                }
            }
        }
    }
}
