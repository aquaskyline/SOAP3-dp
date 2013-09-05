/*
 *
 *    UsageInterface.c
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

#ifndef __INI_PARAM_H__
#define __INI_PARAM_H__

#include "UsageInterface.h"

// This function is to print out the help on the usage of the program
void UIPrintUsageOverview ( char * program_name )
{
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  %s -- Usage Guide  \n", PROJECT_NAME );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  Usage:   %s <command> <other parameters>\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  Command: single              Alignment of single-end reads.\n" );
    fprintf ( stderr, "           single-multi        Alignment of multiple sets of single-end reads.\n" );
    fprintf ( stderr, "           pair                Alignment of paired-end reads.\n" );
    fprintf ( stderr, "           pair-multi          Alignment of multiple sets of paired-end reads.\n" );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  To see the detailed usage for the specific command, please type:\n" );
    fprintf ( stderr, "           %s <command>\n", program_name );
    fprintf ( stderr, "\n" );
}




void UIPrintAlignmentOptions ( int casenum )
{
    // casenum: 1: single; 2: single-multi; 3: pair; 4: pair-multi
    fprintf ( stderr, "      options:\n" );

    switch ( casenum )
    {
        case 1:
            fprintf ( stderr, "               -o <output file prefix> (default: queryFileName/queryBAMFileName)\n" );
            break;

        case 3:
            fprintf ( stderr, "               -u <max value of insert size> (default: 500)\n" );
            fprintf ( stderr, "               -v <min value of insert size> (default: 1)\n" );
            fprintf ( stderr, "               -o <output file prefix> (default: queryFileName1/queryBAMFileName)\n" );
            break;

        case 4:
            break;
    }

    fprintf ( stderr, "               -L <length of the longest read in the input> (default: 120)\n" );
    fprintf ( stderr, "               -h <output option>  1: all valid alignments; \n" );
    fprintf ( stderr, "                                   2: all best alignments (default); \n" );
    fprintf ( stderr, "                                   3: unique best alignments; \n" );
    fprintf ( stderr, "                                   4: random best alignments; \n" );
    fprintf ( stderr, "               -b <output format>  1: Plain; \n" );
    fprintf ( stderr, "                                   2: SAM (default); \n" );
    fprintf ( stderr, "                                   3: BAM \n" );
    fprintf ( stderr, "               -I                  The input is in the Illumina 1.3+ \n" );
    fprintf ( stderr, "                                   FASTQ-like format (i.e. Phred+64)\n" );
    fprintf ( stderr, "                                   (by default, this option is not enabled\n" );
    fprintf ( stderr, "                                   and presumes Phred+33 format).\n" );
    fprintf ( stderr, "               -c <GPU device ID> \n" );
    fprintf ( stderr, "               -p                  Output MD string and NM tag.\n" );

    if ( casenum == 1 || casenum == 3 )
    {
        fprintf ( stderr, "               -A <sample name>    Assign the sample name\n" );
        fprintf ( stderr, "                                   (default: \"default\").\n" );
        fprintf ( stderr, "               -D <read group ID>  Assign the read group ID\n" );
        fprintf ( stderr, "                                   (default: queryFileName1/queryBAMFileName).\n" );
        fprintf ( stderr, "               -R <read group option> Assign the read group option\n" );
        fprintf ( stderr, "                                   Example format: 'LB:...\\tPL:...\\tPU:...'\n" );
    }

    fprintf ( stderr, "               -s <max # of mismatches> (<=4) disable dynamic programming and\n" );
    fprintf ( stderr, "                                   perform alignment with mismatches only.\n" );
    fprintf ( stderr, "                                   If max # of mismatches is not specified,\n" );
    fprintf ( stderr, "                                   default number is 3 for reads of length >= 50;\n" );
    fprintf ( stderr, "                                   2 otherwise. Using this option is not recommended.\n" );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "       Note: For SAM output format, mapping quality value is only available\n" );
    fprintf ( stderr, "             when all-best or all-valid output option is selected.\n" );
    fprintf ( stderr, "\n" );
}

void UIprintUsageSingle ( char * program_name )
{
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  Alignment of single-end reads:\n" );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in FASTA/FASTQ/GZIP format:\n" );
    fprintf ( stderr, "  %s single <bwtCodeIndex> <queryFileName> [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in BAM format:\n" );
    fprintf ( stderr, "  %s single <bwtCodeIndex> <queryBAMFileName> -bam [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "      bwtCodeIndex     : The bwt index built by %s\n", PROJECT_BUILDER_BINARY );
    fprintf ( stderr, "      queryFileName    : The FASTA/FASTQ query file\n" );
    fprintf ( stderr, "                         The read files can be compressed in GZ format.\n\n" );
    fprintf ( stderr, "      queryBAMFileName : The BAM query file\n" );
    fprintf ( stderr, "\n" );
    UIPrintAlignmentOptions ( 1 );
}

void UIprintUsageSingleList ( char * program_name )
{
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  Alignment of multiple sets of single-end reads:\n" );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in FASTA/FASTQ/GZIP format:\n" );
    fprintf ( stderr, "  %s single-multi <bwtCodeIndex> <infoListFile> [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in BAM format:\n" );
    fprintf ( stderr, "  %s single-multi <bwtCodeIndex> <infoListFile> -bam [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "      bwtCodeIndex  : The bwt index built by %s\n", PROJECT_BUILDER_BINARY );
    fprintf ( stderr, "      infoListFile  : Each line inside the file contains information of each set of paired-end reads.\n" );
    fprintf ( stderr, "                      Column 1: first read file name;          Column 2: output prefix;\n" );
    fprintf ( stderr, "                      Column 3 (optional): read group ID;      Column 4 (optional): sample name;\n" );
    fprintf ( stderr, "                      Column 5 (optional): read group options.\n" );
    fprintf ( stderr, "                      Delimiter is a tab character and only the following three cases are allowed:\n" );
    fprintf ( stderr, "                      Case 1: Only contains the first 2 columns.\n" );
    fprintf ( stderr, "                      Case 2: Only contains the first 4 columns.\n" );
    fprintf ( stderr, "                      Case 3: Contains all the 5 columns.\n" );
    fprintf ( stderr, "\n" );
    UIPrintAlignmentOptions ( 2 );
}

void UIprintUsagePair ( char * program_name )
{
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  Alignment of paired-end reads:\n" );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in FASTA/FASTQ/GZIP format:\n" );
    fprintf ( stderr, "  %s pair <bwtCodeIndex> <queryFileName1> <queryFileName2> [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in BAM format:\n" );
    fprintf ( stderr, "  %s pair <bwtCodeIndex> <queryBAMFileName> -bam [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "      bwtCodeIndex    : The bwt index built by %s\n", PROJECT_BUILDER_BINARY );
    fprintf ( stderr, "      queryFileName1  : The FASTA/FASTQ file of the first reads of pairs.\n" );
    fprintf ( stderr, "      queryFileName2  : The FASTA/FASTQ file of the second reads of pairs.\n" );
    fprintf ( stderr, "                        The read files can be compressed in GZ format.\n" );
    fprintf ( stderr, "      queryBAMFileName: The BAM file contains both the first and the second reads of the pairs.\n" );
    fprintf ( stderr, "\n" );
    UIPrintAlignmentOptions ( 3 );
}

void UIprintUsagePairList ( char * program_name )
{
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  Alignment of multiple sets of paired-end reads:\n" );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in FASTA/FASTQ/GZIP format:\n" );
    fprintf ( stderr, "  %s pair-multi <bwtCodeIndex> <infoListFile> [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  For read files in BAM format:\n" );
    fprintf ( stderr, "  %s pair-multi <bwtCodeIndex> <infoListFile> -bam [options]\n", program_name );
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "      bwtCodeIndex    : The bwt index built by %s\n", PROJECT_BUILDER_BINARY );
    fprintf ( stderr, "      infoListFile    : Each line inside the file contains information of each set of paired-end reads.\n" );
    fprintf ( stderr, "                        Column 1: first read file name;          Column 2: second read file name;\n" );
    fprintf ( stderr, "                        Column 3: minimum value of insert size;  Column 4: maximum value of insert size;\n" );
    fprintf ( stderr, "                        Column 5: output prefix;                 Column 6 (optional): read group ID;\n" );
    fprintf ( stderr, "                        Column 7 (optional): sample name;        Column 8 (optional): read group options.\n" );
    fprintf ( stderr, "                        Note: If the read files are in BAM format, the second read file name is not needed,\n" );
    fprintf ( stderr, "                        and thus there are at most seven columans.\n" );
    fprintf ( stderr, "                        Delimiter is a tab character and only the following three cases are allowed:\n" );
    fprintf ( stderr, "                        Case 1: Only contains the first 5 columns (or 4 columns for BAM).\n" );
    fprintf ( stderr, "                        Case 2: Only contains the first 7 columns (or 6 columns for BAM).\n" );
    fprintf ( stderr, "                        Case 3: Contains all the 8 columns (or 7 columns for BAM).\n" );
    fprintf ( stderr, "\n" );
    UIPrintAlignmentOptions ( 4 );
}

#endif
