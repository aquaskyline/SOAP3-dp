/*

   BWTLoad.c        BWTLoad - Loading BWT index into memory

   Copyright (C) 2006, Wong Chi Kwong.

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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "iniparser.h"
#include "Socket.h"
#include "MiscUtilities.h"
#include "MemManager.h"
#include "Timing.h"
#include "BWTLoad.h"

// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();
void PrintIni();
void PrintShortDesc();
void PrintHelp();
void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseLocation, const char *databaseName);

// Query parameters
void ParseQueryParameterFile(char *queryParameterFileName);
void ValidateQueryParameters();
void PrintQueryParameters();

    // Parameters
    int Confirmation = FALSE;
    int UnloadBWTIndex = FALSE;
        
    // Database parameters
    char DatabaseName[MAX_FILENAME_LEN+1] = "";
    char DatabaseLocation[MAX_FILENAME_LEN+1] = "./";

    // DatabaseFiles parameters
    char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.ann";
    char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.amb";
    char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.pac";
    char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.bwt";
    char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.fmv";
    char SaValueFileName[MAX_FILENAME_LEN+1] = "*.sa";
    char SaIndexFileName[MAX_FILENAME_LEN+1] = "*.sai";


int main(int argc, char** argv) {

    // Program input
    dictionary *programInput;

    // Socket
    Socket *tempSocket;
    Socket *bwtSocket;

    // Index
    BWTIndexPointer bwtIndexPointer;
    BWTIndexFileName bwtIndexFileName, requestBwtIndexFileName;

    // Working variables
    int c;
    char filename[MAX_FILENAME_LEN];

    // Performance statistics
    double startTime;
    double elapsedTime;

    // Program input
    programInput = ParseInput(argc, argv);
    PrintShortDesc();

    // Ini
    if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
        *(argv[0] + strlen(argv[0]) - 4) = '\0';
    }
    sprintf(filename, "%s.ini", argv[0]);            // First default ini is <program name>.ini
    ParseIniFile(filename);
    sprintf(filename, "%s.ini", DatabaseName);        // Second default ini is <database name>.ini
    ParseIniFile(filename);
    printf("\n");
    ProcessIni();
    ValidateIni();
    PrintIni();
    iniparser_freedict(programInput);

    memcpy(bwtIndexFileName.AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN+1);
    memcpy(bwtIndexFileName.AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN+1);
    memcpy(bwtIndexFileName.PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN+1);
    memcpy(bwtIndexFileName.BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN+1);
    memcpy(bwtIndexFileName.BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN+1);
    memcpy(bwtIndexFileName.SaValueFileName, SaValueFileName, MAX_FILENAME_LEN+1);
    memcpy(bwtIndexFileName.SaIndexFileName, SaIndexFileName, MAX_FILENAME_LEN+1);

    if (UnloadBWTIndex) {
        printf("Unloading BWT index.\n");
    }
    if (Confirmation == TRUE) {
        printf("Press Y to go or N to cancel. ");
        c = (char)getchar();
        while (c != 'y' && c != 'Y' && c != 'n' && c!= 'N') {
            c = (char)getchar();
        }
        if (c == 'n' || c == 'N') {
            exit(0);
        }
    }

    // Measure searching performance
    startTime = setStartTime();

    // First try to connect to the resident process
    tempSocket = SocketInitiateConnection(BWT_INDEX_SOCKET_TYPE, BWT_INDEX_SOCKET_NAME);
    if (tempSocket != NULL) {
        printf("Communicating with the resident process..");
        if (UnloadBWTIndex) {
            bwtIndexFileName.BWTCodeFileName[0] = '\0';
            if (SocketSend(tempSocket, &bwtIndexFileName, sizeof(BWTIndexFileName)) < 0) {
                fprintf(stderr, "BWTLoad(): Socket communication error!\n");
                exit(1);
            }
        } else {
            if (SocketSend(tempSocket, &bwtIndexFileName, sizeof(BWTIndexFileName)) < 0) {
                fprintf(stderr, "BWTLoad(): Socket communication error!\n");
                exit(1);
            }
        }
        if (SocketReceive(tempSocket, &bwtIndexPointer, sizeof(BWTIndexPointer)) < 0) {
            fprintf(stderr, "BWTLoad(): Error reloading BWT index!\n");
            exit(1);
        }
        SocketEndConnection(tempSocket);
        exit(0);
    } else {
        if (UnloadBWTIndex) {
            fprintf(stderr, "BWTLoad(): Cannot communicate with the resident process!\n");
            exit(1);
        }
    }

    // Load Database
    printf("Loading Database..");
    fflush(stdout);

    BWTBlastLoad(NULL, BWTCodeFileName, BWTOccValueFileName, SaValueFileName, SaIndexFileName,
                 PackedDNAFileName, AnnotationFileName, AmbiguityFileName,
                 &bwtIndexPointer.bwt, &bwtIndexPointer.saIndexRange, &bwtIndexPointer.saIndexRangeNumOfChar,
                 &bwtIndexPointer.packedDNA, &bwtIndexPointer.annotation, &bwtIndexPointer.seqOffset, &bwtIndexPointer.numOfSeq,
                 &bwtIndexPointer.ambiguity, &bwtIndexPointer.numOfAmbiguity, &bwtIndexPointer.minSeqLength);

    printf("\n%u %u %d %u %u %u %d %u %d %u\n", bwtIndexPointer.bwt, bwtIndexPointer.saIndexRange, bwtIndexPointer.saIndexRangeNumOfChar, bwtIndexPointer.packedDNA, bwtIndexPointer.annotation, bwtIndexPointer.seqOffset, bwtIndexPointer.numOfSeq, bwtIndexPointer.ambiguity, bwtIndexPointer.numOfAmbiguity, bwtIndexPointer.minSeqLength);

    // Finished loading
    elapsedTime = getElapsedTime(startTime);

    printf("Finished.\n");
    printf("There are %d sequences in the database.\n\n", bwtIndexPointer.numOfSeq);
    printf("Elapsed time = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
    printf("\n");

    // Create socket
    bwtSocket = SocketCreate(BWT_INDEX_SOCKET_TYPE, BWT_INDEX_SOCKET_NAME);
    if (bwtSocket == NULL) {
        fprintf(stderr, "BWTLoad(): Error creating socket!\n");
        exit(1);
    }

    // Accept connection
    for (;;) {
        if (SocketAcceptConnection(bwtSocket) < 0) {
            fprintf(stderr, "BWTLoad(): Error accepting connection!\n");
            exit(1);
        }
        if (SocketReceive(bwtSocket, &requestBwtIndexFileName, sizeof(BWTIndexFileName)) < 0) {
            fprintf(stderr, "BWTLoad(): Error receiving request!\n");
            exit(1);
        }
        if (requestBwtIndexFileName.BWTCodeFileName[0] == '\0') {
            BWTBlastFree(NULL, bwtIndexPointer.bwt, bwtIndexPointer.saIndexRange, bwtIndexPointer.saIndexRangeNumOfChar,
                         bwtIndexPointer.packedDNA, bwtIndexPointer.seqOffset, bwtIndexPointer.annotation, bwtIndexPointer.numOfSeq, 
                         bwtIndexPointer.ambiguity, bwtIndexPointer.numOfAmbiguity);
            printf("BWTLoad(): Resident processing unloaded.\n");
            exit(0);
        }

        // Check filenames in the request
        if (requestBwtIndexFileName.AnnotationFileName[0] == '-' || strncmp(requestBwtIndexFileName.AnnotationFileName, bwtIndexFileName.AnnotationFileName, MAX_FILENAME_LEN) == 0 ||
            requestBwtIndexFileName.AmbiguityFileName[0] == '-' || strncmp(requestBwtIndexFileName.AmbiguityFileName, bwtIndexFileName.AmbiguityFileName, MAX_FILENAME_LEN) == 0 ||
            requestBwtIndexFileName.PackedDNAFileName[0] == '-' || strncmp(requestBwtIndexFileName.PackedDNAFileName, bwtIndexFileName.PackedDNAFileName, MAX_FILENAME_LEN) == 0 ||
            requestBwtIndexFileName.BWTCodeFileName[0] == '-' || strncmp(requestBwtIndexFileName.BWTCodeFileName, bwtIndexFileName.BWTCodeFileName, MAX_FILENAME_LEN) == 0 ||
            requestBwtIndexFileName.SaValueFileName[0] == '-' || strncmp(requestBwtIndexFileName.SaValueFileName, bwtIndexFileName.SaValueFileName, MAX_FILENAME_LEN) == 0 ||
            requestBwtIndexFileName.SaIndexFileName[0] == '-' || strncmp(requestBwtIndexFileName.SaIndexFileName, bwtIndexFileName.SaIndexFileName, MAX_FILENAME_LEN) == 0) {
        } else {
            memcpy(&bwtIndexFileName, &requestBwtIndexFileName, sizeof(BWTIndexFileName));
            BWTBlastLoad(NULL, BWTCodeFileName, BWTOccValueFileName, SaValueFileName, SaIndexFileName,
                 PackedDNAFileName, AnnotationFileName, AmbiguityFileName,
                 &bwtIndexPointer.bwt, &bwtIndexPointer.saIndexRange, &bwtIndexPointer.saIndexRangeNumOfChar,
                 &bwtIndexPointer.packedDNA, &bwtIndexPointer.annotation, &bwtIndexPointer.seqOffset, &bwtIndexPointer.numOfSeq,
                 &bwtIndexPointer.ambiguity, &bwtIndexPointer.numOfAmbiguity, &bwtIndexPointer.minSeqLength);
        }
        if (SocketSend(bwtSocket, &bwtIndexPointer, sizeof(BWTIndexPointer)) < 0) {
            printf("BWTLoad(): Error replying with BWT index pointers!\n");
            exit(1);
        }
        if (SocketEndConnection(bwtSocket) < 0) {
            printf("BWTLoad(): Error ending connection!\n");
            exit(1);
        }
    }

}

dictionary *ParseInput(int argc, char** argv) {

    dictionary *programInput;
    char t1[3] = "-c";    // specify that this is a boolean type parameter; no following argument
    char t2[3] = "-x";    // specify that this is a boolean type parameter; no following argument
    char *d[2];

    d[0] = t1;
    d[1] = t2;

    programInput = paraparser_load(argc, argv, 2, d);    // 3 parameters are boolean type

    // Whether the user wants to stop the resident process
    UnloadBWTIndex = iniparser_find_entry(programInput, "parameter:-s");

    // Whether confirmation is needed
    Confirmation = iniparser_find_entry(programInput, "parameter:-c");

    if (!UnloadBWTIndex) {

        // Get database and query name
        if (!iniparser_find_entry(programInput, "parameter:-d")) {
            // Database name may be entered through argument
            if (!iniparser_find_entry(programInput, "argument:1")) {
                PrintHelp();
                exit(1);
            }
            iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
        } else {
            iniparser_copystring(programInput, "parameter:-d", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
        }

    }

    return programInput;

}

void ParseIniFile(char *iniFileName) {

    dictionary *ini;

    printf("Loading %s ..", iniFileName);
    ini = iniparser_load(iniFileName, FALSE);
    if (ini == NULL) {
        printf("not found.\n");
        return;
    }
    printf("done.\n");

    // Database parameters
    iniparser_copystring(ini, "Database:Location", DatabaseLocation, DatabaseLocation, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:SaIndexFileName", SaIndexFileName, SaIndexFileName, MAX_FILENAME_LEN);

    iniparser_freedict(ini);

}

void ProcessIni() {

    ProcessFileName(AnnotationFileName, AnnotationFileName, DatabaseLocation, DatabaseName);
    ProcessFileName(AmbiguityFileName, AmbiguityFileName, DatabaseLocation, DatabaseName);
    ProcessFileName(PackedDNAFileName, PackedDNAFileName, DatabaseLocation, DatabaseName);
    ProcessFileName(BWTCodeFileName, BWTCodeFileName, DatabaseLocation, DatabaseName);
    ProcessFileName(BWTOccValueFileName, BWTOccValueFileName, DatabaseLocation, DatabaseName);
    ProcessFileName(SaValueFileName, SaValueFileName, DatabaseLocation, DatabaseName);
    ProcessFileName(SaIndexFileName, SaIndexFileName, DatabaseLocation, DatabaseName);

}

void ValidateIni() {

    if (AnnotationFileName[0] == '\0') {
        fprintf(stderr, "Annotation file is not specified!\n");
        exit(1);
    }
    if (AmbiguityFileName[0] == '\0') {
        fprintf(stderr, "Ambiguity file is not specified!\n");
        exit(1);
    }    
    if (PackedDNAFileName[0] == '\0') {
        fprintf(stderr, "Packed DNA file is not specified!\n");
        exit(1);
    }
    if (BWTCodeFileName[0] == '\0') {
        fprintf(stderr, "BWT Code file is not specified!\n");
        exit(1);
    }
    if (BWTOccValueFileName[0] == '\0') {
        fprintf(stderr, "BWT Occ value file is not specified!\n");
        exit(1);
    }
    if (SaValueFileName[0] == '\0') {
        fprintf(stderr, "SA value file is not specified!\n");
        exit(1);
    }
    if (SaIndexFileName[0] == '\0') {
        fprintf(stderr, "SA index file is not specified!\n");
        exit(1);
    }

}


void PrintIni() {

    char boolean[2];

    boolean[0] = 'N';
    boolean[1] = 'Y';

    printf("Annotation file          : %s\n", AnnotationFileName);
    printf("Ambiguity file           : %s\n", AmbiguityFileName);
    printf("Packed DNA file          : %s\n", PackedDNAFileName);
    printf("BWT Code file            : %s\n", BWTCodeFileName);
    printf("BWT Occ value file       : %s\n", BWTOccValueFileName);
    printf("SA value file            : %s\n", SaValueFileName);
    printf("SA index file            : %s\n", SaIndexFileName);
    printf("\n");

}

void PrintShortDesc() {

    printf("BWTLoad v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
    printf("BWTLoad comes with ABSOLUTELY NO WARRENTY.\n");
    printf("BWTLoad is free software, and you are welcome to\n");
    printf("redistribute it under certain conditions.\n");
    printf("For details type BWTBlast.\n");
    printf("\n");

}

void PrintHelp() {

    printf("BWTLoad v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
    printf("\n");

    printf("This program is free software; you can redistribute it and/or\n");
    printf("modify it under the terms of the GNU General Public License\n");
    printf("as published by the Free Software Foundation; either version 2\n");
    printf("of the License, or (at your option) any later version.\n");
    printf("\n");

    printf("This program is distributed in the hope that it will be useful,\n");
    printf("but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
    printf("GNU General Public License for more details.\n");
    printf("\n");

    printf("You should have received a copy of the GNU General Public License\n");
    printf("along with this program; if not, write to the Free Software\n");
    printf("Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n");
    printf("\n");

    printf("Syntax: BWTLoad <database> [-c confirm] [-x exit/unload]\n");

}

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseLocation,
                     const char *databaseName) {

    char tempChar[MAX_FILENAME_LEN];
    unsigned int i;

    if (inputFileName == NULL || inputFileName[0] == '\0' || inputFileName[0] == ' ' || inputFileName[0] == '-') {
        if (outputFileName != inputFileName) {
            outputFileName[0] = '-';
        }
        return;
    }

    if (strlen(databaseLocation) + strlen(databaseName) + strlen(inputFileName) > MAX_FILENAME_LEN) {
        fprintf(stderr, "File length is too long!\n");
        exit(1);
    }

    strncpy(tempChar, inputFileName, MAX_FILENAME_LEN);

    // locate the *
    for (i=0; i<MAX_FILENAME_LEN; i++) {
        if (tempChar[i] == '*') {
            break;
        }
    }
    if (i<MAX_FILENAME_LEN) {
        tempChar[i] = '\0';
        sprintf(outputFileName, "%s%s%s%s", databaseLocation, tempChar, databaseName, tempChar + i + 1);
    } else {
        sprintf(outputFileName, "%s%s", databaseLocation, tempChar);
    }

}

