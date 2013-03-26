//
//    2BWT-Builder.c
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
/*

   BWTFormatdb.c        Build index for FASTA database

   This program builds index for FASTA database for use of BWTBlastn.

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


   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Thus, changing all references to 2bwt lib to subdirectory.

   Date   : 23rd October 2011
   Author : Edward MK Wu
   Change : Fix a rounding error when building reverse packed sequence.
   
*/

#include <stdio.h>
#include <stdlib.h>

#include "../2bwt-lib/TypeNLimit.h"
#include "../2bwt-lib/BWTConstruct.h"
#include "../2bwt-lib/MiscUtilities.h"
#include "../2bwt-lib/DNACount.h"
#include "../2bwt-lib/TextConverter.h"
#include "../2bwt-lib/MemManager.h"
#include "../2bwt-lib/iniparser.h"
#include "../2bwt-lib/HSP.h"
#include "../2bwt-lib/Timing.h"

#include "LTConstruct.h"
#include "HOCCConstruct.h"

#include "../Release.h"


// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();
void PrintIni();
void PrintShortDesc();
void PrintHelp();

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName);

    // Parameters
    char IniFileName[MAX_FILENAME_LEN+1];
    int Confirmation;
    
    // BuildTasks parameters
    int ParseFASTA = TRUE;
    int BuildBWT = TRUE;
    int BuildSaValue = TRUE;
    int BuildSaIndex = FALSE;
    int BuildLookUp = TRUE;
    int BuildHOT = FALSE;

    // Memory parameters
    unsigned int PoolSize = 2097152;                // 2M  - fixed; not configurable through ini

    // Display parameters
    int ShowProgress = FALSE;

    // Database parameters
    char FASTAFileName[MAX_FILENAME_LEN+1] = "";
    char DatabaseName[MAX_FILENAME_LEN+1] = "";
    char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.index.ann";
    char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.index.amb";
    char TranslateFileName[MAX_FILENAME_LEN+1] = "*.index.tra";
    char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.index.pac";
    char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.index.bwt";
    char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.index.fmv";
    char SaValueFileName[MAX_FILENAME_LEN+1] = "*.index.sa";
    char SaIndexFileName[MAX_FILENAME_LEN+1] = "*.index.sai";
    
    char RevPackedDNAFileName[MAX_FILENAME_LEN+1] = "*.index.rev.pac";
    char RevBWTCodeFileName[MAX_FILENAME_LEN+1] = "*.index.rev.bwt";
    char RevBWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.index.rev.fmv";
    
    char LookupTableFileName[MAX_FILENAME_LEN+1] = "*.index.lkt";
    char RevLookupTableFileName[MAX_FILENAME_LEN+1] = "*.index.rev.lkt";

    char HighOccPatternFileName[MAX_FILENAME_LEN+1] = "*.index.ht";
    char HighOccFileName[MAX_FILENAME_LEN+1] = "*.index.hocc";
    char HighOccPackedFileName[MAX_FILENAME_LEN+1] = "*.index.hpac";
    char HighOccAuxValueFileName[MAX_FILENAME_LEN+1] = "*.index.hv";
    // Parse FASTA parameters
    unsigned int FASTARandomSeed = 0;
    int MaskLowerCase = FALSE;

    // Build BWT parameters
    unsigned int OccValueFreq = 256;
    float TargetNBit = 5;
    unsigned int InitialMaxBuildSize = 10000000;
    unsigned int IncMaxBuildSize = 10000000;

    // Build SA value parameters
    unsigned int SaValueFreq = 1;

    // Build SA index parameters
    unsigned int SaIndexNumOfChar = 12;
    
    //Look Up Table parameters
    unsigned int LookUpTableSize=13;
    
    //High Occurrences Pattern Hash Table parameters
    unsigned int HashPatternLength=35;
    unsigned int HashOccThreshold=4;


void printBinary(unsigned long long seq,int len) {
    char text[64];
    int i,j=63;
    for (i=0;i<len;i++) {
        text[j--]=seq %2;
        seq>>=1;
    }
    for (i=j+1;i<64;i++) {
        if (text[i]==0) printf("0");
        if (text[i]==1) printf("1");
        if ((i-j-1) % 4 ==3) printf(" "); 
    }printf("\n");
}

void BuildReversePacked(const char *inputFileName, unsigned int *textLength, const unsigned int convertToWordPacked, const unsigned int trailerBufferInWord) {

    FILE *inputFile;
    FILE *outputFile;
    unsigned char tempChar[4];
    unsigned int writeBuffer=0;
    unsigned int writeBufferIndex=0;
    unsigned int *packedText;
    unsigned int packedFileLen;
    unsigned char lastByteLength;
    unsigned int wordToProcess;
    long long i;
    int j,k;
    unsigned int trailerBufferIn128;

    trailerBufferIn128 = (trailerBufferInWord + 3) / 4 * 4;

    inputFile = (FILE*)(FILE*)fopen64(inputFileName, "rb");
    outputFile = (FILE*)(FILE*)fopen64(RevPackedDNAFileName, "wb");

    if (inputFile == NULL) {
        fprintf(stderr, "BuildReversePacked() : Cannot open inputFileName!\n");
        exit(1);
    }

    fseek(inputFile, -1, SEEK_END);
    packedFileLen = ftell(inputFile);
    if ((int)packedFileLen < 0) {
        fprintf(stderr, "BuildReversePacked(): Cannot determine file length!\n");
        exit(1);
    }
    fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);

    *textLength = (packedFileLen - 1) * 4 + lastByteLength;

    if (ShowProgress) {
        printf("Packed file size = %u\n",(unsigned int) packedFileLen);
        printf("Text Length = %u\n",*textLength);
    }

    wordToProcess = (*textLength + 64 - 1) / 64 * 4 + trailerBufferIn128 * 4;        // allocate multiple of 128 bit + trailer buffer
    
    packedText = (unsigned int*) MMUnitAllocate(wordToProcess * sizeof(unsigned int));    // allocate 1 more word at end
    for (i=(*textLength)/16; i<wordToProcess; i++) {
        packedText[i] = 0;
    }

    fseek(inputFile, 0, SEEK_SET);
    fread(packedText, 1, packedFileLen, inputFile);
    fclose(inputFile);
    
    
    
    //printf("lastByteLength = %u\n",lastByteLength);
    long long currentLocation = (wordToProcess)*16;
    //printf("currentLocation = %u\n",currentLocation);
    //printf("*textLength = %u\n",*textLength);
    for (i=wordToProcess-1; i>=0; i--) {
        *(unsigned int*)tempChar = packedText[i];
        packedText[i] = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];

        unsigned int mk = packedText[i];
        for (j=0;j<16;j++) {
            if (writeBufferIndex>=16) {
                *(unsigned int*)tempChar = writeBuffer;
                writeBuffer = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];
                fwrite(&writeBuffer,sizeof(unsigned int),1,outputFile);
                writeBufferIndex=0;
            }
            writeBuffer<<=2;
            unsigned char c = mk & 3;
            mk>>=2;
            currentLocation--;
            if (currentLocation>=*textLength) {continue;}
            writeBuffer|=c;writeBufferIndex++;
            
        }
    }
    //*/
    
    
    
    //printf("Finished main loop..\n");
    if (writeBufferIndex>0) {
        //printf("Wiping..%d\n",writeBufferIndex);
        for (k=writeBufferIndex;k<16;k++) writeBuffer<<=2;
        *(unsigned int*)tempChar = writeBuffer;
        if (writeBufferIndex<=4) {
            fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
        } else if (writeBufferIndex<=8) {
            fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
            fwrite(&tempChar[2],sizeof(unsigned char),1,outputFile);
        } else if (writeBufferIndex<=12) {
            fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
            fwrite(&tempChar[2],sizeof(unsigned char),1,outputFile);
            fwrite(&tempChar[1],sizeof(unsigned char),1,outputFile);
        } else if (writeBufferIndex<=16) {
            fwrite(&tempChar[3],sizeof(unsigned char),1,outputFile);
            fwrite(&tempChar[2],sizeof(unsigned char),1,outputFile);
            fwrite(&tempChar[1],sizeof(unsigned char),1,outputFile);
            fwrite(&tempChar[0],sizeof(unsigned char),1,outputFile);
        }
    }
    if (writeBufferIndex % 4 == 0) {
        unsigned char c = 0;
        fwrite(&c, 1, 1, outputFile);
    }
    fwrite(&lastByteLength,sizeof(unsigned char),1,outputFile);
    fclose(outputFile);
    free(packedText);
    return;
}

int main(int argc, char** argv) {
    

    char c;
    MMPool *mmPool;
    dictionary *programInput;
    double startTime;
    double elapsedTime = 0, totalElapsedTime = 0;

    char filename[MAX_FILENAME_LEN+1];
    BWT *bwt = NULL;
    BWT *rev_bwt = NULL;
    HSP *hsp = NULL;
    unsigned int textLength = 0;
    unsigned int numSeq;

    BWTInc *bwtInc = NULL;
    BWTInc *rev_bwtInc = NULL;

    // Program input
    programInput = ParseInput(argc, argv);
    PrintShortDesc();

    // Ini
    if (strcmp(argv[0] + strlen(argv[0]) - 4, ".exe") == 0) {
        *(argv[0] + strlen(argv[0]) - 4) = '\0';
    }
    sprintf(filename, "%s.ini", argv[0]);
    ParseIniFile(filename);
    //printf("\n");
    ProcessIni();
    ValidateIni();
    PrintIni();

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

    startTime = setStartTime();

    MMMasterInitialize(1, 0, FALSE, NULL);
    mmPool = MMPoolCreate(PoolSize);

    // Parse FASTA file to produce packed DNA and annotation file
    if (ParseFASTA == TRUE) {

        printf("Parsing FASTA file..\n");
        numSeq = HSPParseFASTAToPacked(FASTAFileName, AnnotationFileName, PackedDNAFileName, AmbiguityFileName, TranslateFileName, FASTARandomSeed, MaskLowerCase);
        printf("Finished. Parsed %u sequences.\n", numSeq);
        
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");
        
        //Parse packed DNA to construct the packed reversed DNA
        
        printf("Parsing FASTA file reverse..\n");
        unsigned int textLen;
        BuildReversePacked(PackedDNAFileName,&textLen,TRUE, 1);
        //printf("Reversed Packed DNA generated..\n");

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

    }
    
    if (BuildLookUp == TRUE) {
        //Construct Look-Up Table
        printf("Building Look-Up..\n");
        BuildLookupTable(PackedDNAFileName,LookupTableFileName,LookUpTableSize);
        //printf("Look-Up Table is built.\n");
        BuildLookupTable(RevPackedDNAFileName,RevLookupTableFileName,LookUpTableSize);
        //printf("Reversed Look-Up Table is built.\n");
        
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Finished.\nElapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");
    }

    // Construct BWTInc from text
    if (BuildBWT == TRUE) {

        printf("Building BWT..\n");

        bwtInc = BWTIncConstructFromPacked(mmPool, PackedDNAFileName, ShowProgress, 
                                           TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

        printf("Finished constructing BWT in %u iterations.  ", bwtInc->numberOfIterationDone);
        
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        printf("Saving BWT..\n");
        BWTSaveBwtCodeAndOcc(bwtInc->bwt, BWTCodeFileName, BWTOccValueFileName);
        printf("Finished saving BWT.  ");
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        textLength = bwtInc->bwt->textLength;


        //Building Reversed BWT
        printf("Building Reversed BWT..\n");

        rev_bwtInc = BWTIncConstructFromPacked(mmPool, RevPackedDNAFileName, ShowProgress, 
                                           TargetNBit, InitialMaxBuildSize, IncMaxBuildSize);

        printf("Finished constructing Reversed BWT in %u iterations.  ", rev_bwtInc->numberOfIterationDone);
        
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        printf("Saving BWT..\n");
        BWTSaveBwtCodeAndOcc(rev_bwtInc->bwt, RevBWTCodeFileName, RevBWTOccValueFileName);
        printf("Finished saving BWT.  ");
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        textLength = rev_bwtInc->bwt->textLength;

        BWTIncFree(mmPool, bwtInc);
        BWTIncFree(mmPool, rev_bwtInc);
    }

    // Load BWT
    if (BuildSaValue || BuildSaIndex) {

        printf("Loading BWT...\n");

        bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, NULL, NULL, NULL, NULL, 0);
        //Use BWT to build the hash table

        printf("Finished loading BWT.  ");

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        textLength = bwt->textLength;

    }
    

    if (BuildSaValue) {

        printf("Building SA value...\n");
        
        if (ShowProgress) {
            BWTGenerateSaValue(bwt, SaValueFreq, bwt->textLength / SaValueFreq / 10);
        } else {
            BWTGenerateSaValue(bwt, SaValueFreq, 0);
        }
        BWTSaveSaValue(bwt, SaValueFileName);

        printf("Finished building SA value.  ");

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

    }
    

    /*if (BuildSaIndex) {

        printf("Building SA index...\n");
        
        BWTGenerateCachedSaIndex(bwt, SaIndexNumOfChar, SaIndexFileName);

        printf("Finished building SA index.  ");

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

    }*/

    // Free BWT
    if (BuildSaValue || BuildSaIndex) {
        BWTFree(mmPool, bwt, 0);
    }
    
    
    
    
    
    
    
    
    
    if (BuildHOT) {
        printf("Loading BWT and SA Values...\n");

        bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SaValueFileName, NULL, NULL, NULL, 0);
        rev_bwt = BWTLoad(mmPool, RevBWTCodeFileName, RevBWTOccValueFileName, NULL, NULL, NULL, NULL, 0);
        hsp = HSPLoad(mmPool, PackedDNAFileName, AnnotationFileName, AmbiguityFileName, TranslateFileName, 1, 0);
        //Use BWT to build the hash table

        printf("Finished loading BWT and SA Values.  ");

        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        textLength = bwt->textLength;

    }
    //Build High Occurrences Pattern Hash Table
    if (BuildHOT) {
        printf("Building High-Occ Hash Table...\n");
        HOCCBuilder * hoccBuilder = HOCCInitialHighOccWorkGround();
        
        unsigned int l,r;
        unsigned char ec;
        for (ec=0;ec<4;ec++) {
            l = bwt->cumulativeFreq[ec]+1;
            r = bwt->cumulativeFreq[ec+1];
            HOCCExtractionHighOccPattern(hoccBuilder, bwt,rev_bwt,HashPatternLength-1,HashOccThreshold,l,r,l,r);
        }
        
        HOCCPopulateHighOccPatternToBitPattern(hoccBuilder, bwt,hsp,HighOccPatternFileName,HighOccPackedFileName,HighOccAuxValueFileName,HighOccFileName);
        
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Finished.\nElapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");
        HOCCFreeHighOccWorkGround(hoccBuilder);
    }
    if (BuildHOT) {
        
        BWTFree(mmPool, bwt, 0);
        
        BWTFree(mmPool, rev_bwt, 0);
        
        HSPFree(mmPool, hsp, 1, 0);
    }
        

    // Finished all construction tasks
    printf("Index building is completed.\n");
    totalElapsedTime = getElapsedTime(startTime);
    printf("Total elapsed time = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, totalElapsedTime);
    printf("\n");

    //MMMasterPrintReport(stdout, FALSE, FALSE, FALSE);
    if (BuildSaValue) {
        //fprintf(stdout, "Number of char   :  %u\n", textLength);
        //fprintf(stdout, "Bit per char     :  %.2f\n", (float)MMMasterMaxTotalByteDispatched() * BITS_IN_BYTE / textLength);
        //printf("\n");
    }

    MMPoolFree(mmPool);

    iniparser_freedict(programInput);

    return 0;

}

dictionary *ParseInput(int argc, char** argv) {

    dictionary *programInput;
    char t1[3] = "-c";    // specify that this is a boolean type parameter
    char t2[3] = "-U";    // specify that this is a boolean type parameter
    char *d[2];

    d[0] = t1;
    d[1] = t2;
    
    programInput = paraparser_load(argc, argv, 2, d);    // 2 boolean type parameters

    // Get database name
    if (!iniparser_find_entry(programInput, "argument:1")) {
        PrintHelp();
        exit(1);
    }
    iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
    if (strlen(DatabaseName) + 4 > MAX_FILENAME_LEN) {
        PrintHelp();
        exit(1);
    }

    // Get FASTA file name
    iniparser_copystring(programInput, "argument:2", FASTAFileName, DatabaseName, MAX_FILENAME_LEN);
    if (strlen(FASTAFileName) > MAX_FILENAME_LEN) {
        PrintHelp();
        exit(1);
    }


    // Whether confirmation is needed
    Confirmation = iniparser_find_entry(programInput, "parameter:-c");

    MaskLowerCase = iniparser_find_entry(programInput, "parameter:-U");

    return programInput;

}

void ParseIniFile(char *iniFileName) {

    dictionary *ini;

    //printf("Loading %s ..", iniFileName);
    ini = iniparser_load(iniFileName, FALSE);
    if (ini == NULL) {
        printf("[WARNING] Expecting ini file %s but not found. Default values in builder is used.\n",iniFileName);
        printf("[WARNING] For best alignment performance, using builder with ini shipped with %s is advised.\n",PROJECT_NAME);
        return;
    }
    //printf("done.\n");

    // BuildTasks parameters
    ParseFASTA = iniparser_getboolean(ini, "BuildTasks:ParseFASTA", ParseFASTA);
    BuildBWT = iniparser_getboolean(ini, "BuildTasks:BuildBWT", BuildBWT);
    BuildSaValue = iniparser_getboolean(ini, "BuildTasks:BuildSaValue", BuildSaValue);
    BuildLookUp = iniparser_getboolean(ini, "BuildTasks:BuildLookUp", BuildLookUp);
    BuildSaIndex = iniparser_getboolean(ini, "BuildTasks:BuildSaIndex", BuildSaIndex);
    BuildHOT = iniparser_getboolean(ini, "BuildTasks:BuildHOT", BuildHOT);

    // Display parameters
    ShowProgress = iniparser_getboolean(ini, "Display:ShowProgress", ShowProgress);

    // Parse FASTA parameters
    FASTARandomSeed = iniparser_getint(ini, "ParseFASTA:RandomSeed", FASTARandomSeed);
    if (FASTARandomSeed == 0) {
        FASTARandomSeed = getRandomSeed();
    }

    // Build BWT parameters
    OccValueFreq = iniparser_getint(ini, "BuildBWT:OccValueFreq", OccValueFreq);
    TargetNBit = (float)iniparser_getdouble(ini, "BuildBWT:TargetNBit", TargetNBit);
    InitialMaxBuildSize = iniparser_getint(ini, "BuildBWT:InitialMaxBuildSize", InitialMaxBuildSize);
    IncMaxBuildSize = iniparser_getint(ini, "BuildBWT:IncMaxBuildSize", IncMaxBuildSize);

    // Build SA value parameters
    SaValueFreq = iniparser_getint(ini, "BuildSAValue:SaValueFreq", SaValueFreq);

    // Build SA index parameters
    SaIndexNumOfChar = iniparser_getint(ini, "BuildSAIndex:SaIndexNumOfChar", SaIndexNumOfChar);

    // Build Look Up parameters
    LookUpTableSize = iniparser_getint(ini, "BuildLookUp:TableSize", LookUpTableSize);
    

    // Build High Occurrences Pattern Hash Table parameters
    HashPatternLength = iniparser_getint(ini, "BuildHOT:PatternLength", HashPatternLength);
    HashOccThreshold = iniparser_getint(ini, "BuildHOT:OccurrencesThreshold", HashOccThreshold);
    // Database parameters
    iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:TranslateFileName", TranslateFileName, TranslateFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:SaIndexFileName", SaIndexFileName, SaIndexFileName, MAX_FILENAME_LEN);
    
    iniparser_copystring(ini, "Database:RevPackedDNAFileName", RevPackedDNAFileName, RevPackedDNAFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:RevBWTCodeFileName", RevBWTCodeFileName, RevBWTCodeFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:RevBWTOccValueFileName", RevBWTOccValueFileName, RevBWTOccValueFileName, MAX_FILENAME_LEN);
    
    iniparser_copystring(ini, "Database:LookupTableFileName", LookupTableFileName, LookupTableFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:RevLookupTableFileName", RevLookupTableFileName, RevLookupTableFileName, MAX_FILENAME_LEN);
    
    iniparser_copystring(ini, "Database:HighOccPatternFileName", HighOccPatternFileName, HighOccPatternFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:HighOccFileName", HighOccFileName, HighOccFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:HighOccPackedFileName", HighOccPackedFileName, HighOccPackedFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:HighOccAuxValueFileName", HighOccAuxValueFileName, HighOccAuxValueFileName, MAX_FILENAME_LEN);
    
    iniparser_copystring(ini, "Database:TranslateFileName", TranslateFileName, TranslateFileName, MAX_FILENAME_LEN);
    
    iniparser_freedict(ini);

}

void ProcessIni() {

    ProcessFileName(AnnotationFileName, AnnotationFileName, DatabaseName);
    ProcessFileName(AmbiguityFileName, AmbiguityFileName, DatabaseName);
    ProcessFileName(TranslateFileName, TranslateFileName, DatabaseName);
    ProcessFileName(PackedDNAFileName, PackedDNAFileName, DatabaseName);
    ProcessFileName(RevPackedDNAFileName, RevPackedDNAFileName, DatabaseName);
    ProcessFileName(BWTCodeFileName, BWTCodeFileName, DatabaseName);
    ProcessFileName(RevBWTCodeFileName, RevBWTCodeFileName, DatabaseName);
    ProcessFileName(BWTOccValueFileName, BWTOccValueFileName, DatabaseName);
    ProcessFileName(RevBWTOccValueFileName, RevBWTOccValueFileName, DatabaseName);
    ProcessFileName(SaValueFileName, SaValueFileName, DatabaseName);
    ProcessFileName(SaIndexFileName, SaIndexFileName, DatabaseName);
    
    ProcessFileName(LookupTableFileName, LookupTableFileName, DatabaseName); 
    ProcessFileName(RevLookupTableFileName, RevLookupTableFileName, DatabaseName);
    
    ProcessFileName(HighOccPatternFileName, HighOccPatternFileName, DatabaseName);
    ProcessFileName(HighOccFileName, HighOccFileName, DatabaseName);
    ProcessFileName(HighOccPackedFileName, HighOccPackedFileName, DatabaseName);
    ProcessFileName(HighOccAuxValueFileName, HighOccAuxValueFileName, DatabaseName);

    ProcessFileName(TranslateFileName, TranslateFileName, DatabaseName);
}

void ValidateIni() {

    if (!ParseFASTA && !BuildBWT && !BuildSaValue && !BuildSaIndex && !BuildLookUp && !BuildHOT) {
        fprintf(stderr, "No action is specified!\n");
        exit(1);
    }
    if (BuildLookUp) {
        if (PackedDNAFileName[0] == '\0') {
            fprintf(stderr, "Packed DNA file name is not specified!\n");
            exit(1);
        }
        if (RevPackedDNAFileName[0] == '\0') {
            fprintf(stderr, "Reversed Packed DNA file name is not specified!\n");
            exit(1);
        }
    }
    if (ParseFASTA) {
        if (PackedDNAFileName[0] == '\0') {
            fprintf(stderr, "Packed DNA file name is not specified!\n");
            exit(1);
        }
        if (AnnotationFileName[0] == '\0') {
            fprintf(stderr, "Annotation file name is not specified!\n");
            exit(1);
        }
        if (AmbiguityFileName[0] == '\0') {
            fprintf(stderr, "Ambiguity file name is not specified!\n");
            exit(1);
        }
    }
    if (BuildBWT) {
        if (PackedDNAFileName[0] == '\0') {
            fprintf(stderr, "Packed DNA file is not specified!\n");
            exit(1);
        }
        if (BWTCodeFileName[0] == '\0') {
            fprintf(stderr, "BWT code file name is not specified!\n");
            exit(1);
        }
        if (BWTOccValueFileName[0] == '\0') {
            fprintf(stderr, "BWT Occ value file name is not specified!\n");
            exit(1);
        }
        if (TargetNBit < 2.5) {
            fprintf(stderr, "Target NBit should be at least 2.5!\n");
            exit(1);
        }
    }
    if (BuildSaValue) {
        if (BWTCodeFileName[0] == '\0') {
            fprintf(stderr, "BWT code file is not specified!\n");
            exit(1);
        }
        if (BWTOccValueFileName[0] == '\0') {
            fprintf(stderr, "BWT Occ value file is not specified!\n");
            exit(1);
        }
        if (SaValueFileName[0] == '\0') {
            fprintf(stderr, "SA value file name is not specified!\n");
            exit(1);
        }
        if (SaValueFreq <= 0) {
            fprintf(stderr, "SA value frequency must > 0!\n");
            exit(1);
        }
    }

    if (BuildSaIndex) {
        if (BWTCodeFileName[0] == '\0') {
            fprintf(stderr, "BWT code file is not specified!\n");
            exit(1);
        }
        if (BWTOccValueFileName[0] == '\0') {
            fprintf(stderr, "BWT Occ value file is not specified!\n");
            exit(1);
        }
        if (SaIndexFileName[0] == '\0') {
            fprintf(stderr, "SA index file name is not specified!\n");
            exit(1);
        }
        if (SaIndexNumOfChar <= 0) {
            fprintf(stderr, "SA index number of character must > 0!\n");
            exit(1);
        }
        if (SaIndexNumOfChar > 13) {
            fprintf(stderr, "SA index number of character must <= 13!\n");
            exit(1);
        }
    }

    if (BuildHOT) {
        if (BWTCodeFileName[0] == '\0') {
            fprintf(stderr, "BWT code file is not specified!\n");
            exit(1);
        }
        if (BWTOccValueFileName[0] == '\0') {
            fprintf(stderr, "BWT Occ value file is not specified!\n");
            exit(1);
        }
        if (RevBWTCodeFileName[0] == '\0') {
            fprintf(stderr, "Reversed BWT code file is not specified!\n");
            exit(1);
        }
        if (RevBWTOccValueFileName[0] == '\0') {
            fprintf(stderr, "Reversed BWT Occ value file is not specified!\n");
            exit(1);
        }
    }
}


void PrintIni() {

    char boolean[2];

    boolean[0] = 'N';
    boolean[1] = 'Y';

    /*printf("Parse FASTA file    : %c\n", boolean[ParseFASTA]);
    printf("Build BWT           : %c\n", boolean[BuildBWT]);
    printf("Build SA value      : %c\n", boolean[BuildSaValue]);
    printf("Build SA index      : %c\n", boolean[BuildSaIndex]);
    printf("\n");

    printf("Show progress       : %c\n", boolean[ShowProgress]);
    printf("\n");

    if (ParseFASTA) {
        printf("Parse FASTA :\n");
        printf("Mask lower case         : %c\n", boolean[MaskLowerCase]);
        printf("Random seed             : %u\n", FASTARandomSeed);
        printf("\n");
    }

    if (BuildBWT) {
        printf("Build BWT :\n");
        printf("Target N Bits           : %.2f\n", TargetNBit);
        printf("Occ value frequency     : %u\n", OccValueFreq);
        printf("Initial Max Build Size  : %u    Inc Max Build Size : %u\n", 
                InitialMaxBuildSize, IncMaxBuildSize);
        printf("\n");
    }

    if (BuildSaValue) {
        printf("Build SA value :\n");
        printf("SA value frequency      : %u\n", SaValueFreq);
        printf("\n");
    }

    if (BuildSaIndex) {
        printf("Build SA index :\n");
        printf("SA index no. of char    : %u\n", SaIndexNumOfChar);
        printf("\n");
    }

    printf("Annotation file          : %s\n", AnnotationFileName);
    printf("Ambigurity file          : %s\n", AmbiguityFileName);
    printf("Packed DNA file          : %s\n", PackedDNAFileName);
    printf("BWT Code file            : %s\n", BWTCodeFileName);
    printf("BWT Occ value file       : %s\n", BWTOccValueFileName);
    printf("SA value file            : %s\n", SaValueFileName);
    printf("SA index file            : %s\n", SaIndexFileName);
    printf("\n");
    printf("Reversed Packed DNA file          : %s\n", RevPackedDNAFileName);
    printf("Reversed BWT Code file            : %s\n", RevBWTCodeFileName);
    printf("Reversed BWT Occ value file       : %s\n", RevBWTOccValueFileName);
    printf("\n");
    printf("Look-Up Table file                : %s\n", LookupTableFileName);
    printf("Reversed Look-Up Table file       : %s\n", RevLookupTableFileName);
    printf("\n");*/

}

void PrintShortDesc() {

    /*printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
    printf("BWTFormatdb comes with ABSOLUTELY NO WARRENTY.\n");
    printf("BWTFormatdb is free software, and you are welcome to\n");
    printf("redistribute it under certain conditions.\n");
    printf("For details type BWTFormatdb.\n");
    printf("\n");*/

}

void PrintHelp() {

    /*printf("BWTFormatdb v1.0, Copyright (C) 2006, Wong Chi Kwong.\n");
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
    printf("\n");*/

    printf("\n[Main] %s Builder v%d.%d.%d (%s) - Usage guide:\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
    printf("\n");
    printf("Syntax: %s <sequence file>\n",PROJECT_BUILDER_BINARY);

}

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName) {

    char tempChar[MAX_FILENAME_LEN];
    unsigned int i;

    if (inputFileName == NULL) {
        if (outputFileName != inputFileName) {
            outputFileName[0] = '\0';
        }
        return;
    }

    if (strlen(databaseName) + strlen(inputFileName) > MAX_FILENAME_LEN) {
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
        sprintf(outputFileName, "%s%s%s", tempChar, databaseName, tempChar + i + 1);
    } else {
        sprintf(outputFileName, "%s", tempChar);
    }

}

