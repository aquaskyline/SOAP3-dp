/*

   SearchBWT.c        Search BWT - Search pattern by BWT-index

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
#include "MiscUtilities.h"
#include "Timing.h"
#include "MemManager.h"
#include "BWT.h"
#include "TextConverter.h"
#include "iniparser.h"

dictionary *ParseInput(int argc, char** argv);
dictionary *ParseIniFile();
void ProcessIni();
void ValidateIni();
void PrintIni();
void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseLocation, const char *databaseName);
unsigned int ProcessSearchPattern();


#define MAX_SEARCH_PATTERN_LENGTH    65536

    // Parameters

    char DatabaseName[MAX_FILENAME_LEN+1] = "";
    char IniFileName[MAX_FILENAME_LEN+1] = "";
    char PatternFileName[MAX_FILENAME_LEN+1] = "";
    char LogFileName[MAX_FILENAME_LEN+1] = "";
    unsigned int NumberOfSearchPatternFile;

    int Confirmation;
    
    // Action parameters
    int SABinarySearch;
    int BackwardSearch, BackwardSearchCheck;
    int HammingDistSearch, EditDistSearch;
    int SubPatternHammingDistSearch, SubPatternEditDistSearch;
    int FindTextPosition;
    unsigned int MaxNumberOfTextPosition, TextCheckCostFactor;
    unsigned int MaxNumberOfHitGroups;

    // Memory parameters
    int PoolSize;

    // Search pattern parameters
    unsigned int MaxErrorAllowed;

    // Database parameters
    char DatabaseLocation[MAX_FILENAME_LEN+1] = "./";

    // DatabaseFiles parameters
    char AnnotationFileName[MAX_FILENAME_LEN+1] = "*.ann";
    char AmbiguityFileName[MAX_FILENAME_LEN+1] = "*.amb";
    char PackedDNAFileName[MAX_FILENAME_LEN+1] = "*.pac";
    char BWTCodeFileName[MAX_FILENAME_LEN+1] = "*.bwt";
    char BWTOccValueFileName[MAX_FILENAME_LEN+1] = "*.fmv";
    char SaValueFileName[MAX_FILENAME_LEN+1] = "*.sa";
    char SaIndexFileName[MAX_FILENAME_LEN+1] = "*.sai";

    // File location parameters
    char *fileName;

    // Search pattern 
    unsigned char *searchPattern = NULL;
    unsigned int searchPatternAllocated = 0;
    unsigned int *searchPatternPosition = NULL;
    unsigned int searchPatternPositionAllocated = 0;

int main(int argc, char** argv) {

    char c;
    unsigned int i, j, k;
    MMPool *mmPool;
    dictionary *programInput, *ini;
    char argumentText[11] = "argument:0";
    double startTime;
    double elapsedTime = 0, totalElapsedTime = 0;

    BWT *bwt;
    HSP *hsp;

    unsigned char charMap[256];

    unsigned char *convertedKey;
    unsigned int *packedKey;

    unsigned int found;
    unsigned int numberOfPattern, numberOfPatternFound;
    unsigned int saIndexLeft, saIndexRight;
    SaIndexGroupNew *saIndexGroup;
    SaIndexList *tempSaIndex1, *tempSaIndex2;
    unsigned int numberOfSaIndexGroup;
    HitList *hitList;
    unsigned int textPositionMatched, textPositionRetrieved;
    unsigned long long totalTextPositionMatched, totalTextPositionRetrieved;
    BWTSaRetrievalStatistics BWTSaRetrievalStatistics;

    FILE *logFile;
    

    init(textPositionRetrieved);    // to avoid compiler warning only
    init(textPositionMatched);        // to avoid compiler warning only

    programInput = ParseInput(argc, argv);
    ini = ParseIniFile();

    ProcessIni();
    ValidateIni();

    PrintIni();

    if (Confirmation == TRUE) {
        printf("BWT Search. Press Y to go or N to cancel. ");
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

    // Load Database
    printf("Loading Database..");
    fflush(stdout);

    bwt = BWTLoad(mmPool, BWTCodeFileName, BWTOccValueFileName, SaValueFileName, NULL, SaIndexFileName, NULL);
    HSPFillCharMap(charMap);
    hsp = HSPLoad(mmPool, PackedDNAFileName, AnnotationFileName, AmbiguityFileName, MAX_SEARCH_PATTERN_LENGTH / CHAR_PER_WORD);
    if (bwt->textLength != hsp->dnaLength) {
        fprintf(stderr, "BWT-BLAST: Database length inconsistent!\n");
        exit(1);
    }


    if (LogFileName[0] != '\0') {
        logFile = fopen64(LogFileName, "w");
        if (logFile == NULL) {
            fprintf(stderr, "Cannot open log file!\n");
            exit(1);
        }
    } else {
        logFile = NULL;
    }

    packedKey = MMPoolDispatch(mmPool, WordPackedLengthFromText(MAX_SEARCH_PATTERN_LENGTH, BIT_PER_CHAR));
    convertedKey = MMPoolDispatch(mmPool, MAX_SEARCH_PATTERN_LENGTH);

    saIndexGroup = MMUnitAllocate(MaxNumberOfHitGroups * sizeof(SaIndexGroup));
    hitList = MMUnitAllocate(MaxNumberOfTextPosition * sizeof(HitList));
    tempSaIndex1 = MMUnitAllocate(MaxNumberOfTextPosition * sizeof(SaIndexList));
    tempSaIndex2 = MMUnitAllocate(MaxNumberOfTextPosition * sizeof(SaIndexList));

    elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
    printf("Elapsed time = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
    totalElapsedTime += elapsedTime;
    printf("\n");

    // Process search pattern files
    for (i=1; i<=NumberOfSearchPatternFile; i++) {

        argumentText[9] = '0' + (char)(i + 2);
        iniparser_copystring(programInput, argumentText, PatternFileName, PatternFileName, MAX_FILENAME_LEN);

        printf("Loading search pattern : %s\n", PatternFileName);
        numberOfPattern = ProcessSearchPattern();
        printf("Finished loading search pattern.\n");
        elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
        printf("Elapsed time = ");
        printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
        totalElapsedTime += elapsedTime;
        printf("\n");

        if (LogFileName[0] != '\0') {
            fprintf(logFile, "Searching for pattern : %s\n", PatternFileName);
        }

        if (SABinarySearch == TRUE) {

            printf("Start forward search with SA index.\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Forward search with SA index.\n");
            }

            numberOfPatternFound = 0;
            totalTextPositionMatched = 0;
            totalTextPositionRetrieved = 0;
            BWTInitializeSaRetrievalStatistics(&BWTSaRetrievalStatistics);

            for (j=0; j<numberOfPattern; j++) {
                //ConvertTextToWordPacked(searchPattern + searchPatternPosition[j], packedKey, charMap, ALPHABET_SIZE, 
                //                    searchPatternPosition[j+1] - searchPatternPosition[j]);
                //found = BWTForwardSearchSaIndex(packedKey, searchPatternPosition[j+1] - searchPatternPosition[j],
                //                         bwt, hsp->packedDNA, &saIndexLeft, &saIndexRight);

                ConvertTextToCode(searchPattern + searchPatternPosition[j], convertedKey, charMap, 
                                    searchPatternPosition[j+1] - searchPatternPosition[j]);
                found = BWTSaBinarySearch(convertedKey, searchPatternPosition[j+1] - searchPatternPosition[j],
                                         bwt, hsp->packedDNA, &saIndexLeft, &saIndexRight, packedKey);

                if (found) {
                    numberOfPatternFound++;
                    textPositionMatched = saIndexRight - saIndexLeft + 1;
                    totalTextPositionMatched += textPositionMatched;

                    if (FindTextPosition) {
                        saIndexGroup->startSaIndex = saIndexLeft;
                        if (textPositionMatched <= MaxNumberOfTextPosition) {
                            saIndexGroup->numOfMatch = textPositionMatched;
                        } else {
                            saIndexGroup->numOfMatch = MaxNumberOfTextPosition;
                            printf("Max no of text positions reached! Only %u out of %u text positions retrieved.\n", MaxNumberOfTextPosition, textPositionMatched);
                        }
                        saIndexGroup->posQuery = 0;
                        saIndexGroup->info = 0;
                        textPositionRetrieved = BWTTextPosition(bwt, saIndexGroup, 1, hitList, tempSaIndex1, tempSaIndex2, &BWTSaRetrievalStatistics, FALSE);
                        totalTextPositionRetrieved += textPositionRetrieved;
                    }
                }

                if (LogFileName[0] != '\0') {
                    if (found) {
                        fprintf(logFile, "Pattern number %u found. SA Index from %u to %u.  ", 
                                j + 1, saIndexLeft, saIndexRight);
                        if (FindTextPosition) {
                            fprintf(logFile, "%u matches of %u located. ", textPositionRetrieved, textPositionMatched);
                            for (k=0; k<textPositionRetrieved; k++) {
                                fprintf(logFile, "%u ", hitList[k].posText);
                            }
                        }
                    } else {
                        fprintf(logFile, "Pattern number %u not found.", j + 1); 
                    }
                    fprintf(logFile, "\n");
                }
            }

            printf("Finished forward search with SA index.\n");
            printf("%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
            if (FindTextPosition) {
                printf("Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                              BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                              BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
            }
            elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
            printf("Elapsed time = ");
            printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
            totalElapsedTime += elapsedTime;
            printf("\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Finished forward search with SA index.\n");
                fprintf(logFile, "%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
                if (FindTextPosition) {
                    fprintf(logFile, "Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                                  BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                                  BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
                }
                fprintf(logFile, "Elapsed time = ");
                printElapsedTime(logFile, FALSE, FALSE, TRUE, 4, elapsedTime);
                fprintf(logFile, "\n");
            }

        }

        if (BackwardSearch == TRUE) {

            printf("Start backward search.\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Backward search.\n");
            }

            numberOfPatternFound = 0;
            totalTextPositionMatched = 0;
            totalTextPositionRetrieved = 0;
            BWTInitializeSaRetrievalStatistics(&BWTSaRetrievalStatistics);

            for (j=0; j<numberOfPattern; j++) {
                ConvertTextToCode(searchPattern + searchPatternPosition[j], convertedKey, charMap, 
                                    searchPatternPosition[j+1] - searchPatternPosition[j]);
                found = BWTBackwardSearch(convertedKey, searchPatternPosition[j+1] - searchPatternPosition[j],
                                         bwt, &saIndexLeft, &saIndexRight);
                if (found) {
                    numberOfPatternFound++;
                    textPositionMatched = saIndexRight - saIndexLeft + 1;
                    totalTextPositionMatched += textPositionMatched;

                    if (FindTextPosition) {
                        saIndexGroup->startSaIndex = saIndexLeft;
                        if (textPositionMatched <= MaxNumberOfTextPosition) {
                            saIndexGroup->numOfMatch = textPositionMatched;
                        } else {
                            saIndexGroup->numOfMatch = MaxNumberOfTextPosition;
                            printf("Max no of text positions reached! Only %u out of %u text positions retrieved.\n", MaxNumberOfTextPosition, textPositionMatched);
                        }
                        saIndexGroup->posQuery = 0;
                        saIndexGroup->info = 0;
                        textPositionRetrieved = BWTTextPosition(bwt, saIndexGroup, 1, hitList, tempSaIndex1, tempSaIndex2, &BWTSaRetrievalStatistics, FALSE);
                        totalTextPositionRetrieved += textPositionRetrieved;
                    }
                }

                if (LogFileName[0] != '\0') {
                    if (found) {
                        fprintf(logFile, "Pattern number %u found. SA Index from %u to %u. Total %u occurrences ", j + 1, saIndexLeft, saIndexRight, saIndexRight - saIndexLeft + 1);
                        if (FindTextPosition) {
                            fprintf(logFile, "%u matches of %u located. ", textPositionRetrieved, textPositionMatched);
                            for (k=0; k<textPositionRetrieved; k++) {
                                fprintf(logFile, "%u ", hitList[k].posText);
                            }
                        }
                    } else {
                        fprintf(logFile, "Pattern number %u not found.", j + 1); 
                    }
                    fprintf(logFile, "\n");
                }
            }

            printf("Finished backward search.\n");
            printf("%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
            if (FindTextPosition) {
                printf("Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                              BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                              BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
            }
            elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
            printf("Elapsed time = ");
            printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
            totalElapsedTime += elapsedTime;
            printf("\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Finished backward search.\n");
                fprintf(logFile, "%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
                if (FindTextPosition) {
                    fprintf(logFile, "Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                                  BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                                  BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
                }
                fprintf(logFile, "Elapsed time = ");
                printElapsedTime(logFile, FALSE, FALSE, TRUE, 4, elapsedTime);
                fprintf(logFile, "\n");
            }

        }

        if (BackwardSearchCheck == TRUE) {

            printf("Start backward search with text.\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Backward search with text.\n");
            }

            numberOfPatternFound = 0;
            totalTextPositionMatched = 0;
            totalTextPositionRetrieved = 0;
            BWTInitializeSaRetrievalStatistics(&BWTSaRetrievalStatistics);

            for (j=0; j<numberOfPattern; j++) {
                ConvertTextToCode(searchPattern + searchPatternPosition[j], convertedKey, charMap, 
                                    searchPatternPosition[j+1] - searchPatternPosition[j]);
                ConvertTextToWordPacked(searchPattern + searchPatternPosition[j], packedKey, charMap, ALPHABET_SIZE, 
                                    searchPatternPosition[j+1] - searchPatternPosition[j]);
                found = BWTBackwardSearchCheckWithText(convertedKey, packedKey, searchPatternPosition[j+1] - searchPatternPosition[j],
                                                  bwt, hsp->packedDNA, TextCheckCostFactor, MaxNumberOfTextPosition, 
                                                  hitList, &saIndexLeft, &saIndexRight);
                if (found) {
                    numberOfPatternFound++;
                    textPositionMatched = saIndexRight - saIndexLeft + 1;
                    totalTextPositionMatched += textPositionMatched;

                    // Find text position if text check not used
                    if (FindTextPosition && saIndexLeft != 0) {
                        saIndexGroup->startSaIndex = saIndexLeft;
                        if (textPositionMatched <= MaxNumberOfTextPosition) {
                            saIndexGroup->numOfMatch = textPositionMatched;
                        } else {
                            saIndexGroup->numOfMatch = MaxNumberOfTextPosition;
                            printf("Max no of text positions reached! Only %u out of %u text positions retrieved.\n", MaxNumberOfTextPosition, textPositionMatched);
                        }
                        saIndexGroup->posQuery = 0;
                        saIndexGroup->info = 0;
                        textPositionRetrieved = BWTTextPosition(bwt, saIndexGroup, 1, hitList, tempSaIndex1, tempSaIndex2, &BWTSaRetrievalStatistics, FALSE);
                        totalTextPositionRetrieved += textPositionRetrieved;
                    } else {
                        textPositionRetrieved = textPositionMatched;
                        totalTextPositionRetrieved += textPositionRetrieved;
                    }
                }

                if (LogFileName[0] != '\0') {
                    if (found) {
                        fprintf(logFile, "Pattern number %u found. ", j + 1);
                        if (FindTextPosition) {
                            fprintf(logFile, "%u matches of %u located. ", textPositionRetrieved, textPositionMatched);
                            for (k=0; k<textPositionRetrieved; k++) {
                                fprintf(logFile, "%u ", hitList[k].posText);
                            }
                        }
                    } else {
                        fprintf(logFile, "Pattern number %u not found.", j + 1); 
                    }
                    fprintf(logFile, "\n");
                }
            }

            printf("Finished backward search with text.\n");
            printf("%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
            if (FindTextPosition) {
                printf("Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                              BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                              BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
            }
            elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
            printf("Elapsed time = ");
            printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
            totalElapsedTime += elapsedTime;
            printf("\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Finished backward search with text.\n");
                fprintf(logFile, "%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
                if (FindTextPosition) {
                    fprintf(logFile, "Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                                  BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                                  BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
                }
                fprintf(logFile, "Elapsed time = ");
                printElapsedTime(logFile, FALSE, FALSE, TRUE, 4, elapsedTime);
                fprintf(logFile, "\n");
            }

        }

        if (HammingDistSearch == TRUE) {

            printf("Start hamming distance %u approximate match.\n", MaxErrorAllowed);
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Hamming distance %u approximate match.\n", MaxErrorAllowed);
            }

            numberOfPatternFound = 0;
            totalTextPositionMatched = 0;
            totalTextPositionRetrieved = 0;
            BWTInitializeSaRetrievalStatistics(&BWTSaRetrievalStatistics);

            for (j=0; j<numberOfPattern; j++) {
                ConvertTextToCode(searchPattern + searchPatternPosition[j], convertedKey, charMap, 
                                    searchPatternPosition[j+1] - searchPatternPosition[j]);
                numberOfSaIndexGroup = BWTHammingDistMatchOld(convertedKey, searchPatternPosition[j+1] - searchPatternPosition[j],
                                                           bwt, MaxErrorAllowed, saIndexGroup, MaxNumberOfHitGroups, 0, 0);
                if (numberOfSaIndexGroup) {
                    numberOfPatternFound++;
                    textPositionMatched = 0;
                    for (k=0; k<numberOfSaIndexGroup; k++) {
                        textPositionMatched += saIndexGroup[k].numOfMatch;
                    }
                    if (textPositionMatched > MaxNumberOfTextPosition) {
                        textPositionMatched = 0;
                        for (k=0; k<numberOfSaIndexGroup && textPositionMatched < MaxNumberOfTextPosition; k++) {
                            textPositionMatched += saIndexGroup[k].numOfMatch;
                        }
                        numberOfSaIndexGroup = k - 1;
                        saIndexGroup[numberOfSaIndexGroup].numOfMatch -= textPositionMatched - MaxNumberOfTextPosition;
                        for (; k<numberOfSaIndexGroup; k++) {
                            textPositionMatched += saIndexGroup[k].numOfMatch;
                        }
                        printf("Max no of text positions reached! Only %u out of %u text positions retrieved.\n", MaxNumberOfTextPosition, textPositionMatched);
                    }
                    totalTextPositionMatched += textPositionMatched;

                    if (FindTextPosition) {
                        textPositionRetrieved = BWTTextPosition(bwt, saIndexGroup, numberOfSaIndexGroup, hitList, tempSaIndex1, tempSaIndex2, &BWTSaRetrievalStatistics, FALSE);
                        totalTextPositionRetrieved += textPositionRetrieved;
                    }
                }

                if (LogFileName[0] != '\0') {
                    if (numberOfSaIndexGroup) {
                        fprintf(logFile, "Pattern number %u found. %u matches", j + 1, textPositionMatched);
                        if (FindTextPosition) {
                            fprintf(logFile, "%u matches of %u located. ", textPositionRetrieved, textPositionMatched);
                            for (k=0; k<textPositionRetrieved; k++) {
                                fprintf(logFile, "%u ", hitList[k].posText);
                            }
                        }
                    } else {
                        fprintf(logFile, "Pattern number %u not found.", j + 1); 
                    }
                    fprintf(logFile, "\n");
                }
            }

            printf("Finished hamming distance search.\n");
            printf("%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
            if (FindTextPosition) {
                printf("Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                              BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                              BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
            }
            elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
            printf("Elapsed time = ");
            printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
            totalElapsedTime += elapsedTime;
            printf("\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Finished hamming distance search.\n");
                fprintf(logFile, "%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
                if (FindTextPosition) {
                    fprintf(logFile, "Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                                  BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                                  BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
                }
                fprintf(logFile, "Elapsed time = ");
                printElapsedTime(logFile, FALSE, FALSE, TRUE, 4, elapsedTime);
                fprintf(logFile, "\n");
            }

        }

        if (EditDistSearch == TRUE) {

            printf("Start edit distance %u approximate match.\n", MaxErrorAllowed);
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Edit distance %u approximate match.\n", MaxErrorAllowed);
            }

            numberOfPatternFound = 0;
            totalTextPositionMatched = 0;
            totalTextPositionRetrieved = 0;
            BWTInitializeSaRetrievalStatistics(&BWTSaRetrievalStatistics);

            for (j=0; j<numberOfPattern; j++) {
                ConvertTextToCode(searchPattern + searchPatternPosition[j], convertedKey, charMap, 
                                    searchPatternPosition[j+1] - searchPatternPosition[j]);
                numberOfSaIndexGroup = BWTEditDistMatchOld(convertedKey, searchPatternPosition[j+1] - searchPatternPosition[j],
                                                        bwt, MaxErrorAllowed, (SaIndexGroupWithLengthError*)saIndexGroup, MaxNumberOfHitGroups);
                if (numberOfSaIndexGroup > MaxNumberOfHitGroups) {
                    fprintf(stderr, "numberOfSaIndexGroup > MaxNumberOfHitGroups!\n");
                }
                if (numberOfSaIndexGroup) {
                    numberOfPatternFound++;
                    textPositionMatched = 0;
                    for (k=0; k<numberOfSaIndexGroup; k++) {
                        textPositionMatched += saIndexGroup[k].numOfMatch;
                    }
                    if (textPositionMatched > MaxNumberOfTextPosition) {
                        textPositionMatched = 0;
                        for (k=0; k<numberOfSaIndexGroup && textPositionMatched < MaxNumberOfTextPosition; k++) {
                            textPositionMatched += saIndexGroup[k].numOfMatch;
                        }
                        numberOfSaIndexGroup = k - 1;
                        saIndexGroup[numberOfSaIndexGroup].numOfMatch -= textPositionMatched - MaxNumberOfTextPosition;
                        for (; k<numberOfSaIndexGroup; k++) {
                            textPositionMatched += saIndexGroup[k].numOfMatch;
                        }
                        printf("Max no of text positions reached! Only %u out of %u text positions retrieved.\n", MaxNumberOfTextPosition, textPositionMatched);
                    }
                    totalTextPositionMatched += textPositionMatched;

                    if (FindTextPosition) {
                        textPositionRetrieved = BWTTextPosition(bwt, saIndexGroup, numberOfSaIndexGroup, hitList, tempSaIndex1, tempSaIndex2, &BWTSaRetrievalStatistics, FALSE);
                        totalTextPositionRetrieved += textPositionRetrieved;
                    }
                }

                if (LogFileName[0] != '\0') {
                    if (numberOfSaIndexGroup) {
                        fprintf(logFile, "Pattern number %u found. %u matches. ", j + 1, textPositionMatched);
                        if (FindTextPosition) {
                            fprintf(logFile, "%u matches of %u located. ", textPositionRetrieved, textPositionMatched);
                            for (k=0; k<textPositionRetrieved; k++) {
                                fprintf(logFile, "%u ", hitList[k].posText);
                            }
                        }
                    } else {
                        fprintf(logFile, "Pattern number %u not found.", j + 1); 
                    }
                    fprintf(logFile, "\n");
                }
            }

            printf("Finished edit distance search.\n");
            printf("%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
            if (FindTextPosition) {
                printf("Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                              BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                              BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
            }
            elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
            printf("Elapsed time = ");
            printElapsedTime(stdout, FALSE, FALSE, TRUE, 4, elapsedTime);
            totalElapsedTime += elapsedTime;
            printf("\n");
            if (LogFileName[0] != '\0') {
                fprintf(logFile, "Finished edit distance search.\n");
                fprintf(logFile, "%u patterns out of %u found. %llu of %llu matches located.\n", numberOfPatternFound, numberOfPattern, totalTextPositionRetrieved, totalTextPositionMatched);
                if (FindTextPosition) {
                    fprintf(logFile, "Hits: BWT/Cached/Diag/Dup      : %u/%u/%u/%u\n", 
                                  BWTSaRetrievalStatistics.bwtSaRetrieved, BWTSaRetrievalStatistics.cachedSaRetrieved, 
                                  BWTSaRetrievalStatistics.saDiagonalLinked, BWTSaRetrievalStatistics.saDuplicated);
                }
                fprintf(logFile, "Elapsed time = ");
                printElapsedTime(logFile, FALSE, FALSE, TRUE, 4, elapsedTime);
                fprintf(logFile, "\n");
            }

        }

    }

    MMPoolReturn(mmPool, packedKey, WordPackedLengthFromText(MAX_SEARCH_PATTERN_LENGTH, BIT_PER_CHAR));
    MMPoolReturn(mmPool, convertedKey, MAX_SEARCH_PATTERN_LENGTH);

    HSPFree(mmPool, hsp, MAX_SEARCH_PATTERN_LENGTH / CHAR_PER_WORD);
    BWTFree(mmPool, bwt);

    if (searchPattern != NULL) {
        MMUnitFree(searchPattern, searchPatternAllocated);
    }
    if (searchPatternPosition != NULL) {
        MMUnitFree(searchPatternPosition, searchPatternPositionAllocated);
    }
    MMUnitFree(saIndexGroup, MaxNumberOfHitGroups * sizeof(unsigned int));
    MMUnitFree(hitList, MaxNumberOfTextPosition * sizeof(unsigned int));

    MMPoolFree(mmPool);

    iniparser_freedict(programInput);
    iniparser_freedict(ini);

    return 0;

}

dictionary *ParseInput(int argc, char** argv) {

    dictionary *programInput;
    char t1[3] = "-c";    // specify that this is a boolean type parameter
    char *d[1];

    d[0] = t1;

    programInput = paraparser_load(argc, argv, 1, d);

    // Check if ini file is specified
    iniparser_copystring(programInput, "argument:1", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
    if (DatabaseName[0] == '\0') {
        fprintf(stderr, "Syntax: %s <database> <ini> <pattern files> [-c Confirm] -l <log file>\n", argv[0]);
        exit(1);
    }
    iniparser_copystring(programInput, "argument:2", IniFileName, IniFileName, MAX_FILENAME_LEN);
    if (IniFileName[0] == '\0') {
        fprintf(stderr, "Syntax: %s <database> <ini> <pattern files> [-c Confirm] -l <log file>\n", argv[0]);
        exit(1);
    }
    iniparser_copystring(programInput, "argument:3", PatternFileName, PatternFileName, MAX_FILENAME_LEN);
    if (PatternFileName[0] == '\0') {
        fprintf(stderr, "Syntax: %s <database> <ini> <pattern files> [-c Confirm] -l <log file>\n", argv[0]);
        exit(1);
    }

    NumberOfSearchPatternFile = argc - 3;

    // Whether confirmation is needed
    Confirmation = iniparser_find_entry(programInput, "parameter:-c");

    // Log file
    iniparser_copystring(programInput, "-l:3", LogFileName, LogFileName, MAX_FILENAME_LEN);

    return programInput;

}

dictionary *ParseIniFile() {

    dictionary *ini;

    ini = iniparser_load(IniFileName, FALSE);
    if (ini == NULL) {
        fprintf(stderr, "Cannot open ini file\n");
        exit(1);
    }

    // Action parameters
    SABinarySearch = iniparser_getboolean(ini, "Action:SABinarySearch", FALSE);
    BackwardSearch = iniparser_getboolean(ini, "Action:BackwardSearch", FALSE);
    BackwardSearchCheck = iniparser_getboolean(ini, "Action:BackwardSearchCheck", FALSE);
    HammingDistSearch = iniparser_getboolean(ini, "Action:HammingDistSearch", FALSE);
    EditDistSearch = iniparser_getboolean(ini, "Action:EditDistSearch", FALSE);
    FindTextPosition = iniparser_getboolean(ini, "Action:FindTextPosition", FALSE);
    TextCheckCostFactor = iniparser_getint(ini, "Action:TextCheckCostFactor", 20);
    MaxNumberOfTextPosition = iniparser_getint(ini, "Action:MaxNumberOfTextPosition", 100000);
    MaxNumberOfHitGroups = iniparser_getint(ini, "Action:MaxNumberOfHitGroups", 10000);

    // Memory parameters
    PoolSize = iniparser_getint(ini, "Memory:PoolSize", 4096);

    // Search pattern parameters
    MaxErrorAllowed = iniparser_getint(ini, "SearchPattern:MaxErrorAllowed", 0);
    
    // Database parameters
    iniparser_copystring(ini, "Database:Location", DatabaseLocation, DatabaseLocation, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:AnnotationFileName", AnnotationFileName, AnnotationFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:AmbiguityFileName", AmbiguityFileName, AmbiguityFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:PackedDNAFileName", PackedDNAFileName, PackedDNAFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTCodeFileName", BWTCodeFileName, BWTCodeFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:BWTOccValueFileName", BWTOccValueFileName, BWTOccValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:SaValueFileName", SaValueFileName, SaValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:CachedSaIndexFileName", SaIndexFileName, SaIndexFileName, MAX_FILENAME_LEN);


    return ini;

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

    if (SABinarySearch == FALSE && 
        BackwardSearch == FALSE && BackwardSearchCheck == FALSE &&
        HammingDistSearch == FALSE && EditDistSearch == FALSE) {
        fprintf(stderr, "No action is specified!\n");
        exit(1);
    }

    if (BackwardSearchCheck == TRUE && TextCheckCostFactor == 0) {
        fprintf(stderr, "Text checking cost factor must be set!\n");
        exit(1);
    }

    if (FindTextPosition == TRUE && MaxNumberOfTextPosition == 0) {
        fprintf(stderr, "Must set a maximum number of text position to report!\n");
        exit(1);
    }

    if ((HammingDistSearch == TRUE || EditDistSearch == TRUE) &&
        MaxErrorAllowed == 0) {
        fprintf(stderr, "Max error allowed not specified!\n");
        exit(1);
    }

    if (NumberOfSearchPatternFile == 0) {
        fprintf(stderr, "No search pattern file is specified!\n");
        exit(1);
    }
    
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

    printf("\n");
    printf("SA Binary Search              : %c\n", 
        boolean[SABinarySearch]);
    printf("Backward Search               : %c\n", 
        boolean[BackwardSearch]);
    printf("Backward Search Text Check    : %c   Text Check Cost Factor               : %u\n", 
        boolean[BackwardSearchCheck], TextCheckCostFactor);

    printf("Hamming Distance Match        : %c   Edit Distance Match                  : %c\n", 
        boolean[HammingDistSearch], boolean[EditDistSearch]);

    printf("Find Text Position            : %c   Max Number Of Position to Report     : %u\n", 
        boolean[FindTextPosition], MaxNumberOfTextPosition);
    printf("\n");

    if (HammingDistSearch || EditDistSearch) {
        printf("Max Error Allowed        : %u\n", MaxErrorAllowed);
    }

    printf("Pattern File (x%u)        : %s\n", NumberOfSearchPatternFile, PatternFileName);
    printf("\n");

    if (LogFileName[0] != '\0') {
        printf("Log File                 : %s\n", LogFileName);
        printf("\n");
    }

    printf("Annotation file          : %s\n", AnnotationFileName);
    printf("Ambiguity file           : %s\n", AmbiguityFileName);
    printf("Packed DNA file          : %s\n", PackedDNAFileName);
    printf("BWT Code file            : %s\n", BWTCodeFileName);
    printf("BWT Occ value file       : %s\n", BWTOccValueFileName);
    printf("SA value file            : %s\n", SaValueFileName);
    printf("SA index file            : %s\n", SaIndexFileName);
    printf("\n");

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

unsigned int ProcessSearchPattern() {

    FILE *searchPatternFile;
    char c;
    unsigned int i, j, n;

    searchPatternFile = fopen64(PatternFileName, "r");
    if (searchPatternFile == NULL) {
        fprintf(stderr, "Cannot open search pattern file : %s\n", PatternFileName);
        exit(1);
    }

    j = 0;
    n = 0;
    // Find out the size of search pattern
    while (!feof(searchPatternFile)) {
        c = (char)fgetc(searchPatternFile);
        // Find the first comment line
        while (!feof(searchPatternFile) && c != '>') {
            c = (char)fgetc(searchPatternFile);
        }
        while (!feof(searchPatternFile)) {
            // Skip comment line
            while (!feof(searchPatternFile) && c != '\n') {
                c = (char)fgetc(searchPatternFile);
            }
            // Count no. of characters
            while (!feof(searchPatternFile) && c != '>') {
                if (c != '\n') {
                    j++;
                }
                c = (char)fgetc(searchPatternFile);
            }
            // Count no. of patterns
            n++;
        }
    }

    if (j > searchPatternAllocated) {
        if (searchPattern != NULL) {
            MMUnitFree(searchPattern, searchPatternAllocated);
        }
        searchPattern = MMUnitAllocate(j);
        searchPatternAllocated = j;
    }

    if (n + 1 > searchPatternPositionAllocated) {
        if (searchPatternPosition != NULL) {
            MMUnitFree(searchPatternPosition, searchPatternPositionAllocated);
        }
        searchPatternPosition = MMUnitAllocate((n + 1) * sizeof(unsigned int));
        searchPatternPositionAllocated = n + 1;
    }

    fseek(searchPatternFile, 0, SEEK_SET);

    j = 0;
    n = 0;
    i = 0;

    searchPatternPosition[0] = 0;

    // Read pattern into memory
    while (!feof(searchPatternFile)) {
        c = (char)fgetc(searchPatternFile);
        // Find the first comment line
        while (!feof(searchPatternFile) && c != '>') {
            c = (char)fgetc(searchPatternFile);
        }
        while (!feof(searchPatternFile)) {
            // Skip comment line
            while (!feof(searchPatternFile) && c != '\n') {
                c = (char)fgetc(searchPatternFile);
            }
            // Store pattern
            while (!feof(searchPatternFile) && c != '>') {
                if (c != '\n') {
                    searchPattern[j] = c;
                    j++;
                }
                c = (char)fgetc(searchPatternFile);
            }
            // Store pointers to pattern
            n++;
            searchPatternPosition[n] = j;
        }
    }


    fclose(searchPatternFile);

    return n;

}

