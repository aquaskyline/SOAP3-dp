/*

   BlastCompare.c        Compare blast results in table (-m 8) format

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


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MiscUtilities.h"
#include "iniparser.h"

#define DB_SEQ_NAME_ALLOCATION_UNIT        256
#define BLAST_RESULT_ALLOCATION_UNIT    65536
#define MAX_EVALUE    1e+1

double minEvalue = 1e-100;
int maxRank = 100;
int compareByEvalue = TRUE;
int summary = FALSE;
int averageLogEvalue = FALSE;

typedef struct BlastResult {
    int dbSeq;
    unsigned int queryStart;
    unsigned int queryEnd;
    unsigned int textStart;
    unsigned int textEnd;
    double evalue;
    int minRank;
    int maxRank;
    char resultLine[MAX_FILENAME_LEN+1];
} BlastResult;


int BlastResultDbSeq(const void *blastResult, const int index1, const int index2) {

    return ((BlastResult*)blastResult + index1)->dbSeq - ((BlastResult*)blastResult + index2)->dbSeq;

}

int BlastResultEvalueSeq(const void *blastResult, const int index1, const int index2) {

    if (((BlastResult*)blastResult + index1)->evalue > ((BlastResult*)blastResult + index2)->evalue) {
        return 1;
    } else if (((BlastResult*)blastResult + index1)->evalue < ((BlastResult*)blastResult + index2)->evalue) {
        return -1;
    } else {
        return 0;
    }

}


dictionary *ParseInput(int argc, char** argv);
void ValidateParameters();

    char ResultFileName1[MAX_FILENAME_LEN+1] = "";
    char ResultFileName2[MAX_FILENAME_LEN+1] = "";
    char UnmatchedFileName1[MAX_FILENAME_LEN+1] = "";
    char UnmatchedFileName2[MAX_FILENAME_LEN+1] = "";

    char QueryName[MAX_FILENAME_LEN+1] = "";


int main(int argc, char** argv) {

    char **dbSeqName;
    BlastResult *blastResult1, *blastResult2;
    int blastResult1Allocated, blastResult2Allocated;
    int numOfBlastResult1 = 0, numOfBlastResult2 = 0;
    int numOfDbSeqAllocated;
    int numOfDbSeq = 0;

    char temp[MAX_FILENAME_LEN+1] = "";
    char tempResultLine[MAX_FILENAME_LEN+1] = "";
    char tempDbSeqName[MAX_FILENAME_LEN+1] = "";
    unsigned int tempQueryStart, tempQueryEnd;
    unsigned int tempTextStart, tempTextEnd;
    double tempEvalue;

    double *matched1, *matched2;
    double *unmatched1, *unmatched2;
    double *totalLogEvalue1, *totalLogEvalue2;
    int numOfSlot;
    int minSlot, maxSlot;
    double weight, residualWeight;

    FILE *resultFile1, *resultFile2;
    FILE *unmatchedFile1 = NULL, *unmatchedFile2 = NULL;

    dictionary *programInput;

    int i, j, k;
    int startIndex;

    double cumulativeMatched1, cumulativeMatched2;
    double cumulativeUnmatched1, cumulativeUnmatched2;
    double cumulativetotalLogEvalue1, cumulativetotalLogEvalue2;

    if (argc < 4 && (argc < 2 || (argv[1][0] != '-' && argv[1][0] != 's'))) {
        printf("BlastCompare, Copyright (C) 2006, Wong Chi Kwong.\n");
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

        printf("Syntax: blastcompare <result file 1> <result file 2> <query name>\n");
        printf("                     [unmatched file 1] [unmatched file 2]\n");
        printf("                     -e <min E-Value>\n");
        printf("                     -r <max rank>\n");
        printf("                     -i <rank interval>\n");
        exit(1);

    }

    programInput = ParseInput(argc, argv);
    ValidateParameters();

    if (compareByEvalue) {
        numOfSlot = (int)(log10(MAX_EVALUE) - log10(minEvalue)) + 1;
    } else {
        numOfSlot = maxRank + 1;
    }

    matched1 = calloc(numOfSlot, sizeof(double));
    matched2 = calloc(numOfSlot, sizeof(double));
    unmatched1 = calloc(numOfSlot, sizeof(double));
    unmatched2 = calloc(numOfSlot, sizeof(double));
    totalLogEvalue1 = calloc(numOfSlot, sizeof(double));
    totalLogEvalue2 = calloc(numOfSlot, sizeof(double));

    if (summary) {
        fgets(tempResultLine, MAX_FILENAME_LEN, stdin);
        while (!feof(stdin)) {
            for (i=0; i<MAX_FILENAME_LEN; i++) {
                if (tempResultLine[i] == ' ') {
                    tempResultLine[i] = '_';
                }
                if (tempResultLine[i] == '\t') {
                    break;
                }
            }
            if (averageLogEvalue) {
                sscanf(tempResultLine, "%s %lf %lf %lf %lf %lf %lf %lf\n", 
                    temp, &tempEvalue, &cumulativeMatched1, &cumulativeUnmatched1, &cumulativetotalLogEvalue1, &cumulativeMatched2, &cumulativeUnmatched2, &cumulativetotalLogEvalue2);
            } else {
                sscanf(tempResultLine, "%s %lf %lf %lf %lf %lf\n", 
                    temp, &tempEvalue, &cumulativeMatched1, &cumulativeUnmatched1, &cumulativeMatched2, &cumulativeUnmatched2);
            }
            if (compareByEvalue) {
                if (tempEvalue == 0.0) {
                    tempEvalue = minEvalue;
                }
                minSlot = (int)ceil(log10(tempEvalue) - log10(minEvalue));
                if (minSlot < 0) {
                    minSlot = 0;
                }
            } else {
                minSlot = (int)(tempEvalue) - 1;
            }
            if (minSlot < numOfSlot) {
                matched1[minSlot] += cumulativeMatched1;
                unmatched1[minSlot] += cumulativeUnmatched1;
                matched2[minSlot] += cumulativeMatched2;
                unmatched2[minSlot] += cumulativeUnmatched2;
                if (cumulativetotalLogEvalue1 > 0) {
                    totalLogEvalue1[minSlot] += cumulativeMatched1 * log10(cumulativetotalLogEvalue1);
                }
                if (cumulativetotalLogEvalue2 > 0) {
                    totalLogEvalue2[minSlot] += cumulativeMatched2 * log10(cumulativetotalLogEvalue2);
                }
            }
            fgets(tempResultLine, MAX_FILENAME_LEN, stdin);
        }

        for (i=0; i<numOfSlot; i++) {
            if (compareByEvalue) {
                tempEvalue = pow(10, i + (int)log10(minEvalue));
                printf("%.0le", tempEvalue);
            } else {
                printf("%d", i + 1);
            }
            printf("\t%.0lf", matched1[i]);
            if (compareByEvalue) {
                printf("\t%.0lf", unmatched1[i]);
            } else {
                printf("\t%.1lf", unmatched1[i]);
            }
            if (averageLogEvalue) {
                if (matched1[i] > 0) {
                    printf("\t%.0le", pow(10, totalLogEvalue1[i] / matched1[i]));
                } else {
                    printf("\t0");
                }
            }
            printf("\t%.0lf", matched2[i]);
            if (compareByEvalue) {
                printf("\t%.0lf", unmatched2[i]);
            } else {
                printf("\t%.1lf", unmatched2[i]);
            }
            if (averageLogEvalue) {
                if (matched2[i] > 0) {
                    printf("\t%.0le", pow(10, totalLogEvalue2[i] / matched2[i]));
                } else {
                    printf("\t0");
                }
            }
            printf("\n");
        }

        exit(0);
    }

    resultFile1 = fopen64(ResultFileName1, "r");
    if (resultFile1 == NULL) {
        fprintf(stderr, "Cannot open ResultFile1!\n");
        exit(1);
    }
    resultFile2 = fopen64(ResultFileName2, "r");
    if (resultFile2 == NULL) {
        fprintf(stderr, "Cannot open resultFile2!\n");
        exit(1);
    }
    if (UnmatchedFileName1[0] != '\0' && UnmatchedFileName1[0] != ' ' && UnmatchedFileName1[0] != '-') {
        unmatchedFile1 = fopen64(UnmatchedFileName1, "w");
        if (unmatchedFile1 == NULL) {
            fprintf(stderr, "Cannot open UnmatchedFile1!\n");
            exit(1);
        }
    }
    if (UnmatchedFileName2[0] != '\0' && UnmatchedFileName2[0] != ' ' && UnmatchedFileName2[0] != '-') {
        unmatchedFile2 = fopen64(UnmatchedFileName2, "w");
        if (unmatchedFile2 == NULL) {
            fprintf(stderr, "Cannot open UnmatchedFile2!\n");
            exit(1);
        }
    }

    dbSeqName = malloc(DB_SEQ_NAME_ALLOCATION_UNIT * sizeof(char*));
    numOfDbSeqAllocated = DB_SEQ_NAME_ALLOCATION_UNIT;

    blastResult1 = malloc(BLAST_RESULT_ALLOCATION_UNIT * sizeof(BlastResult));
    blastResult1Allocated = BLAST_RESULT_ALLOCATION_UNIT;
    blastResult2 = malloc(BLAST_RESULT_ALLOCATION_UNIT * sizeof(BlastResult));
    blastResult2Allocated = BLAST_RESULT_ALLOCATION_UNIT;

    fgets(tempResultLine, MAX_FILENAME_LEN, resultFile1);

    while(!feof(resultFile1)) {

        if (tempResultLine[0] == '#') {
            fgets(tempResultLine, MAX_FILENAME_LEN, resultFile1);
            continue;
        }

        sscanf(tempResultLine, "%s %s %s %s %s %s %u %u %u %u %lf %s\n", 
            temp, tempDbSeqName + 1, temp, temp, temp, temp,
            &tempQueryStart, &tempQueryEnd, &tempTextStart, &tempTextEnd, &tempEvalue, temp);

        if (tempDbSeqName[1] == '\0' || tempDbSeqName[1] == '\n') {
            break;
        }

        if (tempTextStart < tempTextEnd) {
            tempDbSeqName[0] = '+';
        } else {
            tempDbSeqName[0] = '-';
        }

        for (i=0; i<numOfDbSeq; i++) {
            if (strncmp(tempDbSeqName, dbSeqName[i], MAX_FILENAME_LEN) == 0) {
                break;
            }
        }
        if (i >= numOfDbSeq) {
            if (numOfDbSeq >= numOfDbSeqAllocated) {
                dbSeqName = realloc(dbSeqName, (numOfDbSeqAllocated + DB_SEQ_NAME_ALLOCATION_UNIT) * sizeof(char*));
                numOfDbSeqAllocated += DB_SEQ_NAME_ALLOCATION_UNIT;
            }
            dbSeqName[i] = malloc(MAX_FILENAME_LEN);
            memcpy(dbSeqName[i], tempDbSeqName, MAX_FILENAME_LEN);
            numOfDbSeq++;
        }

        if (numOfBlastResult1 >= blastResult1Allocated) {
            blastResult1 = realloc(blastResult1, (blastResult1Allocated + BLAST_RESULT_ALLOCATION_UNIT) * sizeof(BlastResult));
            blastResult1Allocated += BLAST_RESULT_ALLOCATION_UNIT;
        }
        blastResult1[numOfBlastResult1].dbSeq = i;
        blastResult1[numOfBlastResult1].queryStart = tempQueryStart;
        blastResult1[numOfBlastResult1].queryEnd = tempQueryEnd;
        if (tempTextStart < tempTextEnd) {
            blastResult1[numOfBlastResult1].textStart = tempTextStart;
            blastResult1[numOfBlastResult1].textEnd = tempTextEnd;
        } else {
            blastResult1[numOfBlastResult1].textStart = tempTextEnd;
            blastResult1[numOfBlastResult1].textEnd = tempTextStart;
        }
        blastResult1[numOfBlastResult1].evalue = tempEvalue;
        strncpy(blastResult1[numOfBlastResult1].resultLine, tempResultLine, MAX_FILENAME_LEN);
        numOfBlastResult1++;

        fgets(tempResultLine, MAX_FILENAME_LEN, resultFile1);

    }

    fclose(resultFile1);

    fgets(tempResultLine, MAX_FILENAME_LEN, resultFile2);

    while(!feof(resultFile2)) {

        if (tempResultLine[0] == '#') {
            fgets(tempResultLine, MAX_FILENAME_LEN, resultFile2);
            continue;
        }

        sscanf(tempResultLine, "%s %s %s %s %s %s %u %u %u %u %lf %s\n", 
            temp, tempDbSeqName + 1, temp, temp, temp, temp,
            &tempQueryStart, &tempQueryEnd, &tempTextStart, &tempTextEnd, &tempEvalue, temp);

        if (tempDbSeqName[1] == '\0' || tempDbSeqName[1] == '\n') {
            break;
        }

        if (tempTextStart < tempTextEnd) {
            tempDbSeqName[0] = '+';
        } else {
            tempDbSeqName[0] = '-';
        }

        for (i=0; i<numOfDbSeq; i++) {
            if (strncmp(tempDbSeqName, dbSeqName[i], MAX_FILENAME_LEN) == 0) {
                break;
            }
        }
        if (i >= numOfDbSeq) {
            if (numOfDbSeq >= numOfDbSeqAllocated) {
                dbSeqName = realloc(dbSeqName, (numOfDbSeqAllocated + DB_SEQ_NAME_ALLOCATION_UNIT) * sizeof(char*));
                numOfDbSeqAllocated += DB_SEQ_NAME_ALLOCATION_UNIT;
            }
            dbSeqName[i] = malloc(MAX_FILENAME_LEN);
            memcpy(dbSeqName[i], tempDbSeqName, MAX_FILENAME_LEN);
            numOfDbSeq++;
        }

        if (numOfBlastResult2 >= blastResult2Allocated) {
            blastResult2 = realloc(blastResult2, (blastResult2Allocated + BLAST_RESULT_ALLOCATION_UNIT) * sizeof(BlastResult));
            blastResult2Allocated += BLAST_RESULT_ALLOCATION_UNIT;
        }
        blastResult2[numOfBlastResult2].dbSeq = i;
        blastResult2[numOfBlastResult2].queryStart = tempQueryStart;
        blastResult2[numOfBlastResult2].queryEnd = tempQueryEnd;
        if (tempTextStart < tempTextEnd) {
            blastResult2[numOfBlastResult2].textStart = tempTextStart;
            blastResult2[numOfBlastResult2].textEnd = tempTextEnd;
        } else {
            blastResult2[numOfBlastResult2].textStart = tempTextEnd;
            blastResult2[numOfBlastResult2].textEnd = tempTextStart;
        }
        blastResult2[numOfBlastResult2].evalue = tempEvalue;
        strncpy(blastResult2[numOfBlastResult2].resultLine, tempResultLine, MAX_FILENAME_LEN);
        numOfBlastResult2++;

        fgets(tempResultLine, MAX_FILENAME_LEN, resultFile2);

    }

    fclose(resultFile2);

    // Determine rank
    QSort(blastResult1, numOfBlastResult1, sizeof(BlastResult), BlastResultEvalueSeq);
    QSort(blastResult2, numOfBlastResult2, sizeof(BlastResult), BlastResultEvalueSeq);

    j = 0;
    for (i=0; i<numOfBlastResult1; i++) {
        if (blastResult1[i].evalue == blastResult1[j].evalue) {
            blastResult1[i].minRank = j;
        } else {
            while (j < i) {
                blastResult1[j].maxRank = i - 1;
                j++;
            }
            blastResult1[i].minRank = i;
        }
    }
    while (j < i) {
        blastResult1[j].maxRank = i - 1;
        j++;
    }

    j = 0;
    for (i=0; i<numOfBlastResult2; i++) {
        if (blastResult2[i].evalue == blastResult2[j].evalue) {
            blastResult2[i].minRank = j;
        } else {
            while (j < i) {
                blastResult2[j].maxRank = i - 1;
                j++;
            }
            blastResult2[i].minRank = i;
        }
    }
    while (j < i) {
        blastResult2[j].maxRank = i - 1;
        j++;
    }


    QSort(blastResult1, numOfBlastResult1, sizeof(BlastResult), BlastResultDbSeq);
    QSort(blastResult2, numOfBlastResult2, sizeof(BlastResult), BlastResultDbSeq);

    startIndex = 0;
    for (i=0; i<numOfBlastResult1; i++) {
        while (startIndex < numOfBlastResult2 && blastResult2[startIndex].dbSeq < blastResult1[i].dbSeq) {
            startIndex++;
        }
        for (j=startIndex; j < numOfBlastResult2 && blastResult2[j].dbSeq == blastResult1[i].dbSeq; j++) {
            if (blastResult1[i].queryStart <= blastResult2[j].queryEnd && blastResult2[j].queryStart <= blastResult1[i].queryEnd &&
                blastResult1[i].textStart <= blastResult2[j].textEnd && blastResult2[j].textStart <= blastResult1[i].textEnd) {
                break;
            }
        }
        residualWeight = 0.0;
        if (compareByEvalue) {
            tempEvalue = blastResult1[i].evalue;
            if (tempEvalue == 0.0) {
                tempEvalue = minEvalue;
            }
            minSlot = (int)ceil(log10(tempEvalue) - log10(minEvalue));
            if (minSlot < 0) {
                minSlot = 0;
            }
            if (minSlot >= numOfSlot) {
                minSlot = numOfSlot - 1;
            }
            weight = 1.0;
            maxSlot = minSlot;
        } else {
            weight = (double)1 / (double)(blastResult1[i].maxRank - blastResult1[i].minRank + 1);
            minSlot = blastResult1[i].minRank;
            maxSlot = blastResult1[i].maxRank;
            if (maxSlot > maxRank) {
                if (minSlot > maxRank) {
                    residualWeight = 1;
                    weight = 0;
                    minSlot = maxRank;
                } else {
                    residualWeight = (double)(maxSlot - maxRank) * weight;
                }
                maxSlot = maxRank;
            }
        }
        if (j < numOfBlastResult2 && blastResult2[j].dbSeq == blastResult1[i].dbSeq) {
            for (k=minSlot; k<=maxSlot; k++) {
                matched1[k] += weight;
            }
            matched1[maxRank] += residualWeight;
        } else {
            for (k=minSlot; k<=maxSlot; k++) {
                unmatched1[k] += weight;
            }
            unmatched1[maxRank] += residualWeight;
            if (unmatchedFile1 != NULL) {
                fprintf(unmatchedFile1, "%s", blastResult1[i].resultLine);
            }
        }
        for (k=minSlot; k<=maxSlot; k++) {
            if (blastResult1[i].evalue < 1.0e-180) {
                totalLogEvalue1[k] += log10(1.0e-180) * weight;
            } else {
                totalLogEvalue1[k] += log10(blastResult1[i].evalue) * weight;
            }
        }
        if (blastResult1[i].evalue < 1.0e-180) {
            totalLogEvalue1[maxRank] += log10(1.0e-180) * residualWeight;
        } else {
            totalLogEvalue1[maxRank] += log10(blastResult1[i].evalue) * residualWeight;
        }
    }

    startIndex = 0;
    for (i=0; i<numOfBlastResult2; i++) {
        while (startIndex < numOfBlastResult1 && blastResult1[startIndex].dbSeq < blastResult2[i].dbSeq) {
            startIndex++;
        }
        for (j=startIndex; j < numOfBlastResult1 && blastResult1[j].dbSeq == blastResult2[i].dbSeq; j++) {
            if (blastResult2[i].queryStart <= blastResult1[j].queryEnd && blastResult1[j].queryStart <= blastResult2[i].queryEnd &&
                blastResult2[i].textStart <= blastResult1[j].textEnd && blastResult1[j].textStart <= blastResult2[i].textEnd) {
                break;
            }
        }
        residualWeight = 0.0;
        if (compareByEvalue) {
            tempEvalue = blastResult2[i].evalue;
            if (tempEvalue == 0.0) {
                tempEvalue = minEvalue;
            }
            minSlot = (int)ceil(log10(tempEvalue) - log10(minEvalue));
            if (minSlot < 0) {
                minSlot = 0;
            }
            if (minSlot >= numOfSlot) {
                minSlot = numOfSlot - 1;
            }
            weight = 1.0;
            maxSlot = minSlot;
        } else {
            weight = (double)1 / (double)(blastResult2[i].maxRank - blastResult2[i].minRank + 1);
            minSlot = blastResult2[i].minRank;
            maxSlot = blastResult2[i].maxRank;
            if (maxSlot > maxRank) {
                if (minSlot > maxRank) {
                    residualWeight = 1;
                    weight = 0;
                    minSlot = maxRank;
                } else {
                    residualWeight = (double)(maxSlot - maxRank) * weight;
                }
                maxSlot = maxRank;
            }
        }
        if (j < numOfBlastResult1 && blastResult1[j].dbSeq == blastResult2[i].dbSeq) {
            for (k=minSlot; k<=maxSlot; k++) {
                matched2[k] += weight;
            }
            matched2[maxRank] += residualWeight;
        } else {
            for (k=minSlot; k<=maxSlot; k++) {
                unmatched2[k] += weight;
            }
            unmatched2[maxRank] += residualWeight;
            if (unmatchedFile2 != NULL) {
                fprintf(unmatchedFile2, "%s", blastResult2[i].resultLine);
            }
        }
        for (k=minSlot; k<=maxSlot; k++) {
            if (blastResult2[i].evalue < 1.0e-180) {
                totalLogEvalue2[k] += log10(1.0e-180) * weight;
            } else {
                totalLogEvalue2[k] += log10(blastResult2[i].evalue) * weight;
            }
        }
        if (blastResult2[i].evalue < 1.0e-180) {
            totalLogEvalue2[maxRank] += log10(1.0e-180) * residualWeight;
        } else {
            totalLogEvalue2[maxRank] += log10(blastResult2[i].evalue) * residualWeight;
        }
    }

    cumulativeMatched1 = 0.0;
    cumulativeMatched2 = 0.0;
    cumulativeUnmatched1 = 0.0;
    cumulativeUnmatched2 = 0.0;
    cumulativetotalLogEvalue1 = 0.0;
    cumulativetotalLogEvalue2 = 0.0;

    for (i=0; i<numOfSlot; i++) {
        cumulativeMatched1 += matched1[i];
        cumulativeMatched2 += matched2[i];
        cumulativeUnmatched1 += unmatched1[i];
        cumulativeUnmatched2 += unmatched2[i];
        cumulativetotalLogEvalue1 += totalLogEvalue1[i];
        cumulativetotalLogEvalue2 += totalLogEvalue2[i];
        if (cumulativeMatched1 > 0 || cumulativeMatched2 > 0 || cumulativeUnmatched1 > 0 || cumulativeUnmatched2 > 0) {
            printf("%s", QueryName);
            if (compareByEvalue) {
                tempEvalue = pow(10, i + (int)log10(minEvalue));
                printf("\t%.0le", tempEvalue);
            } else {
                printf("\t%d", i + 1);
            }
            printf("\t%.0lf", cumulativeMatched1 + cumulativeUnmatched1);
            if (compareByEvalue) {
                printf("\t%.0lf", cumulativeUnmatched1);
            } else {
                printf("\t%.1lf", cumulativeUnmatched1);
            }
            if (averageLogEvalue) {
                if (cumulativeMatched1 + cumulativeUnmatched1 > 0) {
                    printf("\t%.0le", pow(10, cumulativetotalLogEvalue1 / (cumulativeMatched1 + cumulativeUnmatched1)));
                } else {
                    printf("\t0");
                }
            }
            printf("\t%.0lf", cumulativeMatched2 + cumulativeUnmatched2);
            if (compareByEvalue) {
                printf("\t%.0lf", cumulativeUnmatched2);
            } else {
                printf("\t%.1lf", cumulativeUnmatched2);
            }
            if (averageLogEvalue) {
                if (cumulativeMatched2 + cumulativeUnmatched2 > 0) {
                    printf("\t%.0le", pow(10, cumulativetotalLogEvalue2 / (cumulativeMatched2 + cumulativeUnmatched2)));
                } else {
                    printf("\t0");
                }
            }
            printf("\n");
        }
    }

    if (unmatchedFile1 != NULL) {
        fclose(unmatchedFile1);
    }
    if (unmatchedFile2 != NULL) {
        fclose(unmatchedFile2);
    }

    for (i=0; i<numOfDbSeq; i++) {
        free(dbSeqName[i]);
    }
    free(dbSeqName);

    free(blastResult1);
    free(blastResult2);

    iniparser_freedict(programInput);

    return 0;

}

dictionary *ParseInput(int argc, char** argv) {

    dictionary *programInput;

    char t1[3] = "-s";    // specify that this is a boolean type parameter; no following argument
    char t2[14] = "-avglogevalue";    // specify that this is a boolean type parameter; no following argument
    char *d[2];

    d[0] = t1;
    d[1] = t2;

    programInput = paraparser_load(argc, argv, 2, d);

    iniparser_copystring(programInput, "argument:1", ResultFileName1, ResultFileName1, MAX_FILENAME_LEN);
    iniparser_copystring(programInput, "argument:2", ResultFileName2, ResultFileName2, MAX_FILENAME_LEN);
    iniparser_copystring(programInput, "argument:3", QueryName, QueryName, MAX_FILENAME_LEN);
    iniparser_copystring(programInput, "argument:4", UnmatchedFileName1, UnmatchedFileName1, MAX_FILENAME_LEN);
    iniparser_copystring(programInput, "argument:5", UnmatchedFileName2, UnmatchedFileName2, MAX_FILENAME_LEN);

    minEvalue = iniparser_getdouble(programInput,"parameter:-e", minEvalue);
    maxRank = iniparser_getint(programInput,"parameter:-r", maxRank);

    if (iniparser_find_entry(programInput, "parameter:-e") && iniparser_find_entry(programInput, "parameter:-r")) {
        fprintf(stderr, "Either E-value or rank can be used!\n");
        exit(1);
    }
    if (iniparser_find_entry(programInput, "parameter:-r")) {
        compareByEvalue = FALSE;
    }

    if (iniparser_find_entry(programInput, "parameter:-s")) {
        summary = TRUE;
    }
    if (iniparser_find_entry(programInput, "parameter:-avglogevalue")) {
        averageLogEvalue = TRUE;
    }


    return programInput;

}

void ValidateParameters() {

    if (summary == FALSE) {

        if (ResultFileName1[0] == '\0' || ResultFileName1[0] == ' ' || ResultFileName1[0] == '-') {
            fprintf(stderr, "ResultFile1 is not entered!\n");
            exit(1);
        }
        if (ResultFileName2[0] == '\0' || ResultFileName2[0] == ' ' || ResultFileName2[0] == '-') {
            fprintf(stderr, "ResultFile2 is not entered!\n");
            exit(1);
        }
        if (QueryName[0] == '\0' || QueryName[0] == ' ' || QueryName[0] == '-') {
            fprintf(stderr, "QueryName is not entered!\n");
            exit(1);
        }

    }

}

