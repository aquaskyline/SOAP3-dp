/* ==========================================================================
 *
 * HSPstatistic.h
 *
 * author: Wong Swee Seong
 * date: 30-03-2006
 *
 * gapped and ungapped extension for string matching
 * 
 * ==========================================================================*/


#ifndef _HSP_statistic_
#define _HSP_statistic_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define DNA_CHAR_SIZE 16

typedef struct stat__score_block
{
  int match, mismatch, gapOpen, gapExt;
} STAT_scoreBlock;

typedef struct stat__Xdropoffs
{
  double ungapXdropoff, gapXdropoff, gapXdropoffFinal;
} STAT_Xdropoffs;

void initializeHSPstatistic(uint64_t databaseSize, uint64_t numOfSequences, uint64_t minSeqLength, double *dbCharProb,
                            uint64_t querySize, uint64_t numOfContext, double *queryCharProb, int scoringMatrix[DNA_CHAR_SIZE][DNA_CHAR_SIZE]);
void printHSPstatistic(FILE * outFile);
int32_t calcUngapCutoffScore();
int32_t calcGapCutoffScore();
int32_t getUngapXdropoff();
int32_t getGapXdropoff();
int32_t getGapXdropoffFinal();

int32_t stat_ungapNormalized2nominal(double normalizedScore);
double stat_ungapNominal2normalized(int32_t nominalScore);
double stat_ungapCalcEvalue(double normalizedScore);
int32_t stat_ungapEvalue2nominal(double evalue);
double stat_gapNominal2normalized(int32_t nominalScore);
int32_t stat_gapNormalized2nominal(double normalizedScore);
double stat_gapCalcEvalue(double normalizedScore);
int32_t stat_gapEvalue2nominal(double evalue);
 
#ifdef __cplusplus
}
#endif

#endif
