/*
  see karlin.c for details.
*/

#ifndef _karlin_
#define _karlin_

/*
wss:
The following definitions are added to remove dependency on blast.h
*/
#include <stdint.h>

extern double BlastKarlin_lambda;
extern double BlastKarlin_K;
extern double BlastKarlin_H;

void BlastKarlinBlkCalc(double* scoreProbabilities, int32_t min, int32_t max);

int32_t BlastComputeLengthAdjustment(double K,
                             double logK,
                             double alpha_d_lambda,
                             double beta,
                             int32_t query_length,
                             uint32_t db_length,
                             int32_t db_num_seqs,
                             int32_t *length_adjustment);
#endif

