/*

   2BWT-Benchmark.c         SAMPLE CODE FOR 2BWT-LIB
   
   Copyright (C) 2011, Edward Wu.
   Website: http://www.cs.hku.hk/~mkewu/2bwt
   
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
   
   Date   : 7th May 2013
   Author : Edward MK Wu
   Change : Sample file to benchmark 2BWT index.
   
            For complete guide on what the 2BWT index library provided
            please check out the file 2BWT-Interface.h for the description
            of all library calls.

*/

////////////////////////////////////
// DEBUG FLAG
////////////////////////////////////
// Define below item to generate a list of 
// index everytime the program is run.
//
// #define DEBUG_GENERATE_LIST 
////////////////////////////////////
// Define below item to generate the output of the BWTDecode
// for manual correctness check capturing the output
//
//#define DEBUG_GENERATE_OUTPUT 
////////////////////////////////////
// Define below item to perform BWTDecode Testing
//
#define DEBUG_PERFORM_DECODE 
////////////////////////////////////
// Define below item to perform BWTDecodeALL Testing
//
#define DEBUG_PERFORM_DECODE_ALL 
////////////////////////////////////

#define BM_NUM_OF_THREAD       4
#define BM_TEST_CHARACTER_CODE 3

#include <stdio.h>
#include <stdlib.h>
#include "../Timing.h"
#include "../2BWT-Interface.h"
#include <time.h>
#include <pthread.h>

typedef struct ThreadArg {
   unsigned int * bucketBody;
   unsigned int bucketStartIdx;
   unsigned int bucketSize;
   BWT * bwt;
} ThreadArg;

void Generate(unsigned int bwtLength) {
    time_t sec;
    unsigned long long i;
    sec = time (NULL);
    srand(sec);
    
    char str[256];
    sprintf(str, "2BWT-Benchmark.Gen-%u",(unsigned int) sec);
    FILE * fp = fopen(str,"w");
    
    unsigned long long testBwtIndexesCount = 244*1000000;
    for (i=0;i<testBwtIndexesCount;i++) {
        fprintf(fp,"%u\n",rand() % bwtLength);
    }

    fclose(fp);
}

void *ThreadFuncAll ( void * _arg ) {
    ThreadArg * arg = (ThreadArg*) _arg;
    unsigned int k;
    unsigned int ALIGN_64 occValue[ALPHABET_SIZE];
    unsigned int startIdx = arg->bucketStartIdx;

    BWT * bwt = arg->bwt;
    unsigned int * body = arg->bucketBody;
    for (k=startIdx;k<startIdx+arg->bucketSize;k++) {
        BWTAllOccValue(bwt,body[k],occValue);
        #ifdef DEBUG_GENERATE_OUTPUT
            printf("DEBUG ALLPARA %llu\t%u\t%u\n",k,body[k],occValue[BM_TEST_CHARACTER_CODE]);
        #endif
    }
}

void *ThreadFunc ( void * _arg ) {
    ThreadArg * arg = (ThreadArg*) _arg;
    unsigned int k;
    unsigned int startIdx = arg->bucketStartIdx;

    BWT * bwt = arg->bwt;
    unsigned int * body = arg->bucketBody;
    for (k=startIdx;k<startIdx+arg->bucketSize;k++) {
        unsigned int result = BWTOccValue(bwt,body[k],BM_TEST_CHARACTER_CODE);
        #ifdef DEBUG_GENERATE_OUTPUT
            printf("DEBUG PARAL %llu\t%u\t%u\n",k,body[k],result);
        #endif
    }
}

int main(int argc,char ** argv) {
    double startTime, lastEventTime, timestamp;
    int j,c;
    unsigned long long i,k,offset;
    ThreadArg threadArgs[BM_NUM_OF_THREAD];
    pthread_t threads[BM_NUM_OF_THREAD];

    /*__m512i a;
    __m512 b;
    a = _mm512_set1_epi32(0xFFFFFFFE);    // Character selection mask for even bits
    b = a;
    unsigned long long ab[8];
    _mm512_storenrngo_ps(ab,b);
    for (i=0;i<8;i++) {printf("%llu %u\n",ab[i]>>1,_mm_countbits_64(ab[i]));}*/

    //Variables for backward and forward search
    unsigned int l,r,rev_l,rev_r;
    
    //Variables for search all sa ranges functions
    unsigned int result_l[ALPHABET_SIZE];
    unsigned int result_r[ALPHABET_SIZE];
    unsigned int result_rev_l[ALPHABET_SIZE];
    unsigned int result_rev_r[ALPHABET_SIZE];
    
    //Variables for result
    int sequenceId;
    unsigned int saCount;
    
    // Load up the index with the below statement
    printf("Loading index ... "); 
    fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT("ncbi.genome1-5.fa.index",".sa");
    printf("DONE\n\n"); 

    
    // Convert the pattern into 2BWT recognised coding scheme
    BWT * bwt = idx2BWT->bwt;

    #ifdef DEBUG_GENERATE_LIST
        Generate(bwt->textLength);
        return 0;
    #endif
    
    startTime = setStartTime();
    lastEventTime = startTime;
    
    // -- Settings for BM_NUM_OF_THREAD = 200;
    unsigned long long testBwtIndexesCount = 10000000;
    unsigned int bucket = testBwtIndexesCount/BM_NUM_OF_THREAD;
    
    // -- Settings for BM_NUM_OF_THREAD = 244;
    //unsigned long long testBwtIndexesCount = 24400000;
    //unsigned int bucket = testBwtIndexesCount/BM_NUM_OF_THREAD;
    
    // -- Settings for BM_NUM_OF_THREAD = 244;
    //unsigned long long testBwtIndexesCount = 244*1000000;
    //unsigned int bucket = testBwtIndexesCount/BM_NUM_OF_THREAD;
    
    FILE * fp;
    if (argc>1) {
        fp = fopen(argv[1],"r");
    }
    unsigned int * testBwtIndexes = (unsigned int*) malloc(sizeof(unsigned int) * testBwtIndexesCount);
    for (i=0;i<testBwtIndexesCount;i++) {
        if (argc>1) {
            fscanf(fp,"%u",&testBwtIndexes[i]);
        } else {
            testBwtIndexes[i] = rand() % bwt->textLength;
        }
    }
    if (argc>1) {
        fclose(fp);
    }
    
    timestamp = getElapsedTime(startTime);
    printf("(Elapsed time (Generate Indexes) : %9.4f seconds)\n\n", timestamp);
    lastEventTime = timestamp;

    printf("Test Data Size = %llu\n",testBwtIndexesCount);
    
#ifdef DEBUG_PERFORM_DECODE
    ////////////////////////////////////////////////////////
    // Version of software that split the buffer to MIC cores
    ////////////////////////////////////////////////////////
    printf("MIC Parallelism %u\n",BM_NUM_OF_THREAD);
    for (j=0;j<BM_NUM_OF_THREAD;j++) {
        threadArgs[j].bwt=bwt;
        threadArgs[j].bucketStartIdx=j*bucket;
        threadArgs[j].bucketBody=testBwtIndexes;
        threadArgs[j].bucketSize=bucket;
    }
    #ifndef DEBUG_GENERATE_OUTPUT
    for (j=0;j<5;j++) {
    #endif
        for (i=0; i<BM_NUM_OF_THREAD; i++){
            if (pthread_create(&(threads[i]), NULL, ThreadFunc, (void*) &(threadArgs[i])))
                fprintf(stderr,"Can't create thread %d!\n",i), exit(1);
        }
        for (i=0; i<BM_NUM_OF_THREAD; i++){
            if (threads[i]!=0 ) {
                if (pthread_join(threads[i],NULL))
                    fprintf(stderr,"Thread %d crashed!\n",i), exit(1);
                threads[i]=0;
            }
        }
        timestamp = getElapsedTime(startTime);
        printf("(Elapsed time (Search Index) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        lastEventTime = timestamp;
    #ifndef DEBUG_GENERATE_OUTPUT
    }
    #endif
    
    ////////////////////////////////////////////////////////
    // Version of software that does it with one core only
    ////////////////////////////////////////////////////////
    printf("Single Core\n");
    #ifndef DEBUG_GENERATE_OUTPUT
    for (j=0;j<5;j++) {
    #endif
        for (i=0;i<testBwtIndexesCount;i++) {
            unsigned int result = BWTOccValue(bwt,testBwtIndexes[i],BM_TEST_CHARACTER_CODE);
            #ifdef DEBUG_GENERATE_OUTPUT
            printf("DEBUG SINGLE %llu\t%u\t%u\n",i,testBwtIndexes[i],result);
            #endif
        }
        timestamp = getElapsedTime(startTime);
        printf("(Elapsed time (Search Index) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        lastEventTime = timestamp;
    #ifndef DEBUG_GENERATE_OUTPUT
    }
    #endif
#endif
    
    
#ifdef DEBUG_PERFORM_DECODE_ALL
    ////////////////////////////////////////////////////////
    // Version of software that split the buffer to MIC cores
    ////////////////////////////////////////////////////////
    printf("MIC ALL Parallelism %u\n",BM_NUM_OF_THREAD);
    #ifndef DEBUG_GENERATE_OUTPUT
    for (j=0;j<5;j++) {
    #endif
        for (i=0; i<BM_NUM_OF_THREAD; i++){
            if (pthread_create(&(threads[i]), NULL, ThreadFuncAll, (void*) &(threadArgs[i])))
                fprintf(stderr,"Can't create thread %d!\n",i), exit(1);
        }
        for (i=0; i<BM_NUM_OF_THREAD; i++){
            if (threads[i]!=0 ) {
                if (pthread_join(threads[i],NULL))
                    fprintf(stderr,"Thread %d crashed!\n",i), exit(1);
                threads[i]=0;
            }
        }
        timestamp = getElapsedTime(startTime);
        printf("(Elapsed time (Search Index) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        lastEventTime = timestamp;
    #ifndef DEBUG_GENERATE_OUTPUT
    }
    #endif
    
    ////////////////////////////////////////////////////////
    // Version of software that does it with one core only
    ////////////////////////////////////////////////////////
    printf("Single ALL Core\n");
    #ifndef DEBUG_GENERATE_OUTPUT
    for (j=0;j<5;j++) {
    #endif
        unsigned int ALIGN_64 occValue[ALPHABET_SIZE];
        for (i=0;i<testBwtIndexesCount;i++) {
            BWTAllOccValue(bwt,testBwtIndexes[i],occValue);
            #ifdef DEBUG_GENERATE_OUTPUT
                printf("DEBUG ALLSINGLE %llu\t%u\t%u\n",i,testBwtIndexes[i],occValue[BM_TEST_CHARACTER_CODE]);
            #endif
        }
        timestamp = getElapsedTime(startTime);
        printf("(Elapsed time (Search Index) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        lastEventTime = timestamp;
    #ifndef DEBUG_GENERATE_OUTPUT
    }
    #endif
#endif


    // Free up the 2BWT index
    printf("\nFree index ... "); 
    fflush(stdout);
    BWTFree2BWT(idx2BWT);
    printf("DONE\n"); 
    free(testBwtIndexes);
    return 0;
}
