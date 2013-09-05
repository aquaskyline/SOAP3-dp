#ifndef _AIO_THREAD_H
#define _AIO_THREAD_H
#include <pthread.h>
#include "QueryParser.h"

typedef unsigned int uint;

enum ReadsFileType { SINGLE_END_TYPE, PAIR_END_TYPE, BAM_TYPE };
typedef struct InputFilePointers
{
    union
    {
        struct  //single read
        {
            gzFile queryFile;
        } single;

        struct  //pair end
        {
            gzFile queryFile;
            gzFile queryFile2;
        } pair;

        struct  //bam file
        {
            bamFile bamQueryFile;
            bam_header_t * bamHeader;
            bam1_t * bam;
        } bam;
    } filePointers;

    enum ReadsFileType fileType;
} InputFilePointers;

InputFilePointers * InputFilePointersCreate ();
void InputFilePointersSetSingle ( InputFilePointers * ifp, gzFile queryFile );
void InputFilePointersSetPair ( InputFilePointers * ifp, gzFile queryFile, gzFile queryFile2 );
void InputFilePointersSetBam ( InputFilePointers * ifp, bamFile bamQueryFile, bam_header_t * bamHeader, bam1_t * bam );
void InputFilePointersFree ( InputFilePointers * ifp );

//These status will be observed by the IO thread and  decide whether the IO thread will fill the buffer
enum BufferStatus { BUF_EMPTY, BUF_FILLING, BUF_FILLED, BUF_PROCESSING, BUF_PROCESSED, BUF_FINISHED}; //EMPTY->FILLING->FILLED->PROCESSING->PROCESSED->EMPTY
typedef struct InputReadsBuffer
{
    enum BufferStatus status;

    size_t filledNum;

    //limits
    uint maxReadLength;
    uint maxNumQueries;
    uint wordPerQuery;
    uint qualityConstant;

    //buffered data
    uint * queries;
    uint * readLengths;
    uint * readIDs;
    char * upkdQualities;
    char * upkdQueryNames;
    char isFastq;
    int maxLenReadName;
} InputReadsBuffer;

InputReadsBuffer * InputReadsBufferCreate ();
InputReadsBuffer * InputReadsBufferFullCreate ( uint maxReadLength, uint maxNumQueries, uint wordPerQuery, uint qualityConstant,
        uint * queries, uint * readLengths, uint * readIDs, char * upkdQualities, char * upkdQueryNames,
        char isFastq, int maxReadNameLen );
bool InputReadsBufferCheckStatus ( InputReadsBuffer * irb, enum BufferStatus status );
void InputReadsBufferSetStatus ( InputReadsBuffer * irb, enum BufferStatus status );
void InputReadsBufferFree ( InputReadsBuffer * irb );

void InputReadsBufferClear ( InputReadsBuffer * irb );

bool InputReadsBufferWaitForStatus ( InputReadsBuffer * irb, enum BufferStatus status );

//These status will be observed by the main thread and IO thread.
//The IO thread changes AIO_BUF_UNFILLED to AIO_BUF_FILLED, and the main thread changes AIO_BUF_FILLED to AIO_BUF_UNFILLED.
//The IO thread waits for AIO_BUF_UNFILLED, and the main thread waits for AIO_BUF_FILLED.
//when AIO_BUF_UNFILLED, the IO thread will change "bufferUnfilled" BUF_PROCESSED->BUF_EMPTY->BUF_FILLING->BUF_FILLED, then asign "bufferUnfilled"
//to "bufferFilled" and mark it AIO_BUF_UNFILLED.
//when AIO_BUF_FILLED, the main thread will change "bufferFilled" BUF_FILLED->BUF_PROCESSING->BUF_PROCESSED.
enum AIOInputBufferStatus {AIO_BUF_INIT, AIO_BUF_UNFILLED, AIO_BUF_FILLED, AIO_BUF_FINISHED};

typedef struct AIOInputBuffer
{
    int threadId;
    enum AIOInputBufferStatus status;
    InputFilePointers * reads;
    InputReadsBuffer * bufferFilled; //point to buffer0 or buffer1
    InputReadsBuffer * bufferUnFilled; //not used

    InputReadsBuffer * buffer0;
    InputReadsBuffer * buffer1;
} AIOInputBuffer;

AIOInputBuffer * AIOInputBufferCreate ( InputReadsBuffer * buffer0, InputReadsBuffer * buffer1 );
void AIOInputBufferClear ( AIOInputBuffer * aiob );
void AIOInputBufferFree ( AIOInputBuffer * aiob );

bool AIOInputBufferWaitForStatus ( AIOInputBuffer * aiob, enum AIOInputBufferStatus status );

bool AIOInputBufferCheckStatus ( AIOInputBuffer * aiob, enum AIOInputBufferStatus status );

/*
//The IO thread's main job
void AIOInputBufferFill(AIOInputBuffer *aiob);

int AIOInputThreadCreate(AIOInputBuffer *aiob);

//main thread's job, return the number of reads loaded
int LoadReadsFromAIOBuffer(AIOInputBuffer *aiob);
*/



void AIOInputBufferFill (
    AIOInputBuffer * aiob, //required
    size_t & bufferSize,    //the following args is for adapting
    uint & bufferIndex,
    unsigned char * charMap,
    char & queryChar,
    char * queryFileBuffer,
    size_t & bufferSize2,
    uint & bufferIndex2,
    char & queryChar2,
    char * queryFileBuffer2
);

void AIOInputThreadCreate (
    AIOInputBuffer * aiob, //required
    size_t & bufferSize,    //the following args is for adapting
    uint & bufferIndex,
    unsigned char * charMap,
    char & queryChar,
    char * queryFileBuffer,
    size_t & bufferSize2,
    uint & bufferIndex2,
    char & queryChar2,
    char * queryFileBuffer2
);


// This function resets the status of buffer with processed reads.
// It should be called before calling LoadReadsFromAIOBuffer (except first calling).
void ResetBufferStatusToUnfilled ( AIOInputBuffer * aiob );

InputReadsBuffer * LoadReadsFromAIOBuffer (
    AIOInputBuffer * aiob //required
);





#endif
