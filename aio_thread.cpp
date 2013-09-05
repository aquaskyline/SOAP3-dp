#include "aio_thread.h"
#include "CPUfunctions.h"

//InputFilePointers
InputFilePointers * InputFilePointersCreate ()
{
    InputFilePointers * t = ( InputFilePointers * ) calloc ( 1, sizeof ( InputFilePointers ) );

    if ( !t )
    {
        fprintf ( stderr, "InputFilePointersCreate failed!\n" );
        exit ( -1 );
    }

    return t;
};

void InputFilePointersSetSingle ( InputFilePointers * ifp, gzFile queryFile )
{
    if ( !ifp ) { return; }

    ifp->fileType = SINGLE_END_TYPE;
    ifp->filePointers.single.queryFile = queryFile;
};
void InputFilePointersSetPair ( InputFilePointers * ifp, gzFile queryFile, gzFile queryFile2 )
{
    if ( !ifp ) { return; }

    ifp->fileType = PAIR_END_TYPE;
    ifp->filePointers.pair.queryFile = queryFile;
    ifp->filePointers.pair.queryFile2 = queryFile2;
};

void InputFilePointersSetBam ( InputFilePointers * ifp, bamFile bamQueryFile, bam_header_t * bamHeader, bam1_t * bam )
{
    if ( !ifp ) { return; }

    ifp->fileType = BAM_TYPE;
    ifp->filePointers.bam.bamQueryFile = bamQueryFile;
    ifp->filePointers.bam.bamHeader = bamHeader;
    ifp->filePointers.bam.bam = bam;
};
void InputFilePointersFree ( InputFilePointers * ifp )
{
    if ( !ifp ) { return; }

    free ( ( void * ) ifp );
};


//InputReadsBuffer
InputReadsBuffer * InputReadsBufferCreate ()
{
    InputReadsBuffer * t = ( InputReadsBuffer * ) calloc ( 1, sizeof ( InputReadsBuffer ) );

    if ( !t )
    {
        fprintf ( stderr, "InputReadsBufferCreate failed!\n" );
        exit ( -1 );
    }

    t->status = BUF_EMPTY;
    t->filledNum = 0;
    return t;
};

InputReadsBuffer * InputReadsBufferFullCreate ( uint maxReadLength, uint maxNumQueries, uint wordPerQuery, uint qualityConstant,
        uint * queries, uint * readLengths, uint * readIDs, char * upkdQualities, char * upkdQueryNames,
        char isFastq, int maxLenReadName )
{
    InputReadsBuffer * t = InputReadsBufferCreate ();
    t->maxReadLength = maxReadLength;
    t->maxNumQueries = maxNumQueries;
    t->wordPerQuery = wordPerQuery;
    t->qualityConstant = qualityConstant;
    t->queries = queries;
    t->readLengths = readLengths;
    t->readIDs = readIDs;
    t->upkdQualities = upkdQualities;
    t->upkdQueryNames = upkdQueryNames;
    t->isFastq = isFastq;
    t->maxLenReadName = maxLenReadName;
    return t;
};
bool InputReadsBufferCheckStatus ( InputReadsBuffer * irb, enum BufferStatus status )
{
    if ( !irb )
    {
        fprintf ( stderr, "NULL pointer exception in InputReadsBufferCheckStatus().\n" );
        exit ( -1 );
    }

    return irb->status == status ? 1 : 0;
};
void InputReadsBufferSetStatus ( InputReadsBuffer * irb, enum BufferStatus status )
{
    if ( !irb ) { return; }

    irb->status = status;
};
void InputReadsBufferFree ( InputReadsBuffer * irb )
{
    if ( !irb ) { return; }

    free ( ( void * ) irb );
};
void InputReadsBufferClear ( InputReadsBuffer * irb )
{
    if ( !irb ) { return; }

    irb->status = BUF_EMPTY;
    irb->filledNum = 0;
};

bool InputReadsBufferWaitForStatus ( InputReadsBuffer * irb, enum BufferStatus status )
{
    if ( !irb )
    {
        fprintf ( stderr, "NULL pointer exception in InputReadsBufferWaitForStatus().\n" );
        exit ( -1 );
    }

    while ( 1 )
    {
        if ( irb->status == status )
        {
            return 1;
        }
        else
        {
            sleep ( 1 );
        }
    }
};


//AIOInputBuffer
AIOInputBuffer * AIOInputBufferCreate ( InputReadsBuffer * buffer0, InputReadsBuffer * buffer1 )
{
    AIOInputBuffer * t = ( AIOInputBuffer * ) calloc ( 1, sizeof ( AIOInputBuffer ) );
    t->status = AIO_BUF_INIT;
    t->buffer0 = buffer0;
    t->buffer1 = buffer1;
    t->bufferUnFilled = t->buffer0;
    return t;
};
void AIOInputBufferClear ( AIOInputBuffer * aiob )
{
    if ( !aiob )
    {
        fprintf ( stderr, "NULL pointer exception in AIOInputBufferClear().\n" );
        exit ( -1 );
    }

    if ( aiob->buffer0 )
    {
        InputReadsBufferClear ( aiob->buffer0 );
    }

    if ( aiob->buffer1 )
    {
        InputReadsBufferClear ( aiob->buffer1 );
    }

    aiob->status = AIO_BUF_INIT;
};
void AIOInputBufferFree ( AIOInputBuffer * aiob )
{
    if ( aiob )
    {
        InputFilePointersFree ( aiob->reads );
        InputReadsBufferFree ( aiob->buffer0 );
        InputReadsBufferFree ( aiob->buffer1 );
        free ( ( void * ) aiob );
    }
};

bool AIOInputBufferWaitForStatus ( AIOInputBuffer * aiob, enum AIOInputBufferStatus status )
{
    if ( !aiob )
    {
        fprintf ( stderr, "NULL pointer exception in AIOInputBufferWaitForStatus().\n" );
        exit ( -1 );
    }

    while ( 1 )
    {
        if ( aiob->status == status )
        {
            return 1;
        }
        else
        {
            sleep ( 1 );
        }
    }
};


bool AIOInputBufferCheckStatus ( AIOInputBuffer * aiob, enum AIOInputBufferStatus status )
{
    if ( !aiob )
    {
        fprintf ( stderr, "NULL pointer exception in AIOInputBufferCheckStatus().\n" );
        exit ( -1 );
    }

    return aiob->status == status ? 1 : 0;
};

//core
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
)
{
    if ( !aiob )
    {
        fprintf ( stderr, "NULL pointer exception in AIOInputBufferFill().\n" );
        exit ( -1 );
    }

#ifdef DEBUG_AIO_THREAD
    static size_t times0, times1;
#endif

    if ( AIOInputBufferCheckStatus ( aiob, AIO_BUF_INIT ) )
    {
        aiob->status = AIO_BUF_UNFILLED;
    }

    while ( 1 )
    {
        InputReadsBuffer * bufferUnFilled = aiob->buffer0; //check buffer0

        if ( InputReadsBufferCheckStatus ( bufferUnFilled, BUF_PROCESSED ) )
        {
            InputReadsBufferClear ( bufferUnFilled );
        }

        if ( InputReadsBufferCheckStatus ( bufferUnFilled, BUF_EMPTY ) ) // begin to fill buffer
        {
#ifdef DEBUG_AIO_THREAD
            fprintf ( stderr, "Filling buffer0 times:%llu\n", ++times0 );
#endif
            InputReadsBufferSetStatus ( bufferUnFilled, BUF_FILLING );
            int num = 0;

            switch ( aiob->reads->fileType )
            {
                case SINGLE_END_TYPE:
                    //call loadSingleReadsGz
                {
                    num = loadSingleReadsGz (
                              aiob->reads->filePointers.single.queryFile,
                              queryFileBuffer,
                              charMap,
                              bufferUnFilled->queries,
                              bufferUnFilled->readLengths,
                              bufferUnFilled->readIDs,
                              bufferUnFilled->upkdQualities,
                              bufferUnFilled->upkdQueryNames,
                              bufferUnFilled->maxReadLength,
                              bufferUnFilled->maxNumQueries,
                              bufferSize,
                              queryChar,
                              bufferIndex,
                              0,//    accumReadNum,
                              bufferUnFilled->wordPerQuery,
                              bufferUnFilled->qualityConstant,
                              bufferUnFilled->isFastq,
                              bufferUnFilled->maxLenReadName );
                    bufferUnFilled->filledNum = num;
                    break;
                }

                case PAIR_END_TYPE:
                    //call loadPairReadsGz2
                {
                    num = loadPairReadsGz2 (
                              aiob->reads->filePointers.pair.queryFile,
                              queryFileBuffer,
                              aiob->reads->filePointers.pair.queryFile2,
                              queryFileBuffer2,
                              charMap,
                              bufferUnFilled->queries,
                              bufferUnFilled->readLengths,
                              bufferUnFilled->readIDs,
                              bufferUnFilled->upkdQualities,
                              bufferUnFilled->upkdQueryNames,
                              bufferUnFilled->maxReadLength,
                              bufferUnFilled->maxNumQueries,
                              bufferSize,
                              queryChar,
                              bufferIndex,
                              bufferSize2,
                              queryChar2,
                              bufferIndex2,
                              0,//accumReadNum,
                              bufferUnFilled->wordPerQuery,
                              bufferUnFilled->qualityConstant,
                              bufferUnFilled->isFastq,
                              bufferUnFilled->maxLenReadName );
                    bufferUnFilled->filledNum = num;
                    break;
                }

                case BAM_TYPE:
                    //call loadBAMReads
                {
                    num = loadBAMReads (
                              aiob->reads->filePointers.bam.bamQueryFile,
                              aiob->reads->filePointers.bam.bamHeader,
                              aiob->reads->filePointers.bam.bam,
                              charMap,
                              bufferUnFilled->queries,
                              bufferUnFilled->readLengths,
                              bufferUnFilled->readIDs,
                              bufferUnFilled->upkdQualities,
                              bufferUnFilled->upkdQueryNames,
                              bufferUnFilled->maxReadLength,
                              bufferUnFilled->maxNumQueries,
                              bufferUnFilled->wordPerQuery,
                              bufferUnFilled->qualityConstant,
                              bufferUnFilled->maxLenReadName );
                    bufferUnFilled->filledNum = num;
                    bufferUnFilled->isFastq = 1;
                    break;
                }

                default:
                    fprintf ( stderr, "ERROR fileType!\n" );
            }

            if ( bufferUnFilled->filledNum > 0 )
            {
                InputReadsBufferSetStatus ( bufferUnFilled, BUF_FILLED );
                AIOInputBufferWaitForStatus ( aiob, AIO_BUF_UNFILLED ); //wait for unfilled
                aiob->bufferFilled = bufferUnFilled;
                aiob->status = AIO_BUF_FILLED;
            }
            else
            {
                InputReadsBufferSetStatus ( bufferUnFilled, BUF_FINISHED );
                AIOInputBufferWaitForStatus ( aiob, AIO_BUF_UNFILLED ); //wait for unfilled
                aiob->bufferFilled = bufferUnFilled;
                aiob->status = AIO_BUF_FINISHED;
                return;
            }
        }

        bufferUnFilled = aiob->buffer1; //check buffer1

        if ( InputReadsBufferCheckStatus ( bufferUnFilled, BUF_PROCESSED ) )
        {
            InputReadsBufferClear ( bufferUnFilled );
        }

        if ( InputReadsBufferCheckStatus ( bufferUnFilled, BUF_EMPTY ) ) // begin to fill buffer
        {
#ifdef DEBUG_AIO_THREAD
            fprintf ( stderr, "Filling buffer1 times:%llu\n", ++times1 );
#endif
            InputReadsBufferSetStatus ( bufferUnFilled, BUF_FILLING );
            int num = 0;

            switch ( aiob->reads->fileType )
            {
                case SINGLE_END_TYPE:
                    //call loadSingleReadsGz
                {
                    num = loadSingleReadsGz (
                              aiob->reads->filePointers.single.queryFile,
                              queryFileBuffer,
                              charMap,
                              bufferUnFilled->queries,
                              bufferUnFilled->readLengths,
                              bufferUnFilled->readIDs,
                              bufferUnFilled->upkdQualities,
                              bufferUnFilled->upkdQueryNames,
                              bufferUnFilled->maxReadLength,
                              bufferUnFilled->maxNumQueries,
                              bufferSize,
                              queryChar,
                              bufferIndex,
                              0,//    accumReadNum,
                              bufferUnFilled->wordPerQuery,
                              bufferUnFilled->qualityConstant,
                              bufferUnFilled->isFastq,
                              bufferUnFilled->maxLenReadName );
                    bufferUnFilled->filledNum = num;
                    break;
                }

                case PAIR_END_TYPE:
                    //call loadPairReadsGz2
                {
                    num = loadPairReadsGz2 (
                              aiob->reads->filePointers.pair.queryFile,
                              queryFileBuffer,
                              aiob->reads->filePointers.pair.queryFile2,
                              queryFileBuffer2,
                              charMap,
                              bufferUnFilled->queries,
                              bufferUnFilled->readLengths,
                              bufferUnFilled->readIDs,
                              bufferUnFilled->upkdQualities,
                              bufferUnFilled->upkdQueryNames,
                              bufferUnFilled->maxReadLength,
                              bufferUnFilled->maxNumQueries,
                              bufferSize,
                              queryChar,
                              bufferIndex,
                              bufferSize2,
                              queryChar2,
                              bufferIndex2,
                              0,//accumReadNum,
                              bufferUnFilled->wordPerQuery,
                              bufferUnFilled->qualityConstant,
                              bufferUnFilled->isFastq,
                              bufferUnFilled->maxLenReadName );
                    bufferUnFilled->filledNum = num;
                    break;
                }

                case BAM_TYPE:
                    //call loadBAMReads
                {
                    num = loadBAMReads (
                              aiob->reads->filePointers.bam.bamQueryFile,
                              aiob->reads->filePointers.bam.bamHeader,
                              aiob->reads->filePointers.bam.bam,
                              charMap,
                              bufferUnFilled->queries,
                              bufferUnFilled->readLengths,
                              bufferUnFilled->readIDs,
                              bufferUnFilled->upkdQualities,
                              bufferUnFilled->upkdQueryNames,
                              bufferUnFilled->maxReadLength,
                              bufferUnFilled->maxNumQueries,
                              bufferUnFilled->wordPerQuery,
                              bufferUnFilled->qualityConstant,
                              bufferUnFilled->maxLenReadName );
                    bufferUnFilled->filledNum = num;
                    bufferUnFilled->isFastq = 1;
                    break;
                }

                default:
                    fprintf ( stderr, "ERROR fileType!\n" );
            }

            if ( bufferUnFilled->filledNum > 0 )
            {
                InputReadsBufferSetStatus ( bufferUnFilled, BUF_FILLED );
                AIOInputBufferWaitForStatus ( aiob, AIO_BUF_UNFILLED ); //wait for unfilled
                aiob->bufferFilled = bufferUnFilled;
                aiob->status = AIO_BUF_FILLED;
            }
            else
            {
                InputReadsBufferSetStatus ( bufferUnFilled, BUF_FINISHED );
                AIOInputBufferWaitForStatus ( aiob, AIO_BUF_UNFILLED ); //wait for unfilled
                aiob->bufferFilled = bufferUnFilled;
                aiob->status = AIO_BUF_FINISHED;
                return;
            }
        }
    }
};


void ResetBufferStatusToUnfilled ( AIOInputBuffer * aiob )
{
    if ( !aiob )
    {
        fprintf ( stderr, "NULL pointer exception in LoadReadsFromAIOBuffer().\n" );
        exit ( -1 );
    }

    //When this function is called, it means that the reads loaded before have been processed.
    //so should mark the buffer previous used BUF_PROCESSED, mark current buffer UNFILLED

    if ( aiob->bufferFilled ) //at the begining, bufferFilled is null.
    {
        if ( InputReadsBufferCheckStatus ( aiob->bufferFilled, BUF_PROCESSING ) )
        {
            InputReadsBufferSetStatus ( aiob->bufferFilled, BUF_PROCESSED );
        }

        if ( aiob->status == AIO_BUF_FILLED )
        {
            aiob->status = AIO_BUF_UNFILLED;
        }
    }
}

InputReadsBuffer * LoadReadsFromAIOBuffer ( AIOInputBuffer * aiob )
{
    /*
       if(AIOInputBufferWaitForStatus(aiob,AIO_BUF_FILLED)){//should be changed to wait AIO_BUF_FILLED or AIO_BUF_FINISHED
       return aiob->bufferFilled;
       }*/
#ifdef DEBUG_AIO_THREAD
    static  size_t times;
    times++;
    fprintf ( stderr, "Enter LoadReadsFromAIOBuffer %llu times\n", times );
#endif

    if ( !aiob )
    {
        fprintf ( stderr, "NULL pointer exception in LoadReadsFromAIOBuffer().\n" );
        exit ( -1 );
    }

    while ( 1 ) //wait AIO_BUF_FILLED or AIO_BUF_FINISHED
    {
        if ( aiob->status == AIO_BUF_FILLED )
        {
            if ( aiob->bufferFilled->status == BUF_FILLED )
            {
                InputReadsBufferSetStatus ( aiob->bufferFilled, BUF_PROCESSING );
            }
            else
            {
                fprintf ( stderr, "ERROR LoadReadsFromAIOBuffer\n" );
            }

#ifdef DEBUG_AIO_THREAD
            fprintf ( stderr, "Exit LoadReadsFromAIOBuffer %llu times\n", times );
#endif
            return aiob->bufferFilled;
        }
        else if ( aiob->status == AIO_BUF_FINISHED )
        {
#ifdef DEBUG_AIO_THREAD
            fprintf ( stderr, "Exit LoadReadsFromAIOBuffer %llu times\n", times );
#endif
            return aiob->bufferFilled;
        }
        else
        {
            sleep ( 1 );
#ifdef DEBUG_AIO_THREAD

            if ( !aiob->bufferFilled )
            {
                fprintf ( stderr, "Sleep status %llu %d \n", times, aiob->status );
            }
            else
            {
                fprintf ( stderr, "Sleep status %llu %d %d\n", times, aiob->status, aiob->bufferFilled->status );
            }

#endif
        }
    }
};

typedef struct IOParams
{
    AIOInputBuffer * aiob;
    size_t * bufferSize;
    uint * bufferIndex;
    unsigned char * charMap;
    char * queryChar;
    char * queryFileBuffer;
    size_t * bufferSize2;
    uint * bufferIndex2;
    char * queryChar2;
    char * queryFileBuffer2;
} IOParams;

void * AIOThreadRun ( void * args )
{
    if ( !args )
    {
        fprintf ( stderr, "AIOThreadRun args NULL\n" );
        exit ( -1 );
    }

    IOParams * param = ( IOParams * ) args;

    if ( !param->aiob )
    {
        fprintf ( stderr, "aiob NULL AIOThreadRun()\n" );
        exit ( -1 );
    }

    AIOInputBufferFill (
        param->aiob,
        * ( param->bufferSize ),
        * ( param->bufferIndex ),
        param->charMap,
        * ( param->queryChar ),
        param->queryFileBuffer,
        * ( param->bufferSize2 ),
        * ( param->bufferIndex2 ),
        * ( param->queryChar2 ),
        param->queryFileBuffer2
    );
    free ( ( void * ) args );
    return NULL;
}
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
)
{
    if ( !aiob )
    {
        fprintf ( stderr, "aiob NULL AIOInputThreadCreate()\n" );
        exit ( -1 );
    }

    IOParams * param = ( IOParams * ) malloc ( sizeof ( IOParams ) );
    param->aiob = aiob;
    param->bufferSize = &bufferSize;
    param->bufferIndex = &bufferIndex;
    param->charMap = charMap;
    param->queryChar = &queryChar;
    param->queryFileBuffer = queryFileBuffer;
    param->bufferSize2 = &bufferSize2;
    param->bufferIndex2 = &bufferIndex2;
    param->queryChar2 = &queryChar2;
    param->queryFileBuffer2 = queryFileBuffer2;
    pthread_t thread_id;
    int t;
    t = pthread_create ( &thread_id, NULL, AIOThreadRun, param );

    // if(!t){
    //  fprintf(stderr,"Create an IO thread.\n");
    // }
    if ( t )
    {
        fprintf ( stderr, "Error in Creating the IO thread.\n" );
    }
};
















