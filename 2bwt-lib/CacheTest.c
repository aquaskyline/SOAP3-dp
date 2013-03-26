/* CacheTest.c        Testing speed of cache

   Copyright 2006, Wong Chi Kwong, all rights reserved.

   This module tests the speed of cache.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/

#include <stdlib.h>
#include <intrin.h>
#include <mmintrin.h>
#include "MiscUtilities.h"
#include "MemManager.h"
#include "iniparser.h"
#include "r250.h"
#include "Timing.h"

dictionary *ParseInput(int argc, char** argv);
void ValidateParameters();

    
    unsigned int MemorySize;
    unsigned int NumberOfAccess;
    char *ReadWrite;
    unsigned int NumberOfIteration;


int main(int argc, char** argv) {

    unsigned int* __restrict address;
    unsigned int* __restrict memory;

    unsigned int i, j;
    unsigned int t=0;
    unsigned int numberOfMemoryLocation;

    dictionary *programInput;
    double startTime;
    double elapsedTime = 0, totalElapsedTime = 0;
    double initializationTime, experimentTime;

    programInput = ParseInput(argc, argv);
    ValidateParameters();

    // Initialize memory manager
    MMMasterInitialize(0, 0, FALSE, NULL);

    numberOfMemoryLocation = MemorySize / sizeof(int);

    // Allocate memory
    address = MMUnitAllocate((NumberOfAccess + 8) * sizeof(unsigned int));
    memory = MMUnitAllocate((numberOfMemoryLocation + 32) * sizeof(unsigned int));

    // Set random seed
    r250_init(getRandomSeed());

    // Set start time
    startTime = setStartTime();

    printf("Initialize memory pointers with random values..");


    for (i=0; i<numberOfMemoryLocation + 32; i++) {
        memory[i] = r250();
    }
    // Initialize address and memory
    for (i=0; i<NumberOfAccess + 8; i++) {
        address[i] = (r250() % numberOfMemoryLocation) / 4 * 4;
    }

    printf("finished.\n");

    elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
    totalElapsedTime += elapsedTime;
    initializationTime = elapsedTime;

    // Test memory speed for read
    if (ReadWrite[0] == 'R' || ReadWrite[0] == 'r') {
        for (i=0; i<NumberOfIteration; i++) {
            for (j=0; j<NumberOfAccess; j++) {


                t += memory[address[j]];

            }
        }
    }
    // Test memory speed for write
    if (ReadWrite[0] == 'W' || ReadWrite[0] == 'w') {
        for (i=0; i<NumberOfIteration; i++) {
            for (j=0; j<NumberOfAccess; j++) {
                memory[address[j]] += address[j];
            }
        }
    }

    elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
    totalElapsedTime += elapsedTime;
    experimentTime = elapsedTime;

    // So that compiler does not remove code for variables t and r
    if (t==0) {
        printf("\n");
    }

    printf("Experiment completed.\n");

    printf("Initialization time   = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, initializationTime);
    printf("Experiment time       = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, experimentTime);
    
    MMUnitFree(address, (NumberOfAccess + 8) * sizeof(unsigned int));
    MMUnitFree(memory, (numberOfMemoryLocation + 32) * sizeof(unsigned int));

    iniparser_freedict(programInput);

    return 0;

}

dictionary *ParseInput(int argc, char** argv) {

    dictionary *programInput;

    programInput = paraparser_load(argc, argv, 0, NULL);

    MemorySize = iniparser_getint(programInput, "argument:1", 0);
    if (!MemorySize) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <Read/Write> <No. of iteration>\n", argv[0]);
        exit(1);
    }
    NumberOfAccess = iniparser_getint(programInput, "argument:2", 0);
    if (!NumberOfAccess) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <Read/Write> <No. of iteration>\n", argv[0]);
        exit(1);
    }
    ReadWrite = iniparser_getstring(programInput, "argument:3", 0);
    if (!NumberOfAccess) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <Read/Write> <No. of iteration>\n", argv[0]);
        exit(1);
    }
    NumberOfIteration = iniparser_getint(programInput, "argument:4", 0);
    if (!NumberOfIteration) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <Read/Write> <No. of iteration>\n", argv[0]);
        exit(1);
    }

    return programInput;

}

void ValidateParameters() {

    if (MemorySize == 0) {
        fprintf(stderr, "Memory Size = 0!\n");
        exit(1);
    }
    if (MemorySize % 4 != 0) {
        fprintf(stderr, "Memory Size must be multiple of 4!\n");
        exit(1);
    }
    if (NumberOfAccess == 0) {
        fprintf(stderr, "Number of Access = 0!\n");
        exit(1);
    }
    if (ReadWrite[0] != 'R' && ReadWrite[0] != 'r' && ReadWrite[0] != 'W' && ReadWrite[0] != 'w') {
        fprintf(stderr, "Read/Write is invalid!\n");
        exit(1);
    }
    if (NumberOfIteration == 0) {
        fprintf(stderr, "Number of Iteration = 0!\n");
        exit(1);
    }

}

