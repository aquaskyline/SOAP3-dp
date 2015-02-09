/* SSETest.c        Testing speed of SSE decode of BWT

   Copyright 2006, Wong Chi Kwong, all rights reserved.

   This module tests the speed of SSE decode of BWT.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/

#include <stdlib.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include "MiscUtilities.h"
#include "MemManager.h"
#include "iniparser.h"
#include "r250.h"
#include "Timing.h"
#include "DNACount.h"

dictionary *ParseInput(int argc, char** argv);
void ValidateParameters();

    
    unsigned int MemorySize;
    unsigned int NumberOfAccess;
    char *Mode;
    unsigned int NumberOfIteration;


void GenerateDNAOccCountTable1(unsigned int *dnaDecodeTable) {

    unsigned int i, j, c, t;

    for (i=0; i<65536; i++) {
        dnaDecodeTable[i] = 0;
        c = i;
        for (j=0; j<8; j++) {
            t = c & 0x00000003;
            dnaDecodeTable[i] += 1 << (t * 8);
            c >>= 2;
        }
    }

}

int main(int argc, char** argv) {

    #define CHAR_PER_128 64
    #define CHAR_PER_256 128
    #define CHAR_PER_WORD 16


    unsigned int* __restrict address;
    unsigned int* __restrict memory;
    unsigned int* __restrict dnaDecodeTable;
    unsigned int* __restrict index1;
    unsigned int* __restrict index2;
    unsigned int* __restrict character;

    unsigned int i, j;
    unsigned int t=0, r=0, r1=0;
    unsigned int temp;
    unsigned int t4[4] = { 0, 0, 0, 0 };
    unsigned int r4[4] = { 0, 0, 0, 0 };
    unsigned int numberOfMemoryLocation;

    dictionary *programInput;
    double startTime;
    double elapsedTime = 0, totalElapsedTime = 0;
    double initializationTime, experimentTime;

    unsigned int numChar1, numChar2, minIndex, maxIndex, minIndex128, maxIndex128;

    const static unsigned int ALIGN_16 partitionOne1[4]  = { 47, 31, 15, 0 };
    const static unsigned int ALIGN_16 partitionOne2[4]  = { 0, 15, 31, 47 };
    const static unsigned int ALIGN_16 partitionZero1[4]  = { 63, 47, 31, 15 };
    const static unsigned int ALIGN_16 partitionZero2[4]  = { 15, 31, 47, 63 };

    // SSE registers
    __m128i r1e, r2e;
    __m128i mcl;

    // Count one alphabet only
    __m128i m0, m1;
    __m128i r1a, r1b, r1c;
    __m128i r2a, r2b, r2c;

    // Count all alphabets only
    __m128i rc, rg, rt;
    __m128i ra1, ra2;
    __m128i rc1, rc2;
    __m128i rg1, rg2;
    __m128i rt1, rt2;

    programInput = ParseInput(argc, argv);
    ValidateParameters();

    // Initialize memory manager
    MMMasterInitialize(0, 0, FALSE, NULL);

    numberOfMemoryLocation = MemorySize / sizeof(int);

    // Allocate memory
    address = MMUnitAllocate((NumberOfAccess + 8) * sizeof(unsigned int));
    memory = MMUnitAllocate((numberOfMemoryLocation + 32) * sizeof(unsigned int));
    index1 = MMUnitAllocate((NumberOfAccess + 8) * sizeof(unsigned int));
    index2 = MMUnitAllocate((NumberOfAccess + 8) * sizeof(unsigned int));
    character = MMUnitAllocate((NumberOfAccess + 8) * sizeof(unsigned int));

    dnaDecodeTable = MMUnitAllocate(65536 * sizeof(unsigned int));

    // Set random seed
    r250_init(getRandomSeed());
//    r250_init(1);

    // Set start time
    startTime = setStartTime();

    printf("Initialize memory pointers with random values..");

    GenerateDNAOccCountTable1(dnaDecodeTable);


    for (i=0; i<numberOfMemoryLocation + 32; i++) {
        memory[i] = r250();
    }
    // Initialize address and memory
    for (i=0; i<NumberOfAccess + 8; i++) {
        address[i] = (r250() % numberOfMemoryLocation) / 8 * 8;    // align to 32 byte
        index1[i] = (r250() % 3) * CHAR_PER_128;
        index2[i] = r250() % 129;
        character[i] = r250() % 4;
    }

    printf("finished.\n");

    elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
    totalElapsedTime += elapsedTime;
    initializationTime = elapsedTime;

    // Test SSE Decode
    if (Mode[0] == 'S' || Mode[0] == 's') {

        for (i=0; i<NumberOfIteration; i++) {
            for (j=0; j<NumberOfAccess; j++) {

    // Ordering of index1 and index2 is not important; this module will handle the ordering
    // index1 and index2 can be on the same aligned 128 bit region or can be on adjacant aligned 128 bit region
    // If index1 and index2 are in the same aligned 128 bit region, one of them must be on the boundary
    // These requirements are to reduce the no. of branches in the program flow

    // Pre-fetch next read
    //_mm_prefetch((char*)(memory + address[j+1]), _MM_HINT_NTA);
    //_mm_prefetch((char*)(memory + address[j+2]), _MM_HINT_NTA);
    //_mm_prefetch((char*)(memory + address[j+3]), _MM_HINT_NTA);
    //_mm_prefetch((char*)(memory + address[j+4]), _MM_HINT_NTA);

    // Sort index1 and index2
    temp = (index1[j] - index2[j]) & -(index1[j] < index2[j]);
    minIndex = index2[j] + temp;
    maxIndex = index1[j] - temp;

    // Locate 128 bit boundary
    minIndex128 = lastAlignedBoundary(minIndex, CHAR_PER_128);
    maxIndex128 = lastAlignedBoundary(maxIndex - (maxIndex - minIndex > CHAR_PER_128), CHAR_PER_128);

    // Determine no.of characters to count
    numChar1 = maxIndex128 - minIndex;
    numChar2 = maxIndex - maxIndex128;

    // Set character extraction masks 
    m0 = _mm_set1_epi32(0xFFFFFFFF + (character[j] & 1));    // Character selection mask for even bits
    m1 = _mm_set1_epi32(0xFFFFFFFF + (character[j] >> 1));    // Character selection mask for odd bits
    mcl = _mm_set1_epi32(0x55555555);                        // Set bit-clearing mask to 0x55555555....(alternate 1-bit)

    // Set counting mask for 2 x 128 bits
    
    r1a = _mm_set1_epi32(numChar1);        // Load number of characters into register
    r2a = _mm_set1_epi32(numChar2);        // Load number of characters into register

    r1b = _mm_load_si128((__m128i*)partitionOne1);    // Load partition into register
    r2b = _mm_load_si128((__m128i*)partitionOne2);    // Load partition into register

    r1c = _mm_load_si128((__m128i*)partitionZero1);    // Load partition into register
    r2c = _mm_load_si128((__m128i*)partitionZero2);    // Load partition into register

    r1b = _mm_cmpgt_epi32(r1a, r1b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    r2b = _mm_cmpgt_epi32(r2a, r2b);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
                
    r1c = _mm_cmpgt_epi32(r1a, r1c);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    r2c = _mm_cmpgt_epi32(r2a, r2c);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros

    r1b = _mm_srli_epi32(r1b, (16 - numChar1 % 16) * 2);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary 
    r2b = _mm_slli_epi32(r2b, (16 - numChar2 % 16) * 2);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary

    r1c = _mm_or_si128(r1b, r1c);    // Combine two masks
    r2c = _mm_or_si128(r2b, r2c);    // Combine two masks

    r1c = _mm_and_si128(r1c, mcl);    // Combine with bit-clearing mask (now = 0x55555555....)
    r2c = _mm_and_si128(r2c, mcl);    // Combine with bit-clearing mask (now = 0x55555555....)

    // Load encoding into register
    r1e = _mm_load_si128((__m128i *)(memory + address[j] + minIndex128 / CHAR_PER_WORD));    // Load encoding into register
    r2e = _mm_load_si128((__m128i *)(memory + address[j] + maxIndex128 / CHAR_PER_WORD));    // Load encoding into register

    // Start counting

    r1b = _mm_srli_epi32(r1e, 1);    // Shift encoding to right by 1 bit
    r2b = _mm_srli_epi32(r2e, 1);    // Shift encoding to right by 1 bit

    r1a = _mm_xor_si128(r1e, m0);    // Check even-bits with mask
    r2a = _mm_xor_si128(r2e, m0);    // Check even-bits with mask

    r1b = _mm_xor_si128(r1b, m1);    // Check odd-bits with mask
    r2b = _mm_xor_si128(r2b, m1);    // Check odd-bits with mask

    r1a = _mm_and_si128(r1a, r1b);    // Combine even and odd bits
    r2a = _mm_and_si128(r2a, r2b);    // Combine even and odd bits

    r1a = _mm_and_si128(r1a, r1c);    // Combine with counting mask, which has been combined with bit-clearing mask of 0x55555555.... 
    r2a = _mm_and_si128(r2a, r2c);    // Combine with counting mask, which has been combined with bit-clearing mask of 0x55555555.... 

    // Combine 2 x 128 bits and continue counting

    r1a = _mm_add_epi32(r1a, r2a);        // Combine 2 x 128 bits by adding them together

    mcl = _mm_set1_epi32(0x33333333);    // Set bit-clearing mask to 0x33333333....(alternate 2-bits)

    r1b = _mm_srli_epi32(r1a, 2);        // Shift intermediate result to right by 2 bit
    r1a = _mm_and_si128(r1a, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1b = _mm_and_si128(r1b, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together

    mcl = _mm_set1_epi32(0x0F0F0F0F);    // Set bit-clearing mask to 0x0F0F0F0F....(alternate 4-bits)
    m0 = _mm_setzero_si128();            // Set an all-zero mask

    r1b = _mm_srli_epi32(r1a, 4);        // Shift intermediate result to right by 2 bit
    r1a = _mm_add_epi32(r1a, r1b);        // Combine shifted and non-shifted intermediate results by adding them together
    r1a = _mm_and_si128(r1a, mcl);        // Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)

    r1a = _mm_sad_epu8(r1a, m0);        // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit

    r = _mm_extract_epi16(r1a, 0) + _mm_extract_epi16(r1a, 4);    // Extract result from register and store into variable

/*
                if (index1[j] < index2[j]) {
                    r1 = ForwardDNAOccCount(memory + address[j] + index1[j] / CHAR_PER_WORD, index2[j] - index1[j], character[j], dnaDecodeTable);
                } else {
                    r1 = BackwardDNAOccCount(memory + address[j] + index1[j] / CHAR_PER_WORD, index1[j] - index2[j], character[j], dnaDecodeTable);
                }
                if (r != r1) {
                    printf("%u\n",j);
                }

*/
                t +=r;

            }
        }
    }

    // Test SSE Decode all alphabets
    if (Mode[0] == 'A' || Mode[0] == 'a') {

        for (i=0; i<NumberOfIteration; i++) {
            for (j=0; j<NumberOfAccess; j++) {

    // Ordering of index1 and index2 is not important; this module will handle the ordering
    // index1 and index2 can be on the same aligned 128 bit region or can be on adjacant aligned 128 bit region
    // If index1 and index2 are in the same aligned 128 bit region, one of them must be on the boundary
    // These requirements are to reduce the no. of branches in the program flow

    // Pre-fetch next read
    //_mm_prefetch((char*)(memory + address[j+1]), _MM_HINT_NTA);
    //_mm_prefetch((char*)(memory + address[j+2]), _MM_HINT_NTA);
    //_mm_prefetch((char*)(memory + address[j+3]), _MM_HINT_NTA);
    //_mm_prefetch((char*)(memory + address[j+4]), _MM_HINT_NTA);

    // Sort index1 and index2
    r = (index1[j] - index2[j]) & -(index1[j] < index2[j]);
    minIndex = index2[j] + r;
    maxIndex = index1[j] - r;

    // Locate 128 bit boundary
    minIndex128 = lastAlignedBoundary(minIndex, CHAR_PER_128);
    maxIndex128 = lastAlignedBoundary(maxIndex - (maxIndex - minIndex > CHAR_PER_128), CHAR_PER_128);

    // Determine no.of characters to count
    numChar1 = maxIndex128 - minIndex;
    numChar2 = maxIndex - maxIndex128;

    // Set character extraction masks 
    mcl = _mm_set1_epi32(0x55555555);                        // Set bit-clearing mask to 0x55555555....(alternate 1-bit)

    // Set counting mask for 2 x 128 bits
    ra1 = _mm_set1_epi32(numChar1);        // Load number of characters into register
    ra2 = _mm_set1_epi32(numChar2);        // Load number of characters into register

    rc1 = _mm_load_si128((__m128i*)partitionOne1);    // Load partition into register
    rc2 = _mm_load_si128((__m128i*)partitionOne2);    // Load partition into register

    rg1 = _mm_load_si128((__m128i*)partitionZero1);    // Load partition into register
    rg2 = _mm_load_si128((__m128i*)partitionZero2);    // Load partition into register

    rc1 = _mm_cmpgt_epi32(ra1, rc1);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones
    rc2 = _mm_cmpgt_epi32(ra2, rc2);                // Compare to generate 4x32 bit mask; the word with counting boundary is all ones

    rg1 = _mm_cmpgt_epi32(ra1, rg1);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros
    rg2 = _mm_cmpgt_epi32(ra2, rg2);                // Compare to generate 4x32 bit mask; the word with counting boundary is all zeros

    rc1 = _mm_srli_epi32(rc1, (16 - numChar1 % 16) * 2);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary 
    rc2 = _mm_slli_epi32(rc2, (16 - numChar2 % 16) * 2);    // Shift bits so that all word comform to the requirement of counting the word with counting boundary

    ra1 = _mm_or_si128(rc1, rg1);    // Combine two masks
    ra2 = _mm_or_si128(rc2, rg2);    // Combine two masks

    // Load encoding into register
    r1e = _mm_load_si128((__m128i *)(memory + address[j] + minIndex128 / CHAR_PER_WORD));    // Load encoding into register
    r2e = _mm_load_si128((__m128i *)(memory + address[j] + maxIndex128 / CHAR_PER_WORD));    // Load encoding into register

    // Start counting
    r1e = _mm_and_si128(r1e, ra1);    // Combine encoding with counting mask
    r2e = _mm_and_si128(r2e, ra2);    // Combine encoding with counting mask

    // ra1, ra2, rc1, rc2, rg1, rg2, rt1, rt2 all retired

    // Shift and combine with character selection mask

    ra1 = _mm_srli_epi32(r1e, 1);    // Shift encoding to right by 1 bit
    ra2 = _mm_srli_epi32(r2e, 1);    // Shift encoding to right by 1 bit

    rt1 = _mm_and_si128(r1e, mcl);    // Check even-bits = '1'
    rt2 = _mm_and_si128(r2e, mcl);    // Check even-bits = '1'

    rg1 = _mm_and_si128(ra1, mcl);    // Check odd-bits = '1'
    rg2 = _mm_and_si128(ra2, mcl);    // Check odd-bits = '1'

    rc1 = _mm_andnot_si128(r1e, mcl);    // Check even-bits = '0'
    rc2 = _mm_andnot_si128(r2e, mcl);    // Check even-bits = '0'

    ra1 = _mm_andnot_si128(ra1, mcl);    // Check odd-bits = '0'
    ra2 = _mm_andnot_si128(ra2, mcl);    // Check odd-bits = '0'

    // r1e, r2e retired

    // Count for 'c' 'g' 't'

    r1e = _mm_and_si128(ra1, rt1);        // Combine even and odd bits
    r2e = _mm_and_si128(ra2, rt2);        // Combine even and odd bits
    ra1 = _mm_and_si128(rg1, rc1);        // Combine even and odd bits
    ra2 = _mm_and_si128(rg2, rc2);        // Combine even and odd bits
    rc1 = _mm_and_si128(rg1, rt1);        // Combine even and odd bits
    rc2 = _mm_and_si128(rg2, rt2);        // Combine even and odd bits

    rc = _mm_add_epi32(r1e, r2e);        // Combine 2 x 128 bits by adding them together
    rg = _mm_add_epi32(ra1, ra2);        // Combine 2 x 128 bits by adding them together
    rt = _mm_add_epi32(rc1, rc2);        // Combine 2 x 128 bits by adding them together

    // All except rc, rg, rt retired

    // Continue counting rc, rg, rt

    mcl = _mm_set1_epi32(0x33333333);    // Set bit-clearing mask to 0x33333333....(alternate 2-bits)

    rc1 = _mm_srli_epi32(rc, 2);        // Shift intermediate result to right by 2 bit
    rg1 = _mm_srli_epi32(rg, 2);        // Shift intermediate result to right by 2 bit
    rt1 = _mm_srli_epi32(rt, 2);        // Shift intermediate result to right by 2 bit

    rc2 = _mm_and_si128(rc, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    rg2 = _mm_and_si128(rg, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    rt2 = _mm_and_si128(rt, mcl);        // Clear alternate 2-bits of intermediate result by combining with bit-clearing mask (now = 0x33333333....)

    rc1 = _mm_and_si128(rc1, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    rg1 = _mm_and_si128(rg1, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)
    rt1 = _mm_and_si128(rt1, mcl);        // Clear alternate 2-bits of shifted intermediate result by combining with bit-clearing mask (now = 0x33333333....)

    rc = _mm_add_epi32(rc1, rc2);        // Combine shifted and non-shifted intermediate results by adding them together
    rg = _mm_add_epi32(rg1, rg2);        // Combine shifted and non-shifted intermediate results by adding them together
    rt = _mm_add_epi32(rt1, rt2);        // Combine shifted and non-shifted intermediate results by adding them together

    mcl = _mm_set1_epi32(0x0F0F0F0F);    // Set bit-clearing mask to 0x0F0F0F0F....(alternate 4-bits)
    r1e = _mm_setzero_si128();            // Set an all-zero mask

    rc1 = _mm_srli_epi32(rc, 4);        // Shift intermediate result to right by 2 bit
    rg1 = _mm_srli_epi32(rg, 4);        // Shift intermediate result to right by 2 bit
    rt1 = _mm_srli_epi32(rt, 4);        // Shift intermediate result to right by 2 bit

    rc2 = _mm_add_epi32(rc, rc1);        // Combine shifted and non-shifted intermediate results by adding them together
    rg2 = _mm_add_epi32(rg, rg1);        // Combine shifted and non-shifted intermediate results by adding them together
    rt2 = _mm_add_epi32(rt, rt1);        // Combine shifted and non-shifted intermediate results by adding them together

    rc = _mm_and_si128(rc2, mcl);        // Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)
    rg = _mm_and_si128(rg2, mcl);        // Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)
    rt = _mm_and_si128(rt2, mcl);        // Clear alternate 4-bits of intermediate result by combining with bit-clearing mask (now = 0xOFOFOFOF....)

    rc = _mm_sad_epu8(rc, r1e);            // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    rg = _mm_sad_epu8(rg, r1e);            // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit
    rt = _mm_sad_epu8(rt, r1e);            // Treating the 128 bit as 16 x 8 bit; summing up the 1st 8 x 8 bit into 1st 64-bit and 2nd 8 x 8 bit into 2nd 64-bit

    r4[1] = _mm_extract_epi16(rc, 0) + _mm_extract_epi16(rc, 4);    // Extract result from register and store into variable
    r4[2] = _mm_extract_epi16(rg, 0) + _mm_extract_epi16(rg, 4);    // Extract result from register and store into variable
    r4[3] = _mm_extract_epi16(rt, 0) + _mm_extract_epi16(rt, 4);    // Extract result from register and store into variable

    r4[0] = maxIndex - minIndex - r4[1] - r4[2] - r4[3];

/*
                if (index1[j] < index2[j]) {
                    ForwardDNAAllOccCount(memory + address[j] + index1[j] / CHAR_PER_WORD, index2[j] - index1[j], t4, dnaDecodeTable);
                } else {
                    BackwardDNAAllOccCount(memory + address[j] + index1[j] / CHAR_PER_WORD, index1[j] - index2[j], t4, dnaDecodeTable);
                }
                if (r4[0] != t4[0] || r4[1] != t4[1] || r4[2] != t4[2] || r4[3] != t4[3]) {
                    printf("%u\n",j);
                }

*/

                t4[0] += r4[0];
                t4[1] += r4[1];
                t4[2] += r4[2];
                t4[3] += r4[3];

            }
        }
    }

    // Test decode by lookup
    if (Mode[0] == 'L' || Mode[0] == 'l') {
        for (i=0; i<NumberOfIteration; i++) {
            for (j=0; j<NumberOfAccess; j++) {

                // Pre-fetch next read
                //_mm_prefetch((char*)(memory + address[j+1]), _MM_HINT_T0);
                //_mm_prefetch((char*)(memory + address[j+2]), _MM_HINT_NTA);
                //_mm_prefetch((char*)(memory + address[j+3]), _MM_HINT_NTA);
                //_mm_prefetch((char*)(memory + address[j+4]), _MM_HINT_NTA);

                if (index1[j] < index2[j]) {
                    r = ForwardDNAOccCount(memory + address[j] + index1[j] / CHAR_PER_WORD, index2[j] - index1[j], character[j], dnaDecodeTable);
                } else {
                    r = BackwardDNAOccCount(memory + address[j] + index1[j] / CHAR_PER_WORD, index1[j] - index2[j], character[j], dnaDecodeTable);
                }

                t += r;

            }
        }
    }

    elapsedTime = getElapsedTime(startTime) - totalElapsedTime;
    totalElapsedTime += elapsedTime;
    experimentTime = elapsedTime;

    // So that compiler does not remove code for variables t and r
    if ((Mode[0] == 'S' || Mode[0] == 's' || Mode[0] == 'L' || Mode[0] == 'l') && t == 0) {
        printf("\n");
    }
    if ((Mode[0] == 'A' || Mode[0] == 'a') && t4[0] == 0 && t4[1] == 0 && t4[2] == 0 && t4[3] == 0) {
        printf("\n");
    }

    printf("Experiment completed.\n");

    printf("Initialization time   = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, initializationTime);
    printf("Experiment time       = ");
    printElapsedTime(stdout, FALSE, FALSE, TRUE, 2, experimentTime);
    
    MMUnitFree(address, (NumberOfAccess + 8) * sizeof(unsigned int));
    MMUnitFree(memory, (numberOfMemoryLocation + 32) * sizeof(unsigned int));
    MMUnitFree(index1, (NumberOfAccess + 8) * sizeof(unsigned int));
    MMUnitFree(index2, (NumberOfAccess + 8) * sizeof(unsigned int));
    MMUnitFree(character, (NumberOfAccess + 8) * sizeof(unsigned int));
    MMUnitFree(dnaDecodeTable, 65536 * sizeof(unsigned int));

    iniparser_freedict(programInput);

    return 0;

}

dictionary *ParseInput(int argc, char** argv) {

    dictionary *programInput;

    programInput = paraparser_load(argc, argv, 0, NULL);

    MemorySize = iniparser_getint(programInput, "argument:1", 0);
    if (!MemorySize) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <SSE/Lookup> <No. of iteration>\n", argv[0]);
        exit(1);
    }
    NumberOfAccess = iniparser_getint(programInput, "argument:2", 0);
    if (!NumberOfAccess) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <SSE/Lookup> <No. of iteration>\n", argv[0]);
        exit(1);
    }
    Mode = iniparser_getstring(programInput, "argument:3", 0);
    if (!NumberOfAccess) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <SSE/Lookup> <No. of iteration>\n", argv[0]);
        exit(1);
    }
    NumberOfIteration = iniparser_getint(programInput, "argument:4", 0);
    if (!NumberOfIteration) {
        fprintf(stderr, "Syntax: %s <Memory size> <Number of Access> <SSE/Lookup> <No. of iteration>\n", argv[0]);
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
    if (Mode[0] != 'S' && Mode[0] != 's' && Mode[0] != 'L' && Mode[0] != 'l' && Mode[0] != 'A' && Mode[0] != 'a') {
        fprintf(stderr, "Testing must be SSE, All alphabet SSE or Lookup!\n");
        exit(1);
    }
    if (NumberOfIteration == 0) {
        fprintf(stderr, "Number of Iteration = 0!\n");
        exit(1);
    }

}

