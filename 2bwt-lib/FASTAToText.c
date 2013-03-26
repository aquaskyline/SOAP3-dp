/* FASTAtoText.c        FASTA to Text

   Copyright 2004, Wong Chi Kwong, all rights reserved.

   This module convert FASTA DNA format to plain text format.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/

#include <stdio.h>
#include "TypeNLimit.h"
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {

    int c;
    FILE* inputFile;
    FILE* outputFile;
    char uppercase[256];
    char lowercase[256];
    char dna[4] = {'A', 'C', 'G', 'T'};

    if (argc != 3 && argc != 4) {
        fprintf(stderr, "%s <input> <output>\n", argv[0]);
        fprintf(stderr, "%s <input> <output> a\n", argv[0]);
        fprintf(stderr, "%s <input> <output> n\n", argv[0]);
        fprintf(stderr, "%s <input> <output> r\n", argv[0]);
        return 1;
    }
    inputFile = (FILE*)fopen64(argv[1], "r");
    if (inputFile == NULL) {
        fprintf(stderr, "Input file not exists!\n");
        return 1;
    }
    outputFile = (FILE*)fopen64(argv[2], "r");
    if (outputFile != NULL) {
        fprintf(stderr, "Output file exists!\n");
        return 1;
    }
    outputFile = (FILE*)fopen64(argv[2], "wb");

    srand(clock());

    for (c=0; c<256; c++) {
        uppercase[c] = (char)c;
    }
    for (c='a'; c<='z'; c++) {
        uppercase[c] = (char)c + 'A' - 'a';
    }
    for (c=0; c<256; c++) {
        lowercase[c] = (char)c;
    }
    for (c='A'; c<='Z'; c++) {
        lowercase[c] = (char)c + 'a' - 'A';
    }

    c = fgetc(inputFile);
    while (c == '\n' || c == '\r') {
        c = fgetc(inputFile);
    }
    while (c == '>') {
        while (c != '\n' && c != '\r' && c != EOF) {
            c = fgetc(inputFile);
        }
        while (c == '\n' || c == '\r') {
            c = fgetc(inputFile);
        }
    }

    while (c != EOF) {
        if (argc == 4 || (c != 'n' && c != 'N')) {
            if ((c == 'n' || c == 'N') && (argc == 4 && argv[3][0] == 'a')) {
                if (c == 'n') {
                    fputc('a', outputFile);
                } else {
                    fputc('A', outputFile);
                }
            } else {
                if ((c == 'n' || c == 'N') && (argc == 4 && argv[3][0] == 'r')) {
                    if (c == 'n') {
                        fputc(lowercase[dna[rand() % 4]], outputFile);
                    } else {
                        fputc(uppercase[dna[rand() % 4]], outputFile);
                    }
                } else {
                    fputc(c, outputFile);
                }
            }
        }
        c = fgetc(inputFile);
        while (c == '\n' || c == '\r') {
            c = fgetc(inputFile);
        }
        while (c == '>') {
            while (c != '\n' && c != '\r' && c != EOF) {
                c = fgetc(inputFile);
            }
            while (c == '\n' || c == '\r') {
                c = fgetc(inputFile);
            }
        }
    }

    fclose(inputFile);
    fclose(outputFile);

}
