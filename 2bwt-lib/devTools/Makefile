####################################################################
#
#   Modification History
#
#   Date   : 24th May 2011
#   Author : Edward MK Wu
#   Change : Rename 2BWT to SOAP
#
#   Date   : 19th June 2011
#   Author : Edward MK Wu
#   Change : Packaging 2BWT library as a separate product.
#            Thus, changing all references to 2bwt lib to subdirectory.
#
####################################################################

CC = g++
DEFINE = 
CFLAGS = -O3 -funroll-loops -march=nocona -maccumulate-outgoing-args -Wno-write-strings -Wno-unused-result 
LIBS = -lm 
THREADLIBS = -lpthread

BWTLIB = ..
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/blast_dust.o $(BWTLIB)/DNACount.o \
             $(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o \
             $(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o \
             $(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o \
             $(BWTLIB)/2BWT-Interface.o

$(BWTLIB) : force_look
	cd $(BWTLIB); $(MAKE)


all: 2BWT-Benchmark

2BWT-Benchmark:	2BWT-Benchmark.o $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) 2BWT-Benchmark.o $(BWTOBJLIBS) -o 2bwt-benchmark

