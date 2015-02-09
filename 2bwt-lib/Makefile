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

all: 2BWT-Builder 2BWT-Interface.o 2BWT-Sample

BWT.o:	BWT.c BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h Makefile

BWTConstruct.o:	BWTConstruct.c BWTConstruct.h TypeNLimit.h BWT.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h QSufSort.h r250.h Makefile

BlastCompare.o:	BlastCompare.c MiscUtilities.h TypeNLimit.h iniparser.h dictionary.h Makefile

dictionary.o:	dictionary.c dictionary.h inistrlib.h dictionary.h Makefile

blast_dust.o:	blast_dust.c blast_dust.h Makefile

DNACount.o:	DNACount.c DNACount.h TypeNLimit.h MiscUtilities.h Makefile

HSP.o:	HSP.c TextConverter.h TypeNLimit.h MemManager.h MiscUtilities.h r250.h HSP.h HSPstatistic.h Makefile

HSPstatistic.o:	HSPstatistic.c karlin.h HSPstatistic.h Makefile

iniparser.o:	iniparser.c iniparser.h dictionary.h inistrlib.h Makefile

inistrlib.o:	inistrlib.c inistrlib.h Makefile

karlin.o:	karlin.c karlin.h Makefile

MemManager.o:	MemManager.c MiscUtilities.h TypeNLimit.h MemManager.h Makefile

MiscUtilities.o:	MiscUtilities.c MiscUtilities.h TypeNLimit.h Makefile

QSufSort.o:	QSufSort.c QSufSort.h TypeNLimit.h MiscUtilities.h Makefile

r250.o:	r250.c r250.h Makefile

Socket.o:	Socket.c Socket.h TypeNLimit.h MemManager.h MiscUtilities.h Makefile

TextConverter.o:	TextConverter.c TextConverter.h TypeNLimit.h MemManager.h MiscUtilities.h r250.h Makefile

Timing.o:	Timing.c Timing.h Makefile

IndexCtl.o:	IndexCtl.c IndexCtl.h Makefile

2BWT-Builder.o:	2BWT-Builder.c Makefile

2BWT-Builder:	2BWT-Builder.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o Timing.o Socket.o HSP.o karlin.o HSPstatistic.o IndexCtl.o Makefile
	$(CC) $(CFLAGS) 2BWT-Builder.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o HSP.o karlin.o HSPstatistic.o Timing.o Socket.o IndexCtl.o -o 2bwt-builder -lm

2BWT-Sample:	2BWT-Sample.o 2BWT-Interface.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o Timing.o Socket.o HSP.o karlin.o HSPstatistic.o IndexCtl.o Makefile
	$(CC) $(CFLAGS) 2BWT-Sample.o 2BWT-Interface.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o HSP.o karlin.o HSPstatistic.o Timing.o Socket.o IndexCtl.o -o 2bwt-sample -lm

2BWT-Interface.o:	2BWT-Interface.c 2BWT-Interface.h Makefile

