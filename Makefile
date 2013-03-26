#
#
#    Makefile
#    Soap3(gpu)
#
#    Copyright (C) 2011, HKU
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 2
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#
####################################################################
#
#   Modification History
#
#   Date   : 6th January 2012
#   Author : Edward MK Wu
#   Change : Standardising SOAP3 package
#
####################################################################

CC = g++
NVCC = /usr/local/cuda/bin/nvcc
CUDAFLAG = -cuda -arch=sm_20 --ptxas-options=-v
LIBFLAG = -L/usr/local/cuda/lib64/ -lcuda -lcudart
CFLAGS = -O3 -funroll-loops -march=core2 -fomit-frame-pointer -maccumulate-outgoing-args -Wno-unused-result -lm -static-libgcc -mpopcnt

BWTLIB = 2bwt-lib
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/DNACount.o $(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o $(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o $(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o
#BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/DNACount.o $(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o

SAMLIB = samtools-0.1.18
#SAMOBJLIBS = $(SAMLIB)/sam.c $(SAMLIB)/bam.c $(SAMLIB)/bgzf.c $(SAMLIB)/kstring.c $(SAMLIB)/bam_import.c $(SAMLIB)/faidx.c $(SAMLIB)/bam_pileup.c $(SAMLIB)/bam_aux.c $(SAMLIB)/sam_header.c $(SAMLIB)/razf.c
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o $(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o

CPULIB = 2bwt-flex
CPUOBJLIB = $(CPULIB)/HOCC.o $(CPULIB)/LT.o $(CPULIB)/LTConstruct.o $(CPULIB)/HOCCConstruct.o $(CPULIB)/SRA2BWTCheckAndExtend.o $(CPULIB)/SRA2BWTMdl.o

all:	SOAP3-DP SOAP3-Builder BGS-Build BGS-View BGS-View-PE Sample 

other:	BGS-View BGS-View-PE Sample

$(BWTLIB) :	force_look
	cd $(BWTLIB); $(MAKE)

$(SAMLIB) :	force_look
	cd $(SAMLIB); $(MAKE)

$(CPULIB) : 	force_look
	cd $(CPULIB); $(MAKE) 


# Why do we parse the source code twice? The reason we perform the *.cu->*.cpp->*o conversion
# is that NVCC does not seem to support -funroll-loops compilation flag. Converting *.cu->*.o
# directly will result in atl least 10% performance reducation.
#
# First Parse: Compiling CUDA into standard CPP
#
.sample.cpp :	sample.cu Makefile
			$(NVCC) $(CUDAFLAG) sample.cu -o .sample.cpp

.soap3-dp-module.cpp :	soap3-dp-module.cu Makefile
			$(NVCC) $(CUDAFLAG) soap3-dp-module.cu -o .soap3-dp-module.cpp

.SOAP3-DP.cpp :	SOAP3-DP.cu Makefile
			$(NVCC) $(CUDAFLAG) SOAP3-DP.cu -o .SOAP3-DP.cpp

.DV-Kernel.cpp :	DV-Kernel.cu Makefile
			$(NVCC) $(CUDAFLAG) DV-Kernel.cu -o .DV-Kernel.cpp

.DV-DPfunctions.cpp :	DV-DPfunctions.cu Makefile
			$(NVCC) $(CUDAFLAG) DV-DPfunctions.cu -o .DV-DPfunctions.cpp

.alignment.cpp :		alignment.cu Makefile
			$(NVCC) $(CUDAFLAG) alignment.cu -o .alignment.cpp

.DV-SemiDP.cpp :		DV-SemiDP.cu Makefile
			$(NVCC) $(CUDAFLAG) DV-SemiDP.cu -o .DV-SemiDP.cpp

.DV-DPForBothUnalign.cpp:	DV-DPForBothUnalign.cu Makefile
			$(NVCC) $(CUDAFLAG) DV-DPForBothUnalign.cu -o .DV-DPForBothUnalign.cpp

.DV-DPForSingleReads.cpp:	DV-DPForSingleReads.cu Makefile
			$(NVCC) $(CUDAFLAG) DV-DPForSingleReads.cu -o .DV-DPForSingleReads.cpp

# Second Parse: Normal CC compilation
#
sample.o:	.sample.cpp Makefile
	$(CC) $(CFLAGS) -c .sample.cpp -o sample.o

soap3-dp-module.o:	.soap3-dp-module.cpp Makefile
	$(CC) $(CFLAGS) -c .soap3-dp-module.cpp -o soap3-dp-module.o

SOAP3-DP.o:	.SOAP3-DP.cpp Makefile
	$(CC) $(CFLAGS) -c .SOAP3-DP.cpp -o SOAP3-DP.o

DV-Kernel.o:	.DV-Kernel.cpp Makefile
	$(CC) $(CFLAGS) -c .DV-Kernel.cpp -o DV-Kernel.o

DV-DPfunctions.o:	.DV-DPfunctions.cpp Makefile
	$(CC) $(CFLAGS) -c .DV-DPfunctions.cpp -o DV-DPfunctions.o

alignment.o:	.alignment.cpp Makefile
	$(CC) $(CFLAGS) -c .alignment.cpp -o alignment.o

DV-SemiDP.o:	.DV-SemiDP.cpp Makefile
	$(CC) $(CFLAGS) -c .DV-SemiDP.cpp -o DV-SemiDP.o

DV-DPForBothUnalign.o:	.DV-DPForBothUnalign.cpp Makefile
	$(CC) $(CFLAGS) -c .DV-DPForBothUnalign.cpp -o DV-DPForBothUnalign.o

DV-DPForSingleReads.o:	.DV-DPForSingleReads.cpp Makefile
	$(CC) $(CFLAGS) -c .DV-DPForSingleReads.cpp -o DV-DPForSingleReads.o

# Other non-CUDA source files
#
BGS-IO.o:	BGS-IO.c Makefile
	$(CC) $(CFLAGS) -c BGS-IO.c -o BGS-IO.o

BGS-HostAlgnmtMdl.o:	BGS-HostAlgnmtMdl.c Makefile
	$(CC) $(CFLAGS) -c BGS-HostAlgnmtMdl.c -o BGS-HostAlgnmtMdl.o

BGS-HostAlgnmtAlgo2.o:	BGS-HostAlgnmtAlgo2.c Makefile
	$(CC) $(CFLAGS) -c BGS-HostAlgnmtAlgo2.c -o BGS-HostAlgnmtAlgo2.o

BGS-HostAlgnmtAlgoSingle.o:	BGS-HostAlgnmtAlgoSingle.c Makefile
	$(CC) $(CFLAGS) -c BGS-HostAlgnmtAlgoSingle.c -o BGS-HostAlgnmtAlgoSingle.o

PE.o:	PE.c Makefile
	$(CC) $(CFLAGS) -c PE.c -o PE.o

SAList.o:	SAList.c Makefile
	$(CC) $(CFLAGS) -c SAList.c -o SAList.o

UsageInterface.o:	UsageInterface.c Makefile
	$(CC) $(CFLAGS) -c UsageInterface.c -o UsageInterface.o

IndexHandler.o:	IndexHandler.c Makefile
	$(CC) $(CFLAGS) -c IndexHandler.c -o IndexHandler.o

IniParam.o:	IniParam.c Makefile
	$(CC) $(CFLAGS) -c IniParam.c -o IniParam.o

FileUtilities.o:	FileUtilities.c Makefile
	$(CC) $(CFLAGS) -c FileUtilities.c -o FileUtilities.o
	
CPUfunctions.o:	CPUfunctions.cpp Makefile
	$(CC) $(CFLAGS) -c CPUfunctions.cpp -o CPUfunctions.o

PEAlgnmt.o:	PEAlgnmt.c Makefile
	$(CC) $(CFLAGS) -c PEAlgnmt.c -o PEAlgnmt.o

SAM.o:	SAM.c Makefile
	$(CC) $(CFLAGS) -c SAM.c -o SAM.o

OutputDPResult.o:	OutputDPResult.cpp Makefile
	$(CC) $(CFLAGS) -c OutputDPResult.cpp -o OutputDPResult.o

AlgnResult.o:	AlgnResult.cpp Makefile
	$(CC) $(CFLAGS) -c AlgnResult.cpp -o AlgnResult.o

global_arrays.o:	global_arrays.cpp Makefile
	$(CC) $(CFLAGS) -c global_arrays.cpp -o global_arrays.o

aio_thread.o:	aio_thread.cpp Makefile
	$(CC) $(CFLAGS) -c aio_thread.cpp -o aio_thread.o

ssse3_popcount.o:	ssse3_popcount.c Makefile
	$(CC) $(CFLAGS) -c ssse3_popcount.c -o ssse3_popcount.o

QueryParser.o:	QueryParser.c Makefile
	$(CC) $(CFLAGS) -c QueryParser.c -o QueryParser.o
	
# Application
#
SOAP3-DP:	SOAP3-DP.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o global_arrays.o aio_thread.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB)
	$(CC) $(CFLAGS) $(LIBFLAG) SOAP3-DP.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o global_arrays.o aio_thread.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB) -o soap3-dp

Sample: sample.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o soap3-dp-module.o global_arrays.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB)
	$(CC) $(CFLAGS) $(LIBFLAG) sample.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o soap3-dp-module.o global_arrays.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB) -o sample

BGS-View: BGS-View.c $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) BGS-View.c $(BWTOBJLIBS) -o BGS-View

BGS-View-PE: BGS-View-PE.c $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) BGS-View-PE.c $(BWTOBJLIBS) -o BGS-View-PE

SOAP3-Builder: $(CPULIB)/2BWT-Builder.c $(BWTLIB)/BWTConstruct.o $(BWTOBJLIBS) $(CPUOBJLIB) Makefile
	$(CC) $(CFLAGS) $(CPULIB)/2BWT-Builder.c $(BWTLIB)/BWTConstruct.o $(BWTOBJLIBS) $(CPUOBJLIB) -o soap3-dp-builder

BGS-Build: BGS-Build.cpp $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) BGS-Build.cpp $(BWTOBJLIBS) -o BGS-Build

cleanALL:
	echo "Cleaning up all object codes, including 2BWT-LIB, SAM and all other library compiled."
	rm -f .*.cpp *.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB)

clean:
	echo "Cleaning up all SOAP3 object code includes CPU and CUDA compiled."
	rm -f .*.cpp *.o

cleanCPU:
	echo "Cleaning up only CPU SOAP3 object code."
	rm -f *.o
