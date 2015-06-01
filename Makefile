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

ARCH := $(shell uname -p)

CC = g++
NVCC = /usr/local/cuda/bin/nvcc
CFLAGS=-O3 -funroll-loops -Wno-unused-result -fopenmp
CUDAFLAG = -cuda --ptxas-options=-v
LIBFLAG = -L/usr/local/cuda/lib64/ -lcuda -lcudart -lpthread -lm -lz
ifeq ($(ARCH), x86_64)
	CUDAFLAG += -arch=sm_35
	CFLAGS += -march=core2 -maccumulate-outgoing-args -mpopcnt -fomit-frame-pointer
else
	CUDAFLAG += -arch=sm_35
	CFLAGS += -D__STDC_LIMIT_MACROS -mcpu=power8 -mtune=power8 -maltivec -fsigned-char

	ifeq ($(USEHUGETLB),yes)
	LIBFLAG += -B/opt/libhugetlbfs/share/libhugetlbfs -Wl,--hugetlbfs-align
	endif
endif

BWTLIB = 2bwt-lib
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/DNACount.o $(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o $(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o $(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o

SAMLIB = samtools-0.1.18
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o $(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o

CPULIB = 2bwt-flex
CPUOBJLIB = $(CPULIB)/HOCC.o $(CPULIB)/LT.o $(CPULIB)/LTConstruct.o $(CPULIB)/HOCCConstruct.o $(CPULIB)/SRA2BWTCheckAndExtend.o $(CPULIB)/SRA2BWTMdl.o

all:	SOAP3-DP SOAP3-Builder BGS-Build BGS-View BGS-View-PE Sample

temp:
	echo $(ARCH)

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

.%.cpp: %.cu Makefile
	$(NVCC) $(CUDAFLAG) $< -o $@

# Second Parse: Normal CC compilation
#

%.o: .%.cpp Makefile
	$(CC) $(CFLAGS) $< -o $@ -c

# Other non-CUDA source files
#
%.o: %.cpp Makefile
	$(CC) $(CFLAGS) $< -o $@ -c

# Application
#
SOAP3-DP:	SOAP3-DP.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o global_arrays.o aio_thread.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB)
	$(CC) $(CFLAGS) SOAP3-DP.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o global_arrays.o aio_thread.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB) $(VPROFLIB) $(LIBFLAG) -o soap3-dp

Sample: sample.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o soap3-dp-module.o global_arrays.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB)
	$(CC) $(CFLAGS) sample.o BGS-IO.o BGS-HostAlgnmtAlgo2.o BGS-HostAlgnmtAlgoSingle.o DV-Kernel.o PE.o SAList.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o DV-DPForBothUnalign.o DV-DPForSingleReads.o soap3-dp-module.o global_arrays.o ssse3_popcount.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB) $(VPROFLIB) $(LIBFLAG) -o sample

BGS-View: BGS-View.cpp $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) BGS-View.cpp $(BWTOBJLIBS) -o BGS-View

BGS-View-PE: BGS-View-PE.cpp $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) BGS-View-PE.cpp $(BWTOBJLIBS) -o BGS-View-PE

SOAP3-Builder: $(CPULIB)/2BWT-Builder.c $(BWTLIB)/BWTConstruct.o $(BWTOBJLIBS) $(CPUOBJLIB) Makefile
	$(CC) $(CFLAGS) $(CPULIB)/2BWT-Builder.c $(BWTLIB)/BWTConstruct.o $(BWTOBJLIBS) $(CPUOBJLIB) -o soap3-dp-builder

BGS-Build: BGS-Build.cpp $(BWTOBJLIBS) Makefile
	$(CC) $(CFLAGS) BGS-Build.cpp $(BWTOBJLIBS) -o BGS-Build

cleanALL:
	echo "Cleaning up all object codes, including 2BWT-LIB, SAM and all other library compiled."
	rm -f .*.cpp *.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB) $(BWTLIB)/BWTConstruct.o

clean:
	echo "Cleaning up all SOAP3 object code includes CPU and CUDA compiled."
	rm -f .*.cpp *.o

cleanCPU:
	echo "Cleaning up only CPU SOAP3 object code."
	rm -f *.o
