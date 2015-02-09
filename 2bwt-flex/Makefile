#
#    Makefile
#
#    SOAP2 / 2BWT
#
#    Copyright (C) 2012, HKU
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
CFLAGS = -O3 -funroll-loops -march=nocona -maccumulate-outgoing-args -Wno-unused-result -Wno-write-strings
LIBS = -lm -lz
THREADLIBS = -lpthread

UTILLIB = utilities
UTILOBJLIBS = $(UTILLIB)/BufferedFileReader.o

BWTLIB = 2bwt-lib
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/blast_dust.o $(BWTLIB)/DNACount.o \
			$(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o \
			$(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o \
			$(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o

PELIB = pair-end
PEOBJLIBS = $(PELIB)/PEAlgnmt.o

SAMLIB = samtools-0.1.18
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o \
				$(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o

SRALIB = .
SRAOBJLIBS = $(SRALIB)/LT.o $(SRALIB)/HOCC.o $(SRALIB)/OutputAnyFormat.o $(SRALIB)/OutputBinaryFormat.o \
			$(SRALIB)/OutputSAMFormat.o $(SRALIB)/SRACore.o $(SRALIB)/SRAOccCollector.o $(SRALIB)/SRA2BWTReport.o \
			$(SRALIB)/SRA2BWTMdl.o $(SRALIB)/SRAConfigParser.o $(SRALIB)/SRA2BWTCheckAndExtend.o $(SRALIB)/SRA2BWTOperations.o \
			$(SRALIB)/SRAQueryParser.o $(SRALIB)/2BWT-SRAAlgnmt.o $(SRALIB)/DPCore.o $(SRALIB)/DPOperations.o $(SRALIB)/DPOccCollector.o \
			$(SRALIB)/SRADPReport.o
PEAOBJLIBS = $(SRALIB)/2BWT-PEAlgnmt.o $(SRALIB)/PE2BWTReport.o $(SRALIB)/PEDPReport.o

all:	SOAP2-Builder SOAP2-CPU SOAP2-CPU-SE SOAP2-CPU-PE SOAP2-Viewer

clean:
	rm $(UTILLIB)/*.o $(BWTLIB)/*.o $(PELIB)/*.o $(SRALIB)/*.o *.o

$(UTILLIB) : force_look
	cd $(UTILLIB); $(MAKE)
$(SAMLIB) : force_look
	cd $(SAMLIB); $(MAKE)
$(BWTLIB) : force_look
	cd $(BWTLIB); $(MAKE)
$(PELIB) : force_look
	cd $(PELIB); $(MAKE)

# SOAP2-CPU
2BWT-Aligner.o:				2BWT-Aligner.c Makefile
SOAP2-CPU: 2BWT-Aligner.o Makefile
	$(CC) $(CFLAGS) $(LIBS) 2BWT-Aligner.o -o soap2-dp

# SOAP2-CPU-SE
LT.o:	LT.c LT.h Makefile
HOCC.o:	HOCC.c HOCC.h Makefile
OutputAnyFormat.o:          OutputAnyFormat.c OutputAnyFormat.h Makefile
OutputBinaryFormat.o:       OutputBinaryFormat.c OutputBinaryFormat.h Makefile
OutputSAMFormat.o:          OutputSAMFormat.c OutputSAMFormat.h Makefile
SRACore.o:	                SRACore.c SRACore.h Makefile
SRAOccCollector.o:	        SRAOccCollector.c SRAOccCollector.h Makefile
SRA2BWTReport.o:	        SRA2BWTReport.c SRA2BWTReport.h Makefile
SRA2BWTMdl.o:	            SRA2BWTMdl.c SRA2BWTMdl.h Makefile
SRAConfigParser.o:		    SRAConfigParser.c SRAConfigParser.h Makefile
SRA2BWTCheckAndExtend.o:    SRA2BWTCheckAndExtend.c SRA2BWTCheckAndExtend.h Makefile
SRA2BWTOperations.o:        SRA2BWTOperations.c SRA2BWTOperations.h Makefile
SRAQueryParser.o:		    SRAQueryParser.c SRAQueryParser.h Makefile
DPCore.o:					DPCore.c DPCore.h Makefile
DPOperations.o:				DPOperations.c DPOperations.h Makefile
DPOccCollector.o:	        DPOccCollector.c DPOccCollector.h Makefile
2BWT-SRAAlgnmt.o:		    2BWT-SRAAlgnmt.c 2BWT-SRAAlgnmt.h Makefile
SRADPReport.o:		    	SRADPReport.c SRADPReport.h Makefile
2BWT-Combined.o:	        2BWT-Combined.c Makefile
SOAP2-CPU-SE: 2BWT-Combined.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) 2BWT-Combined.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o soap2-cpu-se
# SOAP2-Viewer
2BWT-Viewer.o:	2BWT-Viewer.c 2BWT-Viewer.h Makefile
SOAP2-Viewer:	2BWT-Viewer.o $(BWTOBJLIBS) $(UTILOBJLIBS) SRAQueryParser.o OutputAnyFormat.o OutputBinaryFormat.o Makefile
	$(CC) $(CFLAGS) $(LIBS) 2BWT-Viewer.o $(BWTOBJLIBS) $(UTILOBJLIBS) SRAQueryParser.o OutputAnyFormat.o OutputBinaryFormat.o -o soap2-dp-viewer

# SOAP2-Builder
LTConstruct.o:	LTConstruct.c LTConstruct.h Makefile
HOCCConstruct.o:	HOCCConstruct.c HOCCConstruct.h Makefile
2BWT-Builder.o:	2BWT-Builder.c Makefile
SOAP2-Builder:	2BWT-Builder.o $(BWTOBJLIBS) $(BWTLIB)/BWTConstruct.o LTConstruct.o HOCCConstruct.o Makefile
	$(CC) $(CFLAGS) $(LIBS) 2BWT-Builder.o $(BWTOBJLIBS) $(BWTLIB)/BWTConstruct.o LTConstruct.o HOCCConstruct.o -o soap2-dp-builder

# SOAP2-CPU-PE
PE2BWTReport.o:	            PE2BWTReport.c PE2BWTReport.h Makefile
2BWT-PEAlgnmt.o:	        2BWT-PEAlgnmt.c 2BWT-PEAlgnmt.h Makefile
2BWT-PE.o:	                2BWT-PE.c Makefile
PEDPReport.o:	            PEDPReport.c Makefile
SOAP2-CPU-PE: 2BWT-PE.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) 2BWT-PE.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS) -o soap2-cpu-pe
