gcc -c PEAlgnmt.c
export BWTLIB="../2bwt-lib/"

gcc testDrive.c PEAlgnmt.o $BWTLIB/BWT.o $BWTLIB/dictionary.o $BWTLIB/blast_dust.o $BWTLIB/DNACount.o $BWTLIB/HSP.o $BWTLIB/HSPstatistic.o $BWTLIB/iniparser.o $BWTLIB/inistrlib.o $BWTLIB/karlin.o $BWTLIB/MemManager.o $BWTLIB/MiscUtilities.o $BWTLIB/QSufSort.o $BWTLIB/r250.o $BWTLIB/TextConverter.o $BWTLIB/Timing.o $BWTLIB/Socket.o $BWTLIB/2BWT-Interface.o -lm

