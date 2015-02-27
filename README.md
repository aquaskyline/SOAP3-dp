# SOAP3-dp

SOAP3-dp, like its predecessor SOAP3, is a GPU-based software for aligning short reads to a reference sequence. It improves SOAP3 in terms of both speed and sensitivity by skillful exploitation of whole-genome indexing and dynamic programming on a GPU. SOAP3 is limited to find alignments with at most 4 mismatches, while SOAP3-dp can find alignments involving mismatches, INDELs, and small gaps. The number of reads aligned, especially for paired-end data, typically increases 5 to 10 percent from SOAP3 to SOAP3-dp. More interestingly, SOAP3-dp's alignment time is much shorter than SOAP3, as it is found that GPU-based dynamic programming when coupled with indexing can be much more efficient. For example, when aligning length-100 single-end reads with the human genome, SOAP3 typically requires tens of seconds per million reads, while SOAP3-dp takes only a few seconds.

  * [Hardware and Platform](#Section1)
  * [Usage](#Section2)
	  * [Index builder](#Section2.1)
	  * [Aligner](#Section2.2)
	  * [About the output](#Section2.3)
  * [Running SOAP3-dp on a machine with multiple GPU devices](#Section3)
  * [Fine tuning via configuration file soap3-dp.ini](#Section4)

# <a name="Section1"></a>1. Hardware and Platform

To run SOAP3-dp, you need a linux workstation equipped with

  * a multi-core CPU with at least 16GB main memory (default quad-core);
  * a CUDA-enabled GPU with compute capability 2.0 and at least 3GB memory.

SOAP3-dp has been tested with the following GPU: NVIDIA Tesla C2070 (6GB memory), GTX 580 (3.2GB memory) and GTX 680 (4GB memory). And it should also work using Tesla M2050, C2050, C2075 and Quadro 6000.

SOAP3-dp was developed under the 64-bit linux platform and the CUDA Driver version 4.2. SOAP3-dp can also run with CUDA Driver version 4.1 and 4.2 (users are advised to remove the Virtual Memory limit using the command "ulimit -v unlimited").

\* See Section 4 for how to change the default with the INI file soap3-dp.ini.


# <a name="Section2"></a>2. Usage

SOAP3-dp consists of 3 parts: (i) Index (2BWT and GPU-2BWT) builder and (ii) Aligner. There is an auxiliary program to merge the output files.

After downloading the compressed file, run the following commands to extract the different programs:

```bash
  % gunzip soap3-dp-2.1-cudaX-linux.tar.gz
  % tar -xvf soap3-dp-2.1-cudaX-linux.tar
```

You should get file files:
  soap3-dp-builder, soap3-dp-builder.ini, BGS-Build, soap3-dp, soap3-dp.ini.

## <a name="Section2.1"></a>2.1 Index builder

Soap3-builder preprocesses the FASTA reference sequence(s) given in a file [ref sequence file] to produce a few indexes for alignment. The indexes will be stored in files whose names are in the form [ref sequence file].index. Note the following restrictions to the reference sequences:

  * At most 65,000 reference sequences (chromosomes) in one single FASTA file.
  * The total length of all reference sequences is at most 4 billion.
  * All invalid characters (i.e., anything other than A, C, G, T) will be replaced by character G; regions that contain more than 10 invalid characters will not be involved in any alignment.

### Step 1: Build the 2BWT index

Define a sampling rate for SA construction (default: full SA). Using lower sampling rate consumes less host memory but will slightly decrease the speed. In production, if memory allowed, please use full SA. For memory limited environments, please use 1/2 or 1/4 sampling.

Please change option "SaValueFreq" in configuration file "soap3-dp-builder.ini" for following sampling rates:

Full SA, comsumes ~17G host memory:

  > SaValueFreq=1

1/2 SA, comsumes ~11.5G host memory:

  > SaValueFreq=2

1/4 SA, comsumes ~9G host memory:

  > SaValueFreq=4

Then, run the builder:

Syntax:

```bash
  % ./soap3-dp-builder
```

For example:

```bash
  % ./soap3-dp-builder genome.fa
```

### Step 2: Convert the 2BWT index to the GPU-2BWT index

Syntax:

```bash
  % ./BGS-Build .index
```

For example:

```bash
  % ./BGS-Build genome.fa.index
```

Additional files with the filename-prefix "[ref sequence file].index" will be generated.

With all these index files, you are now ready to use soap3-dp to perform alignment for the reads.

## <a name="Section2.2"></a>2.2 Aligner

The aligner supports read files in FASTA format, FASTQ format, compressed gzip format or BAM format. Moreover, the aligner can work on one or multiple input files. Precisely, there are the following few cases:

  1. [One set of paired-end reads in FASTA, FASTQ or GZIP format](#pair)
  2. [One set of paired-end reads in BAM format](#pair-bam)
  3. [Multiple sets of paired-end reads in FASTA, FASTQ or GZIP format](#list)
  4. [Multiple sets of paired-end reads in BAM format](#list-bam)
  5. [One set of single-end reads](#single)
  6. [Multiple sets of single-end reads](#list-single)

### <a name="pair"></a>1. One set of paired-end reads in FASTA, FASTQ or GZIP format

Input reads from two files [read file 1] and [read file 2], and both read files are in FASTA, FASTQ or GZIP format (The aligner can detect these format automatically).

Syntax:

```bash
  % ./soap3-dp pair [ref seq file].index [read file 1] [read file 2] [options]
```

```
Options:
-u [max value of insert size]	                    default: 500;
 	 
-v [min value of insert size]	                    default: 1;
 	 
-L [length of the longest read in the input]	    default: 120; normally in the range [75,200];
 	                                                if < 100, better using GPU with 6GB memory
 	 
-h [output option]	                                1: all valid alignments;
 	                                                2: all best alignment (default);
 	                                                3: unique best alignment;
 	                                                4: random best alignment
 	 
-b [output format]	                                1: Succinct;
 	                                                2: SAM (default);
 	                                                3: BAM
 	 
-o [output prefix]	                                The output file names:
 	                                                [output prefix].dpout.1, [output prefix].gout.1, [output prefix].gout.2, ...
 	                                                default: [read file 1];
 	 
-c [GPU device ID]	                                Specify the GPU device for running SOAP3-dp;
 	 
-I	                                                The input is in the Illumina 1.3+ FASTQ-like format
 	                                                default: this option is not enabled;
 	 
-A [sample name]	                                default: "default";
 	 
-D [read group ID]	                                default: [read file 1];
 	 
-R [read group option]	                            Assign the read group option
 	                                                default: this option is not enabled;

-p                                                  Output MD string and NM tag.
 	 
-s [max # of mismatches]	                        any integer from 0 to 4,
 	                                                disable dynamic programming and
 	                                                perform alignment with mismatches only.
 	                                                If max # of mismatches is not specified,
 	                                                default number is 3 for reads of length >= 50;
 	                                                or 2 otherwise.
 	                                                This option is not recommended.
```

Example: Assume read files query1.fa and query2.fa, and insert size range [200, 500].

```bash
  % ./soap3-dp pair genome.fa.index query1.fa query2.fa -u 500 -v 200
```

See [Section 2.3](#Section2.3) about the details of output format.

### <a name="pair-bam"></a>2. One set of paired end reads

Input reads from a BAM file [read file 1].

Syntax:

```bash
  % ./soap3-dp pair [ref seq file].index [read file 1] -bam [options]
```

[read file 1] is in BAM format and contains both the first and the second reads of the pairs. The first and the second reads of the same pair appear consecutively.

```
Options:
-u [max value of insert size]						default: 500;
 	 
-v [min value of insert size]						default: 1;
 	 
-L [length of the longest read in the input]		default: 120; normally in the range [75,200];
 													if < 100, better using GPU with 6GB memory
 	 
-h [output option]									1: all valid alignments;
 													2: all best alignment (default);
 													3: unique best alignment;
 													4: random best alignment
 	 
-b [output format]									1: Succinct;
 													2: SAM (default);
 													3: BAM
 	 
-o [output prefix]									The output file names:
 													[output prefix].dpout.1, [output prefix].gout.1, [output prefix].gout.2, ...
 													default: [read file 1];
 	 
-c [GPU device ID]									Specify the GPU device for running SOAP3-dp;
 	 
-I													The input is in the Illumina 1.3+ FASTQ-like format
 													default: this option is not enabled;
 	 
-A [sample name]									default: "default";
 	 
-D [read group ID]									default: [read file 1];
 	 
-R [read group option]								Assign the read group option
 													default: this option is not enabled;
 	 
-p                                                  Output MD string and NM tag.

-s [max # of mismatches]							any integer from 0 to 4,
 													disable dynamic programming and
 													perform alignment with mismatches only.
                                                    If max # of mismatches is not specified,
                                                 	default number is 3 for reads of length >= 50;
                                                 	or 2 otherwise.
                                                 	This option is not recommended.
```

### <a name="list"></a>3. Multiple sets of paired-end reads in FASTA, FASTQ or GZIP format.

Syntax:

```bash
  % ./soap3-dp pair-multi [ref seq file].index [info list file] [options]
```

[info list file] contains a list of read information.

```
Format of the file [info list file]:
Each line inside the file contains information of each set of paired-end reads.
 	 
Column 1: first read file name;						Column 2: second read file name;
Column 3: minimum value of insert size;				Column 4: maximum value of insert size;
Column 5: output prefix;							Column 6 (optional): read group ID;
Column 7 (optional): sample name;					Column 8 (optional): read group options.
 	 
Delimiter is a tab character and only the following three cases are allowed:
Case 1: Only contains the first 5 columns.
Case 2: Only contains the first 7 columns.
Case 3: Contains all the 8 columns.
```

```
Options:
-L [length of the longest read in the input]		default: 120; normally in the range [75,200];
 													if < 100, better using GPU with 6GB memory
 	 
-h [output option]									1: all valid alignments;
 													2: all best alignment (default);
 													3: unique best alignment;
 													4: random best alignment
 	 
-b [output format]									1: Succinct;
 													2: SAM (default);
 													3: BAM
 	 
-c [GPU device ID]									Specify the GPU device for running SOAP3-dp;
 	 
-I													The input is in the Illumina 1.3+ FASTQ-like format
 													default: this option is not enabled;
 	 
-A [sample name]									default: "default";
 	 
-D [read group ID]									default: [read file 1];
 	 
-R [read group option]								Assign the read group option
 													default: this option is not enabled;
 	 
-p                                                  Output MD string and NM tag.

-s [max # of mismatches]							any integer from 0 to 4,
 													disable dynamic programming and
 													perform alignment with mismatches only.
 													If max # of mismatches is not specified,
 													default number is 3 for reads of length >= 50;
 													or 2 otherwise.
 													This option is not recommended.

```

### <a name="list-bam"></a>4. Multiple sets of paired-end reads in BAM format.

Syntax:

```bash
  % ./soap3-dp pair-multi [ref seq file].index [info list file] -bam [options]
```

[info list file] contains a list of read information.

```
Each line inside the file contains information of each set of paired-end reads.
 	 
Column 1: BAM read file name;						Column 2: minimum value of insert size;
Column 3: maximum value of insert size;				Column 4: output prefix;
Column 5 (optional): read group ID;					Column 6 (optional): sample name;
Column 7 (optional): read group options.	 
 	 
Delimiter is a tab character and only the following three cases are allowed:
Case 1: Only contains the first 4 columns.
Case 2: Only contains the first 6 columns.
Case 3: Contains all the 7 columns.
```

```
Options:
-L [length of the longest read in the input]		default: 120; normally in the range [75,200];
 													if < 100, better using GPU with 6GB memory
 	 
-h [output option]									1: all valid alignments;
 													2: all best alignment (default);
 													3: unique best alignment;
 													4: random best alignment
 	 
-b [output format]									1: Succinct;
 													2: SAM (default);
 													3: BAM
 	 
-c [GPU device ID]									Specify the GPU device for running SOAP3-dp;
 	 
-I													The input is in the Illumina 1.3+ FASTQ-like format
 													default: this option is not enabled;
 	 
-A [sample name]									default: "default";
 	 
-D [read group ID]									default: [read file 1];
 	 
-R [read group option]								Assign the read group option
 													default: this option is not enabled;
 	 
-p                                                  Output MD string and NM tag.

-s [max # of mismatches]							any integer from 0 to 4,
 													disable dynamic programming and
 													perform alignment with mismatches only.
													If max # of mismatches is not specified,
 													default number is 3 for reads of length >= 50;
 													or 2 otherwise.
 													This option is not recommended.
 
```

### <a name="single"></a>5. One set of single-end reads: Input reads from one file.

Syntax:
For read file [read file 1] in FASTA/FASTQ/GZIP format:

```bash
  % ./soap3-dp single [reference seq index] [read file 1] [options]
```

For read file [read file 1] in BAM format:

```bash
  % ./soap3-dp single [reference seq index] [read file 1] -bam [options]
```

```
Options:
-L [length of the longest read in the input]		default: 120; normally in the range [75,200];
 													if < 100, better using GPU with 6GB memory
 	 
-h [output option]									1: all valid alignments;
 													2: all best alignment (default);
 													3: unique best alignment;
 													4: random best alignment
 	 
-b [output format]									1: Succinct;
 													2: SAM (default);
 													3: BAM
 	 
-o [output prefix]									The output file names:
 													[output prefix].dpout.1, [output prefix].gout.1, [output prefix].gout.2, ...
 													default: [read file 1];
 	 
-c [GPU device ID]									Specify the GPU device for running SOAP3-dp;
 	 
-I													The input is in the Illumina 1.3+ FASTQ-like format
 													default: this option is not enabled;
 	 
-A [sample name]									default: "default";
 	 
-D [read group ID]									default: [read file 1];
 	 
-R [read group option]								Assign the read group option
 													default: this option is not enabled;
 	 
-p                                                  Output MD string and NM tag.

-s [max # of mismatches]							any integer from 0 to 4,
 													disable dynamic programming and
 													perform alignment with mismatches only.
 													If max # of mismatches is not specified,
 													default number is 3 for reads of length >= 50;
 													or 2 otherwise.
 													This option is not recommended.
```

Example 1: The read file query.fq contains a set of length-100 single reads, aligned with default options. All best alignments are reported for each read.

```bash
  % ./soap3-dp single genome.fa.index query.fq -L 100
```

Example 2: Suppose that the read file query.bam is in BAM format and all valid alignments are needed.

```bash
  % ./soap3-dp single genome.fa.index query.bam -bam -h 1
```

### <a name="list-single"></a>6. Multiple sets of single-end reads

Syntax:

For the files in FASTA/FASTQ/GZIP format:

```bash
  % ./soap3-dp single-multi [ref seq file].index [info list file] [options]
```

For the files in BAM format:

```bash
    % ./soap3-dp single-multi [ref seq file].index [info list file] -bam [options]
```

[info list file] contains a list of read information.

```
Format of the file [info list file]:
Each line inside the file contains information of each read.
 	 
Column 1: read file name;							Column 2: output prefix;
Column 3 (optional): read group ID;					Column 4 (optional): sample name;
Column 5 (optional): read group options.	 
 	 
Delimiter is a tab character and only the following three cases are allowed:
Case 1: Only contains the first 2 columns.
Case 2: Only contains the first 4 columns.
Case 3: Contains all the 5 columns.
```

```
Options:
-L [length of the longest read in the input]		default: 120; normally in the range [75,200];
 													if < 100, better using GPU with 6GB memory
 	 
-h [output option]									1: all valid alignments;
 													2: all best alignment (default);
 													3: unique best alignment;
 													4: random best alignment
 	 
-b [output format]									1: Succinct;
 													2: SAM (default);
 													3: BAM
 	 
-c [GPU device ID]									Specify the GPU device for running SOAP3-dp;
 	 
-I													The input is in the Illumina 1.3+ FASTQ-like format
 													default: this option is not enabled;
 	 
-A [sample name]									default: "default";
 	 
-D [read group ID]									default: [read file name];
 	 
-R [read group option]								Assign the read group option
 													default: this option is not enabled;

-p                                                  Output MD string and NM tag.

-s [max # of mismatches]							any integer from 0 to 4,
 													disable dynamic programming and
 													perform alignment with mismatches only.
 													If max # of mismatches is not specified,
 													default number is 3 for reads of length >= 50;
 													or 2 otherwise.
 													This option is not recommended.

```

## <a name="Section2.3"></a>2.3 About the output

In the "default" setting:

  SOAP3-dp divides the reads among 3 threads. Each thread outputs the alignment results in a separate file; the files are named [read file 1].gout.1, [read file 1].gout.2, and [read file 1].gout.3, etc.
  
  The alignment results computed by dynamic programming are in another separate file named [read file 1].dpout.1.
  
  If the input file contains a list of read files, then for each file name in the list, there will be a set of output files.

Succinct plain format:

  By default, each alignment is reported in succinct plain text format with 5 fields: Read #, Chromosome Id#, Offset, Strand, # of Mismatch.
  
  Regarding the alignments in the file [read file 1].dpout.1, there are two types: Type A - reported by SOAP3 with # of mismatch being shown; Type B - reported by dynamic programming with DP score being shown. Each alignment is reported with 7 fields: Read#, Chromosome Id#, Offset, Strand, # of Mismatch or DP Score, [0: Type A; 1: Type B], CIGAR-string**.

SAM format:

  When the option -b 2 is selected, each alignment is reported in SAM format.

In SOAP3-dp version 2.0, the mapping quality score is included in the output SAM format if all-valid and all-best output option is selected.

BAM format:

  In SOAP3-dp version 2.0, when the option -b 3 is selected, each alignment is reported in BAM format.

\* The mapping quality score of an alignment ranges from 0 to 40. The higher the score, the more reliable the alignment is. When a read is found to have no hit (alignment), a null alignment with score zero is outputted. In general, an alignment is given a score only if it is the best alignment; its score is determined by the following three factors: the number of best alignments, the number of mismatches/insertions/deletions, and the base quality values of the mismatched positions. When the best alignment is unique, it achieves the maximum score 40 if it is an exact matching, and around 35 if there is one mismatch; the score gradually decreases to 10 when the number of mismatches/INDELs increases to five. The score of a best but not unique alignment is usually very low (less than 5).

\*\* The CIGAR-string reported is same as defined in SAM format.

# <a name="Section3"></a>3. Running SOAP3-dp on a machine with multiple GPU devices

In SOAP3-dp version 2.0, one can run more than one instance of SOAP3-dp in the same machine if the machine has more than one GPU card. By default, each instance of SOAP3-dp requires one GPU and four CPU threads (i.e. one for controlling GPU and three for outputting the alignment results).

To save the total memory usage, it is highly recommended to set the parameter "ShareIndex" to 1 inside the file "soap3-dp.ini" so that the index in the memory (in CPU side) can be shared, and also it is required to enlarge the "kernel memlock size" of the system to at least 15G. To do so, one may use root account to run the command: "echo '* - memlock unlimited' >> /etc/security/limits.conf" and then re-login to the system again. Note that after using this command, the "kernel memlock size" of the system will be set to 15G permanently.

Then when starting each instance of SOAP3-dp, one needs to specify different GPU device by using an option "-c [GPU device ID]". For example, there are two GPU cards (with device ID 0 and device ID 1) in the machine, and assume that there are at least 8 CPU cores. To start two instances of SOAP3-dp and align two sets of paired-end reads in parallel: (1) [read file 1] and [read file 2], and (2) [read file A] and [read file B]:

```bash
  % ./soap3-dp pair [ref seq file].index [read file 1] [read file 2] -c 0 > out1.log &
  % ./soap3-dp pair [ref seq file].index [read file A] [read file B] -c 1 > out2.log &
```

Note that the files "out1.log" and "out2.log" store the output messages from the first instance and the second instance of SOAP3-dp respectively.

# <a name="Section4"></a>4. Fine tuning via configuration file soap3-dp.ini

The file soap3-dp.ini contains a number of parameters that can fine tune the performance of the aligner program soap3-dp. For single-user systems, soap3-dp.ini and the program soap3-dp are usually put in the same directory. For multi-user systems, the program soap3-dp can be put in a common directory (say, /usr/bin). The following two ways of invoking the alignment program will refer to different ini files: (1) /usr/bin/soap3-dp ... : It refers to the common ini file /usr/bin/soap3-dp.ini. (2) soap3-dp ... : Here it is assumed that the PATH variable has included the common directory /usr/bin (for locating the alignment program) and soap3-dp.ini in the current working directory is the default ini file.\*\*

> NumOfCpuThreads (default 4)  

By default soap3-dp assumes a quad-core CPU and uses three CPU threads for output purpose. The file "soap3-dp.ini" contains a parameter "NumOfCpuThreads" with default value 3. If the CPU has x > 1 core, the number 3 should be changed to x -1.

If there is more than one GPU card installed in a machine, one may run more than one instance of soap3-dp at the same time. Each soap3-dp instance needs one GPU card and (by default) 3 CPU threads for outputting the results.

> StrandArrangement (default +/-):

For paired-end read alignment, by default, soap3-dp only considers +/- as the proper strands of a paired-end alignment (i.e. forward strand for the alignment in left and reverse strand for the alignment in right). One can set the parameter "StrandArrangement" to "+/+" or "-/-".

> MaxOutputPerRead and MaxOutputPerPair (default 1000):

These two parameters are used in the all-valid-alignment mode to limit the number of alignments per read or per pairs. MaxOutputPerRead is for single-end data, and MaxOutputPair is for paired-end data.

> MatchScore (default 1), MismatchScore (default -2), GapOpenScore (default -3), GapExtendScore (default -1):

These four parameters are the scorings of single match, single mismatch, gap opening, and gap extension in the dynamic programming module. MatchScore can be any integer greater than or equal to 1; the other scores can be any integers in the range [-10,-1].

> DPScoreThreshold:

This parameter defines the score threshold when performing dynamic programming. By default, its value is equal to 0.3 x read length. One can specify an integer (no less than 0) for the threshold used in dynamic programming.

> BWALikeScore (default 1):

This parameter defines whether BWA-like mapping quality score is reported. By default, this is enabled (by setting BWALikeScore to 1). This option is useful for some software, like GATK, which was designed for BWA''s MAPQ scoring score. If this option is enabled, for single-end alignment, score will be in the range [0, 37], while for paired-end alignment, score will be in the range [0, 60].

> MinMAPQ and MaxMAPQ:

These two parameters define the minimum value and the maximum value (i.e. the range) of the mapping quality score for each alignment in SAM format, if BWA-like mapping quality score is disabled (by setting BWALikeScore to 0). By default, the mapping quality score is in the range between 1 and 40.

> ShareIndex (default 0):

This parameter defines whether the index in the memory (in CPU side) is shared when more than one instance of SOAP3-dp is running. By default, this option is not enabled. If sharing of the index is enabled (i.e. by setting ShareIndex to 1), it is required to enlarge the kernel memlock size of the system to 15G (one may check the Section 3 in this manual for more information).

> MaxLenReadName (default 64):

This parameter defines the maximum length allowed for the read descriptions inside the input query file. Larger the value, more memory space will be needed. For example, by setting MaxLenReadName=256, the program needs 25G memory in total using full SA. By default, MaxLenReadName=64, and the program needs 20G memory. The too-long read description will be trimmed.

The parameters used for the dynamic programming can be updated according to different requirement. For example, for a very tight alignment of allowing at most 4 mismatches/insertions/deletions, assuming the read length is 100, one may set the match score as 1, mismatch/gap-open/gap-extension score as -10, and the threshold as 55.

\*\* [ Due to historical reason,there is another ini file soap3-dp-builder.ini, which is needed by the index builder program. Users are not recommended to change the parameters in this ini file. Note that soaps3builder.ini and the program soap3-dp-builder can be kept in the same or different directory the same way as soap3-dp.ini and soap3-dp. See above. ]
