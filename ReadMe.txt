
SOAP3-DP version 2.3.167 (Oct 28 2011)
===============================

Old versions of SOAP3: SOAP3 version 91 (April 2011), maintained by BGI http://soap.genomics.org.cn/soap3.html;  
                       SOAP3 version 112 (July 2011), maintained by NIH http://biowulf.nih.gov/apps/soap3.html


1. Introduction
===============

SOAP3 is a GPU-based software for aligning short reads to a reference sequence.  It can find all alignments with k mismatches, where k is chosen from 0 to 4. When compared with its previous version SOAP2, SOAP3 can be up to tens of times faster.  For example, when aligning length-100 reads with the human genome, SOAP3 is the first software that can find all 4-mismatch alignments in tens of seconds per million reads.

The alignment program in this package is optimized to work for multi-millions of short reads each time by running a multi-core CPU and the GPU concurrently.

To exploit the parallelism of the GPU effectively, SOAP3 is using an adapted version of the 2BWT index of SOAP2 (the new index is called the GPU-2BWT). The index and alignment algorithms were developed by the algorithms research group of the University of Hong Kong (T.W. Lam, C.M. Liu, Thomas Wong, Edward Wu and S.M. Yiu), in collaboration with BGI (Beijing Genome Institute).



2. Hardware & Platform
=======================

To run SOAP3, you need a linux workstation equipped with

   (i)   a multi-core CPU with at least 20GB main memory (default quad-core) * and

   (ii)  a CUDA-enabled GPU with at least 3GB memory (default 6GB). *

SOAP3 has been tested with the following GPU: NVIDIA Tesla C2070 (6GB memory), Tesla M2050 (3GB memory), GTX 580 (3GB memory).  
And it should also work using Tesla C1060, C2050, C2070, Quadro 6000, GeForce GTX 590 (3 GB).

SOAP3 was developed under the 64-bit linux platform and the CUDA Driver version 3.2.

* See Section 4 for how to change the default with the INI file soap3-dp.ini.


3. Usage
========

SOAP3 consists of 3 parts: (i) Index (2BWT and GPU-2BWT) builder; (ii) Aligner; (iii) Output viewer.

After downloading the compressed file, run the following commands to extract the different programs:

% gunzip soap3-r146-x64-linux.tar.gz
% tar -xvf soap3-r146-x64-linux.tar

You should get ten files: 

soap3-dp-builder, soap3-dp-builder.ini, BGS-Build, soap3-dp, soap3-dp.ini, BGS-View, BGS-View-PE, make_view.sh, make_view_sam.sh make_view_simple.sh


3.1 Index builder
=================

Soap3-builder preprocesses the FASTA reference sequence(s) given in a file <ref sequence file> to produce a few indexes for alignment. The indexes will be stored in files whose names are in the form <ref sequence file>.index.***. Note the following restrictions to the reference sequences:

(i) At most 256 reference sequences (chromosomes) in one single FASTA file.
(ii) The total length of all reference sequences is at most 4 billion.
(ii) All invalid characters (i.e., anything other than A, C, G, T) will be replaced by character G; regions that contain more than 10 invalid characters will not be involved in any alignment.


Step 1: Build the 2BWT index.

   Syntax:

   % ./soap3-dp-builder <ref sequence file>

   For example:

   % ./soap3-dp-builder genome.fa

 
Step 2: Convert the 2BWT index to the GPU2-BWT index.

   Syntax:

   % ./BGS-Build <ref sequence file>.index

   For example:

   % ./BGS-Build genome.fa.index

Additional files with the filename-prefix "<ref sequence file>.index" will be generated. 

With all these index files, you are now ready to use soap3-dp to perform alignment for the reads.



3.2 Aligner
===========

The aligner supports query files in FASTA format or FASTQ format. It will detect the format automatically.


*  For paired-end reads: Input reads from two files <query file 1> and <query file 2>.   

   Syntax:
   % ./soap3-dp pair <ref seq file>.index <query file 1> <query file 2> -u <max insert size> -v <min insert size> [options]

   <max insert size>, <min insert size> : normally set to the average insert size plus 3 times of the standard deviation, 
   and minus 3 times of the standard deviation, respectively. 

   Options:
         -m <max # of mismatches> (any integer from 0 to 4, 
                                   default: 2 for reads with length < 50; 3 for reads with length >= 50)
         -L <length of the longest read in the input> (default: 120; normally in the range [75,200]; if < 100, better using GPU with 6GB memory))
         -h <alignment type> (1: all valid alignments; 
                              4: random best alignment (default))
         -b <output format>  (0: Succinct [Binary format] (fastest);
                              1: Simple [plain text] (default); 
                              2: SAM v1.4 [plain text]) 

         Example: Assume query files query1.fa and query2.fa, and insert size range [200, 500].

            % ./soap3-dp pair genome.fa.index query1.fa query2.fa -u 500 -v 200 -m 4

         See Section 3.3 about the details of output format.


* For single-end reads: Input reads from one file <query file 1>.  

   Syntax:
   % ./soap3-dp single <reference seq index> <query file 1> [options]

   options:
         -m <max # of mismatches> (any integer from 0 to 4, 
                                   default: 2 for reads with length < 50; 3 for reads with length >= 50)
         -L <length of the longest read in the input> (default: 120; normally in the range [75,200])
         -h <alignment type> (1: all valid alignments; 
                              2: all best alignments;
                              3: unique best alignments;
                              4: random best alignments (default))
         -b <output format>  (0: Succinct [Binary format];
                              1: Simple [plain text] (default); 
                              2: SAM v1.4 [plain text]) 

         Example 1: The query file query.fa contains a set of length-100 single reads, aligned with default options (up to 3 mismatches).
                    A best alignment is reported for each read (break tie randomly).

               % ./soap3-dp single genome.fa.index query.fa -L 100

         Example 2: Suppose that the maximum number of mismatches allowed is 4, and all valid alignments are needed.  

               % ./soap3-dp single genome.fa.index query.fa -m 4 -h 1





3.3 About the output
====================

In the "default" setting,

  * soap3-dp divides the reads among 3 threads. Each thread outputs the alignment results in a separate file;
    the files are named <query file 1>.gout.1, <query file 1>.gout.2, and <query file 1>.gout.3, etc.
  * Each alignment is reported in simple plain text format with 5 fields: Read #, Chromosome Id#, Offset,  Strand, # of Mismatch.

One can use the option -b 2 to output in SAM format, and -b 1 in binary format (the latter is for efficiency
and requires a viewer; see Section 3.3.2 below).


3.3.1 Merging alignment results in simple plain format (default) into one file <query file1>.out.

    The program make_view_simple.sh merges all (default 3) output files 
    <query file 1>.gout.1, <query file 1>.gout.2, and <query file 1>.gout.3
    into one file <query file1>.out.

    Syntax (for both single-end and paired-end input):

      % ./make_view_simple.sh <query file 1>

    For example:
    % ./make_view_simple.sh query1.fa


3.3.2 Merging alignment results in SAM format into one file <query file 1>.out

    If the SAM format is chosen, then soap3-dp outputs the results in the SAM format into 3 (default) files.
    The program make_view_sam.sh merges all output files into one.

     Syntax (for both single-end and paired-end input):

     % ./make_view_sam.sh <query file 1>


3.3.2  Alignment results in binary format (for efficiency when outputting large number of alignments)

    If the binary output option is chosen, soap3-dp outputs in binary format, which is not readable. 
    The program make_view_binary.sh converts the binary format into plain text format and
    merges all the output files into one, namely, <query file 1>.out

   Syntax - singe-end reads:

      % ./make_view_binary.sh single <query file 1>

   Syntax - paired-end read alignment:

      % ./make_view_binary.sh pair <query file 1> <query file 2>

    No matter it is the single-end or paired-end setting, the result are stored in a singe file <query file 1>.out.  


         

4. Fine tuning via configuration file soap3-dp.ini

NumOfCpuThreads (default 3): 
   By default soap3-dp assumes a quad-core CPU and uses three CPU threads for output purpose.  
   The file "soap3-dp.ini" contains a parameter "NumOfCpuThreads" with default value 3.
   If the CPU has x > 1 core, the number 3 should be changed to x -1.   

GPUMemory (default 6):
   If a GPU with less 6GB memory is used, one can modify the value of the parameter "GPUMemory" 
   from 6 to 3, 4 or 5 accordingly.

StrandArrangement (default +/-):
   For paired-end read alignment, by default, soap3-dp only considers +/- as the proper strands of a paired-end alignment 
   (i.e. forward strand for the alignment in left and reverse strand for the alignment in right). 
   One can set the parameter "StrandArrangement" to "+/+" or "-/-". 

MaxOutputPerRead and MaxOutputPerPair (default 1000)
   These two parameters are used in the all-valid-alignment mode to limit the number of alignments per read or per pairs. 
   MaxOutputPerRead is for single-end data, and MaxOutputPair is for MaxOutputPerPair.
   If these parameters are set to a huge value (say, 100K), it is recommended to output in binary format to save time.

                 
5. Performance & Reference

Below is the performance of SOAP3 (version 135) when running on GTX 580 (3.2 GB RAM) to align 
length-100 reads to the human genome with 4 mismatches allowed.

  *  Human reference genome: 37.1
  *  Read data: Sequence read archive (SRA) accession #: SRR211279: volume 25 M (x 2), read length 100, average insert size 300, SD 30    
  *  Index loading time: 134 seconds; use option -m 4 to align with up to 4 mismatches.

     =======================================================================================================
                                           |       Single-end reads        |         Paired-end reads
                4 mismatches (-m 4)        |  Read loading  +  Alignment   |   Read loading  +  Alignment
                                           |     time            time      |      time            time
                                           |  per M reads      per M reads |   per M pairs      per M pairs
     =======================================================================================================
     Random best alignment                 |      ~2  sec   +  ~5 sec      |  ~4 sec +  ~20 sec
        (default output plain text)        |                               |

     ALL valid alignments (default <= 1000)|                               |
        output default plain text          |      ~2  sec   +  ~35 sec     |  ~4 sec +  ~77 sec
        output in binary format  -b 0      |      ~2  sec   +  ~35 sec     |  ~4 sec +  ~74 sec
                                           |                               |
                               
    ALL valid alignments (unlimited)       |                               |  
        output default plain text          |      ~2 sec    +  ~62 sec     |  ~4 sec +  ~84 sec
        output in binary format  -b 0      |      ~2 sec    +  ~38 sec     |  ~4 sec +  ~75 sec

     =======================================================================================================

Reference:
  SOAP3: GPU-based Compressed Indexing and Ultra-fast Parallel Alignment of Short Reads. MASSIVE 2011. http://i.cs.hku.hk/~twlam/soap3.pdf

Contact: 
  BGI and HKU Algorithm group




6. Remark

SAM-tools v0.1.18 is included in SOAP3-GPU package to facilitate outputting alignment result into SAM output format. We have slightly modified the original code of SAM-tools to make it compilable under g++. Please see http://samtools.sourceforge.net/ for details of this package.                                                           




7. Copyright


Copyright (C) 2011, Computer Science Department, The University of Hong Kong

SOAP3 is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
 
SOAP3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should find a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

