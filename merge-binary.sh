#!/bin/bash

#############################################################################################
# This script is to transform the binary alignment output into a readable alignment result. #
#############################################################################################

if [ "$#" -lt 1 ];
then
    echo "For single read alignment:"
    echo "     Usage: $0 single [input query file]"
    echo "     For example: $0 single query.fa"
    echo "."
    echo "For paired-end read alignment:"
    echo "     Usage: $0 pair [input query file 1] [input query file 2]"
    echo "     For example: $0 pair query_1.fa query_2.fa"
    echo "."
    echo "The script will read the corresponding alignment files like [input query file].gout.1"
    echo "and transform them into readable alignment result files."
    exit 1
fi

if [ "$1" == "single" ];
then

		input_query=$2
		
		if [ ! -f $input_query.gout.1 ];
		then
			 echo "The corresponding file $2.gout.1 cannot be found"
		else
			 echo "#Read_Pos	Read_Name	Read	Chromosome	Offset	Read_Length	Strand	Type" > $input_query.out
			 if [ -f $input_query.gout.1 ];
			 then
				  ./BGS-View $input_query $input_query.gout.1 >> $input_query.out
			 fi
			 for i in 2 3 4 5 6
			 do
				  if [ -f $input_query.gout.$i ];
				  then
					   ./BGS-View $input_query $input_query.gout.$i >> $input_query.out
				  fi
			 done
                         if [ -f $input_query.dpout.1 ];
                         then
                                 ./BGS-View $input_query $input_query.dpout.1 >> $input_query.out
                         fi
			 echo "------------------------------------------------------"
			 echo "Alignment file - $input_query.out has been created"
			 echo "------------------------------------------------------"
		fi
		
fi

if [ "$1" == "pair" ];
then
		input_query_1=$2
		input_query_2=$3
		
		if [ ! -f $input_query_1.gout.1 ];
		then
			 echo "The corresponding file $2.gout.1 cannot be found"
		else
			 echo "#Read_Pos	Read_Name	Read	Chromosome	Offset	Read_Length	Strand	Type" > $input_query_1.out
			 if [ -f $input_query_1.gout.1 ];
			 then
				  ./BGS-View-PE $input_query_1 $input_query_2 $input_query_1.gout.1 >> $input_query_1.out
			 fi
			 for i in 2 3 4 5 6
			 do
				  if [ -f $input_query_1.gout.$i ];
				  then
					   ./BGS-View-PE $input_query_1 $input_query_2 $input_query_1.gout.$i >> $input_query_1.out
				  fi
			 done
                         if [ -f $input_query_1.dpout.1 ];
                         then
                                 ./BGS-View-PE $input_query_1 $input_query_2 $input_query_1.dpout.1 >> $input_query_1.out
                         fi
			 echo "------------------------------------------------------"
			 echo "Alignment file - $input_query_1.out has been created"
			 echo "------------------------------------------------------"
		fi
fi
