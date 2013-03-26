#!/bin/bash

#############################################################################################
# This script is to concatenate all alignment outputs of SAM format into one file.          #
#############################################################################################

if [ "$#" -lt 1 ];
then
    echo "For single read alignment:"
    echo "     Usage: $0 [input query file]"
    echo "     For example: $0 query.fa"
    echo "."
    echo "For paired-end read alignment:"
    echo "     Usage: $0 [input query file 1]"
    echo "     For example: $0 query_1.fa"
    echo "."
    echo "The script will read the corresponding alignment files like [input query file].gout.1"
    echo "and combine them into an alignment result file."
    exit 1
fi

input_query=$1
		
if [ ! -f $input_query.gout.1 ];
then
	 echo "The corresponding file $input_query.gout.1 cannot be found"
else
	 cp $input_query.gout.1 $input_query.out
	 
	 for i in 2 3 4 5 6
	 do
		  if [ -f $input_query.gout.$i ];
		  then
                           egrep "^@" -v $input_query.gout.$i >> $input_query.out
		  fi
	 done
         if [ -f $input_query.dpout.1 ];
         then
                  egrep "^@" -v $input_query.dpout.1 >> $input_query.out
         fi
         if [ -f $input_query.unpair ];
         then
                  egrep "^@" -v $input_query.unpair >> $input_query.out
         fi
	 echo "------------------------------------------------------"
	 echo "Alignment file - $input_query.out has been created"
	 echo "------------------------------------------------------"
fi
