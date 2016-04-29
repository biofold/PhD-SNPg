#!/bin/sh
# Search position in sorted bed file with chromosomes in format chrN  
# Run tools/bed_range.sh chromosome position bed_file 

ichr=$1
pos=$2
bed=$3

pre=`echo $ichr |cut -b 1-3`

if [ "$pre" != "chr" ]
then
	ichr="chr"$ichr
fi

grep  -w "^"$ichr $bed 2>/dev/null |awk -v p=$pos 'BEGIN {r="NOT-CODING"} {if (p>=$2 && p<=$3) {r=$2"-"$3; exit}; if (p<$2) {exit} } END {print r}' 
