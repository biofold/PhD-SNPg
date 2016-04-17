#!/bin/bash
path=`dirname $0`
if [[ -z "$var" ]]
then
  hg='all'
else
  hg=$2
fi


if [[ '$hg' == 'all' || '$hg' == 'hg19' ]]
then
   /bin/sh $path/hg19/get_hg19.sh
fi

if [[ '$hg' == 'all' || '$hg' == 'hg38' ]]
then
   /bin/sh $path/hg38/get_hg38.sh
fi
