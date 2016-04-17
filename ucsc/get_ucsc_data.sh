#!/bin/bash
path=`dirname $0`
echo $path
if [[ -z "$var" ]]
then
  hg='all'
else
  hg=$2
fi


if [[ '$hg' == 'all' || '$hg' == 'hg19' ]]
then
  echo '- Download hg19.2bit'
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit -P $path/hg19/hg19.2bit
  echo '- hg19.phyloP46way.primate.bw'
  wget http://snps.biofold.org/phd-snp-g/ucsc/hg19/hg19.phyloP46way.primate.bw -P $path/hg19/hg19.phyloP46way.primate.bw
  echo '- Download hg19.100way.phyloP100way.bw'
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw -P $path/hg19/hg19.100way.phyloP100way.bw
fi

if [[ '$hg' == 'all' || '$hg' == 'hg38' ]]
then
  echo '- Download hg38.2bit'
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit -P $path/hg19/hg38.2bit
  echo '- Download hg38.phyloP7way.bw'
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP7way.bw -P $path/hg38/hg38.phylo7way.bw
  echo '- Download hg38.phyloP100way.bw'
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw -P $path/hg38/hg38.phyloP100way.bw
fi
