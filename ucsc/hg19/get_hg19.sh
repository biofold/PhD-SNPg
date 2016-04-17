#!/bin/bash
# Expected files in this directory 
# - hg19.2bit
# - hg19.phyloP46way.primate.bw
# - hg19.100way.phyloP100way.bw
# from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/

path=`dirname $0`

echo '- Download hg19.2bit'
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit -P $path/hg19.2bit
echo '- hg19.phyloP46way.primate.bw'
wget http://snps.biofold.org/PhD-SNPg/ucsc/hg19/hg19.phyloP46way.primate.bw -P $path/hg19.phyloP46way.primate.bw
echo '- Download hg19.100way.phyloP100way.bw'
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw -P $path/hg19.100way.phyloP100way.bw

