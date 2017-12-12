#!/bin/bash
path=`dirname $0`
# Expected files in this directory
# - hg38.2bit
# - hg38.phyloP7way.bw
# - hg38.phyloP100way.bw
# from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/

echo '- Download hg38.2bit'
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit -O $path/hg38.2bit
echo '- Download hg38.phyloP7way.bw'
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP7way/hg38.phyloP7way.bw -O $path/hg38.phylo7way.bw
echo '- Download hg38.phyloP100way.bw'
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw -O $path/hg38.phyloP100way.bw
