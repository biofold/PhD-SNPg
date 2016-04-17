# PhD-SNPg


INTRODUCTION
      
      Emidio Capriotti, 2016
      Scripts are licensed under a Creative Commons by NC-SA

      PhD-SNPg is a program for the annotation of single nucleotide variants
      that uses data from the UCSC repository.


INSTALLATION

      Minimum requirements:
      wget, zcat, scikit-learn

      Run:
        python setup.py install arch_type

      For Linux 64bit architecture there are two compiled versions:
        - linux.x86_64
        - linux.x86_64.v287

      Installation time depends on the network speed.
      About 30G UCSC files need to be downloaded.

      Test:
        python setup.py test	



MANUAL INSTALLATION

      1) Download phd-snp-g script from github
        - git clone https://github.com/ecapriotti/phd-snp-g

      2) Required python library: scikit-learn
        - git://github.com/scikit-learn/scikit-learn

      3) Required UCSC tools and data:
        - bigWigToBedGraph and twoBitToFa from
          http://hgdownload.cse.ucsc.edu/admin/exe
          in ucsc/exe directory

        - For hg19 based predictions:
          hg19.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
          hg19.phyloP46way.primate.bw http://snps.biofold.org/phd-snp-g/ucsc/hg19/hg19.phyloP46way.primate.bw	
          hg19.100way.phyloP100way.bw: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
          in ucsc/hg19 directory
          Alterantively hg19 bundle is available at http://snps.biofold.org/phd-snp-g/ucsc/hg19.tar.gz		

        - For hg38 based predictions:
          hg38.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
          hg38.phastCons7way.bw http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.bw
          hg38.phastCons100way.bw http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw
          in ucsc/hg38 directory
          Alterantively hg38 bundle is available at http://snps.biofold.org/phd-snp-g/ucsc/hg38.tar.gz


HOW TO RUN
		
      PhD-SNPg can take in input a single variation or a file containing multiple 
      single nucleotide variants.

      - For single variants use the option -c:
        python predict_variants.py chr7,158715219,A,G -g hg19 -c

      - For input file the input can be either: 
	
        plain tab separated file with 4 columns: chr, position, ref, alt
        python predict_variants.py test/test_variants_hg38.tsv -g hg38
       
        vcf file with in the firt 5 columns: chr, position, rsid, ref, alt  
        python predict_variants.py test/test_variants_hg19.vcf.gz --vcf -g hg19
