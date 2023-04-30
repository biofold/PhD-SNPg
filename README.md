# PhD-SNPg


INTRODUCTION
      
      Emidio Capriotti, 2016.
      Scripts are licensed under the Creative Commons by NC-SA license.

      PhD-SNPg is a program for the annotation of single nucleotide variants that uses data from the UCSC repository.


INSTALLATION

      Minimum requirements:
      wget, curl, zcat, scikit-learn.

      Run:
        python2 setup.py install arch_typ

      For Linux 64bit architecture there are two compiled versions:
        - linux.x86_64
        - linux.x86_64.v369
        - linux.x86_64.v385
        - macOSX.arm64
        - macOSX.x86_64

      Installation time depends on the network speed.
      About 30G UCSC files need to be downloaded.

      For light installation:
        python2 setup.py install arch_typ --web

      With the web_install option the PhD-SNPg
      runs without downloading the UCSC data.
      The functionality of the program depends 
      on the network speed.

      Test:
        python2 setup.py test

      For web installation:
        python2 setup.py test --web


MANUAL INSTALLATION

      1) Download PhD-SNPg script from github
        - git clone https://github.com/biofold/PhD-SNPg.git

      2) Required python2 library: scikit-learn-0.17
          It is available in tools directory or at
          https://pypi.python.org/simple/scikit-learn/

	- Untar the scikit-learn-0.17.tar.gz 
          move in the directory and run
          python2 setup.py install --install-lib=..
          
      3) Required UCSC tools and data:
        - bigWigToBedGraph, twoBitToFa and liftOver from
          http://hgdownload.cse.ucsc.edu/admin/exe
          in ucsc/exe directory

        - For hg19 based predictions:
          hg19.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
          hg19.phyloP46way.primate.bw: http://snps.biofold.org/PhD-SNPg/ucsc/hg19/hg19.phyloP46way.primate.bw	
          hg19.100way.phyloP100way.bw: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
          hg19ToHg38.over.chain.gz: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
          in ucsc/hg19 directory
          Alterantively hg19 bundle is available at http://snps.biofold.org/PhD-SNPg/ucsc/hg19.tar.gz		

        - For hg38 based predictions:
          hg38.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
          hg38.phyloP7way.bw:  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP7way/hg38.phyloP7way.bw
          hg38.phyloP100way.bw:  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw
          hg38.phyloP470way.bw:  http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP470way/hg38.phyloP470way.bw
          in ucsc/hg38 directory
          Alterantively hg38 bundle is available at http://snps.biofold.org/PhD-SNPg/ucsc/hg38.tar.gz


HOW TO RUN
		
      PhD-SNPg can take in input a single variation or a file containing multiple single nucleotide variants.
      For web installation append --web at the of all commands.

      - For input file the input can be either: 
	
        plain tab separated file with 4 columns: chr, position, ref, alt
        python2 predict_variants.py test/test_variants_hg38.tsv -g hg38
       
        vcf file with in the firt 5 columns: chr, position, rsid, ref, alt  
        python2 predict_variants.py test/test_variants_hg19.vcf --vcf -g hg19


OUTPUT

      PhD-SNPg returns in output:

        PREDICTION: Pathogenic or Benign
        SCORE: a probabilistic score between 0 and 1. If the score is >0.5 the variants is predicted to be Pathogenic.
        FDR: The false discovery rate associated to higher/lower SCORE.
        PhyloP100: PhyloP100 in the mutated position.
        AvgPhyloP100: Average value of PhyloP100 in a 5-nucleotide window around the mutated position.

        The scores added as extra columns to the input file. An example of output is reported below.

        #CHROM  POS     REF     ALT     CODING  PREDICTION      SCORE   FDR   PhyloP100       AvgPhyloP100
        1       197125161       C       T        Yes    Pathogenic      0.954   0.024   7.361   4.082
        2       31526225        G       A        Yes    Pathogenic      0.792   0.061   1.911   2.799
        3       46998228        T       G        No     Pathogenic      0.990   0.007   7.762   5.511
        4       663761          C       G        No     Benign  0.269   0.068   -0.211  -0.410
        5       74750639        C       T        Yes    Benign  0.006   0.004   -0.074  5.880
        11      121125583       A       G        Yes    Benign  0.001   0.003   -9.071  2.506
