# PhD-SNPg


INTRODUCTION
      
      Emidio Capriotti, 2016.
      Scripts are licensed under the Creative Commons by NC-SA license.

      PhD-SNPg is a program for the annotation of single nucleotide variants that uses data from the UCSC repository.


INSTALLATION

      Minimum requirements:
      wget, curl, zcat, scikit-learn.

      Run:
        python setup.py install arch_typ

      For Linux 64bit architecture there are two compiled versions:
        - linux.x86_64
        - linux.x86_64.v287
        - macOSX.x86_64

      Installation time depends on the network speed.
      About 30G UCSC files need to be downloaded.

      For light installation:
        python setup.py install arch_typ --web

      With the web_install option the PhD-SNPg
      runs without downloading the UCSC data.
      The functionality of the program depends 
      on the network speed.

      Test:
        python setup.py test

      For web installation:
        python setup.py test --web


MANUAL INSTALLATION

      1) Download PhD-SNPg script from github
        - git clone https://github.com/biofold/PhD-SNPg.git

      2) Required python library: scikit-learn-0.17
          It is available in tools directory or at
          https://pypi.python.org/simple/scikit-learn/

	- Untar the scikit-learn-0.17.tar.gz 
          move in the directory and run
          python setup.py install --install-lib=..
          
      3) Required UCSC tools and data:
        - bigWigToBedGraph and twoBitToFa from
          http://hgdownload.cse.ucsc.edu/admin/exe
          in ucsc/exe directory

        - For hg19 based predictions:
          hg19.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
          hg19.phyloP46way.primate.bw http://snps.biofold.org/PhD-SNPg/ucsc/hg19/hg19.phyloP46way.primate.bw	
          hg19.100way.phyloP100way.bw: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
          in ucsc/hg19 directory
          Alterantively hg19 bundle is available at http://snps.biofold.org/PhD-SNPg/ucsc/hg19.tar.gz		

        - For hg38 based predictions:
          hg38.2bit: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
          hg38.phastCons7way.bw http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.bw
          hg38.phastCons100way.bw http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw
          in ucsc/hg38 directory
          Alterantively hg38 bundle is available at http://snps.biofold.org/PhD-SNPg/ucsc/hg38.tar.gz


HOW TO RUN
		
      PhD-SNPg can take in input a single variation or a file containing multiple single nucleotide variants.
      For web installation append --web at the of all commands.

      - For single variants use the option -c:
        python predict_variants.py chr7,158715219,A,G -g hg19 -c

      - For input file the input can be either: 
	
        plain tab separated file with 4 columns: chr, position, ref, alt
        python predict_variants.py test/test_variants_hg38.tsv -g hg38
       
        vcf file with in the firt 5 columns: chr, position, rsid, ref, alt  
        python predict_variants.py test/test_variants_hg19.vcf --vcf -g hg19


OUTPUT

      PhD-SNPg returns in output: 

        PREDICTION: Pathogenic or Benign
        SCORE: a probabilistic score between 0 and 1. If the score is >0.5 the variants is predicted to be Pathogenic.
        FDR: The false discovery rate associated to higher/lower SCORE.
        PhyloP100: PhyloP100 in the mutated position.
        AvgPhyloP100: Average value of PhyloP100 in a 7-nucleotide window around the mutated position.

        The scores added as extra columns to the input file. An example of output is reported below.

        #CHROM  POS     REF     ALT     CODING  PREDICTION      SCORE   FDR   PhyloP100       AvgPhyloP100
        1       45331676        G       A       Yes     Pathogenic      0.990   0.022   7.723   3.313
        1       237634938       G       T       Yes     Benign  0.005   0.011   -4.248  5.939
        2       26461838        G       A       Yes     Benign  0.011   0.022   -0.064  6.116
        2       166009835       A       G       Yes     Benign  0.176   0.072   0.591   4.188
        2       174753570       G       C       Yes     Pathogenic      0.777   0.087   -1.937  5.439
        5       44305045        G       A       Yes     Pathogenic      0.985   0.026   1.655   3.964

