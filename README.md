# PhD-SNPg


INTRODUCTION
      
      Emidio Capriotti, 2016.
      Scripts are licensed under the Creative Commons by NC-SA license.

      PhD-SNPg is a program for the annotation of single nucleotide variants that uses data from the UCSC repository.


INSTALLATION

      Minimum requirements:
      wget, zcat, scikit-learn.

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

      1) Download PhD-SNPg script from github
        - git clone https://github.com/biofold/PhD-SNPg.git

      2) Required python libraries: scikit-learn-0.17
	  They are already available in tools directory.

	- Untar the scikit-learn-0.17.tar.gz directory and run
          python setup.py install --install-lib=../
	  https://pypi.python.org/simple/scikit-learn/

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

      - For single variants use the option -c:
        python predict_variants.py chr7,158715219,A,G -g hg19 -c

      - For input file the input can be either: 
	
        plain tab separated file with 4 columns: chr, position, ref, alt
        python predict_variants.py test/test_variants_hg38.tsv -g hg38
       
        vcf file with in the firt 5 columns: chr, position, rsid, ref, alt  
        python predict_variants.py test/test_variants_hg19.vcf.gz --vcf -g hg19


OUTPUT

      PhD-SNPg returns in output: 

        PREDICTION: Pathogenic or Benign
        SCORE: a probabilistic score between 0 and 1. If the score is >0.5 the variants is predicted to be Pathogenic.
        FDR: The false discovery rate associated to higher/lower SCORE.
        PhyloP100: PhyloP100 in the mutated position.
        AvgPhyloP100: Average value of PhyloP100 in a 7-nucleotide window around the mutated position.

        The scores added as extra columns to the input file. An example of output is reported below.

        #CHROM  POS     REF     ALT     PREDICTION      SCORE   FDR   PhyloP100       AvgPhyloP100
        1       10042376        C       G       Benign  0.160   0.071   -0.159  4.509
        1       197094291       C       T       Pathogenic      0.990   0.486   7.304   3.451
        2       31751295        G       A       Pathogenic      0.901   0.322   1.810   3.205
        2       71797809        C       T       Pathogenic      0.918   0.342   1.181   3.131
        2       179577870       T       C       Benign  0.028   0.023   -6.363  2.406
        5       74046464        C       T       Benign  0.040   0.029   -0.070  4.422

