import os, sys
from commands import getstatusoutput


def test():
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	ucsc_dir = prog_dir+'/ucsc'
	ucsc_exe = ucsc_dir+'/exe'
	prog_cat = 'zcat'
	print >> sys.stderr, 'W'


def get_ucsc_tools(arch_type):
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	ucsc_dir = prog_dir+'/ucsc'
	twobit='twoBitToFa'
	bwg='bigWigToBedGraph'
	ucsc_tool='http://hgdownload.cse.ucsc.edu/admin/exe/'
	prog_get='wget'
	wtwobit=ucsc_tool+'/'+arch_type+'/'+twobit
	cmd=prog_get+' '+wtwobit+' -O '+ucsc_dir+'/exe/'+twobit+'; chmod a+x '+ucsc_dir+'/exe/*'
	print >> sys.stderr,'   Download twoBitToFa'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	wbwg=ucsc_tool+'/'+arch_type+'/'+bwg
	cmd=prog_get+' '+wbwg+' -O '+ucsc_dir+'/exe/'+bwg+'; chmod a+x '+ucsc_dir+'/exe/*'
	print >> sys.stderr,'   Download bigWigToBedGraph'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	print >> sys.stderr,'   Test tools'
	cmd=ucsc_dir+'/exe/twoBitToFa; '+ucsc_dir+'/exe/bigWigToBedGraph'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0 and out[0]!=65280:
		print >> sys.stderr,'  ERROR: Incorrect architecture',arch_type
		cmd='rm '+ucsc_dir+'/exe/*'
		outx=getstatusoutput(cmd)
	else:
		print >> sys.stderr,'   Downloaded UCSC Tools'
	return out
	

def get_ucsc_data(hg,namefile,odir,ucsc_dat='http://hgdownload.cse.ucsc.edu/goldenPath'):
	prog_dir = os.path.dirname(os.path.abspath(__file__))
        ucsc_dir = prog_dir+'/ucsc/'+hg
	prog_get = 'wget'
	print '\n   Download',namefile
	data = ucsc_dat+'/'+hg+'/'+odir+'/'+namefile
	cmd = prog_get+' '+data+' -O '+ucsc_dir+'/'+namefile
	print cmd
	out=getstatusoutput(cmd)
	print cmd
	if not os.path.isfile(ucsc_dir+'/'+namefile): 
		print >> sys.stderr,'ERROR: Not found ',namefile
		sys.exit(1)
	return out
	

def setup(arch_type,hg='all',web=False):
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	prog_get = 'wget'
	print '1) Check wget'
	cmd='which '+prog_get
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0:
		print >> sys.stderr, "ERROR: Command wget not available."
		print sys.exit(1)
	print '\n2) Compile scikit-learn-0.17 and check joblib'
        cmd='cd '+prog_dir+'/tools/; tar -xzvf scikit-learn-0.17.tar.gz;' 
	cmd=cmd+'cd scikit-learn-0.17; python setup.py install --install-lib='+prog_dir+'/tools'
	print cmd
        out=getstatusoutput(cmd)
        print out[1]
        if out[0]!=0:
                print >> sys.stderr, "ERROR: scikit-learn istallation failed."
                sys.exit(1)	
        cmd='python -c \'from sklearn.externals import joblib\''
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0:
		print >> sys.stderr, "ERROR: joblib library not available."
		sys.exit(1)
	print '\n3) Download UCSC Tools'
	out=get_ucsc_tools(arch_type)
	if out[0]!=0 and out[0]!=65280:
		print >> sys.stderr,'ERROR: Incorrect architecture check your system or compile it.'
		print sys.exit(1)
	
	if web: sys.exit(0)

	dcount=0
	print '\n4) Download UCSC Data. It can take several minutes depending on the newtork speed.'
	if (hg=='all' or hg=='hg19'):
		out=get_ucsc_data('hg19','hg19.2bit','bigZips')
		if out[0]==0: dcount+=1
		biofold='http://snps.biofold.org/PhD-SNPg/ucsc'
		out=get_ucsc_data('hg19','hg19.phyloP46way.primate.bw','',biofold)
		if out[0]==0: dcount+=1
		out=get_ucsc_data('hg19','hg19.100way.phyloP100way.bw','phyloP100way')
		if out[0]==0: dcount+=1
		if dcount<3:
			print >> sys.stderr, 'ERROR: Problem in hg19 data downloading.'
			sys.exit(1)
	dcount=0
	if (hg=='all' or hg=='hg38'):
		out=get_ucsc_data('hg38','hg38.2bit','bigZips')
		if out[0]==0: dcount+=1
		out=get_ucsc_data('hg38','hg38.phyloP7way.bw','phyloP7way')
		if out[0]==0: dcount+=1
		out=get_ucsc_data('hg38','hg38.phyloP100way.bw','phyloP100way')
		if out[0]==0: dcount+=1
		if dcount<3:
			print >> sys.stderr, 'ERROR: Problem in hg38 data downloading.'
			sys.exit(1)
	print   '   Downloaded UCSC data'


def test():
	hgs=['hg19','hg38']
	prog_dir = os.path.dirname(os.path.abspath(__file__))
        ucsc_dir = prog_dir+'/ucsc'
	test_dir = prog_dir+'/test'
	ucsc_tool =  ucsc_dir+'/exe'
	print '1) Test python libray scikit-learn'
	cmd='cd '+prog_dir+'/tools; python -c \'import joblib; print joblib.__doc__\''
	print cmd
	out=getstatusoutput(cmd+' |head -n 9')
	print out[1]
	if out[0]!=0:
		print >> sys.stderr,'ERROR: scikit-learn not installed'
		sys.exit(1)
	print '\n2) Test zcat command'
	cmd='zcat -f '+test_dir+'/test_variants_hg19.vcf.gz '
	print cmd
	out=getstatusoutput(cmd+' |grep -A 2 \'#CHROM\'')
	print out[1]
	if out[0]!=0:
		print >> sys.stderr,'ERROR: Command zcat not available'
		sys.exit(1)
	print '\n3) Check hg19 files'
        cmd='cd '+ucsc_dir+'/hg19/; md5sum -c hg19.md5'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0:
		print >> sys.stderr,'WARNING: md5sum failed on hg19'
	#files=['hg19.2bit','hg19.100way.phyloP100way.bw','hg19.phyloP46way.primate.bw']
	#for i in files:
	#	cmd='ls '+ucsc_dir+'/hg19/'+i
	#	print cmd
	#	out=getstatusoutput(cmd)
	#	print out[1]
	#	if out[0]!=0:		
	#		print >> sys.stderr,'WARNING: Not found file',i
	#		if 'hg19' in hgs: hgs.remove('hg19')
	print '\n4) Check hg38 files'
	cmd='cd '+ucsc_dir+'/hg38/; md5sum -c hg38.md5'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0:
		print >> sys.stderr,'WARNING: md5sum failed on hg38'
	#files=['hg38.2bit','hg38.phyloP100way.bw','hg38.phyloP7way.bw']
	#for i in files:
	#	cmd='ls '+ucsc_dir+'/hg38/'+i
	#	print cmd
	#	out=getstatusoutput(cmd)
	#	print out[1]
	#	if out[0]!=0:
	#		print >> sys.stderr,'WARNING: Not found file',i
	#		if 'hg38' in hgs: hgs.remove('hg38')
	if hgs==[]:
		print >> sys.stderr,'ERROR: UCSC data file not correctly downloaded'
		sys.exit()
	print '\n5) Test twoBitToFa command'
	twobit=ucsc_dir+'/'+hgs[0]+'/'+hgs[0]+'.2bit'
	cmd=ucsc_tool+'/twoBitToFa '+twobit+' stdout -seq=chr1 -start=10008 -end=10010'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0:
		print >> sys.stderr,'ERROR: twoBitToFa not working'
		sys.exit(1)
	print '\n6) Test bigWigToBedGraph command'
	if hgs[-1]=='hg19':
		bwg=ucsc_dir+'/'+hgs[-1]+'/'+hgs[-1]+'.100way.phyloP100way.bw'
	else:
		bwg=ucsc_dir+'/'+hgs[-1]+'/'+hgs[-1]+'.phyloP100way.bw'
	cmd=ucsc_tool+'/bigWigToBedGraph '+bwg+' stdout -chrom=chr1 -start=100008 -end=100012'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0:
		print >> sys.stderr,'ERROR: bigWigToBedGraph not working'
		sys.exit(1)
	print '\n7) Test predict_variants.py'
	cmd='python predict_variants.py test/test_variants_'+hgs[0]+'.tsv -g '+hgs[0]
	print cmd
	out=getstatusoutput(cmd+' |head -n 5')
	print out[1]
	if out[0]!=0:
		print >> sys.stderr,'ERROR: predict_variants.py not working'
		sys.exit(1)
	return



if __name__ == '__main__':
	if len(sys.argv)<3:
		print 'python setup.py cmd arch_type [hg] '
		print '- cmd: install, web_install or test'
		print '- arch_type: linux.x86_64, linux.x86_64.v287, macOSX.x86_64, etc'
		print '- hg: all, hg19, hg38'
		sys.exit(0)
	web=False
	opt=sys.argv[1]
	if opt=='web_install':
		web=True
		opt='install'
	if opt=='install':
		arch_type=sys.argv[2]
		hg='all'
		if len(sys.argv)>3: hg=sys.argv[3]
		setup(arch_type,hg,web)
	
	elif opt=='test':
		test()

	else:
		print 'python setup.py cmd arch_type [hg] '
		print '- cmd: install or test'
		print '- arch_type: linux.x86_64, linux.x86_64.v287, macOSX.x86_64, etc'
		print '- hg: all, hg19, hg38' 



