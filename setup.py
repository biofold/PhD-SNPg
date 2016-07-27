import os, sys 
from commands import getstatusoutput
global ucsc_tool, ucsc_path, biofold_path
ucsc_tool='http://hgdownload.cse.ucsc.edu/admin/exe'
ucsc_path='http://hgdownload.cse.ucsc.edu/goldenPath'
biofold_path='http://snps.biofold.org/PhD-SNPg/ucsc'


def get_ucsc_tools(arch_type):
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	ucsc_dir = prog_dir+'/ucsc'
	twobit='twoBitToFa'
	bwg='bigWigToBedGraph'
	prog_get='wget'
	wtwobit=ucsc_tool+'/'+arch_type+'/'+twobit
	cmd=prog_get+' '+wtwobit+' -O '+ucsc_dir+'/exe/'+twobit+'; chmod a+x '+ucsc_dir+'/exe/'+twobit
	print >> sys.stderr,'   Download twoBitToFa'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	wbwg=ucsc_tool+'/'+arch_type+'/'+bwg
	cmd=prog_get+' '+wbwg+' -O '+ucsc_dir+'/exe/'+bwg+'; chmod a+x '+ucsc_dir+'/exe/'+bwg
	print >> sys.stderr,'   Download bigWigToBedGraph'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	print >> sys.stderr,'   Test tools'
	cmd=ucsc_dir+'/exe/twoBitToFa'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0 and out[0]!=65280:
		print >> sys.stderr,'  ERROR: Incorrect architecture',arch_type
		cmd='rm '+ucsc_dir+'/exe/twoBitToFa'
		getstatusoutput(cmd)
		return out
	cmd=ucsc_dir+'/exe/bigWigToBedGraph'
	print cmd
	out=getstatusoutput(cmd)
	print out[1]
	if out[0]!=0 and out[0]!=65280:
		print >> sys.stderr,'  ERROR: Incorrect architecture',arch_type
		cmd='rm '+ucsc_dir+'/exe/bigWigToBedGraph'
		getstatusoutput(cmd)
		return out
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
			print >> sys.stderr, 'ERROR: Problem in downloading hg19 data.'
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
			print >> sys.stderr, 'ERROR: Problem in downloading hg38 data'
			sys.exit(1)
	print   '   Downloaded UCSC data'


def test(hg='all',web=False):
	if hg=='hg19':
		hgs=['hg19']
	elif hg=='hg38':
		hgs=['hg38']
	else:
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

	if web:	
		if hg=='all' or hg=='hg19':
			ihg='hg19'
			print '\n3) Check web hg19 files'
			hg19s=[ucsc_path+'/'+ihg+'/bigZips/'+ihg+'.2bit', \
				biofold_path+'/'+ihg+'/'+ihg+'.phyloP46way.primate.bw',\
				ucsc_path+'/'+ihg+'/phyloP100way/'+ihg+'.100way.phyloP100way.bw']
			for hgfile in hg19s:
				cmd='curl --head -s '+hgfile
				out=getstatusoutput(cmd)
				if out[1].find('HTTP/1.1 200 OK')==-1:
					print >> sys.stderr,'ERROR: Data file',hgfile,'not available.'
					sys.exit(1)	
				else:
					print 'File',hgfile,'available.'

		if hg=='all' or hg=='hg38':
			ihg='hg38'
			print '\n4) Check web hg38 files'	
			hg38s=[ ucsc_path+'/'+ihg+'/bigZips/'+ihg+'.2bit', \
				ucsc_path+'/'+ihg+'/phyloP7way/'+ihg+'.phyloP7way.bw',\
				ucsc_path+'/'+ihg+'/phyloP100way/'+ihg+'.phyloP100way.bw']
			for hgfile in hg38s:
				cmd='curl --head -s '+hgfile
				out=getstatusoutput(cmd)			
				if out[1].find('HTTP/1.1 200 OK')==-1:
					print >> sys.stderr,'ERROR: Data file',hgfile,'not available.'
					sys.exit(1)
				else:
					print 'File',hgfile,'available.'
	else:
		if hg=='all' or hg=='hg19':
			print '\n3) Check local hg19 files. Please wait, md5sum can take few minutes.'
        		cmd='cd '+ucsc_dir+'/hg19/; md5sum -c hg19.md5'
			print cmd
			out=getstatusoutput(cmd)
			print out[1]
			if out[0]!=0:
				print >> sys.stderr,'WARNING: md5sum failed on hg19'
				print >> sys.stderr,'hg19 predictions will run only in web mode '
		
		if hg=='all' or hg=='hg38':
			print '\n4) Check local hg38 files. Please wait, md5sum can take few minutes.'
			cmd='cd '+ucsc_dir+'/hg38/; md5sum -c hg38.md5'
			print cmd
			out=getstatusoutput(cmd)
			print out[1]
			if out[0]!=0:
				print >> sys.stderr,'WARNING: md5sum failed on hg38'
				print >> sys.stderr,'hg38 predictions will run only in web mode '

	#if hgs==[]:
	#	print >> sys.stderr,'ERROR: UCSC data file not correctly downloaded'
	#	sys.exit()

	print '\n5) Test twoBitToFa command'
	for ihg in hgs:	
		if web:
			twobit=ucsc_path+'/'+ihg+'/bigZips/'+ihg+'.2bit'
		else:
			twobit=ucsc_dir+'/'+ihg+'/'+ihg+'.2bit'
		cmd=ucsc_tool+'/twoBitToFa '+twobit+' stdout -seq=chr1 -start=10008 -end=10010'
		print cmd
		out=getstatusoutput(cmd)
		print out[1]
		if out[0]!=0:
			print >> sys.stderr,'ERROR: twoBitToFa not working with',ihs
			sys.exit(1)

	print '\n6) Test bigWigToBedGraph command'
	for ihg in hgs:
		if web:
			if ihg=='hg19':
				bwg=biofold_path+'/'+ihg+'/'+ihg+'.phyloP46way.primate.bw'
			else:
				bwg=ucsc_path+'/'+ihg+'/phyloP100way/'+ihg+'.phyloP100way.bw'
		else:
			if ihg=='hg19':
				bwg=ucsc_dir+'/'+ihg+'/'+ihg+'.100way.phyloP100way.bw'
			else:
				bwg=ucsc_dir+'/'+ihg+'/'+ihg+'.phyloP100way.bw'
		cmd=ucsc_tool+'/bigWigToBedGraph '+bwg+' stdout -chrom=chr1 -start=100008 -end=100012'
		print cmd
		out=getstatusoutput(cmd)
		print out[1]
		if out[0]!=0:
			print >> sys.stderr,'ERROR: bigWigToBedGraph not working with',ihg
			sys.exit(1)

	print '\n7) Test predict_variants.py'
	for ihg in hgs:
		if web:
			cmd='python predict_variants.py test/test_short_variants_'+ihg+'.tsv -g '+ihg+' --web '
		else:
			cmd='python predict_variants.py test/test_short_variants_'+ihg+'.tsv -g '+ihg
		print cmd
		out=getstatusoutput(cmd)
		print out[1]
		if out[0]!=0:
			print >> sys.stderr,'ERROR: predict_variants.py not working with',ihg
			sys.exit(1)
	return



def get_options():
	import optparse
        desc = 'Script for installing and testing PhD-SNPg.'
        parser = optparse.OptionParser("usage: %prog cmd arch_type [-g hg] [--web]", description=desc)
        parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='all', help='Genome version')
	parser.add_option('--web', action='store_true', dest='web', default=False, help='Use UCSC web files')
	(options, args) = parser.parse_args()
	hg='all'
	web=False
	if options.hg: hg=options.hg.lower()
	if hg!='hg19' and hg!='hg38': hg='all'
	if options.web: web=True
	if len(args)<1:
		print 'python setup.py cmd arch_type [-g hg] [--web]'
		print '  cmd: install or test'
		print '  arch_type: linux.x86_64, linux.x86_64.v287, macOSX.x86_64, etc'
		print '  -g = hg: all, hg19, hg38'
		print '  -web = not download file'
		sys.exit(0)
	return args,hg,web



if __name__ == '__main__':
	args,hg,web=get_options()
	opt=args[0]
	if opt=='install' and len(args)>1:
		arch_type=args[2]
		arch_list=['linux.x86_64.v287','linux.x86_64',\
			'macOSX.x86_64','macOSX.i386','macOSX.ppc']
		if arch_type not in arch_list:
			print >> sys.stderr,'ERROR: Incorrect architecture type.'
			print >> sys.stderr,'Available UCSC precomplied tools are only for '+', '.join(arch_list)
			sys.exit(1)
		setup(arch_type,hg,web)
	
	elif opt=='test':
		test(hg,web)
		
	else:
		print 'python setup.py cmd arch_type [-g hg] [--web]'
		print '  cmd: install or test'
		print '  arch_type: linux.x86_64, linux.x86_64.v287, macOSX.x86_64, etc'
		print '  -g = hg: all, hg19, hg38' 
		print '  -web = not download file'



