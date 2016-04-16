#!/usr/bin/env python
import os, sys, subprocess, __builtin__
from commands import getstatusoutput
from sklearn.externals import joblib

def global_vars():
	global tool_dir, prog_dir, ucsc_dir, ucsc_exe, verbose, hg19, hg38, prog_cat
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	tool_dir = prog_dir+'/tools'
	ucsc_dir = prog_dir+'/ucsc'
	ucsc_exe = ucsc_dir+'/exe'
	prog_cat = 'zcat'
	sys.path.append(tool_dir)
	hg19={}
	hg19['fasta']='hg19.2bit'
	hg19['phylop']=['hg19.phyloP46way.primate.bw','hg19.phyloP46way.placental.bw','hg19.100way.phyloP100way.bw']
	hg38={}
	hg38['fasta']='hg38.2bit'
	hg38['phylop']=['hg38.phyloP7way.bw','hg38.phyloP20way.bw','hg38.phyloP100way.bw']
	return



def make_prediction(ichr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,win=3,dbfasta='hg38.2bit',dbpp1='hg38.phyloP7way.bw',dbpp2='hg38.phyloP100way.bw',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpp1,dbpp2,fprog,cprog)
	if seq=='': 
		print >> sys.stderr, 'ERROR: Sequence not found for position',ichr,ipos
		sys.exit(1)
	if seq_input==[]: 
		print >> sys.stderr, 'ERROR: Incorrect nucleotide in position',ichr,ipos
		sys.exit(1)
	#if cons_input1==[] or cons_input2==[]: 
	#Check only P100
	if cons_input2==[]:
		print >> sys.stderr, 'ERROR: Incorrect conservation data in position',ichr,ipos
		sys.exit(1)
	if cons_input1==[]: cons_input1=[0.0 for i in range(2*win+1)]
	model=joblib.load(modfile)
	X=[seq_input + cons_input1+ cons_input2 ]
	y_pred=prediction(X,model)
	print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])



def make_file_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win=3,s='\t',dbfasta='hg38.2bit',dbpp1='hg38.phyloP7way.bw',dbpp2='hg38.phyloP100way.bw',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	model=joblib.load(modfile)
	proc = subprocess.Popen([prog_cat,'-f',namefile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()        
	c=1
	for line in stdout.split('\n'):
		if line == '': continue 
		if line[0]=='#':
			if vcf and line.find('#CHROM')==0: line=line+'\tPREDICTION' 
			print line
			continue 	
		v=line.rstrip().split(s)
		if len(v)<4:
			print line
                        continue
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		if vcf:
			if fpass and len(v)>6 and v[6]!='PASS':
				print line+'\tNA'
				continue
			(ichr,pos,rs,wt,nw)=tuple(v[:5])
		else:
			(ichr,pos,wt,nw)=tuple(v[:4])
		if len(wt)>1 or len(nw)>1 or 'ACGT'.find(wt)==-1 or 'ACGT'.find(nw)==-1 or wt==nw:
			print line+'\tNA'
			#print '\t'.join(str(i) for i in [ichr,pos,wt,nw,'NA'])
			continue
		nchr=ichr
		if nchr.find('chr')==-1: nchr='chr'+ichr
		ipos=int(pos)
		(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(nchr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpp1,dbpp2,fprog,cprog)
		if seq=='': 
			print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos
			print line+'\tNA'
			continue
		if seq_input==[]: 
			print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
			print line+'\tNA'
			continue
		#if cons_input1==[] or cons_input2==[]:
		#Check only P100
		if cons_input2==[]:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
			print line+'\tNA'
			continue
		if cons_input1==[]: cons_input1=[0.0 for i in range(2*win+1)]
		X=[ seq_input + cons_input1+ cons_input2 ]
		y_pred=prediction(X,model)
		print line+'\t'+'%.4f' %y_pred[0]
		#print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])	
	return 


def prediction(X,model):
        y_pred = model.predict_proba(X)[:, 1]
        return y_pred		


def get_options():
	global hg, coord, vcf, fpass, win
	import optparse
	desc = 'Script for scoring single nucleotide variants'
	parser = optparse.OptionParser("usage: %prog variant_file", description=desc)
	parser.add_option('-m','--mod-file', action='store', type='string', dest='mfile', help='Model file')
	parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='hg38', help='Genome version')
	parser.add_option('-v','--verbose', action='store_true', dest='ver', default=False, help='Verbose mode')
	parser.add_option('-c','--coordinate', action='store_true', dest='coord', default=False, help='Coordinate input')
	parser.add_option('--vcf', action='store_true', dest='vcf', default=False, help='VCF file input')
	parser.add_option('--pass', action='store_true', dest='fpass', default=False, help='Predict only PASS variants. Check column 7 in vcf file')
	(options, args) = parser.parse_args()
	outfile = ''
	modfile = prog_dir + '/data/model_w3_p7_100.pkl'
	hg='hg38'
	coord=False
	vcf=False
	fpass=False
	win=3
	if options.mfile: modfile=options.mfile
	if options.hg.lower()=='hg19': hg='hg19'
	if options.coord: coord = True
	if options.fpass: fpass=True
	if options.vcf: vcf = True
	__builtin__.verbose=False
	if options.ver: __builtin__.verbose=True
	if hg=='hg19':
		fasta=hg19['fasta']
		dbpp1=hg19['phylop'][0]
		dbpp2=hg19['phylop'][2]
	else:
		fasta=hg38['fasta']
		dbpp1=hg38['phylop'][0]
		dbpp2=hg38['phylop'][2]
	if not os.path.isfile(modfile):
                print >> sys.stderr,'ERROR: Data model file not found'
		sys.exit(1)
	opts=(outfile,modfile,fasta,dbpp1,dbpp2)
	return args,opts



if __name__ == '__main__':
	global_vars()
	from score_variant import get_phdsnp_input	
	args,opts=get_options()
	(outfile,modfile,fasta,dbpp1,dbpp2)=opts
	ucsc_dbs=ucsc_dir+'/'+hg
	if len(args)>0:
		if coord: 
			(ichr,ipos,wt,nw)=sys.argv[1].split(',')[:4]
			ipos=int(ipos)
			if len(wt)>1 or len(nw)>1 or 'ACGT'.find(wt)==-1 or 'ACGT'.find(nw)==-1 or wt==nw:
				print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nuleotide',wt,nw
				sys.exit(1)
			make_prediction(ichr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,win,fasta,dbpp1,dbpp2)
		else:	
			namefile=args[0]
			if not os.path.isfile(namefile):
				print >> sys.stderr(),'ERROR: Input file not found',namefile
				sys.exit(1)
			make_file_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win,'\t',fasta,dbpp1,dbpp2)
	else:
		print 'predict_variants.py variant_file'
