#!/usr/bin/env python
import os, sys, __builtin__
from commands import getstatusoutput
from sklearn.externals import joblib

def global_vars():
	global tool_dir, prog_dir, ucsc_dir, ucsc_exe, verbose, hg19, hg38
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	tool_dir = prog_dir+'/tools'
	ucsc_dir = prog_dir+'/ucsc'
	ucsc_exe = ucsc_dir+'/exe'
	sys.path.append(tool_dir)
	hg19={}
	hg19['fasta']='hg19.2bit'
	hg19['phylop']=['hg19.phyloP46way.primate.bw','hg19.phyloP46way.placental.bw', \
			'hg19.100way.phyloP100way.bw']
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
	if cons_input1==[] or cons_input2==[]: 
		print >> sys.stderr, 'ERROR: Incorrect conservation data in position',ichr,ipos
		sy.exit(1)
	model=joblib.load(modfile)
	X=[seq_input + cons_input1+ cons_input2 ]
	y_pred=prediction(X,model)
	print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])



def make_file_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win=3,s='\t',dbfasta='hg38.2bit',dbpp1='hg38.phyloP7way.bw',dbpp2='hg38.phyloP100way.bw',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	model=joblib.load(modfile)
	f=open(namefile)
	c=1
	for line in f:	
		v=line.rstrip().split(s)
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		(ichr,pos,wt,nw)=v[:4]
		ipos=int(pos)
		(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpp1,dbpp2,fprog,cprog)
		if seq=='': print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos		
		if seq_input==[]: print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
		if cons_input1==[] or cons_input2==[]: print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
		X=[seq_input + cons_input1+ cons_input2 ]
		y_pred=prediction(X,model)
		print ichr,pos,wt,nw,y_pred	
	return 


def prediction(X,model):
        y_pred = model.predict_proba(X)[:, 1]
        return y_pred		


def get_options():
	import optparse
	desc = 'Script for scoring single nucleotide variants'
	parser = optparse.OptionParser("usage: [-h] [-t var_type] [-o outfile]", description=desc)
	parser.add_option('-o','--output', action='store', type='string', dest='outfile', help='Output file')
	parser.add_option('-m','--mod-file', action='store', type='string', dest='mfile', help='Model file')
	parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='hg38', help='Genome version')
	parser.add_option('-w','--win', action='store', type='int', dest='win', help='Window length')
	parser.add_option('-v','--verbose', action='store_true', dest='ver', default=False, help='Verbose mode')
	(options, args) = parser.parse_args()
	outfile = ''
	modfile = prog_dir + '/data/model_w3_p7_100.pkl'
	hg='hg38'
	win=3
	if options.mfile: modfile=options.mfile
	if options.hg.lower()=='hg19': hg='hg19'
	if options.win>0: win=options.win
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
	opts=(outfile,modfile,win,hg,fasta,dbpp1,dbpp2)
	return args,opts



if __name__ == '__main__':
	global_vars()
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	sys.path.append(prog_dir+'/tools')
	from generate_input import get_phdsnp_input	
	args,opts=get_options()
	(outfile,modfile,win,hg,fasta,dbpp1,dbpp2)=opts
	ucsc_dbs=ucsc_dir+'/'+hg
	if len(args)>1:
		ichr=sys.argv[1]
		ipos=int(sys.argv[2])
		wt=sys.argv[3]
		nw=sys.argv[4]
		make_prediction(ichr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,win,fasta,dbpp1,dbpp2)
				
