#!/usr/bin/env python
import os,sys
from commands import getstatusoutput
from sklearn.externals import joblib


def global_vars():
	global tool_dir, prog_dir, ucsc_dir, ucsc_exe, verbose, hg19, hg38
	verbose = False
	tool_dir = os.path.dirname(os.path.abspath(__file__))
	prog_dir = '/'.join(tool_dir.split('/')[:-1])
	ucsc_dir = prog_dir+'/ucsc'
	ucsc_exe = ucsc_dir+'/exe'
	hg19={}
	hg19['fasta']='hg19.2bit'
	hg19['phylop']=['hg19.phyloP46way.primate.bw','hg19.phyloP46way.placental.bw', \
			'hg19.100way.phyloP100way.bw']
	hg38={}
	hg38['fasta']='hg38.2bit'
	hg38['phylop']=['hg38.phyloP7way.bw','hg38.phyloP20way.bw','hg38.phyloP100way.bw']
	return


def check_seq(seq,slen,right=True):
	n=len(seq)
	if n<slen and right: seq=seq+(slen-n)*'X'
	if n<slen and not right: seq=(slen-n)*'X'+seq
	return seq.upper()


def get_sequence(ichr,ipos,ucsc_exe,ucsc_dbs,win=3,dbname='hg38.2bit',prog='twoBitToFa'):
	s=max(0,ipos-win)
	e=ipos+win+1
	cmd1=ucsc_exe+'/'+prog+' '+ucsc_dbs+'/'+dbname+' stdout -seq='+ichr+' -start='+str(s)+' -end='+str(ipos+1)+' | grep -v "^>" '
	cmd2=ucsc_exe+'/'+prog+' '+ucsc_dbs+'/'+dbname+' stdout -seq='+ichr+' -start='+str(ipos)+' -end='+str(e)+' | grep -v "^>" '
	if verbose: print >> sys.stderr, cmd1
	if verbose: print >> sys.stderr, cmd2
	out1=getstatusoutput(cmd1)
	out2=getstatusoutput(cmd2)
	if out1[0]!=0:  
		print >> sys.stderr,'ERROR: Sequence fetch -', out1[1]
		return ''
	if out2[0]!=0:	
		print >> sys.stderr,'ERROR: Sequence fetch -', out2[1]
		return ''
	seq1=check_seq(out1[1],win+1,True)
	seq2=check_seq(out2[1],win+1,False)
	return seq1[-1],seq1[:-1]+seq2	


def get_seqinput(seq,mut):
	n='ACGT'
	vs=[]
	wt=mut[0]
	pos=int(mut[1:-1])-1
	nw=mut[-1]
	for i in range(len(seq)):
		v=[0,0,0,0,0]	
		if i!=pos:
			v[n.find(seq[i])]=100 
		else:
			if seq[i]!=wt: 
				print >> sys.stderr,'ERROR: Nucleotide not matching -',seq[i],'not',wt
				return []
			else:
				v[n.find(seq[i])]=-100
				v[n.find(nw)]=100 
		vs=vs+v
	return vs


def check_conservation(cons,win=3):
	if cons.count('NA')==2*win+1:
        	print >> sys.stderr, 'ERROR: Profile not available in region.'
		return []
	else:
		cons=[0 if x == "NA" else x for x in cons]
	return cons


def get_conservation(ichr,ipos,ucsc_exe,ucsc_dbs,win=3,dbname='hg38.phyloP100way.bw',prog='bigWigToBedGraph'):
	s=max(0,ipos-win)
        e=ipos+win+1
	v=[]
	for i in range(s,e):
		cmd=ucsc_exe+'/'+prog+' '+ucsc_dbs+'/'+dbname+' stdout -chrom='+ichr+' -start='+str(i)+' -end='+str(i+1)
		if verbose: print >> sys.stderr, cmd
		out=getstatusoutput(cmd)
		if out[0]!=0:
                	print >> sys.stderr,'ERROR: Conservation fetch -', out[1]
                	sys.exit(1)
		if out[1]=='':
			v.append('NA')
		else:
			v.append(float(out[1].split()[-1]))
	cons=check_conservation(v,win)
	return cons


def get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win=3,dbfasta='hg38.2bit',dbpp1='hg38.phyloP7way.bw',dbpp2='hg38.phyloP100way.bw',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	# for 0 starting genome
	ipos=ipos-1
	nuc,seq=get_sequence(ichr,ipos,ucsc_exe,ucsc_dbs,win,dbfasta,fprog)
        seq_input=get_seqinput(seq,wt+str(win+1)+nw)
        cons_input1=get_conservation(ichr,ipos,ucsc_exe,ucsc_dbs,win,dbpp1,cprog)
        cons_input2=get_conservation(ichr,ipos,ucsc_exe,ucsc_dbs,win,dbpp2,cprog)
	return nuc,seq,seq_input,cons_input1,cons_input2


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
		sys.exit(1)
	model=joblib.load(modfile)
	X=[seq_input + cons_input1+ cons_input2 ]
	y_pred=prediction(X,model)
	print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])


def get_file_input(namefile,ucsc_exe,ucsc_dbs,win=3,s='\t',dbfasta='hg38.2bit',dbpp1='hg38.phyloP7way.bw',dbpp2='hg38.phyloP100way.bw',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	vlines=[]
	f=open(namefile)
	c=1
	for line in f:
		v=line.rstrip().split(s)
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		(ichr,pos,wt,nw)=v[:4]
		ipos=int(pos)
		(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpp1,dbpp2,fprog,cprog)
		if seq=='': 
			print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos
		if seq_input==[]:
			print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
		if cons_input1==[] or cons_input2==[]:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
		vlines.append([v,seq,seq_input,cons_input1,cons_input2])
		c=c+1
	return vlines
		

def make_file_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win=3,s='\t',dbfasta='hg38.2bit',dbpp1='hg38.phyloP7way.bw',dbpp2='hg38.phyloP100way.bw',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	model=joblib.load(modfile)
	f=open(namefile)
	c=1
	for line in f:	
		v=line.rstrip().split(s)
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		(ichr,pos,wt,nw)=v[:4]
		if len(wt)>1 and len(nw)>1:
			print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'NA'])
			continue
		ipos=int(pos)
		(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpp1,dbpp2,fprog,cprog)
		if seq=='': print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos		
		if seq_input==[]: print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
		if cons_input1==[] or cons_input2==[]: print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
		X=[seq_input + cons_input1+ cons_input2 ]
		y_pred=prediction(X,model)
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])	
	return 


def prediction(X,model):
        y_pred = model.predict_proba(X)[:, 1]
        return y_pred		


def get_options():
	import optparse
	desc = 'Script for running ContrastRank scoring pipeline'
	parser = optparse.OptionParser("usage: [-h] [-t var_type] [-o outfile]", description=desc)
	parser.add_option('-o','--output', action='store', type='string', dest='outfile', help='Output file')
	parser.add_option('-m','--mod-file', action='store', type='string', dest='mfile', help='Model file')
	parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='hg38', help='Genome version')
	parser.add_option('-w','--win', action='store', type='int', dest='win', help='Window length')
	(options, args) = parser.parse_args()
	outfile=''
	modfile=''
	hg='hg38'
	win=3
	if options.mfile: modfile=options.mfile
	if options.hg.lower()=='hg19': hg='hg19'
	if options.win>0: win=options.win
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
	args,opts=get_options()
	(outfile,modfile,win,hg,fasta,dbpp1,dbpp2)=opts
	ucsc_dbs=ucsc_dir+'/'+hg
	if len(args)>1:
		ichr=sys.argv[1]
		ipos=int(sys.argv[2])
		wt=sys.argv[3]
		nw=sys.argv[4]
		if len(wt)>1 or len(nw)>1 or 'ACGT'.find(wt)==-1 or 'ACGT'.find(nw)==-1 or wt==nw:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nuleotide',wt,nw
			sys.exit()
		if modfile=='':
			(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,fasta,dbpp1,dbpp2)
			print nuc, seq, '\t'.join([str(i) for i in seq_input]), '\t'.join([str(i) for i in cons_input1]), '\t'.join([str(i) for i in cons_input2])
		else:
			make_prediction(ichr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,win,'\t',fasta,dbpp1,dbpp2)
				
