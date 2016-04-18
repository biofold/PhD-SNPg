#!/usr/bin/env python
import os,sys,pickle
from commands import getstatusoutput
from sklearn.externals import joblib


def global_vars():
	global tool_dir, prog_dir, ucsc_dir, ucsc_exe, prog_dat, verbose, hg19, hg38
	verbose = False
	tool_dir = os.path.dirname(os.path.abspath(__file__))
	prog_dir = '/'.join(tool_dir.split('/')[:-1])
	prog_dat = prog_dir + '/data/model'
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
	nuc=''
	seq=''
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
		return seq,nuc
	if out2[0]!=0:	
		print >> sys.stderr,'ERROR: Sequence fetch -', out2[1]
		return seq,nuc
	seq1=check_seq(out1[1],win+1,True)
	seq2=check_seq(out2[1],win+1,False)
	nuc=seq1[-1]
	seq=seq1[:-1]+seq2
	return nuc,seq


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
	nuc=''
	seq=''
	seq_input=[]
	cons_input1=[]
	cons_input2=[]
	ipos=ipos-1
	nuc,seq=get_sequence(ichr,ipos,ucsc_exe,ucsc_dbs,win,dbfasta,fprog)
	if nuc=='' or seq=='': return nuc,seq,seq_input,cons_input1,cons_input2
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
	y_pred,y_fdrs,c_pred=prediction(X,model)
	if y_pred==[]:
		print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
		print '\t'.join([str(i) for i in [ichr,ipos,wt+','+nw] ])+'\tNA\tNA\tNA\tNA\tNA\tNA'
	else:
		print "#CHROM\tPOS\tREF\tALT\tPREDICTION\tSCORE\tFDR\t1-NPV\tPhyloP100\tAvgPhyloP100"
		pp100=cons_input2[win]
		avgpp100=sum(cons_input2)/float(len(cons_input2))
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,c_pred[0],'%.3f' %y_pred[0],'%.3f' %y_fdrs[0][0],'%.3f' %y_fdrs[0][1],pp100,avgpp100])
	return


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
	print "#CHROM\tPOS\tREF\tALT\tPREDICTION\tSCORE\tFDR\t1-NPV\tPhyloP100\tAvgPhyloP100"
	for line in f:	
		v=line.rstrip().split(s)
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		(ichr,pos,wt,nw)=v[:4]
		if len(wt)>1 and len(nw)>1:
			print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'NA','NA','NA','NA','NA','NA'])
			continue
		ipos=int(pos)
		(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpp1,dbpp2,fprog,cprog)
		if seq=='': print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos		
		if seq_input==[]: print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
		if cons_input1==[] or cons_input2==[]: print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
		X=[seq_input + cons_input1+ cons_input2 ]
		y_pred,y_fdrs,c_pred=prediction(X,model)
		if y_pred==[]:
			print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		pp100=cons_input2[win]
		avgpp100=sum(cons_input2)/float(len(cons_input2))
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,c_pred[0],'%.3f' %y_pred[0],'%.3f' %y_fdrs[0][0],'%.3f' %y_fdrs[0][1],pp100,avgpp100])	
	return 


def prediction(X,model,fdr_file='fdr_mean.pkl'):
	y_pred=[]
	y_fdrs=[]
	c_pred=[]
	try:
        	y_pred = model.predict_proba(X)[:, 1]
		c_pred = ['Pathogenic' if i>0.5 else 'Benign' for i in y_pred ]
		if os.path.isfile(prog_dat+'/'+fdr_file):
			fdr_dic=pickle.load(open(prog_dat+'/'+fdr_file))
			y_fdrs=[fdr_dic[round(i,3)] for i in y_pred]
		else:
			y_fdrs=['NA' for i in y_pred]
	except:
		print >> sys.stderr,'WARNING: Prediction errorr check input and scoring models.'
	return y_pred,y_fdrs,c_pred


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
		ochr=ichr
		if ichr.find('chr')==-1: ochr='chr'+ichr
		if ochr=='chrMT': ochr='chrM'
		ipos=int(sys.argv[2])
		wt=sys.argv[3]
		nw=sys.argv[4]
		if len(wt)>1 or len(nw)>1 or 'ACGT'.find(wt)==-1 or 'ACGT'.find(nw)==-1 or wt==nw:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nuleotide',wt,nw
			sys.exit()
		if modfile=='':
			(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,win,fasta,dbpp1,dbpp2)
			if seq_input!=[] and cons_input1!=[] and cons_input2!=[]: 
				chr_data='\t'.join([ichr,str(ipos),wt+','+nw,seq])
				seq_data='\t'.join([str(i) for i in seq_input])
				cons_data1='\t'.join([str(i) for i in cons_input1])
				cons_data2='\t'.join([str(i) for i in cons_input2])
				print '%s\t%s\t%s\t%s' %(chr_data,seq_data,cons_data1,cons_data2)
			else:
				print >> sys.stderr, 'ERROR: Variants',ichr,ipos,wt+','+nw
		else:
			make_prediction(ochr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,win,fasta,dbpp1,dbpp2)
				
