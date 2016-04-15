#!/usr/bin/env python
import os,sys
from commands import getstatusoutput


def global_vars():
	global prog_dir, ucsc_dir, verbose
	verbose = False
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	ucsc_dir = prog_dir+'/ucsc'
	return


def check_seq(seq,slen,right=True):
	n=len(seq)
	if n<slen and right: seq=seq+(slen-n)*'X'
	if n<slen and not right: seq=(slen-n)*'X'+seq
	return seq.upper()


def get_sequence(ichr,ipos,ucsc_dir,win=3,dbname='hg38.2bit',prog='twoBitToFa'):
	s=max(0,ipos-win)
	e=ipos+win+1
	cmd1=ucsc_dir+'/'+prog+' '+ucsc_dir+'/'+dbname+' stdout -seq='+ichr+' -start='+str(s)+' -end='+str(ipos+1)+' | grep -v "^>" '
	cmd2=ucsc_dir+'/'+prog+' '+ucsc_dir+'/'+dbname+' stdout -seq='+ichr+' -start='+str(ipos)+' -end='+str(e)+' | grep -v "^>" '
	if verbose: print cmd1
	if verbose: cmd2
	out1=getstatusoutput(cmd1)
	out2=getstatusoutput(cmd2)
	if out1[0]!=0:  
		print >> sys.stderr,'ERROR: Sequence fetch -', out1[1]
		sys.exit(1)
	if out2[0]!=0:	
		print >> sys.stderr,'ERROR: Sequence fetch -', out2[1]
		sys.exit(1)
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
				sys.exit(1)
			else:
				v[n.find(seq[i])]=-100
				v[n.find(nw)]=100 
		vs.append(v)
	return vs


def check_conservation(cons,win=3):
	if cons.count('NA')==2*win+1:
        	print >> sys.stderr, 'ERROR: Profile not available in region.'
		sys.exit(1)
	else:
		cons=[0 if x == "NA" else x for x in cons]
	return cons


def get_conservation(ichr,ipos,ucsc_dir,win=3,dbname='hg38.phyloP100way.bw',prog='bigWigToBedGraph'):
	s=max(0,ipos-win)
        e=ipos+win+1
	v=[]
	for i in range(s,e):
		cmd=ucsc_dir+'/'+prog+' '+ucsc_dir+'/'+dbname+' stdout -chrom='+ichr+' -start='+str(i)+' -end='+str(i+1)
		if verbose: print cmd
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


def get_phdsnp_input(ichr,ipos,ucsc_dir,win=3,dbfasta='hg38.2bit',dbpp1='hg38.phyloP7way.bw',dbpp2='hg38.phyloP100way.bw',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	nuc,seq=get_sequence(ichr,ipos,ucsc_dir,win,dbfasta,fprog)
        seq_input=get_seqinput(seq,wt+str(win+1)+nw)
        cons_input1=get_conservation(ichr,ipos,ucsc_dir,win,dbpp1,cprog)
        cons_input2=get_conservation(ichr,ipos,ucsc_dir,win,dbpp2,cprog)
	return nuc,seq,seq_input,cons_input1,cons_input2



def get_options():
	import optparse
	desc = 'Script for running ContrastRank scoring pipeline'
	parser = optparse.OptionParser("usage: [-h] [-t var_type] [-o outfile]", description=desc)
        parser.add_option('-t', '--threshold',action='store',type='float',dest='th', default=0.005, help='1000Genomes allele frequency filter')
	parser.add_option('-o','--output', action='store', type='string', dest='outfile', help='Output file')
	parser.add_option('-q', action='store_true', dest='queue', help='Use queue system')
	(options, args) = parser.parse_args()
	return args,options




if __name__ == '__main__':
	global_vars()
	args,opts=get_options()
	if len(args)>1:
		ichr=sys.argv[1]
		ipos=int(sys.argv[2])
		wt=sys.argv[3]
		nw=sys.argv[4]
		win=3
		(nuc,seq,seq_input,cons_input1,cons_input2)=get_phdsnp_input(ichr,ipos,ucsc_dir)
		#nuc,seq=get_sequence(ichr,ipos,ucsc_dir)
		#seq_input=get_seqinput(seq,wt+str(win+1)+nw)
		#cons_input7=get_conservation(ichr,ipos,ucsc_dir,3,'hg38.phyloP7way.bw')
		#cons_input100=get_conservation(ichr,ipos,ucsc_dir,3,'hg38.phyloP100way.bw')
		print nuc, seq, seq_input, cons_input1, cons_input2
