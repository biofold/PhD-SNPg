#!/usr/bin/env python
import os, sys, pickle
from commands import getstatusoutput


def global_vars():
	global tool_dir, prog_dir, ucsc_dir, ucsc_exe, ucsc_web, prog_dat, verbose, hg19, hg38, biofold
	verbose = False
	tool_dir = os.path.dirname(os.path.abspath(__file__))
	prog_dir = '/'.join(tool_dir.split('/')[:-1])
	prog_dat = prog_dir + '/data/model'
	ucsc_dir = prog_dir+'/ucsc'
	ucsc_exe = ucsc_dir+'/exe'
	ucsc = 'http://hgdownload.cse.ucsc.edu/goldenPath'
        biofold = 'http://snps.biofold.org/PhD-SNPg/ucsc'
        ucsc_web = {'hg19.2bit':ucsc+'/hg19/bigZips','hg38.2bit':ucsc+'/hg38/bigZips',\
		'hg19.phyloP46way.primate.bw':biofold+'/hg19','hg19.phyloP46way.placental.bw':biofold+'/hg19','hg19.100way.phyloP100way.bw':ucsc+'/hg19/phyloP100way',\
		'hg38.phyloP7way.bw':ucsc+'/hg38/phyloP7way','hg38.phyloP20way.bw':ucsc+'/hg38/phyloP20way','hg38.phyloP100way.bw':ucsc+'/hg38/phyloP100way'}
	sys.path.insert(0,tool_dir)
	hg19={}
	hg19['fasta']='hg19.2bit'
	hg19['phylop']=['hg19.phyloP46way.primate.bw','hg19.phyloP46way.placental.bw', \
			'hg19.100way.phyloP100way.bw']
	hg19['coding']='hg19_coding.bed'
	hg38={}
	hg38['fasta']='hg38.2bit'
	hg38['phylop']=['hg38.phyloP7way.bw','hg38.phyloP20way.bw','hg38.phyloP100way.bw']
	hg38['phastc']=['hg38.phastCons7way.bw', 'hg38.phastCons20way.bw', 'hg38.phastCons100way.bw']
	hg38['coding']='hg38_coding.bed'
	return


def parse_variants(ichr,pos,wt,nw,ucsc_exe,ucsc_dbs,web=False,dbfasta='hg38.2bit',fprog='twoBitToFa',wseq=100):
	ipos=pos-1
	n_pos=pos
	n_wt=''
	n_nw=''
	if wt=='-': wt=''
	if nw=='-': nw=''
	db_file=ucsc_dbs+'/'+dbfasta
	if web: db_file=ucsc_web[dbfasta]+'/'+dbfasta
	cmd=ucsc_exe+'/'+fprog+' '+db_file+' stdout -seq='+ichr+' -start='+str(ipos)+' -end='+str(ipos+wseq)+' | grep -v "^>" | tr -d "\n" '
	out=getstatusoutput(cmd)
        if out[0]!=0:
                print >> sys.stderr,'ERROR: Sequence fetch -', out[1]
                return n_wt,n_nw,n_pos
	wt_seq=out[1].upper()
	if wt_seq.find(wt)!=0:
		print >> sys.stderr,'ERROR: Sequence fetch -', out[1]
		return n_wt,n_nw,n_pos
	nw_seq=nw+wt_seq[len(wt):]
	lwt=len(wt_seq)
	lnw=len(nw_seq)
	#print wt,nw,pos,wt_seq,nw_seq
	n=min(lnw,lwt)
	v_nuc=['ACGTN'.find(i) for i in nw_seq]
	if -1 in v_nuc:
		print >> sys.stderr,'ERROR: Sequence fetch -', out[1]
		return n_wt,n_nw,n_pos
	for i in range(n):
		#print n_pos,wt_seq[i],nw_seq[i]
		if wt_seq[i]!=nw_seq[i]:
			n_wt=wt_seq[i]
			n_nw=nw_seq[i]
			# Stop at the first nuclotide 
			break
		else:
			n_pos=n_pos+1
	return n_wt,n_nw,n_pos


def check_seq(seq,slen,right=True):
	n=len(seq)
	if n<slen and right: seq=seq+(slen-n)*'X'
	if n<slen and not right: seq=(slen-n)*'X'+seq
	return seq.upper()


def get_sequence(ichr,ipos,ucsc_exe,ucsc_dbs,web=False,win=2,dbname='hg38.2bit',prog='twoBitToFa'):
	nuc=''
	seq=''
	s=max(0,ipos-win)
	e=ipos+win+1
	db_file=ucsc_dbs+'/'+dbname
	if web: db_file=ucsc_web[dbname]+'/'+dbname
	cmd1=ucsc_exe+'/'+prog+' '+db_file+' stdout -seq='+ichr+' -start='+str(s)+' -end='+str(ipos+1)+' | grep -v "^>" | tr -d "\n" '
	cmd2=ucsc_exe+'/'+prog+' '+db_file+' stdout -seq='+ichr+' -start='+str(ipos)+' -end='+str(e)+' | grep -v "^>" | tr -d "\n" '
	if verbose: print >> sys.stderr, cmd1
	if verbose: print >> sys.stderr, cmd2
	out1=getstatusoutput(cmd1)
	out2=getstatusoutput(cmd2)
	if out1[0]!=0:  
		print >> sys.stderr,'ERROR: Sequence fetch -', out1[1]
		sys.exit(1)
		#return seq,nuc
	if out2[0]!=0:	
		print >> sys.stderr,'ERROR: Sequence fetch -', out2[1]
		sys.exit(1)
		#return seq,nuc
	seq1=check_seq(out1[1],win+1,True)
	seq2=check_seq(out2[1],win+1,False)
	nuc=seq1[-1]
	seq=seq1[:-1]+seq2
	return nuc,seq


def get_seqinput(seq,wt,nw,pos):
	n='ACGT'
	vs=[]
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


def check_conservation(cons,win=2):
	if cons.count('NA')==2*win+1:
        	print >> sys.stderr, 'ERROR: Profile not available in region.'
		return []
	else:
		cons=[0 if x == "NA" else x for x in cons]
	return cons


def get_conservation(ichr,ipos,ucsc_exe,ucsc_dbs,web=False,win=2,dbname='hg38.phyloP100way.bw',prog='bigWigToBedGraph'):
	s=max(0,ipos-win)
        e=ipos+win+1
	v=[]
	for i in range(s,e):
		db_file=ucsc_dbs+'/'+dbname
		if web: db_file=ucsc_web[dbname]+'/'+dbname
		cmd=ucsc_exe+'/'+prog+' '+db_file+' stdout -chrom='+ichr+' -start='+str(i)+' -end='+str(i+1)
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


def get_fconservation(ichr,ipos,ucsc_exe,ucsc_dbs,web=False,win=2,dbname='hg38.phyloP100way.bw',prog='bigWigToBedGraph'):
	s=max(0,ipos-win)
	e=ipos+win+1
	v=["NA" for i in range(2*win+1)]
	db_file=ucsc_dbs+'/'+dbname
	if web: db_file=ucsc_web[dbname]+'/'+dbname
	cmd=ucsc_exe+'/'+prog+' '+db_file+' stdout -chrom='+ichr+' -start='+str(s)+' -end='+str(e)
	if verbose: print >> sys.stderr, cmd
	out=getstatusoutput(cmd)
	if out[0]!=0:
		print >> sys.stderr,'ERROR: Conservation fetch -', out[1]
		sys.exit(1)
	else:
		for line in out[1].split('\n'):
			vd=line.split()
			if len(vd)<3: continue
			r1=int(vd[1])
			r2=int(vd[2])	
			for i in range(r1,r2):
				v[i-s]=float(vd[-1])
	cons=check_conservation(v,win)
        return cons


def get_phdsnp_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	# for 0 starting genome
	nuc=''
	seq=''
	seq_input=[]
	cons_input1=[]
	cons_input2=[]
	#Set genome starting position to 0
	ipos=ipos-1
	nuc,seq=get_sequence(ichr,ipos,ucsc_exe,ucsc_dbs,web,win,dbfasta,fprog)
	if nuc=='' or seq=='': return nuc,seq,seq_input,cons_input1,cons_input2
        seq_input=get_seqinput(seq,wt,nw,win)
        cons_input1=get_fconservation(ichr,ipos,ucsc_exe,ucsc_dbs,web,win,dbpps[0],cprog)
        cons_input2=get_fconservation(ichr,ipos,ucsc_exe,ucsc_dbs,web,win,dbpps[1],cprog)
	return nuc,seq,seq_input,cons_input1,cons_input2


def get_snv_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	# for 0 starting genome
	nuc=''
	seq=''
	seq_input=[]
	cons_input=[]
	#Set genome starting position to 0
	ipos=ipos-1
	nuc,seq=get_sequence(ichr,ipos,ucsc_exe,ucsc_dbs,web,win,dbfasta,fprog)
	if nuc=='' or seq=='': return nuc,seq,seq_input,cons_input
	seq_input=get_seqinput(seq,wt,nw,win)
	for dbpp in dbpps:
		cons=get_fconservation(ichr,ipos,ucsc_exe,ucsc_dbs,web,win,dbpp,cprog)
		cons_input.append(cons)
	return nuc,seq,seq_input,cons_input


def get_indel_input(ichr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	# for 0 starting genome        
	nuc=''
	seq=''
	seq_input=[]
	cons_input=[]
	r_cod=[]
	#Set genome starting position to 0        
	ipos=ipos-1
	nuc,seq=get_sequence(ichr,ipos,ucsc_exe,ucsc_dbs,web,win,dbfasta,fprog)
	if nuc=='' or seq=='': return nuc,seq,seq_input,cons_input,r_cod	
	seq_input=get_seqinput(seq,wt,nw,win)
	for dbpp in dbpps:
		cons=get_fconservation(ichr,ipos,ucsc_exe,ucsc_dbs,web,win,dbpp,cprog)
		cons_input.append(cons)
        if pklcod!='': r_cod=get_coding_range(ichr,ipos,ucsc_exe,ucsc_dbs,pklcod)
	return nuc,seq,seq_input,cons_input,r_cod


def get_coding_range(ichr,ipos,ucsc_exe,ucsc_dbs,pklcod):
	r_cod=[]
	if not os.path.isfile(ucsc_dbs+'/'+pklcod):
		print >> sys.stderr,'ERROR: Coding region bed file not found.'
		sys.exit(1)
	cmd=ucsc_exe+'/bed_range.sh '+str(ichr)+' '+str(ipos+1)+' '+ucsc_dbs+'/'+pklcod
	out=getstatusoutput(cmd)
	if out[0]==0:
		if out[1]!="NOT-CODING": r_cod=[int(i) for i in out[1].split('-')]
	else:
		print >> sys.stderr,'WARNING: Coding region check failed. ',out[1]
	return r_cod	


def make_prediction(ichr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	lwt=len(wt)
	lnw=len(nw)
	if wt=='-':
		lwt=1
		lnw+=1
	if nw=='-':
		lwt+=1
		lnw=1
	n_wt,n_nw,n_pos=parse_variants(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web,dbfasta,fprog)
	if n_wt=='' or n_nw=='':
		print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		sys.exit()
	if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
		print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
		sys.exit()
	if pklcod=='':
		r_cod=[]
		(nuc,seq,seq_input,cons_input)=get_snv_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,fprog,cprog)
	else:
		(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
	if seq=='': 
		print >> sys.stderr, 'ERROR: Sequence not found for position',ichr,ipos
		sys.exit(1)
	if seq_input==[]: 
		print >> sys.stderr, 'ERROR: Incorrect nucleotide in position',ichr,ipos
		sys.exit(1)
	if cons_input.count([])>0: 
		print >> sys.stderr, 'ERROR: Incorrect conservation data in position',ichr,ipos
		sys.exit(1)
	if cons_input!=[]:
		cons_input1=cons_input[0]
		cons_input2=cons_input[1]
	else:
		print >> sys.stderr, 'ERROR: Incorrect conservation data in position',ichr,ipos
		sys.exit(1)
	try:
		model=joblib.load(modfile)
	except:
		print >> sys.stderr,'ERROR: Program not able to load modfile. Please check that you have installed a compatible version joblib.'
		sys.exit(1)
	if pklcod=='':
		X=[seq_input + cons_input1+ cons_input2 ]
		y_pred,y_fdrs,c_pred=prediction(X,model)
		v_fdr=[y_fdrs[0][0],y_fdrs[0][1]]
	else:
		p_cod=0
		if r_cod!=[]: p_cod=1
		X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
		y_pred,y_fdrs,c_pred=prediction(X,model)
		v_fdr=[y_fdrs[0][2],y_fdrs[0][3]]
	if y_pred==[]:
		print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
		print '\t'.join([str(i) for i in [ichr,ipos,wt+','+nw] ])+'\tNA\tNA\tNA\tNA\tNA'
	else:
		print "#CHROM\tPOS\tREF\tALT\tPREDICTION\tSCORE\tFDR\tPhyloP100\tAvgPhyloP100"
		pp100=cons_input2[win]
		avgpp100=sum(cons_input2)/float(len(cons_input2))
		if c_pred[0] == "Pathogenic": d_fdr=v_fdr[0]
		if c_pred[0] == "Benign": d_fdr=v_fdr[1]
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,c_pred[0],'%.3f' %y_pred[0],'%.3f' %d_fdr,'%.3f' %pp100,'%.3f' %avgpp100])
	return


def get_file_input(namefile,ucsc_exe,ucsc_dbs,web=False,win=2,s='\t',dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	vlines=[]
	f=open(namefile)
	c=1
	for line in f:
		v=line.rstrip().split(s)
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		(ichr,pos,wt,nw)=v[:4]
		ipos=int(pos)
		lwt=len(wt)
		lnw=len(nw)
		if wt=='-':
			lwt=1
			lnw+=1
		if nw=='-':
			lwt+=1
			lnw=1
		n_wt,n_nw,n_pos=parse_variants(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web,dbfasta,fprog)
		if n_wt=='' or n_nw=='':
			print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
		(nuc,seq,seq_input,cons_input,r_cod)=get_indels_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,pklcod,dbpps,fprog,cprog)
		cons_input1=cons_input[0]
		cons_input2=cons_input[1]
		if seq=='': 
			print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos
		if seq_input==[]:
			print >> sys.stderr, 'WARNING: Incorrect nucleotide in line '+str(c)+'. Genome location:',ichr,pos
		if cons_input1==[] or cons_input2==[]:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
		p_cod=0
		if r_cod!=[]: p_cod=1
		vlines.append([v,seq,seq_input,cons_input1,cons_input2,[lwt,lnw,p_cod]])
		c=c+1
	return vlines
		

def make_file_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web=False,win=2,s='\t',dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	try:
		model1=joblib.load(modfile[0])
		model2=joblib.load(modfile[1])
	except:
		print >> sys.stderr,'ERROR: Program not able to load modfile. Please check that you have installed a compatible version joblib.'
		sys.exit(1)
	f=open(namefile)
	c=1
	print "#CHROM\tPOS\tREF\tALT\tPREDICTION\tSCORE\tFDR\tPhyloP100\tAvgPhyloP100"
	for line in f:	
		v=line.rstrip().split(s)
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		(ichr,pos,wt,nw)=v[:4]
		ipos=int(pos)
		lwt=len(wt)
		lnw=len(nw)
		if wt=='-':
			lwt=1
			lnw+=1
		if nw=='-':
			lwt+=1
			lnw=1
		n_wt,n_nw,n_pos=parse_variants(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web,dbfasta,fprog)
		if n_wt=='' or n_nw=='':
			print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw		
		if wt==nw or nw.find(',')>-1:
			print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'NA','NA','NA','NA','NA'])
			continue
		if len(wt)==1 and len(nw)==1:
			r_cod=[]
			(nuc,seq,seq_input,cons_input)=get_snv_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,fprog,cprog)
		else: 
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
		if seq=='': print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos		
		if seq_input==[]: print >> sys.stderr, 'WARNING: Incorrect nucleotide in line '+str(c)+'. Genome location:',ichr,pos
		if cons_input!=[]:
			cons_input1=cons_input[0]
			cons_input2=cons_input[1]
		else:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
			print line+'\tNA\tNA\tNA\tNA\tNA'
			continue
		if cons_input1==[] or cons_input2==[]: print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
		if len(wt)==1 and len(nw)==1:
			X=[seq_input + cons_input1+ cons_input2 ]
			y_pred,y_fdrs,c_pred=prediction(X,model1)
			v_fdr=[y_fdrs[0][0],y_fdrs[0][1]]
		else:
			p_cod=0
			if r_cod!=[]: p_cod=1
			X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
			y_pred,y_fdrs,c_pred=prediction(X,model2)
			v_fdr=[y_fdrs[0][2],y_fdrs[0][3]]
		if y_pred==[]:
			print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
			print line+'\tNA\tNA\tNA\tNA\tNA'
			continue
		pp100=cons_input2[win]
		avgpp100=sum(cons_input2)/float(len(cons_input2))
		if c_pred[0] == "Pathogenic": d_fdr=v_fdr[0]
		if c_pred[0] == "Benign": d_fdr=v_fdr[1]
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,c_pred[0],'%.3f' %y_pred[0],'%.3f' %d_fdr,pp100,avgpp100])	
	return 


def prediction(X,model,fdr_file='fdr_mean.pkl'):
	y_pred=[]
	y_fdrs=[]
	c_pred=[]
	nf=model.n_features
	if len(X[0]) != nf:
		print >> sys.stderr,'ERROR: Model expecting',nf,'features. Input vector with',len(X[0]),'features. Check the modfile or the window size.'
		sys.exit(1)
	try:
		y_pred = model.predict_proba(X)[:, 1]
		c_pred = ['Pathogenic' if i>0.5 else 'Benign' for i in y_pred ]
		if os.path.isfile(prog_dat+'/'+fdr_file):
			fdr_dic=pickle.load(open(prog_dat+'/'+fdr_file))
			y_fdrs=[fdr_dic[round(i,3)] for i in y_pred]
		else:
			y_fdrs=['NA' for i in y_pred]
	except:
		print >> sys.stderr,'WARNING: Prediction error check input and scoring models.'
	return y_pred,y_fdrs,c_pred


def get_options():
	import optparse
	desc = 'Script for running ContrastRank scoring pipeline'
	parser = optparse.OptionParser("usage: [-h] [-t var_type] [-o outfile]", description=desc)
	parser.add_option('-o','--output', action='store', type='string', dest='outfile', help='Output file')
	parser.add_option('-m','--mod-file', action='store', type='string', dest='mfile', help='Model file')
	parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='hg38', help='Genome version')
	parser.add_option('-w','--win', action='store', type='int', dest='win', help='Window length')
	parser.add_option('--web', action='store_true', dest='web', default=False, help='Use UCSC web files')
	
	(options, args) = parser.parse_args()
	outfile=''
	modfile=''
	hg='hg38'
	web=False
	win=2
	if options.mfile: 
		modfile=options.mfile
		if not os.path.isfile(modfile):
			print >> sys.stderr,'ERROR: Modfile',modfile,'not found.'
			sys.exit(1)
	if options.hg.lower()=='hg19': hg='hg19'
	if options.web: web=True
	if options.win>0: win=options.win
	if hg=='hg19':
		fasta=hg19['fasta']
		dbpps=[hg19['phylop'][0],hg19['phylop'][2]]
		pklcod=cod=hg19['coding']
	else:
		fasta=hg38['fasta']
		dbpps=[hg38['phylop'][0],hg38['phylop'][2]]
		#dbpps=hg38['phylop']+hg38['phastc']
		pklcod=hg38['coding']
	opts=(outfile,modfile,web,win,hg,fasta,dbpps,pklcod)
	return args,opts



if __name__ == '__main__':
	global_vars()
	from sklearn.externals import joblib
	args,opts=get_options()
	(outfile,modfile,web,win,hg,fasta,dbpps,pklcod)=opts
	ucsc_dbs=ucsc_dir+'/'+hg
	if len(args)>1:
		ichr=sys.argv[1]
		pos=sys.argv[2]
		wt=sys.argv[3]
		nw=sys.argv[4]
		ochr=ichr
		if ichr.find('chr')==-1: ochr='chr'+ichr
		if ochr=='chrMT': ochr='chrM'
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect input position',chr,pos
			sys.exit(1)
		if wt==nw:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',ichr,ipos,wt,nw
			sys.exit()
		if modfile=='':
			lwt=len(wt)
			lnw=len(nw)
			if wt=='-':
				lwt=1
				lnw+=1
			if nw=='-':
				lwt+=1
				lnw=1
			n_wt,n_nw,n_pos=parse_variants(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web,fasta)
			if n_wt=='' or n_nw=='':
				print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
				sys.exit()
			if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
				print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
				sys.exit()	
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(ochr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,fasta,dbpps,pklcod)
			if seq_input!=[] and cons_input.count([])==0: 
				chr_data='\t'.join([ichr,str(ipos),str(n_pos),wt+','+nw,n_wt+str(n_pos)+n_nw,seq])
				seq_data='\t'.join([str(i) for i in seq_input])
				cons_data=''
				for cons in cons_input:
					cons_data=cons_data+'\t'.join([str(i) for i in cons])+'\t'
				cons_data=cons_data[:-1]
				p_cod=0
				if r_cod!=[]: p_cod=1
				m_len='%d\t%d\t%d' %(lwt,lnw,p_cod)
				print '%s\t%s\t%s\t%s' %(chr_data,seq_data,cons_data,m_len)
			else:
				print >> sys.stderr, 'ERROR: Variants',ichr,ipos,wt+','+nw
		else:
			if len(wt)==1 and len(nw)==1: pklcod=''	
			make_prediction(ochr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,web,win,fasta,dbpps,pklcod)
	else:
		print "python score_variants.py chromosome position ref_nuc alt_nuc -g hg_version"
				
