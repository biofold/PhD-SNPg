#!/usr/bin/env python
import os, sys, subprocess, pickle
import  __builtin__
from commands import getstatusoutput

def global_vars():
	global tool_dir, prog_dir, prog_dat, ucsc_dir, ucsc_exe, ucsc_web, verbose, hg19, hg38, prog_cat, biofold
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	tool_dir = prog_dir+'/tools'
	prog_dat = prog_dir+'/data/model'
	ucsc_dir = prog_dir+'/ucsc'
	ucsc_exe = ucsc_dir+'/exe'
	prog_cat = 'zcat'
	ucsc = 'http://hgdownload.cse.ucsc.edu/goldenPath'
        biofold = 'http://snps.biofold.org/PhD-SNPg/ucsc'
	ucsc_web = {'hg19.2bit':ucsc+'/hg19/bigZips','hg38.2bit':ucsc+'/hg38/bigZips',\
		'hg19.phyloP46way.primate.bw':biofold+'/hg19','hg19.phyloP46way.placental.bw':biofold+'/hg19','hg19.100way.phyloP100way.bw':ucsc+'/hg19/phyloP100way',\
		'hg38.phyloP7way.bw':ucsc+'/hg38/phyloP7way','hg38.phyloP20way.bw':ucsc+'/hg38/phyloP20way','hg38.phyloP100way.bw':ucsc+'/hg38/phyloP100way','hg38.phyloP470way.bw':ucsc+'/hg38/phyloP470way'}
	__builtin__.ucsc_web=ucsc_web
	sys.path.insert(0,tool_dir)
	hg19={}
	hg19['fasta']='hg19.2bit'
	hg19['phylop']=['hg19.phyloP46way.primate.bw','hg19.phyloP46way.placental.bw','hg19.100way.phyloP100way.bw']
	hg19['coding']='hg19_coding.bed'
	hg38={}
	hg38['fasta']='hg38.2bit'
	hg38['phylop']=['hg38.phyloP7way.bw','hg38.phyloP100way.bw','hg38.phyloP470way.bw']
	hg38['coding']='hg38_coding.bed'
	return



def make_prediction(ichr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	lwt=len(wt)
	lnw=len(nw)
	dat='PhyloP100'
	if wt=='-':
		lwt=1
		lnw+=1
	if nw=='-':
		lwt+=1
		lnw=1
	n_wt,n_nw,n_pos=parse_variants(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web,dbfasta,fprog)
	if n_wt=='' or n_nw=='':
		print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		sys.exit(1)
	if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
		print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
		sys.exit(1)
	if pklcod=='':
		(nuc,seq,seq_input,cons_input,r_cod)=get_snv_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
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
		print >> sys.stderr, 'ERROR: Incorrect conservation data for position',ichr,ipos
		sys.exit(1)
	if cons_input2==[]:
		print >> sys.stderr, 'ERROR: Incorrect conservation data for position',ichr,pos
		sys.exit(1)
	if cons_input1==[]:
		cons_input1=[0.0 for i in range(2*win+1)]
		print >> sys.stderr,'WARNING: '+dat+' data not found for position',ichr,pos
	try:
		model=joblib.load(modfile)
	except:
		print >> sys.stderr,'ERROR: Program not able to load modfile. Please check that you have installed a compatible version joblib.'
		sys.exit(1)
	p_cod=0
	cod='No'
	if r_cod!=[]:
		p_cod=1
		cod='Yes'
	if pklcod=='':
		X=[seq_input + cons_input1+ cons_input2 ]
		y_pred,y_fdrs,c_pred=prediction(X,model)
		v_fdr=[y_fdrs[0][0],y_fdrs[0][1]]
	else:
		X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
		y_pred,y_fdrs,c_pred=prediction(X,model)
		v_fdr=[y_fdrs[0][2],y_fdrs[0][3]]
	if y_pred==[]:
		print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
		print '\t'.join([str(i) for i in [ichr,ipos,wt,nw] ])+'\tNA\tNA\tNA\tNA\tNA\tNA'
	else:
		print '#CHROM\tPOS\tREF\tALT\tCODING\tPREDICTION\tSCORE\tFDR\t'+dat+'\tAvg'+dat
		if vyear=='2017':
			pp100=cons_input2[win]
			avgpp100=sum(cons_input2)/float(len(cons_input2))
		else:
			pp100=cons_input1[win]
			avgpp100=sum(cons_input1)/float(len(cons_input1))
		if c_pred[0] == "Pathogenic": d_fdr=v_fdr[0]
		if c_pred[0] == "Benign": d_fdr=v_fdr[1]
		if ipos!=n_pos: cod=cod+'('+n_wt+str(n_pos)+n_nw+')'
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,cod,c_pred[0],'%.3f' %y_pred[0],'%.3f' %d_fdr,'%.3f' %pp100,'%.3f' %avgpp100])
	return


def make_vcffile_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP100way.bw','hg38.phyloP470way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph',inputfile='prova_input.txt'):
	v_input=[]
	try:
		model1=joblib.load(modfile[0])
		model2=joblib.load(modfile[1])
	except:
		print >> sys.stderr,'ERROR: Program not able to load modfile. Please check that you have installed a compatible version joblib.'
		sys.exit(1)
	proc = subprocess.Popen([prog_cat,'-f',namefile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	dat='PhyloP100'
	c=1
	for line in stdout.split('\n'):
		if line == '': continue
		if line.find('##')==0:
			print line
			continue
		line='\t'.join(line.split())
		if line.find('#')==0:
			if line.find('#CHROM')==0: line=line+'\tCODING\tPREDICTION\tSCORE\tFDR\t'+dat+'\tAvg'+dat
			print line
			continue
		v=line.rstrip().split('\t')
		if len(v)<5:
			print >> sys.stderr,'ERROR: Incorrect line ',c
			print >> sys.stderr, 'ERROR:', line
			continue
		if fpass and len(v)>6 and v[6]!='PASS':
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		(ichr,pos,rs,wt,nw)=tuple(v[:5])
		if wt==nw or nw.find(',')>-1:
			print >> sys.stderr, 'ERROR: Incorrect input line.',ichr,pos,wt,nw
			#print '\t'.join(str(i) for i in [ichr,pos,",",wt,nw,'NA','NA','NA','NA','NA','NA'])
			continue
		nchr=ichr
		if nchr.find('chr')==-1: nchr='chr'+ichr
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect input data. The VCF input file should have al least 5 columns (chr,position,id,ref,alt).'
			print >> sys.stderr, 'ERROR:', line
			continue
		lwt=len(wt)
		lnw=len(nw)
		if wt=='-':
			lwt=1
			lnw+=1
		if nw=='-':
			lwt+=1
			lnw=1
		n_wt,n_nw,n_pos=parse_variants(nchr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web,dbfasta,fprog)
		if n_wt=='' or n_nw=='':
			print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
			print >> sys.stderr, 'ERROR:', line
			continue
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
			print >> sys.stderr, line
			continue
		if len(wt)==1 and len(nw)==1:
			(nuc,seq,seq_input,cons_input,r_cod)=get_snv_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
		else:
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
		if seq=='':
			print >> sys.stderr, 'ERROR: Sequence not found for line',c,ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if seq_input==[]:
			print >> sys.stderr, 'ERROR: Incorrect nucleotide in line '+str(c)+'. Genome location:',ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if cons_input!=[]:
			cons_input1=cons_input[0]
			cons_input2=cons_input[1]
		else:
			print >> sys.stderr, 'ERROR: Incorrect conservation data for line',c,ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		#if cons_input1==[] or cons_input2==[]:
		#Check only P100
		if cons_input2==[]:
			print >> sys.stderr, 'ERROR: Incorrect conservation data for line',c,ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if cons_input1==[]:
			cons_input1=[0.0 for i in range(2*win+1)]
			print >> sys.stderr,'WARNING: '+dat+' data not found in line',c,ichr,pos
		p_cod=0
		cod='No'
		if r_cod!=[]:
			p_cod=1
			cod='Yes'
		if len(wt)==1 and len(nw)==1 and wt!='-' and nw!='-':
			X=[seq_input + cons_input1+ cons_input2 ]
			y_pred,y_fdrs,c_pred=prediction(X,model1)
			v_fdr=[y_fdrs[0][0],y_fdrs[0][1]]
		else:
			X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
			y_pred,y_fdrs,c_pred=prediction(X,model2)
			v_fdr=[y_fdrs[0][2],y_fdrs[0][3]]
		if y_pred==[]:
			print >> sys.stderr,'WARNING: Variant not scored. Check modfile and input'
			print >> sys.stderr, 'WARNING:',line
			#print line+'\tNA\tNA\tNA\tNA\tNA'
			continue
		if vyear=='2017':
			pp100=cons_input2[win]
			avgpp100=sum(cons_input2)/float(len(cons_input2))
		else:
			pp100=cons_input1[win]
                        avgpp100=sum(cons_input1)/float(len(cons_input1))
		#print pp100,avgpp100,cons_input2
		if c_pred[0] == "Pathogenic": d_fdr=v_fdr[0]
		if c_pred[0] == "Benign": d_fdr=v_fdr[1]
		if ipos!=n_pos: cod=cod+'('+n_wt+str(n_pos)+n_nw+')'
		print line+'\t'+'%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f' %(cod,c_pred[0],y_pred[0],d_fdr,pp100,avgpp100)
		v_input.append(line+'\t'+'\t'.join([str(i) for i in X[0]])+'\n')
		#print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])
		c=c+1
	if inputfile!='': open(inputfile,'w').writelines(v_input)
	return


def make_vcffile_multialleles_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP100way.bw','hg38.phyloP470way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph',sep='@'):
	nucs='ACGTN'
	try:
		model1=joblib.load(modfile[0])
		model2=joblib.load(modfile[1])
	except:
		print >> sys.stderr,'ERROR: Program not able to load modfile. Please check that you have installed a compatible version joblib.'
		sys.exit(1)
	list_pred=[]
	proc = subprocess.Popen([prog_cat,'-f',namefile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	dat='PhyloP100'
	c=0
	#print '#CHROM\tPOS\tID\tREF\tALT\tCODING\tPREDICTION\tSCORE\tFDR\t'+dat+'\tAvg'+dat
	for line in stdout.split('\n'):
		list_pred=[]
		if line.find('##')==0:
			print line
			continue
		line='\t'.join(line.split())
		if line.find('#')==0:
			if line.find('#CHROM')==0: line=line+'\tCODING\tPREDICTION\tSCORE\tFDR\t'+dat+'\tAvg'+dat
			print line
			continue
		v=line.rstrip().split('\t')
		c=c+1
		if len(v)<5:
			print >> sys.stderr,'ERROR: Incorrect line ',c
			print >> sys.stderr, 'ERROR:', line
			continue
		if fpass and len(v)>6 and v[6]!='PASS':
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		(ichr,pos,rs,wt,nw)=tuple(v[:5])
		list_nw=nw.split(',')
		nchr=ichr
		if nchr.find('chr')==-1: nchr='chr'+ichr
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect input data. The VCF input file should have al least 5 columns (chr,position,id,ref,alt).'
			print >> sys.stderr, 'ERROR:', line
			continue
		##Option for selecting snv
		#if len(wt)!=1 or len(nw)!=1 or nucs.find(wt)==-1 or nucs.find(nw)==-1 or wt==nw:
		#	print >> sys.stderr, 'ERROR: Not single nucleotide variant',','.join([nchr,pos,wt,nw])
		#	continue
		for inw in list_nw:
			if wt==inw:
				print >> sys.stderr, 'ERROR: Incorrect nucleotide in line '+str(c)+'. Genome location:',ichr,pos
				print >> sys.stderr, 'ERROR:', line
				continue
			lwt=len(wt)
			lnw=len(inw)
			if wt=='-':
				lwt=1
				lnw+=1
			if inw=='-':
				lwt+=1
				lnw=1
			n_wt,n_nw,n_pos=parse_variants(nchr,ipos,wt,inw,ucsc_exe,ucsc_dbs,web,dbfasta,fprog)
			if n_wt=='' or n_nw=='':
				print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,inw
				print >> sys.stderr, 'ERROR:', line
				continue
			if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
				print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,inw
				print >> sys.stderr, 'ERROR:', line
				continue
			if len(wt)==1 and len(inw)==1:
				r_cod=[]
				(nuc,seq,seq_input,cons_input,r_cod)=get_snv_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
			else:
				(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
			if seq=='':
				print >> sys.stderr, 'ERROR: Sequence not found for line '+str(c)+'. Genome location:',ichr,pos
				print >> sys.stderr, 'ERROR:', line
				#list_pred.append(5*['NA'])
				continue
			if seq_input==[]:
				print >> sys.stderr, 'ERROR: Incorrect nucleotide in line '+str(c)+'. Genome location:',ichr,pos
				print >> sys.stderr, 'ERROR:', line
				#list_pred.append(5*['NA'])
				continue
			if cons_input!=[]:
				cons_input1=cons_input[0]
				cons_input2=cons_input[1]
			else:
				print >> sys.stderr, 'ERROR: Incorrect conservation data for line '+str(c)+'. Genome location:',ichr,pos
				print >> sys.stderr, 'ERROR:', line
				#list_pred.append(5*['NA'])
				continue
			#if cons/1_input1==[] or cons_input2==[]:
			#Check only P100
			if cons_input2==[]:
				print >> sys.stderr, 'ERROR: Incorrect conservation data for line '+str(c)+'. Genome location:',ichr,pos
				print >> sys.stderr, 'ERROR:', line
				#list_pred.append(5*['NA'])
				continue
			if cons_input1==[]:
				cons_input1=[0.0 for i in range(2*win+1)]
				print >> sys.stderr,'WARNING: '+dat+' data not found for line '+str(c)+'. Genome location:',ichr,pos
			p_cod=0
			cod='No'
			if r_cod!=[]:
				p_cod=1
				cod='Yes'
			if len(wt)==1 and len(inw)==1 and wt!='-' and inw!='-':
				X=[seq_input + cons_input1+ cons_input2 ]
				y_pred,y_fdrs,c_pred=prediction(X,model1)
				v_fdr=[y_fdrs[0][0],y_fdrs[0][1]]
			else:
				X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
				y_pred,y_fdrs,c_pred=prediction(X,model2)
				v_fdr=[y_fdrs[0][2],y_fdrs[0][3]]
			if y_pred==[]:
				print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
				print >> sys.stderr, 'WARNING:', line
				#list_pred.append(5*['NA'])
				continue
			if vyear=='2017':	
				pp100=cons_input2[win]
				avgpp100=sum(cons_input2)/float(len(cons_input2))
			else:
				pp100=cons_input1[win]
				avgpp100=sum(cons_input1)/float(len(cons_input1))
			if c_pred[0] == "Pathogenic": d_fdr=v_fdr[0]
			if c_pred[0] == "Benign": d_fdr=v_fdr[1]
			if ipos!=n_pos: cod=cod+'('+n_wt+str(n_pos)+n_nw+')'
			list_pred.append(['%s' %cod,'%s' %c_pred[0],'%.3f' %y_pred[0],'%.3f' %d_fdr,'%.3f' %pp100,'%.3f' %avgpp100])
		#print list_pred
		if len(list_pred)!=len(list_nw):
			continue
			#out_data=5*('NA',)
		else:
			out_data=tuple([ ':'.join(single_pred[i]  for single_pred in list_pred)  for i in range(6)])
		#print pp100,avgpp100,cons_input2
		if hg == 'hg19': line=v[5].replace(sep,'\t')
		print line+'\t'+'%s\t%s\t%s\t%s\t%s\t%s' %out_data
		#print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])
	return


def make_tsvfile_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web=False,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP100way.bw','hg38.phyloP470way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph',sep='@'):
	nucs='ACGTN'
	dat='PhyloP100'
	v_input=[]
	try:
		model1=joblib.load(modfile[0])
		model2=joblib.load(modfile[1])
	except:
		print >> sys.stderr,'ERROR: Program not able to load modfile. Please check that you have installed a compatible version joblib.'
		sys.exit(1)
	proc = subprocess.Popen([prog_cat,'-f',namefile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()
	c=0
	print '#CHROM\tPOS\tREF\tALT\tCODING\tPREDICTION\tSCORE\tFDR\t'+dat+'\tAvg'+dat
	for line in stdout.split('\n'):
		if line == '': continue
		v=line.rstrip().split()
		c=c+1
		if len(v)<4:
			print >> sys.stderr,'ERROR: Incorrect line ',c
			print >> sys.stderr, 'ERROR:', line
			continue
		(ichr,pos,wt,nw)=tuple(v[:4])
		if wt==nw or nw.find(',')>-1:
			print >> sys.stderr, 'ERROR: Incorrect input line.',ichr,pos,wt,nw
			print >> sys.stderr, 'ERROR:', line
			continue
		nchr=ichr
		if nchr.find('chr')==-1: nchr='chr'+ichr
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect input data. The tsv input file should have has four columns (chr,position,ref,alt).'
			print >> sys.stderr, 'ERROR:', line
			continue
			sys.exit(1)
		lwt=len(wt)
		lnw=len(nw)
		##Option for selecting snv
		#if lwt!=1 or lnw!=1 or nucs.find(wt)==-1 or nucs.find(nw)==-1 or wt==nw:
		#	print >> sys.stderr, 'ERROR: Not single nucleotide variant',','.join([nchr,pos,wt,nw])
		#	continue
		if wt=='-':
			lwt=1
			lnw+=1
		if nw=='-':
			lwt+=1
			lnw=1
		n_wt,n_nw,n_pos=parse_variants(nchr,ipos,wt,nw,ucsc_exe,ucsc_dbs,web,dbfasta,fprog)
		if n_wt=='' or n_nw=='':
			print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
			print >> sys.stderr, 'ERROR:', line
			continue
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
			print >> sys.stderr, 'ERROR:', line
			continue
		if len(wt)==1 and len(nw)==1:
			r_cod=[]
			(nuc,seq,seq_input,cons_input,r_cod)=get_snv_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
		else:
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
		if seq=='':
			print >> sys.stderr, 'ERROR: Sequence not found for line '+str(c)+'. Genome location:',ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if seq_input==[]:
			print >> sys.stderr, 'ERROR: Incorrect nucleotide in line '+str(c)+'. Genome location:',ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if cons_input!=[]:
			cons_input1=cons_input[0]
			cons_input2=cons_input[1]
		else:
			print >> sys.stderr, 'ERROR: Incorrect conservation data for line '+str(c)+'. Genome location:',ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		#if cons_input1==[] or cons_input2==[]:
		#Check only P100
		if cons_input2==[]:
			print >> sys.stderr, 'ERROR: Incorrect conservation data for line '+str(c)+'. Genome location:',ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if cons_input1==[]:
			cons_input1=[0.0 for i in range(2*win+1)]
			print >> sys.stderr,'WARNING: PhyloP100 data not found for line',c,'mutated site',ichr,pos
		p_cod=0
		cod='No'
		if r_cod!=[]:
			p_cod=1
			cod='Yes'
		if len(wt)==1 and len(nw)==1 and wt!='-' and nw!='-':
			X=[seq_input + cons_input1+ cons_input2 ]
			y_pred,y_fdrs,c_pred=prediction(X,model1)
			v_fdr=[y_fdrs[0][0],y_fdrs[0][1]]
		else:
			X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
			y_pred,y_fdrs,c_pred=prediction(X,model2)
			v_fdr=[y_fdrs[0][2],y_fdrs[0][3]]
		if y_pred==[]:
			print >> sys.stderr,'WARNING: Variant not scored. Check modfile and input'
			print >> sys.stderr, 'WARNING:', line
			continue
		if vyear=='2017':
			pp100=cons_input2[win]
			avgpp100=sum(cons_input2)/float(len(cons_input2))
		else:
			pp100=cons_input1[win]
			avgpp100=sum(cons_input1)/float(len(cons_input1))
		#print pp100,avgpp100,cons_input2
		if c_pred[0] == "Pathogenic": d_fdr=v_fdr[0]
		if c_pred[0] == "Benign": d_fdr=v_fdr[1]
		if ipos!=n_pos: cod=cod+'('+n_wt+str(n_pos)+n_nw+')'
		if hg == 'hg19': line=v[4].replace(sep,'\t')
		print line+'\t'+'%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f' %(cod,c_pred[0],y_pred[0],d_fdr,pp100,avgpp100)
		#print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])
	return



def make_file_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web=False,win=2,s='\t',dbfasta='hg38.2bit',dbpps=['hg38.phyloP100way.bw','hg38.phyloP470way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	try:
		model1=joblib.load(modfile[0])
		model2=joblib.load(modfile[1])
	except:
		print >> sys.stderr,'ERROR: Program not able to load modfile. Please check that you have installed a compatible version joblib.'
		sys.exit(1)
	f=open(namefile)
	c=1
	dat='PhyloP100'
	print "#CHROM\tPOS\tREF\tALT\tCODING\tPREDICTION\tSCORE\tFDR\t'+dat+'\tAvg"+dat
	for line in f:
		v=line.rstrip().split(s)
		if len(v)<4:
			print >> sys.stderr,'ERROR: Incorrect line ',c,line.rstrip()
			print >> sys.stderr, line
			continue
		(ichr,pos,wt,nw)=v[:4]
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect genome location',pos,'. Check your input file'
			print >> sys.stderr, 'ERROR:', line
			continue
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
			print >> sys.stderr, 'ERROR:', line
			continue
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
			print >> sys.stderr, 'ERROR:', line
			continue
		if wt==nw or nw.find(',')>-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
			print >> sys.stderr, 'ERROR:', line
			continue
		if len(wt)==1 and len(nw)==1:
			r_cod=[]
			(nuc,seq,seq_input,cons_input,r_cod)=get_snv_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
		else:
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,web,win,dbfasta,dbpps,pklcod,fprog,cprog)
		if seq=='':
			print >> sys.stderr, 'ERROR: Sequence not found for line',c,ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if seq_input==[]:
			print >> sys.stderr, 'ERROR: Incorrect nucleotide in line '+str(c)+'. Genome location:',ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if cons_input!=[]:
			cons_input1=cons_input[0]
			cons_input2=cons_input[1]
		else:
			print >> sys.stderr, 'ERROR: Incorrect conservation data for line',c,ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if cons_input2==[]:
			print >> sys.stderr, 'ERROR: Incorrect conservation data for line',c,ichr,pos
			print >> sys.stderr, 'ERROR:', line
			continue
		if cons_input1==[]:
			cons_input1=[0.0 for i in range(2*win+1)]
			print >> sys.stderr,'WARNING: PhyloP100 data not found for line',c,ichr,pos
		p_cod=0
		cod='No'
		if r_cod!=[]:
			p_cod=1
			cod='Yes'
		if len(wt)==1 and len(nw)==1 and wt!='-' and nw!='-':
			X=[seq_input + cons_input1+ cons_input2 ]
			y_pred,y_fdrs,c_pred=prediction(X,model1)
			v_fdr=[y_fdrs[0][0],y_fdrs[0][1]]
		else:
			X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
			y_pred,y_fdrs,c_pred=prediction(X,model2)
			v_fdr=[y_fdrs[0][2],y_fdrs[0][3]]
		if y_pred==[]:
			print >> sys.stderr,'WARNING: Variant not scored. Check modfile and input'
			print >> sys.stderr, 'WARNING:', line
			continue
		if vyear=='2017':
			pp100=cons_input2[win]
			avgpp100=sum(cons_input2)/float(len(cons_input2))
		else:
			pp100=cons_input1[win]
			avgpp100=sum(cons_input1)/float(len(cons_input1))
		if c_pred[0] == "Pathogenic": d_fdr=v_fdr[0]
		if c_pred[0] == "Benign": d_fdr=v_fdr[1]
		if ipos!=n_pos: cod=cod+'('+n_wt+str(n_pos)+n_nw+')'
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,cod,c_pred[0],'%.3f' %y_pred[0],'%.3f' %d_fdr,pp100,avgpp100])
	return




def prediction(X,model,fdr_file='fdr_mean'):
	y_pred=[]
	y_fdrs=[]
	c_pred=[]
	fdr_file=fdr_file+'_'+vyear+'.pkl'
	nf=model.n_features
	if len(X[0]) != nf:
		print >> sys.stderr,'ERROR: Model expecting',nf,'features. Input vector with',len(X[0]),'features. Check the modfile or the window size.'
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


def liftover(filein,fileout,sep='@'):
	vhead=[]
	ftmp=filein+'.inlift'
	ferr=filein+'.erlift'
	f=open(filein)
	fw=open(ftmp,'w')
	for line in f:
		if line.find('#')==0:
			vhead.append(line)
			continue
		line=line.replace(',','\t')
		v=line.rstrip().split()
		chrom=v[0]
		if chrom.find('chr')==-1: chrom='chr'+v[0]
		fw.write(chrom+'\t'+v[1]+'\t'+str(int(v[1])+1)+'\t'+sep.join(v)+'\n')
	fw.close()
	out=getstatusoutput(ucsc_dir+'/exe/liftOver '+ftmp+' '+ucsc_dir+'/hg19/hg19ToHg38.over.chain.gz '+fileout+' /dev/stderr')
	unmap=out[1].split('\n')
	if len(unmap)>1: unmap=unmap[2:]
	if len(unmap)>0:
		for l in unmap:
			if l.find('#')==0: continue
			print >>sys.stderr,'ERROR: Incorrect mapping. Check input line',' '.join(l.split()[-1].split(sep))
		#open(ferr,'w').write('\n'.join())
	f=open(fileout)
	vout=[]
	for line in f:
		v=line.rstrip().split()
		t=v[3].split(sep)
		if len(t)<4:
			print >> sys.stderr,"ERROR: Incorrect mapping of",v[3]
			continue
		if line.find('chr')==0: v[0]=v[0][3:]
		if vcf:
			vout.append('\t'.join(v[:2])+'\t.\t'+t[3]+'\t'+t[4]+'\t'+v[3]+'\n')
		else:
			vout.append('\t'.join(v[:2])+'\t'+t[2]+'\t'+t[3]+'\t'+v[3]+'\n')
	if vcf: vout=vhead+vout
	open(fileout,'w').writelines(vout)
	return unmap


def get_options():
	global hg, coord, vcf, fpass, web, win, vyear
	import optparse
	desc = 'Script for scoring single nucleotide variants'
	parser = optparse.OptionParser("usage: %prog variant_file", description=desc)
	parser.add_option('-m','--mod-file', action='store', type='string', dest='mfile', help='Model file')
	parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='hg38', help='Genome version')
	parser.add_option('-v','--verbose', action='store_true', dest='ver', default=False, help='Verbose mode')
	parser.add_option('-c','--coordinate', action='store_true', dest='coord', default=False, help='Coordinate input')
	parser.add_option('--vcf', action='store_true', dest='vcf', default=False, help='VCF file input')
	parser.add_option('--web', action='store_true', dest='web', default=False, help='Use UCSC web files')
	parser.add_option('--pass', action='store_true', dest='fpass', default=False, help='Predict only PASS variants. Check column 7 in vcf file')
	parser.add_option('-y','--version', action='store', type='string', dest='vyear', default='2022', help='Version year')
	(options, args) = parser.parse_args()
	outfile = ''
	#modfile = [prog_dir + '/data/model/snv_model_w5_p7_500.pkl',prog_dir + '/data/model/indel_model_w5_p7_500.pkl']
	hg='hg38'
	coord=False
	vcf=False
	fpass=False
	web=False
	win=2
	vyear='2022'
	if options.mfile: modfile=options.mfile
	if options.hg.lower()=='hg19': hg='hg19'
	if options.coord: coord = True
	if options.fpass: fpass=True
	if options.web: web=True
	if options.vyear=='2017': vyear='2017'
	if options.vcf: vcf = True
	__builtin__.verbose=False
	if options.ver: __builtin__.verbose=True
	fasta=hg38['fasta']
	if vyear=='2017':
		dbpps=[hg38['phylop'][0],hg38['phylop'][1]]
	else:
		dbpps=[hg38['phylop'][1],hg38['phylop'][2]]
	#dbpps=hg38['phylop']+hg38['phastc']
	pklcod=hg38['coding']
	modfile = [prog_dat + '/snv_model_'+vyear+'_w5_p7_500_hg38.pkl',prog_dat + '/indel_model_'+vyear+'_w5_p7_500_hg38.pkl']
	if not os.path.isfile(modfile[0]) or not os.path.isfile(modfile[1]):
                print >> sys.stderr,'ERROR: Data model files not found',modfile
		sys.exit(1)
	opts=(outfile,modfile,fasta,dbpps,vyear,pklcod)
	return args,opts



if __name__ == '__main__':
	global_vars()
	from sklearn.externals import joblib
	from score_variants import parse_variants, get_snv_input, get_indel_input
	args,opts=get_options()
	(outfile,modfile,fasta,dbpps,vyear,pklcod)=opts
	ucsc_dbs=ucsc_dir+'/hg38'
	if len(args)>0:
		namefile=args[0]
		if not os.path.isfile(namefile):
			print >> sys.stderr,'ERROR: Input file not found',namefile
			sys.exit(1)
		if hg=='hg19':
			input_hg38=namefile+'.hg38'
			olift=liftover(namefile,input_hg38)
			namefile=input_hg38
			#hg='hg38'
			#vcf=True
		if vcf:
			make_vcffile_multialleles_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web,win,fasta,dbpps,pklcod)
				#make_vcffile_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web,win,fasta,dbpps,pklcod)
		else:
			make_tsvfile_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,web,win,fasta,dbpps,pklcod)
	else:
		print 'predict_variants.py variant_file -g hg_version'
