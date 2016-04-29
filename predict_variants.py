#!/usr/bin/env python
import os, sys, subprocess, pickle
import  __builtin__
from commands import getstatusoutput
from sklearn.externals import joblib


def global_vars():
	global tool_dir, prog_dir, prog_dat, ucsc_dir, ucsc_exe, verbose, hg19, hg38, prog_cat
	prog_dir = os.path.dirname(os.path.abspath(__file__))
	tool_dir = prog_dir+'/tools'
	prog_dat = prog_dir+'/data/model'
	ucsc_dir = prog_dir+'/ucsc'
	ucsc_exe = ucsc_dir+'/exe'
	prog_cat = 'zcat'
	sys.path.append(tool_dir)
	hg19={}
	hg19['fasta']='hg19.2bit'
	hg19['phylop']=['hg19.phyloP46way.primate.bw','hg19.phyloP46way.placental.bw','hg19.100way.phyloP100way.bw']
	hg19['coding']='hg19_coding.bed'
	hg38={}
	hg38['fasta']='hg38.2bit'
	hg38['phylop']=['hg38.phyloP7way.bw','hg38.phyloP20way.bw','hg38.phyloP100way.bw']
	hg38['coding']='hg38_coding.bed'
	return



def make_prediction(ichr,ipos,wt,nw,modfile,ucsc_exe,ucsc_dbs,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	lwt=len(wt)
	lnw=len(nw)
	if wt=='-':
		lwt=1
		lnw+=1
	if nw=='-':
		lwt+=1
		lnw=1
	n_wt,n_nw,n_pos=parse_variants(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,dbfasta,fprog)
	if n_wt=='' or n_nw=='':
		print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		sys.exit()
	if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
		print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
		sys.exit()
	if pklcod=='':
		r_cod=[]
		(nuc,seq,seq_input,cons_input)=get_snv_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,fprog,cprog)
	else:
		(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,fprog,cprog)
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
	model=joblib.load(modfile)
	if pklcod=='':
		X=[seq_input + cons_input1+ cons_input2 ]
	else:
		p_cod=0
		if r_cod!=[]: p_cod=1
		X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
	y_pred,y_fdrs,c_pred=prediction(X,model)
	if y_pred==[]:
		print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
		print '\t'.join([str(i) for i in [ichr,ipos,wt+','+nw] ])+'\tNA\tNA\tNA\tNA\tNA\tNA'
	else:
		print "#CHROM\tPOS\tREF\tALT\tPREDICTION\tSCORE\tFDR\t1-NPV\tPhyloP100\tAvgPhyloP100"
		pp100=cons_input2[win]
		avgpp100=sum(cons_input2)/float(len(cons_input2))
		print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,c_pred[0],'%.3f' %y_pred[0],'%.3f' %y_fdrs[0][0],'%.3f' %y_fdrs[0][1],'%.3f' %pp100,'%.3f' %avgpp100])
	return


def make_vcffile_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	model1=joblib.load(modfile[0])
	model2=joblib.load(modfile[1])	
	proc = subprocess.Popen([prog_cat,'-f',namefile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()        
	c=1
	for line in stdout.split('\n'):
		if line == '': continue 
		if line[0]=='#':
			if line.find('#CHROM')==0: line=line+'\tPREDICTION\tSCORE\tFDR\t1-NPV\tPhyloP100\tAvgPhyloP100'
			print line
			continue 	
		v=line.rstrip().split('\t')
		if len(v)<5:
			print >> sys.stderr,'WARNING: Incorrect line ',c	
			print line
			continue
		if fpass and len(v)>6 and v[6]!='PASS':
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		(ichr,pos,rs,wt,nw)=tuple(v[:5])
		if wt==nw or nw.find(',')>-1:
			print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'NA','NA','NA','NA','NA','NA'])
			continue
		nchr=ichr
		if nchr.find('chr')==-1: nchr='chr'+ichr
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect input data. The VCF input file should have al least 5 columns (chr,position,id,ref,alt).'
			sys.exit(1)
		lwt=len(wt)
		lnw=len(nw)
		if wt=='-':
			lwt=1
			lnw+=1
		if nw=='-':
			lwt+=1
			lnw=1
		n_wt,n_nw,n_pos=parse_variants(nchr,ipos,wt,nw,ucsc_exe,ucsc_dbs,dbfasta,fprog)
		if n_wt=='' or n_nw=='':
			print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw		
		if len(wt)==1 and len(nw)==1:
			r_cod=[]
			(nuc,seq,seq_input,cons_input)=get_snv_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,fprog,cprog)
		else: 
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,fprog,cprog)
		if seq=='': 
			print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		if seq_input==[]: 
			print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		if cons_input!=[]:
			cons_input1=cons_input[0]
			cons_input2=cons_input[1]
		else:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		#if cons_input1==[] or cons_input2==[]:
		#Check only P100
		if cons_input2==[]:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		if cons_input1==[]: cons_input1=[0.0 for i in range(2*win+1)]
		if len(wt)==1 and len(nw)==1:
			X=[seq_input + cons_input1+ cons_input2 ]
			y_pred,y_fdrs,c_pred=prediction(X,model1)
		else:
			p_cod=0
			if r_cod!=[]: p_cod=1
			X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
			y_pred,y_fdrs,c_pred=prediction(X,model2)		
		if y_pred==[]:
			print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		pp100=cons_input2[win]
		avgpp100=sum(cons_input2)/float(len(cons_input2))	
		#print pp100,avgpp100,cons_input2
		print line+'\t'+'%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' %(c_pred[0],y_pred[0],y_fdrs[0][0],y_fdrs[0][1],pp100,avgpp100)
		#print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])	
		c=c+1
	return 


def make_vcffile_multialleles_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	model1=joblib.load(modfile[0])
	model2=joblib.load(modfile[1])
	list_pred=[]	
	proc = subprocess.Popen([prog_cat,'-f',namefile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()        
	c=1
	for line in stdout.split('\n'):
		list_pred=[]
		if line == '': continue 
		if line[0]=='#':
			if line.find('#CHROM')==0: line=line+'\tPREDICTION\tSCORE\tFDR\t1-NPV\tPhyloP100\tAvgPhyloP100'
			print line
			continue 	
		v=line.rstrip().split('\t')
		if len(v)<5:
			print >> sys.stderr,'WARNING: Incorrect line ',c	
			print line
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
			sys.exit(1)
		for inw in list_nw:
			if wt==inw:
				print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
				print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
				continue
			lwt=len(wt)
			lnw=len(inw)
			if wt=='-':
				lwt=1
				lnw+=1
			if inw=='-':
				lwt+=1
				lnw=1
			n_wt,n_nw,n_pos=parse_variants(nchr,ipos,wt,inw,ucsc_exe,ucsc_dbs,dbfasta,fprog)
			if n_wt=='' or n_nw=='':
				print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,inw
			if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
				print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,inw		
			if len(wt)==1 and len(inw)==1:
				r_cod=[]
				(nuc,seq,seq_input,cons_input)=get_snv_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,fprog,cprog)
			else: 
				(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,fprog,cprog)
			if seq=='': 
				print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos
				list_pred.append(6*['NA'])
				continue
			if seq_input==[]: 
				print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
				list_pred.append(6*['NA'])
				continue
			if cons_input!=[]:
				cons_input1=cons_input[0]
				cons_input2=cons_input[1]
			else:
				print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
				list_pred.append(6*['NA'])
				continue
			#if cons_input1==[] or cons_input2==[]:
			#Check only P100
			if cons_input2==[]:
				print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
				list_pred.append(6*['NA'])
				continue
			if cons_input1==[]: cons_input1=[0.0 for i in range(2*win+1)]
			if len(wt)==1 and len(inw)==1:
				X=[seq_input + cons_input1+ cons_input2 ]
				y_pred,y_fdrs,c_pred=prediction(X,model1)
			else:
				p_cod=0
				if r_cod!=[]: p_cod=1
				X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
				y_pred,y_fdrs,c_pred=prediction(X,model2)		
			if y_pred==[]:
				print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
				list_pred.append(6*['NA'])
				continue
			pp100=cons_input2[win]
			avgpp100=sum(cons_input2)/float(len(cons_input2))
			list_pred.append(['%s' %c_pred[0],'%.3f' %y_pred[0],'%.3f' %y_fdrs[0][0],'%.3f' %y_fdrs[0][1],'%.3f' %pp100,'%.3f' %avgpp100])
		#print list_pred
		if list_pred==[]:
			out_data=6*('NA',)
		else:
			out_data=tuple([ ':'.join(single_pred[i]  for single_pred in list_pred)  for i in range(6)])	
		#print pp100,avgpp100,cons_input2
		print line+'\t'+'%s\t%s\t%s\t%s\t%s\t%s' %out_data		
		#print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])	
		c=c+1
	return 


def make_tsvfile_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win=2,dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	model1=joblib.load(modfile[0])
	model2=joblib.load(modfile[1])	
	proc = subprocess.Popen([prog_cat,'-f',namefile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = proc.communicate()        
	c=1
	print "#CHROM\tPOS\tREF\tALT\tPREDICTION\tSCORE\tFDR\t1-NPV\tPhyloP100\tAvgPhyloP100"
	for line in stdout.split('\n'):
		if line == '': continue 
		v=line.rstrip().split()
		if len(v)<4:
			print >> sys.stderr,'WARNING: Incorrect line ',c	
			print line
			continue
		(ichr,pos,wt,nw)=tuple(v[:4])
		if wt==nw or nw.find(',')>-1:
			print >> sys.stderr, 'WARNING: Incorrect input line.',ichr,pos,wt,nw
			#print '\t'.join(str(i) for i in [ichr,pos,wt,nw,'NA','NA','NA','NA','NA','NA'])
			continue
		nchr=ichr
		if nchr.find('chr')==-1: nchr='chr'+ichr
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect input data. The tsv input file should have has four columns (chr,position,ref,alt).'
			sys.exit(1)
		lwt=len(wt)
		lnw=len(nw)
		if wt=='-':
			lwt=1
			lnw+=1
		if nw=='-':
			lwt+=1
			lnw=1
		n_wt,n_nw,n_pos=parse_variants(nchr,ipos,wt,nw,ucsc_exe,ucsc_dbs,dbfasta,fprog)
		if n_wt=='' or n_nw=='':
			print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw		
		if len(wt)==1 and len(nw)==1:
			r_cod=[]
			(nuc,seq,seq_input,cons_input)=get_snv_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,fprog,cprog)
		else: 
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(nchr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,fprog,cprog)
		if seq=='': 
			print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos
			#print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		if seq_input==[]: 
			print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
			#print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		if cons_input!=[]:
			cons_input1=cons_input[0]
			cons_input2=cons_input[1]
		else:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
			#print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		#if cons_input1==[] or cons_input2==[]:
		#Check only P100
		if cons_input2==[]:
			print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
			#print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		if cons_input1==[]: cons_input1=[0.0 for i in range(2*win+1)]
		if len(wt)==1 and len(nw)==1:
			X=[seq_input + cons_input1+ cons_input2 ]
			y_pred,y_fdrs,c_pred=prediction(X,model1)
		else:
			p_cod=0
			if r_cod!=[]: p_cod=1
			X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
			y_pred,y_fdrs,c_pred=prediction(X,model2)		
		if y_pred==[]:
			print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
			#print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		pp100=cons_input2[win]
		avgpp100=sum(cons_input2)/float(len(cons_input2))	
		#print pp100,avgpp100,cons_input2
		print line+'\t'+'%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' %(c_pred[0],y_pred[0],y_fdrs[0][0],y_fdrs[0][1],pp100,avgpp100)
		#print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'%.4f' %y_pred[0]])	
		c=c+1
	return 



def make_file_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win=2,s='\t',dbfasta='hg38.2bit',dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw'],pklcod='hg38_coding.pkl',fprog='twoBitToFa',cprog='bigWigToBedGraph'):
	model1=joblib.load(modfile[0])
	model2=joblib.load(modfile[1])
	f=open(namefile)
	c=1
	print "#CHROM\tPOS\tREF\tALT\tPREDICTION\tSCORE\tFDR\t1-NPV\tPhyloP100\tAvgPhyloP100"
	for line in f:	
		v=line.rstrip().split(s)
		if len(v)<4: print >> sys.stderr,'WARNING: Incorrect line ',c,line.rstrip()
		(ichr,pos,wt,nw)=v[:4]
		try:
			ipos=int(pos)
		except:
			print >> sys.stderr,'ERROR: Incorrect genome location',pos,'. Check your input file'
			sys.exit(1)
		lwt=len(wt)
		lnw=len(nw)
		if wt=='-':
			lwt=1
			lnw+=1
		if nw=='-':
			lwt+=1
			lnw=1
		n_wt,n_nw,n_pos=parse_variants(ochr,ipos,wt,nw,ucsc_exe,ucsc_dbs,dbfasta,fprog)
		if n_wt=='' or n_nw=='':
			print >> sys.stderr, 'ERROR: Incorrect mutation mapping. Check position',ichr,ipos,wt,nw
		if 'ACGTN'.find(n_wt)==-1 or 'ACGTN'.find(n_nw)==-1:
			print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw		
		if wt==nw or nw.find(',')>-1:
			print '\t'.join(str(i) for i in [ichr,ipos,wt,nw,'NA','NA','NA','NA','NA','NA'])
			continue
		if len(wt)==1 and len(nw)==1:
			r_cod=[]
			(nuc,seq,seq_input,cons_input)=get_snv_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,fprog,cprog)
		else: 
			(nuc,seq,seq_input,cons_input,r_cod)=get_indel_input(ichr,n_pos,n_wt,n_nw,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,fprog,cprog)
		if seq=='': print >> sys.stderr, 'WARNING: Sequence not found for line',c,ichr,pos		
		if seq_input==[]: print >> sys.stderr, 'WARNING: Incorrect nucleotide in line',c,ichr,pos
		if cons_input!=[]:
			cons_input1=cons_input[0]
			cons_input2=cons_input[1]
		else:
			print >> sys.stderr,'WARNING: Variants not scored. Check modfile and input'
			print line+'\tNA\tNA\tNA\tNA\tNA\tNA'
			continue
		if cons_input1==[] or cons_input2==[]: print >> sys.stderr, 'WARNING: Incorrect conservation data in line',c,ichr,pos
		if len(wt)==1 and len(nw)==1:
			X=[seq_input + cons_input1+ cons_input2 ]
			y_pred,y_fdrs,c_pred=prediction(X,model1)
		else:
			p_cod=0
			if r_cod!=[]: p_cod=1
			X=[seq_input + cons_input1+ cons_input2 + [lwt, lnw, p_cod]]
			y_pred,y_fdrs,c_pred=prediction(X,model2)
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
	modfile = [prog_dir + '/data/model/snv_model_w5_p7_500.pkl',prog_dir + '/data/model/indel_model_w5_p7_500.pkl']
	hg='hg38'
	coord=False
	vcf=False
	fpass=False
	win=2
	if options.mfile: modfile=options.mfile
	if options.hg.lower()=='hg19': hg='hg19'
	if options.coord: coord = True
	if options.fpass: fpass=True
	if options.vcf: vcf = True
	__builtin__.verbose=False
	if options.ver: __builtin__.verbose=True
	if hg=='hg19':
		fasta=hg19['fasta']
		dbpps=[hg19['phylop'][0],hg19['phylop'][2]]
		pklcod=cod=hg19['coding']
	else:
		fasta=hg38['fasta']
		dbpps=[hg38['phylop'][0],hg38['phylop'][2]]
		#dbpps=hg38['phylop']+hg38['phastc']
		pklcod=hg38['coding']
	if not os.path.isfile(modfile[0]) or not os.path.isfile(modfile[1]):
                print >> sys.stderr,'ERROR: Data model files not found'
		sys.exit(1)
	opts=(outfile,modfile,fasta,dbpps,pklcod)
	return args,opts



if __name__ == '__main__':
	global_vars()
	from score_variants import parse_variants, get_snv_input, get_indel_input	
	args,opts=get_options()
	(outfile,modfile,fasta,dbpps,pklcod)=opts
	ucsc_dbs=ucsc_dir+'/'+hg
	if len(args)>0:
		if coord: 
			(ichr,ipos,wt,nw)=sys.argv[1].split(',')[:4]
			ochr=ichr
			if ichr.find('chr')==-1: ochr='chr'+ichr
			if ochr=='chrMT': ochr='chrM'
			ipos=int(ipos)
			if wt==nw:
				print >> sys.stderr, 'ERROR: Incorrect wild-type or mutant nucleotide',wt,nw
				sys.exit(1)
			if len(wt)==1 and len(nw)==1:
				pklcod=''
				pred_model=modfile[0]
			else:
				pred_model=modfile[1]
			make_prediction(ochr,ipos,wt,nw,pred_model,ucsc_exe,ucsc_dbs,win,fasta,dbpps,pklcod)
		else:	
			namefile=args[0]
			if not os.path.isfile(namefile):
				print >> sys.stderr,'ERROR: Input file not found',namefile
				sys.exit(1)
			if vcf:
				make_vcffile_multialleles_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win,fasta,dbpps,pklcod)
			else:
				make_tsvfile_predictions(namefile,modfile,ucsc_exe,ucsc_dbs,win,fasta,dbpps,pklcod)
	else:
		print 'predict_variants.py variant_file'
