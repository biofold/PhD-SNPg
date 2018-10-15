#!/usr/bin/env python
import sys, os
import numpy as np
import __builtin__
global pdir
ptools=os.path.dirname(os.path.abspath(__file__))
sys.path=[ptools]+sys.path
pdir='/'.join(ptools.split('/')[:-1])
sys.path=[pdir]+sys.path
from sklearn.externals import joblib
from score_variants import global_vars, get_file_input, prediction



def global_vars(refgen):
	global ucsc_exe, ucsc_dbs, dbfasta, pklcod, dbpps, modfiles
	__builtin__.verbose=False
	ucsc_exe=pdir+'/ucsc/exe'
	ucsc_dbs=pdir+'/ucsc/'+refgen
	prog_dat=pdir+'/data/model'
	__builtin__.prog_dat=prog_dat
	dbfasta=refgen+'.2bit'
	pklcod=refgen+'_coding.pkl'
	if refgen=='hg38':                
		dbpps=['hg38.phyloP7way.bw','hg38.phyloP100way.bw']
	elif refgen=='hg19':
		dbpps=['hg19.phyloP46way.primate.bw','hg19.100way.phyloP100way.bw']
	else:
		print >> sys.stderr,'Incorrect reference genome. PhD-SNPg uses only hg19 and hg38.'
		sys.exit(1)
	modfiles=[ prog_dat+ '/snv_model_w5_p7_500_'+refgen+'.pkl',prog_dat + '/indel_model_w5_p7_500_'+refgen+'.pkl']


def get_input_data(namefile,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod):
	idata=[]
	sdata=[]
	input_list=get_file_input(namefile,ucsc_exe,ucsc_dbs,False,win,'\t',dbfasta,dbpps,pklcod)
	for i in input_list:
		if len(i[2])==0:
			print >> sys.stderr, 'ERROR: Incorrect variant. Genome location:',i[0][0],i[0][1]
			print >> sys.stderr, 'ERROR:','\t'.join(i[0])
			continue
		elif len(i[4])==0:
			print >> sys.stderr, 'ERROR: Conservation score not available. Genome location:',i[0][0],i[0][1]
			print >> sys.stderr, 'ERROR:','\t'.join(i[0])
			continue
		elif len(i[3])==0:
			i[3]=np.zeros(2*win+1).tolist()
		else:
			pass
		if i[5][0]==1 and i[5][1]==1:
			sdata.append(i)
		else:
			idata.append(i)
	return sdata,idata


def pred_sdata(data,model,ofile):
	out=[]
	model=joblib.load(model)
	X=[i[2]+i[3]+i[4] for i in data]
	p=prediction(X,model)
	for i in range(len(data)):
		if p[0][i]>0.5: 
			fdr=p[1][i][0]
		else:
			fdr=p[1][i][1]
		#out.append(data[i][0]+[p[2][i],p[0][i],fdr,data[i][5][2],sum(data[i][5])/5])
		out=[p[0][i],fdr,data[i][4][2],sum(data[i][4])/5]
		ofile.write('\t'.join(data[i][0]+[p[2][i]])+'\t'+'\t'.join(['%.3f' %round(j,3) for j in out])+'\n')
	return 


def pred_idata(data,model,ofile):
	out=[]
	model=joblib.load(model)
	X=[i[2]+i[3]+i[4]+i[5] for i in data]
	p=prediction(X,model)
	for i in range(len(data)):
		if p[0][i]>0.5:
			fdr=p[1][i][2]
		else:
			fdr=p[1][i][3]
		#out.append(data[i][0]+[p[2][i],p[0][i],fdr,data[i][5][2],sum(data[i][5])/5])
		out=[p[0][i],fdr,data[i][4][2],sum(data[i][4])/5]
		ofile.write('\t'.join(data[i][0]+[p[2][i]])+'\t'+'\t'.join(['%.3f' %round(j,3) for j in out])+'\n')
	return


def get_options():
        import optparse        
	desc = 'Script for running scoring pipeline'        
	parser = optparse.OptionParser("usage: [-h] [-t var_type] [-o outfile]", description=desc)
	parser.add_option('-o','--output', action='store', type='string', dest='outfile', help='Output file')        
	parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='hg38', help='Genome version')
	(options, args) = parser.parse_args()
	win=2
	hg='hg38'
	outfile=''
	if options.hg:	hg=options.hg.lower()
	if options.outfile: outfile=options.outfile
	return args, hg, outfile, win
	


if __name__  == '__main__':
	args, hg, outfile, win=get_options()
	namefile=args[0]
	refgen=hg
	global_vars(refgen)
	if outfile=='':
		fout=sys.stdout
	else:
		fout=open(outfile,'w')
	if dbpps!=[]:
		sdata,idata=get_input_data(namefile,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod)
		if len(idata)==0 and len(sdata)==0:
			print >> sys.stderr,'WARNING: No mutation data found. Please check your input.'
			sys.exit(1)
		if len(sdata)>0: ps=pred_sdata(sdata,modfiles[0],fout)
		if len(idata)>0: pi=pred_idata(idata,modfiles[1],fout)	
	fout.close()

