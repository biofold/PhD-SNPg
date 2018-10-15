#!/usr/bin/env python
import sys, os
import __builtin__
pdir=os.path.dirname(os.path.abspath(__file__))
sys.path.append(pdir)
sys.path=[pdir+'/tools']+sys.path
from sklearn.externals import joblib
from score_variants import global_vars, get_file_input, prediction


def get_input_data(namefile,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,vcf=False):
	idata=[]
	sdata=[]
	input_list=get_file_input(namefile,ucsc_exe,ucsc_dbs,False,win,'\t',dbfasta,dbpps,pklcod,vcf=vcf)
	for i in input_list:
		if i[5][0]==1 and i[5][1]==1:
			sdata.append(i)
		else:
			idata.append(i)
	return sdata,idata


def pred_sdata(data,model):
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
		print '\t'.join(data[i][0]+[p[2][i]])+'\t'+'\t'.join(['%.3f' %round(j,3) for j in out])
	return 

def pred_idata(data,model):
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
		print '\t'.join(data[i][0]+[p[2][i]])+'\t'+'\t'.join(['%.3f' %round(j,3) for j in out])
	return


def get_options():
        import optparse        
	desc = 'Script for running scoring pipeline'        
	parser = optparse.OptionParser("usage: [-h] [-t var_type] [-o outfile]", description=desc)
	#parser.add_option('-o','--output', action='store', type='string', dest='outfile', help='Output file')        
	parser.add_option('-g','--genome', action='store', type='string', dest='hg', default='hg38', help='Genome version')
	parser.add_option('-v','--vcf', action='store_true',dest='vcf',default=False, help='VCF input file')
	parser.add_option('--verbose', action='store_true', dest='ver', default=False, help='Verbose mode')
	(options, args) = parser.parse_args()
	win=2
	hg='hg38'
	vcf=False
	if options.hg:	hg=options.hg.lower()
	if options.vcf: vcf=True
	if options.ver: __builtin__.verbose=True
	if len(args)==0: 
		parser.print_help()
		sys,exit(1)
	return args, hg, win,vcf
	


if __name__  == '__main__':
	__builtin__.verbose=False
	args, hg, win, vcf=get_options()
	namefile=args[0]
	refgen=hg
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
		dbpps=[]
	modfiles=[ prog_dat+ '/snv_model_w5_p7_500_'+refgen+'.pkl',prog_dat + '/indel_model_w5_p7_500_'+refgen+'.pkl']
	if dbpps!=[]:
		sdata,idata=get_input_data(namefile,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,vcf)
		sdata,idata=get_input_data(namefile,ucsc_exe,ucsc_dbs,win,dbfasta,dbpps,pklcod,vcf)
		if len(idata)==0 and len(sdata)==0:
			print sys.stderr(),'WARNING: No mutation data found. Please check your input.'
			sys.exit(1)
		if len(sdata)>0: ps=pred_sdata(sdata,modfiles[0])
		if len(idata)>0: pi=pred_idata(idata,modfiles[1])	
		#input_list=get_file_input(namefile,ucsc_exe,ucsc_dbs,False,2,'\t',dbfasta,dbpps,pklcod)
		#print input_list[:9]

