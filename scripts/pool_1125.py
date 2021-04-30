# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	File Name:	  compute_frequency
#	Description:	The script is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered
#					on-target efficiency.
#	Author:		 Lei Huang
#	Date:		   2020/11/18
#	E-mail:		 huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

import time
print('Start at:')
print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
import argparse
import os
import subprocess
import logging
import regex
import yaml
import textwrap
import multiprocessing
import pandas as pd
import glob
from Bio import pairwise2
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))

description = '''
------------------------------------------------------------------------------------------------------------------------ 

The script, supporting both paired-end and single-end reads, is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered on-target efficiency.

This is an example for input file with YAML format. Please note that the pcr sequence and guide RNA sequence must be one same strand.
# test.yaml
pcr_seq: TCTGTCTGAAACGGTCCCTGGCTAAACTCCACCCATGGGTTGGCCAGCCTTGCCTTGACCAATAGCCTTGACAAGGCAAACTTGACCAATAGTCTTAGAGTATCCAGTGAGGCCAGGGGCC
guide_rna:
	target: GCCAGCCTTGCCTTGACCAATAG
	PAM: TGGGTTG
adapter:
	seq_1: GCCAGCCTTGCCTTGACCAATAG
	seq_2: GCCAGCCTTGCCTTGACCAATAG
	seq_3: GCCAGCCTTGCCTTGACCAATAG
	seq_4: GCCAGCCTTGCCTTGACCAATAG
	
------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--threads", help="the number of CPU you want to use", default=3, required=True, type=int)
parse.add_argument("--background_read_1", help="the FASTQ file of read 1 of background for paired-end sequencing, default: NA", default="NA", required=False)
parse.add_argument("--background_read_2", help="the FASTQ file of read 2 of background for paired-end sequencing, default: NA", default="NA", required=False)
parse.add_argument("--background_read", help="the FASTQ file of background for single-end sequencing, default: NA", default="NA", required=False)
parse.add_argument("--experiment_read_1", help="the FASTQ file of read 1 of experiment for paired-end sequencing, default: NA", default="NA", required=False)
parse.add_argument("--experiment_read_2", help="the FASTQ file of read 2 of experiment for paired-end sequencing, default: NA", default="NA", required=False)
parse.add_argument("--experiment_read", help="the FASTQ file of experiment for single-end sequencing, default: NA", default="NA", required=False)
parse.add_argument("--output_path", help="the output path", required=True)
parse.add_argument("--adapter_file_path", help="the adapter file path", required=True)
parse.add_argument("--pam_position", help="the PAM position, 5' end or 3' end", choices=['5_end', '3_end'], required=True)
parse.add_argument("--sequencing_type", help="sequencing strategy, paired-end or single-end", choices=['pe', 'se'], required=True)
parse.add_argument("--cleavage_window", help="the number of flanking base pairs at both sides of gRNA without PAM, the mutations occured outside the region will be excluded and regarded as unmodified, default: 5", default=5, type=int)
#parse.add_argument("--minimum_alignment_score", help="the minimum alignment score for merged reads to amplicon, default: 100", default=100, type=int)
parse.add_argument("--genome_path", help="reference genome fasta's path.", required=False)
parse.add_argument("--assembly", help="genome assembly version.", required=False)
parse.add_argument("--annovar_path",help="path to annovar.", required=False)
parse.add_argument("--minimum_percentage_support_snp", help="a paticular mutation type (only for substitution), in background sample, with over #minimum_percentage_support_snp would be considered as a     inherent SNP. Note: you can set a big value to disable this function, like 101, default: 90", default=90, type=int)
parse.add_argument("--heterozygous_threshold", help="background sample merged aligned reads within cleavege window region, over #heterozygous_threshold would be considered as heterozygous", default=40, type=int)
parse.add_argument("--species",help="Species",required=False)
parse.add_argument("--projectid", help="Project ID for uploading to the gRNAdb",required=False)
parse.add_argument("--python_path", help="Path to python, default: python",required=False,default="python")
parse.add_argument("--configure_csv_path", help="path to thr configure csv file.", required=True)

args = parse.parse_args()

pam_position = args.pam_position

def rc(x):
	return x.upper()[::-1].replace('A','t').replace('T','A').replace('C','g').replace('G','C').replace('t','T').replace('g','G')

if not os.path.exists(args.output_path):
	os.mkdir(args.output_path)

sample = pd.read_csv(args.configure_csv_path,sep=',')
sample.columns =  ['Site type','pcr_seq','target','pam']

for i in sample['Site type'].values:
	pcr_seq = sample[sample['Site type'] == i]['pcr_seq'].values[0].upper().strip()
	target = sample[sample['Site type'] == i]['target'].values[0].upper().strip()
	pam = sample[sample['Site type'] == i]['pam'].values[0].upper().strip()
	if not os.path.exists(os.path.join(args.output_path,i)):
		os.mkdir(os.path.join(args.output_path,i))
	with open(os.path.join(args.output_path,i+'/'+i+'_configure.yaml'),'w') as f:
		if pam_position == '3_end':
			guide = target+pam
		elif pam_position == '5_end':
			guide = pam+target
		else:
			print("Please input pam_position as 3_end or 5_end")
		if guide in pcr_seq:
			f.write('pcr_seq: '+pcr_seq+'\nguide_rna:\n    target: '+target+'\n    PAM: '+pam+'\nadapter:\n')
		elif guide in rc(pcr_seq):
			f.write('pcr_seq: '+rc(pcr_seq)+'\nguide_rna:\n    target: '+target+'\n    PAM: '+pam+'\nadapter:\n')
		else:
			with open(os.path.join(args.output_path,'out_error.txt'),'w') as fe:
				fe.write("Guide RNA not on the PCR sequences on site "+i+'\n')
			sys.exit()				
	
		adapter_list = []
		with open(args.adapter_file_path,'r') as ad:
			for k in  ad.readlines():
				adapter_list.append(k.upper().strip())
			for ai in adapter_list:
				if rc(ai) not in adapter_list:
					adapter_list.append(rc(ai))
			ac = 0
			for ai in adapter_list:
				ac += 1
				f.write('    seq_'+str(ac)+': '+ai+'\n')

def pool_run(t):
	if args.sequencing_type == 'pe':
		if args.background_read_1 != None:
			if t == 0:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --background_read_1 '+args.background_read_1+' --background_read_2 '+args.background_read_2+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --background_read_1 '+args.background_read_1+' --background_read_2 '+args.background_read_2+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
			else:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(t))+' --background_read_1 '+args.background_read_1+' --background_read_2 '+args.background_read_2+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(t))+' --background_read_1 '+args.background_read_1+' --background_read_2 '+args.background_read_2+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
		else:
			if t == 0:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
			else:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(t))+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+ ' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(t))+' --experiment_read_1 '+args.experiment_read_1+' --experiment_read_2 '+args.experiment_read_2+' --sequencing_type pe --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
	elif args.sequencing_type == 'se':
		if args.background_read != None:
			if t == 0:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --background_read '+args.background_read+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --background_read '+args.background_read+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
			else:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(t))+' --background_read '+args.background_read+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(t))+' --background_read '+args.background_read+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
		else:
			if t == 0:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'on/on_configure.yaml')+' --output_path '+os.path.join(args.output_path,'on/')+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'on/nohup.log.txt')+' 2>&1')
			else:
				if args.assembly != None:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(i))+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' --annovar_path '+args.annovar_path+' --genome_path '+args.genome_path+' --assembly '+args.assembly+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
		
				else:
					os.system(args.python_path+' ' +  script_dir + '/scripts_1124.py --yaml '+os.path.join(args.output_path,'off'+str(t)+'/off'+str(t)+'_configure.yaml')+' --output_path '+os.path.join(args.output_path,'off'+str(i))+' --experiment_read '+args.experiment_read+' --sequencing_type se --pam_position '+pam_position+' --threads '+str(args.threads)+' --minimum_percentage_support_snp '+str(args.minimum_percentage_support_snp)+' --heterozygous_threshold '+str(args.heterozygous_threshold)+' --cleavage_window '+str(args.cleavage_window)+' > '+os.path.join(args.output_path,'off'+str(t)+'/nohup.log.txt')+' 2>&1')
	else:
		print('Please select sequencing type between se and pe\n')

if __name__ == "__main__":

	process = locals()
	for i in range(len(sample)):
		process['q'+str(i)] = multiprocessing.Process(target=pool_run, args=(i,))
	for i in range(len(sample)):
		process['q'+str(i)].start()
	for i in range(len(sample)):
		process['q'+str(i)].join()			

	def frameshift(x,y):
		if x.split(' ')[0] == 'frameshift':
			return y
		else:
			return 0

	def nonsynonymous(x,y):
		if x.split(' ')[0] == 'nonsynonymous':
			return y
		else:
			return 0

	def stopgain(x,y):
		if x == 'stopgain':
			return y
		else:
			return 0

	def startloss(x,y):
		if x == 'startloss':
			return y
		else:
			return 0

	def lpg(x,y):
		if x == 'Pathogenic' or (x == 'Likely_pathogenic') or (x == 'Pathogenic/Likely_pathogenic'):
			return y
		else:
			return 0
	
	
	os.chdir(args.output_path)
	if not os.path.exists('./comparison/'):
		os.mkdir('comparison')

	with open('comparison/temp.csv','w') as g:
		g.write('Site type,Clean reads,Mapped reads,Matched reads,Mutated reads,Mutation ratio,Indel ratio\n')
		for i in os.listdir('./'):
			s = ''
			if i[0] != '.' and (i[0] != 'c') and (os.path.isdir(i)):
				if os.path.exists(i+'/03.analysis/Statistics'):
					with open(i+'/03.analysis/Statistics','r') as f:
						for k in f.readlines():
							s += k.split('\t')[-1].strip()+','
						g.write(i+','+s[:-1]+'\n')
						crn = s.split(',')[0]
				else:
					with open(i+'/02.bowtie2/experiment/alignment/mapped.seq','r') as f:
						g.write(i+','+crn+(','+str(len(f.readlines())-1))*2+',0,0,0\n')

	c = pd.read_csv('comparison/temp.csv')
	#c = c[c['Mapped reads'] != 0]
	#c = c.dropna(subset=['Mapped reads'])

	def sort_sites(x):
		if x == 'on':
			return 0
		else:
			return int(x[3:])

	c['order'] = c['Site type'].apply(sort_sites)
	c = c.sort_values(by='order')
	c = c.drop('order',axis=1)	

	if args.assembly != None:

		def fs(x,y):
			if os.path.exists(x+'/04.outdir/out_multianno.xlsx'):
				df = pd.read_excel(x+'/04.outdir/out_multianno.xlsx')
				df['frameshift'] = df.apply(lambda row:frameshift(row['ExonicFunc.refGene'],row['Count']),axis=1)
				return df['frameshift'].sum()/y
			else:
				return 0

		def ns(x,y):
			if os.path.exists(x+'/04.outdir/out_multianno.xlsx'):
				df = pd.read_excel(x+'/04.outdir/out_multianno.xlsx')
				df['nonsynonymous'] = df.apply(lambda row:nonsynonymous(row['ExonicFunc.refGene'],row['Count']),axis=1)
				return df['nonsynonymous'].sum()/y
			else:
				return 0

		def sls(x,y):
			if os.path.exists(x+'/04.outdir/out_multianno.xlsx'):
				df = pd.read_excel(x+'/04.outdir/out_multianno.xlsx')
				df['startloss'] = df.apply(lambda row:startloss(row['ExonicFunc.refGene'],row['Count']),axis=1)
				return df['startloss'].sum()/y
			else:
				return 0

		def sgs(x,y):
			if os.path.exists(x+'/04.outdir/out_multianno.xlsx'):
				df = pd.read_excel(x+'/04.outdir/out_multianno.xlsx')
				df['stopgain'] = df.apply(lambda row:stopgain(row['ExonicFunc.refGene'],row['Count']),axis=1)
				return df['stopgain'].sum()/y
			else:
				return 0

		def pg(x,y):
			if os.path.exists(x+'/04.outdir/out_multianno.xlsx'):
				df = pd.read_excel(x+'/04.outdir/out_multianno.xlsx')
				df['Pathogenic/Likely_pathogenic ratio'] = df.apply(lambda row:lpg(row['CLNSIG'],row['Count']),axis=1)
				return df['Pathogenic/Likely_pathogenic ratio'].sum()/y
			else:
				return 0
		
		c['Frameshift ratio'] = c.apply(lambda row:fs(row['Site type'],row['Matched reads']),axis=1)
		c['Nonsynonymous ratio'] = c.apply(lambda row:ns(row['Site type'],row['Matched reads']),axis=1)
		c['Startloss ratio'] = c.apply(lambda row:sls(row['Site type'],row['Matched reads']),axis=1)
		c['Stopgain ratio'] = c.apply(lambda row:sgs(row['Site type'],row['Matched reads']),axis=1)
		
		if args.assembly != 'GRCm38.p6':
			c['Pathogenic/Likely_pathogenic ratio'] = c.apply(lambda row:pg(row['Site type'],row['Matched reads']),axis=1)
	c.to_excel('comparison/out_comparison.xlsx',index=None)

	with PdfPages('comparison/out_comparison_figures1.pdf') as pdf:
		sns.set_style("ticks")
		sns.set(style='white',color_codes=True,palette='RdBu',rc={"figure.figsize": (8, 6)})
		plt.figure(figsize=(4*len(os.listdir('.')), 6))
		Y = c['Mutation ratio']
		Y1 = c['Indel ratio']
		X = np.arange(len(Y))
		bar_width = 0.25
		tick_label = c['Site type']

		for x,y in zip(X,Y):
			plt.text(x+0.005,y+0.0005,'%.3f' %y, ha='center',va='bottom')
		for x,y1 in zip(X,Y1):
			plt.text(x+0.24,y1+0.0005,'%.3f' %y1, ha='center',va='bottom')

		plt.bar(X, Y, bar_width, align="center", color="aquamarine", label="Editing efficiency")
		plt.bar(X+bar_width, Y1, bar_width, color="dodgerblue", align="center", label="Indel ratio")
		plt.xlabel("Site type")
		plt.ylabel("Ratio")
		plt.title('Camparison')
		plt.xticks(X+bar_width/2, tick_label)
		plt.legend()
		pdf.savefig()
		plt.close()
	with open('comparison/out_comparison_table1.csv','w') as f:
		f.write('Site type,Mutation ratio,Indel ratio\n')
		for i in range(len(c)):	
			f.write(c['Site type'].values[i]+','+str(c['Mutation ratio'].values[i])+','+str(c['Indel ratio'].values[i])+'\n')

	if args.assembly != None:
		with PdfPages('comparison/out_comparison_figures2.pdf') as pdf:
			plt.figure(figsize=(4*len(os.listdir('.')), 6))
			Y = c['Nonsynonymous ratio']
			Y1 = c['Frameshift ratio']
			Y2 = c['Startloss ratio']
			Y3 = c['Stopgain ratio']
			X = np.arange(len(Y))
			bar_width = 0.25

			for x,y in zip(X,Y):
				plt.text(x+0.005,y+0.0005,'%.3f' %y, ha='center',va='bottom')
			for x,y1 in zip(X,Y1):
				plt.text(x+0.24,y1+0.0005,'%.3f' %y1, ha='center',va='bottom')
			for x,y2 in zip(X,Y2):
				plt.text(x+0.48,y2+0.0005,'%.3f' %y2, ha='center',va='bottom')
			for x,y3 in zip(X,Y3):
				plt.text(x+0.72,y3+0.0005,'%.3f' %y3, ha='center',va='bottom')

			plt.bar(X,Y,bar_width,align="center", color="aquamarine", label="Nonsynonymous ratio")
			plt.bar(X+bar_width, Y1,bar_width, color="dodgerblue", align="center", label="Frameshift ratio")
			plt.bar(X+2*bar_width,Y2,bar_width, color="lime",align="center", label="Startloss ratio")
			plt.bar(X+3*bar_width,Y3,bar_width, color="gold",align="center", label="Stopgain ratio")
			plt.xlabel("Site type")
			plt.ylabel("Ratio")
			plt.title('Camparison')
			plt.xticks(X+1.5*bar_width, tick_label)
			plt.legend()
			pdf.savefig()
			plt.close(i)
		with open('comparison/out_comparison_table2.csv','w') as f:
			f.write('Site type,Frameshift ratio,Nonsynonymous ratio,Startloss ratio,Stopgain ratio\n')
			for i in range(len(c)):
				f.write(c['Site type'].values[i]+','+str(c['Frameshift ratio'].values[i])+','+str(c['Nonsynonymous ratio'].values[i])+','+str(c['Startloss ratio'].values[i])+','+str(c['Stopgain ratio'].values[i])+'\n')
		
		if args.assembly != 'GRCm38.p6': 
			with PdfPages('comparison/out_comparison_figures3.pdf') as pdf:
				plt.figure(figsize=(4*len(os.listdir('.')), 6))
				Y = c['Pathogenic/Likely_pathogenic ratio']
				X = np.arange(len(Y))
				for x,y in zip(X,Y):
					plt.text(x+0.005,y+0.00005,'%.3f' %y, ha='center',va='bottom')

				plt.bar(X, Y, bar_width, align="center", color="dodgerblue", label="Pathogenic/Likely_pathogenic ratio")
				plt.xlabel("Site type")
				plt.ylabel("Ratio")
				plt.title('Camparison')
				plt.xticks(X,tick_label)
				plt.legend()
				pdf.savefig()
				plt.close()
			with open('comparison/out_comparison_table3.csv','w') as f:
				f.write('Site type,Pathogenic/Likely_pathogenic ratio\n')
				for i in range(len(c)):
					f.write(c['Site type'].values[i]+','+str(c['Pathogenic/Likely_pathogenic ratio'].values[i])+'\n')

	if args.projectid != None:
		def gene(x):
			if os.path.exists(x+'/04.outdir/out_multianno.xlsx'):
				df = pd.read_excel(x+'/04.outdir/out_multianno.xlsx')
			else:
				tmpx=glob.glob(os.path.join(x+'/04.outdir','anno*csv'))[0]
				df = pd.read_csv(tmpx)
			return df['Gene.refGene'].values[0]

		c['Gene'] = c['Site type'].apply(gene)
		c['Species'] = args.species
		c['Assembly'] = args.assembly

		def info_0(x):
			with open(x+'/'+x+'_configure.yaml','r') as f:
				for i in f.readlines():
					if i.strip().split(':')[0] == 'PAM':
						PAM = i.strip().split(':')[1].strip()
					elif i.strip().split(':')[0] == 'target':
						target = i.strip().split(':')[1].strip()
					elif i.strip().split(':')[0] == 'pcr_seq':
						pcr = i.strip().split(':')[1].strip()

				if args.pam_position == '3_end':
					gr = target+PAM
					c['PAM orientation'] = '3\' PAM'
				else:
					gr = PAM+target
					c['PAM orientation'] = '5\' PAM'

			with open(x+'/04.outdir/pcr_seq_out.txt','r') as f:
            			#chr9	minus	95108134	95107961
				line = f.readline()
				if line.split('\t')[1] == 'minus':

					start = int(line.split('\t')[-1].strip())+rc(pcr).find(rc(gr))
				else:
					start = int(line.split('\t')[2])+pcr.find(gr)
				return target+','+PAM+','+line.split('\t')[0]+':'+str(start)+'-'+str(start+len(gr)+1)

		c['info_0'] = c['Site type'].apply(info_0)
		c['Project'] = args.projectid
		c['Target sequence'] = c['info_0'].apply(lambda x:x.split(',')[0])
		c['PAM'] = c['info_0'].apply(lambda x:x.split(',')[1])
		c['Position'] = c['info_0'].apply(lambda x:x.split(',')[2])
		if args.assembly != 'GRCm38.p6':
			c = c[['Project','Site type','Species','Assembly','Gene','Position','Target sequence','PAM','PAM orientation','Mutation ratio', 'Indel ratio', 'Frameshift ratio','Nonsynonymous ratio', 'Startloss ratio', 'Stopgain ratio','Pathogenic/Likely_pathogenic ratio']]
		else:
			c = c[['Project','Site type','Species','Assembly','Gene','Position','Target sequence','PAM','PAM orientation','Mutation ratio', 'Indel ratio', 'Frameshift ratio','Nonsynonymous ratio', 'Startloss ratio', 'Stopgain ratio']]
		c['Date'] = time.strftime("%Y-%m-%d",time.localtime())

		c.to_excel('comparison/out_db.xlsx',index=None)
		print('End at:')
		print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))	
