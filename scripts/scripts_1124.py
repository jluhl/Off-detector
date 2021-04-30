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

import sys
import time
print('Start at:')
print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
import argparse
import itertools
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

def parase_yaml(yaml_file, pam_position, cleavage_window=5):
	with open(yaml_file, 'r') as f:
		items = yaml.load(f.read(), Loader=yaml.FullLoader)
		pcr_seq = items["pcr_seq"].upper()
		target_seq = items["guide_rna"]["target"].upper()
		pam_seq = items["guide_rna"]["PAM"].upper()
		# define the pam_position is the end position for 5'-end PAM, or the start position for 3'-end PAM.
		if pam_position == '5_end':
			guide_seq = pam_seq + target_seq	
			start_index = pcr_seq.find(guide_seq)
			if start_index == -1:
				logging.error("The guide RNA and PCR sequence are not on the same strand, please check!")
				os._exit(0)
			else:
				pam_end_position = start_index + len(pam_seq)
				end_index = start_index + len(guide_seq)
		else:
			guide_seq = target_seq + pam_seq
			start_index = pcr_seq.find(guide_seq)
			if start_index == -1:
				logging.error("The guide RNA and PCR sequence are not on the same strand, please check!")
				os._exit(0)
			else:
				pam_end_position = start_index + len(target_seq)
				end_index = start_index + len(guide_seq)
		adapter = items["adapter"]
		search_window_start = (pcr_seq.find(guide_seq) - cleavage_window) if (pcr_seq.find(guide_seq) - cleavage_window) >=0 else 0
		search_window_end = (pcr_seq.find(guide_seq) + len(guide_seq) + cleavage_window) if (pcr_seq.find(guide_seq) + len(guide_seq) + cleavage_window) <= len(pcr_seq) else len(pcr_seq)
	return pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter, search_window_start, search_window_end


def fastp(read_1, read_2, read, yaml_file, pam_position, sequencing_type, cleavage_window=5):
	pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter, search_window_start, search_window_end = parase_yaml(yaml_file, pam_position, cleavage_window)
	with open("adapter.fa", 'w') as f:
		for k, j in adapter.items():
			print(">", k, file=f)
			print(j, file=f)
	command_pe = "fastp -i " + read_1 + " -o unmerged_read_1.fq.gz -q 28 -u 10 -e 28 -I " + read_2 + " -O unmerged_read_2.fq.gz --thread 6 --merge " \
			"--merged_out merged_reads.fq.gz --adapter_fasta adapter.fa --length_required 50 > log.fastp 2>&1"
	command_se = "fastp -i " + read + " -o clean_read.fq.gz  -q 28 -u 10 -e 28 --thread 6 --adapter_fasta adapter.fa --length_required 50 > log.fastp 2>&1"
	if sequencing_type == "pe":
		subprocess.check_call(command_pe, shell=True)
	elif sequencing_type == "se":
		subprocess.check_call(command_se, shell=True)

def bowtie2(out_path, fastp_path, yaml_file, pam_position, sequencing_type, cleavage_window=5):
	pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter, search_window_start, search_window_end = parase_yaml(yaml_file, pam_position, cleavage_window)
	index_path = os.path.join(out_path, "index")
	alignment_path = os.path.join(out_path, "alignment")
	os.mkdir(alignment_path)
	os.mkdir(index_path)
	os.chdir(index_path)
	with open("pcr_product.fa", 'w') as f:
		print(">PCR_seq", file=f)
		print(pcr_seq, file=f)
	command_index = 'bowtie2-build pcr_product.fa index > log.index 2>&1'
	subprocess.check_call(command_index, shell=True)
	os.chdir(alignment_path)
	merged_reads_path = os.path.join(fastp_path, "merged_reads.fq.gz")
	single_reads_path = os.path.join(fastp_path, "clean_read.fq.gz")
	command_alignment_pe = "bowtie2 --sensitive-local --threads 8 -k 1 -x ../index/index -U " + merged_reads_path + " -S align.sam > log.bowtie2 2>&1"
	command_alignment_se = "bowtie2 --sensitive-local --threads 8 -k 1 -x ../index/index -U " + single_reads_path + " -S align.sam > log.bowtie2 2>&1"
	if sequencing_type == 'pe':
		subprocess.check_call(command_alignment_pe, shell=True)
	elif sequencing_type == 'se':
		subprocess.check_call(command_alignment_se, shell=True)



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
parse.add_argument("--yaml", help="the configured YAML file", required=True)
parse.add_argument("--pam_position", help="the PAM position, 5' end or 3' end", choices=['5_end', '3_end'], required=True)
parse.add_argument("--sequencing_type", help="sequencing strategy, paired-end or single-end", choices=['pe', 'se'], required=True)
parse.add_argument("--cleavage_window", help="the number of flanking base pairs at both sides of gRNA without PAM, the mutations occured outside the region will be excluded and regarded as unmodified, default: 5", default=5, type=int)
#parse.add_argument("--minimum_alignment_score", help="the minimum alignment score for merged reads to amplicon, default: 100", default=100, type=int)
parse.add_argument("--minimum_ref_length", help="a read with over minimum_ref_length is then considered as a qualified one to be aligned against reference sequence, default: 50", default=50, type=int )
parse.add_argument("--minimum_percentage_support_snp", help="a paticular mutation type (only for substitution), in background sample, with over #minimum_percentage_support_snp would be considered as a inherent SNP. Note: you can set a big value to disable this function, like 101, default: 90", default=90, type=int)
parse.add_argument("--genome_path", help="reference genome fasta's path.", required=False)
parse.add_argument("--assembly", help="genome assembly version.", required=False)
parse.add_argument("--annovar_path",help="path to annovar.", required=False)
parse.add_argument("--heterozygous_threshold", help="background sample merged aligned reads within cleavege window region, over #heterozygous_threshold would be considered as heterozygous", default=40, type=int)

args = parse.parse_args()

	# define the format of log
log_file_path = os.path.join(args.output_path, 'my.log')
LOG_FORMAT = "[%(asctime)s][%(levelname)s][%(module)s] %(message)s"
DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"
logging.basicConfig(filename=log_file_path, level=logging.DEBUG, format=LOG_FORMAT, datefmt=DATE_FORMAT)
output_path = os.path.abspath(args.output_path)
outdir_path = os.path.join(output_path, "04.outdir")
logging.info("Start at : "+time.strftime("%Y-%m-%d  %H:%M",time.localtime()))

if not os.path.exists(outdir_path):
	os.mkdir(outdir_path)

yaml_path = os.path.abspath(args.yaml)
yaml_file = yaml_path
pam_position = args.pam_position
cleavage_window = args.cleavage_window
#align_score_threshold = args.minimum_alignment_score
pcr_seq, guide_seq, start_index, end_index, pam_end_position, adapter, search_window_start, search_window_end = parase_yaml(yaml_file, pam_position, cleavage_window)
	################# experiment group only ################

with open(os.path.join(outdir_path, 'pcr_seq.fa'),'w') as f:
	f.write('>pcr_seq\n'+pcr_seq+'\n')

os.system('blastn -db '+args.genome_path+' -query '+os.path.join(outdir_path, 'pcr_seq.fa')+' -out '+os.path.join(outdir_path, 'pcr_seq_out.txt')+' -outfmt \"6 saccver sstrand sstart send\" -max_target_seqs 1')

with open(os.path.join(outdir_path,'pcr_seq_out.txt'),'r') as f:
	line = f.readline().split('\t')
	chr = line[0]
	sstrand = line[1]
	if sstrand == 'minus':
		sstart = line[3].strip()
	else:
		sstart = line[2]

if args.background_read_1 == "NA" and args.background_read == "NA":
	# run fastp
	logging.info("Step1: Running fastp...")
	fastp_path = os.path.join(output_path, "01.fastp")
	if not os.path.exists(fastp_path):
		os.mkdir(fastp_path)
	os.chdir(fastp_path)
	if (not os.path.exists(os.path.join(fastp_path, "merged_reads.fq.gz"))) and (not os.path.exists(os.path.join(fastp_path, "clean_read.fq.gz"))):
		if args.sequencing_type == 'pe':
			read_1_path = os.path.abspath(args.experiment_read_1)
			read_2_path = os.path.abspath(args.experiment_read_2)
			fastp(read_1_path, read_2_path, "NA", yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
		elif args.sequencing_type == 'se':
			read_path = os.path.abspath(args.experiment_read)
			fastp("NA", "NA", read_path, yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
	else:
		logging.warning("The mergered FASTQ has existed, skip fastp...")

		# run bowtie2
	logging.info("Step2: Running bowtie2...")
	bowtie2_path = os.path.join(output_path, "02.bowtie2")
	alignment_path = os.path.join(bowtie2_path, "alignment")
	sam_path = os.path.join(alignment_path, "align.sam")
	if not os.path.exists(bowtie2_path):
		os.mkdir(bowtie2_path)
	os.chdir(bowtie2_path)
	if not os.path.exists(sam_path):
		bowtie2(bowtie2_path, fastp_path, yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
	else:
		logging.warning("The SAM file has existed, skip bowtie2...")

		# run analysis
		#logging.info("Parasing SAM...")

		#parase_experiment_sam(output_path, yaml_path, args.pam_position, cleavage_window=args.cleavage_window, indels=args.indels, substitution=args.substitution, total_mutation_detail=None, snp_candidates=None, minimum_ref_length=args.minimum_ref_length)

		#logging.info("Finish analyszing...")

	################# background + experiment ################
elif (args.background_read_1 != "NA") or (args.background_read != "NA"):
		# run fastp
	logging.info("Step1: Running fastp...")
	fastp_path = os.path.join(output_path, "01.fastp")
	fastp_path_background = os.path.join(fastp_path, "background")
	fastp_path_experiment = os.path.join(fastp_path, "experiment")
	if not os.path.exists(fastp_path):
		os.mkdir(fastp_path)
	os.chdir(fastp_path)
		## background
	if not os.path.exists(fastp_path_background):
		os.mkdir(fastp_path_background)
	os.chdir(fastp_path_background)
	if (not os.path.exists(os.path.join(fastp_path_background, "merged_reads.fq.gz"))) and (not os.path.exists(os.path.join(fastp_path_background, "clean_read.fq.gz"))):
		if args.sequencing_type == 'pe':
			background_read_1_path = os.path.abspath(args.background_read_1)
			background_read_2_path = os.path.abspath(args.background_read_2)
			fastp(background_read_1_path, background_read_2_path, "NA", yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
		elif args.sequencing_type == 'se':
			background_read_path = args.background_read
			fastp("NA", "NA", background_read_path, yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
	else:
		logging.warning("The QC of background has been completed, skip...")

		## experiment
	if not os.path.exists(fastp_path_experiment):
		os.mkdir(fastp_path_experiment)
	os.chdir(fastp_path_experiment)
	if (not os.path.exists(os.path.join(fastp_path_experiment, "merged_reads.fq.gz"))) and (not os.path.exists(os.path.join(fastp_path_experiment, "clean_read.fq.gz"))):
		if args.sequencing_type == 'pe':
			experiment_read_1_path = os.path.abspath(args.experiment_read_1)
			experiment_read_2_path = os.path.abspath(args.experiment_read_2)
			fastp(experiment_read_1_path, experiment_read_2_path, "NA", yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
		elif args.sequencing_type == 'se':
			experiment_read_path = args.experiment_read
			fastp("NA", "NA", experiment_read_path, yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
	else:
		logging.warning("The QC of experiment has been completed, skip...")

		# run bowtie2
	logging.info("Step2: Running bowtie2...")
	bowtie2_path = os.path.join(output_path, "02.bowtie2")
	background_bowtie2_path = os.path.join(bowtie2_path, "background")
	experiment_bowtie2_path = os.path.join(bowtie2_path, "experiment")

	alignment_path = os.path.join(bowtie2_path, "alignment")
	background_sam_path = os.path.join(background_bowtie2_path, "alignment/align.sam")
	experiment_sam_path = os.path.join(experiment_bowtie2_path, "alignment/align.sam")
	if not os.path.exists(bowtie2_path):
		os.mkdir(bowtie2_path)
	os.chdir(bowtie2_path)
		## background
	if not os.path.exists(background_sam_path):
		if not os.path.exists(background_bowtie2_path):
			os.mkdir(background_bowtie2_path)
		os.chdir(background_bowtie2_path)
		bowtie2(background_bowtie2_path, fastp_path_background, yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
	else:
		logging.warning("The SAM file of background has existed, skip...")
	## experiment
	if not os.path.exists(experiment_sam_path):
		if not os.path.exists(experiment_bowtie2_path):
			os.mkdir(experiment_bowtie2_path)
		os.chdir(experiment_bowtie2_path)
		bowtie2(experiment_bowtie2_path, fastp_path_experiment, yaml_path, args.pam_position, args.sequencing_type, args.cleavage_window)
	else:
		logging.warning("The SAM file of background has existed, skip...")

		# run analysis
	logging.info("Step3: Parasing SAM...")

if 'experiment_bowtie2_path' not in dir():
	experiment_bowtie2_path = bowtie2_path

if 'experiment_sam_path' not in dir():
	experiment_sam_path = os.path.join(bowtie2_path, "alignment/align.sam")


os.system('grep -i PCR_seq ' + experiment_sam_path + ' | awk -F \'\t\' \'{print $10}\' > ' + os.path.join(experiment_bowtie2_path,'alignment/mapped.seq'))
#os.system('sed -i \'1d\' '+ os.path.join(experiment_sam_path,'alignment/mapped.seq'))

if 'background_bowtie2_path' in dir():
	logging.info('Step4.0: Analysising control group reads...')
	os.system('grep -i PCR_seq ' + background_sam_path + ' | awk -F \'\t\' \'{print $10}\' > ' + os.path.join(background_bowtie2_path,'alignment/mapped.seq'))
	#os.system('sed -i \'1d\' '+ os.path.join(background_sam_path + 'alignment/mapped.seq'))
	with open(os.path.join(background_bowtie2_path, "alignment/mapped.seq"),'r') as f:
		msf = f.readlines()
		if len(msf) == 1:
			logging.info("Error: No reads mapped to control group data on this site!")
			with open(os.path.join(outdir_path,'gene.annovar'),'w') as g:
				g.write(chr+'\t'+sstart+'\t'+sstart+'\tA\tT\thet\t.\t.\n')
			if args.assembly == 'GRCh37.p13':
				os.system('perl '+os.path.join(args.annovar_path,'table_annovar.pl')+' '+os.path.join(outdir_path, 'gene.annovar')+' '+os.path.join(args.annovar_path,'humandb/')+' -buildver hg19 -out '+os.path.join(outdir_path, 'anno')+' -remove -protocol refGene -operation g -nastring . -csvout')
			elif args.assembly == 'GRCh38.p13':
				os.system('perl '+os.path.join(args.annovar_path,'table_annovar.pl')+' '+os.path.join(outdir_path, 'gene.annovar')+' '+os.path.join(args.annovar_path,'humandb/')+' -buildver hg38 -out '+os.path.join(outdir_path, 'anno')+' -remove -protocol refGene -operation g -nastring . -csvout')
			elif args.assembly == 'GRCm38.p6':
				os.system('perl '+os.path.join(args.annovar_path,'table_annovar.pl')+' '+os.path.join(outdir_path, 'gene.annovar')+' '+os.path.join(args.annovar_path,'mousedb/')+' -buildver mm10 -out '+os.path.join(outdir_path, 'anno')+' -remove -protocol refGene -operation g -nastring . -csvout')
			sys.exit()		
		else:
			with open(os.path.join(background_bowtie2_path, "alignment/mapped.csv"),'w') as g:
				g.write('Seq')
				for i in msf:
					g.write(i.strip()+'\n')

	seqdf0 = pd.read_csv(os.path.join(background_bowtie2_path, "alignment/mapped.csv"),sep=',')
	uni0 =  seqdf0['Seq'].unique()
	len_uni0 = len(uni0)//args.threads

	def para0(x):
		with open(os.path.join(background_bowtie2_path, 'alignment/t_Unique_seq_100_'+str(x)+'.csv'),'w') as f:
			if x + 1 < args.threads:
				for i in uni0[x*len_uni0:(x+1)*len_uni0]:
					temp = bioalign(pcr_seq,i)
					if temp != '':
						f.write(i+','+str(len(seqdf0[seqdf0['Seq'] == i]))+','+temp[0]+','+temp[1]+','+str(temp[2])+'\n')
			else:
				for i in uni0[x*len_uni0:]:
					temp = bioalign(pcr_seq,i)
					if temp != '':
						f.write(i+','+str(len(seqdf0[seqdf0['Seq'] == i]))+','+temp[0]+','+temp[1]+','+str(temp[2])+'\n')

def rc(x):
	return x.upper()[::-1].replace('A','t').replace('T','A').replace('C','g').replace('G','C').replace('t','T').replace('g','G')

#align_score_threshold
ast = 0.8*len(pcr_seq)

def bioalign(x,y):
	a_ori = pairwise2.align.globalms(x,y,2,-1,-3,-1)[0]
	a_rc = pairwise2.align.globalms(x,rc(y),2,-1,-3,-1)[0]
	if a_ori[2] >= a_rc[2] and (a_ori[2] >= ast):
		return a_ori
	elif a_ori[2] < a_rc[2] and (a_rc[2] >= ast):
		return a_rc
	else:
		return ''

def mt_global(ref,alt):
	s_global = ''
	for i in range(len(ref)):
		if ref[i] == '-':
			s_global += 'I'
		elif ref[i] == alt[i]:
			s_global += 'M'
		elif alt[i] == '-':
			s_global += 'D'
		else:
			s_global += 'S'
	return s_global

    #Deletion size classification
    
def d_size(x):
	c = 0
	str_temp = ''
	if len(x) == x.count('D'):
		return str(len(x))
	else:
		for i in x:
			if i == 'D':
				c += 1
			else:
				if c != 0:
					str_temp += (str(c)+',')
				c = 0    
	return str_temp

# Insertion size classification
    
def i_size(x):
	c = 0
	str_temp = ''
	if len(x) == x.count('I'):
		return str(len(x))
	else:
		for i in x:
			if i == 'I':
				c += 1
			else:
				if c != 0:
					str_temp += (str(c)+',')
				c = 0   
	return str_temp

    
# Deletion position distribution
    
def d_pos(x):
	str_temp = ''    
	start_pos_y = 0    
	for i in x.replace('I',''):
		start_pos_y += 1
		if i == 'D':
			str_temp += str(start_pos_y)+','
	return str_temp

   
#Insertion position distribution
def imerge(x):
	while 'II' in x:
		x = x.replace('II','I')
	return x
def i_pos(x):
	str_temp = ''
	start_pos_y = 0
	for i in imerge(x):
		start_pos_y += 1
		if i == 'I':
			str_temp += str(start_pos_y)+','
	return str_temp

#Deletion position distribution
def s_pos(x):
	str_temp = ''
	start_pos_y = 0
	for i in x.replace('I',''):
		start_pos_y += 1
		if i == 'S':
			str_temp += str(start_pos_y)+','
	return str_temp

# grep -c '>' v2/02.bowtie2/experiment/alignment/align.fa
# grep -c PCR_seq align.sam >> 
with open(os.path.join(experiment_bowtie2_path, "alignment/mapped.seq"),'r') as f:
	with open(os.path.join(experiment_bowtie2_path, "alignment/mapped.csv"),'w') as g:
		g.write('Seq')
		for i in f.readlines():
			g.write(i.strip()+'\n')

seqdf = pd.read_csv(os.path.join(experiment_bowtie2_path, "alignment/mapped.csv"),sep=',')
analysis_path = os.path.join(output_path, "03.analysis")
if not os.path.exists(analysis_path):
	os.mkdir(analysis_path)
os.chdir(analysis_path)


uni = seqdf['Seq'].unique()
len_uni = len(uni)//args.threads

print('Amount of unique mapped sequence:')
print(len(uni))
	
def para(x):
	with open('t_Unique_seq_100_'+str(x)+'.csv','w') as f:
		if x + 1 < args.threads:
			for i in uni[x*len_uni:(x+1)*len_uni]:
				temp = bioalign(pcr_seq,i)       
				#if temp[2] >= align_score_threshold:
				if temp != '':
					f.write(i+','+str(len(seqdf[seqdf['Seq'] == i]))+','+temp[0]+','+temp[1]+','+str(temp[2])+'\n')
		else:
			for i in uni[x*len_uni:]:
				temp = bioalign(pcr_seq,i)
				#if temp[2] >= align_score_threshold:
				if temp != '':
					f.write(i+','+str(len(seqdf[seqdf['Seq'] == i]))+','+temp[0]+','+temp[1]+','+str(temp[2])+'\n')	
					#print('para run : ' + temp[0]+'\n')
	
if __name__ == "__main__":
        
	
	cseq = pcr_seq[pcr_seq.find(guide_seq) - args.cleavage_window : pcr_seq.find(guide_seq) + len(guide_seq) + args.cleavage_window]
	print(cseq)
	print(pcr_seq)
	print(guide_seq)

	def bioalign2(x,y):
		p_list = []
		for p0 in pairwise2.align.globalms(x,y,2,-1,-3,-2):		
      			p_list.append(len(p0[1].strip('-')))
		for p0 in pairwise2.align.globalms(x,y,2,-1,-3,-2):
			if len(p0[1].strip('-')) == min(p_list):
				return p0

	def window_seq0(r,a): 
		temp = bioalign2(r,cseq)[1]
		x1 = temp.find(cseq[0])
		x2 = len(temp) - temp[::-1].find(cseq[-1])

		if x1 != 0 and (x2 != len(temp)) and (a[:x1].replace('-','') != '') and  (a[x2:].replace('-','') != ''):
			# cseq = string_return
			string_return = bioalign2(a[x1:x2],cseq)[1]
			return string_return + ',' + a[x1:x2]

		elif x1 == 0 and (a[x2:].replace('-','') != ''):
			#print(temp[x1:x2] + ',' + a[x1:x2]+'\n')
			string_return = bioalign2(a[x1:x2],cseq)[1]
			return string_return + ',' + a[x1:x2]

		elif x2 == len(temp) and (a[:x1].replace('-','') != ''):
			string_return = bioalign2(a[x1:x2],cseq)[1]
			return string_return + ',' + a[x1:x2] 
		else:
			return ''	
	
	bdic = {('A',):'A',('C',):'C',('G',):'G',('T',):'T',('A','G'):'R',('C','T'):'Y',('A','C'):'M',('G','T'):'K',('C','G'):'S',('A','T'):'W',('A','C','T'):'H',('C','G','T'):'B',('A','C','G     '):'V',('A','G','T'):'D',('A','C','G','T'):'N'}
	ambig = {}
	for i in bdic.keys():
		ambig[bdic[i]] = i			
	
	if 'background_bowtie2_path' in dir():
		from weblogo import *		
		process0 = locals()
		for i in range(args.threads):
			process0['p0'+str(i)] = multiprocessing.Process(target=para0, args=(i,))
		for i in range(args.threads):
			process0['p0'+str(i)].start()
		for i in range(args.threads):
			process0['p0'+str(i)].join()

		csv_list = glob.glob(os.path.join(background_bowtie2_path, 'alignment/t*.csv'))
		for i in csv_list:
			fr = open(i,'rb').read()
			with open(os.path.join(background_bowtie2_path, 'alignment/Unique_seq.csv'),'ab') as f:
				f.write(fr)
		
		dfin0 = pd.read_csv(os.path.join(background_bowtie2_path, 'alignment/Unique_seq.csv'),sep=',',header=None)
		dfin0.columns = ['Seq','Count','Ref','Alt','Score']
		
		sudic = sum(dfin0['Count']) 
		if sudic // 10000 > 10:
			sudicx = sudic//10000
			dfin0['Count'] = dfin0['Count'].apply(lambda x:x//sudicx)
			dfin0 = dfin0[dfin0['Count'] != 0]
 
		dfin0['window_seq'] = dfin0.apply(lambda row:window_seq0(row['Ref'],row['Alt']),axis=1)
		dfin0 = dfin0[dfin0['window_seq'] != '']
		print(len(dfin0))
		dfin0['window_merged_reads'] = dfin0['window_seq'].apply(lambda x:x.split(',')[1])
		dfin0['window_cseq'] = dfin0['window_seq'].apply(lambda x:x.split(',')[0])
		#dfin0['window_merged_reads','Count'].to_csv(os.path.join(background_bowtie2_path, 'alignment/seqlogo_seq.csv'))
		# Count the cleavage area base pairs distribution of background aligned merged reads
		
		logging.info('Step4.1: weblogo...')
		dic_cw0 = {}
		for i in range(len(dfin0)):
			thr = dfin0[['window_merged_reads','Count','window_cseq']].values[i]
			c = -1
			if len(thr[2]) == len(thr[0]):
				for j in range(len(thr[2])):
					if thr[2][j] != '-':
						c += 1
						if thr[0][j] in ['A','G','C','T']:
							if c in dic_cw0.keys():
								if thr[0][j] not in dic_cw0[c].keys():
									dic_cw0[c][thr[0][j]] = thr[1]
								else:
									dic_cw0[c][thr[0][j]] += thr[1]
							else:
								dic_cw0[c] = {}
								dic_cw0[c][thr[0][j]] = thr[1]
	
		dic_cw = {}
		for k1 in sorted(dic_cw0.keys()):
			dic_cw[k1] = dic_cw0[k1]
		
		dic_cw_ratio = {}
		for i in dic_cw.keys():
			for j in dic_cw[i].keys():
				if j in ['A','G','C','T']:
					if i not in dic_cw_ratio.keys():
						dic_cw_ratio[i] = {}
						dic_cw_ratio[i][j] = dic_cw[i][j]/sum(dic_cw[i].values())
            
					else:
						dic_cw_ratio[i][j] = dic_cw[i][j]/sum(dic_cw[i].values())
		vdlist = []
		for i in dic_cw.keys():
			vdlist.append(sum(dic_cw[i].values()))

		with open(os.path.join(outdir_path,'sequence_logo.txt'),'w') as f: 
			tsum = max(vdlist)
			for j in dic_cw.keys():
				tstr = ''
				for t in dic_cw[j].keys():
					tstr += dic_cw[j][t]*(t+',')
				f.write((tstr+(tsum-sum(dic_cw[j].values()))*('-,'))[:-1]+'\n')	
		xlp = pd.read_csv(os.path.join(outdir_path,'sequence_logo.txt'),sep=',',header=None)
		with open(os.path.join(outdir_path,'sequence_logo_in.txt'),'w') as g:
			for i in range(len(xlp.columns)):
				str_v = ''
				for v in xlp[i].values:
					str_v += v
				g.write(str_v+'\n')
		f = open(os.path.join(outdir_path,'sequence_logo_in.txt'))
		seqs = read_seq_data(f)
		f.close()
		data = LogoData.from_seqs(seqs)
		options = LogoOptions()
		options.fineprint =  'Sequence logo'
		format = LogoFormat(data, options)
		eps = eps_formatter(data, format)
		with open(os.path.join(outdir_path,'out_weblogo.eps'),'wb') as f:
			f.write(eps)


		threshold_het = args.heterozygous_threshold/100
		threshold_wt = args.minimum_percentage_support_snp/100
		dic_ref = {}
		for i in dic_cw_ratio.keys():
			dic_ref[i] = []
		for i in dic_cw_ratio.keys():
			for j in dic_cw_ratio[i].keys():
				if dic_cw_ratio[i][j] >= threshold_het:
					dic_ref[i].append(j)
			for j in dic_cw_ratio[i].keys():
				if len(dic_ref[i]) == 0:
					dic_ref[i] = [cseq[i]]
				elif len(dic_ref[i]) == 1:
					if dic_cw_ratio[i][j] < threshold_wt and (dic_cw_ratio[i][j] >= threshold_het): 
						if dic_ref[i][0] != cseq[i]:
							dic_ref[i] = [cseq[i]]
							
				else:
					pass

			dic_ref[i] = bdic[tuple(sorted(dic_ref[i]))]
		
		with open(os.path.join(outdir_path,'out_ref_summary.txt'),'w') as fw:
			ncseq = ''
			for i in dic_ref.values():
				ncseq += i
			fw.write('According to the above figure, the sample-specific reference sequence in the quantification window (Ref<sub>ss</sub>) is:\n'+ncseq)
			if ncseq == cseq:
				fw.write(' (same as the conventional reference sequence).')
			else:
				for j in range(len(cseq)):
					if cseq[j] != ncseq[j]:		
						fw.write(" (with Position "+str(j+1)+" flipped from \""+cseq[j]+"\" to \""+ncseq[j]+"\" in the conventional reference sequence).")
				if cseq != ncseq:			
					cseq = ncseq
			fw.write('\nAll the mutations in this amplicon are called based on Ref<sub>ss</sub>.\n')
	
	logging.info("Step4.2 analysising treated group data...")
	process = locals()
	pcr_seq_ori = pcr_seq
	pcr_seq = pcr_seq.replace(pcr_seq[pcr_seq.find(guide_seq) - args.cleavage_window : pcr_seq.find(guide_seq) + len(guide_seq) + args.cleavage_window],cseq)	
	#print('before run : '+pcr_seq+'/n')
	groups = itertools.groupby(cseq, lambda char:char not in ambig)
	splits = []
	
	for b,group in groups:
		if b:
			splits.extend([[g] for g in group])
		else:
			for nuc in group:
				splits.append(ambig[nuc])
	answer = [''.join(p) for p in itertools.product(*splits)]

	for i in range(args.threads):
		process['p'+str(i)] = multiprocessing.Process(target=para, args=(i,))
	for i in range(args.threads):	
		process['p'+str(i)].start()
	for i in range(args.threads):
		process['p'+str(i)].join()

	csv_list = glob.glob('t*.csv')
	for i in csv_list:
		fr = open(i,'rb').read()  
		with open('Unique_seq.csv','ab') as f:
			f.write(fr) 

	dfin = pd.read_csv('Unique_seq.csv',header=None)
	dfin.columns = ['Seq','Count','Ref','Alt','Score']

	def bioalign2(x,y):
		p_list = []
		for p0 in pairwise2.align.globalms(x,y,2,-1,-3,-2):		
      			p_list.append(len(p0[1].strip('-')))
		for p0 in pairwise2.align.globalms(x,y,2,-1,-3,-2):
			if len(p0[1].strip('-')) == min(p_list):
				return p0

	def window_seq(r,a): 
		temp = bioalign2(r,cseq)[1]
		x1 = temp.find(cseq[0])
		x2 = len(temp) - temp[::-1].find(cseq[-1])

		max_score = []
		if x1 != 0 and (x2 != len(temp)) and (a[:x1].replace('-','') != '') and  (a[x2:].replace('-','') != ''):
			#print(temp[x1:x2] + ',' + a[x1:x2]+'\n')
			if len(answer) != 1:
				for best_cseq in answer:
					max_score.append(bioalign2(a[x1:x2],best_cseq)[2])
					#print(str(bioalign2(r,best_cseq)[2])+",")               
				max_score_max = max(max_score)
				for best_cseq in answer:             
					if  bioalign2(a[x1:x2],best_cseq)[2] == max_score_max:
						return bioalign2(a[x1:x2],best_cseq)[1] + ',' + a[x1:x2] + ',' + bioalign2(a[x1:x2],cseq)[1]
						break
			else:
				string_return = bioalign2(a[x1:x2],cseq)[1]
				return string_return + ',' + a[x1:x2] + ',' + string_return

		elif x1 == 0 and (a[x2:].replace('-','') != ''):
			#print(temp[x1:x2] + ',' + a[x1:x2]+'\n')
			if len(answer) != 1:
				for best_cseq in answer:
					max_score.append(bioalign2(a[x1:x2],best_cseq)[2])
				#print(str(bioalign2(r,best_cseq)[2])+",")               
				max_score_max = max(max_score)
				for best_cseq in answer:
				             
					if  bioalign2(a[x1:x2],best_cseq)[2] == max_score_max:
						return bioalign2(a[x1:x2],best_cseq)[1] + ',' + a[x1:x2] + ',' + bioalign2(a[x1:x2],cseq)[1]
						break   
			else:
				string_return = bioalign2(a[x1:x2],cseq)[1]
				return string_return + ',' + a[x1:x2] + ',' + string_return

		elif x2 == len(temp) and (a[:x1].replace('-','') != ''):
			if len(answer) != 1:
			#print(temp[x1:x2] + ',' + a[x1:x2]+'\n')
				for best_cseq in answer:
					max_score.append(bioalign2(a[x1:x2],best_cseq)[2])
				#print(str(bioalign2(r,best_cseq)[2])+",")               
				max_score_max = max(max_score)
				for best_cseq in answer:             
					if  bioalign2(a[x1:x2],best_cseq)[2] == max_score_max:
						return bioalign2(a[x1:x2],best_cseq)[1] + ',' + a[x1:x2] + ',' + bioalign2(a[x1:x2],cseq)[1]
						break               
			else:
				string_return = bioalign2(a[x1:x2],cseq)[1]
				return string_return + ',' + a[x1:x2] + ',' + string_return
		else:
			return ''	

	dfin['window_seq'] = dfin.apply(lambda row:window_seq(row['Ref'],row['Alt']),axis=1)
	dfin = dfin[dfin['window_seq'] != '']
	print(len(dfin))
	dfin['window_seq_ref'] = dfin['window_seq'].apply(lambda x:x.split(',')[0])
	dfin['window_seq_alt'] = dfin['window_seq'].apply(lambda x:x.split(',')[1])
	dfin['window_seq_ref_tab'] = dfin['window_seq'].apply(lambda x:x.split(',')[2])

	def aboutn(r,a):
		sfix = ''
		for i in range(len(a)):
			if a[i] == 'N':
				if r[i] != '-':
					sfix += r[i]
				else:
					sfix += 'N'
			else:
				sfix += a[i]
		return sfix

	dfin['window_seq_alt'] = dfin.apply(lambda row:aboutn(row['window_seq_ref'],row['window_seq_alt']),axis=1)


	with open('window_variants.csv','w') as f:
		c_same = 0
		
		f.write('Window_seq_ref_tab,Window_seq_ref,Window_seq_alt,Count\n')
		for i in dfin['window_seq_alt'].unique():
			if dfin[dfin['window_seq_alt']==i]['window_seq_ref'].values[0] != i:
	
				f.write(dfin[dfin['window_seq_alt']==i]['window_seq_ref_tab'].values[0]+','+dfin[dfin['window_seq_alt']==i]['window_seq_ref'].values[0]+','+i+','+str(dfin[dfin['window_seq_alt']==i]['Count'].sum())+'\n')	
			else:
				c_same += dfin[dfin['window_seq_alt']==i]['Count'].sum()	
	
	csv_in = pd.read_csv('window_variants.csv',sep=',')
	   
	csv_in['museq'] = csv_in.apply(lambda row:mt_global(row['Window_seq_ref'],row['Window_seq_alt']),axis=1)
	csv_in['d_size'] = csv_in['museq'].apply(d_size)
	csv_in['i_size'] = csv_in['museq'].apply(i_size)
	csv_in['s_pos'] = csv_in.apply(lambda row:s_pos(row['museq']),axis=1)
	csv_in['d_pos'] = csv_in.apply(lambda row:d_pos(row['museq']),axis=1)
	csv_in['i_pos'] = csv_in.apply(lambda row:i_pos(row['museq']),axis=1)
	csv_in['Ref_given'] = pcr_seq_ori[pcr_seq_ori.find(guide_seq) - args.cleavage_window : pcr_seq_ori.find(guide_seq) + len(guide_seq) + args.cleavage_window]

	#csv_in[['Ref_given','Window_seq_ref','Window_seq_alt','Count','museq','d_size','i_size','s_pos','d_pos','i_pos']]

	def indel_c(x):
		if 'I' in x or ('D' in x):
			return True
		else:
			return False
	csv_in['Indel'] = csv_in['museq'].apply(indel_c)
	
	with open(os.path.join(analysis_path,'Statistics'),'w') as fw:
		with open(experiment_sam_path,'r') as f:
			fw.write('Total _clean_reads:\t'+str(len(f.readlines())-3)+'\n')
		with open(os.path.join(experiment_bowtie2_path,'alignment/mapped.seq'),'r') as f:
			fw.write('mapped_reads:\t'+str(len(f.readlines())-1)+'\n')
		mrc = csv_in['Count'].sum()
		fw.write('match_reads:\t'+str(mrc+c_same)+'\n')
		fw.write('mutated_reads:\t'+str(mrc)+'\n')
		fw.write('mutation_rate [mutated_reads/match_reads]:\t'+str(mrc/(mrc+c_same))+'\n')
		fw.write('indel_rate :\t'+str(csv_in[csv_in['Indel'] == True]['Count'].sum()/(mrc+c_same))+'\n')

	cleavage_start = 1
	cleavage_end = len(cseq)
	d_dic = {}
	
	for i in csv_in[['Count','d_size']].values:
		if i[1] != '':
			for j in i[1].split(','):
				if j != '':
					if j not in d_dic.keys():
						d_dic[j] = i[0]
					else:
						d_dic[j] += i[0]
				else:
					pass

	i_dic = {}
	for i in csv_in[['Count','i_size']].values:
		if i[1] != '':
			for j in i[1].split(','):
				if j != '':
					if j not in i_dic.keys():
						i_dic[j] = i[0]
					else:
						i_dic[j] += i[0]
				else:
					pass
	sp_dic = {}
	for i in csv_in[['Count','s_pos']].values:
		if i[1] != '':
			for j in i[1].split(','):
				if j != '':
					if int(j) not in sp_dic.keys():
						sp_dic[int(j)] = i[0]
					else:
						sp_dic[int(j)] += i[0]
				else:
					pass

	scw_dic = {}
	for i in range(cleavage_start,cleavage_end+1):
		if i in sp_dic.keys():
			scw_dic[i] = sp_dic[i]
		else:
			scw_dic[i] = 0

	sp_keys = []
	sp_values = []
	ctemp0 = -1
	for i in sorted(scw_dic):
		ctemp0 += 1
		sp_keys.append(str(i)+' '+cseq[ctemp0])
		sp_values.append(scw_dic[i])



	dp_dic = {}
	for i in csv_in[['Count','d_pos']].values:
		if i[1] != '':
			for j in i[1].split(','):
				if j != '':
					if int(j) not in dp_dic.keys():
						dp_dic[int(j)] = i[0]
					else:
						dp_dic[int(j)] += i[0]
				else:
					pass

	dcw_dic = {}
	for i in range(cleavage_start,cleavage_end+1):
		if i in dp_dic.keys():
			dcw_dic[i] = dp_dic[i]
		else:
			dcw_dic[i] = 0

	dp_keys = []
	dp_values = []
	ctemp = -1

	for i in sorted(dcw_dic):
		ctemp += 1
		dp_keys.append(str(i)+' '+cseq[ctemp])
		dp_values.append(dcw_dic[i])


	ip_dic = {}
	for i in csv_in[['Count','i_pos']].values:
		if i[1] != '':
			for j in i[1].split(','):
				if j != '':
					if int(j) not in ip_dic.keys():
						ip_dic[int(j)] = i[0]
					else:
						ip_dic[int(j)] += i[0]
				else:
					pass
				
	print(ip_dic)
	icw_dic = {}
	for i in range(cleavage_start,cleavage_end+1):
		if i in ip_dic.keys():
			icw_dic[i] = ip_dic[i]
		else:
			icw_dic[i] = 0

	ip_keys = []
	ip_values = []
	ctemp2 = -1
	for i in sorted(icw_dic):
		ctemp2 += 1
		ip_keys.append(str(i)+' '+cseq[ctemp2])
		ip_values.append(icw_dic[i])

	logging.info('Step5: Plotting...')
	with open(os.path.join(outdir_path,'out_barplot_dataframe.csv'),'w') as f:
		f.write('pos,seq,subsitution,insertion,deletion')
		c = -1
		for i in range(cleavage_start,cleavage_end+1):
			c += 1
			f.write(str(i)+','+cseq[c]+','+str(scw_dic[i])+','+str(icw_dic[i])+','+str(dcw_dic[i])+'\n')

	sort_d_dic_values = []
	sort_d_dic_keys = []
	for i in d_dic.keys():
		sort_d_dic_keys.append(int(i))
	for i in sorted(sort_d_dic_keys):
		sort_d_dic_values.append(d_dic[str(i)])

	sort_i_dic_values = []
	sort_i_dic_keys = []
	for i in i_dic.keys():
		sort_i_dic_keys.append(int(i))
	for i in sorted(sort_i_dic_keys):
		sort_i_dic_values.append(i_dic[str(i)])

	sns.set_style("ticks")
	sns.set(style='whitegrid',color_codes=True,palette='RdBu',rc={"figure.figsize": (8, 6)})
	sns.set_context("paper",font_scale=2, rc={"lines.linewidth": 5})
	if len(i_dic) != 0:
		with PdfPages(os.path.join(outdir_path,'out_barplot_figures1.pdf')) as pdf:
			plt.figure(figsize=(len(icw_dic), 20))
			plt.title('Insertion Sizes')
			plt.xlabel('Size')
			plt.ylabel('Count')
			sns.barplot(sorted(sort_i_dic_keys),sort_i_dic_values)
			pdf.savefig()
			plt.close()
		with open(os.path.join(outdir_path,'out_barplot_table1.csv'),'w') as f:
			f.write('Insertion Size,Count\n')
			for i in range(len(sort_i_dic_keys)):
				f.write(str(sorted(sort_i_dic_keys)[i])+','+str(sort_i_dic_values[i])+'\n')
	else:
		print('No insertion detected')

	if len(d_dic) != 0:
		with PdfPages(os.path.join(outdir_path,'out_barplot_figures2.pdf')) as pdf:
			plt.figure(figsize=(len(icw_dic), 20))
			plt.title('Deletion Sizes')
			plt.xlabel('Size')
			plt.ylabel('Count')
			sns.barplot(sorted(sort_d_dic_keys),sort_d_dic_values)
			pdf.savefig()
			plt.close()
		with open(os.path.join(outdir_path,'out_barplot_table2.csv'),'w') as f:
			f.write('Deletion Size,Count\n')
			for i in range(len(sort_d_dic_keys)):
				f.write(str(sorted(sort_d_dic_keys)[i])+','+str(sort_d_dic_values[i])+'\n')
	else:
		print('No deletion detected')

	if len(sp_dic) != 0:
		with PdfPages(os.path.join(outdir_path,'out_barplot_figures3.pdf')) as pdf:
			plt.figure(figsize=(len(scw_dic), 20))
			plt.title('Substitution Position')
			plt.xlabel('Coordinate on the Chromosome')
			plt.ylabel('Count')
			plt.xticks(rotation=90)
			sns.barplot(sp_keys,sp_values)
			pdf.savefig()
			plt.close()
		with open(os.path.join(outdir_path,'out_barplot_table3.csv'),'w') as f:
			f.write('Substitution Position,Count\n')
			for i in range(len(sp_keys)):
				f.write(str(sp_keys[i])+','+str(sp_values[i])+'\n')
	else:
		print('No substitution detected')

	if len(i_dic) != 0:
		with PdfPages(os.path.join(outdir_path,'out_barplot_figures4.pdf')) as pdf:
			plt.figure(figsize=(len(icw_dic), 20))
			plt.title('Insertion Position')
			plt.xlabel('Coordinate on the Chromosome')
			plt.ylabel('Count')
			plt.xticks(rotation=90)
			sns.barplot(ip_keys,ip_values)
			pdf.savefig()
			plt.close()
		with open(os.path.join(outdir_path,'out_barplot_table4.csv'),'w') as f:
			f.write('Insertion Position,Count\n')
			for i in range(len(ip_keys)):
				f.write(str(ip_keys[i])+','+str(ip_values[i])+'\n')

	if len(d_dic) != 0:
		with PdfPages(os.path.join(outdir_path,'out_barplot_figures5.pdf')) as pdf:
			plt.figure(figsize=(len(dcw_dic), 20))
			plt.title('Deletion Position')
			plt.xlabel('Coordinate on the Chromosome')
			plt.ylabel('Count')
			plt.xticks(rotation=90)
			sns.barplot(dp_keys,dp_values)
			pdf.savefig()
			plt.close()
		with open(os.path.join(outdir_path,'out_barplot_table5.csv'),'w') as f:
			f.write('Deletion Position,Count\n')
			for i in range(len(dp_keys)):
				f.write(str(dp_keys[i])+','+str(dp_values[i])+'\n')

	with open(os.path.join(outdir_path,'out_barplot_dataframe.csv'),'w') as f:
		f.write('pos,seq,subsitution,insertion,deletion\n')
		c = -1
		for i in range(cleavage_start,cleavage_end+1):
			c += 1
			f.write(str(i)+','+cseq[c]+','+str(scw_dic[i])+','+str(icw_dic[i])+','+str(dcw_dic[i])+'\n')
	
	csv_in = csv_in[['Ref_given','Window_seq_ref_tab','Window_seq_ref','Window_seq_alt','Count','museq','Indel','d_size','i_size','s_pos','d_pos','i_pos']] 
	csv_in.columns = ['Conventional reference','Refwt','Real Ref','Window Alt','Count','Maker','Indel','Deletion Size','Insertion Size','Subsititution Position','Deletion Position','Insertion Position']
	def removedot(x):
		if x != '':
			if str(x)[-1] == ',':
				return str(x)[:-1]
			else:
				return x
		else:
			return ''
	
	csv_in['Insertion Size'] = csv_in['Insertion Size'].apply(removedot)
	csv_in['Deletion Size'] = csv_in['Deletion Size'].apply(removedot)
	csv_in = csv_in.sort_values(by="Count",ascending = False)
	csv_in.to_excel('window_variants.xlsx',index=None)
	csv_in[['Conventional reference','Refwt','Window Alt','Count','Indel','Deletion Size','Insertion Size']].to_excel('out_window_variants.xlsx',index=None)

	print('Cleavage window region unique mutation sequence:')
	print(c)

	print('End at:')
	print(time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
	if args.assembly != None:
		from pyfaidx import Fasta
		
		dfin = pd.read_excel('window_variants.xlsx',index_col=None)
		dfin.columns = ['Conventional reference','Refwt','Real Ref','Window Alt','Count','Maker','Indel','Deletion Size','Insertion Size','Subsititution Position','Deletion Position','Insertion Position']
		cseq = dfin['Real Ref'].values[0].replace('-','').strip()
		dfin['align_seq'] = dfin['Window Alt'].apply(lambda x:x.replace('-','').strip())
	
	
		## Extract DNA seq chr*:*-*
		genome_fasta = Fasta(args.genome_path)
		def genome_seq(x,y,z,m):
			if m == '-':
				return '-'
			else:
				return genome_fasta[x][y-1:z].seq
	
		def eqstr(x,y):
			if x == y:
				return 'T'
			else:
				return 'F'
	
		if sstrand == 'minus':
			print(pcr_seq)
			print(dfin['Conventional reference'].values[0])
			print(len(dfin['Conventional reference'].values[0]))	  #cleavage_start = int(sstart) - pcr_seq.find(dfin['rcw'].values[0].strip()) - len(dfin['reference_cleavage_window_with_SNP'].values[0].strip())
			cleavage_start = int(sstart) + rc(pcr_seq_ori).find(rc(dfin['Conventional reference'].values[0]))
		elif sstrand == 'plus':
			cleavage_start = pcr_seq_ori.find(dfin['Conventional reference'].values[0]) + int(sstart)
		else:
			print('Error : input sstrand')
	
		cleavage_end = cleavage_start+len(cseq)-1
		print(cleavage_end)
		print(cleavage_start)
		
		def end(x):
			for i in range(len(x)):
				if x[::-1][i] == 'M':
					pass
				else:
					return cleavage_end-i
		
		def start(x,y):
			return x-len(y.replace('-',''))+1

		def hx(x):
			if x.strip('-') == '':
				return '-'
			else:
				return x.replace('-','')

		def seq_g(x,y):
			for i in range(len(x)):
				if x[i] != y[i]:
					i1 = i
					break
			for i in range(len(x)):
				if x[::-1][i] != y[::-1][i]:
					i2 = i
					break
			return x[i1:len(x)-i2]+','+y[i1:len(x)-i2]

		if sstrand == 'minus':
			#csv_in[['Conventional reference','Refwt','Window Alt','Count','Maker','Indel','Deletion Size','Insertion Size']]
			csv_in['r2'] = csv_in['Refwt'].apply(rc)
			csv_in['a2'] = csv_in['Window Alt'].apply(rc)
			
			csv_in['m2'] = csv_in.apply(lambda row:mt_global(row['r2'],row['a2']),axis=1)
			csv_in['e'] = csv_in['m2'].apply(end)
			csv_in['seq_short'] = csv_in.apply(lambda row:seq_g(row['r2'],row['a2']),axis=1)
			csv_in['seq_short_r'] = csv_in['seq_short'].apply(lambda x:x.split(',')[0])
			csv_in['seq_short_a'] = csv_in['seq_short'].apply(lambda x:x.split(',')[1])
			csv_in['s'] = csv_in.apply(lambda row:start(row['e'],row['seq_short_r']),axis=1)

			#chr12	25245345	25245345	C	T	het	.	.
			with open(os.path.join(outdir_path, 'sample.annovar'),'w') as f:
				for i in range(len(csv_in)):
					sher = csv_in['s'].values[i]
					eher = csv_in['e'].values[i]
					ques = hx(csv_in['seq_short_r'].values[i])
					if ques != '-':
						f.write(chr+'\t'+str(sher)+'\t'+str(eher)+'\t'+hx(csv_in['seq_short_r'].values[i])+'\t'+hx(csv_in['seq_short_a'].values[i])+'\thet\t.\t.\n')
					else:
						f.write(chr+'\t'+str(eher)+'\t'+str(eher)+'\t'+hx(csv_in['seq_short_r'].values[i])+'\t'+hx(csv_in['seq_short_a'].values[i])+'\thet\t.\t.\n')
		else:
			csv_in['e'] = csv_in['Maker'].apply(end)
			#csv_in['e'] = csv_in['Maker'].apply(end)
			csv_in['seq_short'] = csv_in.apply(lambda row:seq_g(row['Refwt'],row['Window Alt']),axis=1)
			csv_in['seq_short_r'] = csv_in['seq_short'].apply(lambda x:x.split(',')[0])
			csv_in['seq_short_a'] = csv_in['seq_short'].apply(lambda x:x.split(',')[1])
			csv_in['s'] = csv_in.apply(lambda row:start(row['e'],row['seq_short_r']),axis=1)
			with open(os.path.join(outdir_path, 'sample.annovar'),'w') as f:
				for i in range(len(csv_in)):
					sher = csv_in['s'].values[i]
					eher = csv_in['e'].values[i]
					ques = hx(csv_in['seq_short_r'].values[i])
					if ques != '-':
						f.write(chr+'\t'+str(sher)+'\t'+str(eher)+'\t'+hx(csv_in['seq_short_r'].values[i])+'\t'+hx(csv_in['seq_short_a'].values[i])+'\thet\t.\t.\n')
					else:
						f.write(chr+'\t'+str(eher)+'\t'+str(eher)+'\t'+hx(csv_in['seq_short_r'].values[i])+'\t'+hx(csv_in['seq_short_a'].values[i])+'\thet\t.\t.\n')
						
		logging.info("Step6: running ANNOVAR...")
		if args.assembly == 'GRCh38.p13':
			os.system('perl '+os.path.join(args.annovar_path,'table_annovar.pl')+' '+os.path.join(outdir_path, 'sample.annovar')+' '+os.path.join(args.annovar_path,'humandb/')+' -buildver hg38 -out '+os.path.join(outdir_path, 'anno')+' -remove -protocol refGene,clinvar_20190305 -operation g,f -nastring . -csvout')
			csv_in = pd.read_csv(os.path.join(outdir_path,'anno.hg38_multianno.csv'),keep_default_na=False)
		elif args.assembly == 'GRCh37.p13':
			os.system('perl '+os.path.join(args.annovar_path,'table_annovar.pl')+' '+os.path.join(outdir_path, 'sample.annovar')+' '+os.path.join(args.annovar_path,'humandb/')+' -buildver hg19 -out '+os.path.join(outdir_path, 'anno')+' -remove -protocol refGene,clinvar_20190305 -operation g,f -nastring . -csvout')
			csv_in = pd.read_csv(os.path.join(outdir_path,'anno.hg19_multianno.csv'),keep_default_na=False)
		elif args.assembly == 'GRCm38.p6':
			os.system('perl '+os.path.join(args.annovar_path,'table_annovar.pl')+' '+os.path.join(outdir_path, 'sample.annovar')+' '+os.path.join(args.annovar_path,'mousedb/')+' -buildver mm10 -out '+os.path.join(outdir_path, 'anno')+' -remove -protocol refGene -operation g -nastring . -csvout')
			csv_in = pd.read_csv(os.path.join(outdir_path,'anno.mm10_multianno.csv'),keep_default_na=False)		
	   	#['Conventional reference','Real Ref','Alt','Count','Maker','Indel','Deletion Size','Insertion Size','Subsititution Position','Deletion Position','Insertion Position']
		csv_in = pd.merge(csv_in,dfin[['Conventional reference','Refwt','Window Alt','Count','Indel','Deletion Size','Insertion Size']],left_index=True,right_index=True)
		#csv_in = csv_in.join(temp[['Reads_Seqs','Reads_Count','Reads_Length','Mutation_Type']])
		csv_in['Alt'] = csv_in['Alt'].apply(lambda x:x.replace('het','N'))
		#ExonicFunc.refGene
		#frameshift insertion
		#Func.refGene
		def func_fill(x):
			if x == '.':
				for i in csv_in['Func.refGene'].unique():
					if i != '.':
						return i
			else:
				return x
		csv_in['Func.refGene'] = csv_in['Func.refGene'].apply(func_fill)
	
		def gene_fill(x):
			if x == '.':
				for i in csv_in['Gene.refGene'].unique():
					if i != '.':
						return i
			else:
				return x
		csv_in['Gene.refGene'] = csv_in['Gene.refGene'].apply(gene_fill)
	
		def frame_fill(x,y,z):
			if x == 'exonic':
				if y == 'N':
					return 'frameshift insertion'
				else:
					return z
			else:
				return z
		csv_in['ExonicFunc.refGene'] = csv_in.apply(lambda row:frame_fill(row['Func.refGene'],row['Alt'],row['ExonicFunc.refGene']),axis=1)
	
		if sstrand == 'minus':
			csv_in['Ref'] = csv_in['Ref'].apply(rc)
			csv_in['Alt'] = csv_in['Alt'].apply(rc)
			csv_in['Strand'] = '-'
			if csv_in['Conventional reference'].values[0] != csv_in['Refwt'].values[0]:
				csv_in['genome_ref'] = csv_in.apply(lambda row:genome_seq(row['Chr'],row['Start'],row['End'],row['Ref']),axis=1)
				csv_in['genome_ref'] = csv_in['genome_ref'].apply(rc)
				csv_in['Ref_unchaged'] = csv_in.apply(lambda row:eqstr(row['Ref'],row['genome_ref']),axis=1)
				if args.assembly != 'GRCm38.p6':
					csv_in = csv_in[['Conventional reference','Refwt','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','genome_ref','Ref_unchaged','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG']]
				else:
					csv_in = csv_in[['Conventional reference','Refwt','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','genome_ref','Ref_unchaged','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene']]
			else:
				if args.assembly != 'GRCm38.p6':
					csv_in = csv_in[['Conventional reference','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG']]
				else:
					csv_in = csv_in[['Conventional reference','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene']]			
			print(csv_in[csv_in['Alt'] != '0']['Count'].sum())
			csv_in = csv_in[csv_in['Alt'] != '0']

		else:
			csv_in['Strand'] = '+'
			if csv_in['Conventional reference'].values[0] != csv_in['Refwt'].values[0]:
				csv_in['genome_ref'] = csv_in.apply(lambda row:genome_seq(row['Chr'],row['Start'],row['End'],row['Ref']),axis=1)
				csv_in['Ref_unchaged'] = csv_in.apply(lambda row:eqstr(row['Ref'],row['genome_ref']),axis=1)
				if args.assembly != 'GRCm38.p6':
					csv_in = csv_in[['Conventional reference','Refwt','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','genome_ref','Ref_unchaged','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG']]
				else:
					csv_in = csv_in[['Conventional reference','Refwt','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','genome_ref','Ref_unchaged','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene']]
			else:
				if args.assembly != 'GRCm38.p6':
					csv_in = csv_in[['Conventional reference','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','CLNALLELEID','CLNDN','CLNDISDB','CLNREVSTAT','CLNSIG']]
				else:
					csv_in = csv_in[['Conventional reference','Window Alt','Count','Indel','Deletion Size','Insertion Size','Chr','Strand','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene']]

		csv_in.to_excel(os.path.join(outdir_path,'out_multianno.xlsx'),index=None)
		logging.info("Finish analysing...")
		logging.info('End at : '+time.strftime("%Y-%m-%d  %H:%M",time.localtime()))
		
