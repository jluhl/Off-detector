# Off-detector

## Instruction
Here we propose our Off-Detector to facilitate the CRISPR-edited amplicon sequencing data analysis, with functions that existing tools are not able to provide. Off-Detector brings the following five key innovations: first, sample-specific reference sequences are calculated and used in mutant calling; second, functional annotations of CRISPR/Cas-induced variants are added, especially pathogenic variants; third, a batch mode is designed for comparing editing outcomes of multiple amplicons; fourth, a CRISPR/Cas-genome-editing-outcome database is provided, sharing editing outcomes of on-target sites and their corresponding off-target sites with the consent of data uploaders; fifth, processing memory and time are optimized allowing for hundreds of ampliconsâ€™ sequencing data.

## System requirements

## Installation Guide 
Download and install Anaconda from https://www.anaconda.com/products/individual  
Python and packages    
pip install pyfaidx==0.5.9  
pip install biopython  
pip install weblogo  
conda install -y bowtie2  

Install fastp  
wget http://opengene.org/fastp/fastp  
chmod a+x ./fastp  

Download and install Annovar from http://download.openbioinformatics.org/annovar_download_form.php  
perl annotate_variation.pl -downdb -webfrom annovar avdblist humandb/ -buildver hg38  
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/  
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20190305 humandb/  

Install blast

Install samtools

## Demo
Demo data could be downloaded in the directory Off-detector/data/ including treated group sequencing data (e.1.fq.gz and e.2.fq.gz , paired-end sequencing) and control group sequencing data (b.1.fq.gz and b.2.fq.gz, paired-end sequencing). Demo configure file could be downloaded in the directory Off-detector/amplicons.

## Usage
usage: Off-detector/scripts/pool_1125.py   
--experiment_read_1 --experiment_read_2 --background_read_1 --background_read_2   
--sequencing_type --configure_csv_path --pam_position --adapter_file_path   
--cleavage_window --minimum_percentage_support_snp --heterozygous_threshold --projectid   
--species --assembly --threads --output_path --genome_path --annovar_path --python_path  

### Arguments
--cleavage_window : We extend 5â€™ and 3â€™ of the â€?target site+PAMâ€™ region by this parameter.The extended region is called quantification window.Variants outside the quantification window will not be analyzed and included in our results.  
--minimum_percentage_support_snp : Minimum allele frequency for single allele (homozygous SNP site)  
--heterozygous_threshold : Minimum allele frequency for multiple alleles (heterozygous SNP site)  
