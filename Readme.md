#### Detecting Intra Assay contamination using Contest
```
############################################################################################################
##	Intra assay contamination occurs when a patient sample designated for a test orderable is contaminated with patient(s) who came in for the same test orderable. This results in erroneous reports being generated due to variant data in the patient sample being misrepresented. Detecting this is crucial in improving patient care.
##
##	This script has been designed in order to detect these intra assay contamination events. This script uses a tool developed by Broad institute named ContEst. The script runs ContEst on the CLC samples and detects any abnormality in the patient data.
##
##	-------------------------
##	USAGE:
##	-------------------------
##
##	python Dlmp_ContEst.py -r <RunFolder> -c <ConfigFile> -o <OutputDir>
##                 -r CLC_Run_Folder # REQUIRED
##                 -c configfile     # REQUIRED Contest config File
##                -o outdir         # REQUIRED outdir name
##               
##	---------------------------
##	INPUT DATA STRUCTURE:
##	---------------------------
##
##   <RUN_NAME>/samples/<sampleA>
##          /<sampleB> 
##
##   under each sampleX, it will looks for these files
##      	<sampleX>.bam
##		<sampleX>_cmb.vcf
##
##	----------------------------
##	CONFIG FILE:
##	----------------------------
##	A config file is needed. 
##
##	#Scripts and tools
##	PERL=/usr/local/biotools/perl/5.10.1/bin/perl
##	MULTIALLELE=/dlmp/sandbox/cgslIS/codes/CONTAMINATION/MultiAllele_VCFsplit.pl
##	CONTEST=/dlmp/sandbox/cgslIS/codes/CONTAMINATION/contest-1.0.24530-bin/ContEst.jar
##	JAVA=/usr/local/biotools/java/jdk1.6.0_33/bin/java
##	QSUB=/usr/local/biotools/oge/ge2011.11/bin/linux-x64/qsub
##
##	#Variant Min Cvg
##	CVG=10
##	
##	#Contest Reference
##	HAPMAP_VCF=/dlmp/sandbox/cgslIS/db/HapMap_hg19.vcf
##	REFERENCE=/dlmp/sandbox/cgslIS/db/hg19_ordered.fa
##
##	#Memory Specs
##	JVM_MEMORY=5G
##	JOB_MEMORY=10G
##	SUMMARY_MEM=1G
##	QUEUE=sandbox.q  
##
##	-----------------------------------
##	DEPENDENCIES
##	-----------------------------------
##	TOOLS
##		1. Python
##		2. Perl
##		3. Java
##		4. ContEst - Located in "contest-1.0.24530-bin" folder in Git repo
##	SCRIPTS
##		1. MultiAllele_VCFsplit.pl
##	
##	-----------------------------------
##	MAIN SUB STEPS:
##	-----------------------------------
##
##	1.	The program processes the config files for the tools and memory for submitting the jobs.
##	2.	It then scans each sample directory in the run folder for the BAM and VCF File
##	3.	The program processes each VCF to correct for following
##		Multiple alleles are split to separate lines 
##			SCRIPT - /dlmp/sandbox/cgslIS/codes/CONTAMINATION/MultiAllele_VCFsplit.pl
##		Remove all Chromosome X and Y variants
##		Substitute genotype info with 0 when “.” Is present, otherwise ContEst does not work.
##	4.	Create Job Files to submit ContEst tool for each sample
##			SCRIPT - /dlmp/sandbox/cgslIS/codes/CONTAMINATION/contest-1.0.24530-bin/ContEst.jar
##	5.	Submit the ContEst job for each sample and submit one job for the overall report. 
##	6.	Process output from ContEst and create a run level summary file
##
##	-------------------------------------
##	RUN EXAMPLE:
##	-------------------------------------
##
##	python /dlmp/sandbox/cgslIS/codes/CONTAMINATION/Dlmp_ContEst.py -r /dlmp/prod/runs/NGSHM/NGSHM_20160406_SSXT140_MS274A_AN0WH/ -c /dlmp/sandbox/cgslIS/codes/CONTAMINATION/config.txt -o /dlmp/sandbox/cgslIS/Jag/CONTAMINATION/SQA/ContEst
##
##
##	Your job 3632653 ("CONTEST_NGSHM_20160406_SSXT140_MS274A_AN0WH_0") has been submitted
##	Your job 3632654 ("CONTEST_NGSHM_20160406_SSXT140_MS274A_AN0WH_1") has been submitted
##	Your job 3632655 ("CONTEST_NGSHM_20160406_SSXT140_MS274A_AN0WH_2") has been submitted
##	Your job 3632656 ("CONTEST_NGSHM_20160406_SSXT140_MS274A_AN0WH_3") has been submitted
##	Your job 3632657 ("CONTEST_NGSHM_20160406_SSXT140_MS274A_AN0WH_4") has been submitted
##	Your job 3632658 ("CONTEST_NGSHM_20160406_SSXT140_MS274A_AN0WH_5") has been submitted
##	Your job 3632659 ("ContReport") has been submitted
##
##	-----------------------------------
##	Output file example
##	-----------------------------------
##	When jobs finish running, there will be text file report generated in the output folder named “ContEstReport.txt”:
##
##	LOCATION: /dlmp/sandbox/cgslIS/Jag/CONTAMINATION/SQA/ContEst/ContEstReport.txt
##
##	10068475373 4.2
##	10068502986 3.6
##	10068516876 1.4
##	10068518196 0.5
##	10068518568 1.4
##	
##	Each line has a sample name and then a % of contamination estimate for the sample.
```