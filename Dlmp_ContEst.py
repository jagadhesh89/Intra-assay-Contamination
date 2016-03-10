#!/home/m139105/Python-2.7.11/python
import sys
import re
import random
from decimal import *
import sys
import imp
import os
import subprocess
from optparse import OptionParser
from random import randint
from subprocess import call

def GetArgs():
	usage = "python Dlmp_ContEst.py -r <RunFolder> -c <ConfigFile> -o <OutputDir>"
	parser = OptionParser(usage = usage)
	parser.add_option("-c", "--reqConf", dest="config", help="The config file")
	parser.add_option("-r", "--runFolder", dest="run", help="The run folder")
	parser.add_option("-o", "--outFolder", dest="output", help="The output folder")
	(options, args) = parser.parse_args()
	if(len(sys.argv) < 7) :
		print("python Dlmp_ContEst.py -r <RunFolder> -c <ConfigFile> -o <OutputDir>")
		sys.exit()

	ConfigFile = os.path.abspath(options.config)
	OutputPath = os.path.abspath(options.output)
	RunName = os.path.abspath(options.run)+"/"
	SamplePath = RunName+"/samples/"
	if os.path.isdir(RunName+"samples/") == False:
		SamplePath = RunName
	SampleDir = os.listdir(SamplePath)
	if os.path.isdir(RunName+"logs") == False:
		Cmd = 'mkdir ' + RunName+"logs"
		subprocess.call([Cmd], shell=True)
	if os.path.isdir(OutputPath) == False:
		Cmd = 'mkdir ' + OutputPath
		subprocess.call([Cmd], shell=True)
	
	
	return (ConfigFile,SampleDir,SamplePath,OutputPath+"/")
	
def main():
	ConfigFile,SampleDir,SamplePath,OutputPath = GetArgs()
	Perl,SplitAllele,ContEst,Java,Cvg,Hapmap,Ref,Jvm_mem,Job_mem,Summary_mem,Qsub,Queue=ProcessConfigFile(ConfigFile)
	ProcessedVCFs,BAMNames = FindVCF_BAMFiles(SampleDir,SamplePath,SplitAllele,Perl,OutputPath)
	JobIds = ""
	for i in range(0,len(BAMNames)):
		if i == 0:
			JobIds = RunContEst(Java,ContEst,ProcessedVCFs[i],BAMNames[i],OutputPath,Hapmap,Ref,Cvg,i,Jvm_mem,Job_mem,Qsub,Queue)
		else:
			JobIds = JobIds + "," + RunContEst(Java,ContEst,ProcessedVCFs[i],BAMNames[i],OutputPath,Hapmap,Ref,Cvg,i,Jvm_mem,Job_mem,Qsub,Queue)
	PrepRunlevelFile(OutputPath,JobIds,Summary_mem,Qsub,Queue)
	return 1
	
def RunContEst(Java,ContEst,ProcessedVCF,BAMName,OutputPath,Hapmap,Ref,Cvg,Number,Jvm_mem,Job_mem,Qsub,Queue):
	SampleName = BAMName.split("/")[-1]
	if os.path.isdir(OutputPath+"tmp/Jobs") == False:
		Cmd = 'mkdir ' +OutputPath+"tmp/Jobs"
		subprocess.call([Cmd], shell=True)
	Cmd = Java + " -Xmx" + Jvm_mem + " -jar " + ContEst + " -T Contamination -et NO_ET -B:genotypes,vcf "+ProcessedVCF + " -BTI genotypes -B:pop,vcf "+Hapmap+ " -I " + BAMName + " -R " + Ref + " -mbc " + Cvg + " -o " +OutputPath+"tmp/"+SampleName.replace(".bam",".out")
	JobFile = open(OutputPath+"tmp/Jobs/"+SampleName.replace(".bam",".sh"),"w")
	JobFile.write("#!/bin/bash\n"+Cmd)
	JobFile.close()
	Cmd = Qsub +" -q " + Queue + " -l h_vmem="+ Job_mem +" -b y -N CONTEST_"+str(Number) + " -e " + OutputPath+"tmp/Jobs/"+SampleName.replace(".bam",".error") + " -o " + OutputPath+"tmp/Jobs/"+SampleName.replace(".bam",".out") + " " + OutputPath+"tmp/Jobs/"+SampleName.replace(".bam",".sh")
	subprocess.call([Cmd], shell=True)
	return "CONTEST_"+str(Number)

def ProcessVCF(VCFName,Perl,SplitAllele,SamplePath,SampleName,OutputPath):
	if os.path.isdir(OutputPath+"tmp") == False:
		Cmd = 'mkdir ' + OutputPath+"tmp"
		subprocess.call([Cmd], shell=True)
	Cmd = 'cat '+VCFName + " | " + Perl + " " + SplitAllele + " > " + OutputPath+"tmp/"+SampleName+"_MultiLine.vcf" 
	subprocess.call([Cmd], shell=True)
	Cmd = "sed -i 's/\.\//0\//g' " + OutputPath+"tmp/"+SampleName+"_MultiLine.vcf"
	subprocess.call([Cmd], shell=True)
	Cmd = "sed -i 's/0\/0/0/g' " + OutputPath+"tmp/"+SampleName+"_MultiLine.vcf"
	subprocess.call([Cmd], shell=True)
	Cmd = "sed -i '/chrX/Q' " + OutputPath+"tmp/"+SampleName+"_MultiLine.vcf"
	subprocess.call([Cmd], shell=True)
	return OutputPath+"tmp/"+SampleName+"_MultiLine.vcf"
	
def FindVCF_BAMFiles(SampleDir,SamplePath,SplitAllele,Perl,OutputPath):
	BAMNames = []
	ProcessedVCFs = []
	for eachSampleDir in SampleDir:
		if "lane" in eachSampleDir:
			continue
		eachSamplePath = SamplePath+eachSampleDir+"/"
		SampleName = eachSampleDir
		VCFname = eachSamplePath+SampleName+"_cmb.vcf"
		BAMNames.append(eachSamplePath+SampleName+".bam")
		ProcessedVCFs.append(ProcessVCF(VCFname,Perl,SplitAllele,eachSamplePath,SampleName,OutputPath))
	return (ProcessedVCFs,BAMNames)

def ProcessConfigFile(ConfigFile):
	ConfigFile = open(ConfigFile,"r")
	ConfigDict = {}
	for eachLine in ConfigFile:
		if "=" not in eachLine:
			continue
		eachLine_split = eachLine.split("=")
		ConfigDict[eachLine_split[0]] = eachLine_split[1].strip()
	ConfigFile.close()
	return (ConfigDict["PERL"],ConfigDict["MULTIALLELE"],ConfigDict["CONTEST"],ConfigDict["JAVA"],ConfigDict["CVG"],ConfigDict["HAPMAP_VCF"],ConfigDict["REFERENCE"],ConfigDict["JVM_MEMORY"],ConfigDict["JOB_MEMORY"],ConfigDict["SUMMARY_MEM"],ConfigDict["QSUB"],ConfigDict["QUEUE"])
	

def PrepRunlevelFile(OutputPath,JobIds,Summary_mem,Qsub,Queue):
	JobFolder = OutputPath+"tmp/Jobs"
	ContEstTmpFolder = OutputPath+"tmp/"
	Cmd = "for eachSample in "+ ContEstTmpFolder+"*.out; do  FileName=`basename $eachSample .out`;  ContEst=`tail -n 1 $eachSample | cut -f4`;  echo $FileName' '$ContEst >>"+ContEstTmpFolder+"../ContEstReport.txt; done"
	JobFile = open(OutputPath+"tmp/Jobs/ContEstReport.sh","w")
	JobFile.write("#!/bin/bash\n"+Cmd)
	JobFile.close()
	Cmd = Qsub+" -q "+ Queue+ " -l h_vmem="+Summary_mem +" -b y -N ContReport -hold_jid " + JobIds+  " -e " + OutputPath+"tmp/Jobs/ContestReport.error" + " -o " + OutputPath+"tmp/Jobs/ContEstReport.out" + " " + OutputPath+"tmp/Jobs/ContEstReport.sh"
	subprocess.call([Cmd], shell=True)

if __name__ == '__main__':
    main()







