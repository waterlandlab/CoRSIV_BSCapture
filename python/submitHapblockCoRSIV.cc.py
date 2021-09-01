#!/usr/bin/env python
__author__ = "Cristian Coarfa"
__version__ = "1.0"

import os, sys, argparse, re, glob
import datetime
from argparse import RawTextHelpFormatter

from CCLabUtils.simpleTime import SimpleTime

class SubmitEMatrixQTL:
    DEBUG_PROGRESS = None
    DEBUG_CLUSTER_AVERAGE = True
    
    def __init__(self, myArgs):
        self.myArgs = myArgs
           
    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
            """\
Utility %s version %s
Utility that streamlines haplotype block/CoRSIV methylation analysis
            """%(sys.argv[0], __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-x','--corsivMethylationMatrix',   help='data matrix with beta values for CORSiVs',    required=True)
        parser.add_argument('-b','--corsivBedFile',         help='corsiv bed file',         required=True)
        parser.add_argument('-s','--snvFile',               help='GTEx phased SNV file',    required=True)
        parser.add_argument('-c','--chromosome',            help='single chromosome',       required=True)
        parser.add_argument('-S','--tissueMap',             help='map for samples/tissues/capture files',   required=True)
        parser.add_argument('-H','--haplotypeBlocks',       help='haplotype block plink.blocks.det file',   required=True)
        parser.add_argument('-d','--corsivBlockDistance',   help='maximum haplotype block/corsiv distance (bp)',   required=False, default=0)
        parser.add_argument('-o','--outputRoot',            help='output file root',        required=True)
        
        parser.add_argument('-C','--clusterType',       help='batch queuing type (pbs, slurm)',         required=True)
        parser.add_argument('-q','--queue',             help='cluster queue',                           required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def generateClusterJobPreamble(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitEMatrixQTL::generateClusterJobPreamble\n"%SimpleTime.now())
            
        self.myArgs.clusterType = self.myArgs.clusterType.upper()
        self.mappingFolder = os.path.dirname(os.path.abspath(self.myArgs.outputRoot))
        self.dName = os.path.basename(self.myArgs.outputRoot)
        
        os.system("mkdir -p "+self.mappingFolder + " ")

                        
        if (self.myArgs.clusterType == "PBS" ):
            self.clusterJobFile = "%s/hap.%s.%s.pbs"%(self.mappingFolder, self.dName, os.getpid())
            self.clusterJobFileWriter = open(self.clusterJobFile, "w")
            self.clusterJobFileWriter.write("#!/bin/bash\n")
            self.clusterJobFileWriter.write("#PBS -q "+self.myArgs.queue+"\n")
            self.clusterJobFileWriter.write("#PBS -l nodes=1:ppn=2\n")
            self.clusterJobFileWriter.write("#PBS -l vmem=16gb\n")
            self.clusterJobFileWriter.write("#PBS -l pmem=16gb\n")
            self.clusterJobFileWriter.write("#PBS -l walltime=48:00:00\n")
            self.clusterJobFileWriter.write("#PBS -l cput=48:00:00\n")
            # self.clusterJobFileWriter.write("#PBS -m ea\n")
            self.clusterJobFileWriter.write("#PBS -N "+self.dName+"."+str(os.getpid())+"\n")
            self.clusterJobFileWriter.write("#PBS -o "+self.mappingFolder+"/"+self.dName+".hap."+str(os.getpid())+".o\n")
            self.clusterJobFileWriter.write("#PBS -e "+self.mappingFolder+"/"+self.dName+".hap."+str(os.getpid())+".e\n")
            
        elif(self.myArgs.clusterType == "SLURM"):
            self.clusterJobFile = "%s/hap.%s.%s.slurm"%(self.mappingFolder, self.dName, os.getpid())
            self.clusterJobFile = self.mappingFolder +"/hap."+self.dName+"."+os.getpid()+".slurm"
            self.clusterJobFileWriter = open(self.clusterJobFile, "w")
            
            self.clusterJobFileWriter.write("#!/bin/bash\n")
            self.clusterJobFileWriter.write("#SBATCH -p "+self.myArgs.queue+"\n")
            self.clusterJobFileWriter.write("#SBATCH -N 1\n")
            self.clusterJobFileWriter.write("#SBATCH -n 2\n")
            self.clusterJobFileWriter.write("#SBATCH --mem 16gb \n")
            self.clusterJobFileWriter.write("#SBATCH -t 1-0:00 \n") # 1 day, 0 hours, 0 minutes
            self.clusterJobFileWriter.write("#SBATCH -o "+self.mappingFolder +"/" +self.dName + ".hap."+str(os.getpid())+".%N.%J.o\n")
            self.clusterJobFileWriter.write("#SBATCH -e "+self.mappingFolder +"/" +self.dName + ".hap."+str(os.getpid())+".%N.%J.e\n")
            self.clusterJobFileWriter.write("#SBATCH -J hap."+self.dName+"."+ str(os.getpid())+"\n")
            
        else:
            sys.stderr.write("Unknown cluster type !\n")
            sys.exit(4)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitEMatrixQTL::generateClusterJobPreamble\n"%SimpleTime.now())
            
    
    def generateJobCommands(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitEMatrixQTL::generateJobCommands\n"%SimpleTime.now())
            
        self.clusterJobFileWriter.write("hostname\n")
        self.clusterJobFileWriter.write("env\n")
        self.clusterJobFileWriter.write("ulimit -a\n")
           
        self.clusterJobFileWriter.write("source $HOME/.bashrc\n")
        self.clusterJobFileWriter.write("conda activate base\n")
        self.clusterJobFileWriter.write("module load R\n")
        
        tmpDir = "${TMPDIR}"
            
        self.clusterJobFileWriter.write("df "+ tmpDir +"  -h \n")
        self.fullScratch =tmpDir+ "/scratch.hap."+self.dName+"."+str(os.getpid())
        self.clusterJobFileWriter.write("mkdir -p "+self.fullScratch+"\n")
        self.clusterJobFileWriter.write("cd  "+self.fullScratch+"; sleep 1\n")

        input_list = [ os.path.abspath(self.myArgs.corsivMethylationMatrix),
                       os.path.abspath(self.myArgs.corsivBedFile),
                       os.path.abspath(self.myArgs.snvFile),
                       os.path.abspath(self.myArgs.tissueMap),
                       os.path.abspath(self.myArgs.haplotypeBlocks) ]
                       
                       
        self.clusterJobFileWriter.write("rsync -avL %s . \n"%" ".join(input_list))

        haplotype_block_vs_corsiv_command = "hapblockCoRSIV.cc.py -b %s -x %s -s %s -S %s -c %s -o %s -H %s &> log.%s.txt"%(os.path.basename(self.myArgs.corsivBedFile),
                os.path.basename(self.myArgs.corsivMethylationMatrix),
                os.path.basename(self.myArgs.snvFile),
                os.path.basename(self.myArgs.tissueMap),
                self.myArgs.chromosome,
                self.dName,
                os.path.basename(self.myArgs.haplotypeBlocks),
                self.dName)
                
        
        self.clusterJobFileWriter.write("%s\n"%haplotype_block_vs_corsiv_command)
        self.clusterJobFileWriter.write("sleep 1; rsync -av %s* log.* %s/\n"%(self.dName, self.mappingFolder))
        self.clusterJobFileWriter.write("cd  "+self.mappingFolder+"; sleep 1\n")
        
        self.clusterJobFileWriter.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitEMatrixQTL::generateJobCommands\n"%SimpleTime.now())
    
    def submitClusterJob(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitEMatrixQTL::submitClusterJob\n"%SimpleTime.now())
            
        if (self.myArgs.clusterType == "PBS" ):
            os.system("qsub "+self.clusterJobFile)
        elif(self.myArgs.clusterType == "SLURM"):
            os.system("sbatch "+self.clusterJobFile)
        else:
            sys.stderr.write("Unknown cluster type !\n")
            sys.exit(4)
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitEMatrixQTL::submitClusterJob\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitEMatrixQTL::work\n"%SimpleTime.now())
            
        self.generateClusterJobPreamble()
        self.generateJobCommands()
        self.submitClusterJob()  
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitEMatrixQTL::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line: %s\n"%" ".join(sys.argv))
        myArgs = SubmitEMatrixQTL.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = SubmitEMatrixQTL(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
