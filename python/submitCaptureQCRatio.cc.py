#!/usr/bin/env python
__author__ = "Cristian Coarfa"
__version__ = "1.0"

import os, sys, argparse, re, glob
import datetime
from argparse import RawTextHelpFormatter

from CCLabUtils.simpleTime import SimpleTime

class SubmitCaptureQC:
    DEBUG_PROGRESS = None
    DEBUG_CLUSTER_AVERAGE = True
    
    def __init__(self, myArgs):
        self.myArgs = myArgs
           
    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
            """\
Utility %s version %s
Utility that streamlines on-target reads computation
TODO: play with tmp dir location; for now assume ${TMPDIR}
            """%(sys.argv[0], __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-b','--bamFilePattern',    help='BAM file(s) pattern',                                 required=True)
        parser.add_argument('-C','--captureTarget',     help='BED file with capture targets',                       required=True)
        parser.add_argument('-r','--ratio',             help='read ratio overlap w/ capture targets (0,1]',         required=True)
        parser.add_argument('-e','--experiment',        help='experiment label',                                    required=True)
        parser.add_argument('-o','--outputRoot',        help='output file root (includes folder name if needed)',   required=True)
        
        parser.add_argument('-x','--tmpDir',            help='temporary directory',                     required=False)
        parser.add_argument('-c','--clusterType',       help='batch queuing type (pbs, slurm)',         required=True)
        parser.add_argument('-q','--queue',             help='cluster queue',                           required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def generateClusterJobPreamble(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitCaptureQC::generateClusterJobPreamble\n"%SimpleTime.now())
            
        self.myArgs.clusterType = self.myArgs.clusterType.upper()
        self.mappingFolder = os.path.dirname(os.path.abspath(self.myArgs.outputRoot))
        self.dName = os.path.basename(self.myArgs.outputRoot)
        
        os.system("mkdir -p "+self.mappingFolder + " ")

                        
        if (self.myArgs.clusterType == "PBS" ):
            self.clusterJobFile = "%s/captureQCRatio.%s.%s.pbs"%(self.mappingFolder, self.dName, os.getpid())
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
            self.clusterJobFileWriter.write("#PBS -o "+self.mappingFolder+"/"+self.dName+".captureQCRatio."+str(os.getpid())+".o\n")
            self.clusterJobFileWriter.write("#PBS -e "+self.mappingFolder+"/"+self.dName+".captureQCRatio."+str(os.getpid())+".e\n")
            
        elif(self.myArgs.clusterType == "SLURM"):
            self.clusterJobFile = "%s/captureQCRatio.%s.%s.slurm"%(self.mappingFolder, self.dName, os.getpid())
            self.clusterJobFile = self.mappingFolder +"/captureQCRatio."+self.dName+"."+os.getpid()+".slurm"
            self.clusterJobFileWriter = open(self.clusterJobFile, "w")
            
            self.clusterJobFileWriter.write("#!/bin/bash\n")
            self.clusterJobFileWriter.write("#SBATCH -p "+self.myArgs.queue+"\n")
            self.clusterJobFileWriter.write("#SBATCH -N 1\n")
            self.clusterJobFileWriter.write("#SBATCH -n 2\n")
            self.clusterJobFileWriter.write("#SBATCH --mem 16gb \n")
            self.clusterJobFileWriter.write("#SBATCH -t 1-0:00 \n") # 1 day, 0 hours, 0 minutes
            self.clusterJobFileWriter.write("#SBATCH -o "+self.mappingFolder +"/" +self.dName + ".captureQCRatio."+str(os.getpid())+".%N.%J.o\n")
            self.clusterJobFileWriter.write("#SBATCH -e "+self.mappingFolder +"/" +self.dName + ".captureQCRatio."+str(os.getpid())+".%N.%J.e\n")
            self.clusterJobFileWriter.write("#SBATCH -J captureQCRatio"+self.dName+"."+ str(os.getpid())+"\n")
            
        else:
            sys.stderr.write("Unknown cluster type !\n")
            sys.exit(4)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitCaptureQC::generateClusterJobPreamble\n"%SimpleTime.now())
            
    
    def generateJobCommands(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitCaptureQC::generateJobCommands\n"%SimpleTime.now())
            
        self.clusterJobFileWriter.write("hostname\n")
        self.clusterJobFileWriter.write("env\n")
        self.clusterJobFileWriter.write("ulimit -a\n")
           
        self.clusterJobFileWriter.write("source $HOME/.bashrc\n")
        self.clusterJobFileWriter.write("conda activate base\n")
        
        tmpDir = "${TMPDIR}"
        if (self.myArgs.tmpDir):
            tmpDir = os.path.abspath(self.myArgs.tmpDir)
            
        self.clusterJobFileWriter.write("df "+ tmpDir +"  -h \n")
        self.fullScratch =tmpDir+ "/scratch.bismark."+self.dName+"."+str(os.getpid())
        self.clusterJobFileWriter.write("mkdir -p "+self.fullScratch+"\n")
        self.clusterJobFileWriter.write("cd  "+self.fullScratch+"; sleep 1\n")

        captureQCRatio_command = "captureQCRatio.cc.py -b \"%s\" -c %s -r %s -e %s -o %s -k "%(os.path.abspath(self.myArgs.bamFilePattern), os.path.abspath(self.myArgs.captureTarget), self.myArgs.ratio, self.myArgs.experiment, os.path.basename(self.myArgs.outputRoot))
        self.clusterJobFileWriter.write("%s\n"%captureQCRatio_command)
        self.clusterJobFileWriter.write("sleep 1; rsync -av %s.*.xls %s/\n"%(os.path.basename(self.myArgs.outputRoot), self.mappingFolder)) 
        self.clusterJobFileWriter.write("cd  "+self.mappingFolder+"; sleep 1\n")
        
        self.clusterJobFileWriter.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitCaptureQC::generateJobCommands\n"%SimpleTime.now())
    
    def submitClusterJob(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitCaptureQC::submitClusterJob\n"%SimpleTime.now())
            
        if (self.myArgs.clusterType == "PBS" ):
            os.system("qsub "+self.clusterJobFile)
        elif(self.myArgs.clusterType == "SLURM"):
            os.system("sbatch "+self.clusterJobFile)
        else:
            sys.stderr.write("Unknown cluster type !\n")
            sys.exit(4)
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitCaptureQC::submitClusterJob\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubmitCaptureQC::work\n"%SimpleTime.now())
            
        self.generateClusterJobPreamble()
        self.generateJobCommands()
        self.submitClusterJob()  
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubmitCaptureQC::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line: %s\n"%" ".join(sys.argv))
        myArgs = SubmitCaptureQC.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = SubmitCaptureQC(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
