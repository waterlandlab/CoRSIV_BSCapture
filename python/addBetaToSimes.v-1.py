#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
import glob as glob
import subprocess

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleStats import SimpleStats

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class Annotate_CoRSIV_mQTL:
    DEBUG_PROGRESS                 = True
    DEBUG_PROCESS_CORSIV           = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

run a set of corsivs against SNVs
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-m','--simesFile',         help='simes summary',    required=True)
        parser.add_argument('-b','--mQTLFolder',        help='Pearson based mQTL folder',             required=True)
        parser.add_argument('-c','--corsivIndex',       help='column w/ corsiv info (0-based)',             required=True)
        parser.add_argument('-o','--outputRoot',        help='output file root',            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::work\n"%SimpleTime.now())
        
        simesReader = open(self.myArgs.simesFile, "rt")
        simes_with_beta_file = "%s.with_beta.csv"%self.myArgs.outputRoot
        simes_with_beta_writer = open(simes_with_beta_file, "wt")
        corsiv_index = int(self.myArgs.corsivIndex)
        
        line_idx = -1
        for line in simesReader:
            line_idx +=1
            
            if (line_idx ==0):
                simes_with_beta_writer.write("%s,betaCoeffMeth\n"%line.strip())
                continue
        
            ff = line.strip().split(',')
            corsiv_short = ff[corsiv_index]
            snv_info = ff[1]
            glob_string = "%s%smQTL*%s*xls"%(self.myArgs.mQTLFolder, os.path.sep, corsiv_short)
            sys.stderr.write("Checking glob string %s\n"%glob_string)
            corsiv_mqtl_pearson_file_list = glob.glob(glob_string)
            corsiv_mqtl_pearson_file = corsiv_mqtl_pearson_file_list[0]
            
            if (self.DEBUG_PROCESS_CORSIV):
                sys.stderr.write("For corsiv_short %s found mQTL file %s\n"%(corsiv_short, corsiv_mqtl_pearson_file))

            mQTL_line = subprocess.check_output("cat %s | grep %s"%(corsiv_mqtl_pearson_file, snv_info), shell=True, text=True)
            gg = mQTL_line.strip().split('\t')
            betaCoeffMeth = gg[2]
            
            simes_with_beta_writer.write("%s,%s\n"%(line.strip(), betaCoeffMeth))
            
            if (self.DEBUG_PROCESS_CORSIV):
                sys.stderr.write("grep output %s\n"%mQTL_line.strip())
            
            
        simes_with_beta_writer.close() 
        simesReader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_mQTL::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = Annotate_CoRSIV_mQTL.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = Annotate_CoRSIV_mQTL(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
