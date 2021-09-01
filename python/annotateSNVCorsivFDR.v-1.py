#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
import glob as glob

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleStats import SimpleStats

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class Annotate_CoRSIV_mQTL:
    DEBUG_PROGRESS                 = True
    DEBUG_LOAD_CORSIV_DEF          = True
    DEBUG_PROCESS_ONE_MQTL         = False
    DEBUG_PROCESS_ONE_MQTL         = False
    DEBUG_BETA_BY_DISTANCE         = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

run a set of corsivs against SNVs
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-m','--mQTLFile',          help='eMatrixQTL within 1mb for one CORSiVs',    required=True)
        parser.add_argument('-b','--corsivBedFile',     help='corsiv bed file',             required=True)
        parser.add_argument('-F','--FDRlimit',          help='FDR limit (default 0.05)',    default ="0.05", required=False)
        parser.add_argument('-o','--outputRoot',        help='output file root',            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    def loadCoRSIVDef(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::loadCoRSIVDef %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))

        self.corsiv_def_hash = {}
        
        corsiv_reader = open(self.myArgs.corsivBedFile, "rt")

        for line in corsiv_reader:
            
            ff = line.strip().split('\t')
            corsiv_short = ff[3]
            chrom = ff[0]
            chrom_start = int(ff[1])
            chrom_stop  = int(ff[2])
           
            corsiv_info = CStruct(corsiv_short = corsiv_short, chrom = chrom, chrom_start = chrom_start, chrom_stop = chrom_stop)
            self.corsiv_def_hash[corsiv_short]=corsiv_info
            
            if (self.DEBUG_LOAD_CORSIV_DEF):
                sys.stderr.write("Loaded corsiv definition  short: %s chrom/start/stop: %s/%s/%s \n"%(corsiv_short, self.corsiv_def_hash[corsiv_short].chrom, self.corsiv_def_hash[corsiv_short].chrom_start, self.corsiv_def_hash[corsiv_short].chrom_stop))
                
        corsiv_reader.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_mQTL::loadCoRSIVDef %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))
    
    def processOneCoRSIV(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::processOneCoRSIV %s"%(SimpleTime.now(), self.myArgs.mQTLFile))

        corsiv_mQTL_reader = open(self.myArgs.mQTLFile, "rt")
    
        annotated_corsiv = "%s.%s.all.csv"%(self.myArgs.outputRoot, os.path.basename(self.myArgs.mQTLFile))
        annotated_corsiv_writer = open(annotated_corsiv, "wt")
        
        annotated_corsiv_significant = "%s.%s.significant.csv"%(self.myArgs.outputRoot, os.path.basename(self.myArgs.mQTLFile))
        annotated_corsiv_significant_writer = open(annotated_corsiv_significant, "wt")
        
        header = ["SNV", "CoRSIV", "distance", "beta", "t-stat", "pvalue", "PerCoRSIVFDR"]
        
        annotated_corsiv_writer.write("%s\n"%"\t".join(header))
        annotated_corsiv_significant_writer.write("%s\n"%"\t".join(header))
                                
        # step 1: collect snv/corsiv/beta/pvalue/distance for SNV<=1mb
        cis_corsiv_list = []
        cis_pvalue_list = []
        
        line_idx =-1        
        for line in corsiv_mQTL_reader:
            line_idx += 1
            if (line_idx==0):
                continue
            
            ff = line.strip().split('\t')
           
            snv_info = ff[0]
            gg = snv_info.split('_')
            snv_chrom = gg[0]
            snv_chrom_pos = int(gg[1])
        
            corsiv_short = ff[1]
            corsiv_info = self.corsiv_def_hash[corsiv_short]
            
            if (snv_chrom_pos < corsiv_info.chrom_start):
                snv_corsiv_distance = corsiv_info.chrom_start - snv_chrom_pos
                pos = "U"
            elif (snv_chrom_pos > corsiv_info.chrom_stop):
                snv_corsiv_distance = snv_chrom_pos - corsiv_info.chrom_stop
                pos = "D"
            else:
                snv_corsiv_distance = 0
                pos = "O"
            
            if (snv_corsiv_distance<=1000000):
                buffer  = [ff[0], corsiv_short, snv_corsiv_distance, ff[2], ff[3], ff[4]]
                cis_corsiv_list.append(buffer)
                cis_pvalue_list.append(float(ff[4]))
                
            
        corsiv_mQTL_reader.close()

        # step 2: compute FDR
        statsHelper = SimpleStats()
        cis_qvalue_list = statsHelper.p_adjust_bh(cis_pvalue_list)
        
        q_value_limit = float(self.myArgs.FDRlimit)
        
        # step 3: write out corsiv and updated FDR
        for idx in range(len(cis_corsiv_list)):
            buffer = cis_corsiv_list[idx]
            q_value = cis_qvalue_list[idx]
            annotated_corsiv_writer.write("%s,%s\n"%(",".join([str(x) for x in buffer]),  q_value))
            if (q_value<q_value_limit):
                annotated_corsiv_significant_writer.write("%s,%s\n"%(",".join([str(x) for x in buffer]),  q_value))
                
        annotated_corsiv_significant_writer.close()        
        annotated_corsiv_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_mQTL::processOneCoRSIV %s\n"%(SimpleTime.now(), self.myArgs.mQTLFile))

    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::work\n"%SimpleTime.now())
        
        self.loadCoRSIVDef()
        self.processOneCoRSIV()
        
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
