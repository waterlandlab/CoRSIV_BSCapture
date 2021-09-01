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
        
class Select_SNP_CoRSIV_by_FDR:
    DEBUG_PROGRESS                 = True
    DEBUG_LOAD_CORSIV_DEF          = True
    DEBUG_LOAD_ALL_CORSIVs         = True
    DEBUG_PROCESS_ONE_MQTL         = True
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
         
        parser.add_argument('-m','--simesOutputFile',   help='eMatrixQTL within 1mb for one CORSiVs',    required=True)
        parser.add_argument('-b','--corsivBedFile',     help='corsiv bed file',             required=True)
        parser.add_argument('-F','--FDRlimit',          help='FDR limit (default 0.05)',    default ="0.05", required=False)
        parser.add_argument('-a','--allCoRSIVFile',     help='corsiv/snp/simes files, with FDR field added ',  required=True)
        parser.add_argument('-o','--outputRoot',        help='output file root',            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    def loadCoRSIVDef(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_SNP_CoRSIV_by_FDR::loadCoRSIVDef %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))

        self.corsiv_def_hash = {}
        
        corsiv_reader = open(self.myArgs.corsivBedFile, "rt")

        for line in corsiv_reader:
            
            ff = line.strip().split('\t')
            corsiv_short = ff[3].replace('"','')
            chrom = ff[0]
            chrom_start = int(ff[1])
            chrom_stop  = int(ff[2])
           
            corsiv_info = CStruct(corsiv_short = corsiv_short, chrom = chrom, chrom_start = chrom_start, chrom_stop = chrom_stop)
            self.corsiv_def_hash[corsiv_short]=corsiv_info
            
            if (self.DEBUG_LOAD_CORSIV_DEF):
                sys.stderr.write("Loaded corsiv definition  short: %s chrom/start/stop: %s/%s/%s \n"%(corsiv_short, self.corsiv_def_hash[corsiv_short].chrom, self.corsiv_def_hash[corsiv_short].chrom_start, self.corsiv_def_hash[corsiv_short].chrom_stop))
                
        corsiv_reader.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_SNP_CoRSIV_by_FDR::loadCoRSIVDef %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))
    
    def processSimesCoRSIV(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_SNP_CoRSIV_by_FDR::processOneCoRSIV %s"%(SimpleTime.now(), self.myArgs.simesOutputFile))

    
        significant_corsiv_snp_file = "%s.%s.all.csv"%(self.myArgs.outputRoot, os.path.basename(self.myArgs.simesOutputFile))
        significant_corsiv_snp_writer = open(significant_corsiv_snp_file, "wt")
        
        
        header = ["SNV", "CoRSIV", "distance", "pvalue", "simes_qvalue"]
        
        significant_corsiv_snp_writer.write("%s\n"%",".join(header))
                                
        # step 1: collect snv/corsiv/beta/pvalue/distance for SNV<=1mb
        line_idx =-1
        
        simes_output_reader = open(self.myArgs.simesOutputFile, "rt")
        for line in simes_output_reader:
            line_idx += 1
            if (line_idx==0):
                continue
            
            ff = line.strip().split(',')
           
            snv_info = ff[2]
            gg = snv_info.split('_')
            snv_chrom = gg[0]
            snv_chrom_pos = int(gg[1])
        
            corsiv_short = ff[3].replace('"','')
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
                current_pvalue = float(ff[4])
                current_simes = float(ff[5])
                
                if (current_simes<=self.max_simes_under_limit):
                    if (self.DEBUG_PROCESS_ONE_MQTL):
                        sys.stderr.write("Passing corsiv/snp: %s simes %s vs max_simes %s\n"%(line.strip(), current_simes, self.max_simes_under_limit))
                    buffer = [snv_info, corsiv_short, snv_corsiv_distance, current_pvalue, current_simes]
                    significant_corsiv_snp_writer.write("%s\n"%",".join([str(x) for x in buffer]))
                    
        simes_output_reader.close()
        significant_corsiv_snp_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_SNP_CoRSIV_by_FDR::processOneCoRSIV %s\n"%(SimpleTime.now(), self.myArgs.simesOutputFile))

    def loadAllCoRSIVs(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_SNP_CoRSIV_by_FDR::loadAllCoRSIVs\n"%SimpleTime.now())
    
        all_corsiv_reader = open(self.myArgs.allCoRSIVFile, "rt")
        
        self.max_simes_under_limit = 0
        self.max_fdr_value_under_limit = 0
        
        self.fdr_limit = float(self.myArgs.FDRlimit)
        
        line_idx = -1
        for line in all_corsiv_reader:
            line_idx  += 1
            if (line_idx==0):
                continue
            ff = line.strip().split(',')
            line_fdr = float(ff[len(ff)-1])
            line_simes = float(ff[len(ff)-3])
            
            if ((line_fdr<self.fdr_limit) and (self.max_fdr_value_under_limit<line_fdr)):
                self.max_fdr_value_under_limit = line_fdr
                self.max_simes_under_limit = line_simes
                if (self.DEBUG_LOAD_ALL_CORSIVs):
                    sys.stderr.write("LoadAllCoRSIV: line %s line_fdr %s bumped max_fdr_under limit to %s simes to %s\n"%(line.strip(), line_fdr, self.max_fdr_value_under_limit, self.max_simes_under_limit))
            
        all_corsiv_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("LoadAllCoRSIV: max fdr under limit %s max simes under limit %s\n"%(self.max_fdr_value_under_limit, self.max_simes_under_limit))
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_SNP_CoRSIV_by_FDR::work\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_SNP_CoRSIV_by_FDR::work\n"%SimpleTime.now())
        
        # load all corsivs def
        self.loadCoRSIVDef()
        
        # load all corsivs / SNP / Simes file; determine limit of Simes for significance
        self.loadAllCoRSIVs()
        
        # find significant SNPs; add distance to the significant SNPs in the process
        self.processSimesCoRSIV()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_SNP_CoRSIV_by_FDR::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = Select_SNP_CoRSIV_by_FDR.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = Select_SNP_CoRSIV_by_FDR(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
