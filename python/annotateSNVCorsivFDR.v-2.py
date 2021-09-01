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
        
        
        header = ["SNV", "CoRSIV", "distance", "beta", "t-stat", "pvalue", "PerCoRSIVFDR"]
        
        annotated_corsiv_writer.write("%s\n"%"\t".join(header))
                                
        # step 1: collect snv/corsiv/beta/pvalue/distance for SNV<=1mb
        cis_corsiv_list = []
        cis_pvalue_list = []
        
        line_idx =-1
        min_pvalue_beta = False
        min_pvalue = 1
        
        tmp_simes_input = "%s.tmp_simes_input.csv"%self.myArgs.outputRoot
        tmp_simes_output = "%s.tmp_simes_output.csv"%self.myArgs.outputRoot
        r_simes_log = "%s.log_simes.txt"%self.myArgs.outputRoot
        
        tmp_simes_header = ["Link_index", "snv", "corsiv_short", "pValue"]
        tmp_simes_input_writer = open(tmp_simes_input, "wt")
        tmp_simes_input_writer.write("%s\n"%",".join(tmp_simes_header))
        
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
                current_beta = float(ff[2])
                current_pvalue = float(ff[4])
                
                simes_input_buffer = ["link_%s"%line_idx, snv_info, corsiv_short, current_pvalue]
                
                tmp_simes_input_writer.write("%s\n"%",".join([str(x) for x in simes_input_buffer]))
                
                if (min_pvalue>current_pvalue):
                    min_pvalue = current_pvalue
                    min_pvalue_beta = current_beta
        
        tmp_simes_input_writer.close()    
        corsiv_mQTL_reader.close()

        # step 2: compute Simes
        r_cmd_simes = "simes_correction.cc.R %s %s &> %s"%(tmp_simes_input, tmp_simes_output, r_simes_log)
        os.system(r_cmd_simes)
        
        simes_output_reader = open(tmp_simes_output, "rt")
        
        min_simes_qvalue = 1
        min_simes_snp = ""
        line_idx = -1
        for line in simes_output_reader:
            line_idx += 1
            if (line_idx ==0):
                continue
            ff = line.strip().split(',')
            q_value =float(ff[5])
            if (q_value<min_simes_qvalue):
                min_simes_snp = ff[2]
                min_simes_qvalue = q_value
                min_simes_pvalue = ff[4]
        
        simes_output_reader.close()
        
        # step 3: write out corsiv and updated min Simes p-value
                
        annotated_corsiv_significant = "%s.%s.simes.csv"%(self.myArgs.outputRoot, os.path.basename(self.myArgs.mQTLFile))
        annotated_corsiv_significant_writer = open(annotated_corsiv_significant, "wt")
        annotated_corsiv_significant_writer.write("CoRSIV_short,snp,min_pValue,simes_qvalue\n")
        significant_simes_data = [corsiv_short, min_simes_snp, min_simes_pvalue, min_simes_qvalue]
        annotated_corsiv_significant_writer.write("%s\n"%",".join([str(x) for x in significant_simes_data]))
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
