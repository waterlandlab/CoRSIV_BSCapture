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
import math

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleStats import SimpleStats

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class DistillCoRSIVSNPMap:
    DEBUG_PROGRESS                      = True
    DEBUG_SETUP                         = True
    DEBUG_SUMMARIZE_CORSIV_SNP_MAP      = False

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

summarize corsiv snp associations in a distance heatmap
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-m','--twoStepmQTLFile',   help='final two-step mQTL file',    required=True)
        parser.add_argument('-w','--windowSize',        help='window size',                 required=True)
        parser.add_argument('-d','--distanceCap',       help='distance cap',                required=True)
        parser.add_argument('-o','--outputRoot',        help='output file root',            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    def setupAnalysis(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START DistillCoRSIVSNPMap::setupAnalysis\n"%SimpleTime.now())
    
        self.distanceCap = abs(int(self.myArgs.distanceCap))
        self.windowSize = abs(int(self.myArgs.windowSize))
        self.numberOfWindows = 2*(self.distanceCap/self.windowSize)
        if (self.DEBUG_SETUP):
            sys.stderr.write("distance cap %s window size %s number of windows %s\n"%(self.distanceCap, self.windowSize, self.numberOfWindows))
        
        self.corsiv_bin_struct = {}
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START DistillCoRSIVSNPMap::setupAnalysis\n"%SimpleTime.now())


    
    def summarizeCorsivSNPMap(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START DistillCoRSIVSNPMap::summarizeCorsivSNPMap\n"%SimpleTime.now())
    
        corsiv_reader = open(self.myArgs.twoStepmQTLFile, "rt")
        
        line_idx = -1
        for line in corsiv_reader:
            line_idx += 1
            if (line_idx==0):
                continue
            
            ff = line.strip().split(',')
            corsiv = ff[1]
            snp_corsiv_distance = int(ff[2])
            normalized_distance = snp_corsiv_distance + self.distanceCap
            pvalue = float(ff[3])
            if (pvalue<10**-200):
                pvalue = 10**-200
            minus_log10_pvalue = -math.log(pvalue)/math.log(10)
            beta_coeff = float(ff[5])
            
            if (normalized_distance<0):
                normalized_distance = 0
            if (normalized_distance>=2*self.distanceCap):
                normalized_distance = 2*self.distanceCap-1
            bin_index = normalized_distance/self.windowSize
            if (self.DEBUG_SUMMARIZE_CORSIV_SNP_MAP):
                sys.stderr.write("Line %s: corsiv %s distance %s bin_index %s -log10(pvalue) %s beta_coeff %s\n"%(line.strip(), corsiv, snp_corsiv_distance, bin_index, minus_log10_pvalue, beta_coeff))
            
            bin_address = "bin_%s"%bin_index
            if not (corsiv in self.corsiv_bin_struct):
                self.corsiv_bin_struct[corsiv]={}
            if not (bin_address in self.corsiv_bin_struct[corsiv]):
                bin_info = CStruct(corsiv = corsiv, bin_address = bin_address, max_abs_beta=0, max_minus_log10_pvalue=0, snp_count = 0, bin_index = bin_index)
                self.corsiv_bin_struct[corsiv][bin_address]=bin_info
            bin_info = self.corsiv_bin_struct[corsiv][bin_address]
            bin_info.snp_count += 1
            if (abs(beta_coeff)>abs(bin_info.max_abs_beta)):
                if (self.DEBUG_SUMMARIZE_CORSIV_SNP_MAP):
                    sys.stderr.write("bump up max_abs_beta from %s to %s\n"%(bin_info.max_abs_beta, beta_coeff))
                bin_info.max_abs_beta = beta_coeff
            if (minus_log10_pvalue>bin_info.max_minus_log10_pvalue):
                if (self.DEBUG_SUMMARIZE_CORSIV_SNP_MAP):
                    sys.stderr.write("bump up minus_log10_pvalue from %s to %s\n"%(bin_info.max_minus_log10_pvalue, minus_log10_pvalue))
                bin_info.max_minus_log10_pvalue = minus_log10_pvalue
        
        corsiv_reader.close()
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP DistillCoRSIVSNPMap::summarizeCorsivSNPMap\n"%SimpleTime.now())



    def outputCorsivSNPMap(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START DistillCoRSIVSNPMap::outputCorsivSNPMap\n"%SimpleTime.now())
    
        
        output_report_counts    = "%s.counts.xls"%self.myArgs.outputRoot    
        output_report_beta      = "%s.beta.xls"%self.myArgs.outputRoot    
        output_report_pvalue    = "%s.minusLog10pvalue.xls"%self.myArgs.outputRoot
        
        output_report_counts_writer = open(output_report_counts, "wt")
        
        output_report_beta_writer = open(output_report_beta, "wt")
        
        output_report_pvalue_writer = open(output_report_pvalue, "wt")
        
        header =["corsiv"]
        for bin_index in range(self.numberOfWindows):
            header.append("bin_%s"%bin_index)
        
        output_report_counts_writer.write("%s\n"%",".join(header))
        output_report_beta_writer.write("%s\n"%",".join(header))
        output_report_pvalue_writer.write("%s\n"%",".join(header))
        
        for corsiv_name in self.corsiv_bin_struct:
            corsiv_hash = self.corsiv_bin_struct[corsiv_name]
            
            count_vector = np.zeros(self.numberOfWindows)
            beta_vector = np.zeros(self.numberOfWindows)
            pvalue_vector = np.zeros(self.numberOfWindows)
            
            for bin_address in corsiv_hash:
                bin_info = corsiv_hash[bin_address]
                bin_index = bin_info.bin_index
                count_vector[bin_index]=bin_info.snp_count
                beta_vector[bin_index]=bin_info.max_abs_beta
                pvalue_vector[bin_index]=bin_info.max_minus_log10_pvalue
        
            output_report_counts_writer.write("%s,%s\n"%(corsiv_name, ",".join([str(x) for x in count_vector])))
            output_report_beta_writer.write("%s,%s\n"%(corsiv_name, ",".join([str(x) for x in beta_vector])))
            output_report_pvalue_writer.write("%s,%s\n"%(corsiv_name, ",".join([str(x) for x in pvalue_vector])))
        
        output_report_pvalue_writer.close()
        output_report_beta_writer.close()
        output_report_counts_writer.close()
    
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP DistillCoRSIVSNPMap::outputCorsivSNPMap\n"%SimpleTime.now())

    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START DistillCoRSIVSNPMap::work\n"%SimpleTime.now())
        
        
        self.setupAnalysis() 
        self.summarizeCorsivSNPMap()
        self.outputCorsivSNPMap()

        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP DistillCoRSIVSNPMap::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = DistillCoRSIVSNPMap.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = DistillCoRSIVSNPMap(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
