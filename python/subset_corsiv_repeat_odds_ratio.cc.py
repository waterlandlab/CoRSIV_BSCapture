#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
import glob as glob
import gzip
import subprocess
import math
import scipy.stats

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleStats import SimpleStats
from CCLabUtils.simpleData import SimpleData

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class SubsetCorsivRepeatsOddsRatios:
    DEBUG_PROGRESS              = True
    DEBUG_LOAD_ODDS_RATIOS      = True
    DEBUG_SUBSET_OR      = True
    DEBUG_ADJUST_P      = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Analyze combined repeats odds-ratio & p-values by applying an FDR then subsetting only significant odds-ratios
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-r','--repeatOddsRatioFile',   help='repeats odds-ratio & p-values file',  required=True)
        parser.add_argument('-f','--fdrLimit',              help='fdr limit (default 0.05)',      required=False, type=float, default=0.05)
        parser.add_argument('-o','--outputRoot',            help='output file root',  required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
    
   
    def loadOddsRatiosAndPvalues(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubsetCorsivRepeatsOddsRatios::loadOddsRatiosAndPvalues\n"%SimpleTime.now())

        self.or_pv_df = pd.read_csv(self.myArgs.repeatOddsRatioFile, sep="\t", index_col=0, header=0)
       
        self.columns = list(self.or_pv_df.columns)
        self.number_of_radii = int(len(self.columns)/2)
        if (self.DEBUG_LOAD_ODDS_RATIOS):
            sys.stderr.write("Data frame from %s has %s columns and %s radii\n"%(self.myArgs.repeatOddsRatioFile, len(self.columns), self.number_of_radii))
        
        self.radius_stats_hash = {}
        for idx in range(self.number_of_radii):
            if (self.DEBUG_LOAD_ODDS_RATIOS):
                sys.stderr.write("\nProcessing p-values and q-values for column %s %s\n"%(idx, self.columns[idx]))
            p_values = self.or_pv_df.iloc[:,idx+self.number_of_radii]
            statsHelper = SimpleStats()
            q_values = statsHelper.p_adjust_bh(p_values)
            radius_info = CStruct(idx = idx, p_values = p_values, q_values = q_values)
            if (self.DEBUG_ADJUST_P):
                sys.stderr.write("pvalues %s\n"%"\t".join([str(x) for x in list(p_values)]))
                sys.stderr.write("qvalues %s\n"%"\t".join([str(x) for x in list(q_values)]))
                
            self.radius_stats_hash[idx]=radius_info
            
        if (self.DEBUG_LOAD_ODDS_RATIOS):
                sys.stderr.write("Done Processing p-values and q-values\n")
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubsetCorsivRepeatsOddsRatios::loadOddsRatiosAndPvalues\n"%SimpleTime.now())
            
    def subsetCoRSIVRepeats(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubsetCorsivRepeatsOddsRatios::subsetCoRSIVRepeats\n"%SimpleTime.now())
            
        subset_columns = self.columns[:self.number_of_radii]
        self.subset_or = pd.DataFrame(0, index=self.or_pv_df.index,
                                      columns = subset_columns)
        index  = list(self.subset_or.index)
        
        for col_idx in range(self.number_of_radii):
            if (self.DEBUG_SUBSET_OR):
                sys.stderr.write("column %s --> %s\n"%(col_idx, self.columns[col_idx]))
                
            radius_info = self.radius_stats_hash[col_idx]
            
            for row_idx in range(len(index)):
                q_value =radius_info.q_values[row_idx] 
                if (q_value<self.fdr_limit) :
                    current_or = self.or_pv_df.iat[row_idx, col_idx]
                    if (current_or<2**-20):
                        current_or = 2**-20
                    if (current_or>2**20):
                        current_or = 2**20
                    log2_current_or = math.log(current_or)/math.log(2)
                    self.subset_or.iloc[row_idx, col_idx] = log2_current_or
                    if (self.DEBUG_SUBSET_OR):
                        sys.stderr.write("[%s : %s ,%s : %s] q=%s < %s PASS linear %s log2 %s %s\n"%(row_idx, index[row_idx], col_idx, subset_columns[col_idx], q_value, self.fdr_limit, current_or, log2_current_or, self.subset_or.iat[row_idx, col_idx]))
                else:
                    if (self.DEBUG_SUBSET_OR):
                        sys.stderr.write("[%s : %s ,%s : %s] q=%s < %s FAIL\n"%(row_idx, index[row_idx], col_idx, subset_columns[col_idx], q_value, self.fdr_limit))

        output_file = "%s.fdr_limit.%s.xls"%(self.myArgs.outputRoot, self.myArgs.fdrLimit)
        self.subset_or.to_csv(output_file, sep="\t", index_label="Repeats")
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubsetCorsivRepeatsOddsRatios::subsetCoRSIVRepeats\n"%SimpleTime.now())
                      
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubsetCorsivRepeatsOddsRatios::work\n"%SimpleTime.now())
        
        self.fdr_limit = float(self.myArgs.fdrLimit)
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("fdrLimit %s\n"%self.fdr_limit)
        
        self.loadOddsRatiosAndPvalues()
        self.subsetCoRSIVRepeats()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubsetCorsivRepeatsOddsRatios::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = SubsetCorsivRepeatsOddsRatios.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = SubsetCorsivRepeatsOddsRatios(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
