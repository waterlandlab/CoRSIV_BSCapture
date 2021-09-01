#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "2.0"

import os, sys, argparse, re
import glob
import datetime
import math
from argparse import RawTextHelpFormatter
import gzip
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sn

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class CombineCorsivWIRITC:
    DEBUG_PROGRESS                = True
    DEBUG_COMBINE                 = True
    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.
            
combine CORSiV maxWIR with ITC and IIR
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-i','--itc',           help='corsiv ITC and IIR',                  required=True)
        parser.add_argument('-m','--maxWIR',        help='max within individual range file',    required=True)
        parser.add_argument('-o','--outputRoot',    help='output file root',                    required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args


    def loadITC(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineCorsivWIRITC::loadITC\n"%SimpleTime.now())
            
        self.data_itc = pd.read_csv(self.myArgs.itc, sep=",", index_col=0)
        
        if (self.DEBUG_COMBINE):
            sys.stderr.write("CORSiV ITC %s x %s\n"%(self.data_itc.shape[0], self.data_itc.shape[1]))
        
        self.corsiv_itc_index_hash = {}
        self.corsiv_itc_index = list(self.data_itc.index)
        
        for corsiv_itc_idx in range(len(self.corsiv_itc_index)):
            self.corsiv_itc_index_hash[self.corsiv_itc_index[corsiv_itc_idx]]=corsiv_itc_idx
            if (self.DEBUG_COMBINE):
                sys.stderr.write("Short corsiv %s\n "%self.corsiv_itc_index[corsiv_itc_idx])
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineCorsivWIRITC::loadITC\n"%SimpleTime.now())
            
    def loadWIR(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineCorsivWIRITC::loadWIR\n"%SimpleTime.now())
            
        self.data_wir = pd.read_csv(self.myArgs.maxWIR, sep="\t", index_col = 0)
        
        if (self.DEBUG_COMBINE):
            sys.stderr.write("CORSiV WIR %s x %s\n"%(self.data_wir.shape[0], self.data_wir.shape[1]))
                    
        # redo the index
        long_name_index = list(self.data_wir.index)
        short_id_list = []
        
        keep_corsiv_list = []
        
        short_re = re.compile('idx.((\d+)_(\d+))_')
        for idx in range(len(long_name_index)):
            corsiv_long_name = long_name_index[idx]
            short_match = re.search(short_re, corsiv_long_name)
            if (short_match):
                corsiv_short_name = short_match.group(1)
                short_id_list.append(corsiv_short_name)
                keep_corsiv_list.append(idx)
                if (self.DEBUG_COMBINE):
                    sys.stderr.write("long corsiv name %s short corsiv name %s\n"%(corsiv_long_name, corsiv_short_name))
                    
            else:
                short_id_list.append(corsiv_long_name)
                if (self.DEBUG_COMBINE):
                    sys.stderr.write("long corsiv name %s is control\n"%corsiv_long_name)
                    
                
        self.data_wir.index = short_id_list
        self.data_wir = self.data_wir.iloc[keep_corsiv_list,:]
            
        # merge wir dataframe with itc data frame
        merged_df = pd.merge(self.data_wir, self.data_itc.drop_duplicates(), left_index = True, right_index=True, how='left')
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("merge wir %s x %s itc %s x %s --> %s x %s\n"%(self.data_wir.shape[0], self.data_wir.shape[1], self.data_itc.shape[0], self.data_itc.shape[1], merged_df.shape[0], merged_df.shape[1]))

        merged_data_csv = "%s.merged.xls"%self.myArgs.outputRoot
        
        merged_df.to_csv(merged_data_csv, sep="\t", index_label="Short_CORSiV")
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineCorsivWIRITC::loadWIR\n"%SimpleTime.now())
            
    
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineCorsivWIRITC::work\n"%SimpleTime.now())
            
        self.loadITC()
        self.loadWIR()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineCorsivWIRITC::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = CombineCorsivWIRITC.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = CombineCorsivWIRITC(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
