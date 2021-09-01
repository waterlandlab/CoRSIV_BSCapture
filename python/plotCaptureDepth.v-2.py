#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "2.0"

import os, sys, argparse, re
import glob
import datetime
import math
import scipy.stats
from argparse import RawTextHelpFormatter
import random
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
        
class AssessVersion2CoverageCoRSIV:
    DEBUG_PROGRESS                  = True
    DEBUG_LOAD_COVERAGE_COUNTS          = True
    DEBUG_VIOLIN_PLOTS              = True
    DEBUG_SETUP_DELTA_DF            = True

    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Assess basepair counts for each sample and each corsiv
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-c','--corsivCoverageCounts',      help='data matrix with read counts for CoRSiVs',    required=True)
        parser.add_argument('-o','--outputRoot',                help='output file root',                            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def loadCoverageValues(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CoverageCoRSIV::loadCoverageValues\n"%SimpleTime.now())
    
        self.df_coverage = pd.read_csv(self.myArgs.corsivCoverageCounts, sep="\t", index_col=0, header=0)
        if (self.DEBUG_LOAD_COVERAGE_COUNTS):
            sys.stderr.write("Loaded coverage counts values dim %s \n"%str(self.df_coverage.shape))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CoverageCoRSIV::setupTissueDeltaDataFrame\n"%SimpleTime.now())
            
        corsiv_index = list(self.df_coverage.index)
        data_frame_list_of_lists = []
        corsiv_columns = list(self.df_coverage.columns)
        
        for row_idx in range(len(corsiv_index)):
            corsiv_name = corsiv_index[row_idx]
            bits = corsiv_name.split('_')
            bp_length = int(bits[2])-int(bits[1])+1
            
            if (self.DEBUG_LOAD_COVERAGE_COUNTS):
                sys.stderr.write("Setting up BP coverage for corsiv %s bounds %s %s length %s\n"%(corsiv_name, bits[1], bits[2], bp_length))
            
            for col_idx in range(2, len(corsiv_columns)):
                tissue = corsiv_columns[col_idx]
                tissue_bits = tissue.split('_')
                counts = self.df_coverage.iat[row_idx, col_idx]
                bp_coverage = float(counts)*150/float(bp_length)
                if (self.DEBUG_LOAD_COVERAGE_COUNTS):
                    sys.stderr.write("Col [%s] %s reads %s bp_coverage  %s\n"%(col_idx, tissue, counts, bp_coverage))
                
                data_frame_list_of_lists.append([tissue_bits[0], bp_coverage])
                    
        # setup and load delta beta data frame
        columns_bp_coverage = ["Tissue", "BPCoverage"]
        self.bp_coverage_data_frame = pd.DataFrame( data_frame_list_of_lists, columns = columns_bp_coverage)
        
        if (self.DEBUG_LOAD_COVERAGE_COUNTS):
            sys.stderr.write("BP coverage matrix shape %s columns %s\n"%(str(self.bp_coverage_data_frame.shape), ";".join(list(self.bp_coverage_data_frame.columns))))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessVersion2CoverageCoRSIV::loadCoverageValues\n"%SimpleTime.now())
    
    def outputViolinPlots(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CoverageCoRSIV::outputViolinPlots\n"%SimpleTime.now())
        
        violin_plot_pdf = "%s.violin_plot.pdf"%self.myArgs.outputRoot
        
        sn.set(style="whitegrid")
        ax = sn.violinplot(x="Tissue", y="BPCoverage", data=self.bp_coverage_data_frame)
        #ax.set_ylim(-1,1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)
    
        violin_plot_pdf = "%s.violin_plot_0_100x.pdf"%self.myArgs.outputRoot
        ax = sn.violinplot(x="Tissue", y="BPCoverage", data=self.bp_coverage_data_frame)
        ax.set_ylim(0, 100)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)

        violin_plot_pdf = "%s.violin_plot_0_200x.pdf"%self.myArgs.outputRoot
        ax = sn.violinplot(x="Tissue", y="BPCoverage", data=self.bp_coverage_data_frame)
        ax.set_ylim(0, 200)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)
                    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessVersion2CoverageCoRSIV::outputViolinPlots\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CoverageCoRSIV::work\n"%SimpleTime.now())
    
            
        self.loadCoverageValues()
        self.outputViolinPlots()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessVersion2CoverageCoRSIV::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = AssessVersion2CoverageCoRSIV.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = AssessVersion2CoverageCoRSIV(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
