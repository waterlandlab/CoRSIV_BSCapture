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
        
class AssessVersion2CaptureCoRSIV:
    DEBUG_PROGRESS                  = True
    DEBUG_LOAD_BETA_VALUES          = True
    DEBUG_VIOLIN_PLOTS              = True
    DEBUG_SETUP_DELTA_DF            = True

    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Assess differences between methylation values in v1 and v2 of the CoRSIV capture array and plot violin plots
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-x','--corsivMethylationMatrix',   help='data matrix with beta values for CoRSiVs',    required=True)
        parser.add_argument('-o','--outputRoot',                help='output file root',                            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

 
    def boxPlotSeries(self, series_hash, output_file_root, boxplot_title):
        pdf_file = "%s.pdf"%output_file_root
        jpg_file = "%s.jpg"%output_file_root

        fig, ax = plt.subplots(1)
        box_data = []
        labels = []
        for series_key in series_hash:
            box_data.append(series_hash[series_key])
            labels.append(series_key)

        box = plt.boxplot(box_data,patch_artist=True,whis="range")
        ax.set_ylim(0, 1)
        plt.title(boxplot_title)
        plt.xticks(range(1,len(box_data)+1),labels,rotation=90)

        plt.savefig(pdf_file, bbox_inches='tight',dpi=300)
        plt.savefig(jpg_file, bbox_inches='tight',dpi=300)

    def setupTissueDeltaDataFrame(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CaptureCoRSIV::setupTissueDeltaDataFrame\n"%SimpleTime.now())
            
        self.tissue_delta_hash = {}
        corsiv_index = list(self.df_beta_values.index)
        
        number_of_valid_deltas = 0
        data_frame_list_of_lists = []
        for tissue in self.tissue_version_hash:
            tissue_info = self.tissue_version_hash[tissue]
            if (self.DEBUG_SETUP_DELTA_DF):
                sys.stderr.write("Setting up Delta DF for tissue %s %s : v1=%s v2=%s\n"%(tissue, tissue_info.tissue, tissue_info.v1, tissue_info.v2))
            if (tissue_info.v1>=0 and tissue_info.v2>=0):
                if (self.DEBUG_SETUP_DELTA_DF):
                    sys.stderr.write("v1=%s v2=%s both present\n"%(tissue_info.v1, tissue_info.v2))
                    
                self.tissue_delta_hash[tissue]=[]
                for row_idx in range(len(corsiv_index)):
                    beta_value_v1 = self.df_beta_values.iat[row_idx, tissue_info.v1]
                    beta_value_v2 = self.df_beta_values.iat[row_idx, tissue_info.v2]
                    if (np.isnan(beta_value_v1) or np.isnan(beta_value_v2)):
                        if (self.DEBUG_SETUP_DELTA_DF):
                            sys.stderr.write("Eek v1=%s v2=%s \n"%( beta_value_v1, beta_value_v2))
                        continue
                    else:
                        delta_beta = beta_value_v2 - beta_value_v1
                        if (self.DEBUG_SETUP_DELTA_DF):
                            sys.stderr.write("Yay v1=%s v2=%s deltaBeta %s \n"%( beta_value_v1, beta_value_v2, delta_beta))
                        self.tissue_delta_hash[tissue].append(delta_beta)
                        data_frame_list_of_lists.append([tissue, delta_beta])
                
                number_of_valid_deltas += len(self.tissue_delta_hash[tissue])
                if (self.DEBUG_SETUP_DELTA_DF):
                    sys.stderr.write("Tissue %s contributes %s delta betas\n"%(tissue, len(self.tissue_delta_hash[tissue])))
                    
        # setup and load delta beta data frame
        columns_delta_beta = ["Tissue", "DeltaBeta"]
        self.delta_beta_data_frame = pd.DataFrame( data_frame_list_of_lists, columns = columns_delta_beta)
        
        if (self.DEBUG_SETUP_DELTA_DF):
            sys.stderr.write("DeltaBeta matrix shape %s columns %s\n"%(str(self.delta_beta_data_frame.shape), ";".join(list(self.delta_beta_data_frame.columns))))
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessVersion2CaptureCoRSIV::setupTissueDeltaDataFrame\n"%SimpleTime.now())
    
    def loadBetaValues(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CaptureCoRSIV::loadBetaValues\n"%SimpleTime.now())
    
        self.df_beta_values = pd.read_csv(self.myArgs.corsivMethylationMatrix, sep="\t", index_col=0, header=0)
        if (self.DEBUG_LOAD_BETA_VALUES):
            sys.stderr.write("Loaded beta values for 2 rounds: dim %s \n"%str(self.df_beta_values.shape))
        
        # setup v1/v2 index
        self.tissue_hash = {}
        column_list = list(self.df_beta_values.columns)
        self.tissue_version_hash ={}
        
        for col_idx in range(len(column_list)):
            column_name = column_list[col_idx]
            bits = column_name.split('_')
            tissue = bits[0]
            version = column_name[-2:]
            
            if (self.DEBUG_LOAD_BETA_VALUES):
                sys.stderr.write("Parsing column  [%s] %s --> tissue %s version %s \n"%(col_idx, column_name, tissue, version))
                
            if not (tissue in self.tissue_version_hash):
                tissue_info = CStruct(tissue = tissue, v1=-1, v2=-1)
                self.tissue_version_hash[tissue]=tissue_info
            
            tissue_info = self.tissue_version_hash[tissue]
            if (version =="v1"):
                tissue_info.v1 = col_idx
            else:
                tissue_info.v2 = col_idx
        
        self.setupTissueDeltaDataFrame()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessVersion2CaptureCoRSIV::loadBetaValues\n"%SimpleTime.now())
    
    def outputViolinPlots(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CaptureCoRSIV::outputViolinPlots\n"%SimpleTime.now())
        
        violin_plot_pdf = "%s.violin_plot.pdf"%self.myArgs.outputRoot
        
        sn.set(style="whitegrid")
        ax = sn.violinplot(x="Tissue", y="DeltaBeta", data=self.delta_beta_data_frame)
        ax.set_ylim(-1,1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessVersion2CaptureCoRSIV::outputViolinPlots\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessVersion2CaptureCoRSIV::work\n"%SimpleTime.now())
    
            
        self.loadBetaValues()
        self.outputViolinPlots()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessVersion2CaptureCoRSIV::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = AssessVersion2CaptureCoRSIV.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = AssessVersion2CaptureCoRSIV(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
