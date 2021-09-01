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
        
class CoRSIVCoverageReport:
    DEBUG_PROGRESS                  = True
    DEBUG_LOAD_COVERAGE_COUNTS      = False
    DEBUG_VIOLIN_PLOTS              = True
    DEBUG_LOAD_CORSIV               = False
    DEBUG_BINARIZE_COVERAGE         = False

    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Assess basepair counts for each sample and each corsiv
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-c','--corsivCoverageCounts',   help='data matrix with read counts for CoRSiVs',   required=True)
        parser.add_argument('-r','--readLength',             help='read length',                                required=True)
        parser.add_argument('-d','--corsivDefinition',       help='corsiv definition',                          required=True)
        parser.add_argument('-o','--outputRoot',             help='output file root',                           required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def loadCoverageValues(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CoRSIVCoverageReport::loadCoverageValues\n"%SimpleTime.now())
    
        self.df_read_counts = pd.read_csv(self.myArgs.corsivCoverageCounts, sep="\t", index_col=0, header=0)
        if (self.DEBUG_LOAD_COVERAGE_COUNTS):
            sys.stderr.write("Loaded coverage counts values dim %s \n"%str(self.df_read_counts.shape))
        
        self.corsiv_index = list(self.df_read_counts.index)
        self.corsiv_info_index = []
        for row_idx in range(len(self.corsiv_index)):
            corsiv_coordinates = self.corsiv_index[row_idx]
            self.corsiv_info_index.append(self.corsiv_genomic_coord_to_name_hash[corsiv_coordinates])
            
        self.corsiv_columns_plus_2 = list(self.df_read_counts.columns)
        self.corsiv_columns_only = self.corsiv_columns_plus_2[2:]
        self.read_length = float(self.myArgs.readLength)
        
        self.df_bp_coverage = pd.DataFrame(0.0, index=self.corsiv_info_index, columns=self.corsiv_columns_only)
        
        # compute bp coverage
        for row_idx in range(len(self.corsiv_index)):
            corsiv_coordinates = self.corsiv_index[row_idx]
            corsiv_info = self.corsiv_genomic_coord_to_name_hash[corsiv_coordinates]
            bits = corsiv_coordinates.split('_')
            bp_length = int(bits[2])-int(bits[1])+1
            
            if (self.DEBUG_LOAD_COVERAGE_COUNTS):
                sys.stderr.write("Setting up BP coverage for corsiv %s info %s === %s bounds %s %s length %s\n"%(corsiv_coordinates, corsiv_info, self.corsiv_info_index[row_idx],  bits[1], bits[2], bp_length))
            
            for col_idx in range(2, len(self.corsiv_columns_plus_2)):
                tissue = self.corsiv_columns_plus_2[col_idx]
                counts = self.df_read_counts.iat[row_idx, col_idx]
                bp_coverage = float(counts)*self.read_length/float(bp_length)
                
                if (self.DEBUG_LOAD_COVERAGE_COUNTS):
                    sys.stderr.write("Col [%s] %s reads %s bp_coverage  %s\n"%(col_idx, tissue, counts, bp_coverage))
                
                self.df_bp_coverage.iat[row_idx, col_idx-2]=bp_coverage
                    
        if (self.DEBUG_LOAD_COVERAGE_COUNTS):
            sys.stderr.write("BP coverage matrix shape %s\n"%str(self.df_bp_coverage.shape))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CoRSIVCoverageReport::loadCoverageValues\n"%SimpleTime.now())

        
    def outputBinaryHeatmaps(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CoRSIVCoverageReport::outputBinaryHeatmaps\n"%SimpleTime.now())
        
        full_bp_coverage_matrix_file = "%s.full_bp_coverage.xls"%self.myArgs.outputRoot
        self.df_bp_coverage.to_csv(full_bp_coverage_matrix_file, sep='\t', index_label="CoRSIV")

        self.limits = [30, 50, 100, 150, 200]
        
        for limit in self.limits:
            sys.stderr.write("[%s] START binary heatmap at %s coverage\n"%(SimpleTime.now(), limit))
            binary_heatmap_file = "%s.binary_heatmap_limit_%sX.xls"%(self.myArgs.outputRoot, limit)
            df_binarized = pd.DataFrame(0, index=self.corsiv_info_index, columns=self.corsiv_columns_only)
            for row_idx in range(len(self.corsiv_info_index)):
                for col_idx in range(len(self.corsiv_columns_only)):
                    bp_coverage = self.df_bp_coverage.iat[row_idx, col_idx]
                    if (bp_coverage<limit):
                        df_binarized.iat[row_idx, col_idx]=0
                    else:
                        df_binarized.iat[row_idx, col_idx]=1
            df_binarized.to_csv(binary_heatmap_file, sep='\t', index_label="CoRSIV")
            sys.stderr.write("[%s] STOP binary heatmap at %s coverage\n"%(SimpleTime.now(), limit))
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CoRSIVCoverageReport::outputBinaryHeatmaps\n"%SimpleTime.now())
    
    def binarizeCoverages(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CoRSIVCoverageReport::binarizeCoverages\n"%SimpleTime.now())
        
        self.coverage_limits = [5, 10, 20, 30, 50, 100, 150, 200]
        
        data_frame_list_of_lists = []
        for column_idx in range(len(self.corsiv_columns_only)):
            sample_corsivs_array = self.df_bp_coverage.iloc[:,column_idx]
            for coverage_limit in self.coverage_limits:
                count_above_limit = np.count_nonzero(sample_corsivs_array>=coverage_limit)
                ratio_count_above_limit = float(count_above_limit)/float(len(sample_corsivs_array))
                
                if (self.DEBUG_BINARIZE_COVERAGE):
                    sys.stderr.write("[%s] Sample %s has %s (count) %s (ratio) CoRSIVs above %s\n"%(column_idx, self.corsiv_columns_only[column_idx], count_above_limit, ratio_count_above_limit, coverage_limit))
                    
                data_frame_list_of_lists.append([coverage_limit, ratio_count_above_limit])
                
        columns_binarized_coverage = ["Limit", "CoRSIV_Ratio"]
        self.binarized_coverage_data_frame = pd.DataFrame( data_frame_list_of_lists, columns = columns_binarized_coverage)
        
        binarized_coverage_data_file = "%s.binarized_coverage_data.xls"%self.myArgs.outputRoot
        self.binarized_coverage_data_frame.to_csv(binarized_coverage_data_file, sep='\t', index_label="Limit")
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("Binarized BP coverage matrix shape %s columns %s\n"%(str(self.binarized_coverage_data_frame.shape), ";".join(list(self.binarized_coverage_data_frame.columns))))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CoRSIVCoverageReport::binarizeCoverages\n"%SimpleTime.now())
        
            
    def outputViolinPlots(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CoRSIVCoverageReport::outputViolinPlots\n"%SimpleTime.now())
        
        violin_plot_pdf = "%s.violin_plot.pdf"%self.myArgs.outputRoot
        
        sn.set(style="white")
        ax = sn.violinplot(x="Limit", y="CoRSIV_Ratio", data=self.binarized_coverage_data_frame, width=1.2, linewidth=0.5)
        ax.set_ylim(0,1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=300)
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CoRSIVCoverageReport::outputViolinPlots\n"%SimpleTime.now())
    
    def loadCorsivDef(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CoRSIVCoverageReport::loadCorsivDef\n"%SimpleTime.now())
        
        self.corsiv_genomic_coord_to_name_hash={}
        
        corsiv_def_reader = open(self.myArgs.corsivDefinition, "rt")
        
        for line in corsiv_def_reader:
            ff = line.strip().split('\t')
            corsiv_name = ff[3]
            corsiv_key = "%s_%s_%s"%(ff[0], ff[1], ff[2])
            buffer = [ff[3], ff[0], ff[1], ff[2], ff[4]]
            corsiv_info = "%s"%";".join([str(x) for x in buffer])
            self.corsiv_genomic_coord_to_name_hash[corsiv_key]=corsiv_info
            if (self.DEBUG_LOAD_CORSIV):
                sys.stderr.write("[%s] CoRSIV key %s --> CoRSIV info %s %s\n"%(SimpleTime.now(), corsiv_key, corsiv_info, self.corsiv_genomic_coord_to_name_hash[corsiv_key]))
        
        corsiv_def_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CoRSIVCoverageReport::loadCorsivDef\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CoRSIVCoverageReport::work\n"%SimpleTime.now())
        
        self.loadCorsivDef()
        self.loadCoverageValues()
        
        self.binarizeCoverages()
        self.outputViolinPlots()
        
        self.outputBinaryHeatmaps()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CoRSIVCoverageReport::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = CoRSIVCoverageReport.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = CoRSIVCoverageReport(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
