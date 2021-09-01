#!/usr/bin/env python
__author__ = "Cristian Coarfa & Jianrong Dong & Kimal Rajapakshe"
__version__ = 3

import os, sys, argparse, re
import numpy as np
import scipy.stats
import pandas as pd
import seaborn as sn

import matplotlib as mpl
import matplotlib.pyplot as plt

from argparse import RawTextHelpFormatter

from CCLabUtils.simpleTime import SimpleTime

class CaptureCorsivORHeatmaps:
    DEBUG_PROGRESS = True
    DEBUGzs = False
    DEBUG_LOAD_SIGNATURES = False
    DEBUG_SETUP_ANALYSIS = False
    DEBUG_SETUP_ANALYSIS_VERBOSE = False
    DEBUG_SIG_CORRELATIONS = False
    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
  
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s


"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-d','--dataFolder', help='fdr corrected OR folder',required=True)
        parser.add_argument('-H','--heatmapHeight', help='[optional] heatmap height',required=False, default=20)
        parser.add_argument('-W','--heatmapWidth', help='[optional] heatmap width', required=False, default=10)
        parser.add_argument('-f','--rowFont', help='[optional] heatmap width', required=False, default=10)
        parser.add_argument('-o','--outputRoot', help='root of the output plots',required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        
        return args

    def doOneHeatmap(self, support_data):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCorsivORHeatmaps::doOneHeatmap %s\n"%(SimpleTime.now(), support_data))

        data_file_full = "%s/%s"%(self.myArgs.dataFolder, support_data)   
        data_df = pd.read_csv(data_file_full, header = 0, index_col=0, sep='\t')
         
        heatmap_pdf  = "%s%s%s.h_%s.w_%s_rowfont_%s.pdf"%(self.myArgs.outputRoot, os.path.sep, support_data, self.heatmap_height, self.heatmap_width, self.row_font)
        heatmap_jpg  = "%s%s%s.h_%s.w_%s_rowfont_%s.jpg"%(self.myArgs.outputRoot, os.path.sep, support_data, self.heatmap_height, self.heatmap_width, self.row_font)

        sys.stderr.write("data df %s heatmap_pdf %s heatmap_jpg %s\n"%(data_file_full, heatmap_pdf, heatmap_jpg))
        
        sn.set()
        # Draw the full plot
        g = sn.clustermap(data_df, center=0, cmap="bwr", col_cluster=False, row_cluster=False,  vmin=-3, vmax=3, linecolor='black',
            linewidths=.05, figsize=(self.heatmap_width, self.heatmap_height))
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = self.row_font)

        plt.savefig(heatmap_jpg, dpi=150)
        plt.savefig(heatmap_pdf, dpi=150)
        
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCorsivORHeatmaps::doOneHeatmap %s\n"%(SimpleTime.now(), support_data))
        
                     
                     
    def doHeatmaps(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCorsivORHeatmaps::doHeatmaps\n"%SimpleTime.now())

        self.heatmap_height = int(self.myArgs.heatmapHeight)        
        self.heatmap_width  = int(self.myArgs.heatmapWidth)        
        self.row_font  = int(self.myArgs.rowFont)        
        # heatmap of all significant correlations
        
               
        self.doOneHeatmap("odds-ratio.all-corsivs-over-all-controls.cumulative.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.cumulative.corsivs-over-subset10k-450k.20210211.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.cumulative.report-all-corsiv-vs-all-tdmr.01272021.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.cumulative.report-genic-corsiv-vs-genic-controls.02112021.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.cumulative.report-genic-corsiv-vs-genic-subset-4641-450bk-200bp.20210211.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.cumulative.report-genic-corsiv-vs-genic-tdmr.02112021.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.cumulative.report-genic-corsiv-vs-nongenic-corsivs.02112021.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.haplotype.neg_R2_over_pos_R2.cumulative.fdr_limit.0.1.xls")
        self.doOneHeatmap("odds-ratio.cumulative.report-negative-mQTL-over-positive-mQTL.02112021.fdr_limit.0.1.xls")
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCorsivORHeatmaps::doHeatmaps\n"%SimpleTime.now())
        
        
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCorsivORHeatmaps::work\n"%SimpleTime.now())
            
        self.doHeatmaps()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCorsivORHeatmaps::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
	try:
		sys.stderr.write("Command line: %s\n"%" ".join(sys.argv))
		myArgs = CaptureCorsivORHeatmaps.processArguments()
		if (myArgs is None):
			pass
		else:
			bp = CaptureCorsivORHeatmaps(myArgs)
			bp.work()
	except:
		sys.stderr.write("An unknown error occurred.\n")
		raise
