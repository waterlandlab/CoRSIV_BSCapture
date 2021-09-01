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
         
        parser.add_argument('-m','--mQTLFilePattern',   help='data matrix with beta values for CORSiVs',    required=True)
        parser.add_argument('-b','--corsivBedFile',     help='corsiv bed file',         required=True)
        parser.add_argument('-o','--outputRoot',        help='output file root',        required=True)
        
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
    
    def processOneCoRSIV(self, one_mQTL_file):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::processOneCoRSIV %s"%(SimpleTime.now(), one_mQTL_file))

        corsiv_mQTL_reader = open(one_mQTL_file, "rt")
        
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
            
            buffer  = [ff[0], corsiv_short, snv_corsiv_distance, ff[2], ff[3], ff[4] ]
            self.mQTL_annotated_report_writer.write("%s\n"%"\t".join([str(x) for x in buffer]))
            focused_buffer = ["link_%s"%self.link_index, abs(snv_corsiv_distance), ff[2], abs(float(ff[2]))]
            self.link_index += 1
            self.mQTL_focused_dataFrame_writer.write("%s\n"%",".join([str(x) for x in focused_buffer]))
            
        corsiv_mQTL_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_mQTL::processOneCoRSIV %s\n"%(SimpleTime.now(), one_mQTL_file))

    def processAllCoRSIVs(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::processAllCoRSIVs\n"%SimpleTime.now())
        
        self.mQTL_file_list = glob.glob(self.myArgs.mQTLFilePattern)
        self.mQTL_annotated_report = "%s.collated_annotated.xls"%self.myArgs.outputRoot
        self.mQTL_focused_dataFrame = "%s.small_df.xls"%self.myArgs.outputRoot
        self.mQTL_focused_dataFrame_writer = open(self.mQTL_focused_dataFrame, "wt")
        
        
        self.mQTL_annotated_report_writer = open(self.mQTL_annotated_report, "wt")
        header = ["SNP", "CoRSIV", "Distance", "beta", "tstat", "pvalue"]
        self.mQTL_annotated_report_writer.write("%s\n"%"\t".join(header))
        
        focused_header = ["CoRSIV_SNV_Link", "AbsDistance", "beta", "absBeta"]
        self.mQTL_focused_dataFrame_writer.write("%s\n"%",".join(focused_header))

        self.link_index = 0        
        for one_mQTL_file in self.mQTL_file_list:
            self.processOneCoRSIV(one_mQTL_file)        
        
        self.mQTL_focused_dataFrame_writer.close()
        self.mQTL_annotated_report_writer.close()
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::processAllCoRSIVs\n"%SimpleTime.now())
        
    
    def betaCoefficientByDistance(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::betaCoefficientByDistance\n"%SimpleTime.now())
        
        focused_data_frame = pd.read_csv(self.mQTL_focused_dataFrame, header=0, index_col=0)
        
        mean_median_beta_distance_report = "%s.betaCoeff_by_distance.xls"%self.myArgs.outputRoot
        mean_median_beta_distance_report_writer = open(mean_median_beta_distance_report, "wt")
        
        summary_header = ["Dummy_Index", "Window", "Min_Distance", "Max_Distance", "beta_avg", "beta_median", "absBeta_avg", "absBeta_median", "links"]
        mean_median_beta_distance_report_writer.write("%s\n"%",".join(summary_header))
        
        for window_index in range(1000):
            min_distance = window_index * 1000
            max_distance = window_index * 1000 + 999
                
            link_slice = focused_data_frame[(focused_data_frame.AbsDistance >= min_distance) & (focused_data_frame.AbsDistance<=max_distance)]
            
            if (self.DEBUG_BETA_BY_DISTANCE):
                sys.stderr.write("W %s min_d %s max_d %s slice shape %s\n"%(window_index, min_distance, max_distance, link_slice.shape))
                
            if (link_slice.shape[0]>0):
                beta_array = link_slice.iloc[:,1]
                abs_beta_array = link_slice.iloc[:,2]
                beta_avg = np.mean(beta_array)
                beta_median = np.median(beta_array)
                abs_beta_avg = np.mean(abs_beta_array)
                abs_beta_median = np.median(abs_beta_array)
                
                summary_buffer = ["w_%s"%window_index, window_index, min_distance, max_distance, beta_avg, beta_median, abs_beta_avg, abs_beta_median, link_slice.shape[0] ]

                mean_median_beta_distance_report_writer.write("%s\n"%",".join([str(x) for x in summary_buffer]))                
            
        #    meth_slice = chrom_specific_df[(chrom_specific_df.index>=chrom_start) & (chrom_specific_df.index<=chrom_stop) & (chrom_specific_df.N>=self.min_read_coverage)]
        
        
        mean_median_beta_distance_report_writer.close()
        
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_mQTL::betaCoefficientByDistance\n"%SimpleTime.now())
        
    
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::work\n"%SimpleTime.now())
        
        self.loadCoRSIVDef()
        self.processAllCoRSIVs()
        self.betaCoefficientByDistance()
        
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
