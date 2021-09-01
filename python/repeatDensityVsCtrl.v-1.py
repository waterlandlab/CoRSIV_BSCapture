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

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleStats import SimpleStats
from CCLabUtils.simpleData import SimpleData

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class FeatureRepeatOverlap:
    DEBUG_PROGRESS     = True
    DEBUG_OVERLAP      = True
    DEBUG_BACKGROUND   = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Analyze repeat overlaps within large windows or within a tiled window progression for one repeat, one control, and many features
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-r','--repeatFile',            help='UCSC repeat mask BED file (gzipped)',  required=True)
        parser.add_argument('-l','--repeatLabel',           help='repeat label',  required=True)
        parser.add_argument('-c','--controlBedFile',        help='control bed file',      required=True)
        parser.add_argument('-b','--bedFilePattern',        help='bed file(s) pattern to analyze',      required=True)
        parser.add_argument('-w','--genomicIncrement',      help='window increment size (should be divider of 100,000)',      required=True)
        parser.add_argument('-o','--outputRoot',            help='output file root',            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
    
    def getOverlap(self, file1, file2, radius):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::getOverlap\n"%SimpleTime.now())
        
        command = "bedtools window -u -w %s -a %s -b %s | wc -l"%(radius, file1, file2)
       
        overlap = subprocess.getoutput(command)

        if (self.DEBUG_OVERLAP):
            sys.stderr.write("Overlap %s within %s bp from %s is %s\n"%(file1, radius, file2, overlap))
                             
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::getOverlap\n"%SimpleTime.now())
            
        return float(overlap)


    def determineBedOverlaps(self, bed_file_info):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::determineBedOverlaps\n"%SimpleTime.now())


        if (self.DEBUG_BACKGROUND):
            sys.stderr.write("Background for %s and %s within 100k\n"%(bed_file_info.bed_file, self.myArgs.repeatFile))
        
        overlap_100k = self.getOverlap(self.myArgs.repeatFile, bed_file_info.bed_file, 100000)
        bed_file_info.cumulative_overlap_hash[100000]=overlap_100k
        
        for window_index in range(self.window_number+1):
            radius = window_index * self.window_size
            if (self.DEBUG_BACKGROUND):
                sys.stderr.write("Window idx [%s] radius %s\n"%(window_index, radius))
                
            crt_radius_overlap = self.getOverlap(self.myArgs.repeatFile, bed_file_info.bed_file, radius)
            bed_file_info.cumulative_overlap_hash[radius]=crt_radius_overlap
            if (window_index==0):
                bed_file_info.incremental_overlap_hash[radius]=crt_radius_overlap
            else:
                prev_radius = (window_index-1)*self.window_size
                prev_overlap = bed_file_info.cumulative_overlap_hash[prev_radius]
                if (self.DEBUG_BACKGROUND):
                    sys.stderr.write("Window [%s] radius %s prev radius %s prev overlap %s\n"%(window_index, radius, prev_radius, prev_overlap))
                
                bed_file_info.incremental_overlap_hash[radius]= crt_radius_overlap - prev_overlap
                
            if (self.DEBUG_BACKGROUND):
                sys.stderr.write("Window [%s] radius %s crt overlap %s incremental overlap %s\n"%(window_index, radius,
                                   crt_radius_overlap, bed_file_info.incremental_overlap_hash[radius]))


        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::determineBedOverlaps\n"%SimpleTime.now())
    
    def establishBackground(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::establishBackground\n"%SimpleTime.now())
        
        self.window_number = int(100000/self.window_size)
       
        self.control_size = SimpleData.getTextFileSize(self.myArgs.controlBedFile)
        
        self.bed_files = glob.glob(self.myArgs.bedFilePattern)
        if (len(self.bed_files)==0):
            sys.stderr.write("No bed files specified\n")
            sys.exit(4)

        self.bed_file_hash = {}        
        for bed_file in self.bed_files:
            bed_file_base = os.path.basename(bed_file)
            bed_file_size = SimpleData.getTextFileSize(bed_file)
            bed_file_info = CStruct(bed_file = bed_file, bed_file_base = bed_file_base, bed_file_size = bed_file_size,
                                    cumulative_overlap_hash = {}, incremental_overlap_hash={})
            self.bed_file_hash[bed_file_base] = bed_file_info
        
        self.control_file_base = os.path.basename(self.myArgs.controlBedFile)
        self.control_file_info = CStruct(bed_file = self.myArgs.controlBedFile, bed_file_base = self.control_file_base, bed_file_size = self.control_size, cumulative_overlap_hash = {}, incremental_overlap_hash={})
        self.determineBedOverlaps(self.control_file_info)
        
        self.control_100k_overlap = self.control_file_info.cumulative_overlap_hash[100000]
        self.control_0k_overlap = self.control_file_info.cumulative_overlap_hash[0]
        
        if (self.DEBUG_BACKGROUND):
            sys.stderr.write("Window size %s windows number %s\n"%(self.window_size, self.window_number))
        
        for bed_file_base in self.bed_file_hash:
            bed_file_info = self.bed_file_hash[bed_file_base]
            
            self.determineBedOverlaps(bed_file_info)
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::establishBackground\n"%SimpleTime.now())

    def reportBackground(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::reportBackground\n"%SimpleTime.now())
    
        background_report_0_100 = "%s.0k.100k.odds_ratios.xls"%self.myArgs.outputRoot
        background_report_0_100_writer = open(background_report_0_100, "wt")
        
        background_report_0_100_log2 = "%s.0k.100k.log2_odds_ratios.xls"%self.myArgs.outputRoot
        background_report_0_100_log2_writer = open(background_report_0_100_log2, "wt")
        
        # header = ["Repeat"]
        # buffer = [self.myArgs.repeatLabel]
       
        sorted_bed_file_bases = sorted(self.bed_file_hash.keys())
        
        background_report_0_100_writer.write("Repeat\t%s\t%s\n"%("\t".join(sorted_bed_file_bases), "\t".join(sorted_bed_file_bases)))
        background_report_0_100_log2_writer.write("Repeat\t%s\t%s\n"%("\t".join(sorted_bed_file_bases), "\t".join(sorted_bed_file_bases)))
        
        buffer_0k = []
        buffer_100k = []
        
        buffer_0k_log2 = []
        buffer_100k_log2 = []
        
        
        for bed_file_base in sorted_bed_file_bases:
            bed_file_info = self.bed_file_hash[bed_file_base]
            
            bed_file_size = float(bed_file_info.bed_file_size)
            
            overlap_0 = float(bed_file_info.cumulative_overlap_hash[0])
            odds_ratio_0k = (overlap_0/self.control_0k_overlap) / (bed_file_size/self.control_size)
            
            overlap_100k = float(bed_file_info.cumulative_overlap_hash[100000])
            odds_ratio_100k = (overlap_100k/self.control_100k_overlap) / (bed_file_size/self.control_size)
            
            buffer_0k.append(str(odds_ratio_0k))
            buffer_100k.append(str(odds_ratio_100k))
            
            if (odds_ratio_0k<2**-10):
                odds_ratio_0k_log2 = -10
            else:
                odds_ratio_0k_log2 = math.log(odds_ratio_0k)/math.log(2)
                
            buffer_0k_log2.append(str(odds_ratio_0k_log2))
            
            if (odds_ratio_100k<2**-10):
                odds_ratio_100k_log2 = -100
            else:
                odds_ratio_100k_log2 = math.log(odds_ratio_100k)/math.log(2)
                
            buffer_100k_log2.append(str(odds_ratio_100k_log2))
                       
        
        background_report_0_100_writer.write("%s\t%s\t%s\n"%(self.myArgs.repeatLabel, "\t".join(buffer_0k), "\t".join(buffer_100k)))
        background_report_0_100_writer.close()
        
        background_report_0_100_log2_writer.write("%s\t%s\t%s\n"%(self.myArgs.repeatLabel, "\t".join(buffer_0k_log2), "\t".join(buffer_100k_log2)))
        background_report_0_100_log2_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::reportBackground\n"%SimpleTime.now())
       
    
    def reportOddsRatioCumulative(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::computeIncrementalOverlaps\n"%SimpleTime.now())
        
        cumulative_report_linear = "%s.0k_100k_cumulative_odds_ratios.xls"%self.myArgs.outputRoot
        cumulative_report_linear_writer = open(cumulative_report_linear, "wt")
        cumulative_report_log2 = "%s.0k_100k_cumulative_log2_odds_ratios.xls"%self.myArgs.outputRoot
        cumulative_report_log2_writer = open(cumulative_report_log2, "wt")
        
        header = ["Repeat"]
        radius_list = []
        for window_idx in range(self.window_number+1):
            radius=int(window_idx * self.window_size)
            radius_list.append(radius)
            
        sorted_bed_file_bases = sorted(self.bed_file_hash.keys())
       
        for bed_file_base in sorted_bed_file_bases:
            for radius in radius_list:
                radius_kb = int(radius/1000)
                header.append("%s_%sKB"%(bed_file_base, radius_kb))
            
        cumulative_report_linear_writer.write("%s\n"%"\t".join(header))    
        cumulative_report_log2_writer.write("%s\n"%"\t".join(header))    
            
        buffer_log2 = [self.myArgs.repeatLabel]
        buffer_linear = [self.myArgs.repeatLabel]
        
        for bed_file_base in sorted_bed_file_bases:
            bed_file_info = self.bed_file_hash[bed_file_base]
            bed_file_size = float(bed_file_info.bed_file_size)
            
            for radius in radius_list:        
                control_radius_overlap = float(self.control_file_info.cumulative_overlap_hash[radius])
                radius_overlap = float(bed_file_info.cumulative_overlap_hash[radius])
                radius_odds_ratio = (radius_overlap/control_radius_overlap) / (bed_file_size/self.control_size)
                buffer_linear.append(str(radius_odds_ratio))
            
                
                if (radius_odds_ratio<2**-10):
                    radius_odds_ratio_log2 = -10
                else:
                    radius_odds_ratio_log2 = math.log(radius_odds_ratio)/math.log(2)
                
                buffer_log2.append(str(radius_odds_ratio_log2))
    
        cumulative_report_linear_writer.write("%s\n"%"\t".join(buffer_linear))    
        cumulative_report_log2_writer.write("%s\n"%"\t".join(buffer_log2))    
        
        cumulative_report_log2_writer.close()
        cumulative_report_linear_writer.close()
        
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::computeIncrementalOverlaps\n"%SimpleTime.now())
    
    def reportWindowDensities(self):
        pass
    
    def reportOverlaps(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::computeIncrementalOverlaps\n"%SimpleTime.now())
        
        self.reportBackground()
        self.reportOddsRatioCumulative()
        self.reportWindowDensities()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::computeIncrementalOverlaps\n"%SimpleTime.now())
            
              
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::work\n"%SimpleTime.now())
        self.window_size = int(self.myArgs.genomicIncrement)
        if (100000%self.window_size!=0):
            sys.stderr.write("Window size should divide 100,000\n")
            sys.exit(4)
        
        self.establishBackground()
        self.reportOverlaps()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = FeatureRepeatOverlap.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = FeatureRepeatOverlap(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
