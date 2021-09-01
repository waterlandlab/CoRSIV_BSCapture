#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

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
        
class GenomicFeatureDensityPlot:
    DEBUG_PROGRESS                  = True
    DEBUG_WINDOW_PREP               = True
    DEBUG_RUN_WINDOW                = True
    DEBUG_LOAD_CORSIV               = False
    DEBUG_BINARIZE_COVERAGE         = False

    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Generate a density plot of neighboring features around two BED files and plot comparative density
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-b','--bedFile1',              help='bed file 1',   required=True)
        parser.add_argument('-B','--bedFile2',              help='bed file 2',   required=True)
        parser.add_argument('-l','--labelFile1',            help='bed file 1',   required=True)
        parser.add_argument('-L','--labelFile2',            help='bed file 2',   required=True)
        parser.add_argument('-g','--genomicFeaturesBed',    help='genomic features BED file', required=True)
        parser.add_argument('-G','--labelGenomicFeatures',  help='genomic features BED file', required=True)
        parser.add_argument('-C','--chromosomeMap',         help='TAB delimited chromosome name and size', required=True)
        parser.add_argument('-r','--radius',                help='radius around bed file',    required=False, default=100000)
        parser.add_argument('-w','--windowCount',           help='number of windows to split radius in',    required=False, default=100)
        parser.add_argument('-t','--plotTitle',             help='feature description',                     required=True)
        parser.add_argument('-k','--keepTemporaryFiles',    help='[optional] keep intermediate files (default=false)',    required=False, action="store_true")
        parser.add_argument('-o','--outputRoot',            help='output file root',                        required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args


    def prepareWindowRadiusBedFile(self, bed_file):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START GenomicFeatureDensityPlot::prepareWindowRadiusBedFile\n"%SimpleTime.now())
    
        # window_radius_file_name = "%s.radius_%s_windowCount_%s_around_%s.bed"%(self.myArgs.outputRoot, self.myArgs.radius, self.myArgs.windowCount, os.path.basename(bed_file))
        window_radius_file_name = "%s.radius_%s"%(self.myArgs.outputRoot, os.path.basename(bed_file))
        self.cleanup_list.append(window_radius_file_name)
       
        if (self.DEBUG_WINDOW_PREP):
            sys.stderr.write("window_radius_file_name %s\n"%window_radius_file_name)
        
        radius_size = int(self.myArgs.radius)
        window_count = int(self.myArgs.windowCount)
        self.window_size = radius_size/window_count
        window_hash = {}
        
        window_radius_file_writer = open(window_radius_file_name, "wt")
        
        bed_reader = open(bed_file, "rt")
        bed_idx = -1
        for line in bed_reader:
            bed_idx += 1
            ff = line.strip().split("\t")
            chrom = ff[0]
            start = int(ff[1])
            stop = int(ff[2])
            mid_point = int((start+stop)/2)
            if (self.DEBUG_WINDOW_PREP):
                sys.stderr.write("[bed entry %s] chrom %s start %s stop %s midpoint %s\n"%(bed_idx, chrom, start, stop, mid_point))
            
            for window_idx in range(0, window_count):
                in_radius_window_plus_index =  window_count + window_idx
                window_plus_start = int(mid_point + window_idx*self.window_size)
                window_plus_stop = int(mid_point + window_idx*self.window_size+self.window_size-1)
                window_radius_file_writer.write("%s\t%s\t%s\tfeature_%s.window_%s\t1\t+\n"%(chrom, window_plus_start, window_plus_stop, bed_idx,in_radius_window_plus_index))
                
                key = "%s_%s_%s"%(chrom, window_plus_start, window_plus_stop)
                window_hash[key]=in_radius_window_plus_index    
                
                in_radius_window_minus_index = window_count-1-window_idx
                window_minus_stop = int(mid_point - window_idx*self.window_size)
                window_minus_start = int(mid_point - window_idx*self.window_size  - self.window_size + 1)
                if (window_minus_start<1 or window_minus_stop<1):
                    continue
                window_radius_file_writer.write("%s\t%s\t%s\tfeature_%s.window_%s\t1\t+\n"%(chrom, window_minus_start, window_minus_stop, bed_idx, in_radius_window_minus_index))
                
                key = "%s_%s_%s"%(chrom, window_minus_start, window_minus_stop)
                window_hash[key]=in_radius_window_minus_index    
        
        bed_reader.close()
        window_radius_file_writer.close()
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP GenomicFeatureDensityPlot::prepareWindowRadiusBedFile\n"%SimpleTime.now())
        
        return window_radius_file_name, window_hash
    
    def getBedFileSize(self, bed_file):
        bed_reader = open(bed_file, "rt")
        bed_count = 0
        for line in bed_reader:
            bed_count +=1
            
        return bed_count
        
    def runWindowRadiusCoverage(self, window_radius_file, bed_file, window_hash):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START GenomicFeatureDensityPlot::runWindowRadiusCoverage\n"%SimpleTime.now())
    
        # eRNA_command = k
        feature_density = 13
        bed_count = self.getBedFileSize(bed_file)
        if (self.DEBUG_RUN_WINDOW):
            sys.stderr.write("bed_count for %s is %s \n"%(bed_file, bed_count))
        
        tmp_folder = "%s_tmp_%s_eRNA"%(self.myArgs.outputRoot, os.path.basename(window_radius_file))
        
#        eRNA_scores = "%s_eRNA_%s_counts"%(self.myArgs.outputRoot, os.path.basename(window_radius_file))
        eRNA_scores = "%s_%s_cnt"%(self.myArgs.outputRoot, os.path.basename(window_radius_file))
       
       
        eRNA_log = "%s_le_%s.txt"%(self.myArgs.outputRoot, os.path.basename(window_radius_file))
        self.cleanup_list.append(eRNA_log)
        
        eRNA_command = "eRNAScore.cc.py -r %s -e %s -x %s -p 2 -C %s -o %s &> %s"%(self.myArgs.genomicFeaturesBed,
                window_radius_file, tmp_folder, self.myArgs.chromosomeMap, eRNA_scores, eRNA_log)
        
        
        if (self.DEBUG_RUN_WINDOW):
            sys.stderr.write("eRNA command %s\n"%eRNA_command)
        
        os.system(eRNA_command)
        
        # summarize coverage
        window_count = int(self.myArgs.windowCount)
        coverage_sum = np.zeros(2*window_count)
        coverage_avg = np.zeros(2*window_count)
        
        coverage_info_xls = "%s.coverageInfo.xls"%eRNA_scores
        self.cleanup_list.append(coverage_info_xls)
        
        eRNA_reader = open(coverage_info_xls, "rt")
        
        line_idx = -1
        for line in eRNA_reader:
            line_idx += 1
            if (line_idx==0):
                continue
            
            ff = line.strip().split("\t")
            coverage = float(ff[2])
            window_key = ff[0]
            window_idx_in_radius = window_hash[window_key]
            coverage_sum[window_idx_in_radius] += coverage
            
        
        eRNA_reader.close()
        
        for idx in range(2*window_count):
            coverage_avg[idx] = coverage_sum[idx]/float(bed_count)
            
        if not self.myArgs.keepTemporaryFiles:
            rm_command = "rm -rf %s"%tmp_folder
            sys.stderr.write("rm command : %s\n"%rm_command)
            os.system(rm_command)
       
        feature_density = CStruct(bed_regions = bed_count, coverage_sum = coverage_sum, coverage_avg = coverage_avg)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP GenomicFeatureDensityPlot::runWindowRadiusCoverage\n"%SimpleTime.now())
        
        return feature_density
    
    def computeFeatureDensity(self, bed_file):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START GenomicFeatureDensityPlot::computeFeatureDensity\n"%SimpleTime.now())
        
        window_radius_file, window_hash = self.prepareWindowRadiusBedFile(bed_file)
        feature_density = self.runWindowRadiusCoverage(window_radius_file, bed_file, window_hash)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP GenomicFeatureDensityPlot::computeFeatureDensity\n"%SimpleTime.now())
        
        return feature_density        
    
    def outputCombinedFeatureDensityPlot(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START GenomicFeatureDensityPlot::outputCombinedFeatureDensityPlot\n"%SimpleTime.now())
        
       
        support_data = "%s.support_data.xls"%self.myArgs.outputRoot
        support_data_writer = open(support_data, "wt")
        
        support_data_writer.write("Distance from BED center\t%s_%s\t%s_%s\n"%(self.myArgs.labelGenomicFeatures, self.myArgs.labelFile1,
                                                                              self.myArgs.labelGenomicFeatures, self.myArgs.labelFile2))
        
        window_count = int(self.myArgs.windowCount)
        x_axis = np.zeros(2*window_count)
        for idx in range(2*window_count):
            x_axis[idx]=-window_count*self.window_size+idx*self.window_size
            support_data_writer.write("%s\t%s\t%s\n"%(x_axis[idx], self.density_around_bed_file1.coverage_avg[idx], self.density_around_bed_file2.coverage_avg[idx]))
            
        support_data_writer.close()
        
        support_data_writer.close()
        
        simple_plot_pdf = "%s.line_plot.pdf"%self.myArgs.outputRoot
        plt.figure()
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        ax.set_xlabel("Distance from BED midpoint") 
        ax.set_ylabel("Feature Density")
        ax.set_title("Density of %s within a %s-bp radius"%(self.myArgs.plotTitle, self.myArgs.radius))
        
        ax.plot(x_axis, self.density_around_bed_file1.coverage_avg, '-', color='b')
        ax.plot(x_axis, self.density_around_bed_file2.coverage_avg, '-', color='r')
        
        plt.savefig(simple_plot_pdf, dpi=150) 
        
        # two panel plot
        two_panel_plot_pdf = "%s.two_panel_barplot.pdf"%self.myArgs.outputRoot
        two_panel_plot_jpg = "%s.two_panel_barplot.jpg"%self.myArgs.outputRoot
        
       
        max_coverage_density_1 = max(self.density_around_bed_file1.coverage_avg)
        max_coverage_density_2 = max(self.density_around_bed_file2.coverage_avg)
        max_12 = max(max_coverage_density_1, max_coverage_density_2)
        
        sys.stderr.write("Max1 %s max2 %s max_12 %s\n"%(max_coverage_density_1, max_coverage_density_2, max_12))
        
        plt.figure()
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        
        xl = np.arange(2*window_count)
        
        ax1.bar(xl, self.density_around_bed_file1.coverage_avg, width=0.9, color='white', edgecolor='black')
        # ax1.set_xticks([xl[0], xl[window_count-1], xl[2*window_count-1]])
        # ax1.set_xticklabels([x_axis[0], x_axis[window_count-1], x_axis[2*window_count-1]])
        # ax1.set_xticks(xl)
        ax1.set_xticks([xl[0], xl[window_count-1], xl[2*window_count-1]])
        ax1.set_xticklabels([-100000, 0, 100000])
        
        ax1.set_xlabel("Distance from BED midpoint") 
        ax1.set_ylabel("Feature Density")
#        ax1.set_title("Density of %s within a %s-bp radius"%(self.myArgs.plotTitle, self.myArgs.radius))
        ax1.set_title("Density of %s"%self.myArgs.plotTitle)
        ax1.set_ylim(0, max_12)
        
        ax2.bar(xl, self.density_around_bed_file2.coverage_avg, width=0.9, color='white', edgecolor='black')
        ax2.set_xticks([xl[0], xl[window_count-1], xl[2*window_count-1]])
        ax2.set_xticklabels([-100000, 0, 100000])
        
        ax2.set_xlabel("Distance from BED midpoint") 
        ax2.set_ylabel("Feature Density")
        ax2.set_title("Density of %s"%self.myArgs.plotTitle)
#        ax2.set_title("Density of %s within a %s-bp radius"%(self.myArgs.plotTitle, self.myArgs.radius))
        ax2.set_ylim(0, max_12)
        
        plt.savefig(two_panel_plot_jpg, dpi=150)
        plt.savefig(two_panel_plot_pdf, dpi=150)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START GenomicFeatureDensityPlot::outputCombinedFeatureDensityPlot\n"%SimpleTime.now())

    def cleanup(self):
        if not self.myArgs.keepTemporaryFiles:
            for tmp_file in self.cleanup_list:
                sys.stderr.write("cleaning up %s\n"%tmp_file)
                os.unlink(tmp_file)
                
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START GenomicFeatureDensityPlot::work\n"%SimpleTime.now())
        
        self.cleanup_list = []
        self.density_around_bed_file1 = self.computeFeatureDensity(self.myArgs.bedFile1)
        self.density_around_bed_file2 = self.computeFeatureDensity(self.myArgs.bedFile2)
        
        self.outputCombinedFeatureDensityPlot()
        self.cleanup()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP GenomicFeatureDensityPlot::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = GenomicFeatureDensityPlot.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = GenomicFeatureDensityPlot(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
