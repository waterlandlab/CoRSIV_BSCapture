#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "4.0"

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

class FeatureRepeatOverlap:
    DEBUG_PROGRESS     = True
    DEBUG_OVERLAP      = True
    DEBUG_BACKGROUND   = True
    DEBUG_div0   = True
    DEBUG_ODDS_RATIO_DIRECTION  = True
    DEBUG_COMPARISON_DIRECTION  = True

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
        parser.add_argument('-W','--genomicWindow',         help='max genomic window around corsivs (default 100,000)',  type=int,   required=False, default=100000)
        parser.add_argument('-w','--genomicIncrement',      help='window increment size (should be divider of the genomic window size)',      required=True)
        parser.add_argument('-C','--oddsRatioComparisons',  help='file containing pairs of comparisons',   required=False)
        parser.add_argument('-k','--keepTempFiles',         help='keep temporary files (note: that is a lot of files)',   required=False, action="store_true")
        parser.add_argument('-o','--outputRoot',            help='output file root',  required=True)

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

    def getUpstreamOverlap(self, file1, file2, radius):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::getUpstreamOverlap %s %s %s\n"%(SimpleTime.now(), file1, file2, radius))

        # establish a bait file upstream with radius
        tmp_upstream_file = "%s.tmp_upstream.%s.by.%s.bed"%(self.myArgs.outputRoot, os.path.basename(file2), radius)
        tmp_upstream_file_writer = open(tmp_upstream_file, "wt")
        with open(file2, "rt") as file2_reader:
            for line in file2_reader:
                ff = line.strip().split('\t')
                bed_start = int(ff[1])
                upstream_start = bed_start - radius
                if (upstream_start<1):
                    upstream_start = 1
                ff[1] = str(upstream_start)
                tmp_upstream_file_writer.write("%s\n"%"\t".join(ff))

        tmp_upstream_file_writer.close()

        command = "bedtools intersect -u  -a %s -b %s | wc -l"%(file1, tmp_upstream_file)

        overlap = subprocess.getoutput(command)

        if (self.DEBUG_OVERLAP):
            sys.stderr.write("Upstream overlap %s within %s bp from %s is %s\n"%(file1, radius, file2, overlap))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::getUpstreamOverlap %s %s %s\n"%(SimpleTime.now(), file1, file2, radius))

        if not self.myArgs.keepTempFiles:
            os.unlink(tmp_upstream_file)

        return float(overlap)

    def getDownstreamOverlap(self, file1, file2, radius):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::getDownstreamOverlap %s %s %s\n"%(SimpleTime.now(), file1, file2, radius))

        # establish a bait file upstream with radius
        tmp_downstream_file = "%s.tmp_downstream.%s.by.%s.bed"%(self.myArgs.outputRoot, os.path.basename(file2), radius)
        tmp_downstream_file_writer = open(tmp_downstream_file, "wt")
        with open(file2, "rt") as file2_reader:
            for line in file2_reader:
                ff = line.strip().split('\t')
                bed_stop = int(ff[2])
                downstream_stop = bed_stop + radius
                ff[2] = str(downstream_stop)
                tmp_downstream_file_writer.write("%s\n"%"\t".join(ff))

        tmp_downstream_file_writer.close()

        command = "bedtools intersect -u  -a %s -b %s | wc -l"%(file1, tmp_downstream_file)

        overlap = subprocess.getoutput(command)

        if (self.DEBUG_OVERLAP):
            sys.stderr.write("Downstream overlap %s within %s bp from %s is %s\n"%(file1, radius, file2, overlap))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::getDownstreamOverlap %s %s %s\n"%(SimpleTime.now(), file1, file2, radius))

        if not self.myArgs.keepTempFiles:
            os.unlink(tmp_downstream_file)

        return float(overlap)

    def determineBedOverlaps(self, bed_file_info):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::determineBedOverlaps\n"%SimpleTime.now())


        if (self.DEBUG_BACKGROUND):
            sys.stderr.write("Background for %s and %s within 100k\n"%(bed_file_info.bed_file, self.myArgs.repeatFile))

        overlap_genomic_window = self.getOverlap(self.myArgs.repeatFile, bed_file_info.bed_file, self.myArgs.genomicWindow)
        bed_file_info.cumulative_overlap_hash[self.genomic_window]=overlap_genomic_window

        for window_index in range(self.window_number+1):
            radius = window_index * self.genomic_increment
            if (self.DEBUG_BACKGROUND):
                sys.stderr.write("Window idx [%s] radius %s\n"%(window_index, radius))

            crt_radius_overlap = self.getOverlap(self.myArgs.repeatFile, bed_file_info.bed_file, radius)
            bed_file_info.cumulative_overlap_hash[radius]=crt_radius_overlap
            if (window_index==0):
                bed_file_info.incremental_overlap_hash[radius]=crt_radius_overlap
                bed_file_info.directional_overlap_hash[radius]=crt_radius_overlap
            else:
                # determine the upstream overlap
                bed_file_info.directional_overlap_hash[-radius]=self.getUpstreamOverlap(self.myArgs.repeatFile, bed_file_info.bed_file, radius)
                # determine the downstream overlap
                bed_file_info.directional_overlap_hash[radius]=self.getDownstreamOverlap(self.myArgs.repeatFile, bed_file_info.bed_file, radius)
                prev_radius = (window_index-1)*self.genomic_increment
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

        self.window_number = int(self.genomic_window/self.genomic_increment)

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
                                    cumulative_overlap_hash = {}, incremental_overlap_hash={}, directional_overlap_hash={})
            self.bed_file_hash[bed_file_base] = bed_file_info

        self.control_file_base = os.path.basename(self.myArgs.controlBedFile)
        self.control_file_info = CStruct(bed_file = self.myArgs.controlBedFile, bed_file_base = self.control_file_base,
                            bed_file_size = self.control_size, cumulative_overlap_hash = {}, incremental_overlap_hash={}, directional_overlap_hash={})
        self.determineBedOverlaps(self.control_file_info)

        self.control_genomic_window_overlap = self.control_file_info.cumulative_overlap_hash[self.genomic_window]
        self.control_0k_overlap = self.control_file_info.cumulative_overlap_hash[0]

        if (self.DEBUG_BACKGROUND):
            sys.stderr.write("Window size %s windows number %s\n"%(self.genomic_increment, self.window_number))

        for bed_file_base in self.bed_file_hash:
            bed_file_info = self.bed_file_hash[bed_file_base]

            self.determineBedOverlaps(bed_file_info)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::establishBackground\n"%SimpleTime.now())

    def reportBackground(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::reportBackground\n"%SimpleTime.now())

        background_report_0_100 = "%s.0k.%sk.odds_ratios.xls"%(self.myArgs.outputRoot, int(self.genomic_window/1000))
        background_report_0_100_writer = open(background_report_0_100, "wt")

        background_report_0_100_log2 = "%s.0k.%sk.log2_odds_ratios.xls"%(self.myArgs.outputRoot, int(self.genomic_window/1000))
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
            if (self.DEBUG_div0):
            	sys.stderr.write("bed_file_size %s self.control_size %s overlap_0 %s self.control_0k_overlap %s\n"%(bed_file_size, self.control_size, overlap_0, self.control_0k_overlap))

            if (self.control_0k_overlap>0):
                odds_ratio_0k = (overlap_0/self.control_0k_overlap) / (bed_file_size/self.control_size)
            else:
                if (overlap_0>0):
                    odds_ratio_0k = 2**10
                else:
                    odds_ratio_0k = 1


            overlap_genomic_window = float(bed_file_info.cumulative_overlap_hash[self.genomic_window])
            odds_ratio_100k = (overlap_genomic_window/self.control_genomic_window_overlap) / (bed_file_size/self.control_size)

            buffer_0k.append(str(odds_ratio_0k))
            buffer_100k.append(str(odds_ratio_100k))

            if (odds_ratio_0k<2**-10):
                odds_ratio_0k_log2 = -10
            else:
                odds_ratio_0k_log2 = math.log(odds_ratio_0k)/math.log(2)

            buffer_0k_log2.append(str(odds_ratio_0k_log2))

            if (odds_ratio_100k<2**-10):
                odds_ratio_100k_log2 = -10
            else:
                odds_ratio_100k_log2 = math.log(odds_ratio_100k)/math.log(2)

            buffer_100k_log2.append(str(odds_ratio_100k_log2))


        background_report_0_100_writer.write("%s\t%s\t%s\n"%(self.myArgs.repeatLabel, "\t".join(buffer_0k), "\t".join(buffer_100k)))
        background_report_0_100_writer.close()

        background_report_0_100_log2_writer.write("%s\t%s\t%s\n"%(self.myArgs.repeatLabel, "\t".join(buffer_0k_log2), "\t".join(buffer_100k_log2)))
        background_report_0_100_log2_writer.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::reportBackground\n"%SimpleTime.now())


    def reportOddsRatioCumulative(self, cumulativeFlag):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::computeIncrementalOverlaps cumulativeFlag %s\n"%(SimpleTime.now(), cumulativeFlag))

        if (cumulativeFlag):
            string_label = "cumulative"
        else:
            string_label = "incremental"

        cumulative_report_linear = "%s.0k_%sk_%s_odds_ratios.xls"%(self.myArgs.outputRoot, int(self.genomic_window/1000), string_label)
        cumulative_report_linear_writer = open(cumulative_report_linear, "wt")
        cumulative_report_log2 = "%s.0k_%sk_%s_log2_odds_ratios.xls"%(self.myArgs.outputRoot, int(self.genomic_window/1000), string_label)
        cumulative_report_log2_writer = open(cumulative_report_log2, "wt")

        header = ["Repeat"]
        radius_list = []
        for window_idx in range(self.window_number+1):
            radius=int(window_idx * self.genomic_increment)
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

            # also produce a per-file overlap with log2 odds ratio and p-value
            individual_bed_overlap_report = "%s.individual_overlap_wrt_control.%s.%s.xls"%(self.myArgs.outputRoot, bed_file_base, string_label)
            individual_bed_overlap_report_writer = open(individual_bed_overlap_report, "wt")
            individual_bed_overlap_report_header = ["Repeat"]

            individual_header_dist = []
            individual_header_pvalue = []

            for radius in radius_list:
                radius_kb = int(radius/1000)
                individual_header_dist.append("%s_%sKB_OR"%(bed_file_base, radius_kb))
                individual_header_pvalue.append("%s_%sKB_Pval"%(bed_file_base, radius_kb))

            individual_bed_overlap_report_writer.write("Repeat\t%s\t%s\n"%("\t".join(individual_header_dist), "\t".join(individual_header_pvalue)))

            ind_buffer_or = []
            ind_buffer_pv = []

            for radius in radius_list:
                if (cumulativeFlag):
                    control_radius_overlap = float(self.control_file_info.cumulative_overlap_hash[radius])
                    radius_overlap = float(bed_file_info.cumulative_overlap_hash[radius])
                else:
                    control_radius_overlap = float(self.control_file_info.incremental_overlap_hash[radius])
                    radius_overlap = float(bed_file_info.incremental_overlap_hash[radius])



                if (control_radius_overlap>0):
                    radius_odds_ratio = (radius_overlap/control_radius_overlap) / (bed_file_size/self.control_size)
                else:
                    if (radius_overlap>0):
                        radius_odds_ratio = 2**10
                    else:
                        radius_odds_ratio = 1

#                radius_odds_ratio = (radius_overlap/control_radius_overlap) / (bed_file_size/self.control_size)
                buffer_linear.append(str(radius_odds_ratio))


                if (radius_odds_ratio<2**-10):
                    radius_odds_ratio_log2 = -10
                else:
                    radius_odds_ratio_log2 = math.log(radius_odds_ratio)/math.log(2)

                oddsRatio, pValue = scipy.stats.fisher_exact( [[ radius_overlap, control_radius_overlap], [bed_file_size, self.control_size]])
                ind_buffer_or.append(oddsRatio)
                ind_buffer_pv.append(pValue)

                buffer_log2.append(str(radius_odds_ratio_log2))

            individual_bed_overlap_report_writer.write("%s\t%s\t%s"%(self.myArgs.repeatLabel, "\t".join([str(x) for x in ind_buffer_or]),
                                                                     "\t".join([str(x) for x in ind_buffer_pv])))
            individual_bed_overlap_report_writer.close()

        cumulative_report_linear_writer.write("%s\n"%"\t".join(buffer_linear))
        cumulative_report_log2_writer.write("%s\n"%"\t".join(buffer_log2))

        cumulative_report_log2_writer.close()
        cumulative_report_linear_writer.close()


        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::computeIncrementalOverlaps cumulativeFlag %s\n"%(SimpleTime.now(), cumulativeFlag))

    def reportOddsRatioWithDirection(self): # report only cumulative overlaps
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::reportOddsRatioWithDirection \n"%SimpleTime.now())

        radius_list = []
        for window_idx in range(-self.window_number,0):
            radius=int(window_idx * self.genomic_increment)
            radius_list.append(radius)

        for window_idx in range(self.window_number+1):
            radius=int(window_idx * self.genomic_increment)
            radius_list.append(radius)

        if (self.DEBUG_ODDS_RATIO_DIRECTION):
            sys.stderr.write("radius list: %s\n"%' '.join([str(x) for x in radius_list]))

        sorted_bed_file_bases = sorted(self.bed_file_hash.keys())

        for bed_file_base in sorted_bed_file_bases:
            if (self.DEBUG_ODDS_RATIO_DIRECTION):
                sys.stderr.write("computing cumulative odds-ratio/p-value for %s vs control %s\n"%(bed_file_base, self.control_file_base))

            bed_file_info = self.bed_file_hash[bed_file_base]
            bed_file_size = float(bed_file_info.bed_file_size)

            # also produce a per-file overlap with linear odds ratio (computed by fisher test) and p-value
            individual_bed_overlap_report = "%s.individual_overlap_wrt_control.%s.direction.cumulative.xls"%(self.myArgs.outputRoot,bed_file_base)
            individual_bed_overlap_report_writer = open(individual_bed_overlap_report, "wt")
            individual_bed_overlap_report_header = ["Repeat"]

            individual_header_dist = []
            individual_header_pvalue = []

            for radius in radius_list:
                radius_kb = int(radius/1000)
                individual_header_dist.append("%s_%sKB_OR"%(bed_file_base, radius_kb))
                individual_header_pvalue.append("%s_%sKB_Pval"%(bed_file_base, radius_kb))

            individual_bed_overlap_report_writer.write("Repeat\t%s\t%s\n"%("\t".join(individual_header_dist), "\t".join(individual_header_pvalue)))

            ind_buffer_or = []
            ind_buffer_pv = []

            for radius in radius_list:
                # if (cumulativeFlag): only do cumulative
                control_radius_overlap = float(self.control_file_info.directional_overlap_hash[radius])
                radius_overlap = float(bed_file_info.directional_overlap_hash[radius])

                oddsRatio, pValue = scipy.stats.fisher_exact( [[ radius_overlap, control_radius_overlap], [bed_file_size, self.control_size]])
                ind_buffer_or.append(oddsRatio)
                ind_buffer_pv.append(pValue)

            individual_bed_overlap_report_writer.write("%s\t%s\t%s"%(self.myArgs.repeatLabel,
                                    "\t".join([str(x) for x in ind_buffer_or]),
                                    "\t".join([str(x) for x in ind_buffer_pv])))

            individual_bed_overlap_report_writer.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::reportOddsRatioWithDirection \n"%SimpleTime.now())


    def performComparisons(self, cumulativeFlag):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::performComparisons cumulativeFlag %s\n"%(SimpleTime.now(), cumulativeFlag))

        if (cumulativeFlag):
            string_label = "cumulative"
        else:
            string_label = "incremental"

        radius_list = []
        for window_idx in range(self.window_number+1):
            radius=int(window_idx * self.genomic_increment)
            radius_list.append(radius)

        # read comparison pairs
        comparison_reader = open(self.myArgs.oddsRatioComparisons, "rt")

        for line in comparison_reader:
            ff = line.strip().split()

            comparison_label = ff[0]
            bed_file1 = ff[1]
            bed_file2 = ff[2]

            # also produce a per-file overlap with log2 odds ratio and p-value
            one_comparison_report = "%s.comparison.%s.%s.xls"%(self.myArgs.outputRoot, comparison_label, string_label)
            one_comparison_report_writer = open(one_comparison_report, "wt")

            individual_header_dist = []
            individual_header_pvalue = []

            for radius in radius_list:
                radius_kb = int(radius/1000)
                individual_header_dist.append("%s_%sKB_OR"%(comparison_label, radius_kb))
                individual_header_pvalue.append("%s_%sKB_Pval"%(comparison_label, radius_kb))

            one_comparison_report_writer.write("Repeat\t%s\t%s\n"%("\t".join(individual_header_dist), "\t".join(individual_header_pvalue)))

            ind_buffer_or = []
            ind_buffer_pv = []

            group1_info = self.bed_file_hash[bed_file1]
            group2_info = self.bed_file_hash[bed_file2]

            bed_file1_size = float(group1_info.bed_file_size)
            bed_file2_size = float(group2_info.bed_file_size)

            for radius in radius_list:
                if (cumulativeFlag):
                    bed_file1_overlap = float(group1_info.cumulative_overlap_hash[radius])
                    bed_file2_overlap = float(group2_info.cumulative_overlap_hash[radius])
                else:
                    bed_file1_overlap = float(group1_info.incremental_overlap_hash[radius])
                    bed_file2_overlap = float(group2_info.incremental_overlap_hash[radius])

                oddsRatio, pValue = scipy.stats.fisher_exact( [[ bed_file2_overlap, bed_file1_overlap],
                                                               [ bed_file2_size, bed_file1_size ]])
                # adjust for "infinite" odds ratio
                if (self.DEBUG_div0):
                    sys.stderr.write("comparison: bed_file2_overlap %s bed_file2_size %s bed_file1_overlap %s bed_file1_size %s or %s p %s\n"%
                                (bed_file2_overlap, bed_file2_size, bed_file1_overlap, bed_file1_size, oddsRatio, pValue))

                ind_buffer_or.append(oddsRatio)
                ind_buffer_pv.append(pValue)


            one_comparison_report_writer.write("%s\t%s\t%s"%(self.myArgs.repeatLabel, "\t".join([str(x) for x in ind_buffer_or]),
                                                                     "\t".join([str(x) for x in ind_buffer_pv])))
            one_comparison_report_writer.close()

        comparison_reader.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::performComparisons cumulativeFlag %s\n"%(SimpleTime.now(), cumulativeFlag))

    def performComparisonsWithDirection(self): # only do cumulative
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::performComparisonsWithDirection \n"%SimpleTime.now())

        radius_list = []
        for window_idx in range(-self.window_number,0):
            radius=int(window_idx * self.genomic_increment)
            radius_list.append(radius)

        for window_idx in range(self.window_number+1):
            radius=int(window_idx * self.genomic_increment)
            radius_list.append(radius)

        if (self.DEBUG_ODDS_RATIO_DIRECTION):
            sys.stderr.write("radius list: %s\n"%' '.join([str(x) for x in radius_list]))

        # read comparison pairs
        comparison_reader = open(self.myArgs.oddsRatioComparisons, "rt")

        for line in comparison_reader:
            ff = line.strip().split()

            comparison_label = ff[0]
            bed_file1 = ff[1]
            bed_file2 = ff[2]

            if (self.DEBUG_COMPARISON_DIRECTION):
                sys.stderr.write("comparison cumulative odds-ratio/p-value for %s: %s vs %s\n"%(comparison_label, bed_file1, bed_file2))

            # also produce a per-file overlap with log2 odds ratio and p-value
            one_comparison_report = "%s.comparison.%s.with_direction.cumulative.xls"%(self.myArgs.outputRoot, comparison_label)
            one_comparison_report_writer = open(one_comparison_report, "wt")

            individual_header_dist = []
            individual_header_pvalue = []

            for radius in radius_list:
                radius_kb = int(radius/1000)
                individual_header_dist.append("%s_%sKB_OR"%(comparison_label, radius_kb))
                individual_header_pvalue.append("%s_%sKB_Pval"%(comparison_label, radius_kb))

            one_comparison_report_writer.write("Repeat\t%s\t%s\n"%("\t".join(individual_header_dist), "\t".join(individual_header_pvalue)))

            ind_buffer_or = []
            ind_buffer_pv = []

            group1_info = self.bed_file_hash[bed_file1]
            group2_info = self.bed_file_hash[bed_file2]

            bed_file1_size = float(group1_info.bed_file_size)
            bed_file2_size = float(group2_info.bed_file_size)

            for radius in radius_list:
                # if (cumulativeFlag): only do cumulative
                bed_file1_overlap = float(group1_info.directional_overlap_hash[radius])
                bed_file2_overlap = float(group2_info.directional_overlap_hash[radius])

                oddsRatio, pValue = scipy.stats.fisher_exact( [[ bed_file2_overlap, bed_file1_overlap],
                                                               [ bed_file2_size, bed_file1_size ]])
                # adjust for "infinite" odds ratio
                if (self.DEBUG_COMPARISON_DIRECTION):
                    sys.stderr.write("direction x comparison: bed_file2_overlap %s bed_file2_size %s bed_file1_overlap %s bed_file1_size %s or %s p %s\n"%
                                (bed_file2_overlap, bed_file2_size, bed_file1_overlap, bed_file1_size, oddsRatio, pValue))

                ind_buffer_or.append(oddsRatio)
                ind_buffer_pv.append(pValue)


            one_comparison_report_writer.write("%s\t%s\t%s"%(self.myArgs.repeatLabel, "\t".join([str(x) for x in ind_buffer_or]),
                                                                     "\t".join([str(x) for x in ind_buffer_pv])))
            one_comparison_report_writer.close()

        comparison_reader.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::performComparisonsWithDirection \n"%SimpleTime.now())

    def reportOverlaps(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::computeIncrementalOverlaps\n"%SimpleTime.now())

        self.reportBackground()
        self.reportOddsRatioCumulative(True)
        # self.reportOddsRatioCumulative(False)
        self.performComparisons(True)
        # self.performComparisons(False)

        self.reportOddsRatioWithDirection() # only cumulative report
        self.performComparisonsWithDirection() # only cumulative report

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FeatureRepeatOverlap::computeIncrementalOverlaps\n"%SimpleTime.now())


    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FeatureRepeatOverlap::work\n"%SimpleTime.now())
        self.genomic_window = int(self.myArgs.genomicWindow)
        self.genomic_increment = int(self.myArgs.genomicIncrement)

        if (self.genomic_window%self.genomic_increment!=0):
            sys.stderr.write("Genomic window size %s should be divisible with %s\n"%(self.genomic_window, self.genomic_increment))
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
