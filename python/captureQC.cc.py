#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import glob
import datetime
import math
from argparse import RawTextHelpFormatter
import gzip
import random

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleData import SimpleData

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class BamCaptureQC:
    DEBUG_PROGRESS              = True
    DEBUG_CONVERT_BAM_TO_BED    = True
    DEBUG_SETUP                 = True
    DEBUG_OVERLAP               = True
    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.
            
Compute QC metrics for BAM files resulting from DNA capture
TODO: pool of processes
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-b','--bamFilePattern',    help='BAM file(s) pattern',                                 required=True)
        parser.add_argument('-c','--captureTarget',     help='BED file with capture targets',                       required=True)
        parser.add_argument('-r','--radius',            help='radius around capture targets (in basepairs)',        required=True)
        parser.add_argument('-e','--experiment',        help='experiment label',                                    required=True)
        parser.add_argument('-k','--keepTmp',           help='[optional] keep tmp files (default False)',           action="store_true")
        parser.add_argument('-o','--outputRoot',        help='output file root (includes folder name if needed)',   required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    
    def setupAnalysis(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BamCaptureQC::setupAnalysis\n"%SimpleTime.now())
            
        self.tmp_folder = "%s.tmp.%s"%(os.path.abspath(self.myArgs.outputRoot), os.getpid())
        
        if (self.DEBUG_SETUP):
            sys.stderr.write("Setting up tmp folder %s\n"%self.tmp_folder)
        
        os.system("mkdir -p %s"%self.tmp_folder)
        
        # collect bam file list
        self.bam_file_list = glob.glob(self.myArgs.bamFilePattern)
        
        if (self.DEBUG_SETUP):
            sys.stderr.write("bam files: %s\n"%"\n\t".join(self.bam_file_list))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BamCaptureQC::setupAnalysis\n"%SimpleTime.now())
            
    
    def convertBedToBam(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BamCaptureQC::convertBedToBam\n"%SimpleTime.now())
        
        self.full_bed_list = [] 
        self.overlap_bed_list = []
        
        self.mapped_reads_count = 0 
        
        bam_idx = -1
        for bam_file in self.bam_file_list:
            bam_idx += 1
            full_bed_file = "%s/full_bed.%s.idx.%s.bed.gz"%(self.tmp_folder, os.path.basename(bam_file), bam_idx)
            self.full_bed_list.append(full_bed_file)
            
            bam_to_bed_command = "bamToBed -i %s | awk '$5>=10' | pigz -p 2 -c > %s"%(bam_file, full_bed_file)
            
            if (self.DEBUG_CONVERT_BAM_TO_BED):
                sys.stderr.write("[%s] converting BAM file %s to BED file %s\nCommand: %s\n"%(SimpleTime.now(), bam_file, full_bed_file, bam_to_bed_command))
                
            os.system(bam_to_bed_command)
            current_bed_count = SimpleData.getGZFileSize(full_bed_file)
            self.mapped_reads_count += current_bed_count
            
            if (self.DEBUG_CONVERT_BAM_TO_BED):
                sys.stderr.write("[%s] got BED file %s of size %s overall %s reads \n"%(SimpleTime.now(), full_bed_file, current_bed_count, self.mapped_reads_count))
                
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BamCaptureQC::convertBedToBam\n"%SimpleTime.now())
            
    def determineCaptureOverlap(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BamCaptureQC::determineCaptureOverlap\n"%SimpleTime.now())
        
        self.overlap_bed_list = []
        self.overlap_read_count = 0
        
        bed_idx = -1
        for full_bed in self.full_bed_list:
            bed_idx += 1
            overlap_bed_file = "%s.overlap.%s.idx.%s.bed.gz"%(full_bed, os.path.basename(self.myArgs.captureTarget), bed_idx)
            bedtool_window_command = "bedtools window -w %s -u -a %s -b %s | pigz -p 2 -c > %s"%(self.myArgs.radius, full_bed, self.myArgs.captureTarget, overlap_bed_file)
            
            if (self.DEBUG_OVERLAP):
                sys.stderr.write("bedtools window command %s\n"%bedtool_window_command)
                
            os.system(bedtool_window_command)
            
            current_overlap_bed_count = SimpleData.getGZFileSize(overlap_bed_file)
            self.overlap_read_count += current_overlap_bed_count
            
            if (self.DEBUG_CONVERT_BAM_TO_BED):
                sys.stderr.write("[%s] got overlap BED file %s of size %s overall %s reads \n"%(SimpleTime.now(), overlap_bed_file, current_overlap_bed_count, self.overlap_read_count))
                         
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BamCaptureQC::determineCaptureOverlap\n"%SimpleTime.now())
                
    
    def reportQC(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BamCaptureQC::reportQC\n"%SimpleTime.now())
            
        report_file = "%s.report.xls"%self.myArgs.outputRoot
        report_file_writer = open(report_file, "wt")
        report_file_writer.write("Experiment\tMapped_reads\tOverlap_reads\tOverlap_reads_ratio\n")
        
        report_file_writer.write("%s\t%s\t%s\t%s\n"%(self.myArgs.experiment, self.mapped_reads_count, self.overlap_read_count, float(self.overlap_read_count)/float(self.mapped_reads_count)))
        
        report_file_writer.close()
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BamCaptureQC::reportQC\n"%SimpleTime.now())
                
    def cleanup(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BamCaptureQC::reportQC\n"%SimpleTime.now())
        
        if not (self.myArgs.keepTmp):
            cleanup_command = "rm -rf %s"%self.tmp_folder
            sys.stderr.write("Cleanup command %s\n"%cleanup_command)
            os.system(cleanup_command)
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BamCaptureQC::reportQC\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BamCaptureQC::work\n"%SimpleTime.now())
            
        self.setupAnalysis()
        self.convertBedToBam()
        self.determineCaptureOverlap()
        self.reportQC()
        self.cleanup()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BamCaptureQC::work\n"%SimpleTime.now())
            

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = BamCaptureQC.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = BamCaptureQC(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
