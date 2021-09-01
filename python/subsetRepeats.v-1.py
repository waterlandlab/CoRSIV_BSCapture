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

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleStats import SimpleStats

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class SubsetRepeats:
    DEBUG_PROGRESS                 = True
    DEBUG_LOAD_CORSIV_DEF          = True
    DEBUG_LOAD_ALL_CORSIVs         = True
    DEBUG_PROCESS_ONE_MQTL         = True
    DEBUG_BETA_BY_DISTANCE         = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Preprocess UCSC repeats and select those with a count above user specified threshold (minimal 10,000)
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-r','--repeatMaskFile',        help='UCSC repeat mask file (gzipped)',  required=True)
        parser.add_argument('-t','--repeatCountThreshold',  help='repeat count threshold',      required=False, default=10000)
        parser.add_argument('-o','--outputRoot',            help='output file root',            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    def countRepeats(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubsetRepeats::work\n"%SimpleTime.now())

        self.repeat_count_hash = {}
        
        repeat_reader = gzip.open(self.myArgs.repeatMaskFile, "rt")
        
        for line in repeat_reader:
            ff = line.strip().split('\t')
            l1 = ff[11]
            l2 = ff[12]
            l3 = ff[10]
            
            key_1 = l1
            key_2 = "%s_%s"%(l1, l2)
            key_3 = "%s_%s_%s"%(l1, l2, l3)
            
            if not (key_1 in self.repeat_count_hash):
                self.repeat_count_hash[key_1]=0
            self.repeat_count_hash[key_1] += 1
            
            if not (key_2 in self.repeat_count_hash):
                self.repeat_count_hash[key_2]=0
            self.repeat_count_hash[key_2] += 1
            
            if not (key_3 in self.repeat_count_hash):
                self.repeat_count_hash[key_3]=0
            self.repeat_count_hash[key_3] += 1
        
        repeat_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubsetRepeats::work\n"%SimpleTime.now())
   
    def subsetRepeats(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubsetRepeats::work\n"%SimpleTime.now())

        repeat_count_threshold = 10000
        
        repeat_count_report = "%s.repeat_count.xls"%self.myArgs.outputRoot
        passing_repeat_count_report = "%s.passing_repeat_count.xls"%self.myArgs.outputRoot
        repeat_count_report_writer = open(repeat_count_report, "wt") 
        passing_repeat_count_writer = open(passing_repeat_count_report, "wt")
        
        repeat_count_report_writer.write("RepeatKey\tCount\n")
        passing_repeat_count_writer.write("RepeatKey\tCount\n")
        
        fd_hash = {}
        for repeat_key in self.repeat_count_hash:
            repeat_count = self.repeat_count_hash[repeat_key]
        
            repeat_count_report_writer.write("%s\t%s\n"%(repeat_key, repeat_count))     
            
            if (repeat_count>repeat_count_threshold):
                rmsk_file = "%s.%s.bed.gz"%(self.myArgs.outputRoot, repeat_key)
                fd_writer = gzip.open(rmsk_file, "wt")
                fd_hash[repeat_key] = fd_writer
                passing_repeat_count_writer.write("%s\t%s\n"%(repeat_key, repeat_count))     
       
        repeat_count_report_writer.close() 
        passing_repeat_count_writer.close()
        
        repeat_reader = gzip.open(self.myArgs.repeatMaskFile, "rt")
        for line in repeat_reader:
            ff = line.strip().split('\t')
            l1 = ff[11]
            l2 = ff[12]
            l3 = ff[10]
            
            key_1 = l1
            key_2 = "%s_%s"%(l1, l2)
            key_3 = "%s_%s_%s"%(l1, l2, l3)
            
            chrom       = ff[5]
            chrom_start = ff[6]
            chrom_stop  = ff[7]
            buffer = [chrom, chrom_start, chrom_stop]
            
            if (self.repeat_count_hash[key_1] >= repeat_count_threshold):
                fd = fd_hash[key_1]
                fd.write("%s\t%s\t1\t+\n"%("\t".join(buffer), key_1))
            
            if (self.repeat_count_hash[key_2] >= repeat_count_threshold):
                fd = fd_hash[key_2]
                fd.write("%s\t%s\t1\t+\n"%("\t".join(buffer), key_2))
            
            if (self.repeat_count_hash[key_3] >= repeat_count_threshold):
                fd = fd_hash[key_3]
                fd.write("%s\t%s\t1\t+\n"%("\t".join(buffer), key_3))
        
        for repeat_key in fd_hash:
            fd = fd_hash[repeat_key]
            fd.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubsetRepeats::work\n"%SimpleTime.now())
   
       
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SubsetRepeats::work\n"%SimpleTime.now())
        
        # load and count repeats 1 level, 2 level, 3 levels
        self.countRepeats()
        
        # select repeats and write them in individual bed.gz files
        self.subsetRepeats()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SubsetRepeats::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = SubsetRepeats.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = SubsetRepeats(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
