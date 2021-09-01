#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
import glob as glob
import subprocess
import math
from operator import itemgetter

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class Select_mQTL_CoRSIV_Quartiles:
    DEBUG_PROGRESS              = True
    DEBUG_LOAD_CORSIV_DEF       = True
    DEBUG_LOAD_R2               = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Given the result of mQTL goodness of fit R2 correlation, select
*  Q1 and Q4
*  Q1+Q2 vs Q3+Q4 eg bottom half and top half
*  and generate bed files for other downstream analyses, including repeat overlaps.

"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)

        parser.add_argument('-m','--mQTLR2CSVFile', help='TAB-delimited text file containing the mQTL R2 results',     required=True)
        parser.add_argument('-b','--corsivBedDef',  help='bed file with corsiv definition',         required=True)
        parser.add_argument('-c','--indexCorsiv',   help='0-based index for corsiv name column',    required=True)
        parser.add_argument('-r','--indexR2',       help='0-based index for R2 column',             required=True)
        parser.add_argument('-o','--outputRoot',    help='output file root',                        required=True)

        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def loadCoRSIVDef(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_mQTL_CoRSIV_Quartiles::loadCoRSIVDef\n"%SimpleTime.now())

        self.corsiv_def_hash={}

        with open(self.myArgs.corsivBedDef, "rt") as corsiv_bed_reader:
            for line in corsiv_bed_reader:
                ff = line.strip().split('\t')
                corsiv_id = ff[3]
                self.corsiv_def_hash[corsiv_id]=line
                if (self.DEBUG_LOAD_CORSIV_DEF):
                    sys.stderr.write("Corsiv %s ==> def ==> %s"%(corsiv_id, line))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_mQTL_CoRSIV_Quartiles::loadCoRSIVDef\n"%SimpleTime.now())

    def loadMQTLR2(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_mQTL_CoRSIV_Quartiles::loadMQTLR2\n"%SimpleTime.now())

        self.corsiv_r2_list = []
        self.corsiv_index = int(self.myArgs.indexCorsiv)
        self.r2_index = int(self.myArgs.indexR2)

        with open(self.myArgs.mQTLR2CSVFile, "rt") as mqtl_reader:
            line_idx = -1
            for line in mqtl_reader:
                line_idx += 1
                if (line_idx==0):
                    continue
                else:
                    ff = line.strip().split('\t')
                    corsiv_id = ff[self.corsiv_index]
                    corsiv_r2 = ff[self.r2_index]

                    if (self.DEBUG_LOAD_R2):
                        sys.stderr.write("Corsiv %s R2 %s \n"%(corsiv_id, corsiv_r2))
                    self.corsiv_r2_list.append( (corsiv_id, corsiv_r2) )

        self.sorted_corsiv_r2_list = sorted(self.corsiv_r2_list, key=itemgetter(1))

        if (self.DEBUG_LOAD_R2):
            sys.stderr.write("CorsivIdx CoRSIV R2\n")
            for corsiv_idx in range(len(self.sorted_corsiv_r2_list)):
                sys.stderr.write("[%s] %s %s\n"%(corsiv_idx, self.sorted_corsiv_r2_list[corsiv_idx][0],
                    self.sorted_corsiv_r2_list[corsiv_idx][1]))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_mQTL_CoRSIV_Quartiles::loadMQTLR2\n"%SimpleTime.now())

    def outputBedSlice(self, start_index, stop_index, output_file):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_mQTL_CoRSIV_Quartiles::outputBedSlice [%s, %s] -> %s\n"%(SimpleTime.now(),
                        start_index, stop_index, output_file))

        with open(output_file, "wt") as bed_slice_writer:
            for idx in range(start_index, stop_index+1):
                corsiv_id = self.sorted_corsiv_r2_list[idx][0]
                corsiv_def_line =self.corsiv_def_hash[corsiv_id]
                bed_slice_writer.write("%s"%corsiv_def_line)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_mQTL_CoRSIV_Quartiles::outputBedSlice [%s, %s] -> %s\n"%(SimpleTime.now(),
                        start_index, stop_index, output_file))

    def outputResults(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_mQTL_CoRSIV_Quartiles::outputResults\n"%SimpleTime.now())

        corsiv_count = len(self.sorted_corsiv_r2_list)

        q1_index = int(corsiv_count/4)-1
        q2_index = int(corsiv_count/2)-1
        q3_index = corsiv_count - q1_index -1

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("count %s q1 %s q2 %s q3 %s\n"%(corsiv_count, q1_index, q2_index, q3_index))

        q1_file = "%s.Q1.bed"%self.myArgs.outputRoot
        q1q2_file = "%s.Q1_Q2.bed"%self.myArgs.outputRoot
        q3q4_file = "%s.Q3_Q4.bed"%self.myArgs.outputRoot
        q4_file = "%s.Q4.bed"%self.myArgs.outputRoot

        self.outputBedSlice(0, q1_index, q1_file)
        self.outputBedSlice(0, q2_index, q1q2_file)
        self.outputBedSlice(q2_index+1, corsiv_count-1, q3q4_file)
        self.outputBedSlice(q3_index, corsiv_count-1, q4_file)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_mQTL_CoRSIV_Quartiles::outputResults\n"%SimpleTime.now())


    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Select_mQTL_CoRSIV_Quartiles::work\n"%SimpleTime.now())

        self.loadCoRSIVDef()
        self.loadMQTLR2()
        self.outputResults()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Select_mQTL_CoRSIV_Quartiles::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = Select_mQTL_CoRSIV_Quartiles.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = Select_mQTL_CoRSIV_Quartiles(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
