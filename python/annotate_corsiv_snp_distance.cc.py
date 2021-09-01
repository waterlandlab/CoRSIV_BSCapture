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
    DEBUG_PROCESS_ONE_MQTL_VERBOSE = True
    DEBUG_BETA_BY_DISTANCE         = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

annotate corsiv-snv distance on a csv file
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)

        parser.add_argument('-m','--corsivSNVFilePattern',      help='CSV file(s) pattern containing CORSiVs and nearby SNVs',    required=True)
        parser.add_argument('-b','--corsivBedFile',             help='corsiv bed file',         required=True)
        parser.add_argument('-c','--corsivIndex',               help='field index [0-based] containing the corsiv',         required=True)
        parser.add_argument('-s','--snpIndex',                  help='field index [0-based] containing the SNV',         required=True)
        parser.add_argument('-o','--outputRoot',                help='output file root',        required=True)

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

    def process_one_corsiv_snp_file(self, one_corsiv_snv_file):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::process_one_corsiv_snp_file %s"%(SimpleTime.now(), one_corsiv_snv_file))

        corsiv_snp_annotated_file = "%s.dist-anno.%s"%(self.myArgs.outputRoot, os.path.basename(one_corsiv_snv_file))
        corsiv_snp_annotated_writer = open(corsiv_snp_annotated_file, "wt")

        corsiv_snp_reader = open(one_corsiv_snv_file, "rt")
        corsiv_index = int(self.myArgs.corsivIndex)
        snp_index = int(self.myArgs.snpIndex)

        line_idx =-1
        for line in corsiv_snp_reader:
            line_idx += 1
            ff = line.strip().split(',')
            if (self.DEBUG_PROCESS_ONE_MQTL_VERBOSE):
                sys.stderr.write("Bug [%s, %s]: %s\n"%(one_corsiv_snv_file, line_idx, line.strip()))

            if (line_idx==0):
                corsiv_snp_annotated_writer.write("%s,Distance\n"%",".join(ff))
            else:
                snv_info = ff[snp_index].replace('"','')
                gg = snv_info.split('_')
                snv_chrom = gg[0]
                snv_chrom_pos = int(gg[1])

                corsiv_short = ff[corsiv_index].replace('"','')
                corsiv_info = self.corsiv_def_hash[corsiv_short]

                if (snv_chrom_pos < corsiv_info.chrom_start):
                    snv_corsiv_distance = snv_chrom_pos - corsiv_info.chrom_start
                elif (snv_chrom_pos > corsiv_info.chrom_stop):
                    snv_corsiv_distance = snv_chrom_pos - corsiv_info.chrom_stop
                else:
                    snv_corsiv_distance = 0

                corsiv_snp_annotated_writer.write("%s,%s\n"%(",".join(ff), snv_corsiv_distance))

        corsiv_snp_reader.close()

        corsiv_snp_annotated_writer.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_mQTL::process_one_corsiv_snp_file %s"%(SimpleTime.now(), one_corsiv_snv_file))

    def process_all_corsiv_snp_files(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::process_all_corsiv_snp_files\n"%SimpleTime.now())

        self.corsiv_snp_file_list = glob.glob(self.myArgs.corsivSNVFilePattern)
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("Found %s corsiv / snv files\n"%len(self.corsiv_snp_file_list))

        for one_corsiv_snv_file in self.corsiv_snp_file_list:
            self.process_one_corsiv_snp_file(one_corsiv_snv_file)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_mQTL::process_all_corsiv_snp_files\n"%SimpleTime.now())


    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_mQTL::work\n"%SimpleTime.now())

        self.loadCoRSIVDef()
        self.process_all_corsiv_snp_files()

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
