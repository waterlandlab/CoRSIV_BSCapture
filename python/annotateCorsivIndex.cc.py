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
from CCLabUtils.simpleStats import SimpleStats

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class AnnotateCoRSIVIndex:
    DEBUG_PROGRESS           = True
    DEBUG_LOAD_CORSIV_DEF    = True
    DEBUG_ANNOTATE_CORSIV    = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Sort and index corsivs by chromosome (1 through  22) and start position
Load and annotate corsiv-based TAB delimited reports
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-b','--corsivBedFile',     help='corsiv bed file',                 required=True)
        parser.add_argument('-r','--corsivReportCSV',   help='comma separated variables corsiv-based report',  required=True)
        parser.add_argument('-c','--corsivColumn',      help='0-based column for the CoRSIV',   required=True)
        parser.add_argument('-o','--outputRoot',        help='output file root',                required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    # load corsiv info chrom as string, as number (1-22), startPos, endPos, name
    # sort by chrom index number and start pos
    # assing a corsiv index to each corsiv
    def loadAndIndexCoRSIVs(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AnnotateCoRSIVIndex::loadAndIndexCoRSIVs %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))

        self.corsiv_def_hash = {}
        self.corsiv_list = []
        
        corsiv_reader = open(self.myArgs.corsivBedFile, "rt")

        for line in corsiv_reader:
            
            ff = line.strip().split('\t')
            corsiv_short = ff[3]
            chrom = ff[0]
            chrom_number = int(ff[0][3:])
            chrom_start = int(ff[1])
            chrom_stop  = int(ff[2])
           
            corsiv_info = CStruct(corsiv_short = corsiv_short, chrom = chrom, chrom_number = chrom_number,
                                  chrom_start = chrom_start, chrom_stop = chrom_stop, corsiv_index = -1)
            self.corsiv_def_hash[corsiv_short]=corsiv_info
            self.corsiv_list.append(corsiv_info)
            
            if (self.DEBUG_LOAD_CORSIV_DEF):
                sys.stderr.write("Loaded corsiv definition  short: %s chrom number %s chrom/start/stop: %s/%s/%s \n"%(corsiv_short, corsiv_info.chrom_number, self.corsiv_def_hash[corsiv_short].chrom, self.corsiv_def_hash[corsiv_short].chrom_start, self.corsiv_def_hash[corsiv_short].chrom_stop))
                
        corsiv_reader.close()

        corsiv_index_file = "%s.corsiv_index.xls"%self.myArgs.outputRoot
        corsiv_index_file_writer = open(corsiv_index_file, "wt")
        corsiv_index_file_writer.write("CoRSIV_Index\tChrom\tChromStart\tChromStop\tCorSIV\n")
        
        # sort all corsivs based on chrom and start
        self.corsiv_list.sort(key=lambda x:(x.chrom_number, x.chrom_start))
        
        for idx in range(len(self.corsiv_list)):
            corsiv_info = self.corsiv_list[idx]
            corsiv_info.corsiv_index = idx
            buffer = [corsiv_info.corsiv_index, corsiv_info.chrom, corsiv_info.chrom_start, corsiv_info.chrom_stop, corsiv_info.corsiv_short]
            corsiv_index_file_writer.write("%s\n"%"\t".join([str(x) for x in buffer]))
        
        corsiv_index_file_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AnnotateCoRSIVIndex::loadAndIndexCoRSIVs %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))
    
    def annotateCoRSIVIndex(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AnnotateCoRSIVIndex::annotateCoRSIVIndex\n"%SimpleTime.now())

        corsiv_report_reader = open(self.myArgs.corsivReportCSV, "rt")
    
        annotated_corsiv = "%s.corsiv_index.csv"%self.myArgs.outputRoot
        annotated_corsiv_writer = open(annotated_corsiv, "wt")
        
        
        line_idx = -1
        corsiv_report_column_idx = int(self.myArgs.corsivColumn)
        for line in corsiv_report_reader:
            ff = line.strip().replace('"','').split(',')
            line_idx += 1
            if (line_idx==0):
                ff[0]="DummyIndex"
                annotated_corsiv_writer.write("CorsivIndex1,CorsivIndex,Chrom,ChromStart,ChromStop,%s\n"%",".join(ff))
            else:
                corsiv_name  = ff[corsiv_report_column_idx]
                if (corsiv_name in self.corsiv_def_hash):
                    corsiv_info = self.corsiv_def_hash[corsiv_name]
                    annotated_corsiv_writer.write("%s,%s,%s,%s,%s,%s\n"%(corsiv_info.corsiv_index, corsiv_info.corsiv_index,
                                                corsiv_info.chrom, corsiv_info.chrom_start, corsiv_info.chrom_stop, ','.join(ff)))
                    if (self.DEBUG_ANNOTATE_CORSIV):
                        sys.stderr.write("Found corsiv %s in hash index %s\n"%(corsiv_name, corsiv_info.corsiv_index))
                                         
        annotated_corsiv_writer.close()
        corsiv_report_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AnnotateCoRSIVIndex::annotateCoRSIVIndex\n"%SimpleTime.now())

    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AnnotateCoRSIVIndex::work\n"%SimpleTime.now())
        
        self.loadAndIndexCoRSIVs()
        self.annotateCoRSIVIndex()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AnnotateCoRSIVIndex::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = AnnotateCoRSIVIndex.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = AnnotateCoRSIVIndex(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
