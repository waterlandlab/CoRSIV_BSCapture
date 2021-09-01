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

class ComputeMedianMQTLCandidate:
    DEBUG_PROGRESS       = True
    DEBUG_PARSE_MQTL     = False 

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Determine median mQTL candidate distance
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)

        parser.add_argument('-m','--mQTLCandidates',  help='eMatrixQTL result',             required=True)
        parser.add_argument('-p','--pValueLimit',     help='pValue limit',                  required=True)
        parser.add_argument('-o','--outputRoot',      help='output file root',              required=True)

        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def compileMQTLCAndidates(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START ComputeMedianMQTLCandidate::compileMQTLCAndidates\n"%SimpleTime.now())

        self.passing_distance_list = []
        self.p_value_limit =float(self.myArgs.pValueLimit)

        with open(self.myArgs.mQTLCandidates, "rt") as mqtl_reader:
            line_idx = -1
            for line in mqtl_reader:
                line_idx +=1
                if (line_idx==0):
                    continue
                ff = line.strip().split("\t")
                p_value = float(ff[5])
                distance = float(ff[2])
                if (self.DEBUG_PARSE_MQTL):
                    sys.stderr.write("distance %s pvalue %s\n"%(distance, p_value))
                if (p_value<self.p_value_limit):
                    self.passing_distance_list.append(distance)
                    if (self.DEBUG_PARSE_MQTL):
                        sys.stderr.write("pass distance %s pvalue %s\n"%(distance, p_value))

        output_report = "%s.median_mqtl.xls"%self.myArgs.outputRoot
        median_distance = np.median(self.passing_distance_list)

        with open(output_report, "wt") as report_writer:
            report_writer.write("%s\t%s\n"%(os.path.basename(self.myArgs.mQTLCandidates), median_distance))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP ComputeMedianMQTLCandidate::compileMQTLCAndidates\n"%SimpleTime.now())

    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START ComputeMedianMQTLCandidate::work\n"%SimpleTime.now())

        self.compileMQTLCAndidates()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP ComputeMedianMQTLCandidate::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = ComputeMedianMQTLCandidate.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = ComputeMedianMQTLCandidate(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
