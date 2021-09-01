#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import glob as glob
import numpy as np
import scipy.stats

from CCLabUtils.simpleTime import SimpleTime
from CCLabUtils.simpleStats import SimpleStats

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

# Approach
# * for each best mQTL file create a hash corsiv/SNP
# * add the corsiv/SNP as info to a list
# * sort list by by SNP chromosome and chromosome pos
# * then output the heatmap by traversing the corsiv/SNP list in sort order
class SummarizeMQTLStats:
    DEBUG_PROGRESS              = True
    DEBUG_SELECT_CORSIV_MQTL    = True
    DEBUG_LOAD_METHYLATION      = True
    DEBUG_LOAD_SNP              = True


    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

compute R2 for a CORSIV and mQTL
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)

        parser.add_argument('-x','--mQTLCandidateFolder',   help='mQTL candidates folder',  required=True)
        parser.add_argument('-t','--tissue',                help='target tissue',  required=True)
        parser.add_argument('-b','--globalMQTLFile',        help='tissue wide mQTL',  required=True)
        parser.add_argument('-s','--mQTLFileSnpIndex',      help='column id (0-based) in mQTL file for the SNP id',  required=True)
        parser.add_argument('-o','--outputRoot',            help='output file root',   required=True)

        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def selectCoRSIVMQTL(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::selectCoRSIVMQTL\n"%SimpleTime.now())

        self.corsiv_folder_basename = os.path.basename(self.myArgs.mQTLCandidateFolder)
        if (self.DEBUG_SELECT_CORSIV_MQTL):
            sys.stderr.write("mQTL candidate folder %s -> %s\n"%(os.path.abspath(self.myArgs.mQTLCandidateFolder), self.corsiv_folder_basename))

        match_corsiv_id = re.search(r'(chr\d+.me_idx..*)', self.corsiv_folder_basename)
        self.corsiv_id = match_corsiv_id.group(1)

        if (self.DEBUG_SELECT_CORSIV_MQTL):
            sys.stderr.write("mQTL candidate folder  %s corsiv id %s\n"%( self.corsiv_folder_basename, self.corsiv_id))

        # now go through the mQTL summary
        self.snp_column_idx = int(self.myArgs.mQTLFileSnpIndex)
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("mQTL SNV idx %s\n"%self.snp_column_idx)


        mQTL_summary_reader = open(self.myArgs.globalMQTLFile, "rt")
        for line in mQTL_summary_reader:
            if (line.find(self.corsiv_id)>0):
                ff = line.strip().split(',')
                self.snp_id = ff[self.snp_column_idx].replace('"','')
                self.mqtl_fdr = ff[len(ff)-2]
                self.beta_coefficient_methylation = float(ff[len(ff)-3])
                check_corsiv_id = ff[7].replace('"','')
                if (self.DEBUG_SELECT_CORSIV_MQTL):
                    sys.stderr.write("Found mQTL %s check corsiv %s %s beta %s fdr %s line %s \n"%(self.snp_id, check_corsiv_id, self.corsiv_id, self.beta_coefficient_methylation, self.mqtl_fdr, line.strip()))

        mQTL_summary_reader.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::selectCoRSIVMQTL\n"%SimpleTime.now())

    def loadMethylation(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::loadMethylation\n"%SimpleTime.now())

        # find the right methylation file
        file_pattern = "%s%smethylation*%s*"%(self.myArgs.mQTLCandidateFolder, os.path.sep, self.myArgs.tissue)

        if (self.DEBUG_LOAD_METHYLATION):
            sys.stderr.write("file pattern %s\n"%file_pattern)

        methylation_list = glob.glob(file_pattern)

        methylation_file = methylation_list[0]

        if (self.DEBUG_LOAD_METHYLATION):
            sys.stderr.write("methylation info file %s\n"%methylation_file)

        methylation_file_reader = open(methylation_file, "rt")
        line_idx = -1
        for line in methylation_file_reader:
            line_idx += 1
            ff = line.strip().split('\t')
            if (line_idx==0):
                self.methylation_samples = ff[1:]
                if (self.DEBUG_LOAD_METHYLATION):
                    sys.stderr.write("Methylation samples %s\n"%'\t'.join(self.methylation_samples))
            else:
                self.methylation_array = np.asfarray(ff[1:])
                if (self.DEBUG_LOAD_METHYLATION):
                    sys.stderr.write("Methylation array %s\n"%'\t'.join([str(x) for x in self.methylation_array]))

        methylation_file_reader.close()

        if (self.DEBUG_LOAD_METHYLATION):
            sys.stderr.write("methylation file %s\n"%methylation_file)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::loadMethylation\n"%SimpleTime.now())

    def loadGeneticInfo(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::loadGeneticInfo\n"%SimpleTime.now())

        # find the right methylation file
        file_pattern = "%s%ssnv*%s*"%(self.myArgs.mQTLCandidateFolder, os.path.sep, self.myArgs.tissue)

        if (self.DEBUG_LOAD_SNP):
            sys.stderr.write("snp file pattern %s\n"%file_pattern)

        snv_list = glob.glob(file_pattern)

        snv_file = snv_list[0]

        if (self.DEBUG_LOAD_SNP):
            sys.stderr.write("SNV info file %s\n"%snv_file)


        snv_file_reader = open(snv_file, "rt")
        line_idx = -1
        # here we need to find the right snp
        for line in snv_file_reader:
            if (self.DEBUG_LOAD_SNP):
                sys.stderr.write("Snp file line %s\n"%line_idx)

            line_idx += 1
            ff = line.strip().split('\t')
            if (line_idx==0):
                self.snv_samples = ff[1:]
                if (self.DEBUG_LOAD_SNP):
                    sys.stderr.write("SNv samples %s\n"%'\t'.join(self.snv_samples))
                if (self.methylation_samples != self.snv_samples):
                    sys.stderr.write("Eek: methylation and snv samples differ\n")
            else:
                if (self.DEBUG_LOAD_SNP):
                    sys.stderr.write("Snp %s vs snp_id %s \n"%(ff[0], self.snp_id))

                if (ff[0]==self.snp_id):
                    for sample_idx in range(1,len(ff)):
                        if (ff[sample_idx]=="NA"):
                            ff[sample_idx]=0
                    self.snp_array = np.asfarray(ff[1:])
                    if (self.DEBUG_LOAD_SNP):
                        sys.stderr.write("Snp %s %s Raw snp array %s\n"%(ff[0], self.snp_id, '\t'.join([str(x) for x in self.snp_array])))
                    self.adjusted_snp_array = self.beta_coefficient_methylation * self.snp_array
                    if (self.DEBUG_LOAD_SNP):
                        sys.stderr.write("Snp %s %s adjusted snp array %s\n"%(ff[0], self.snp_id, '\t'.join([str(x) for x in self.adjusted_snp_array])))

        snv_file_reader.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::loadGeneticInfo\n"%SimpleTime.now())


    def computeR2(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::computeR2\n"%SimpleTime.now())

        rho_value, rho_pvalue = scipy.stats.pearsonr(self.methylation_array, self.adjusted_snp_array)
        r2 = rho_value * rho_value
        dnam_min = min(self.methylation_array)
        dnam_max = max(self.methylation_array)
        dnam_range = dnam_max - dnam_min

        output_report = "%s.output_%s_%s.txt"%(self.myArgs.outputRoot, self.corsiv_id, self.myArgs.tissue)
        output_report_writer = open(output_report, "wt")
        out_header = ["Corsiv", "mQTL_SNP", "beta", "FDR", "R", "R2", "DNA_Methylation_Range", "DNA_Methylation_Min", "DNA_Methylation_Max"]
        out_values = [self.corsiv_id, self.snp_id, self.beta_coefficient_methylation, self.mqtl_fdr, rho_value, r2, dnam_range, dnam_min, dnam_max]
        output_report_writer.write("%s\n"%"\t".join(out_header))
        output_report_writer.write("%s\n"%"\t".join([str(x) for x in out_values]))
        output_report_writer.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::computeR2\n"%SimpleTime.now())

    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummarizeMQTLStats::work\n"%SimpleTime.now())

        # determine corsiv ID, and also SNP id
        self.selectCoRSIVMQTL()
        self.loadMethylation()
        self.loadGeneticInfo()
        self.computeR2()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummarizeMQTLStats::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = SummarizeMQTLStats.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = SummarizeMQTLStats(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
