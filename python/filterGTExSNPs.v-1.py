#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
from scipy.stats import rankdata

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class FilterGTExSNPsBySampleList:
    DEBUG_PROGRESS                      = True
    DEBUG_SETUP                         = True
    DEBUG_LOAD_SAMPLE_MAP               = True
    DEBUG_LOAD_SNV_SCHEMA               = True
    DEBUG_ANNOTATE_CORSIVs              = True
    DEBUG_CORSIV_X_SAMPLE_SNV           = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

determine SNPs that overlap with a set of bed regions specified by the user and a set of samples
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)

        parser.add_argument('-b','--corsivBedFile',         help='corsiv bed file',         required=True)
        parser.add_argument('-s','--snvFile',               help='GTEx phased SNV file',    required=True)
        parser.add_argument('-c','--chromosome',            help='single chromosome',       required=True)
        parser.add_argument('-t','--tissueMap',             help='file containing the capture sample description',   required=True)
        parser.add_argument('-H','--vcfHeader',             help='file containing the samples',   required=True)
        parser.add_argument('-o','--outputRoot',            help='output file root',        required=True)
        parser.add_argument('-k','--keepTempFiles',         help='[optional] keep temorary files',  action="store_true")

        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def loadGTExDonors(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::loadGTExDonors \n"%SimpleTime.now())

        self.gtex_donor_set = set()
        sample_map_reader = open(self.myArgs.tissueMap)

        line_idx = -1
        for line in sample_map_reader:
            line_idx +=1
            if (line_idx==0):
                continue
            else:
                ff = line.strip().split("\t")
                seq_file_id = ff[0]
                gtex_donor_id = ff[1]
                gtex_tissue_id = ff[2]
                gender = ff[3]
                tissue = "_".join(ff[4].split())

                if (self.DEBUG_LOAD_SAMPLE_MAP):
                    sys.stderr.write("Loading %s seq id --> (%s, %s, %s, %s)\n"%(seq_file_id, gtex_donor_id, gtex_tissue_id, gender, tissue))

                self.gtex_donor_set.add(gtex_donor_id)

        sample_map_reader.close()

        if (self.DEBUG_LOAD_SAMPLE_MAP):
            sys.stderr.write("Loaded %s distinct GTEx donors\n"%len(self.gtex_donor_set))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::loadGTExDonors \n"%SimpleTime.now())


    def loadSNVSchema(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::loadSNVSchema\n"%SimpleTime.now())

        self.gtex_snv_to_index_hash = {}

        with open(self.myArgs.vcfHeader, "rt") as snv_reader:
            line = snv_reader.readline()

        self.header_ff = line.strip().split('\t')

        for snv_sample_index in range(9,len(self.header_ff)):
            sample_name = self.header_ff[snv_sample_index]
            self.gtex_snv_to_index_hash[sample_name]=snv_sample_index

            if (self.DEBUG_LOAD_SNV_SCHEMA):
                sys.stderr.write("[%s] SNV sample %s\n"%(snv_sample_index, sample_name))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::loadSNVSchema\n"%SimpleTime.now())

    # build a hash of GTEx sample to CoRSIV column
    # build a hash of GTEx sample to SNV column
    # find common samples between corsivs and SNVs
    def setupCorsivsAndSNVs(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::setupCorsivsAndSNVs\n"%SimpleTime.now())

        self.loadGTExDonors()
        self.loadSNVSchema()

        self.gtex_genotype_donors = set(self.gtex_snv_to_index_hash.keys())
        self.gtex_donors_capture_genotype = self.gtex_genotype_donors.intersection(self.gtex_donor_set)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("Found %s GTEx donors with capture and genotypes\n"%len(self.gtex_donors_capture_genotype))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::setupCorsivsAndSNVs\n"%SimpleTime.now())

    # convert VCF to BED
    def convertVCFToBedFile(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::convertVCFToBedFile \n"%SimpleTime.now())

        self.snv_bed_file = "%s.vcf_as.bed"%self.myArgs.outputRoot
        snv_bed_file_writer = open(self.snv_bed_file, "wt")

        with open(self.myArgs.snvFile, "rt") as vcf_file_reader:
            for line in vcf_file_reader:
                if (line.find("#")==0): # skip VCF comment lines
                    continue
                ff = line.strip().split("\t")
                chrom = ff[0]
                if (chrom != self.myArgs.chromosome):
                    continue

                chrom_start = int(ff[1])
                chrom_stop = chrom_start + 1
                snv_name = ff[2]
                snv_bed_buffer = [chrom, chrom_start, chrom_stop, snv_name, "1", "+"]
                snv_bed_file_writer.write("%s\n"%"\t".join([str(x) for x in snv_bed_buffer]))

        snv_bed_file_writer.close()

        self.cleanup_list.append(self.snv_bed_file)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::convertVCFToBedFile \n"%SimpleTime.now())

    # annotate corsivs w/ SNPs (hash or corsivs with set of SNVs)
    # encode corsiv as chrom-start-stop, so the same code works for other genomic regions (all corsivs, naive controls, random controls, 450k controls)
    def annotateCoRSIVsWithVCF(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::annotateCoRSIVsWithVCF \n"%SimpleTime.now())


        self.corsiv_to_snv_hash = {}
        # setup only corsivs on current chromosome
        with open(self.myArgs.corsivBedFile, "rt") as corsiv_def_reader:
            for line in corsiv_def_reader:
                ff = line.strip().split('\t')
                chrom = ff[0]
                if (chrom != self.myArgs.chromosome):
                    continue
                genomic_key = "_".join(ff[0:4])
                genomic_name = ff[3]
                genomic_region_info = CStruct(genomic_key=genomic_key, genomic_name=genomic_name, snv_set=set(), sample_x_snv={})
                for gtex_donor in self.gtex_donors_capture_genotype:
                    genomic_region_info.sample_x_snv[gtex_donor]=0

                self.corsiv_to_snv_hash[genomic_key]=genomic_region_info

        if (self.DEBUG_ANNOTATE_CORSIVs):
            sys.stderr.write("setup info for %s regions on chrom %s\n"%(len(self.corsiv_to_snv_hash), self.myArgs.chromosome))

        # run bedtools window command
        tmp_connect_snv_corsiv_file = "%s.temp_connect_snv_corsiv.txt"%self.myArgs.outputRoot
        self.cleanup_list.append(tmp_connect_snv_corsiv_file)

        bedtools_window = "bedtools window -w 0 -a %s -b %s > %s"%(self.snv_bed_file, self.myArgs.corsivBedFile, tmp_connect_snv_corsiv_file)
        if (self.DEBUG_ANNOTATE_CORSIVs):
            sys.stderr.write("bedtools window command %s\n"%bedtools_window)
        os.system(bedtools_window)

        self.snv_to_corsiv_hash = {}

        with open(tmp_connect_snv_corsiv_file, "rt") as tmp_connect_snv_corsiv_file_reader:
            for line in tmp_connect_snv_corsiv_file_reader:
                ff = line.strip().split('\t')
                snv_name = ff[3]
                genomic_key = "_".join(ff[6:10])
                if not genomic_key in self.corsiv_to_snv_hash:
                    sys.stderr.write("Eek %s : %s !!\n"%(line.strip(), genomic_key))
                else:
                    genomic_region_info = self.corsiv_to_snv_hash[genomic_key]
                    genomic_region_info.snv_set.add(snv_name)
                    # there should not be conflicts of SNVs with multiple regions
                    if (snv_name in self.snv_to_corsiv_hash):
                        sys.stderr.write("Gee whiz overlapping genomic regions %s %s\n"%(genomic_key, self.snv_to_corsiv_hash[snv_name]))
                    self.snv_to_corsiv_hash[snv_name]=genomic_key

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::annotateCoRSIVsWithVCF \n"%SimpleTime.now())

    # traverse vcf file
    # for each SNP associated with a corsiv
    # for each GTEx donor w/ capture & SNV
    #    count the corsiv x sample SNV in a hash
    def determineCoRSIV_x_SNP_Sample(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::determineCoRSIV_x_SNP_Sample \n"%SimpleTime.now())

        self.corsiv_x_samples_snv_hash = {}

        with open(self.myArgs.snvFile, "rt") as snv_reader:
            line_idx = -1
            for line in snv_reader:
                line_idx += 1
                ff = line.strip().split('\t')
                snv_name = ff[2]
                if not (snv_name in self.snv_to_corsiv_hash):
                    if (self.DEBUG_CORSIV_X_SAMPLE_SNV):
                        sys.stderr.write("[%s] Skip SNV snv_name %s\n"%(line_idx, snv_name))
                    continue
                else:
                    genomic_key = self.snv_to_corsiv_hash[snv_name]
                    genomic_region_info = self.corsiv_to_snv_hash[genomic_key]
                    # traverse all donors w/ capture & SNVs
                    for gtex_donor in self.gtex_donors_capture_genotype:
                        ff_idx = self.gtex_snv_to_index_hash[gtex_donor]
                        snv_string = ff[ff_idx]
                        if (snv_string=="0|0"):
                            pass
                        else:
                            genomic_region_info.sample_x_snv[gtex_donor]+=1

                        if (self.DEBUG_CORSIV_X_SAMPLE_SNV):
                            sys.stderr.write("[%s, %s] corsiv %s donor [%s] %s %s -> snv_string %s snv_count %s\n"%(line_idx, snv_name,
                                        genomic_key,
                                        ff_idx, gtex_donor, self.header_ff[ff_idx], snv_string,
                                        genomic_region_info.sample_x_snv[gtex_donor]))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::determineCoRSIV_x_SNP_Sample\n"%SimpleTime.now())


    # output the matrix of corsiv x samples
    def outputCoRSIV_x_SNP_Sample(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::outputCoRSIV_x_SNP_Sample \n"%SimpleTime.now())

        report = "%s.region_x_snp_reports.xls"%self.myArgs.outputRoot
        with open(report, "wt") as report_writer:
            sorted_donors = sorted(list(self.gtex_donors_capture_genotype))
            report_writer.write("Genomic_Region\t%s\n"%"\t".join(sorted_donors))

            for genomic_key in self.corsiv_to_snv_hash:
                genomic_info = self.corsiv_to_snv_hash[genomic_key]
                snv_buffer = []
                for gtex_donor in sorted_donors:
                    snv_buffer.append(genomic_info.sample_x_snv[gtex_donor])
                report_writer.write("%s\t%s\n"%(genomic_key, "\t".join([str(x) for x in snv_buffer])))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::outputCoRSIV_x_SNP_Sample\n"%SimpleTime.now())

    def computeSampleLevelCoRSIVSNVs(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::computeSampleLevelCoRSIVSNVs \n"%SimpleTime.now())

        # convert VCF to BED
        self.convertVCFToBedFile()
        # annotate corsivs w/ SNPs (hash or corsivs with set of SNVs)
        # encode corsiv as chrom-start-stop, so the same code works for other genomic regions (all corsivs, naive controls, random controls, 450k controls)
        self.annotateCoRSIVsWithVCF()
        # traverse vcf file
            # for each SNP associated with a corsiv
            # for each GTEx donor w/ capture & SNV
                # count the corsiv x sample SNV in a hash
        self.determineCoRSIV_x_SNP_Sample()
        # output the matrix of corsiv x samples
        self.outputCoRSIV_x_SNP_Sample()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::computeSampleLevelCoRSIVSNVs\n"%SimpleTime.now())

    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START FilterGTExSNPsBySampleList::work\n"%SimpleTime.now())

        self.cleanup_list=[]

        self.setupCorsivsAndSNVs()
        self.computeSampleLevelCoRSIVSNVs()

        # cleanup temporary files
        if not self.myArgs.keepTempFiles:
            for file_name in self.cleanup_list:
                os.unlink(file_name)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP FilterGTExSNPsBySampleList::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = FilterGTExSNPsBySampleList.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = FilterGTExSNPsBySampleList(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
