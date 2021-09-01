#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import glob as glob

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
class CombineBestmQTLSNPs:
    DEBUG_PROGRESS              = True
    DEBUG_LOAD_BEST_mQTL        = True
    DEBUG_LOAD_ONE_mQTL         = True
    DEBUG_OUTPUT_COMBINED_MQTL  = True
    DEBUG_OUTPUT_BLOCK  = True

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Combine best mQTL results across a set of CoRSIV/mQTL results
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-m','--mQTLfiles',        help='mQTL file pattern',  required=True)
        parser.add_argument('-b','--snpToHapBlock',    help='file for SNPs to haplotype block',  required=True)
        parser.add_argument('-o','--outputRoot',       help='output file root',   required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    def loadOneTissuemQTL(self, mQTL_csv):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineBestmQTLSNPs::loadOneTissuemQTL\n"%SimpleTime.now())
        
        tissue_mQTL_SNP_CoRSIV_hash = {}
        tissue_mQTL_hapblock_CoRSIV_hash = {}
        
        mQTL_reader = open(mQTL_csv, "rt")
       
        line_index = -1
        for line in mQTL_reader:
            line_index += 1
            if (line_index==0):
                continue
            
            ff = line.strip().split(',')
            corsiv_index = int(ff[0])
            chrom_string = ff[2]
            chrom_number = int(chrom_string.replace("chr",""))
            snp_string = ff[8]
            gg = ff[8].split('_')
            snp_pos = int(gg[1])
            beta = float(ff[len(ff)-3])
            if (beta>0):
                beta_sign = 1
            else:
                beta_sign = -1
                
            if (self.DEBUG_LOAD_ONE_mQTL):
                sys.stderr.write("corsiv index %s chrom [%s] %s snp_string %s snp_pos %s full line %s\n"%(corsiv_index, chrom_string, chrom_number, snp_string, snp_pos, line.strip()))
            
            snp_corsiv_key = "%s\t%s\t%s"%(chrom_number, snp_pos, corsiv_index)
            tissue_mQTL_SNP_CoRSIV_hash[snp_corsiv_key]=beta_sign
            
            # lookup hapblock
            if (snp_string in self.snp_to_hapblock_hash):
                hapblock_info = self.snp_to_hapblock_hash[snp_string]
                
                if (self.DEBUG_LOAD_ONE_mQTL):
                    sys.stderr.write("snp string %s hapblock_info %s\n"%(snp_string, hapblock_info))
                    
                hapblock_start = hapblock_info[2]
                block_corsiv_key = "%s\t%s\t%s"%(chrom_number, hapblock_start, corsiv_index)
                tissue_mQTL_hapblock_CoRSIV_hash[block_corsiv_key]=beta_sign
                
                if (self.DEBUG_LOAD_ONE_mQTL):
                    sys.stderr.write("corsiv index %s chrom [%s] %s snp_string %s snp_pos %s hapblock %s\n"%(corsiv_index, chrom_string, chrom_number, snp_string, snp_pos, hapblock_start))
            
        mQTL_reader.close()
             
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineBestmQTLSNPs::loadOneTissuemQTL\n"%SimpleTime.now())
        
        return tissue_mQTL_SNP_CoRSIV_hash, tissue_mQTL_hapblock_CoRSIV_hash
   
    def loadmQTLFiles(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineBestmQTLSNPs::loadmQTL\n"%SimpleTime.now())
        
        self.mQTL_csv_list = glob.glob(self.myArgs.mQTLfiles)
        self.mQTL_file_basename_list = []
        
        if (self.DEBUG_LOAD_BEST_mQTL):
            sys.stderr.write("mQTL files: %s\n"%";".join(self.mQTL_csv_list))
       
        self.mQTL_SNP_CoRSIV_hash = {}
        for mQTL_csv in self.mQTL_csv_list:
            if (self.DEBUG_LOAD_BEST_mQTL):
                sys.stderr.write("Loading mQTL %s\n"%mQTL_csv)
            tissue_basename = os.path.basename(mQTL_csv)
            self.mQTL_file_basename_list.append(tissue_basename)
            self.mQTL_SNP_CoRSIV_hash[tissue_basename]=self.loadOneTissuemQTL(mQTL_csv) # remember that we get a pair here
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineBestmQTLSNPs::loadmQTL\n"%SimpleTime.now())
        
    def outputCombinedmQTLProfiles(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineBestmQTLSNPs::outputCombinedmQTLProfiles\n"%SimpleTime.now())
        
        # collect all keys and sort them in genomic coordinates order; use corsiv index as last key element
        self.snp_corsiv_key_set = set()
        self.snp_corsiv_key_hash = {}
        
        for tissue in self.mQTL_SNP_CoRSIV_hash:
            if (self.DEBUG_OUTPUT_COMBINED_MQTL):
                sys.stderr.write("collecting keys from tissue %s\n"%tissue)
            
            tissue_hash = self.mQTL_SNP_CoRSIV_hash[tissue]
            for snp_key in tissue_hash:
                gg = snp_key.split('\t')
                chrom_number = int(gg[0])
                snp_pos = int(gg[1])
                corsiv_index = int(gg[2])
                tuple_key = (chrom_number, snp_pos, corsiv_index)
                self.snp_corsiv_key_set.add( tuple_key)
                self.snp_corsiv_key_hash[tuple_key] = snp_key
        
        # then output combined report
        self.snp_corsiv_key_list = list(self.snp_corsiv_key_set)
        self.snp_corsiv_key_list.sort()
        
        output_report = "%s.output_report.genomic_coord_corsiv_index.xls"%self.myArgs.outputRoot
        output_report_writer = open(output_report, "wt")
        
        output_report_writer.write("Feature\t%s\n"%"\t".join(self.mQTL_file_basename_list))
        
        for idx in range(len(self.snp_corsiv_key_list)):
            tuple_key = self.snp_corsiv_key_list[idx]
            snp_key = self.snp_corsiv_key_hash[tuple_key]
            heatmap_key = "chr%s_%s_corsiv_%s"%(tuple_key[0], tuple_key[1], tuple_key[2])
            snp_buffer = [heatmap_key]
            
            for tissue in self.mQTL_file_basename_list:
                tissue_hash = self.mQTL_SNP_CoRSIV_hash[tissue]
                if (snp_key in tissue_hash):
                    snp_buffer.append(tissue_hash[snp_key])
                else:
                    snp_buffer.append(0)
            
#            output_report_writer.write("%s\n"%heatmap_key)
            output_report_writer.write("%s\n"%"\t".join([str(x) for x in snp_buffer]))
            
        output_report_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineBestmQTLSNPs::outputCombinedmQTLProfiles\n"%SimpleTime.now())
            
    def outputCombinedmQTLProfiles2(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineBestmQTLSNPs::outputCombinedmQTLProfiles2\n"%SimpleTime.now())
        
        # collect all keys and sort them in genomic coordinates order; use corsiv index as last key element
        self.snp_corsiv_key_set = set()
        self.snp_corsiv_key_hash = {}
        
        # only output hapblock here (eg skip SNPs outside of haplotype blocks)
        self.block_corsiv_key_set = set()
        self.block_corsiv_key_hash = {}
        
        # collect both snp keys and hapblock keys
        for tissue in self.mQTL_SNP_CoRSIV_hash:
            if (self.DEBUG_OUTPUT_COMBINED_MQTL):
                sys.stderr.write("collecting keys from tissue %s\n"%tissue)
            
            tissue_hash_snp, tissue_hash_block = self.mQTL_SNP_CoRSIV_hash[tissue]
            
            for snp_key in tissue_hash_snp:
                gg = snp_key.split('\t')
                chrom_number = int(gg[0])
                snp_pos = int(gg[1])
                corsiv_index = int(gg[2])
                snp_tuple_key = (chrom_number, snp_pos, corsiv_index)
                self.snp_corsiv_key_set.add( snp_tuple_key)
                self.snp_corsiv_key_hash[snp_tuple_key] = snp_key
                
            if (self.DEBUG_OUTPUT_BLOCK):
                sys.stderr.write("collecting blocks from tissue %s\n"%tissue)
            
            for block_key in tissue_hash_block:
                gg = block_key.split('\t')
                chrom_number = int(gg[0])
                block_pos = int(gg[1])
                corsiv_index = int(gg[2])
                block_tuple_key = (chrom_number, block_pos, corsiv_index)
                if not(block_tuple_key in self.block_corsiv_key_set):
                    self.block_corsiv_key_set.add(block_tuple_key)
                    self.block_corsiv_key_hash[block_tuple_key] = block_key
            
                if (self.DEBUG_OUTPUT_BLOCK):
                    sys.stderr.write("tissue %s block_key %s block_tuple  %s \n"%(tissue, block_key, block_tuple_key))
                
        # then output combined report
                
        self.snp_corsiv_key_list = list(self.snp_corsiv_key_set)
        self.snp_corsiv_key_list.sort()
        
        output_report = "%s.snp_report.genomic_coord_corsiv_index.xls"%self.myArgs.outputRoot
        output_report_writer = open(output_report, "wt")
        output_report_writer.write("SNP\t%s\n"%"\t".join(self.mQTL_file_basename_list))
        
        for idx in range(len(self.snp_corsiv_key_list)):
            tuple_key = self.snp_corsiv_key_list[idx]
            snp_key = self.snp_corsiv_key_hash[tuple_key]
            heatmap_key = "chr%s_%s_corsiv_%s"%(tuple_key[0], tuple_key[1], tuple_key[2])
            snp_buffer = [heatmap_key]
            
            for tissue in self.mQTL_file_basename_list:
                tissue_hash_snp, tissue_hash_block = self.mQTL_SNP_CoRSIV_hash[tissue]
                if (snp_key in tissue_hash_snp):
                    snp_buffer.append(tissue_hash_snp[snp_key])
                else:
                    snp_buffer.append(0)
            
            output_report_writer.write("%s\n"%"\t".join([str(x) for x in snp_buffer]))
            
        output_report_writer.close()

        # block level report
        self.block_corsiv_key_list = list(self.block_corsiv_key_set)
        self.block_corsiv_key_list.sort()
        
        output_report = "%s.block_report.genomic_coord_corsiv_index.xls"%self.myArgs.outputRoot
        output_report_writer = open(output_report, "wt")
        output_report_writer.write("Block\t%s\n"%"\t".join(self.mQTL_file_basename_list))
        
        for idx in range(len(self.block_corsiv_key_list)):
            tuple_key = self.block_corsiv_key_list[idx]
            block_key = self.block_corsiv_key_hash[tuple_key]
            heatmap_key = "chr%s_%s_corsiv_%s"%(tuple_key[0], tuple_key[1], tuple_key[2])
            block_buffer = [heatmap_key]
            
            for tissue in self.mQTL_file_basename_list:
                tissue_hash_snp, tissue_hash_block = self.mQTL_SNP_CoRSIV_hash[tissue]
                if (block_key in tissue_hash_block):
                    block_buffer.append(tissue_hash_block[block_key])
                else:
                    block_buffer.append(0)
            
            output_report_writer.write("%s\n"%"\t".join([str(x) for x in block_buffer]))
            
        output_report_writer.close()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineBestmQTLSNPs::outputCombinedmQTLProfiles2\n"%SimpleTime.now())
            
    def loadHapBlockInfo(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineBestmQTLSNPs::work\n"%SimpleTime.now())
        
        self.snp_to_hapblock_hash = {}
        hap_block_reader = open(self.myArgs.snpToHapBlock, "rt")
        
        for line in hap_block_reader:
            ff = line.strip().split('\t')
            chrom = ff[0]
            snp_id = ff[3]
            gg = snp_id.split('_')
            hap_block_start = int(ff[7])
            hap_block_stop  = int(ff[8])
            self.snp_to_hapblock_hash[snp_id]=(snp_id, chrom, hap_block_start, hap_block_stop)
        
        hap_block_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineBestmQTLSNPs::work\n"%SimpleTime.now())
        
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CombineBestmQTLSNPs::work\n"%SimpleTime.now())
        
        self.loadHapBlockInfo()
        # load and count repeats 1 level, 2 level, 3 levels
        self.loadmQTLFiles()
        
        # select repeats and write them in individual bed.gz files
        # self.outputCombinedmQTLProfiles()
        self.outputCombinedmQTLProfiles2()
        # note: we will have a separate output only for SNPs in haplotype blocks
        
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CombineBestmQTLSNPs::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = CombineBestmQTLSNPs.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = CombineBestmQTLSNPs(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
