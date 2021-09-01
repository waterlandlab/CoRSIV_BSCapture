#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
from scipy.stats import rankdata

import matplotlib as mpl
import pandas as pd
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sn


from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class Run_CoRSIV_mQTL:
    DEBUG_PROGRESS                      = True
    DEBUG_SETUP                         = True
    DEBUG_LOAD_SAMPLE_MAP               = True
    DEBUG_TISSUE_COVERAGE               = True
    DEBUG_CLEANUP_COLUMNS               = True
    DEBUG_CORSIV_SHORT_HASH             = True
    DEBUG_LOAD_SNV_SCHEMA               = True
    DEBUG_PROCESS_ONE_CORSIV            = True
    DEBUG_PROCESS_ONE_CORSIV_VERBOSE    = False

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

run a set of corsivs against SNVs
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
         
        parser.add_argument('-x','--corsivMethylationMatrix',   help='data matrix with beta values for CORSiVs',    required=True)
        parser.add_argument('-b','--corsivBedFile',         help='corsiv bed file',         required=True)
        parser.add_argument('-s','--snvFile',               help='GTEx phased SNV file',    required=True)
        parser.add_argument('-c','--chromosome',            help='single chromosome',       required=True)
        parser.add_argument('-S','--tissueMap',             help='map for samples/tissues/capture files',   required=True)
        parser.add_argument('-H','--haplotypeBlocks',       help='haplotype block plink.blocks.det file',   required=True)
        parser.add_argument('-d','--corsivBlockDistance',   help='maximum haplotype block/corsiv distance (bp)',   required=False, default=1000)
        parser.add_argument('-o','--outputRoot',            help='output file root',        required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args
   
    def loadTissueMap(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::loadTissueMap\n"%SimpleTime.now())

        self.tissue_hash = {}
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

                tissue_info = CStruct(seq_file_id = seq_file_id, gtex_donor_id=gtex_donor_id, gtex_tissue_id=tissue, gender = gender, tissue = tissue)
                self.tissue_hash[seq_file_id] = tissue_info

        sample_map_reader.close()
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::loadTissueMap\n"%SimpleTime.now())
            
            
    def cleanupColumns(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::cleanupColumns\n"%SimpleTime.now())
        
        # clean up tissue names
        column_list = list(self.coverage_df.columns)
        
        cleanup_column_list = []
        for col_idx in range(len(column_list)):
            col_header = column_list[col_idx]
            ff = col_header.split("_")
            seq_file_id = ff[0]
            if (seq_file_id in self.tissue_hash):
                new_column_name = "%s %s %s %s"%(seq_file_id, self.tissue_hash[seq_file_id].gtex_donor_id, self.tissue_hash[seq_file_id].gender, self.tissue_hash[seq_file_id].tissue)
                cleanup_column_list.append(new_column_name)
                if (self.DEBUG_CLEANUP_COLUMNS):
                    sys.stderr.write("Extract Seq ID: [%s] %s --> new column name %s\n"%(col_idx, col_header, new_column_name))
            else:
                sys.stderr.write("Could not find info on %s --> %s\n"%(col_header, seq_file_id))
                sys.exit(4)
                
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::cleanupColumns\n"%SimpleTime.now())
            
        return cleanup_column_list          
            
    def findTissueCoverage(self, column_list):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::findTissueCoverage\n"%SimpleTime.now())
            
        tissue_count_hash = {}
        max_coverage = 0
        tissue_to_donor_hash = {}
        
        for col_idx in range(len(column_list)):
            col_name = column_list[col_idx]
            ff = col_name.split()
            seq_file_id = ff[0]
            
            gtex_donor_id = ff[1]
            if not (gtex_donor_id in tissue_count_hash):
                tissue_count_hash[gtex_donor_id]=set()
                
            # mark the tissue type for this donor 
            gtex_tissue_id = self.tissue_hash[seq_file_id].gtex_tissue_id
            if not (gtex_tissue_id in tissue_to_donor_hash):
                tissue_to_donor_hash[gtex_tissue_id]={}
            tissue_to_donor_hash[gtex_tissue_id][gtex_donor_id]=col_idx
            
            if (self.DEBUG_TISSUE_COVERAGE):
                sys.stderr.write("Added to tissue %s (%s --> %s : %s)\n"%(gtex_tissue_id, gtex_donor_id, col_idx, col_name))
                
            tissue_count_hash[gtex_donor_id].add(col_idx)
            new_set_length = len(tissue_count_hash[gtex_donor_id])
            if (new_set_length>max_coverage):
                max_coverage = new_set_length
                if (self.DEBUG_TISSUE_COVERAGE):
                    sys.stderr.write("bumped max tissue coverage to %s for %s : %s\n"%(max_coverage, gtex_donor_id, seq_file_id))
          
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::findTissueCoverage\n"%SimpleTime.now())
            
        return tissue_count_hash, max_coverage, tissue_to_donor_hash
    
    def loadBetaValues(self):    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::loadBetaValues\n"%SimpleTime.now())
            
        self.coverage_df = pd.read_csv(self.myArgs.corsivMethylationMatrix, sep="\t", index_col=0, header=0)
        
        cleanup_column_list = self.cleanupColumns()
        self.coverage_df.columns = cleanup_column_list
        
        self.multi_tissue_hash, self.max_coverage, self.tissue_to_donor_hash = self.findTissueCoverage(cleanup_column_list)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::loadBetaValues\n"%SimpleTime.now())
            
    def setupCoRSIVToMethylationIndexHash(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::setupCoRSIVToMethylationIndexHash\n"%SimpleTime.now())
        
        self.coverage_index = list(self.coverage_df.index)
        
        self.corsiv_short_to_methylation_index_hash = {}
        
        for corsiv_index in range(len(self.coverage_index)):
            corsiv_long = self.coverage_index[corsiv_index]
            ff = corsiv_long.split(';')
            gg = ff[3].split('.')
            corsiv_short = '.'.join(gg[0:3])
            self.corsiv_short_to_methylation_index_hash[corsiv_short]=corsiv_index
            if (self.DEBUG_CORSIV_SHORT_HASH):
                sys.stderr.write("[%s] Corsiv long %s short %s corsiv index %s\n"%(corsiv_index, corsiv_long, corsiv_short, self.corsiv_short_to_methylation_index_hash[corsiv_short]))
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::setupCoRSIVToMethylationIndexHash\n"%SimpleTime.now())
    
    def loadSNVSchema(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::loadSNVSchema\n"%SimpleTime.now())
        
        self.gtex_snv_to_index_hash = {}
        
        snv_reader = open(self.myArgs.snvFile, "rt")

        line = snv_reader.readline()
        ff = line.strip().split('\t')
    
        for snv_sample_index in range(9,len(ff)):
            sample_name = ff[snv_sample_index]
            self.gtex_snv_to_index_hash[sample_name]=snv_sample_index
            
            if (self.DEBUG_LOAD_SNV_SCHEMA):
                sys.stderr.write("[%s] SNV sample %s\n"%(snv_sample_index, sample_name))
        
        snv_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::loadSNVSchema\n"%SimpleTime.now())
        
    # build a hash of GTEx sample to CoRSIV column
    # build a hash of GTEx sample to SNV column
    # find common samples between corsivs and SNVs
    def setupCorsivsAndSNVs(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::setupCorsivsAndSNVs\n"%SimpleTime.now())
        
        self.loadTissueMap()
        self.loadBetaValues()
        self.setupCoRSIVToMethylationIndexHash()
        self.loadSNVSchema()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::setupCorsivsAndSNVs\n"%SimpleTime.now())
        
    
    # traverse corsiv file and process one by one
    def processAllCoRSIVs(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processAllCoRSIVs\n"%SimpleTime.now())

        corsiv_reader = open(self.myArgs.corsivBedFile, "rt")
        
        for line in corsiv_reader:
            self.processOneCoRSIV(line)
        
        corsiv_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::processAllCoRSIVs\n"%SimpleTime.now())
    
    def filterSNVByHapBlock(self, hap_line, hap_block_name, corsiv_short, out_corsiv_folder):
        
        one_line_hap_bed = "%s%sone_line_%s.%s.bed"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name)
        
        one_line_writer = open(one_line_hap_bed, "wt")
        one_line_writer.write(hap_line)
        one_line_writer.close()
        
        temp_SNV_file = "%s%stmp_SNV.1bp.%s.%s.vcf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name)
        bedtools_command = "bedtools window -u -w 1 -a %s -b %s > %s"%(self.myArgs.snvFile, one_line_hap_bed, temp_SNV_file)
        
        if (self.DEBUG_PROCESS_ONE_CORSIV):
            sys.stderr.write("out corsiv folder %s tmp snv file %s one-line hap bed %s command %s\n"%(out_corsiv_folder, temp_SNV_file, one_line_hap_bed, bedtools_command))
        
        os.system(bedtools_command)
        
        hap_block_snvs_set = self.hapblock_to_snv_hash[hap_block_name] 
        # finally filter by SNVs in the block
        temp_SNV_file_block = "%s%stmp_SNV.1bp.%s.%s.hap_snvs.vcf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name)
        temp_SNV_file_block_writer = open(temp_SNV_file_block, "wt")
        temp_SNV_file_reader = open(temp_SNV_file, "rt")
        for line in temp_SNV_file_reader:
            ff = line.strip().split('\t')
            snv_id = ff[2]
            if (snv_id in hap_block_snvs_set):
                temp_SNV_file_block_writer.write(line)
                
        temp_SNV_file_reader.close()
        
        temp_SNV_file_block_writer.close()
        return temp_SNV_file_block
    
    # for each tissue
    # find common samples between tissue/capture and SNV that also have methylation
    # setup SNP and DATA for eMatrixQTL
    # invoke script
    def processOneCoRSIVOneTissue(self, gtex_tissue, filtered_SNV_file, out_corsiv_folder, corsiv_short, hap_block_name):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processOneCoRSIVOneTissue corsiv %s tissue %s"%(SimpleTime.now(), gtex_tissue, corsiv_short))
    
        corsiv_index = self.corsiv_short_to_methylation_index_hash[corsiv_short]
        
        donor_hash = self.tissue_to_donor_hash[gtex_tissue]
        
        if (self.DEBUG_PROCESS_ONE_CORSIV):
                sys.stderr.write("Process tissue %s with %s donors\n"%(gtex_tissue, len(donor_hash)))
        
        donors_with_methylation_and_SNV_list = []
        methylation_donors_col_idx = []
        snv_donor_col_idx = []
        
        for gtex_donor in donor_hash:
            if (self.DEBUG_PROCESS_ONE_CORSIV_VERBOSE):
                sys.stderr.write("Checking donor %s"%gtex_donor)
            if gtex_donor in self.gtex_snv_to_index_hash:
                gtex_donor_methylation_col_idx = donor_hash[gtex_donor]
                gtex_donor_snv_col_idx = self.gtex_snv_to_index_hash[gtex_donor]
                if (self.DEBUG_PROCESS_ONE_CORSIV):
                    sys.stderr.write("Yay! capture donor %s --> (capture  %s ,SNV %s) \n"%(gtex_donor, donor_hash[gtex_donor], gtex_donor_snv_col_idx))
                
                # check if methylation present
                meth_is_na_flag = np.isnan(self.coverage_df.iat[corsiv_index, gtex_donor_methylation_col_idx])
                if not meth_is_na_flag:
                    if (self.DEBUG_PROCESS_ONE_CORSIV_VERBOSE):
                        sys.stderr.write("TwoYays! capture donor %s --> (capture  %s ,SNV %s) methylation %s \n"%(gtex_donor, donor_hash[gtex_donor], self.gtex_snv_to_index_hash[gtex_donor], self.coverage_df.iat[corsiv_index, gtex_donor_methylation_col_idx] ))
        
                    donors_with_methylation_and_SNV_list.append(gtex_donor)                
                    methylation_donors_col_idx.append(gtex_donor_methylation_col_idx)
                    snv_donor_col_idx.append(gtex_donor_snv_col_idx)
            else:
                if (self.DEBUG_PROCESS_ONE_CORSIV_VERBOSE):
                    sys.stderr.write("Sigh! capture donor %s not in SNV\n"%gtex_donor)
        
        if (self.DEBUG_PROCESS_ONE_CORSIV):
            sys.stderr.write("For corsiv short %s and tissue %s found %s %s donors with methSNV\n"%(corsiv_short, gtex_tissue, len(methylation_donors_col_idx), len(snv_donor_col_idx)))
            
        if len(methylation_donors_col_idx)>=20:
            if (self.DEBUG_PROCESS_ONE_CORSIV):
                sys.stderr.write("methSNV: For corsiv short %s and tissue %s found %s %s donors with methSNV: eMatrixQTL\n"%(corsiv_short, gtex_tissue, len(methylation_donors_col_idx), len(snv_donor_col_idx)))
                
            support_haplotype_methylation_data_file = "%s%smethylation_hap_corsiv.%s.hap_block.%s.tissue.%s.txt"%(out_corsiv_folder, os.path.sep,
                                                                                                                 corsiv_short, hap_block_name, gtex_tissue)
            support_haplotype_methylation_data_writer = open(support_haplotype_methylation_data_file, "wt")
            support_haplotype_methylation_data_writer.write("CoRSIV\t%s\n"%"\t".join(donors_with_methylation_and_SNV_list))
            
            
            methylation_tissue_data = self.coverage_df.iloc[corsiv_index, methylation_donors_col_idx]
            support_haplotype_methylation_data_writer.write("%s_beta\t%s\n"%(corsiv_short, "\t".join([str(x) for x in methylation_tissue_data])))
            rank_methylation_tissue_data = rankdata(methylation_tissue_data)
            support_haplotype_methylation_data_writer.write("%s\t%s\n"%(corsiv_short, "\t".join([str(x) for x in rank_methylation_tissue_data])))
                
        
            # setup haplotype SNV data

            filtered_SNV_reader = open(filtered_SNV_file, "rt")
            
            haplotype_sum = np.zeros(len(snv_donor_col_idx))
            for snv_line in filtered_SNV_reader:
                snv_ff = snv_line.strip().split('\t')
                snv_id = snv_ff[2]
                buffer = [snv_id]
                for sample_idx_idx in range(len(snv_donor_col_idx)):
                    snv_sample_idx = snv_donor_col_idx[sample_idx_idx]
                    snv_string = snv_ff[snv_sample_idx]
                    if (snv_string=="0|0"):
                        haplotype_sum[sample_idx_idx] += 0 
                    elif (snv_string=="0|1") or (snv_string=="1|0"):
                        haplotype_sum[sample_idx_idx] += 1 
                    elif (snv_string=="1|1"):
                        haplotype_sum[sample_idx_idx] += 2 
                    else:
                        sys.stderr.write("Eek: snv %s snv_string %s sample_idx %s sample %s \n"%(snv_id, snv_string, snv_sample_idx, donors_with_methylation_and_SNV_list[sample_idx_idx]))
                        buffer.append("NA")
                        
            filtered_SNV_reader.close()
            
            support_haplotype_methylation_data_writer.write("%s_hap\t%s\n"%(corsiv_short, "\t".join([str(x) for x in haplotype_sum])))
            
            
            list_of_lists = []
            for idx in range(len(snv_donor_col_idx)):
                hap_sum = int(haplotype_sum[idx])
                meth_ratio = methylation_tissue_data[idx]
                list_of_lists.append([hap_sum, meth_ratio])
            
            columns_hap_meth = ["Haplotype_group", "Methylation_ratio"]
            hap_meth_data_frame = pd.DataFrame(list_of_lists,  columns = columns_hap_meth)
            
            plt.figure()
            sn.set(style="white")
            pdf_boxplot = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.boxplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
            ax = sn.boxplot(x="Haplotype_group", y="Methylation_ratio", data=hap_meth_data_frame, showfliers = False)
            plt.savefig(pdf_boxplot, bbox_inches='tight', dpi=150)

            plt.figure()
            sn.set(style="white")
            pdf_scatterplot =  "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.scatterplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
            ax = sn.scatterplot(x="Haplotype_group", y="Methylation_ratio", data=hap_meth_data_frame)
            plt.savefig(pdf_scatterplot, bbox_inches='tight', dpi=150)
            
        else:
            if (self.DEBUG_PROCESS_ONE_CORSIV):
                sys.stderr.write("methSNVFail: For corsiv short %s and tissue %s found only %s %s donors with methSNV\n"%(corsiv_short, gtex_tissue, len(methylation_donors_col_idx), len(snv_donor_col_idx)))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::processOneCoRSIVOneTissue corsiv %s tissue %s"%(SimpleTime.now(), gtex_tissue, corsiv_short))
            
    
    # Filter SNV file by 1MB from CoRSIV
    # Find samples with capture
    # Filter by those in capture
    # For each tissue
    # * Setup SNP and DATA inputs for eMatrixQTL
    # * minimum 20 donors per tissue w/ SNV
    # * Run eMatrixQTL
    def processOneCoRSIV(self, corsiv_line):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processOneCoRSIV %s"%(SimpleTime.now(), corsiv_line))
        
        ff = corsiv_line.strip().split('\t')
        chrom = ff[0]
        corsiv_short = ff[3]
        
        if (self.DEBUG_PROCESS_ONE_CORSIV):
            sys.stderr.write("Processing CoRSIV %s: chrom %s corsiv_short %s\n"%(corsiv_line.strip(), chrom, corsiv_short))
            
        if (chrom!=self.myArgs.chromosome):
            sys.stderr.write("Skipping line %s on chromosome %s different from user specified %s\n"%(corsiv_line.strip(), chrom, self.myArgs.chromosome))
            return
        
        out_corsiv_folder = "%s.out_corsiv_folder.%s"%(self.myArgs.outputRoot, corsiv_short)
        os.system("mkdir -p %s"%out_corsiv_folder)
        
        one_line_corsiv_bed = "%s%sone_line_%s.bed"%(out_corsiv_folder, os.path.sep, corsiv_short)
        one_line_writer = open(one_line_corsiv_bed, "wt")
        one_line_writer.write(corsiv_line)
        one_line_writer.close()
        
        # run bedtools w/ all hapblocks def
        hapblock_corsiv_overlap_bed = "%s%shap_corsiv_overlap.%s.bed"%(out_corsiv_folder, os.path.sep, corsiv_short)
        
        hap_block_intersect = "bedtools window -u -w %s -a %s -b %s > %s"%(int(self.myArgs.corsivBlockDistance)+1,
                                          self.haplotype_blocks_bed, one_line_corsiv_bed, hapblock_corsiv_overlap_bed)
        if (self.DEBUG_PROCESS_ONE_CORSIV):
            sys.stderr.write("hapblock/corsiv overlap command %s\n"%hap_block_intersect)
        
        os.system(hap_block_intersect)
       
        hap_corsiv_overlap_reader =  open(hapblock_corsiv_overlap_bed)
        for hap_line in hap_corsiv_overlap_reader:
            ff = hap_line.strip().split("\t")
            hap_block_name = ff[3]
            if (self.DEBUG_PROCESS_ONE_CORSIV):
                sys.stderr.write("Hap block line %s name %s\n"%(hap_line.strip(), hap_block_name))
                
            filtered_SNV_file =  self.filterSNVByHapBlock(hap_line, hap_block_name, corsiv_short, out_corsiv_folder)
         
            for gtex_tissue in self.tissue_to_donor_hash:
                self.processOneCoRSIVOneTissue(gtex_tissue, filtered_SNV_file, out_corsiv_folder, corsiv_short, hap_block_name)
    
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::processOneCoRSIV %s"%(SimpleTime.now(), corsiv_line))
    
    def setupHaplotypeBlocks(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::setupHaplotypeBlocks\n"%SimpleTime.now())

        self.hapblock_to_snv_hash = {}        

        self.haplotype_blocks_bed = "%s.hap_blocks.bed"%self.myArgs.outputRoot
        self.haplotype_blocks_hash = {}
        hap_bed_writer = open(self.haplotype_blocks_bed, "wt")
        
        hap_reader = open(self.myArgs.haplotypeBlocks, "rt")
        line_idx = -1
        for line in hap_reader:
            line_idx += 1
            if (line_idx ==0):
                continue
           
            ff = line.strip().split()
            hap_block_name = "_".join(ff[0:3])
            hap_block_nSNVs = ff[4]
            buffer = ["chr%s"%ff[0], ff[1], ff[2], hap_block_name, str(hap_block_nSNVs), "+"]
            hapblock_line = "\t".join(buffer)
            hap_bed_writer.write("%s\n"%hapblock_line)
            self.haplotype_blocks_hash[hap_block_name]=hapblock_line
            hap_block_SNVs = set(ff[5].split('|'))
            self.hapblock_to_snv_hash[hap_block_name]=hap_block_SNVs
            
        hap_reader.close()
        
        hap_bed_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::setupHaplotypeBlocks\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::work\n"%SimpleTime.now())
        
        self.setupCorsivsAndSNVs()
        self.setupHaplotypeBlocks()
        
        self.processAllCoRSIVs()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = Run_CoRSIV_mQTL.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = Run_CoRSIV_mQTL(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
