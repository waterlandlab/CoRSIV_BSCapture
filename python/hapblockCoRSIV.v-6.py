#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np

import scipy.stats
from scipy.stats import rankdata
from scipy.stats import pearsonr
from scipy.stats import spearmanr

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
    DEBUG_PROCESS_ONE_CORSIV_VERBOSE    = True
    DEBUG_PROCESS_HAPLOTYPE_STRINGS     = True
    DEBUG_PROCESS_HAPLOTYPE_STRINGS_GLM = True
    DEBUG_PROCESS_HAPLOTYPE_INFO        = True

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
        parser.add_argument('-d','--corsivBlockDistance',   help='maximum haplotype block/corsiv distance (bp)',   required=False, default=0)
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
    
    
    def processHaplotypeAlleleSum(self, out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue,
                                           haplotype_sum, methylation_tissue_data):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processHaplotypeAlleleSum corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(), gtex_tissue, corsiv_short, hap_block_name))
   
        # compute Pearson
        # compute Spearman
        pcc_value, pcc_pvalue = scipy.stats.pearsonr(haplotype_sum, methylation_tissue_data)
        rho_value, rho_pvalue = scipy.stats.spearmanr(haplotype_sum, methylation_tissue_data)
        
        correlation_report = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.correlations"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        correlation_report_writer = open(correlation_report, "wt")
        
        correlation_report_writer.write("CoRSIV\tTissue\tDonors\tHap_block\tPCC\tPCC_pvalue\tRho\tRho_pvalue\n")
        buffer = [corsiv_short, gtex_tissue, len(haplotype_sum), hap_block_name, pcc_value, pcc_pvalue, rho_value, rho_pvalue]
        correlation_report_writer.write("%s\n"%"\t".join([str(x) for x in buffer]))
        correlation_report_writer.close()
        
        
        list_of_lists = []
        for idx in range(len(haplotype_sum)):
            hap_sum = int(haplotype_sum[idx])
            meth_ratio = methylation_tissue_data[idx]
            list_of_lists.append([hap_sum, meth_ratio])
        
        columns_hap_meth = ["Haplotype allele sum", "CoRSIV Methylation"]
        hap_meth_data_frame = pd.DataFrame(list_of_lists,  columns = columns_hap_meth)
        
        plt.figure()
        sn.set(style="white")
        pdf_boxplot = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.boxplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.boxplot(x="Haplotype allele sum", y="CoRSIV Methylation", data=hap_meth_data_frame, showfliers = False)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        plt.savefig(pdf_boxplot, bbox_inches='tight', dpi=150)

        plt.figure()
        sn.set(style="white")
        pdf_scatterplot =  "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.scatterplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.scatterplot(x="Haplotype allele sum", y="CoRSIV Methylation", data=hap_meth_data_frame)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        plt.savefig(pdf_scatterplot, bbox_inches='tight', dpi=150)
        
        plt.figure()
        sn.set(style="white")
        pdf_swarmplot = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.swarmplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.swarmplot(x="Haplotype allele sum", y="CoRSIV Methylation", data=hap_meth_data_frame)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        plt.savefig(pdf_swarmplot, bbox_inches='tight', dpi=150)

        plt.figure()
        sn.set(style="white")
        pdf_violinplot = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.violinplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.violinplot(x="Haplotype allele sum", y="CoRSIV Methylation", data=hap_meth_data_frame, cut=0)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        plt.savefig(pdf_violinplot, bbox_inches='tight', dpi=150)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::processHaplotypeAlleleSum corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(), gtex_tissue, corsiv_short, hap_block_name))


    def processHaplotypeStringByGLM(self, out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue,
                                         hap_string_to_hap_sum, hap_string_hash):
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processHaplotypeStringByGLM corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(), corsiv_short, gtex_tissue, hap_block_name))
        
        # set up an array of allele sum/lexicographic for hap blocks w/ at least 2 values
        # output them in a txt file
        # run the glm driver
        # load the results and output a glm file output that we can then integrate over all the corsivs
        # plot -log10(p-value) vs glm R2
        # plot hap allele sum R and R2 vs glm R2
        
      
        threshold = 2
        hap_block_list = []
        for hap_string in hap_string_hash:
            hap_methylation_list = hap_string_hash[hap_string]
            hap_sum = hap_string_to_hap_sum[hap_string]
            
            if (len(hap_methylation_list)>=threshold):
                if (self.DEBUG_PROCESS_HAPLOTYPE_STRINGS_GLM):
                    sys.stderr.write("Hap string %s sample count %s passes GLM threshold %s\n"%(hap_string, len(hap_methylation_list), threshold))
                hap_block_list.append( (hap_sum, hap_string) )
            else:
                if (self.DEBUG_PROCESS_HAPLOTYPE_STRINGS):
                    sys.stderr.write("Hap string %s  sample count %s misses GLM threshold %s\n"%(hap_string, len(hap_methylation_list), threshold))
       
        sorted_hap_block_list = sorted(hap_block_list)
        
        
        glm_input = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.glm_input.txt"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        
        glm_input_writer = open(glm_input, "wt")
        glm_input_writer.write("Sample_Index\tHapAlleleSum\tHapBlockIndex\tMethylation\n")
        
        sample_index = 0
        for hap_index in range(len(sorted_hap_block_list)): 
            hap_sum, hap_string = sorted_hap_block_list[hap_index]
            hap_methylation_list = hap_string_hash[hap_string]
            if (self.DEBUG_PROCESS_HAPLOTYPE_STRINGS_GLM):
                sys.stderr.write("hap block index %s hap_sum %s hap_string %s with %s samples \n"%(hap_index, hap_sum, hap_string, len(hap_methylation_list)))
            for methylation_value in hap_methylation_list:
                sample_index += 1
                glm_input_writer.write("Sample_%s\t%s\t%s\t%s\n"%(sample_index, hap_sum, hap_index, methylation_value))
            
        glm_input_writer.close()
        
        glm_output = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.glm_output.txt"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        
        glm_command = "glm_corsiv_hapsum.cc.R %s %s"%(glm_input, glm_output)
        
        if (self.DEBUG_PROCESS_HAPLOTYPE_STRINGS_GLM):
            sys.stderr.write("GLM command %s\n"%glm_command)
        os.system(glm_command)
        
        # extract the results and output a GLM model results
        glm_formal_record = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.glm"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        glm_formal_record_writer = open(glm_formal_record, "wt")
       
        r = open(glm_output, "rt")
        line = r.readline()
        line = r.readline()
        r.close()
        ff = line.strip().split(',')
        glm_formal_record_writer.write("CoRSIV\tHapName\tTissue\tGLM_pvalue\tGLM_R2\n")
        buffer = [corsiv_short, hap_block_name, gtex_tissue, ff[1], ff[2]]
        glm_formal_record_writer.write("%s\n"%"\t".join(buffer))
        
        glm_formal_record_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processHaplotypeStrings corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(),  corsiv_short, gtex_tissue, hap_block_name))

    
    def processHaplotypeInfo(self, out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue, methylation_tissue_data, sample_haplotype_info_list):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processHaplotypeInfo corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(), gtex_tissue, corsiv_short, hap_block_name))

        # fill hapString1.hapString2 and find the haplotypes with at least two instances
        hap_allele_string_hash = {}
        two_allele_hash = {}
        hap_reflexive_allele_hash = {}
        hap_012_allele_hash = {}
        
        for idx in range(len(sample_haplotype_info_list)):
            sample_haplotype_info = sample_haplotype_info_list[idx]
            allele1_string = sample_haplotype_info.haplotype_allele1_string
            allele1_sum = sample_haplotype_info.haplotype_allele1_sum
            if not (allele1_string in hap_allele_string_hash):
                hap_allele_info = CStruct(hap_allele_sum=allele1_sum, hap_allele_string = allele1_string, allele_string_idx=-1)
                hap_allele_string_hash[allele1_string]=hap_allele_info
                
            allele2_string = sample_haplotype_info.haplotype_allele2_string
            allele2_sum = sample_haplotype_info.haplotype_allele2_sum
            if not (allele2_string in hap_allele_string_hash):
                hap_allele_info = CStruct(hap_allele_sum=allele2_sum, hap_allele_string = allele2_string, allele_string_idx=-1)
                hap_allele_string_hash[allele2_string]=hap_allele_info
            
            methylation_ratio = methylation_tissue_data[idx]
            
            two_allele_string = "%s\n%s"%(allele1_string, allele2_string)
            if not (two_allele_string in two_allele_hash):
                two_allele_info = CStruct(two_allele_string=two_allele_string, count=0, two_allele_sum =allele1_sum + allele2_sum, meth_set=set())
                two_allele_hash[two_allele_string]  = two_allele_info 
                
            two_allele_hash[two_allele_string].count += 1
            two_allele_hash[two_allele_string].meth_set.add(methylation_ratio)
            
            # set up reflexive two allele haplotype
            hap_reflexive = "%s\n%s"%(allele1_string, allele2_string)
            if (allele2_string<allele1_string):
                hap_reflexive = "%s\n%s"%(allele2_string, allele1_string)
            sample_haplotype_info.reflexive_haplotype = hap_reflexive
            
            if not (hap_reflexive in hap_reflexive_allele_hash):
                hap_reflexive_info = CStruct(hap_reflexive = hap_reflexive, allele1_string=allele1_string, allele2_string = allele2_string,
                                             two_allele_sum= allele1_sum + allele2_sum, count=0, meth_set=set(), hap_reflexive_index=-1)
                hap_reflexive_allele_hash[hap_reflexive]=hap_reflexive_info
            
            hap_reflexive_allele_hash[hap_reflexive].count += 1
            hap_reflexive_allele_hash[hap_reflexive].meth_set.add(methylation_ratio)
            
            # setup the 012 haplotype hash
            haplotype_012_string = sample_haplotype_info.haplotype_012_string
            if not (haplotype_012_string in hap_012_allele_hash):
                hap_012_info = CStruct(haplotype_012_string = haplotype_012_string,
                                       allele1_string =allele1_string, allele2_string = allele2_string,
                                       two_allele_sum = allele1_sum + allele2_sum, count=0, meth_set=set(), hap_012_index =-1)
                hap_012_allele_hash[haplotype_012_string]=hap_012_info
            
            hap_012_allele_hash[haplotype_012_string].count += 1
            hap_012_allele_hash[haplotype_012_string].meth_set.add(methylation_ratio)
        
        if (self.DEBUG_PROCESS_HAPLOTYPE_INFO):
            sys.stderr.write("Single haplotype allele strings: %s\n"%" ".join(list(hap_allele_string_hash.keys())))
            sys.stderr.write("Reflexive haplotype allele strings: %s\n"%" ".join(list(hap_reflexive_allele_hash.keys())))
            sys.stderr.write("012 haplotype allele strings: %s\n"%" ".join(list(hap_012_allele_hash.keys())))
            
        # sort the allele string by allele sum and lexicographic order
        hap_allele_string_list = []
        for hap_allele_string in hap_allele_string_hash:
            hap_allele_info = hap_allele_string_hash[hap_allele_string]
            hap_allele_string_list.append( (hap_allele_info.hap_allele_sum, hap_allele_info.hap_allele_string))
        
        hap_allele_string_list.sort()
        for hap_allele_idx in range(len(hap_allele_string_list)):
            hap_allele_sum, hap_allele_string = hap_allele_string_list[hap_allele_idx]
            hap_allele_info = hap_allele_string_hash[hap_allele_string]
            hap_allele_info.allele_string_idx = hap_allele_idx
            
        if (self.DEBUG_PROCESS_HAPLOTYPE_INFO):
            sys.stderr.write("Sorted Single haplotype allele strings: %s\n"%" ".join( [str(x) for x in hap_allele_string_list]))

        two_allele_string_list = []        
        for two_allele_string in two_allele_hash:
            two_allele_info = two_allele_hash[two_allele_string]
            two_allele_string_list.append( (two_allele_info.two_allele_sum, two_allele_string))
            
        two_allele_string_list.sort()
        if (self.DEBUG_PROCESS_HAPLOTYPE_INFO):
            sys.stderr.write("Sorted two allele strings: %s\n"%" ".join( [str(x) for x in hap_allele_string_list]))
            
        # sort by allele sum, lexicographic the reflexive allele list
        hap_reflexive_string_list = []        
        for hap_reflexive in hap_reflexive_allele_hash:
            hap_reflexive_info = hap_reflexive_allele_hash[hap_reflexive]
            hap_reflexive_string_list.append( (hap_reflexive_info.two_allele_sum, hap_reflexive))
            
        hap_reflexive_string_list.sort()
        if (self.DEBUG_PROCESS_HAPLOTYPE_INFO):
            sys.stderr.write("Sorted reflexive haplotype strings: %s\n"%" ".join( [str(x) for x in hap_reflexive_string_list]))
        
        for hap_reflexive_index in range(len(hap_reflexive_string_list)):
            two_allele_sum, hap_reflexive = hap_reflexive_string_list[hap_reflexive_index]
            hap_reflexive_info = hap_reflexive_allele_hash[hap_reflexive]
            hap_reflexive_info.hap_reflexive_index = hap_reflexive_index
                    
        # sort by allele sum, lexicographic the 012 allele list
        hap_012_string_list = []        
        for haplotype_012_string in hap_012_allele_hash:
            hap_012_info = hap_012_allele_hash[haplotype_012_string]
            hap_012_string_list.append( (hap_012_info.two_allele_sum, haplotype_012_string) )
            
        hap_012_string_list.sort()
        if (self.DEBUG_PROCESS_HAPLOTYPE_INFO):
            sys.stderr.write("Sorted reflexive haplotype strings: %s\n"%" ".join( [str(x) for x in hap_012_string_list]))
                    
        for hap_012_index in range(len(hap_012_string_list)):
            two_allele_sum, haplotype_012_string = hap_012_string_list[hap_012_index]
            hap_012_info = hap_012_allele_hash[haplotype_012_string]
            hap_012_info.hap_012_index = hap_012_index
            
        # now setup R
        haplotype_info_glm_input = "%s%shap_info_glm_input.corsiv.%s.block.%s.tissue.%s.txt"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        
        haplotype_info_glm_input_writer = open(haplotype_info_glm_input, "wt")
        header = ["Sample", "allele1_index", "allele1_sum", "allele2_index", "allele2_sum", "allele_sum", "methylation", "reflexiveHaplotype_index", "haplotype_012_index"]
        haplotype_info_glm_input_writer.write("%s\n"%"\t".join(header))
        
        # setup also list of lists to print swarmplot
        for idx in range(len(sample_haplotype_info_list)):
            sample_haplotype_info = sample_haplotype_info_list[idx]
            allele1_string = sample_haplotype_info.haplotype_allele1_string
            allele1_info =hap_allele_string_hash[allele1_string]
            allele1_index = allele1_info.allele_string_idx
            allele1_sum = sample_haplotype_info.haplotype_allele1_sum
            
            allele2_string = sample_haplotype_info.haplotype_allele2_string
            allele2_info =hap_allele_string_hash[allele2_string]
            allele2_index = allele2_info.allele_string_idx
            allele2_sum = sample_haplotype_info.haplotype_allele2_sum
            
            two_allele_string = "%s\n%s"%(allele1_string, allele2_string)
            methylation_ratio = methylation_tissue_data[idx]
            
            reflexive_haplotype = sample_haplotype_info.reflexive_haplotype
            hap_reflexive_index = hap_reflexive_allele_hash[reflexive_haplotype].hap_reflexive_index
            
            haplotype_012_string = sample_haplotype_info.haplotype_012_string
            hap_012_index = hap_012_allele_hash[haplotype_012_string].hap_012_index
            
            if (hap_reflexive_allele_hash[reflexive_haplotype].count>1):
                buffer = ["Sample_%s"%idx, allele1_index, allele1_sum, allele2_index, allele2_sum, allele1_sum + allele2_sum, methylation_ratio, hap_reflexive_index, hap_012_index]
                haplotype_info_glm_input_writer.write("%s\n"%"\t".join([str(x) for x in buffer]))
                    
        haplotype_info_glm_input_writer.close()
        haplotype_info_glm_output = "%s%shap_info_glm_output.corsiv.%s.block.%s.tissue.%s.txt"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        
        glm_command = "glm_corsiv_allele.cc.R %s %s"%(haplotype_info_glm_input, haplotype_info_glm_output)
        if (self.DEBUG_PROCESS_HAPLOTYPE_INFO):
            sys.stderr.write("glm_command %s\n"%glm_command)
            
        os.system(glm_command)
        
        haplotype_info_glm_output_reader=open(haplotype_info_glm_output, "rt")
        line = haplotype_info_glm_output_reader.readline()
        line = haplotype_info_glm_output_reader.readline()
        haplotype_info_glm_output_reader.close()
        ff = line.strip().split(',')
        haplotype_info_glm_output_reader.close()
        
        haplotype_info_glm_report =  "%s%shap_info_glm_report.corsiv.%s.block.%s.tissue.%s.glminfo"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        haplotype_info_glm_report_writer = open(haplotype_info_glm_report, "wt")
       
        report_header = ["CoRSIV", "HapName", "Tissue", "glm0_R2", "glm0_p", "glm1_R2", "glm1_p", "glm2_R2", "glm2_p", "glm3_R2", "glm3_p", "glm4_R2", "glm4_p", "glm5_R2", "glm5_p"]
        haplotype_info_glm_report_writer.write("%s\n"%"\t".join(report_header))
        haplotype_info_glm_report_writer.write("%s\t%s\t%s\t%s\n"%(corsiv_short, hap_block_name, gtex_tissue, "\t".join(ff[1:])))
        haplotype_info_glm_report_writer.close()
        
        # plot swarmplot of two-allele plot
        list_of_lists = []
        for two_allele_pair in two_allele_string_list:
            two_allele_sum, two_allele_string = two_allele_pair
            two_allele_info = two_allele_hash[two_allele_string]
            if (two_allele_info.count>1):
                for methylation_ratio in two_allele_info.meth_set:
                    list_of_lists.append([two_allele_string, methylation_ratio])
        
        columns_hap_allele1_allele2_meth = ["Haplotype allele1 allele2", "CoRSIV Methylation"]
        hap_allele1_allele2_meth_data_frame = pd.DataFrame(list_of_lists,  columns = columns_hap_allele1_allele2_meth)
        
        plt.figure(figsize=(20,10))
        sn.set(style="white")
        pdf_swarmplot = "%s%shap_allele1_allele2_info.corsiv.%s.block.%s.tissue.%s.swarmplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.swarmplot(x="Haplotype allele1 allele2", y="CoRSIV Methylation", data=hap_allele1_allele2_meth_data_frame)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        plt.xticks(rotation=90)
        plt.savefig(pdf_swarmplot, bbox_inches='tight', dpi=150)

        # plot swarmplot of reflexive allele plot
        list_of_lists = []
        for hap_reflexive_pair in hap_reflexive_string_list:
            two_allele_sum, hap_reflexive = hap_reflexive_pair
            hap_reflexive_info = hap_reflexive_allele_hash[hap_reflexive]
            if (hap_reflexive_info.count>1):
                num_samples = len(hap_reflexive_info.meth_set)
                for methylation_ratio in hap_reflexive_info.meth_set:
#                    list_of_lists.append(["%s\n%s"%(str(two_allele_sum).zfill(4), hap_reflexive), methylation_ratio])
                    list_of_lists.append(["all_sum=%s n=%s\n%s"%(two_allele_sum, num_samples, hap_reflexive), methylation_ratio])
        
        columns_hap_reflexive_meth = ["Haplotype reflexive allele1 allele2", "CoRSIV Methylation"]
        hap_reflexive_meth_data_frame = pd.DataFrame(list_of_lists,  columns = columns_hap_reflexive_meth)
        
        plt.figure(figsize=(20,10))
        sn.set(style="white")
        pdf_swarmplot = "%s%shap_reflexive.corsiv.%s.block.%s.tissue.%s.swarmplot.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.swarmplot(x="Haplotype reflexive allele1 allele2", y="CoRSIV Methylation", data=hap_reflexive_meth_data_frame)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        ax.set(ylim=(0, 1))
        plt.xticks(rotation=90)
        plt.savefig(pdf_swarmplot, bbox_inches='tight', dpi=150)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::processHaplotypeInfo corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(), gtex_tissue, corsiv_short, hap_block_name))
    
    def processHaplotypeStrings(self, out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue,
                                         haplotype_sum, methylation_tissue_data, haplotype_string_list):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Run_CoRSIV_mQTL::processHaplotypeStrings corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(), gtex_tissue, corsiv_short, hap_block_name))

        # organize data and select SNP patterns w/ at least 5% of the data
        hap_string_hash = {}
        hap_string_to_hap_sum = {}
        for idx in range(len(haplotype_string_list)):
            hap_string = haplotype_string_list[idx]
            corsiv_methylation = methylation_tissue_data[idx]
            if not (hap_string in hap_string_hash):
                hap_string_hash[hap_string]=[]
            hap_string_hash[hap_string].append(corsiv_methylation)
            hap_sum = haplotype_sum[idx]
            hap_string_to_hap_sum[hap_string]=hap_sum
        
        threshold = int(0.02*len(haplotype_string_list))
        if (self.DEBUG_PROCESS_HAPLOTYPE_STRINGS):
            sys.stderr.write("Threshold: %s\n"%threshold)
        
        # add the glm processing
        self.processHaplotypeStringByGLM(out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue,
                                         hap_string_to_hap_sum, hap_string_hash)    
        
        list_of_lists = []
        range_list_of_lists = []
        
        hap_sum_array = []
        hap_range_array = []
        
        for hap_string in hap_string_hash:
            hap_methylation_list = hap_string_hash[hap_string]
            hap_sum = hap_string_to_hap_sum[hap_string]
            hap_sum_string = str(int(hap_string_to_hap_sum[hap_string])).zfill(4)
            hap_label = "%s_%s"%(hap_sum_string, hap_string)
            
            if (len(hap_methylation_list)>=threshold):
                if (self.DEBUG_PROCESS_HAPLOTYPE_STRINGS):
                    sys.stderr.write("Hap string %s  sample count %s passes threshold %s\n"%(hap_string, len(hap_methylation_list), threshold))
                for corsiv_methylation in hap_methylation_list:
                    list_of_lists.append([hap_label, corsiv_methylation])
                hap_string_range = max(hap_methylation_list)-min(hap_methylation_list)
                range_list_of_lists.append([hap_sum, hap_string_range])
                # prep range correlation 
                hap_sum_array.append(hap_sum)
                hap_range_array.append(hap_string_range)
            else:
                if (self.DEBUG_PROCESS_HAPLOTYPE_STRINGS):
                    sys.stderr.write("Hap string %s  sample count %s misses threshold %s\n"%(hap_string, len(hap_methylation_list), threshold))
                    
        columns_hap_meth = ["Haplotype snps", "CoRSIV Methylation"]
        hap_meth_data_frame = pd.DataFrame(list_of_lists,  columns = columns_hap_meth)
        
        columns_hap_meth_range = ["Haplotype sum", "CoRSIV Methylation Range"]
        hap_meth_range_data_frame = pd.DataFrame(range_list_of_lists,  columns = columns_hap_meth_range)
        
        # then compute correlations and plot swarmplots by string label (hopefully it sorts out alphabetically)
        # also save suporting data for correlation
        # finally collect range vs allele sum; plot as swarmplot
        plt.figure()
        sn.set(style="white")
        pdf_swarmplot = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.hap_swarmstring.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.swarmplot(x="Haplotype snps", y="CoRSIV Methylation", data=hap_meth_data_frame)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        plt.savefig(pdf_swarmplot, bbox_inches='tight', dpi=150)
        
        plt.figure()
        sn.set(style="white")
        pdf_swarmplot = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.hap_swarmrange.pdf"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        ax = sn.swarmplot(x="Haplotype sum", y="CoRSIV Methylation Range", data=hap_meth_range_data_frame)
        ax.set_title("%s %s %s"%(corsiv_short, gtex_tissue, hap_block_name))
        plt.savefig(pdf_swarmplot, bbox_inches='tight', dpi=150)
        
        pcc_value, pcc_pvalue = scipy.stats.pearsonr(hap_sum_array, hap_range_array)
        rho_value, rho_pvalue = scipy.stats.spearmanr(hap_sum_array, hap_range_array)
        
        correlation_report = "%s%smethylation_hap_corsiv.%s.block.%s.tissue.%s.correlstring"%(out_corsiv_folder, os.path.sep, corsiv_short, hap_block_name, gtex_tissue)
        correlation_report_writer = open(correlation_report, "wt")
        correlation_report_writer.write("CoRSIV\tTissue\tDonors\tHap_block\tPCC\tPCC_pvalue\tRho\tRho_pvalue\n")
        buffer = [corsiv_short, gtex_tissue, len(haplotype_sum), hap_block_name, pcc_value, pcc_pvalue, rho_value, rho_pvalue]
        correlation_report_writer.write("%s\n"%"\t".join([str(x) for x in buffer]))
        correlation_report_writer.close()
        
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Run_CoRSIV_mQTL::processHaplotypeStrings corsiv %s tissue %s hapblock %s\n"%(SimpleTime.now(), gtex_tissue, corsiv_short, hap_block_name))
    
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
            haplotype_string_list = ['']*len(snv_donor_col_idx)
            
            sample_haplotype_info_list = []
            
            for sample_idx_idx in range(len(snv_donor_col_idx)):
                snv_sample_idx = snv_donor_col_idx[sample_idx_idx]
                sample_haplotype_info = CStruct(sample_idx = snv_sample_idx,
                                          haplotype_allele1_sum = 0, haplotype_allele2_sum=0, haplotype_allele12_sum = 0,
                                          haplotype_allele1_string = "", haplotype_allele2_string="",
                                          haplotype_allele12_string="", haplotype_012_string="", reflexive_haplotype ="" )
                sample_haplotype_info_list.append(sample_haplotype_info)
                
            for snv_line in filtered_SNV_reader:
                snv_ff = snv_line.strip().split('\t')
                snv_id = snv_ff[2]
                buffer = [snv_id]
                for sample_idx_idx in range(len(snv_donor_col_idx)):
                    sample_info = sample_haplotype_info_list[sample_idx_idx]
                    
                    snv_sample_idx = snv_donor_col_idx[sample_idx_idx]
                    snv_string = snv_ff[snv_sample_idx]
                    if (snv_string=="0|0"):
                        haplotype_sum[sample_idx_idx] += 0
                        haplotype_string_list[sample_idx_idx] += "0"
                        
                        sample_info.haplotype_allele1_sum += 0
                        sample_info.haplotype_allele2_sum += 0
                        sample_info.haplotype_allele12_sum += 0
                        
                        sample_info.haplotype_allele1_string += "0"
                        sample_info.haplotype_allele2_string += "0"
                        sample_info.haplotype_012_string += "0"
                        
                    elif (snv_string=="0|1") or (snv_string=="1|0"):
                        haplotype_sum[sample_idx_idx] += 1 
                        haplotype_string_list[sample_idx_idx] += "1"
                        
                        if (snv_string=="0|1"):
                            sample_info.haplotype_allele1_sum += 0
                            sample_info.haplotype_allele2_sum += 1
                            sample_info.haplotype_allele1_string += "0"
                            sample_info.haplotype_allele2_string += "1"
                        else:
                            sample_info.haplotype_allele1_sum += 1
                            sample_info.haplotype_allele2_sum += 0
                            sample_info.haplotype_allele1_string += "1"
                            sample_info.haplotype_allele2_string += "0"
                            
                        sample_info.haplotype_allele12_sum += 1
                        sample_info.haplotype_012_string += "1"                        
                        
                    elif (snv_string=="1|1"):
                        haplotype_sum[sample_idx_idx] += 2 
                        haplotype_string_list[sample_idx_idx] += "2"
                    
                        sample_info.haplotype_allele1_sum += 1
                        sample_info.haplotype_allele2_sum += 1
                        sample_info.haplotype_allele12_sum += 2
                        
                        sample_info.haplotype_allele1_string += "1"
                        sample_info.haplotype_allele2_string += "1"
                        sample_info.haplotype_012_string += "2"
                    
                    else:
                        sys.stderr.write("Eek: snv %s snv_string %s sample_idx %s sample %s \n"%(snv_id, snv_string, snv_sample_idx, donors_with_methylation_and_SNV_list[sample_idx_idx]))
                        buffer.append("NA")
                        
            filtered_SNV_reader.close()
            
            
            support_haplotype_methylation_data_writer.write("%s_hap\t%s\n"%(corsiv_short, "\t".join([str(x) for x in haplotype_sum])))
            support_haplotype_methylation_data_writer.close()
            
            self.processHaplotypeAlleleSum(out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue,
                                           haplotype_sum, methylation_tissue_data)  
            self.processHaplotypeStrings(out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue,
                                         haplotype_sum, methylation_tissue_data, haplotype_string_list)  
            
            self.processHaplotypeInfo(out_corsiv_folder, corsiv_short, hap_block_name, gtex_tissue,
                                      methylation_tissue_data, sample_haplotype_info_list)  
                        
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
