#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "2.0"

import os, sys, argparse, re
import glob
import datetime
import math
import scipy.stats
from argparse import RawTextHelpFormatter
import random
import gzip
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sn
from numpy.polynomial.polynomial import polyfit

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class CaptureCORSiV_vs_GeneExpression:
    DEBUG_PROGRESS                  = True
    DEBUG_LOAD_SAMPLE_MAP           = False
    DEBUG_CLEANUP_COLUMNS           = False
    DEBUG_TISSUE_COVERAGE           = False
    DEBUG_MULTI_TISSUE_HEATMAPS     = False
    DEBUG_MAX_WIR_VERBOSE           = False
    DEBUG_MAX_WIR_SHORT             = False
    DEBUG_TISSUE_PAIR_RHO           = False
    DEBUG_TISSUE_PAIR_RHO_VERBOSE   = False
    DEBUG_TISSUE_PAIR_RHO_SUPER_VERBOSE   = False
    DEBUG_TISSUE_PAIR_SCATTERPLOT   = False
    DEBUG_TISSUE_PAIR_CORSIV_REPORT = False
    DEBUG_CORSIV_RANGE              = False
    DEBUG_CORSIV_GENE_PAIRS         = False
    DEBUG_LOAD_GTEX_SAMPLE_MAP      = False
    DEBUG_LOAD_GTEX_GENE_EXPRESSION = False
    DEBUG_PCC_RNA_METH              = False
    DEBUG_PCA_RNA_METH_CORSIV       = False
    DEBUG_PCC_RNA_array             = False
    DEBUG_PCC_METH_array            = False
    DEBUG_REPORT_METH_RNA_CORR      = False
    DEBUG_PCC_RNA_METH_PROGRESS     = True
    DEBUG_REPORT_GENE_NOT_FOUND     = False
    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.
            
determine methylation in GTEx capture vs GTEx gene expression correlation  
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-x','--corsivMethylationMatrix',   help='data matrix with beta values for CORSiVs',    required=True)
        parser.add_argument('-S','--tissueMap',                 help='map for samples/tissues/capture files',       required=True)
        parser.add_argument('-d','--minDonorsCorrelation',      help='minimum donors for correlation',              required=True)
        parser.add_argument('-g','--gtexGeneExpression',        help='GTEx gene expression',                        required=True)
        parser.add_argument('-s','--gtexSampleDef',             help='GTEx seq sample description',                 required=True)
        parser.add_argument('-C','--corsivDef',                 help='capture corsiv definition',                   required=True)
        parser.add_argument('-G','--geneModel',                 help='GENCODE gene model',                        required=True)
        parser.add_argument('-R','--geneModelRadius',           help='radius around GENCODE gene model',          required=True)
        parser.add_argument('-o','--outputRoot',                help='output file root',                          required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def annotatedCorsivsWithGenes(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::annotatedCorsivsWithGenes\n"%SimpleTime.now())
        
        self.corsiv_to_genes_hash={}
        
        corsiv_gene_neighbors_file = "%s.corsiv_gene_pairs.bed"%self.myArgs.outputRoot
        bedtools_command = "bedtools window -w %s -a %s -b %s > %s"%(self.myArgs.geneModelRadius, self.myArgs.corsivDef, self.myArgs.geneModel, corsiv_gene_neighbors_file)
        if (self.DEBUG_CORSIV_GENE_PAIRS):
            sys.stderr.write("bedtools command %s\n"%bedtools_command)
       
        os.system(bedtools_command)        
        corsiv_gene_neighbors_file_reader = open(corsiv_gene_neighbors_file, "rt")
        for line in corsiv_gene_neighbors_file_reader:
            ff = line.strip().split('\t')
            corsiv = ff[3]
            corsiv_gene = ff[10]
            if not (corsiv in self.corsiv_to_genes_hash):
                self.corsiv_to_genes_hash[corsiv]=set()
                
            self.corsiv_to_genes_hash[corsiv].add(corsiv_gene)
            
            if (self.DEBUG_CORSIV_GENE_PAIRS):
                sys.stderr.write("CG Pair: %s -> %s\n"%(corsiv, corsiv_gene))
        
        corsiv_gene_neighbors_file_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::annotatedCorsivsWithGenes\n"%SimpleTime.now())
        
        
    def loadTissueMap(self):    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::loadTissueMap\n"%SimpleTime.now())
        
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
    
  
   
                    
    def loadGTExSampleMap(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::loadGTExSampleMap\n"%SimpleTime.now())
    
        self.gene_tissue_hash = {}
        
        
        rna_map_reader = open(self.myArgs.gtexSampleDef, "rt")
        
        for line in rna_map_reader:
            if (line.find("RNASEQ")<0):
                continue
            
            ff = line.strip().split('\t')
            rna_sample_id = ff[0]
            gtex_tissue_id = ff[6]
            gg = rna_sample_id.split('-')
            gtex_donor_id = "-".join(gg[:2])
            
            if (self.DEBUG_LOAD_GTEX_SAMPLE_MAP):
                sys.stderr.write("rna_sample_id %s donor %s tissue %s\n"%(rna_sample_id, gtex_donor_id, gtex_tissue_id))
            
            tissue_info = CStruct(rna_sample_id = rna_sample_id, gtex_donor_id= gtex_donor_id, gtex_tissue_id=gtex_tissue_id)
            
            if not (gtex_tissue_id in self.gene_tissue_hash):
                self.gene_tissue_hash[gtex_tissue_id] = set()
            
            self.gene_tissue_hash[gtex_tissue_id].add(tissue_info)
        
        rna_map_reader.close()
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::loadGTExSampleMap\n"%SimpleTime.now())
        
    def loadGTExGeneExpression(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::loadGTExGeneExpression\n"%SimpleTime.now())

        self.gtex_gene_expression_df = pd.read_csv(self.myArgs.gtexGeneExpression, sep='\t', index_col=0, header=0)
        
        self.gtex_gene_columns_hash = {}
        self.gtex_gene_columns = list(self.gtex_gene_expression_df.columns)
        
        for col_idx in range(len(self.gtex_gene_columns)):
            rna_sample = self.gtex_gene_columns[col_idx]
            self.gtex_gene_columns_hash[rna_sample] = col_idx
            if (self.DEBUG_LOAD_GTEX_GENE_EXPRESSION):
                sys.stderr.write("rna col_idx %s rna sample %s -> %s\n"%(col_idx, rna_sample, self.gtex_gene_columns_hash[rna_sample]))
           
        self.gtex_gene_symbols_to_rowidx_hash = {}
        self.gtex_gene_index = list(self.gtex_gene_expression_df.index)
        
        for gene_idx in range(len(self.gtex_gene_index)):
            gene_symbol = self.gtex_gene_expression_df.iat[gene_idx, 0]
            self.gtex_gene_symbols_to_rowidx_hash[gene_symbol] = gene_idx
            if (self.DEBUG_LOAD_GTEX_GENE_EXPRESSION):
                sys.stderr.write("rna row_idx %s ensid %s symbol %s --> %s\n"%(gene_idx, self.gtex_gene_index[gene_idx], gene_symbol, self.gtex_gene_symbols_to_rowidx_hash[gene_symbol]))
                   
        # now setup for rna tissue to donor map
        
        self.rna_tissue_to_donor_hash = {}
        for rna_tissue in self.gene_tissue_hash:
            self.rna_tissue_to_donor_hash[rna_tissue]={}
            if (self.DEBUG_LOAD_GTEX_GENE_EXPRESSION):
                sys.stderr.write("building donor & rna tissue sequenced map for rna tissue %s\n"%rna_tissue)
                
            for donor_info in self.gene_tissue_hash[rna_tissue]:
                donor_id = donor_info.gtex_donor_id
                rna_sample = donor_info.rna_sample_id
                
                if (rna_sample in self.gtex_gene_columns_hash):
                    self.rna_tissue_to_donor_hash[rna_tissue][donor_id]=self.gtex_gene_columns_hash[rna_sample]
                    if (self.DEBUG_LOAD_GTEX_GENE_EXPRESSION):
                        sys.stderr.write("donor %s rna_sample %s with coverage: column %i\n"%(donor_id, rna_sample, self.gtex_gene_columns_hash[rna_sample]))
            
            if (self.DEBUG_LOAD_GTEX_GENE_EXPRESSION):
                sys.stderr.write("finished donor & rna tissue sequenced map for rna tissue %s with %s pairs\n"%(rna_tissue, len(self.rna_tissue_to_donor_hash[rna_tissue])))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::loadGTExGeneExpression\n"%SimpleTime.now())
    
    def build_list_idx_hash(self, string_list):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::build_list_idx_hash\n"%SimpleTime.now())
            
        list_idx_hash = {}
        
        for idx in range(len(string_list)):
            current_string = string_list[idx]
            list_idx_hash[current_string]=idx
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::build_list_idx_hash\n"%SimpleTime.now())
        return list_idx_hash
    
    def computeMethylationGeneExpCorrelation(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::computeMethylationGeneExpCorrelation\n"%SimpleTime.now())
        
        self.min_donors_correlation = int(self.myArgs.minDonorsCorrelation)
        if (self.min_donors_correlation<10):
            sys.stderr.write("Need at least 10 donors to compute correlations !\nExiting ...\n")
            sys.exit(4)
            
            
        self.meth_tissues = len(self.tissue_to_donor_hash)
        self.rna_tissues  = len(self.gene_tissue_hash)
        
        
        coverage_index = list(self.coverage_df.index)
        if (self.DEBUG_PCC_RNA_METH):
            sys.stderr.write("Top 5 meth corsivs %s\n"%"\n".join(coverage_index[:5]))
            
        self.capture_corsiv_to_row_idx_hash = self.build_list_idx_hash(coverage_index)

        self.capture_corsiv_to_short_corsiv_hash = {}
        
        for capture_corsiv in coverage_index:
            corsiv_split_1 = capture_corsiv.split(';')
            corsiv_split_2 = corsiv_split_1[3].split('.')
            short_corsiv = '.'.join(corsiv_split_2[:3])
            if (short_corsiv in self.corsiv_to_genes_hash):
                genes = ";".join(list(self.corsiv_to_genes_hash[short_corsiv]))
                self.capture_corsiv_to_short_corsiv_hash[capture_corsiv]=(short_corsiv, self.corsiv_to_genes_hash[short_corsiv])
            else:
                genes = None
                
            if (self.DEBUG_PCC_RNA_METH):
                sys.stderr.write("capture_corsivs %s corsiv_split1 %s corsiv_split2 %s short corsiv %s genes %s \n"%(capture_corsiv, corsiv_split_1[3], ";".join(corsiv_split_2), short_corsiv, genes))
            
        self.meth_rna_pair_corr_hash = {}
        
        meth_tissue_list = sorted( list(self.tissue_to_donor_hash.keys()))
        rna_tissue_list = sorted(list(self.gene_tissue_hash.keys()))
        
        self.final_corsiv_set = set()
        
        if (self.DEBUG_PCC_RNA_METH):
            sys.stderr.write("meth tissues: %s: %s\n rna tissues %s : %s\n"%(self.meth_tissues, ";".join(meth_tissue_list), self.rna_tissues, ";".join(rna_tissue_list)))
            
        # iterate over methylation tissues
        for meth_tissue_idx in range(len(meth_tissue_list)):
            meth_tissue_id = meth_tissue_list[meth_tissue_idx]
            # find all donors for each tissue
            meth_tissue_donors = set(self.tissue_to_donor_hash[meth_tissue_id].keys())
            if (self.DEBUG_PCC_RNA_METH):
                sys.stderr.write("Processing meth tissue %s %s with %s donors\n"%(meth_tissue_idx, meth_tissue_id, len(meth_tissue_donors)))
                
            # iterate over gene expression tissue
            for rna_tissue_idx in range(len(rna_tissue_list)):
                rna_tissue_id = rna_tissue_list[rna_tissue_idx]
                rna_tissue_donors = set(self.rna_tissue_to_donor_hash[rna_tissue_id].keys())
                
                # find common donors
                common_meth_rna_donors = meth_tissue_donors.intersection(rna_tissue_donors)
                
                if (self.DEBUG_PCC_RNA_METH):
                    sys.stderr.write("Assesing pair meth %s %s with %s donors vs rna %s %s with %s donors common %s donors\n"%(meth_tissue_idx, meth_tissue_id, len(meth_tissue_donors), rna_tissue_idx, rna_tissue_id, len(rna_tissue_donors), len(common_meth_rna_donors)))
                    
                    
                if (common_meth_rna_donors<self.min_donors_correlation):
                    if (self.DEBUG_PCC_RNA_METH_PROGRESS):
                        sys.stderr.write("Skipping pair meth tissue %s rna tissue %s with fewer than %s donors: %s\n"%(meth_tissue_id, rna_tissue_id, self.min_donors_correlation, len(common_meth_rna_donors)))
                    continue
                else:
                    if (self.DEBUG_PCC_RNA_METH_PROGRESS):
                        sys.stderr.write("Processing pair meth tissue %s rna tissue %s with more than %s donors: %s\n"%(meth_tissue_id, rna_tissue_id, self.min_donors_correlation, len(common_meth_rna_donors)))
                tissue_pair = "METH_%s_RNA_%s"%(meth_tissue_id, rna_tissue_id) 
                self.meth_rna_pair_corr_hash[tissue_pair]={}
                
                # setup columns for common donors in beta values df and in gene expression df
                common_meth_donors_idx = []
                common_rna_donors_idx = []
                
                # find their column id in each dataset
                for meth_and_rna_donor in common_meth_rna_donors:
                    meth_df_col_idx = self.tissue_to_donor_hash[meth_tissue_id][meth_and_rna_donor]
                    common_meth_donors_idx.append(meth_df_col_idx)
                    
                    rna_df_col_idx = self.rna_tissue_to_donor_hash[rna_tissue_id][meth_and_rna_donor]
                    common_rna_donors_idx.append(rna_df_col_idx)
   
                # iterate over corsivs
                for capture_corsiv in self.capture_corsiv_to_short_corsiv_hash:
                    short_corsiv, gene_list = self.capture_corsiv_to_short_corsiv_hash[capture_corsiv]
                    for gene_symbol in gene_list:
                        if (self.DEBUG_PCA_RNA_METH_CORSIV):
                            sys.stderr.write("Trio: M %s R %s C %s gene %s\n"%(meth_tissue_id, rna_tissue_id, capture_corsiv, gene_symbol))
                        
                        if not (gene_symbol in self.gtex_gene_symbols_to_rowidx_hash):
                            if (self.DEBUG_REPORT_GENE_NOT_FOUND):
                                sys.stderr.write("Corsiv gene %s not found in gene expression\n"%gene_symbol)
                            continue
                        
                        final_corsiv_id = "corsiv_%s_gene_%s"%(short_corsiv, gene_symbol)
                        self.final_corsiv_set.add(final_corsiv_id)
                        # count corsivs w/ methylation                    
                        meth_corsiv_row_idx = self.capture_corsiv_to_row_idx_hash[capture_corsiv]
                        rna_gene_row_idx = self.gtex_gene_symbols_to_rowidx_hash[gene_symbol]
                        
                        meth_value_meth_col_idx = []
                        meth_value_rna_col_idx = []
                        
                        for running_common_donor_idx in range(len(common_meth_donors_idx)):
                            meth_col_idx = common_meth_donors_idx[running_common_donor_idx]
                            meth_is_na_flag = np.isnan(self.coverage_df.iat[meth_corsiv_row_idx, meth_col_idx])
                            
                            if  not meth_is_na_flag:
                                meth_value_meth_col_idx.append(meth_col_idx)
                                meth_value_rna_col_idx.append(common_rna_donors_idx[running_common_donor_idx])
                            
                        if (len(meth_value_meth_col_idx)<self.min_donors_correlation):
                            # if < 20: -2
                            self.meth_rna_pair_corr_hash[tissue_pair][final_corsiv_id]=-2
                            if (self.DEBUG_PCA_RNA_METH_CORSIV):
                                sys.stderr.write("ShortTrio: M %s R %s C %s gene %s values %s\n"%(meth_tissue_id, rna_tissue_id, capture_corsiv, gene_symbol, len(meth_value_meth_col_idx)))
                        else:
                        # select appropriate arrays
                            if (self.DEBUG_PCA_RNA_METH_CORSIV):
                                sys.stderr.write("PassTrio: M %s R %s C %s gene %s values %s\n"%(meth_tissue_id, rna_tissue_id, capture_corsiv, gene_symbol, len(meth_value_meth_col_idx)))
                       
                            if (self.DEBUG_PCC_RNA_array):
                                sys.stderr.write("Gene cols: %s\n"%";".join([str(x) for x in sorted(meth_value_rna_col_idx)]))
                            
                            if (self.DEBUG_PCC_METH_array):
                                sys.stderr.write("METH cols: %s\n"%";".join([str(x) for x in sorted(meth_value_meth_col_idx)]))
                                
                            rna_values_array = self.gtex_gene_expression_df.iloc[rna_gene_row_idx, meth_value_rna_col_idx]
                            meth_values_array = self.coverage_df.iloc[meth_corsiv_row_idx, meth_value_meth_col_idx]
                            # compute correlation
                            pcc_value, pcc_pvalue = scipy.stats.pearsonr(rna_values_array, meth_values_array)
                            if (pcc_pvalue<0.05):
                                # if significant > yes
                                self.meth_rna_pair_corr_hash[tissue_pair][final_corsiv_id]=pcc_value
                            else:
                                # if not significant 0
                                self.meth_rna_pair_corr_hash[tissue_pair][final_corsiv_id]=0
                                
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::computeMethylationGeneExpCorrelation\n"%SimpleTime.now())
        
    def reportAndVisualizeMethRNACorrelation(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::report_and_visualize_rna_meth_corr\n"%SimpleTime.now())

        # write data frame
        meth_rna_corr_file = "%s.methylation_rna_correlation.xls"%self.myArgs.outputRoot
        
        meth_rna_corr_file_writer = open(meth_rna_corr_file, "wt")
        
        corsiv_header = "Corsiv/Gene"
        meth_rna_pair_list = sorted(self.meth_rna_pair_corr_hash.keys())
        final_corsiv_list = sorted(list(self.final_corsiv_set))
        meth_rna_corr_file_writer.write("%s\t%s\n"%(corsiv_header,'\t'.join(meth_rna_pair_list)))
                                        
        for final_corsiv in final_corsiv_list:
            if (self.DEBUG_REPORT_METH_RNA_CORR):
                sys.stderr.write("outputting final corsiv %s\n"%final_corsiv)
            buffer = [final_corsiv]
            for meth_rna_pair in meth_rna_pair_list:
                final_corsiv_hash = self.meth_rna_pair_corr_hash[meth_rna_pair]
                buffer.append(final_corsiv_hash[final_corsiv])
                
            meth_rna_corr_file_writer.write("%s\n"%'\t'.join([str(x) for x in buffer]))
                
        meth_rna_corr_file_writer.close()
        
        # load data frame
        meth_rna_corr_df = pd.read_csv(meth_rna_corr_file, header=0, index_col=0, sep='\t')
        sys.stderr.write("meth_rna_corr_df shape %s\n"%str(meth_rna_corr_df.shape))
        # plot rna/meth corelation using cmap, euclidian and correlation
        pdf_heatmap_file = "%s.meth_rna_corr_heatmap.pdf"%self.myArgs.outputRoot
        jpg_heatmap_file = "%s.meth_rna_corr_heatmap.jpg"%self.myArgs.outputRoot
        
        plt.figure()
        sn.set()
#        sn.clustermap(meth_rna_corr_df,  cmap='vlag', col_colors=None, row_cluster=True, col_cluster=True, metric="euclidean", linewidths=0.0, rasterized=True, figsize=(30,30))
        sn.clustermap(meth_rna_corr_df,  cmap='vlag', col_colors=None, row_cluster=True, col_cluster=True, metric="correlation", linewidths=0.0, rasterized=True, figsize=(30,30))
        plt.savefig(pdf_heatmap_file, dpi=150)
        plt.savefig(jpg_heatmap_file, dpi=150)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::report_and_visualize_rna_meth_corr\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START CaptureCORSiV_vs_GeneExpression::work\n"%SimpleTime.now())
    
        self.annotatedCorsivsWithGenes() # bedtools window your friend
        
        self.loadTissueMap()
        self.loadBetaValues()
        
        self.loadGTExSampleMap()
        self.loadGTExGeneExpression()
        
        self.computeMethylationGeneExpCorrelation()
        self.reportAndVisualizeMethRNACorrelation()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP CaptureCORSiV_vs_GeneExpression::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = CaptureCORSiV_vs_GeneExpression.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = CaptureCORSiV_vs_GeneExpression(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
