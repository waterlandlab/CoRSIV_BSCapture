#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

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

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class SummaryCorsivBetaGeneExpression:
    DEBUG_PROGRESS          = True
    DEBUG_LOAD_CORR         = True
    DEBUG_BUILD_SUMMARY     = True
    DEBUG_OUTPUT_SUMMARY    = True
    DEBUG_SETUP_TISSUE_CENTRIC_REPORT   = True
    DEBUG_DO_TISSUE_CENTRIC_REPORT      = True
    DEBUG_OUTPUT_TISSUE_CENTRIC_REPORT  = True
    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

Summarize DNA methylation vs RNA expression correlation at corsivs
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-c','--correlationMatrix',     help='data matrix with  CoRSiVs methylation/expression correlation',    required=True)
        parser.add_argument('-R','--minimumR2',             help='minimum R^2 consider',    required=True)
        parser.add_argument('-o','--outputRoot',            help='output file root',                            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def loadCoverageValues(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::loadCoverageValues\n"%SimpleTime.now())
    
        self.df_coverage = pd.read_csv(self.myArgs.corsivCoverageCounts, sep="\t", index_col=0, header=0)
        if (self.DEBUG_LOAD_COVERAGE_COUNTS):
            sys.stderr.write("Loaded coverage counts values dim %s \n"%str(self.df_coverage.shape))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::setupTissueDeltaDataFrame\n"%SimpleTime.now())
            
        corsiv_index = list(self.df_coverage.index)
        data_frame_list_of_lists = []
        corsiv_columns = list(self.df_coverage.columns)
        
        for row_idx in range(len(corsiv_index)):
            corsiv_name = corsiv_index[row_idx]
            bits = corsiv_name.split('_')
            bp_length = int(bits[2])-int(bits[1])+1
            
            if (self.DEBUG_LOAD_COVERAGE_COUNTS):
                sys.stderr.write("Setting up BP coverage for corsiv %s bounds %s %s length %s\n"%(corsiv_name, bits[1], bits[2], bp_length))
            
            for col_idx in range(2, len(corsiv_columns)):
                tissue = corsiv_columns[col_idx]
                tissue_bits = tissue.split('_')
                counts = self.df_coverage.iat[row_idx, col_idx]
                bp_coverage = float(counts)*150/float(bp_length)
                if (self.DEBUG_LOAD_COVERAGE_COUNTS):
                    sys.stderr.write("Col [%s] %s reads %s bp_coverage  %s\n"%(col_idx, tissue, counts, bp_coverage))
                
                data_frame_list_of_lists.append([tissue_bits[0], bp_coverage])
                    
        # setup and load delta beta data frame
        columns_bp_coverage = ["Tissue", "BPCoverage"]
        self.bp_coverage_data_frame = pd.DataFrame( data_frame_list_of_lists, columns = columns_bp_coverage)
        
        if (self.DEBUG_LOAD_COVERAGE_COUNTS):
            sys.stderr.write("BP coverage matrix shape %s columns %s\n"%(str(self.bp_coverage_data_frame.shape), ";".join(list(self.bp_coverage_data_frame.columns))))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::loadCoverageValues\n"%SimpleTime.now())
    
    def outputViolinPlots(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::outputViolinPlots\n"%SimpleTime.now())
        
        violin_plot_pdf = "%s.violin_plot.pdf"%self.myArgs.outputRoot
        
        sn.set(style="whitegrid")
        ax = sn.violinplot(x="Tissue", y="BPCoverage", data=self.bp_coverage_data_frame)
        #ax.set_ylim(-1,1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)
    
        violin_plot_pdf = "%s.violin_plot_0_100x.pdf"%self.myArgs.outputRoot
        ax = sn.violinplot(x="Tissue", y="BPCoverage", data=self.bp_coverage_data_frame)
        ax.set_ylim(0, 100)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)

        violin_plot_pdf = "%s.violin_plot_0_200x.pdf"%self.myArgs.outputRoot
        ax = sn.violinplot(x="Tissue", y="BPCoverage", data=self.bp_coverage_data_frame)
        ax.set_ylim(0, 200)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)
                    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::outputViolinPlots\n"%SimpleTime.now())
    
    
    def loadCorrelationMatrix(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::loadCorrelationMatrix\n"%SimpleTime.now())
    
        self.df_corr = pd.read_csv(self.myArgs.correlationMatrix, sep="\t", index_col=0, header=0)
        
        if (self.DEBUG_LOAD_CORR):
            sys.stderr.write("Loaded correlation values dim %s \n"%str(self.df_corr.shape))
        
        self.corr_columns = list(self.df_corr.columns)
        self.corr_index = list(self.df_corr.index)
       
        self.col_idx_to_rna_array  = []
        self.rna_labels = {}
        self.corr_df_col_to_idx_hash = {}
        # setup col to rna tissue array
        
        for col_idx in range(len(self.corr_columns)):
            current_column = self.corr_columns[col_idx]
            rna_start = current_column.find("RNA_")
            rna_string = current_column[rna_start+4:]
            if (self.DEBUG_LOAD_CORR):
                sys.stderr.write("Column [%s] : %s RNA_start %s RNA_string %s\n"%(col_idx, current_column, rna_start, rna_string))
                
            self.col_idx_to_rna_array.append(rna_string)
            self.rna_labels[rna_string]=1
            
            self.corr_df_col_to_idx_hash[current_column]=col_idx
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::loadCorrelationMatrix\n"%SimpleTime.now())
        
    def summarizeCounts(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::summarizeCounts\n"%SimpleTime.now())
        
        
        self.hash_counts = {}
        self.hash_counts["ALL_positive_corr"]=[]
        self.hash_counts["ALL_negative_corr"]=[]
        for rna_label in self.rna_labels:
            rna_pos_string = "%s_positive_corr"%rna_label
            rna_neg_string = "%s_negative_corr"%rna_label
            self.hash_counts[rna_pos_string]=[]
            self.hash_counts[rna_neg_string]=[]
        
        for corsiv_idx in range(len(self.corr_index)):
            current_corsiv = self.corr_index[corsiv_idx]
            
            corsiv_corr_count = {}
            for rna_key in self.hash_counts.keys():
                corsiv_corr_count[rna_key]=0
            
            for col_idx in range(len(self.corr_columns)):
                current_column = self.corr_columns[col_idx]
                rna_string = self.col_idx_to_rna_array[col_idx]
                corr_value = float(self.df_corr.iat[corsiv_idx, col_idx])
                
                if (corr_value>=-1 and corr_value<-self.min_abs_pcc):
                    corsiv_corr_count["ALL_negative_corr"] +=1
                elif (corr_value>self.min_abs_pcc and corr_value<=1):
                    corsiv_corr_count["ALL_positive_corr"] +=1
            
            if (self.DEBUG_BUILD_SUMMARY):
                if (corsiv_corr_count["ALL_positive_corr"]>0 or corsiv_corr_count["ALL_negative_corr"]>0):
                    sys.stderr.write("Yay: Corsiv [%s]: %s  ALL_positive_corr %s ALL_negative_corr %s\n"%(corsiv_idx, current_corsiv, corsiv_corr_count["ALL_positive_corr"], corsiv_corr_count["ALL_negative_corr"]))
                else:
                    sys.stderr.write("Nope: Corsiv [%s]: %s  ALL_positive_corr %s ALL_negative_corr %s\n"%(corsiv_idx, current_corsiv, corsiv_corr_count["ALL_positive_corr"], corsiv_corr_count["ALL_negative_corr"]))
            
            for rna_key in self.hash_counts.keys():
                if (corsiv_corr_count[rna_key]>0):
                    self.hash_counts[rna_key].append(corsiv_corr_count[rna_key])
                    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::summarizeCounts\n"%SimpleTime.now())
        
    def outputSummary(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::work\n"%SimpleTime.now())
        
        for rna_key in self.hash_counts:
            corr_counts = self.hash_counts[rna_key]
            if (len(corr_counts)>0):
                sys.stderr.write("rna_key %s count>0 %s\n"%(rna_key, len(corr_counts)))
           
            plt.figure()
            sn.set()
            hist_plot_pdf = "%s.hist_plot.%s.pdf"%(self.myArgs.outputRoot, rna_key)
            plt.hist(np.asfarray(corr_counts), bins=50)
            plt.savefig(hist_plot_pdf, bbox_inches='tight', dpi=150)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::work\n"%SimpleTime.now())
    
    def setupTargetSurrogateCalc(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::setupTargetSurrogateCalc\n"%SimpleTime.now())
    
        self.tissue_centric_corr_hash = {}
        self.surrogate_tissues = ["Whole_Blood", "Skin"]
        
        self.labels_hash = {}
        for surrogate_tissue in self.surrogate_tissues:
            self.labels_hash[surrogate_tissue]={}
            
        self.report_header = []
        positive_corr_target_label= "Positive RNA/Meth correlation"
        self.report_header.append(positive_corr_target_label)
        self.labels_hash["positive_corr_target_label"]=positive_corr_target_label
        
        for surrogate_tissue in self.surrogate_tissues:
            positive_corr_target_positive_corr_surrogate_label = "Surrogate %s positive_corr/target positive_corr/surrogate"%surrogate_tissue
            self.labels_hash[surrogate_tissue]["positive_corr_target_positive_corr_surrogate_label"]=positive_corr_target_positive_corr_surrogate_label
            self.report_header.append(positive_corr_target_positive_corr_surrogate_label)
            
            positive_corr_target_negative_corr_surrogate_label = "Surrogate %s positive_corr/target negative_corr/surrogate"%surrogate_tissue
            self.labels_hash[surrogate_tissue]["positive_corr_target_negative_corr_surrogate_label"]=positive_corr_target_negative_corr_surrogate_label
            self.report_header.append(positive_corr_target_negative_corr_surrogate_label)
        
        negative_corr_target_label= "Negative RNA/Meth correlation"
        self.report_header.append(negative_corr_target_label)
        self.labels_hash["negative_corr_target_label"]=negative_corr_target_label
                       
        for surrogate_tissue in self.surrogate_tissues:
            negative_corr_target_positive_corr_surrogate_label = "Surrogate %s negative_corr/target positive_corr/surrogate"%surrogate_tissue
            self.labels_hash[surrogate_tissue]["negative_corr_target_positive_corr_surrogate_label"]=negative_corr_target_positive_corr_surrogate_label
            self.report_header.append(negative_corr_target_positive_corr_surrogate_label)
            
            negative_corr_target_negative_corr_surrogate_label = "Surrogate %s negative_corr/target negative_corr/surrogate"%surrogate_tissue
            self.labels_hash[surrogate_tissue]["negative_corr_target_negative_corr_surrogate_label"]=negative_corr_target_negative_corr_surrogate_label
            self.report_header.append(negative_corr_target_negative_corr_surrogate_label)
        
        self.venn_hash = {}
        self.venn_labels = ["Target", "Whole_Blood", "Skin","Two_Surrogates"]
            
        for rna_tissue in self.rna_labels:
            pos_corr_label = "PosCorr_%s"%rna_tissue
            neg_corr_label = "NegCorr_%s"%rna_tissue
            self.venn_hash[pos_corr_label]={}
            self.venn_hash[neg_corr_label]={}
            for venn_label in self.venn_labels:
                self.venn_hash[pos_corr_label][venn_label]=0
                self.venn_hash[neg_corr_label][venn_label]=0
                
            self.tissue_centric_corr_hash[rna_tissue] = {}
            for label in self.report_header:
                self.tissue_centric_corr_hash[rna_tissue][label]=0
                if (self.DEBUG_SETUP_TISSUE_CENTRIC_REPORT):
                    sys.stderr.write("Setup report for %s x %s\n"%(rna_tissue, label))
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::tissueCentricCounts\n"%SimpleTime.now())
    
    def doTargetSurrogateReport(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::doTargetSurrogateReport\n"%SimpleTime.now())
            
        # iterate over each corsiv/gene
        for corsiv_idx in range(len(self.corr_index)):
            current_corsiv = self.corr_index[corsiv_idx]
            
            if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                sys.stderr.write("Processing corsiv %s\n"%current_corsiv)
            
            # for each rna tissue
            for rna_tissue in self.rna_labels:
                
                pos_corr_label = "PosCorr_%s"%rna_tissue
                neg_corr_label = "NegCorr_%s"%rna_tissue
                
                target_meth_rna_corr_string = "METH_%s_RNA_%s"%(rna_tissue, rna_tissue)
                target_meth_rna_corr_col_idx = self.corr_df_col_to_idx_hash[target_meth_rna_corr_string]
                target_meth_rna_corr = self.df_corr.iat[corsiv_idx, target_meth_rna_corr_col_idx]
                
                if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                    sys.stderr.write("TT_corr: Corsiv %s Target/RNA %s pair_string %s col_idx %s check column %s corr %s\n"%(current_corsiv, rna_tissue, target_meth_rna_corr_string, target_meth_rna_corr_col_idx, self.corr_columns[target_meth_rna_corr_col_idx], target_meth_rna_corr))
                
                # check target rna/meth correlation
                if (target_meth_rna_corr>self.min_abs_pcc) and (target_meth_rna_corr<=1):
                    report_label = self.labels_hash["positive_corr_target_label"]
                    self.tissue_centric_corr_hash[rna_tissue][report_label]+=1
                
                    if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                        sys.stderr.write("Positive_target_corr\n")
                
                    venn_support_set = set()
                   
                    for surrogate_tissue in self.surrogate_tissues:
                        surrogate_meth_target_rna_corr_string = "METH_%s_RNA_%s"%(surrogate_tissue, rna_tissue)
                        surrogate_meth_target_rna_corr_col_idx = self.corr_df_col_to_idx_hash[surrogate_meth_target_rna_corr_string]
                        surrogate_meth_target_rna_corr = self.df_corr.iat[corsiv_idx, surrogate_meth_target_rna_corr_col_idx]
                
                        if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                            sys.stderr.write("TS_corr: Corsiv %s Target %s surrogate %s pair_string %s col_idx %s check column %s corr %s\n"%(current_corsiv, rna_tissue, surrogate_tissue, surrogate_meth_target_rna_corr_string, surrogate_meth_target_rna_corr_col_idx, self.corr_columns[surrogate_meth_target_rna_corr_col_idx], surrogate_meth_target_rna_corr))
                
                        # check surrogate rna/meth correlation
                        if (surrogate_meth_target_rna_corr>self.min_abs_pcc) and (surrogate_meth_target_rna_corr<=1):
                            report_label = self.labels_hash[surrogate_tissue]["positive_corr_target_positive_corr_surrogate_label"]
                            self.tissue_centric_corr_hash[rna_tissue][report_label]+=1
                            venn_support_set.add(surrogate_tissue)
                            
                            if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                                sys.stderr.write("pos_target_pos_surrogate\n")
                        
                        elif (surrogate_meth_target_rna_corr>=-1) and (surrogate_meth_target_rna_corr<-self.min_abs_pcc):
                            report_label = self.labels_hash[surrogate_tissue]["positive_corr_target_negative_corr_surrogate_label"]
                            self.tissue_centric_corr_hash[rna_tissue][report_label]+=1
                            
                            if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                                sys.stderr.write("pos_target_neg_surrogate\n")
                    
                    if (len(venn_support_set)==0):
                        self.venn_hash[pos_corr_label]["Target"]+=1
                    elif (len(venn_support_set)==2):
                        self.venn_hash[pos_corr_label]["Two_Surrogates"]+=1
                    else:
                        surrogate_supporting = list(venn_support_set)[0]
                        self.venn_hash[pos_corr_label][surrogate_supporting]+=1
                            
                            
                # check target rna/meth correlation
                elif (target_meth_rna_corr>=-1) and (target_meth_rna_corr<-self.min_abs_pcc):
                    report_label = self.labels_hash["negative_corr_target_label"]
                    self.tissue_centric_corr_hash[rna_tissue][report_label]+=1
                    
                    venn_support_set = set()
                    
                    for surrogate_tissue in self.surrogate_tissues:
                        surrogate_meth_target_rna_corr_string = "METH_%s_RNA_%s"%(surrogate_tissue, rna_tissue)
                        surrogate_meth_target_rna_corr_col_idx = self.corr_df_col_to_idx_hash[surrogate_meth_target_rna_corr_string]
                        surrogate_meth_target_rna_corr = self.df_corr.iat[corsiv_idx, surrogate_meth_target_rna_corr_col_idx]
                
                        if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                            sys.stderr.write("TS_corr: Corsiv %s Target %s surrogate %s pair_string %s col_idx %s check column %s corr %s\n"%(current_corsiv, rna_tissue, surrogate_tissue, surrogate_meth_target_rna_corr_string, surrogate_meth_target_rna_corr_col_idx, self.corr_columns[surrogate_meth_target_rna_corr_col_idx], surrogate_meth_target_rna_corr))
                
                        # check surrogate rna/meth correlation
                        if (surrogate_meth_target_rna_corr>self.min_abs_pcc) and (surrogate_meth_target_rna_corr<=1):
                            report_label = self.labels_hash[surrogate_tissue]["negative_corr_target_positive_corr_surrogate_label"]
                            self.tissue_centric_corr_hash[rna_tissue][report_label]+=1
                            
                            if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                                sys.stderr.write("neg_target_pos_surrogate\n")
                        
                        elif (surrogate_meth_target_rna_corr>=-1) and (surrogate_meth_target_rna_corr<-self.min_abs_pcc):
                            report_label = self.labels_hash[surrogate_tissue]["negative_corr_target_negative_corr_surrogate_label"]
                            self.tissue_centric_corr_hash[rna_tissue][report_label]+=1
                            venn_support_set.add(surrogate_tissue)
                            
                            if (self.DEBUG_DO_TISSUE_CENTRIC_REPORT):
                                sys.stderr.write("neg_target_neg_surrogate\n")                                                
                        
                    if (len(venn_support_set)==0):
                        self.venn_hash[neg_corr_label]["Target"]+=1
                    elif (len(venn_support_set)==2):
                        self.venn_hash[neg_corr_label]["Two_Surrogates"]+=1
                    else:
                        surrogate_supporting = list(venn_support_set)[0]
                        self.venn_hash[neg_corr_label][surrogate_supporting]+=1
                                            
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::doTargetSurrogateReport\n"%SimpleTime.now())
    
    def outputTargetSurrogateReport(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::outputTargetSurrogateReport\n"%SimpleTime.now())

        target_surrogate_report = "%s.target_surrogate_report.xls"%self.myArgs.outputRoot
        target_surrogate_report_writer = open(target_surrogate_report, "wt")
        
        target_surrogate_report_writer.write("Tissue\t%s\n"%"\t".join(self.report_header))
        
        for rna_tissue in self.rna_labels:
            buffer = [rna_tissue]
            
            for label in self.report_header:
                buffer.append(self.tissue_centric_corr_hash[rna_tissue][label])
            
            target_surrogate_report_writer.write("%s\n"%"\t".join([str(x) for x in buffer]))    
        
        target_surrogate_report_writer.close()
        
        
        target_surrogate_venn_report = "%s.target_surrogate_venn_report.xls"%self.myArgs.outputRoot
        target_surrogate_venn_report_writer = open(target_surrogate_venn_report, "wt")
        target_surrogate_venn_report_writer.write("Tissue\t%s\n"%"\t".join(self.venn_labels))
        
        for rna_tissue in self.rna_labels:
            pos_corr_label = "PosCorr_%s"%rna_tissue
            neg_corr_label = "NegCorr_%s"%rna_tissue
            
            buffer_pos = [pos_corr_label]
            buffer_neg = [neg_corr_label]
                
            
            for label in self.venn_labels:
                buffer_pos.append(self.venn_hash[pos_corr_label][label])
                buffer_neg.append(self.venn_hash[neg_corr_label][label])
            
            target_surrogate_venn_report_writer.write("%s\n"%"\t".join([str(x) for x in buffer_pos]))    
            target_surrogate_venn_report_writer.write("%s\n"%"\t".join([str(x) for x in buffer_neg]))    
        
        target_surrogate_venn_report_writer.close()
                 
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::outputTargetSurrogateReport\n"%SimpleTime.now())
            
    def tissueCentricCounts(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::tissueCentricCounts\n"%SimpleTime.now())
    
        self.setupTargetSurrogateCalc()    
        self.doTargetSurrogateReport()
        
        self.outputTargetSurrogateReport()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::tissueCentricCounts\n"%SimpleTime.now())
    
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START SummaryCorsivBetaGeneExpression::work\n"%SimpleTime.now())
        
        self.min_abs_pcc = math.sqrt(float(self.myArgs.minimumR2))
        sys.stderr.write("Minimum abs PCC %s\n"%self.min_abs_pcc)
    
        self.loadCorrelationMatrix()
        # self.summarizeCounts()
        # self.outputSummary()
        self.tissueCentricCounts()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP SummaryCorsivBetaGeneExpression::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = SummaryCorsivBetaGeneExpression.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = SummaryCorsivBetaGeneExpression(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
