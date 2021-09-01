#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "8.0"

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
        
class AssessCaptureCORSiV:
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
    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.
            
annotate capture samples by individual
select only individuals with at least 2,3,4,5 ... max possible tissues
TODO: pool of processes
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-x','--corsivMethylationMatrix',   help='data matrix with beta values for CORSiVs',    required=True)
        parser.add_argument('-S','--tissueMap',                 help='map for samples/tissues/capture files',       required=True)
        parser.add_argument('-m','--maxWIR',                    help='max within individual range threshold',       required=True)
        parser.add_argument('-d','--minDonorsITC',              help='minimum donors for ITC',                      required=True)
        parser.add_argument('-C','--plotHeatmap',               help='[flag] plot heatmap (default false)',       required=False, action="store_true")
        parser.add_argument('-P','--plotITCCoRSIVList',         help='[optional] provide a file with a list of CoRSIVs for which to generate ITC scatterplots',  required=False)
        parser.add_argument('-o','--outputRoot',                help='output file root',                            required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    
    def loadTissueMap(self):    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::loadTissueMap\n"%SimpleTime.now())
        
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
            sys.stderr.write("[%s] START AssessCaptureCORSiV::loadTissueMap\n"%SimpleTime.now())


    def cleanupColumns(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::cleanupColumns\n"%SimpleTime.now())
        
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
            sys.stderr.write("[%s] STOP AssessCaptureCORSiV::cleanupColumns\n"%SimpleTime.now())
            
        return cleanup_column_list  
       
    def findTissueCoverage(self, column_list):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::findTissueCoverage\n"%SimpleTime.now())
            
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
            sys.stderr.write("[%s] STOP AssessCaptureCORSiV::findTissueCoverage\n"%SimpleTime.now())
            
        return tissue_count_hash, max_coverage, tissue_to_donor_hash
     
    def loadCoverage(self):    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::loadCoverage\n"%SimpleTime.now())
            
        self.coverage_df = pd.read_csv(self.myArgs.corsivMethylationMatrix, sep="\t", index_col=0, header=0)
        
        cleanup_column_list = self.cleanupColumns()
        self.coverage_df.columns = cleanup_column_list
        
        self.multi_tissue_hash, self.max_coverage, self.tissue_to_donor_hash = self.findTissueCoverage(cleanup_column_list)
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::loadCoverage\n"%SimpleTime.now())
    
    def mapColorsToDonors(self, number_of_donors):
        actual_donors = min(50, number_of_donors)
        
        ret = []
        r = int(random.random() * 256)
        g = int(random.random() * 256)
        b = int(random.random() * 256)
        step = 256 / actual_donors
        for i in range(actual_donors):
          r += step
          g += step
          b += step
          r = int(r) % 256
          g = int(g) % 256
          b = int(b) % 256
          ret.append((r,g,b))
          
        full_colors = []
        ret_idx = -1
        for idx in range(number_of_donors):
            ret_idx = (ret_idx +1)%actual_donors
              
    def generateMultiTissueHeatmaps(self):    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::generateMultiTissueHeatmaps\n"%SimpleTime.now())
        
        # identify individuals with multiple tissues
        self.min_coverage_to_tissue_hash = {}
        
        number_of_donors = len(self.multi_tissue_hash)
        donor_color_map = self.mapColorsToDonors(number_of_donors)
        tissue_color_hash = {}
        
        for min_tissue_coverage in range(5, self.max_coverage+1):
            column_ids_set = set()
            for gtex_donor_id in self.multi_tissue_hash:
                donor_tissue_count = len(self.multi_tissue_hash[gtex_donor_id])
                if (donor_tissue_count>=min_tissue_coverage):
                    column_ids_set = column_ids_set.union(self.multi_tissue_hash[gtex_donor_id])
                    
                    if (self.DEBUG_MULTI_TISSUE_HEATMAPS):
                        sys.stderr.write("Donor %s with %s tissues >= %s --> tally %s tissues\n"%(gtex_donor_id, donor_tissue_count, min_tissue_coverage, len(column_ids_set)))
                    
            self.min_coverage_to_tissue_hash[min_tissue_coverage] = column_ids_set
            # TODO: here select samples and perfom heatmap
            slice_df = self.coverage_df.iloc[:,list(column_ids_set)]
            slice_data_output_tab = "%s.at_least_%s_tissues.xls"%(self.myArgs.outputRoot, min_tissue_coverage)
            
            slice_df.to_csv(slice_data_output_tab, sep="\t", index_label="CORSiVs")
            
            slice_data_output_pdf = "%s.at_least_%s_tissues.pdf"%(self.myArgs.outputRoot, min_tissue_coverage)
            slice_data_output_jpg = "%s.at_least_%s_tissues.jpg"%(self.myArgs.outputRoot, min_tissue_coverage)
            
            plt.figure()
            sn.set()
            sn.clustermap(slice_df,  cmap='bwr', col_colors=None, row_cluster=True, col_cluster=True, metric="euclidean", linewidths=0.0, rasterized=True, figsize=(int(float(slice_df.shape[1])*0.3),30))
            plt.savefig(slice_data_output_jpg, dpi=300)
            plt.savefig(slice_data_output_pdf, dpi=300)
            
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::generateMultiTissueHeatmaps\n"%SimpleTime.now())
    
    def analyzeWir(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::analyzeWir; max_WIR = %s\n"%(SimpleTime.now(), self.myArgs.maxWIR))


        # set up CORSiv maxWIR hash
        max_wir_threshold = float(self.myArgs.maxWIR)
        self.max_WIR_hash = {}
        self.passing_max_WIR_list = []
        
        self.corsiv_name_list = list(self.coverage_df.index)

        max_WIR_file = "%s.max_WIR.xls"%self.myArgs.outputRoot
        max_WIR_file_writer = open(max_WIR_file, "wt")
        max_WIR_file_writer.write("CORSiV\tIndividuals\tMax WIR\n")
        
        for corsiv_idx in range(len(self.coverage_df)):
            corsiv_name = self.corsiv_name_list[corsiv_idx]
            if (self.DEBUG_MAX_WIR_SHORT):
                sys.stderr.write("Analyzing maxWIR in corsiv %s\n"%corsiv_name)

            wir_list = []
            
            for gtex_donor_id in self.multi_tissue_hash:
                col_idx_list = list(self.multi_tissue_hash[gtex_donor_id])
                if (self.DEBUG_MAX_WIR_VERBOSE):
                    sys.stderr.write("Corsiv %s donor %s tissue count %s\n"%(corsiv_name, gtex_donor_id, len(col_idx_list)))
                    
                if (len(col_idx_list)>1):
                    beta_values_within_individual = self.coverage_df.iloc[corsiv_idx, col_idx_list]
                    min_beta = np.min(beta_values_within_individual)
                    max_beta = np.max(beta_values_within_individual)
                    
                    wir = max_beta - min_beta
                    if (self.DEBUG_MAX_WIR_VERBOSE):
                        sys.stderr.write("Corsiv %s donor %s tissue count %s min_beta %s max_beta %s WIR: %s\n"%(corsiv_name, gtex_donor_id, len(col_idx_list), min_beta, max_beta, wir))
                        
                    wir_list.append(wir)
            
            wir_float_array = np.asfarray(wir_list)
            max_WIR = np.max(wir_float_array)
            self.max_WIR_hash[corsiv_name]=max_WIR
            max_WIR_file_writer.write("%s\t%s\t%s\n"%(corsiv_name, len(col_idx_list), max_WIR))
            if (max_WIR<=max_wir_threshold):
                self.passing_max_WIR_list.append(corsiv_idx)
            if (self.DEBUG_MAX_WIR_SHORT):
                sys.stderr.write("Corsiv %s ge2 individuals %s maxWIR: %s\n"%(corsiv_name, len(wir_list), max_WIR))


        max_WIR_file_writer.close()
        
       
        # further subset corsivs based on wir 
        self.coverage_df = self.coverage_df.iloc[self.passing_max_WIR_list, :]
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessCaptureCORSiV::analyzeWir\n"%SimpleTime.now())
    
    def corrPairPlotRhoDistribution(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::corrPairPlotRhoDistribution\n"%SimpleTime.now())

        tissue_pair_count_hash = {}
        self.min_donors_itc = int(self.myArgs.minDonorsITC)
        if (self.min_donors_itc<10):
            sys.stderr.write("Need at least 10 donors to compute correlations !\nExiting ...\n")
            sys.exit(4)
        
        self.number_of_tissues = len(self.tissue_to_donor_hash)
        if (self.DEBUG_TISSUE_PAIR_RHO):
            sys.stderr.write("Have overall %s tissues\n"%self.number_of_tissues)
        
        coverage_index = list(self.coverage_df.index)
        
        self.corsiv_set = set(coverage_index)
        
        self.itc_pair_hash = {}
        self.wir_pair_hash = {}
        self.itc_pair_hash_per_corsiv = {}
        
        # set up subplots grid
        self.corsiv_range_hash = {}
        
        tissue_list = sorted( list(self.tissue_to_donor_hash.keys()))
         
        corr_distribution_pair_plot_pdf = "%s.pair_corr.pdf"%self.myArgs.outputRoot  
        plt.figure()
        fig, ax = plt.subplots(self.number_of_tissues, self.number_of_tissues, figsize=(20,20))
    
        for tissue_idx_i in range(self.number_of_tissues):
            gtex_tissue_id_i = tissue_list[tissue_idx_i]
            tissue_i_donors = set(self.tissue_to_donor_hash[gtex_tissue_id_i].keys())
            if (self.DEBUG_TISSUE_PAIR_RHO):
                sys.stderr.write("Processing tissue i %s with %s donors\t"%(gtex_tissue_id_i, len(tissue_i_donors)))
            
            # compute per tissue range
            tissue_i_donors_list = list(tissue_i_donors)
            
            for corsiv_idx in range(len(coverage_index)):
                corsiv_name = coverage_index[corsiv_idx]
                if (self.DEBUG_CORSIV_RANGE):
                    sys.stderr.write("Evaluating range for tissue %s corsiv %s\n"%(gtex_tissue_id_i, corsiv_name))

                if not (corsiv_name in self.corsiv_range_hash):
                    corsiv_range_info = CStruct(corsiv_name=corsiv_name, min_beta=2, max_beta=-1, range_all = "NA", tissue_range_hash={})
                    self.corsiv_range_hash[corsiv_name]=corsiv_range_info
                   
                corsiv_range_info =  self.corsiv_range_hash[corsiv_name]
                
                tissue_i_donors_col_idx = []
                    
                for gtex_donor_id in tissue_i_donors:
                    tissue_i_col_idx = self.tissue_to_donor_hash[gtex_tissue_id_i][gtex_donor_id]
                    i_nan_flag = np.isnan(self.coverage_df.iat[corsiv_idx, tissue_i_col_idx])
                    if not i_nan_flag:
                        tissue_i_donors_col_idx.append(tissue_i_col_idx)
                
                number_of_i_donors = len(tissue_i_donors_col_idx)
                if (number_of_i_donors<2):
                    if (self.DEBUG_CORSIV_RANGE):
                        sys.stderr.write("FailRange: tissue %s corsiv %s donors=%s<2\n"%(gtex_tissue_id_i, corsiv_name, number_of_i_donors))
                else:
                    tissue_i_corsiv_beta_array = self.coverage_df.iloc[corsiv_idx, tissue_i_donors_col_idx]
                    min_beta = np.min(tissue_i_corsiv_beta_array)
                    max_beta = np.max(tissue_i_corsiv_beta_array)
                    if (self.DEBUG_CORSIV_RANGE):
                        sys.stderr.write("PassRange1: Corsiv %s Tissue %s corsiv.min_beta=%s corsiv.max_beta=%s min_beta=%s max_beta=%s\n"%(corsiv_name, gtex_tissue_id_i, corsiv_range_info.min_beta, corsiv_range_info.max_beta, min_beta, max_beta))
                    corsiv_range_info.min_beta = min(corsiv_range_info.min_beta, min_beta)
                    corsiv_range_info.max_beta = max(corsiv_range_info.max_beta, max_beta)
                    corsiv_x_tissue_range = max_beta - min_beta
                    corsiv_range_info.tissue_range_hash[gtex_tissue_id_i]=corsiv_x_tissue_range
                    
                    if (self.DEBUG_CORSIV_RANGE):
                        sys.stderr.write("PassRange2: Corsiv %s Tissue %s corsiv.min_beta=%s corsiv.max_beta=%s min_beta=%s max_beta=%s range %s =  %s\n"%(corsiv_name, gtex_tissue_id_i, corsiv_range_info.min_beta, corsiv_range_info.max_beta, min_beta, max_beta, corsiv_x_tissue_range, corsiv_range_info.tissue_range_hash[gtex_tissue_id_i]))
                                          
            for tissue_idx_j in range(tissue_idx_i+1, self.number_of_tissues):
                gtex_tissue_id_j = tissue_list[tissue_idx_j]
                tissue_j_donors = set(self.tissue_to_donor_hash[gtex_tissue_id_j].keys())
                if (self.DEBUG_TISSUE_PAIR_RHO):
                    sys.stderr.write("Processing tissue j %s with %s donors\t"%(gtex_tissue_id_j, len(tissue_j_donors)))
 
                tissue_pair = "%s==>%s"%(gtex_tissue_id_i, gtex_tissue_id_j)
                if (gtex_tissue_id_i>gtex_tissue_id_j):
                    tissue_pair = "%s==>%s"%(gtex_tissue_id_j, gtex_tissue_id_i)
                    
                common_donors = tissue_i_donors.intersection(tissue_j_donors)
                number_of_common_donors = len(common_donors)
                tissue_pair_count_hash[tissue_pair] = number_of_common_donors

                if (self.DEBUG_TISSUE_PAIR_RHO):
                    sys.stderr.write("Intersection ! [%s = %s] x [%s = %s] ---> %s donors\n"%(gtex_tissue_id_i, len(tissue_i_donors), gtex_tissue_id_j, len(tissue_j_donors), len(common_donors)))
                
                                          
            for tissue_idx_j in range(tissue_idx_i+1, self.number_of_tissues):
                gtex_tissue_id_j = tissue_list[tissue_idx_j]
                tissue_j_donors = set(self.tissue_to_donor_hash[gtex_tissue_id_j].keys())
                if (self.DEBUG_TISSUE_PAIR_RHO):
                    sys.stderr.write("Processing tissue j %s with %s donors\t"%(gtex_tissue_id_j, len(tissue_j_donors)))
 
                tissue_pair = "%s==>%s"%(gtex_tissue_id_i, gtex_tissue_id_j)
                if (gtex_tissue_id_i>gtex_tissue_id_j):
                    tissue_pair = "%s==>%s"%(gtex_tissue_id_j, gtex_tissue_id_i)
                    
                common_donors = tissue_i_donors.intersection(tissue_j_donors)
                number_of_common_donors = len(common_donors)
                tissue_pair_count_hash[tissue_pair] = number_of_common_donors

                if (self.DEBUG_TISSUE_PAIR_RHO):
                    sys.stderr.write("Intersection ! [%s = %s] x [%s = %s] ---> %s donors\n"%(gtex_tissue_id_i, len(tissue_i_donors), gtex_tissue_id_j, len(tissue_j_donors), len(common_donors)))
                
                if (number_of_common_donors>=self.min_donors_itc):
                    tissue_i_common_donors_col_idx = []
                    tissue_j_common_donors_col_idx = []
                    
                    for gtex_donor_id in common_donors:
                        tissue_i_col_idx = self.tissue_to_donor_hash[gtex_tissue_id_i][gtex_donor_id]
                        tissue_j_col_idx = self.tissue_to_donor_hash[gtex_tissue_id_j][gtex_donor_id]
                        tissue_i_common_donors_col_idx.append(tissue_i_col_idx)
                        tissue_j_common_donors_col_idx.append(tissue_j_col_idx)
                        

                    inter_tissue_corr_list = []
                    inter_tissue_range_list = []
                    
                    self.itc_pair_hash_per_corsiv[tissue_pair] = {}
                    
                    for corsiv_idx in range(len(coverage_index)):
                        corsiv_name = coverage_index[corsiv_idx]
                            
                        if (self.DEBUG_TISSUE_PAIR_RHO_VERBOSE):
                            sys.stderr.write("Checking for NAs tissue pair [%s]x[%s]: corsiv [%s] %s\n"%(gtex_tissue_id_i, gtex_tissue_id_j, corsiv_idx, corsiv_name))
                        
                        final_tissue_i_common_donors_col_idx = []
                        final_tissue_j_common_donors_col_idx = []
                        
                        for running_common_donor_idx in range(len(tissue_i_common_donors_col_idx)):
                            i_col_idx = tissue_i_common_donors_col_idx[running_common_donor_idx]
                            j_col_idx = tissue_j_common_donors_col_idx[running_common_donor_idx]
                            
                            value_tissue_i = self.coverage_df.iat[corsiv_idx, i_col_idx]
                            value_tissue_j = self.coverage_df.iat[corsiv_idx, j_col_idx]
                            i_NaN_flag = np.isnan(value_tissue_i)
                            j_NaN_flag = np.isnan(value_tissue_j)
                            
                            if (self.DEBUG_TISSUE_PAIR_RHO_SUPER_VERBOSE):
                                sys.stderr.write("value pair: %s naflag %s and  %s naflag %s\n"%(value_tissue_i,  i_NaN_flag, value_tissue_j, j_NaN_flag))
                                
                            if (not i_NaN_flag) and (not j_NaN_flag):
                                final_tissue_i_common_donors_col_idx.append(i_col_idx)
                                final_tissue_j_common_donors_col_idx.append(j_col_idx)
                                if (self.DEBUG_TISSUE_PAIR_RHO_SUPER_VERBOSE):
                                    sys.stderr.write("Bump common donors to %s\n"%len(final_tissue_i_common_donors_col_idx))
                            else:
                                if (self.DEBUG_TISSUE_PAIR_RHO_SUPER_VERBOSE):
                                    sys.stderr.write("NoBump common donors stays %s\n"%len(final_tissue_i_common_donors_col_idx))
                        
                        if len(final_tissue_i_common_donors_col_idx)>=self.min_donors_itc:    
                            # check at coverage level every indidividual
                            tissue_i_beta_array = self.coverage_df.iloc[corsiv_idx, final_tissue_i_common_donors_col_idx]
                            tissue_j_beta_array = self.coverage_df.iloc[corsiv_idx, final_tissue_j_common_donors_col_idx]
                            rho_value, rho_pvalue = scipy.stats.pearsonr(tissue_i_beta_array, tissue_j_beta_array)
                            
                            for range_idx in range(len(final_tissue_i_common_donors_col_idx)):
                                inter_tissue_range_list.append(abs(tissue_i_beta_array[range_idx]-tissue_j_beta_array[range_idx]))
                                
                            if (self.DEBUG_TISSUE_PAIR_RHO_VERBOSE):
                                sys.stderr.write("Tissues %s x %s: corsiv [%s]: %s --> rho %.4f\n"%(gtex_tissue_id_i, gtex_tissue_id_j, corsiv_idx, corsiv_name, rho_value))
                            
                            inter_tissue_corr_list.append(rho_value)
                            
                            self.itc_pair_hash_per_corsiv[tissue_pair][corsiv_name]=(tissue_i_beta_array, tissue_j_beta_array, rho_value, len(final_tissue_i_common_donors_col_idx))
                        else:
                            if (self.DEBUG_TISSUE_PAIR_RHO_VERBOSE):
                                sys.stderr.write("FailPair: Tissues %s x %s: corsiv [%s]: %s --> only %s common donors\n"%(gtex_tissue_id_i, gtex_tissue_id_j, corsiv_idx, corsiv_name, len(final_tissue_i_common_donors_col_idx)))
                                                       
                # plot the damn thing
                ax[tissue_idx_i, tissue_idx_j].hist(inter_tissue_corr_list, bins=20)
                ax[tissue_idx_i, tissue_idx_j].set_xlim(-1,1)
                
                self.itc_pair_hash[tissue_pair] = inter_tissue_corr_list
                
                ax[tissue_idx_j, tissue_idx_i].hist(inter_tissue_range_list, bins=20)
                ax[tissue_idx_j, tissue_idx_i].set_xlim(0,1)
                
                self.wir_pair_hash[tissue_pair] = inter_tissue_range_list


        # compute overall range per corsiv
        for corsiv_name in  self.corsiv_set:
            corsiv_range_info = self.corsiv_range_hash[corsiv_name]
            if (corsiv_range_info.max_beta>=0) and (corsiv_range_info.min_beta<=1):
                corsiv_range_info.range_all = corsiv_range_info.max_beta - corsiv_range_info.min_beta
                    
            if (self.DEBUG_CORSIV_RANGE):
                sys.stderr.write("RangeAll: Corsiv %s min_beta %s max_beta %s range %s\n"%(corsiv_name, corsiv_range_info.min_beta, corsiv_range_info.max_beta, corsiv_range_info.range_all))
        
        for tissue_idx_j in range(self.number_of_tissues):
            gtex_tissue_id_j = tissue_list[tissue_idx_j]
            ax[self.number_of_tissues-1, tissue_idx_j].set_xlabel(gtex_tissue_id_j)
            
        for tissue_idx_i in range(self.number_of_tissues):
            gtex_tissue_id_i = tissue_list[tissue_idx_i]
            ax[tissue_idx_i, 0].set_ylabel(gtex_tissue_id_i)
                    
            
            
        plt.savefig(corr_distribution_pair_plot_pdf, dpi=300)
        
        # TODO: report also number of donors with both tissues, and conditionally color in Excel
        pair_tissue_count_file = "%s.pair_count.xls"%self.myArgs.outputRoot
        pair_tissue_count_file_writer = open(pair_tissue_count_file, "wt")
        pair_tissue_count_file_writer.write("Tissue Pair\t%s\n"%"\t".join(tissue_list))
        
        for tissue_idx_i in range(self.number_of_tissues):
            gtex_tissue_id_i = tissue_list[tissue_idx_i]
            buffer = [""]*(1+tissue_idx_i)
            for tissue_idx_j in range(tissue_idx_i+1, self.number_of_tissues):
                gtex_tissue_id_j = tissue_list[tissue_idx_j]
                tissue_pair = "%s==>%s"%(gtex_tissue_id_i, gtex_tissue_id_j)
                if (gtex_tissue_id_i>gtex_tissue_id_j):
                    tissue_pair = "%s==>%s"%(gtex_tissue_id_j, gtex_tissue_id_i)
                overlap_size = tissue_pair_count_hash[tissue_pair]
                if (self.DEBUG_TISSUE_PAIR_RHO):
                    sys.stderr.write("x [%s, %s] vs [%s, %s] --> %s common donors\n"%(tissue_idx_i, gtex_tissue_id_i, tissue_idx_j, gtex_tissue_id_j, overlap_size))
                buffer.append(str(overlap_size))
            pair_tissue_count_file_writer.write("%s\t%s\n"%(gtex_tissue_id_i, "\t".join(buffer)))
                
        pair_tissue_count_file_writer.close()
        
        
        # finally plot as long boxplot the itc and the wir
        self.boxPlotSeries(self.itc_pair_hash, "%s.itc_tissue_pairs"%self.myArgs.outputRoot, "ITC over tissue pairs")        
        # self.boxPlotSeries(self.wir_pair_hash, "%s.wir_tissue_pairs"%self.myArgs.outputRoot, "Within Individual Range over tissue pairs")
        
        self.outputCorsivTissueCorrelationMap()
        
   #     self.reportPerCoRSIVCorrelations()
        
        if (self.myArgs.plotITCCoRSIVList):
            self.plotPerCoRSIVScatterPlots()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessCaptureCORSiV::corrPairPlotRhoDistribution\n"%SimpleTime.now())
       
   
    def outputCorsivTissueCorrelationMap(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::outputCorsivTissueCorrelationMap\n"%SimpleTime.now())
            
        full_corsiv_tissuepair_correlation_file = "%s.corsiv_tissuepair_correlation.xls"%self.myArgs.outputRoot
        full_corsiv_tissuepair_correlation_file_writer = open(full_corsiv_tissuepair_correlation_file, "wt")
        
        full_corsiv_tissuepair_correlation_binary_file = "%s.corsiv_tissuepair_correlation_binary.xls"%self.myArgs.outputRoot
        full_corsiv_tissuepair_correlation_binary_file_writer = open(full_corsiv_tissuepair_correlation_binary_file, "wt") 
        
        
        tissue_pair_list = self.itc_pair_hash_per_corsiv.keys()
        full_corsiv_tissuepair_correlation_binary_file_writer.write("CoRSIV\t%s\n"%"\t".join(tissue_pair_list))
        full_corsiv_tissuepair_correlation_file_writer.write("CoRSIV\t%s\n"%"\t".join(tissue_pair_list))

        for corsiv_name in  self.corsiv_set:
            buffer_full = [corsiv_name]
            buffer_binary = [corsiv_name]
            
            for tissue_pair in tissue_pair_list:
                if corsiv_name in self.itc_pair_hash_per_corsiv[tissue_pair]:
                    _, _, rho, _ = self.itc_pair_hash_per_corsiv[tissue_pair][corsiv_name]
                else:
                    rho=0
                buffer_full.append(rho)
                if (rho>=math.sqrt(0.5)):
                    buffer_binary.append(1)
                else:
                    buffer_binary.append(0)
            full_corsiv_tissuepair_correlation_file_writer.write("%s\n"%"\t".join([str(x) for x in buffer_full]))
            full_corsiv_tissuepair_correlation_binary_file_writer.write("%s\n"%"\t".join([str(x) for x in buffer_binary]))
                    
        full_corsiv_tissuepair_correlation_binary_file_writer.close()
        full_corsiv_tissuepair_correlation_file_writer.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessCaptureCORSiV::outputCorsivTissueCorrelationMap\n"%SimpleTime.now())
            
    def reportPerCoRSIVCorrelations(self):

        corsiv_itc_report = "%s.corsiv_itc_per_tissue_pair.xls"%self.myArgs.outputRoot
        corsiv_itc_report_writer =  open(corsiv_itc_report, "wt")
        
        tissue_list = sorted( list(self.tissue_to_donor_hash.keys()))   
        
        tissue_pairs  = sorted(list(self.itc_pair_hash_per_corsiv.keys()))
        
        header = ["CoRSIVs", "Tissues_covered", "min_ITC"]
        for pair in tissue_pairs:
            header.append("%s rho"%pair)
            header.append("%s donors"%pair)
        
        header.append("Range")
        header.append("min_BetaValue")
        header.append("max_BetaValue")
        
        for tissue in tissue_list:
            header.append("Range %s"%tissue)
            
        corsiv_itc_report_writer.write("%s\n"%"\t".join(header))
            
        for corsiv_name in  sorted(list(self.corsiv_set)):
            if (self.DEBUG_TISSUE_PAIR_CORSIV_REPORT):
                sys.stderr.write("Generating  tissue pair report line for for %s\n"%corsiv_name)
            
            min_ITC = 2
            
            buffer = [corsiv_name, "-1", "NA"]
            for tissue_pair in tissue_pairs:
                if (corsiv_name in self.itc_pair_hash_per_corsiv[tissue_pair]):
                    i_beta_values, j_beta_values, rho_value, donor_count = self.itc_pair_hash_per_corsiv[tissue_pair][corsiv_name]
                    buffer.append(str(rho_value))        
                    buffer.append(str(donor_count))
                    if (min_ITC>rho_value):
                        min_ITC = rho_value
                else:
                    buffer.append("NA")
                    buffer.append("NA")
            if (min_ITC<1):
                buffer[2]=str(min_ITC)
                
            if (corsiv_name in self.corsiv_range_hash):
                corsiv_range_info = self.corsiv_range_hash[corsiv_name]
                buffer.append(str(corsiv_range_info.range_all))
                buffer.append(str(corsiv_range_info.min_beta))
                buffer.append(str(corsiv_range_info.max_beta))
                for tissue in tissue_list:
                    if (tissue in corsiv_range_info.tissue_range_hash):
                        buffer.append(str(corsiv_range_info.tissue_range_hash[tissue]))
                    else:
                        buffer.append("NA")
            else:
                buffer.append("NA")
                buffer.append("NA")
                buffer.append("NA")
                for tissue in tissue_list:
                    buffer.append("NA")
            
            corsiv_itc_report_writer.write("%s\n"%"\t".join(buffer))
                            
        corsiv_itc_report_writer.close()
                    
    def plotPerCoRSIVScatterPlots(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::plotPerCoRSIVScatterPlots\n"%SimpleTime.now())
            
        os.system("mkdir -p %s.corsiv_itc"%self.myArgs.outputRoot)
        
        # set up subplots grid
        tissue_list = sorted( list(self.tissue_to_donor_hash.keys()))
        
        
        corsiv_itc_plot_list = []
        corsiv_list_reader = open(self.myArgs.plotITCCoRSIVList, "rt")
        for line in corsiv_list_reader:
            corsiv_itc_plot_list.append(line.strip())
        corsiv_list_reader.close()
            
        plot_idx = 0
        for corsiv_name in  corsiv_itc_plot_list:
            if (self.DEBUG_TISSUE_PAIR_SCATTERPLOT):
                sys.stderr.write("Generating scatterplot for %s\n"%corsiv_name)
                
            cleaned_up_corsiv_name = corsiv_name.replace("N/A", "NA").replace(";+","").replace(" ","_").replace(";","_")
                
            corr_distribution_pair_plot_pdf = "%s.corsiv_itc/%s.pair_corr.pdf"%(self.myArgs.outputRoot, cleaned_up_corsiv_name)
#            corr_distribution_pair_plot_jpg = "%s.corsiv_itc/%s.pair_corr.jpg"%(self.myArgs.outputRoot, cleaned_up_corsiv_name)
            
            plt.figure()
            fig, ax = plt.subplots(self.number_of_tissues-1, self.number_of_tissues-1, figsize=(20,20))
            
            for tissue_idx_i in range(0, self.number_of_tissues-1):
                for tissue_idx_j in range(tissue_idx_i+1, self.number_of_tissues-1):
                    ax[tissue_idx_i, tissue_idx_j].axis('off')
            
            for tissue_idx_i in range(self.number_of_tissues):
                gtex_tissue_id_i = tissue_list[tissue_idx_i]
                if (self.DEBUG_TISSUE_PAIR_SCATTERPLOT):
                    sys.stderr.write("Tissue i %s\t"%gtex_tissue_id_i)
                                      
                for tissue_idx_j in range(tissue_idx_i+1, self.number_of_tissues):
                    gtex_tissue_id_j = tissue_list[tissue_idx_j]
                    if (self.DEBUG_TISSUE_PAIR_RHO):
                        sys.stderr.write("Tissue j %s\n"%gtex_tissue_id_j)
     
                    tissue_pair = "%s==>%s"%(gtex_tissue_id_i, gtex_tissue_id_j)
                    if (gtex_tissue_id_i>gtex_tissue_id_j):
                        tissue_pair = "%s==>%s"%(gtex_tissue_id_j, gtex_tissue_id_i)
                    
                    if (tissue_pair in self.itc_pair_hash_per_corsiv):
                        if (self.DEBUG_TISSUE_PAIR_RHO):
                            sys.stderr.write("OK: %s x %s --> %s\n"%(gtex_tissue_id_i, gtex_tissue_id_j, corsiv_name))
                            
                        if (corsiv_name in self.itc_pair_hash_per_corsiv[tissue_pair]):
                            i_beta_values, j_beta_values, rho_value, donor_count = self.itc_pair_hash_per_corsiv[tissue_pair][corsiv_name]
                            ax[tissue_idx_j-1, tissue_idx_i].scatter(i_beta_values, j_beta_values, color='b')
                            b, m = polyfit(i_beta_values, j_beta_values, 1)
                            
                            ax[tissue_idx_j-1, tissue_idx_i].plot(i_beta_values, b+m*i_beta_values, '-', color='k')
                            ax[tissue_idx_j-1, tissue_idx_i].set_xlim(0,1)
                            ax[tissue_idx_j-1, tissue_idx_i].set_ylim(0,1)
                            ax[tissue_idx_j-1, tissue_idx_i].text(0.2,0.2, "PCC R=%.4f N=%s"%(rho_value, donor_count))
                            
                        #    ax[tissue_idx_i, tissue_idx_j-1].plot(j_beta_values, b+m*j_beta_values, '-')
                        #    ax[tissue_idx_i, tissue_idx_j-1].set_xlim(0,1)
                        #    ax[tissue_idx_i, tissue_idx_j-1].set_ylim(0,1)
                        #    ax[tissue_idx_i, tissue_idx_j-1].text(0.2,0.2, "PCC R=%.4f N=%s"%(rho_value, donor_count))
                                        
                        
            
            for tissue_idx_j in range(self.number_of_tissues-1):
                gtex_tissue_id_j = tissue_list[tissue_idx_j+1]
                ax[tissue_idx_j, 0].set_ylabel(gtex_tissue_id_j)
#                ax[self.number_of_tissues-2, tissue_idx_j].set_xlabel(gtex_tissue_id_j)
            
            for tissue_idx_i in range(self.number_of_tissues-1):
                gtex_tissue_id_i = tissue_list[tissue_idx_i]
                ax[self.number_of_tissues-2, tissue_idx_i].set_xlabel(gtex_tissue_id_i)
#                ax[tissue_idx_i, 0].set_ylabel(gtex_tissue_id_i)
                
            plt.savefig(corr_distribution_pair_plot_pdf, dpi=300)
            
            plot_idx += 1
        
    def boxPlotSeries(self, series_hash, output_file_root, boxplot_title):
        swarmplot_pdf_file = "%s.swarmplot.pdf"%output_file_root
        swarmplot_jpg_file = "%s.swarmplot.jpg"%output_file_root

        list_of_lists = []
        for series_key in series_hash:
            for value in series_hash[series_key]:
                list_of_lists.append([series_key, value])
        df_for_swarmplot = pd.DataFrame(list_of_lists, columns = ["ITC_pair", "PCC"]) 
        
        plt.figure()
        sn.set(style='white')
        ax = sn.swarmplot(x="ITC_pair", y="PCC", data=df_for_swarmplot, size=1)
        ax.set_ylim(-1, 1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.title(boxplot_title)

        plt.savefig(swarmplot_pdf_file, bbox_inches='tight',dpi=300)
        plt.savefig(swarmplot_jpg_file, bbox_inches='tight',dpi=300)

        boxplot_pdf_file = "%s.pdf"%output_file_root
        boxplot_jpg_file = "%s.jpg"%output_file_root
    
        violin_plot_pdf = "%s.violin_plot.pdf"%output_file_root
        violin_plot_jpg = "%s.violin_plot.jpg"%output_file_root
        plt.figure()
        sn.set(style="white")
        ax = sn.violinplot(x="ITC_pair", y="PCC", data=df_for_swarmplot, linewidth=0.5, cut=0)
#        ax = sn.violinplot(x="ITC_pair", y="PCC", data=df_for_swarmplot, width=1.2, linewidth=0.5, cut=0)
        ax.set_ylim(-1,1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=300)
        plt.savefig(violin_plot_jpg, bbox_inches='tight', dpi=300)
        
        fig, ax = plt.subplots(1)
        box_data = []
        labels = []
        for series_key in series_hash:
            box_data.append(series_hash[series_key])
            labels.append(series_key)
    
        box = plt.boxplot(box_data,patch_artist=True,whis="range")
        ax.set_ylim(-1, 1)
        plt.title(boxplot_title)
        plt.xticks(range(1,len(box_data)+1),labels,rotation=90)
    
        plt.savefig(boxplot_pdf_file, bbox_inches='tight',dpi=300)
        plt.savefig(boxplot_jpg_file, bbox_inches='tight',dpi=300)
        
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::work\n"%SimpleTime.now())
    
            
        self.loadTissueMap()
        self.loadCoverage()
        self.analyzeWir()
        if (self.myArgs.plotHeatmap): 
            self.generateMultiTissueHeatmaps()
            
        self.corrPairPlotRhoDistribution()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessCaptureCORSiV::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = AssessCaptureCORSiV.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = AssessCaptureCORSiV(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
