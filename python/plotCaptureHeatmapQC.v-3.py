#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "2.0"

import os, sys, argparse, re
import glob
import datetime
import math
import scipy.stats
from argparse import RawTextHelpFormatter
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
    DEBUG_TISSUE_PAIR_SCATTERPLOT   = False
    
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
    
    def generateMultiTissueHeatmaps(self):    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::generateMultiTissueHeatmaps\n"%SimpleTime.now())
        
        # identify individuals with multiple tissues
        self.min_coverage_to_tissue_hash = {}
        
        for min_tissue_coverage in range(3, self.max_coverage+1):
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
            sn.clustermap(slice_df,  cmap='bwr', col_colors=None, row_cluster=True, col_cluster=True, metric="euclidean", linewidths=0.0, rasterized=True, figsize=(int(float(slice_df.shape[1])*0.5),30))
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
        
        self.number_of_tissues = len(self.tissue_to_donor_hash)
        if (self.DEBUG_TISSUE_PAIR_RHO):
            sys.stderr.write("Have overall %s tissues\n"%self.number_of_tissues)
        
        coverage_index = list(self.coverage_df.index)
        
        self.corsiv_set = set(coverage_index)
        
        self.itc_pair_hash = {}
        self.wir_pair_hash = {}
        self.itc_pair_hash_per_corsiv = {}
        
        # set up subplots grid
        tissue_list = sorted( list(self.tissue_to_donor_hash.keys()))
         
        corr_distribution_pair_plot_pdf = "%s.pair_corr.pdf"%self.myArgs.outputRoot  
        plt.figure()
        fig, ax = plt.subplots(self.number_of_tissues, self.number_of_tissues, figsize=(20,20))
    
        for tissue_idx_i in range(self.number_of_tissues):
            gtex_tissue_id_i = tissue_list[tissue_idx_i]
            tissue_i_donors = set(self.tissue_to_donor_hash[gtex_tissue_id_i].keys())
            if (self.DEBUG_TISSUE_PAIR_RHO):
                sys.stderr.write("Processing tissue i %s with %s donors\t"%(gtex_tissue_id_i, len(tissue_i_donors)))
                                  
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
                    sys.stderr.write("Intersection ! [%s = %s] x [%s = %s] ---> %s donors\t"%(gtex_tissue_id_i, len(tissue_i_donors), gtex_tissue_id_j, len(tissue_j_donors), len(common_donors)))
                
                if (number_of_common_donors>=10):
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
                        tissue_i_beta_array = self.coverage_df.iloc[corsiv_idx, tissue_i_common_donors_col_idx]
                        tissue_j_beta_array = self.coverage_df.iloc[corsiv_idx, tissue_j_common_donors_col_idx]
                        rho_value, rho_pvalue = scipy.stats.pearsonr(tissue_i_beta_array, tissue_j_beta_array)
                        for range_idx in range(number_of_common_donors):
                            inter_tissue_range_list.append(abs(tissue_i_beta_array[range_idx]-tissue_j_beta_array[range_idx]))
                            
                        if (self.DEBUG_TISSUE_PAIR_RHO_VERBOSE):
                            sys.stderr.write("Tissues %s x %s: corsiv [%s]: %s --> rho %.4f\n"%(gtex_tissue_id_i, gtex_tissue_id_j, corsiv_idx, corsiv_name, rho_value))
                        
                        inter_tissue_corr_list.append(rho_value)
                        
                        self.itc_pair_hash_per_corsiv[tissue_pair][corsiv_name]=(tissue_i_beta_array, tissue_j_beta_array, rho_value)
                    
                    # plot the damn thing
                    ax[tissue_idx_i, tissue_idx_j].hist(inter_tissue_corr_list, bins=20)
                    ax[tissue_idx_i, tissue_idx_j].set_xlim(-1,1)
                    
                    self.itc_pair_hash[tissue_pair] = inter_tissue_corr_list
                    
                    ax[tissue_idx_j, tissue_idx_i].hist(inter_tissue_range_list, bins=20)
                    ax[tissue_idx_j, tissue_idx_i].set_xlim(0,1)
                    
                    self.wir_pair_hash[tissue_pair] = inter_tissue_range_list

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
                    sys.stderr.write("x [%s, %s] vs [%s, %s] --> %s common tissues\n"%(tissue_idx_i, gtex_tissue_id_i, tissue_idx_j, gtex_tissue_id_j, overlap_size))
                buffer.append(str(overlap_size))
            pair_tissue_count_file_writer.write("%s\t%s\n"%(gtex_tissue_id_i, "\t".join(buffer)))
                
        pair_tissue_count_file_writer.close()
        
        
        # finally plot as long boxplot the itc and the wir
        self.boxPlotSeries(self.itc_pair_hash, "%s.itc_tissue_pairs"%self.myArgs.outputRoot, "ITC over tissue pairs")        
        self.boxPlotSeries(self.wir_pair_hash, "%s.wir_tissue_pairs"%self.myArgs.outputRoot, "Within Individual Range over tissue pairs")
        
        self.plotPerCoRSIVScatterPlots()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP AssessCaptureCORSiV::corrPairPlotRhoDistribution\n"%SimpleTime.now())
       
    def plotPerCoRSIVScatterPlots(self):
        os.system("mkdir -p %s.corsiv_itc"%self.myArgs.outputRoot)
        
        # set up subplots grid
        tissue_list = sorted( list(self.tissue_to_donor_hash.keys()))
        
        plot_idx = 0
        for corsiv_name in  self.corsiv_set:
            if (self.DEBUG_TISSUE_PAIR_SCATTERPLOT):
                sys.stderr.write("Generating scatterplot for %s\n"%corsiv_name)
                
            cleaned_up_corsiv_name = corsiv_name.replace("N/A", "NA").replace(";+","").replace(" ","_").replace(";","_")
                
            corr_distribution_pair_plot_pdf = "%s.corsiv_itc/%s.pair_corr.pdf"%(self.myArgs.outputRoot, cleaned_up_corsiv_name)
#            corr_distribution_pair_plot_jpg = "%s.corsiv_itc/%s.pair_corr.jpg"%(self.myArgs.outputRoot, cleaned_up_corsiv_name)
            
            plt.figure()
            fig, ax = plt.subplots(self.number_of_tissues, self.number_of_tissues, figsize=(20,20))
            
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
                            
                        i_beta_values, j_beta_values, rho_value = self.itc_pair_hash_per_corsiv[tissue_pair][corsiv_name]
                        ax[tissue_idx_i, tissue_idx_j].scatter(j_beta_values, i_beta_values, color='g')
                        b, m = polyfit(j_beta_values, i_beta_values, 1)
                        
                        
                        ax[tissue_idx_i, tissue_idx_j].plot(j_beta_values, b+m*j_beta_values, '-')
                        ax[tissue_idx_i, tissue_idx_j].set_xlim(0,1)
                        ax[tissue_idx_i, tissue_idx_j].set_ylim(0,1)
                        ax[tissue_idx_i, tissue_idx_j].text(0.2,0.2, "PCC R=%.4f"%rho_value)
                    
                        
            
            for tissue_idx_j in range(self.number_of_tissues):
                gtex_tissue_id_j = tissue_list[tissue_idx_j]
                ax[self.number_of_tissues-1, tissue_idx_j].set_xlabel(gtex_tissue_id_j)
            
            for tissue_idx_i in range(self.number_of_tissues):
                gtex_tissue_id_i = tissue_list[tissue_idx_i]
                ax[tissue_idx_i, 0].set_ylabel(gtex_tissue_id_i)
                
            plt.savefig(corr_distribution_pair_plot_pdf, dpi=300)
#            plt.savefig(corr_distribution_pair_plot_jpg, dpi=300)
            
            plot_idx += 1
        
    def boxPlotSeries(self, series_hash, output_file_root, boxplot_title):
        pdf_file = "%s.pdf"%output_file_root
        jpg_file = "%s.jpg"%output_file_root

        fig, ax = plt.subplots(1)
        box_data = []
        labels = []
        for series_key in series_hash:
            box_data.append(series_hash[series_key])
            labels.append(series_key)

        box = plt.boxplot(box_data,patch_artist=True,whis="range")
        ax.set_ylim(0, 1)
        plt.title(boxplot_title)
        plt.xticks(range(1,len(box_data)+1),labels,rotation=90)

        plt.savefig(pdf_file, bbox_inches='tight',dpi=300)
        plt.savefig(jpg_file, bbox_inches='tight',dpi=300)

        
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START AssessCaptureCORSiV::work\n"%SimpleTime.now())
    
            
        self.loadTissueMap()
        self.loadCoverage()
        self.analyzeWir()
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