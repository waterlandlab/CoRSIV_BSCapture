#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "2.0"

import os, sys, argparse, re
import glob
import datetime
import math
from argparse import RawTextHelpFormatter
import gzip
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class PlotCaptureDepth:
    DEBUG_PROGRESS              = True
    DEBUG_LOAD_CAPTURE_DEF      = False
    DEBUG_LOAD_OTHER_TARGETS    = False
    DEBUG_COMPUTE_DISTRIBUTION  = False
    DEBUG_PLOT_DISTRIBUTION     = False
    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.
            
Plot distribution of capture sequence depth over multiple target regions

TODO: think of something
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-c','--coverageFiles'  , help='bedtools coverage files pattern',       required=True)
        parser.add_argument('-d','--captureDef'     , help='capture BED definition',                required=True)
        parser.add_argument('-t','--targetFiles'    , help='additional BED target files pattern',   required=False)
        parser.add_argument('-o','--outputRoot'     , help='output files root',                     required=True)
        try:
            args = parser.parse_args()
        except:
            args = None
        return args


    def loadCaptureDefinition(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START PlotCaptureDepth::loadCaptureDefinition\n"%SimpleTime.now())

        self.capture_regions = set()
        
        capture_def_reader = open(self.myArgs.captureDef)
        
        for line in capture_def_reader:
            ff = line.split("\t")
            region_name = ff[3]
            self.capture_regions.add(region_name)
            if (self.DEBUG_LOAD_CAPTURE_DEF):
                sys.stderr.write("Loading capture def %s\n"%region_name)
                
        capture_def_reader.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START PlotCaptureDepth::loadCaptureDefinition loaded %s regions\n"%(SimpleTime.now(), len(self.capture_regions)))
    
    
    def prepareAdditionalTargets(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START PlotCaptureDepth::prepareAdditionalTargets\n"%SimpleTime.now())
        self.other_targets_hash = {}
        
        # add all targets as obvious subset
        self.other_targets_hash["all_targets"] = self.capture_regions
        
        other_targets_list = list(glob.glob(self.myArgs.targetFiles))
        
        if (self.DEBUG_LOAD_OTHER_TARGETS):
            sys.stderr.write("Other BED targets: %s\n"%";".join(other_targets_list))
            
        if (len(other_targets_list)>0):
            for target_bed in other_targets_list:
                short_name = os.path.basename(target_bed)
                
                if (self.DEBUG_LOAD_OTHER_TARGETS):
                    sys.stderr.write("Determining overlap with %s short name %s\n"%(target_bed, short_name))
                    
                # run bedtools intersect
                overlap = "%s.capture_overlap_with_%s.bed"%(self.myArgs.captureDef, short_name)
                bedtools_intersect_command = "bedtools intersect -u -a %s -b %s > %s"%(self.myArgs.captureDef, target_bed, overlap)
                
                if (self.DEBUG_LOAD_OTHER_TARGETS):
                    sys.stderr.write("bedtools command %s\n"%bedtools_intersect_command)
                    
                os.system(bedtools_intersect_command)
                
                # load results
                overlap_region_set = set()
                overlap_reader = open(overlap)
                for line in overlap_reader:
                    ff = line.strip().split("\t")
                    region_name = ff[3]
                    overlap_region_set.add(region_name)
                
                overlap_reader.close()
                
                # if overlap count > 0; add common regions
                if (len(overlap_region_set)>0):
                    if (self.DEBUG_LOAD_OTHER_TARGETS):
                        sys.stderr.write("overlap of capture def %s with %s yielded %s regions\n"%(self.myArgs.captureDef, target_bed, len(overlap_region_set)))
                        
                    self.other_targets_hash[short_name]=overlap_region_set
                                    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP PlotCaptureDepth::prepareAdditionalTargets have %s target sets \n"%(SimpleTime.now(), len(self.other_targets_hash)))
    
    
    def updateOneCoverage(self, one_coverage):
        if (self.DEBUG_COMPUTE_DISTRIBUTION):
                sys.stderr.write("[%s] START Updating depth distribution from %s\n"%(SimpleTime.now(), one_coverage))
            
        count_hash = {}    
        for key in self.other_targets_hash:
            count_hash[key]={}
            for depth_target in self.depth_targets:
                count_hash[key][depth_target]=0
            
        coverage_reader = open(one_coverage)
        
        for line in coverage_reader:
            ff = line.strip().split("\t")
            chrom_start = float(ff[1])
            chrom_stop  = float(ff[2])
            region_size = chrom_stop - chrom_start + 1
            region_name = ff[3]
            read_count = float(ff[6])
            read_depth = read_count*150/region_size
            
            if (self.DEBUG_COMPUTE_DISTRIBUTION):
                sys.stderr.write("Coverage %s line %s --> %sx\n"%(one_coverage, line.strip(), read_depth))
            
            # now go through all count_hash and update as warranted
            for key in self.other_targets_hash:
                regions_set = self.other_targets_hash[key]
                if (region_name in regions_set):
                    # if (self.DEBUG_COMPUTE_DISTRIBUTION):
                    #    sys.stderr.write("[%s] region %s found in targets set %s\n"%(one_coverage, region_name, key))
                    
                    for depth_target in self.depth_targets:
                        if (read_depth>=depth_target):
                            if (self.DEBUG_COMPUTE_DISTRIBUTION):
                                sys.stderr.write("[%s][%s][%s] has count %s\n"%(one_coverage, key, depth_target, count_hash[key][depth_target]))
                            
                            count_hash[key][depth_target] += 1
                            
                            if (self.DEBUG_COMPUTE_DISTRIBUTION):
                                sys.stderr.write("[%s][%s][%s] region %s depth %s exceeds %s --> count %s\n"%(one_coverage, key, depth_target, region_name, read_depth, depth_target, count_hash[key][depth_target]))
    
        coverage_reader.close()
        
        # collect fraction for each threshold
        for key in self.other_targets_hash:
            regions_in_target = len(self.other_targets_hash[key])
            if (self.DEBUG_COMPUTE_DISTRIBUTION):
                sys.stderr.write("Updating targets %s with %s regions\n"%(key, regions_in_target))
            
            for depth_target in self.depth_targets:
                count_exceeding_depth = count_hash[key][depth_target]
                fraction_regions_exceeding_depth = float(count_exceeding_depth)/float(regions_in_target)
                self.exceed_coverage_hash[key][depth_target].append(fraction_regions_exceeding_depth)
                
                if (self.DEBUG_COMPUTE_DISTRIBUTION):
                    sys.stderr.write("Updating target (%s, %s regions, depth %s) with count %s fraction %s\n"%(key, regions_in_target, depth_target, count_exceeding_depth, fraction_regions_exceeding_depth))
                
        if (self.DEBUG_COMPUTE_DISTRIBUTION):
                sys.stderr.write("[%s] STOP Updating depth distribution from %s\n"%(SimpleTime.now(), one_coverage))
        
    def computeDistributions(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START PlotCaptureDepth::computeDistributions\n"%SimpleTime.now())
            
        self.depth_targets = [1, 5, 10, 20, 30, 40, 50, 100, 150, 200]
    
        self.exceed_coverage_hash = {}
        for key in self.other_targets_hash:
            self.exceed_coverage_hash[key]={}
            for depth_target in self.depth_targets:
                self.exceed_coverage_hash[key][depth_target]=[]
        
        coverage_files = list(glob.glob(self.myArgs.coverageFiles))
        
        if (self.DEBUG_COMPUTE_DISTRIBUTION):
            sys.stderr.write("[%s] found %s coverage files\n"%(SimpleTime.now(), len(coverage_files)))
            
        for one_coverage in coverage_files:
                
            self.updateOneCoverage(one_coverage)
                
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START PlotCaptureDepth::computeDistributions from %s coverage files\n"%(SimpleTime.now(), len(coverage_files)))
    
    
    def plotOneTargetSet(self, key):
        pdf_file = "%s.overlap_%s.pdf"%(self.myArgs.outputRoot, key) 
        jpg_file = "%s.overlap_%s.pdf"%(self.myArgs.outputRoot, key) 
        
        fig, ax = plt.subplots(1)
        box_data = []
        labels = []
        for depth_target in self.depth_targets:
            box_data.append(self.exceed_coverage_hash[key][depth_target])
            labels.append(str(depth_target))
        
        box = plt.boxplot(box_data,patch_artist=True,whis="range")    
        ax.set_ylim(0, 1)
        plt.title("Sequencing depth for %s"%key)
        plt.xticks(range(1,len(box_data)+1),labels,rotation=90)
            
        plt.savefig(pdf_file, bbox_inches='tight',dpi=300)
        plt.savefig(jpg_file, bbox_inches='tight',dpi=300)
                
    def plotDistributions(self):
        for key in self.exceed_coverage_hash:
            self.plotOneTargetSet(key)
            
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START PlotCaptureDepth::work\n"%SimpleTime.now())

        self.loadCaptureDefinition()    
        self.prepareAdditionalTargets()
        self.computeDistributions()
        self.plotDistributions()
    
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP PlotCaptureDepth::work\n"%SimpleTime.now())
        

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = PlotCaptureDepth.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = PlotCaptureDepth(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
