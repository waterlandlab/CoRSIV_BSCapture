#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import pandas as pd
import numpy as np
import glob as glob

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class Annotate_CoRSIV_nearest_gene:
    DEBUG_PROGRESS          = True
    DEBUG_LOAD_CORSIV_DEF   = True
    DEBUG_LOAD_GENE_DEF     = True
    DEBUG_NEAREST_GENE      = True
    DEBUG_NEAREST_GENE_VERBOSE  = False 

    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

annotate nearest gene for corsivs
TODO: implement something smarter than brute force gene search
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)

        parser.add_argument('-b','--corsivBedFile',   help='corsiv bed file',               required=True)
        parser.add_argument('-g','--geneBedFile',     help='gene body definition (BED)',    required=True)
        parser.add_argument('-c','--chromosome',      help='chromosome to process',         required=True)
        parser.add_argument('-o','--outputRoot',      help='output file root',              required=True)

        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def loadCoRSIVDef(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_nearest_gene::loadCoRSIVDef %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))

        self.corsiv_def_hash = {}

        with open(self.myArgs.corsivBedFile, "rt") as corsiv_reader:
            for line in corsiv_reader:

                ff = line.strip().split('\t')
                chrom = ff[0]
                if (chrom != self.myArgs.chromosome):
                    continue

                corsiv_short = ff[3]
                chrom_start = int(ff[1])
                chrom_stop  = int(ff[2])

                corsiv_info = CStruct(corsiv_short = corsiv_short, chrom = chrom, chrom_start = chrom_start, chrom_stop = chrom_stop,
                            distance_left="NA", gene_left="", distance_right="NA", gene_right="",
                            nearest_abs_distance="NA", nearest_abs_gene="" )
                self.corsiv_def_hash[corsiv_short]=corsiv_info

                if (self.DEBUG_LOAD_CORSIV_DEF):
                    sys.stderr.write("Loaded corsiv definition  short: %s chrom/start/stop: %s/%s/%s \n"%(corsiv_short,
                            self.corsiv_def_hash[corsiv_short].chrom,
                            self.corsiv_def_hash[corsiv_short].chrom_start, self.corsiv_def_hash[corsiv_short].chrom_stop))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_nearest_gene::loadCoRSIVDef %s"%(SimpleTime.now(), self.myArgs.corsivBedFile))

    def loadGeneDef(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_nearest_gene::loadGeneDef %s"%(SimpleTime.now(), self.myArgs.geneBedFile))

        self.gene_def_hash = {}

        with open(self.myArgs.geneBedFile, "rt") as gene_reader:
            for line in gene_reader:

                ff = line.strip().split('\t')
                chrom = ff[0]
                if (chrom != self.myArgs.chromosome):
                    continue

                gene_symbol = ff[4]
                gene_ens = ff[3]
                chrom_start = int(ff[1])
                chrom_stop  = int(ff[2])

                gene_info = CStruct(gene_ens = gene_ens, gene_symbol=gene_symbol,
                            chrom = chrom, chrom_start = chrom_start, chrom_stop = chrom_stop)

                self.gene_def_hash[gene_ens]=gene_info

                if (self.DEBUG_LOAD_GENE_DEF):
                    sys.stderr.write("Loaded gene  ens %s symbol %s chrom/start/stop: %s/%s/%s \n"%(gene_ens, gene_symbol,
                            self.gene_def_hash[gene_ens].chrom, self.gene_def_hash[gene_ens].chrom_start, self.gene_def_hash[gene_ens].chrom_stop))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_nearest_gene::loadGeneDef %s"%(SimpleTime.now(), self.myArgs.geneBedFile))

    def determineNearestGeneOneCorsiv(self, corsiv_id):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_nearest_gene::determineNearestGeneOneCorsiv %s\n"%(SimpleTime.now(), corsiv_id))

        corsiv_info = self.corsiv_def_hash[corsiv_id]

        for gene_ens in self.gene_def_hash:
            gene_info = self.gene_def_hash[gene_ens]

            if (gene_info.chrom_start > corsiv_info.chrom_stop):
                distance_right = gene_info.chrom_start -corsiv_info.chrom_stop
                if (corsiv_info.gene_right ==""):
                    corsiv_info.distance_right = distance_right
                    corsiv_info.gene_right = gene_info.gene_ens
                    if (self.DEBUG_NEAREST_GENE_VERBOSE):
                        sys.stderr.write("Initialize_right corsiv %s right gene %s dist %s\n"%(corsiv_id, corsiv_info.gene_right, corsiv_info.distance_right))
                elif (abs(distance_right)<abs(corsiv_info.distance_right)):
                    corsiv_info.distance_right = distance_right
                    corsiv_info.gene_right = gene_info.gene_ens
                    if (self.DEBUG_NEAREST_GENE_VERBOSE):
                        sys.stderr.write("Update_right corsiv %s right gene %s dist %s\n"%(corsiv_id, corsiv_info.gene_right, corsiv_info.distance_right))
            elif (gene_info.chrom_stop<corsiv_info.chrom_start):
                distance_left = gene_info.chrom_stop -corsiv_info.chrom_start
                if (corsiv_info.gene_left ==""):
                    corsiv_info.distance_left = distance_left
                    corsiv_info.gene_left = gene_info.gene_ens
                    if (self.DEBUG_NEAREST_GENE_VERBOSE):
                        sys.stderr.write("Initialize_left corsiv %s left gene %s dist %s\n"%(corsiv_id, corsiv_info.gene_left, corsiv_info.distance_left))
                elif (abs(distance_left)<abs(corsiv_info.distance_left)):
                    corsiv_info.distance_left = distance_left
                    corsiv_info.gene_left = gene_info.gene_ens
                    if (self.DEBUG_NEAREST_GENE_VERBOSE):
                        sys.stderr.write("Update_left corsiv %s left gene %s dist %s\n"%(corsiv_id, corsiv_info.gene_left, corsiv_info.distance_left))
            else:
                # overlap !!
                corsiv_info.distance_left = 0
                corsiv_info.gene_left = gene_info.gene_ens
                corsiv_info.distance_right = 0
                corsiv_info.gene_right = gene_info.gene_ens
                if (self.DEBUG_NEAREST_GENE_VERBOSE):
                    sys.stderr.write("Update_overlap corsiv %s left gene %s left dist %s right gene %s right dist %s\n"%(corsiv_id,
                        corsiv_info.gene_left, corsiv_info.distance_left, corsiv_info.gene_right, corsiv_info.distance_right))

        # determine also the best distance/best gene
        # note: we expect to find at least one of the gene left and gene right
        if (corsiv_info.gene_left=="") and (corsiv_info.gene_right==""):
            sys.stderr.write("Eeek !!!! corsiv %s gene left and right null\n"%corsiv_id)
            sys.exit(4)
        elif (corsiv_info.gene_left==""):
            corsiv_info.nearest_abs_gene = corsiv_info.gene_right
            corsiv_info.nearest_abs_distance = corsiv_info.distance_right

            if (self.DEBUG_NEAREST_GENE):
                sys.stderr.write("Best_right_default %s %s %s\n"%(corsiv_id, corsiv_info.nearest_abs_gene, corsiv_info.nearest_abs_distance))

        elif (corsiv_info.gene_right==""):
            corsiv_info.nearest_abs_gene = corsiv_info.gene_left
            corsiv_info.nearest_abs_distance = corsiv_info.distance_left

            if (self.DEBUG_NEAREST_GENE):
                sys.stderr.write("Best_left_default %s %s %s\n"%(corsiv_id, corsiv_info.nearest_abs_gene, corsiv_info.nearest_abs_distance ))

        else: # have both right and left distance
            if (abs(corsiv_info.distance_left)<abs(corsiv_info.distance_right)):
                corsiv_info.nearest_abs_gene = corsiv_info.gene_left
                corsiv_info.nearest_abs_distance = corsiv_info.distance_left
                if (self.DEBUG_NEAREST_GENE):
                    sys.stderr.write("Best_left_distance %s: %s %s : proof left %s %s vs right %s %s\n"%(corsiv_id,
                        corsiv_info.nearest_abs_gene, corsiv_info.nearest_abs_distance,
                        corsiv_info.distance_left, corsiv_info.gene_left, corsiv_info.distance_right, corsiv_info.gene_right))
            else:
                corsiv_info.nearest_abs_gene = corsiv_info.gene_right
                corsiv_info.nearest_abs_distance = corsiv_info.distance_right
                if (self.DEBUG_NEAREST_GENE):
                    sys.stderr.write("Best_right_distance %s: %s %s : proof left %s %s vs right %s %s\n"%(corsiv_id,
                        corsiv_info.nearest_abs_gene, corsiv_info.nearest_abs_distance,
                        corsiv_info.distance_left, corsiv_info.gene_left, corsiv_info.distance_right, corsiv_info.gene_right))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_nearest_gene::determineNearestGeneOneCorsiv %s\n"%(SimpleTime.now(), corsiv_id))

    def determineNearestGene(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_nearest_gene::determineNearestGene\n"%SimpleTime.now())

        for corsiv_id in self.corsiv_def_hash:
            self.determineNearestGeneOneCorsiv(corsiv_id)

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_nearest_gene::determineNearestGene\n"%SimpleTime.now())


    def outputNearestGene(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_nearest_gene::outputNearestGene\n"%SimpleTime.now())

        corsiv_distance_report = "%s.nearest_gene.xls"%self.myArgs.outputRoot

        with open(corsiv_distance_report, "wt") as corsiv_gene_writer:
            corsiv_gene_writer.write("CoRSIV\tChrom\tStart\tStop\tNearest Left Gene\tNearest Left Distance\tNearest Right Gene\tNearest Right Distance\tNearest Gene\tNearest Distance\n")

            for corsiv_id in self.corsiv_def_hash:
                corsiv_info = self.corsiv_def_hash[corsiv_id]
                corsiv_buffer = [corsiv_id, corsiv_info.chrom, corsiv_info.chrom_start, corsiv_info.chrom_stop,
                                corsiv_info.distance_left, corsiv_info.gene_left,
                                corsiv_info.distance_right, corsiv_info.gene_right,
                                corsiv_info.nearest_abs_distance, corsiv_info.nearest_abs_gene]

                corsiv_gene_writer.write("%s\n"%"\t".join([str(x) for x in corsiv_buffer]))

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_nearest_gene::outputNearestGene\n"%SimpleTime.now())

    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START Annotate_CoRSIV_nearest_gene::work\n"%SimpleTime.now())

        self.loadCoRSIVDef()
        self.loadGeneDef()
        self.determineNearestGene()
        self.outputNearestGene()

        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP Annotate_CoRSIV_nearest_gene::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = Annotate_CoRSIV_nearest_gene.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = Annotate_CoRSIV_nearest_gene(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
