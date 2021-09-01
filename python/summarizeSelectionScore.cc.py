#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter
import glob as glob
import numpy as np
import scipy.stats
# import pandas as pd
import seaborn as sn

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

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
class SummarizeSelectionScore:
	DEBUG_PROGRESS              = True
	DEBUG_LOAD_BED_FILE    		= True
	DEBUG_LOAD_BED_FILE_SCORES  = True
	DEBUG_OUTPUT = True


	def __init__(self, myArgs):
		self.myArgs = myArgs

	@staticmethod
	def processArguments():
		parser = argparse.ArgumentParser(description=
	"""\
	Utility %s version %s.

	compute R2 for a CORSIV and mQTL
	"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)

		parser.add_argument('-b','--bedFilePattern',   		help='bed file(s) pattern',  required=True)
		parser.add_argument('-x','--selectionScoreFile',    help='target tissue',  required=True)
		parser.add_argument('-o','--outputRoot',            help='output file root',   required=True)

		try:
			args = parser.parse_args()
		except:
			args = None
		return args

	def loadSelectionScores(self, bed_file_long, bed_file_short):
		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] START SummarizeSelectionScore::loadSelectionScores for %s\n"%(SimpleTime.now(), bed_file_long))

		# setup bedtools intersect command to obtain the selection scores
		score_file = "%s.tmp.%s.scores.txt"%(self.myArgs.outputRoot, bed_file_short)
		bedtool_command = "bedtools intersect -wa -wb -a %s -b %s > %s"%(bed_file_long, self.myArgs.selectionScoreFile, score_file)
		if (self.DEBUG_LOAD_BED_FILE_SCORES):
			sys.stderr.write("bedtool command %s\n"%bedtool_command)
		os.system(bedtool_command)

		score_set = []
		score_reader = open(score_file, "rt")
		for line in score_reader:
			ff = line.strip().split()
			score = float(ff[10])
			if (self.DEBUG_LOAD_BED_FILE_SCORES):
				sys.stderr.write("Line |%s| led to score %s\n"%(line.strip(), score))
			score_set.append(score)

		score_reader.close()

		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] STOP SummarizeSelectionScore::loadSelectionScores for %s\n"%(SimpleTime.now(), bed_file_long))

		return score_set

	def loadBedFileAndScores(self):
		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] START SummarizeSelectionScore::loadBedFileAndScores\n"%SimpleTime.now())

		self.bed_file_list = glob.glob(self.myArgs.bedFilePattern)
		self.bed_file_score_hash = {}

		for bed_file_long in self.bed_file_list:
			bed_file_base = os.path.basename(bed_file_long)
			if (self.DEBUG_LOAD_BED_FILE):
				sys.stderr.write("Processing bed file long %s short %s\n"%(bed_file_long, bed_file_base))

			self.bed_file_score_hash[bed_file_base]=self.loadSelectionScores(bed_file_long, bed_file_base)

		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] STOP SummarizeSelectionScore::loadBedFileAndScores\n"%SimpleTime.now())

	def outputResults(self):
		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] START SummarizeSelectionScore::outputResults\n"%SimpleTime.now())

		data_frame_list_of_lists = []
		sorted_beds = sorted(self.bed_file_score_hash.keys())
		for bed_file_short in sorted_beds:
			bed_file_score_list = self.bed_file_score_hash[bed_file_short]

			for score in bed_file_score_list:
				data_frame_list_of_lists.append([bed_file_short, score])

			columns_bp_coverage = ["BedFile", "TajimaDScore"]
			self.score_data_frame = pd.DataFrame( data_frame_list_of_lists, columns = columns_bp_coverage)


#		violin_plot_pdf = "%s.violin_plot.pdf"%self.myArgs.outputRoot
#		sn.set(style="whitegrid")
#		ax = sn.violinplot(x="BedFile", y="TajimaDScore", data=self.score_data_frame)
#		#ax.set_ylim(-1,1)
#		ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
#		plt.savefig(violin_plot_pdf, bbox_inches='tight', dpi=150)

		plt.figure()
		box_plot_pdf = "%s.box_plot.pdf"%self.myArgs.outputRoot
		ax = sn.boxplot(x="BedFile", y="TajimaDScore", data=self.score_data_frame)
		#ax.set_ylim(-1,1)
		ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
		plt.savefig(box_plot_pdf, bbox_inches='tight', dpi=150)

		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] START SummarizeSelectionScore::outputResults\n"%SimpleTime.now())



	def work(self):
		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] START SummarizeSelectionScore::work\n"%SimpleTime.now())

		self.loadBedFileAndScores()
		self.outputResults()

		if (self.DEBUG_PROGRESS):
			sys.stderr.write("[%s] STOP SummarizeSelectionScore::work\n"%SimpleTime.now())

########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = SummarizeSelectionScore.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = SummarizeSelectionScore(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
