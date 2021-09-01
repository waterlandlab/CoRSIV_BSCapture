#!/usr/bin/env python
__author__ = "Cristian Coarfa"

__version__ = "1.0"

import os, sys, argparse, re
import datetime
from argparse import RawTextHelpFormatter

from CCLabUtils.simpleTime import SimpleTime

class CStruct(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
class BedSplit_ChromosomeWindows:
    DEBUG_PROGRESS     = True
    DEBUG_SPLIT        = True

    
    def __init__(self, myArgs):
        self.myArgs = myArgs

    @staticmethod
    def processArguments():
        parser = argparse.ArgumentParser(description=
"""\
Utility %s version %s.

split a bed file in equal size parts per chromosome
assume small number of chromosomes (<2048 or so)
"""%(os.path.basename(sys.argv[0]), __version__), formatter_class=RawTextHelpFormatter)
        
        parser.add_argument('-b','--bedFile',              help='bed file',   required=True)
        parser.add_argument('-s','--splitLines',            help='number of records in split bed files',   required=True)
        parser.add_argument('-o','--outputRoot',            help='output file root',         required=True)
        
        try:
            args = parser.parse_args()
        except:
            args = None
        return args

    def splitBedFile(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BedSplit_ChromosomeWindows::splitBedFile\n"%SimpleTime.now())
        
        self.chrom_fd_hash = {} 
        max_chunk_size = int(self.myArgs.splitLines)
        if (max_chunk_size<1):
            sys.stderr.write("Split Chunk size should be greater then or equal to 1\n.Exiting ...\n")
            sys.exit(4)

        if (self.DEBUG_SPLIT):
            sys.stderr.write("max chunk size %s\n"%max_chunk_size)
            
        bed_reader = open(self.myArgs.bedFile)
        
        for line in bed_reader:
            ff = line.strip().split('\t')
            chrom = ff[0]
            if not (chrom in self.chrom_fd_hash):
                chunk_index = 0
                chunk_name = "%s.chrom.%s.chunk.%s.bed"%(self.myArgs.outputRoot, chrom, chunk_index)
                fd = open(chunk_name, "wt")
                chrom_info = CStruct(chrom = chrom, chunk_index = 0, chunk_size = 0, fd = fd, chunk_name = chunk_name)
                self.chrom_fd_hash[chrom]=chrom_info
                if (self.DEBUG_SPLIT):
                    sys.stderr.write("Add new chunk for chrom %s chunk_index %s chunk_size %s chunk_file %s \n"%(chrom_info.chrom, chrom_info.chunk_index, chrom_info.chunk_size, chrom_info.chunk_name))
                
            chrom_info = self.chrom_fd_hash[chrom]
            chrom_info.fd.write(line)
            chrom_info.chunk_size += 1
            
            if (self.DEBUG_SPLIT):
                    sys.stderr.write("W: chrom %s chunk_index %s crt chunk_size %s vs max_chunk_size %s\n"%(chrom_info.chrom, chrom_info.chunk_index, chrom_info.chunk_size, max_chunk_size))
                    
            if (chrom_info.chunk_size==max_chunk_size):
                chrom_info.fd.close()
                chrom_info.chunk_index +=1
                chunk_name = "%s.chrom.%s.chunk.%s.bed"%(self.myArgs.outputRoot, chrom, chrom_info.chunk_index)
                chrom_info.fd = open(chunk_name, "wt")
                chrom_info.chunk_size = 0
                chrom_info.chunk_name = chunk_name
                if (self.DEBUG_SPLIT):
                    sys.stderr.write("Updated chunk for chrom %s chunk_index %s chunk_size %s chunk_file %s \n"%(chrom_info.chrom, chrom_info.chunk_index, chrom_info.chunk_size, chrom_info.chunk_name))
                
        bed_reader.close()
        
        for chrom in self.chrom_fd_hash:
            chrom_info = self.chrom_fd_hash[chrom]
            chrom_info.fd.close()
        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BedSplit_ChromosomeWindows::splitBedFile\n"%SimpleTime.now())
    
    def work(self):
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] START BedSplit_ChromosomeWindows::work\n"%SimpleTime.now())
        
        self.splitBedFile()        
        if (self.DEBUG_PROGRESS):
            sys.stderr.write("[%s] STOP BedSplit_ChromosomeWindows::work\n"%SimpleTime.now())
    
########################################################################################
# MAIN
########################################################################################

# Process command line options
## Instantiate analyzer using the program arguments
## Analyze this !

if __name__ == '__main__':
    try:
        sys.stderr.write("Command line %s\n"%" ".join(sys.argv))
        myArgs = BedSplit_ChromosomeWindows.processArguments()
        if (myArgs is None):
            pass
        else:
            bp = BedSplit_ChromosomeWindows(myArgs)
            bp.work()
    except:
        sys.stderr.write("An unknown error occurred.\n")
        raise
