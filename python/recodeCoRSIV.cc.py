#!/usr/bin/env python

import os, sys

DEBUG_XFORM=True

sys.stderr.write("command: %s\n"%" ".join(sys.argv))

if (len(sys.argv) != 3):
	sys.exit(4)

input_bed = sys.argv[1]
output_bed = sys.argv[2]
sys.stderr.write("in file %s out file %s\n"%(input_bed, output_bed))

in_reader = open(input_bed, "rt")
out_writer = open(output_bed, "wt")

for line in in_reader:
	if (DEBUG_XFORM):
		sys.stderr.write("Processing line %s\n"%line.strip())
	ff = line.strip().split("\t")
	g = ff[3].split(".")
	m = g[2].split("_")
	if (DEBUG_XFORM):
		sys.stderr.write("Gene: %s Cpg Count: %s\n"%(ff[4], m[2]))
	buffer = [ff[0], ff[1], ff[2], "%s.gene_%s"%(ff[3], ff[4]), m[2], ff[5]]
	out_writer.write("%s\n"%"\t".join(buffer))

out_writer.close()
in_reader.close()
