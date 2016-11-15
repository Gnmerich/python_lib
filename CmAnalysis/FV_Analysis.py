#!/bin/python

from CmAnalysis import CmAnalysis
import CmDataStrucs as CD
import Parsing_Utils as PU
import argparse, sys

#Command line Argument Parsing
p=argparse.ArgumentParser(description='Analyze Flavivirus CMsearch runs')
p.add_argument('TB', help='CMsearch TBLOUT file')
p.add_argument('FA', help='Source seqence file')
#p.add_argument('--threshold', type=float, default=0.001, help='Minimum E-Value (default 0.001)')
#p.add_argument('-o', type=str, help='Output file name', default='vecs_summary.out')
args = p.parse_args()

#Build Analysis Object
analysis = CmAnalysis(args.TB, args.FA)

#Iterate over all cmsearch-hits and print them in BED-format to STDOUT:
for hit in analysis.cm_hits:
    print hit.BedStr()

#Iterate over all unique (compound-) loci and print information to STDOUT
for unique_locus in analysis.unique_loci:
    print 'Sequence Name:', unique_locus.seqname
    print 'Covar Models:', unique_locus.cms
    print 'Range:', unique_locus.range
    print 'Strand:', unique_locus.strand
    #Print all hits of this locus in BED-format to STDOUT:
    print unique_locus.BedStr()
