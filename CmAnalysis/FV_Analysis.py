#!/bin/python

from CmAnalysis import CmAnalysis
import CmDataStrucs as CD
import Parsing_Utils as PU
import argparse, sys

#Command line Argument Parsing
p=argparse.ArgumentParser(description='Analyze Flavivirus CMsearch runs')
p.add_argument('TB', help='CMsearch TBLOUT file')
p.add_argument('FA', help='Source seqence file')
p.add_argument('--threshold', type=float, default=0.001, help='Minimum E-Value (default 0.001)')
#p.add_argument('-o', type=str, help='Output file name', default='vecs_summary.out')
args = p.parse_args()

analysis = CmAnalysis(args.TB, args.FA)

for seqname in analysis.MapLoci('seqname'):
    print seqname
    for locus in analysis.MapLoci('seqname')[seqname]:
        print locus.range, locus.cms
        print locus.evalues

sys.exit()
# for compound in analysis.unique_loci:
#     print compound.seqname
#     print compound.cms
#     print compound.range
#     print compound.ranges

#Print Vectors
#analysis.cmVectors('MVV_all.vecs')
