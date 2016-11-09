#!/bin/python

from CmAnalysis import CmAnalysis
import CmDataStrucs as CD
import Parsing_Utils as PU
import argparse

#Command line Argument Parsing
p=argparse.ArgumentParser(description='Analyze Flavivirus CMsearch runs')
p.add_argument('TB', help='CMsearch TBLOUT file')
p.add_argument('FA', help='Source seqence file')
p.add_argument('--threshold', type=float, default=0.001, help='Minimum E-Value (default 0.001)')
#p.add_argument('-o', type=str, help='Output file name', default='vecs_summary.out')
args = p.parse_args()

analysis = CmAnalysis(args.TB, args.FA)

#loci ausprinten fuer jedes genom
#analysis.cm_hits
#analysis.unique_loci

#
