#!/bin/python

import argparse
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from Bio import Phylo

p = argparse.ArgumentParser(description='Calculate Newick Tree from Distance Matrix', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('DISTMAT', help='File containing Distance Matrix (lower triangular)')
p.add_argument('--method', help='Tree construction algorithm', default='nj', choices=['nj','upgma'])
p.add_argument('--out', help='Output files name (tree+mapping)', default='tmp_tree')
p.add_argument('--draw', help='Draw a picture of the calculated tree', action='store_true')
p.add_argument('--outfmt', help='Output Tree file format', choices=['newick', 'phyloxml', 'nexus'], default='phyloxml')
p.add_argument('--distout', help='Print Distance Matrix to STDOUT', action='store_true')
args = p.parse_args()

def main():
    dist_mat = ParseMatrix(args.DISTMAT)

    if args.distout is True:
        print 'Distance Matrix:'
        print dist_mat

    tree_constructor = DistanceTreeConstructor()
    if args.method == 'nj':
        tree = tree_constructor.nj(dist_mat)
    elif args.method == 'upgma':
        tree = tree_constructor.upgma(dist_mat)

    if args.draw is True:
        Phylo.draw(tree)

    #Write NEWICK file
    Phylo.write(tree, args.out+'.tree', args.outfmt)


### FUNCTIONS ###
def ParseMatrix(filename):
    mat_names = [] # FASTA headers
    mat_names_num = []
    lt_matrix = [] #lower triangular matrix
    with open(filename, 'rU' ) as MAT:
        for l in MAT:
            l = l.strip('\n')
            if len(l) == 0:
                continue
            elif l[0] == '>':
                mat_names.append(l[1:])
            else:
                lt_matrix.append([float(i) for i in l.split()])

    #Switch to integer headers and print Mapping file
    with open(args.out+'.map', 'w') as MAP:
        for index, name in enumerate(mat_names):
            MAP.write(str(index)+'\t'+name+'\n')
            mat_names_num.append(str(index))

    #Fill into Biopython distmat Data Structures
    dist_matrix = _DistanceMatrix(names=mat_names_num, matrix=lt_matrix)
    return dist_matrix


if __name__ == '__main__':
    main()
