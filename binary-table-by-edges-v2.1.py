#!/usr/bin/python
# coding=utf-8
#
# C.D. Dunn†, B.A. Akpınar, and V. Sharma†. An unusual amino acid substitution within the hummingbird
# cytochrome c oxidase alters a key proton-conducting channel. bioRxiv. doi: 10.1101/610915.

# This script calculates fluctuations for all internal edges [node-edge-node] in a given tree
# for each amino acid position in a given alignment.

# The input files are tree and alignment files, with internal node labels and amino acid predictions, such as
# PAGAN-generated ancestral alignment and tree files.
# In the output, '1' indicates that the corresponding edge fluctuates for the given position,
# and '0' indicates otherwise.

# Script usage: $python binary-table-by-edges-v2.py -t [tree_file] -a [alignment_file]

# Dependencies: Bio.SeqIO and Bio.Phylo

import argparse
from Bio import Phylo
from Bio import SeqIO
import time

parser = argparse.ArgumentParser(prog='Binary table by edges',
                                 usage='$python binary-table-by-edges-v2.1.py -t [tree_file] -a [alignment_file]')

parser.add_argument('--treefile', '-t',
                    help='Please specify a NWK tree file with internal node labels')
parser.add_argument('--alignmentfile', '-a',
                    help='Please specify the alignment file used to generate the tree file')


args = parser.parse_args()

tree_file = args.treefile
alignment_file = args.alignmentfile

##########################################################################################
#################### Functions to be used in the script ##################################

# To load the FASTA alignment file


def read_fasta(alignment):  # To read the FASTA alignment file

    aa_dict = {}
    with open(alignment, mode='r') as handle:

        # Using Biopython's parse function to reduce memory footprint
        for record in SeqIO.parse(handle, 'fasta'):

            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = record.seq
            aa_dict[identifier] = sequence

    return aa_dict

# To define all edges in a given tree


def all_edges(tree):

    alledges = []
    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            alledges.append(parent.name + '*' + child.name)

    return alledges

##########################################################################################


# Read the NWK tree file and retrieve all edges in the tree
my_tree = Phylo.read(tree_file, 'newick')
edges_all = sorted(all_edges(my_tree))

print "The total number of edges in the tree is: %d" % len(edges_all)

# Read the FASTA file and record amino acid values for each amino acid position
all_values = read_fasta(alignment_file)
prot_length = len(all_values[all_values.keys()[0]])
print "The number of positions to be analyzed is: %d" % prot_length
print "..."

# Open a new file to output results
try:
    results = open('BinaryF-table_%s.txt' % (tree_file.split('.')[0]), 'w')
except:
    results = open('BinaryF-table_%s.txt' % tree_file, 'w')

header = '\t' + '\t'.join(edges_all) + '\n'
results.write(header)

# For all-identical positions, all edges are assigned 0's
data_line_identical = '\t'.join(['0']*len(edges_all))

start = time.time()
for pos in range(0, prot_length):  # Process one position at a time

    # Retrieve all values for the given position for all nodes and leaves
    pos_values = {}
    for sp in all_values.keys():
        pos_values[sp] = all_values[sp][pos].upper()

    # Decide if the given position is all identical across all nodes and leaves
    if len(list(set(pos_values.values()))) == 1:
        results.write('Pos%r' % (pos + 1) + '\t' + data_line_identical + '\n')

    # If not, compute fluctuations
    else:
        mutated_edges = []

        for item in edges_all:
            node1 = item.split('*')[0]
            node2 = item.split('*')[1]
            if pos_values[node1] != pos_values[node2]:
                mutated_edges.append('1')  # Fluctuating edge

            else:
                mutated_edges.append('0')  # Non-fluctuating edge

        results.write('Pos%r' % (pos + 1) + '\t' + '\t'.join(mutated_edges) + '\n')

end = time.time()

print "Computation took %r seconds." % (end - start)

results.close()

print "Analysis complete."
