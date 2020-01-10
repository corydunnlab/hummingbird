#!/usr/bin/python
# coding=utf-8
#
# C.D. Dunn†, B.A. Akpınar, and V. Sharma†. An unusual amino acid substitution within the hummingbird
# cytochrome c oxidase alters a key proton-conducting channel. bioRxiv. doi: 10.1101/610915.

# This script calculates substitutions for all edges of a given tree for each amino acid position in a given alignment.

# The input files are FASTA alignment and NWK tree files, with internal node labels and amino acid predictions, such as
# PAGAN-generated ancestral alignment and tree files.
# In the output, '1' indicates a substitution in the corresponding edge for the given position,
# and '0' indicates otherwise.

# UPDATE: v3 adds convention for a species of interest automatically, if specified (optional).

# Script usage:
# $python binary-table-by-edges-v3.py -t [tree_file] -a [alignment_file] -conv <species_of_interest> (optional)

# Dependencies: Bio.SeqIO and Bio.Phylo

import argparse
from Bio import Phylo
from Bio import SeqIO
import time

parser = argparse.ArgumentParser(prog='binary-table-by-edges-v3',
                                 usage='$python binary-table-by-edges-v3.py -t [tree_file] -a [alignment_file] '
                                       '-conv <species_of_interest> (optional)')

parser.add_argument('--treefile', '-t',
                    help='Please specify a NWK tree file with internal node labels')
parser.add_argument('--alignmentfile', '-a',
                    help='Please specify the alignment file used to generate the tree file')
parser.add_argument('--convention', '-conv', default='NoSpecies',
                    help='(optional)Please specify a species of interest for which the convention will be added')

args = parser.parse_args()

tree_file = args.treefile
alignment_file = args.alignmentfile
convention_info = args.convention

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

# Read the alignment FASTA file and record values for all alignment positions
all_values = read_fasta(alignment_file)
alignment_length = len(all_values[all_values.keys()[0]])
print "The number of positions to be analyzed is: %d" % alignment_length

# Prepare convention information if specified
correspondences = {}
if convention_info == 'NoSpecies':
    results = open('BinaryF-table_%s.txt' % (tree_file.split('.')[0]), 'w')
    header = '\t' + '\t'.join(edges_all) + '\n'
    results.write(header)
else:  # If a species convention is requested
    results = open('BinaryF-table_%s_%s_convention.txt' % (tree_file.split('.')[0], convention_info), 'w')
    header = '\t\t' + '\t'.join(edges_all) + '\n'
    results.write(header)

    species_sequence = all_values[convention_info]
    j = 0
    for i in range(0, len(species_sequence)):
        if species_sequence[i] != '-':
            correspondences['Pos' + str(i + 1)] = species_sequence[i] + str(j + 1)
            j += 1
        else:
            correspondences['Pos' + str(i + 1)] = '-'
###

data_line_identical = '\t'.join(['0']*len(edges_all))  # Edges with identical values are assigned 0's

# Mapping substitutions
start = time.time()
for pos in range(0, alignment_length):  # Process one position at a time

    # Retrieve all values for the given position for all nodes and leaves
    pos_values = {}
    for sp in all_values.keys():
        pos_values[sp] = all_values[sp][pos].upper()

    # Decide if the given position is all identical across all nodes and leaves
    if len(list(set(pos_values.values()))) == 1:
        data_out = 'Pos%r' % (pos + 1) + '\t' + data_line_identical + '\n'
    # If not, compute substitutions
    else:
        mutated_edges = []

        for item in edges_all:
            node1 = item.split('*')[0]
            node2 = item.split('*')[1]
            if pos_values[node1] != pos_values[node2]:
                mutated_edges.append('1')   # Mutating edge

            else:
                mutated_edges.append('0')  # Non-mutating edge

        data_out = 'Pos%r' % (pos + 1) + '\t' + '\t'.join(mutated_edges) + '\n'

    if convention_info == 'NoSpecies':
        results.write(data_out)
    else:
        results.write(correspondences['Pos%r' % (pos + 1)] + '\t' + data_out)

end = time.time()

print "Computation took %r seconds." % (end - start)

results.close()

print "Analysis complete."
