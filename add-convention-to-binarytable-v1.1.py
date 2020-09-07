#!/usr/bin/python
# coding=utf-8
#
# C.D. Dunn†, B.A. Akpınar, and V. Sharma†. 2019. An unusual amino acid substitution within the hummingbird
# cytochrome c oxidase alters a key proton-conducting channel. G3: Genes, Genomes, Genetics, 10:7, 2477-2485
# doi: 10.1534/g3.120.401312

# This script adds amino acid position [within the full length protein] correspondences of the alignment positions
# based on a user-defined species.

# The output file specifies amino acid position [within the full length protein] as AAi,
# and alignment position as Posi, for the i'th position.

# Script usage: $python add-convention-to-binarytable.py -b [binary_table] -a [alignment_file] -conv <species>

# Dependencies: Bio.SeqIO

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(prog='Add convention to binary table',
                                 usage='$python add-convention-to-binarytable-v1.1.py -b [binary_table]'
                                       '-a [alignment_file] -conv <species>')

parser.add_argument('--binary', '-b',
                    help='Please specify a binary table')
parser.add_argument('--alignmentfile', '-a',
                    help='Please specify an ancestral alignment file')
parser.add_argument('--convention', '-conv',
                    help='Please specify a species of interest for which the convention will be added')


args = parser.parse_args()

bin_table = args.binary
alignment_file = args.alignmentfile
species = args.convention

##########################################################################################
#################### Functions to be used in the script ##################################


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

##########################################################################################
##########################################################################################

# Open a new file to output results

try:
    results = open(bin_table.split('.')[0] + '_%s_convention.txt' % species, 'w')
except:
    results = open(bin_table + '_%s_convention.txt' % species, 'w')

# Read the FASTA file and record all values for each position for the species of interest
all_values = read_fasta(alignment_file)
species_sequence = all_values[species]
all_values.clear()

# Make a correspondence dictionary between amino acid positions and alignment positions for the given species
correspondences = {}
j = 0
for i in range(0, len(species_sequence)):
    if species_sequence[i] != '-':
        correspondences['Pos' + str(i + 1)] = species_sequence[i] + str(j + 1)
        j += 1
    else:
        correspondences['Pos' + str(i + 1)] = '-'

# Read the input binary table and add conventions to all alignment positions, based on the dictionary generated above
with open(bin_table) as f:
    for i, line in enumerate(f):
        if i == 0:
            header = 'AAi\tPosi' + line
            results.write(header)
        else:
            elems = line.split('\t')
            new_line = [correspondences[elems[0]]] + elems
            results.write('\t'.join(new_line))

results.close()

print "Analysis complete."
