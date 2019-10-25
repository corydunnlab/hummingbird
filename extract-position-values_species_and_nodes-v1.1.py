#!/usr/bin/python
# coding=utf-8
#
# C.D. Dunn†, B.A. Akpınar, and V. Sharma†. An unusual amino acid substitution within the hummingbird
# cytochrome c oxidase alters a key proton-conducting channel. bioRxiv. doi: 10.1101/610915.

# This script lists all possible values for each position in a given alignment.
# If a species is given, a species convention is also added to the output file.

# Script usage: $python extract-position-values.py -a [alignment_file] -conv <species_convention> (optional)

# Dependencies: Bio.SeqIO

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(prog='Extract position values',
                                 usage='$python extract-position-values-v1.1.py -a [alignment_file] '
                                       '-conv <species_convention> (optional)')


parser.add_argument('--alignmentfile', '-a',
                    help='Please specify an alignment file used to generate the tree file')
parser.add_argument('--convention', '-conv', default='nospecies',
                    help='Please specify a species to add a species convention to the alignment position')

args = parser.parse_args()

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

# Open a new file to output results
try:
    results = open(alignment_file.split('.')[0] + '_PosValues.txt', 'w')
except:
    results = open(alignment_file + '_PosValues.txt', 'w')

# Read the alignment file
all_values = read_fasta(alignment_file)
alignment_length = len(all_values[all_values.keys()[0]])
print 'Total number of positions in the alignment is: %r\n' % alignment_length

if species != 'nospecies':  # If a species convention is requested
    species_sequence = all_values[species]
    correspondences = {}
    j = 0
    for i in range(0, len(species_sequence)):
        if species_sequence[i] != '-':
            correspondences['Pos' + str(i + 1)] = species_sequence[i] + str(j + 1)
            j += 1
        else:
            correspondences['Pos' + str(i + 1)] = '-'

    header = 'Alignment position\tValues\t%s convention\n' % species
    results.write(header)

    for pos in range(0, alignment_length):  # Process one position at a time
        print 'Processing position: %r' % (pos + 1)

        alignment_pos = 'Pos' + str(pos + 1)

        # Retrieve all values for the given position for all nodes and leaves
        pos_values = []
        for sp in all_values.keys():
            pos_values.append(all_values[sp][pos].upper())

        data_line = [alignment_pos, ', '.join(set(list(pos_values))), correspondences[alignment_pos]]
        results.write('\t'.join(data_line) + '\n')

else:  # If a species convention is NOT requested
    header = 'Alignment position\tValues\n'
    results.write(header)

    for pos in range(0, alignment_length):  # Process one position at a time
        print 'Processing position: %r' % (pos + 1)

        alignment_pos = 'Pos' + str(pos + 1)

        # Retrieve all values for the given position for all nodes and leaves
        pos_values = []
        for sp in all_values.keys():
            pos_values.append(all_values[sp][pos].upper())

        data_line = [alignment_pos, ', '.join(set(list(pos_values)))]
        results.write('\t'.join(data_line) + '\n')

results.close()

print "Analysis complete."
