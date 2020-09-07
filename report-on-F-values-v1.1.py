#!/usr/bin/python
# coding=utf-8
#
# C.D. Dunn†, B.A. Akpınar, and V. Sharma†. 2019. An unusual amino acid substitution within the hummingbird
# cytochrome c oxidase alters a key proton-conducting channel. G3: Genes, Genomes, Genetics, 10:7, 2477-2485
# doi: 10.1534/g3.120.401312

# This script reports node values of fluctuating edges.
# The script accepts two types of query files:
# (1) edge: A comma-separated text file containing [Edge,position] pairs. E.g. Node1*Node2,101
# (2) binary: A binaryF table, the direct output of "binary-table-by-edges.py"

# Script usage is: $python report-on-F-values.py -t [tree_file] -a [alignment_file] -q [query_file] -tp <edge/binary>
# Types for -tp (--type) parameter are 'edge' or 'binary'. Default is 'edge'.

# Dependencies: Bio.SeqIO and Bio.Phylo

import argparse
from Bio import Phylo
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(prog='report on F',
                                 usage='$python report-on-F-values-v.1.1.py -t [tree_file] -a [alignment_file] '
                                       '-q [query_file] -tp <edge/binary>')

parser.add_argument('--treefile', '-t',
                    help='Please specify a NWK tree file with internal node labels')
parser.add_argument('--alignmentfile', '-a',
                    help='Please specify the alignment file used to generate the tree file')
parser.add_argument('--queryfile', '-q',
                    help='Please specify a list of query file: (1) binary, or (2) edge')
parser.add_argument('--type', '-tp', default='edge',
                    help='Please specify the type of the query file: binary/edge')

args = parser.parse_args()

tree_file = args.treefile
alignment_file = args.alignmentfile
query_file = args.queryfile
query_type = args.type

##########################################################################################


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


# To retrieve clades by name

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names


##########################################################################################

# Read the NWK tree
my_tree = Phylo.read(tree_file, 'newick')
clades_dict = lookup_by_names(my_tree)

# Read the FASTA file and record amino acid values for each position
all_values = read_fasta(alignment_file)

if query_type == 'edge':  # If an edge*edge list is given as input

    try:
        results = open('%s_Fvalues.txt' % query_file.split('.')[0], 'w')
    except:
        results = open('%s_Fvalues.txt' % query_file, 'w')

    header = 'Edge name\tAncestral node\tDescendant node/species\tAlignment position\t' \
             'Ancestral character\tDescendant character\n'
    results.write(header)

    with open(query_file) as f:
        for line in f:
            # Define the parent and the child nodes
            if '*' in line:
                node1 = line.split(',')[0].split('*')[0]
                node2 = line.split(',')[0].split('*')[1]
                if clades_dict[node1] in clades_dict[node2]:
                    parent = node2
                    child = node1
                elif clades_dict[node2] in clades_dict[node1]:
                    parent = node1
                    child = node2
                else:
                    print "Are the input nodes connected: %s?" % (node1 + '*' + node2)
                    print "Please check your input."
                    sys.exit()

            # Define the alignment position
            position = line.strip('\n').split(',')[1]
            position = int(position.strip('\r'))

            parent_F = all_values[parent][position - 1].upper()
            child_F = all_values[child][position - 1].upper()

            data_line = [node1 + '*' + node2, parent, child, str(position), parent_F, child_F]
            results.write('\t'.join(data_line) + '\n')
            print 'Query edge: %s is processed.' % (node1 + '*' + node2)

else:  # If a binaryF table is given as input

    try:
        results = open('%s_Fluctuating_Edges.txt' % query_file.split('.')[0], 'w')
    except:
        results = open('%s_Fluctuating_Edges.txt' % query_file, 'w')

    header = 'Alignment position\tEdge name\tAncestral node\tDescendant node/species\tAncestral character' \
             '\tDescendant character\tTotal no. of fluctuations for the alignment position\n'
    results.write(header)

    edges_dict = {}
    pos_edge_dict = {}
    fluct_counter = {}
    with open(query_file) as f:
        for i, line in enumerate(f):
            elems = line.strip('\n').strip('\r').split('\t')
            if i == 0:  # First, get all edge*edge pairs from the header of binaryF table
                for j, elem in enumerate(elems):
                    edges_dict[j] = elem
            else:
                for j, elem in enumerate(elems):
                    if j > 0:
                        if int(elem) == 1:  # Collect information for all fluctuating edges that have 1's in binaryF
                            position_no = int(elems[0].strip('Pos'))
                            if position_no not in pos_edge_dict.keys():
                                pos_edge_dict[position_no] = [edges_dict[j]]
                                fluct_counter[position_no] = 1
                            else:
                                pos_edge_dict[position_no].append(edges_dict[j])
                                fluct_counter[position_no] += 1
        positions_sorted = sorted(pos_edge_dict.keys())
        for item_pos in positions_sorted:
            for item_edge in sorted(pos_edge_dict[item_pos]):  # Get values for the fluctuating edges for each position
                if '*' in item_edge:
                    node1 = item_edge.split('*')[0]
                    node2 = item_edge.split('*')[1]
                    if clades_dict[node1] in clades_dict[node2]:
                        parent = node2
                        child = node1
                    elif clades_dict[node2] in clades_dict[node1]:
                        parent = node1
                        child = node2
                    else:
                        print "Are the input nodes connected: %s?" % (node1 + '*' + node2)
                        print "Please check your input."
                        sys.exit()

                parent_F = all_values[parent][item_pos - 1].upper()
                child_F = all_values[child][item_pos - 1].upper()

                data_line = ['Pos' + str(item_pos), item_edge, parent, child, parent_F, child_F,
                             str(fluct_counter[item_pos])]
                results.write('\t'.join(data_line) + '\n')
                print 'Query edge: %s is processed.' % (node1 + '*' + node2)

print '\nAll queries are processed.'
