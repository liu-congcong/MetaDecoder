#! /usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

from metadecoder.fasta_utility import read_fasta_file


def __init__(parameters):
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'Generate fasta formatted clusters on the basis of an assembly and a cluster file.',
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    parser.add_argument(
        '-a', '--assembly', type = str, required = True, metavar = 'str',
        help = 'The input fasta formatted assembly file.'
    )
    parser.add_argument(
        '-c', '--cluster', type = str, required = True, metavar = 'str',
        help = 'The input cluster file.'
    )
    parser.add_argument(
        '-o', '--output_prefix', type = str, required = True, metavar = 'str',
        help = 'The prefix of output files.'
    )
    return parser.parse_args(parameters)


def read_cluster_file(input_file):
    cluster_id2sequence_id = dict()
    open_file = open(input_file, 'r')
    open_file.readline()
    for line in open_file:
        sequence_id, cluster_id = line.strip('\n').split('\t')
        cluster_id2sequence_id.setdefault(cluster_id, list()).append(sequence_id)
    open_file.close()
    return cluster_id2sequence_id


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])
    sequence_id2sequence = dict((sequence_id, sequence) for sequence_id, sequence in read_fasta_file(parameters.assembly))
    cluster_id2sequence_id = read_cluster_file(parameters.cluster)
    for cluster_id, sequence_ids in cluster_id2sequence_id.items():
        open_file = open(parameters.output_prefix + '.' + cluster_id + '.fasta', 'w')
        for sequence_id in sequence_ids:
            open_file.write('>' + sequence_id + '\n')
            open_file.write(sequence_id2sequence[sequence_id] + '\n')
        open_file.close()