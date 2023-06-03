#! /usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

from metadecoder.fasta_utility import read_fasta_file


def __init__(parameters):
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'Construct the mappings between sequence ID and cluster ID.',
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    parser.add_argument(
        'cluster', type = str, nargs = '+', metavar = 'cluster',
        help = 'The fasta formatted cluster files.\n The file name needs to match the format as follows:\n"program.cluster_id.fasta" e.g. "MetaDecoder.1.fasta".'
    )
    parser.add_argument(
        '-amber', '--amber', default = False, action = 'store_true',
        help = 'Output CAMI headers.\nDefault: False.'
    )
    parser.add_argument(
        '-s', '--sample', default = 'MetaDecoder', metavar = 'str',
        help = 'Sample ID used in CAMI headers.\nDefault: MetaDecoder.'
    )
    parser.add_argument(
        '-nh', '--no_header', default = False, action = 'store_true',
        help = 'Do not output headers.\nDefault: False.'
    )
    return parser.parse_args(parameters)


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])
    if parameters.amber:
        print('@Version:0.9.0', flush = True)
        print('@SampleID:{0}'.format(parameters.sample, ), flush = True)
        print('@@SEQUENCEID', 'BINID', sep = '\t', flush = True)
    elif not parameters.no_header:
        print('Sequence ID', 'Cluster ID', sep = '\t', flush = True)
    for cluster_file in parameters.cluster:
        cluster_id = os.path.splitext(os.path.basename(cluster_file))[0].split('.')[-1]
        for sequence_id, sequence in read_fasta_file(cluster_file):
            print(sequence_id, cluster_id, sep = '\t', flush = True)