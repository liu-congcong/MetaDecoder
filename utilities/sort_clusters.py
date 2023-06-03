#! /usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter


def __init__(parameters):
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = '.',
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    parser.add_argument(
        'sample', type = str, metavar = 'id'
    )
    parser.add_argument(
        'program', type = str, choices = ('concoct', 'maxbin2', 'metabat2', 'metadecoder', 'vamb'), metavar = 'concoct|maxbin2|metabat2|metadecoder|vamb'
    )
    parser.add_argument(
        'cluster', type = str, nargs = '+', metavar = 'fasta'
    )
    return parser.parse_args(parameters)


def rename_concoct_cluster(sample, cluster):
    cluster = os.path.abspath(cluster)
    cluster_id = os.path.basename(cluster).split('.')[0]
    os.rename(cluster, sample + '.concoct.' + str(int(cluster_id) + 1) + '.fasta')
    return None


def rename_maxbin2_cluster(sample, cluster):
    cluster = os.path.abspath(cluster)
    cluster_id = os.path.basename(cluster).split('.')[2]
    os.rename(cluster, sample + '.maxbin2.' + str(int(cluster_id)) + '.fasta')
    return None


def rename_metabat2_cluster(sample, cluster):
    cluster = os.path.abspath(cluster)
    cluster_id = os.path.basename(cluster).split('.')[1]
    os.rename(cluster, sample + '.metabat2.' + cluster_id + '.fasta')
    return None


def rename_metadecoder_cluster(sample, cluster):
    cluster = os.path.abspath(cluster)
    cluster_id = os.path.basename(cluster).split('.')[1]
    os.rename(cluster, sample + '.metadecoder.' + cluster_id + '.fasta')
    return None


def rename_vamb_cluster(sample, cluster):
    cluster = os.path.abspath(cluster)
    cluster_id = os.path.basename(cluster).split('.')[0]
    os.rename(cluster, sample + '.vamb.' + cluster_id + '.fasta')
    return None


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])
    for cluster in parameters.cluster:
        if parameters.program == 'concoct':
            rename_concoct_cluster(parameters.sample, cluster)
        elif parameters.program == 'maxbin2':
            rename_maxbin2_cluster(parameters.sample, cluster)
        elif parameters.program == 'metabat2':
            rename_metabat2_cluster(parameters.sample, cluster)
        elif parameters.program == 'metadecoder':
            rename_metadecoder_cluster(parameters.sample, cluster)
        elif parameters.program == 'vamb':
            rename_vamb_cluster(parameters.sample, cluster)