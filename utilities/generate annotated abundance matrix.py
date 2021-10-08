#!/usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

import numpy


def __init__(parameters):
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'Generate annotated abundance matrix.',
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    parser.add_argument(
        '-a', '--abundance', type = str, nargs = '+', required = True, metavar = 'file',
        help = 'The profiles of clusters reported by CheckM profile.'
    )
    parser.add_argument(
        '-g', '--gtdb', type = str, required = True, metavar = 'file',
        help = 'The annotitons of clusters reported by GTDB.'
    )
    parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = 'file'
    )
    return parser.parse_args(parameters)


def read_abundance_file(cluster_id2abundance, input_file, file_index):
    open_file = open(input_file, 'r')
    open_file.readline()
    for line in open_file:
        lines = line.rstrip('\n').rsplit('\t')
        if lines[0] != 'unbinned':
            cluster_id2abundance[lines[0]] = (file_index, float(lines[5]))
    open_file.close()
    return cluster_id2abundance


def read_gtdb_file(input_file):
    k2p2c2o2f2g2s2cluster_ids = dict()
    open_file = open(input_file, 'r')
    open_file.readline()
    for line in open_file:
        cluster_id, classifications, _ = line.rstrip('\n').split('\t', maxsplit = 2)
        k, p, c, o, f, g, s = classifications.split(';')
        k2p2c2o2f2g2s2cluster_ids.setdefault(k[3 : ], dict()).setdefault(p[3 : ], dict()).setdefault(c[3 : ], dict()).setdefault(o[3 : ], dict()).setdefault(f[3 : ], dict()).setdefault(g[3 : ], dict()).setdefault(s[3 : ], list()).append(cluster_id)
    open_file.close()
    return k2p2c2o2f2g2s2cluster_ids


def main(k2p2c2o2f2g2s2cluster_ids, cluster_id2abundance, column_names, output_file):
    columns = len(column_names)
    classification_k_abundances = numpy.zeros(shape = columns, dtype = numpy.float64)
    classification_p_abundances = numpy.zeros(shape = columns, dtype = numpy.float64)
    classification_c_abundances = numpy.zeros(shape = columns, dtype = numpy.float64)
    classification_o_abundances = numpy.zeros(shape = columns, dtype = numpy.float64)
    classification_f_abundances = numpy.zeros(shape = columns, dtype = numpy.float64)
    classification_g_abundances = numpy.zeros(shape = columns, dtype = numpy.float64)
    classification_s_abundances = numpy.zeros(shape = columns, dtype = numpy.float64)
    open_file = open(output_file, 'w')
    open_file.write('Classification\t' + '\t'.join(column_names) + '\n')
    for k, p2c2o2f2g2s2cluster_ids in k2p2c2o2f2g2s2cluster_ids.items():
        for p, c2o2f2g2s2cluster_ids in p2c2o2f2g2s2cluster_ids.items():
            for c, o2f2g2s2cluster_ids in c2o2f2g2s2cluster_ids.items():
                for o, f2g2s2cluster_ids in o2f2g2s2cluster_ids.items():
                    for f, g2s2cluster_ids in f2g2s2cluster_ids.items():
                        for g, s2cluster_ids in g2s2cluster_ids.items():
                            for s, cluster_ids in s2cluster_ids.items():
                                for cluster_id in cluster_ids:
                                    abundance_index, abundance = cluster_id2abundance[cluster_id]
                                    classification_s_abundances[abundance_index] += abundance
                                if s != '':
                                    open_file.write('|'.join([k, p, c, o, f, g, s]) + '\t' + '\t'.join(classification_s_abundances.astype(numpy.str_).tolist()) + '\n')
                                classification_g_abundances += classification_s_abundances
                                classification_s_abundances *= 0.0
                            if g != '':
                                open_file.write('|'.join([k, p, c, o, f, g]) + '\t' + '\t'.join(classification_g_abundances.astype(numpy.str_).tolist()) + '\n')
                            classification_f_abundances += classification_g_abundances
                            classification_g_abundances *= 0.0
                        if f != '':
                            open_file.write('|'.join([k, p, c, o, f]) + '\t' + '\t'.join(classification_f_abundances.astype(numpy.str_).tolist()) + '\n')
                        classification_o_abundances += classification_f_abundances
                        classification_f_abundances *= 0.0
                    if o != '':
                        open_file.write('|'.join([k, p, c, o]) + '\t' + '\t'.join(classification_o_abundances.astype(numpy.str_).tolist()) + '\n')
                    classification_c_abundances += classification_o_abundances
                    classification_o_abundances *= 0.0
                if c != '':
                    open_file.write('|'.join([k, p, c]) + '\t' + '\t'.join(classification_c_abundances.astype(numpy.str_).tolist()) + '\n')
                classification_p_abundances += classification_c_abundances
                classification_c_abundances *= 0.0
            if p != '':
                open_file.write('|'.join([k, p]) + '\t' + '\t'.join(classification_p_abundances.astype(numpy.str_).tolist()) + '\n')
            classification_k_abundances += classification_p_abundances
            classification_p_abundances *= 0.0
        if k != '':
            open_file.write(k + '\t' + '\t'.join(classification_k_abundances.astype(numpy.str_).tolist()) + '\n')
        classification_k_abundances *= 0.0
    open_file.close()
    return None


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])
    cluster_id2abundance = dict()
    for abundance_index, abundance_file in enumerate(parameters.abundance):
        read_abundance_file(cluster_id2abundance, abundance_file, abundance_index)
    k2p2c2o2f2g2s2cluster_ids = read_gtdb_file(parameters.gtdb)
    main(k2p2c2o2f2g2s2cluster_ids, cluster_id2abundance, parameters.abundance, parameters.output)