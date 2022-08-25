#!/usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from re import compile
from subprocess import PIPE, run


def __init__(parameters):
    parser = ArgumentParser(
        prog = sys.argv[0],
        formatter_class = RawTextHelpFormatter,
        description = 'Add gtdb annotations to tree.',
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    parser.add_argument(
        '-g', '--gtdb', type = str, required = True, metavar = 'file',
        help = 'The input gtdbtk summary file.'
    )
    parser.add_argument(
        '-t', '--tree', type = str, required = True, metavar = 'file',
        help = 'The input tree file.'
    )
    return parser.parse_args(parameters)


def read_tree_file(input_file):
    nodes = list()
    template = compile('[,(]([^:,(]+)')
    open_file = open(input_file, 'r')
    content = open_file.read()
    open_file.close()
    for node in template.finditer(content):
        if node:
            nodes.append(node.group(1))
    return nodes


def read_gtdb_file(input_file, cluster_ids):
    cluster_id2classification = dict()
    novel_species = list()
    open_file = open(input_file, 'r')
    open_file.readline()
    for line in open_file:
        lines = line.rstrip('\n').split('\t', maxsplit = 6)
        if lines[0] in cluster_ids:
            phylum = lines[1].split(';', maxsplit = 2)[1][3 : ]
            if not phylum:
                phylum = 'unknown'
            cluster_id2classification[lines[0]] = phylum
            if lines[5] == 'N/A' or float(lines[5]) < 95.0:
                novel_species.append(lines[0])
    open_file.close()
    return (cluster_id2classification, novel_species)


def get_colors(classifications):
    classification2color = dict()
    run_process = run(['Rscript', '-e', "library(scales);cat(hue_pal()({0}))".format(len(classifications))], stdout = PIPE)
    assert not run_process.returncode, "Run R failed!"
    for classification, color in zip(classifications, run_process.stdout.decode().split()):
        classification2color[classification] = color
    return classification2color


def write_itol_color_file(cluster_id2classification, classification2color, output_file):
    open_file = open(output_file, 'w')
    open_file.write('TREE_COLORS\nSEPARATOR TAB\nDATA\n')
    for cluster_id, classification in cluster_id2classification.items():
        classification = cluster_id2classification[cluster_id]
        open_file.write('\t'.join([cluster_id, 'range', classification2color[classification], classification]) + '\n')
    open_file.close()


def write_itol_binary_file(total_cluster_ids, novel_cluster_ids, output_file):
    open_file = open(output_file, 'w')
    open_file.write('DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tANI\nCOLOR\t#ff0000\nFIELD_SHAPES\t2\nFIELD_LABELS\t≤95\nDATA\n')
    for cluster_id in total_cluster_ids:
        binary_label = '-1'
        open_file.write('\t'.join([cluster_id, ('1' if cluster_id in novel_cluster_ids else '-1')]) + '\n')
    open_file.close()


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])
    cluster_ids = read_tree_file(parameters.tree)
    cluster_id2classification, novel_species = read_gtdb_file(parameters.gtdb, set(cluster_ids))
    classification2color = get_colors(sorted(set(cluster_id2classification.values())))
    output_prefix = os.path.basename(parameters.tree)
    write_itol_color_file(cluster_id2classification, classification2color, output_prefix + '.itol.color.txt')
    write_itol_binary_file(cluster_ids, set(novel_species), output_prefix + '.itol.binary.txt')