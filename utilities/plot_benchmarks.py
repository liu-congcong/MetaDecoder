#! /usr/bin/env python3
import os
import sys
from argparse import ArgumentParser,RawTextHelpFormatter
from math import ceil, sqrt
import matplotlib.pyplot


def __init__(parameters):
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'Plot benchmarks assessed by Amber or CheckM.\n\nPlot Amber:\npython3 {0} --input AMBER_OUTPUT|CHECKM_OUTPUT --output OUPUT'.format(sys.argv[0]),
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    parser.add_argument(
        '-i', '--input', type = str, required = True, metavar = 'file|directory',
        help = 'This parameter is compatible with 2 types of input:\n1. TSV formatted CheckM\'s output [FILE].\n2. Amber\'s output [DIRECTORY].'
    )
    parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = 'prefix'
    )
    parser.add_argument(
        '-p', '--precision', default = 0.95, type = float, required = False, metavar = 'float',
        help = 'Threshold of precision to plot (0.0 - 1.0).\nDefault: 0.95.'
    )
    return parser.parse_args(parameters)


def read_amber_metrics_file(file, precision, recall2count):
    # sample_id	Bin ID	Most abundant genome	Purity (bp)	Completeness (bp)	Bin size (bp) #
    open_file = open(file, 'r')
    open_file.readline() # Remove header line. #
    for line in open_file:
        lines = line.rstrip('\n').split('\t', maxsplit = 5)
        if float(lines[3]) >= precision:
            recall = float(lines[4])
            for recall_threshold in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95):
                if recall >= recall_threshold:
                    recall2count[recall_threshold] = recall2count.get(recall_threshold, 0) + 1
                else:
                    break
    open_file.close()
    return None


def read_amber_directory(amber, precision, sample2program2recall2count):
    genome = os.path.join(amber, 'genome')
    if os.access(genome, os.F_OK):
        sample2program2recall2count['AMBER'] = dict()
        program2recall2count = sample2program2recall2count['AMBER']
        for entry in os.scandir(genome):
            if entry.is_dir() and entry.name.lower() != 'gold standard':
                program2recall2count[entry.name] = dict()
                metrics = os.path.join(entry.path, 'metrics_per_bin.tsv')
                assert os.access(metrics, os.R_OK), 'Can not find "metrics_per_bin.tsv" in {0}.'.format(entry.path)
                read_amber_metrics_file(metrics, precision, program2recall2count[entry.name])
    else:
        for entry in os.scandir(amber):
            if entry.is_file() and entry.name.endswith('.amber_metrics.tsv'):
                sample, program, _ = entry.name.split('.', maxsplit = 2)
                sample2program2recall2count.setdefault(sample, dict()).setdefault(program, dict())
                read_amber_metrics_file(entry.path, precision, sample2program2recall2count[sample][program])
    return None


def read_checkm_file(checkm, precision, sample2program2recall2count):
    open_file = open(checkm, 'r')
    open_file.readline() # Remove header line. #
    # Bin Id	Marker lineage	# genomes	# markers	# marker sets	0	1	2	3	4	5+	Completeness	Contamination	Strain heterogeneity #
    # Bin Id: sample.program.cluster #
    for line in open_file:
        lines = line.strip('\n').split('\t')
        sample, program, _ = lines[0].split('.', maxsplit = 2)
        if 1 - float(lines[12]) / 100 >= precision:
            recall = float(lines[11]) / 100
            for recall_threshold in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95):
                if recall >= recall_threshold:
                    sample2program2recall2count[sample][program][recall_threshold] = sample2program2recall2count.setdefault(sample, dict()).setdefault(program, dict()).get(recall_threshold, 0) + 1
                else:
                    break
    open_file.close()
    return None


def plot(sample2program2recall2count, precision, output):
    samples = len(sample2program2recall2count)
    columns = ceil(sqrt(samples))
    rows = ceil(samples / columns)
    figure, subplots = matplotlib.pyplot.subplots(nrows = rows, ncols = columns, figsize = (3.6, 3.6) if samples == 1 else (columns * 1.8, rows * 1.8), constrained_layout = True, squeeze = False)
    subplots = subplots.flatten()
    for subplot in subplots[samples : ]:
        subplot.remove()
    recalls = list()
    counts = list()
    unique_programs = dict()
    for subplot, (sample, program2recall2count) in zip(subplots, sorted(sample2program2recall2count.items())):
        programs = sorted(program2recall2count.keys())
        for program in programs:
            unique_programs[program] = 0
            for recall, count in program2recall2count[program].items():
                recalls.append(recall)
                counts.append(count)
            subplot.plot(recalls, counts, '*-', markersize = 2, alpha = 0.8, linewidth = 1, label = program)
            subplot.set_title(sample, fontdict = {'fontsize': 10})
            subplot.set_xlim(0.97, 0.48)
            subplot.set_xticks((0.9, 0.8, 0.7, 0.6, 0.5))
            subplot.tick_params(labelsize = 8)
            recalls.clear()
            counts.clear()
    figure.legend(sorted(unique_programs.keys()), title = 'Program', bbox_to_anchor = (1, 0.5), loc = 'center left', fontsize = 8, borderaxespad = 0.0)
    figure.set_constrained_layout_pads(w_pad = 0.05, h_pad = 0.05, wspace = 0.05, hspace = 0.05)
    figure.supxlabel('Recall (Precision>={0})'.format(precision), fontsize = 9)
    figure.supylabel('Count', fontsize = 9)
    figure.savefig(output, dpi = 300, bbox_inches = 'tight')
    return None


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])
    assert parameters.precision >= 0.0 and parameters.precision <= 1.0, 'The value of "--precision" must be a float from 0.0 to 1.0.'
    counts = list()
    sample2program2recall2count = dict()
    if os.path.isdir(parameters.input):
        read_amber_directory(parameters.input, parameters.precision, sample2program2recall2count)
    else:
        read_checkm_file(parameters.input, parameters.precision, sample2program2recall2count)
    open_file = open(parameters.output + '.tsv', 'w')
    open_file.write('Sample\tProgram\tPrecision\tRecall\tCount\n')
    print('Sample', 'Program', 'Precision', 'Recall≥0.95', 'Recall≥0.90', 'Recall≥0.85', 'Recall≥0.80', 'Recall≥0.75', 'Recall≥0.70', 'Recall≥0.65', 'Recall≥0.60', 'Recall≥0.55', 'Recall≥0.50', sep = '\t', flush = True)
    for sample, program2recall2count in sample2program2recall2count.items():
        for program, recall2count in program2recall2count.items():
            for recall in (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95):
                open_file.write('\t'.join([sample, program, str(parameters.precision), str(recall), str(recall2count.get(recall, 0))]) + '\n')
            for recall in (0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5):
                counts.append(str(recall2count.get(recall, 0)))
            print(sample, program, '≥{0}'.format(parameters.precision), '\t'.join(counts), sep = '\t', flush = True)
            counts.clear()
    open_file.close()
    plot(sample2program2recall2count, parameters.precision, parameters.output + '.pdf')