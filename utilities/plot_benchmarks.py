#! /usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from math import ceil, cos, pi, pow, sin, sqrt

import matplotlib.pyplot
import numpy
from matplotlib.ticker import AutoMinorLocator


class HCL:

    def __init__(self, min_hue = 15.0, max_hue = 375.0, chroma = 100.0, luminance = 65.0):
        assert luminance > 0.0, '"luminance" should be a positive float.'

        self.__min_hue = min_hue
        self.__max_hue = max_hue
        self.__chroma = chroma
        self.__luminance = luminance
        self.__white_x = 95.047
        self.__white_y = 100.0
        self.__white_z = 108.883
        self.__gamma = 2.4
        self.__eps = 216.0 / 24389.0 # 0.008856452 #
        self.__kappa = 24389.0 / 27.0 # 903.2963 #
        self.__xyz2rgb_coefficient = ((3.2404542, -1.5371385, -0.4985314), (-0.9692660, 1.8760108, 0.0415560), (0.0556434, -0.2040259, 1.0572252))
        return None


    def __xyz2rgb(self, x, y, z, coefficient_index):
        coefficient = self.__xyz2rgb_coefficient[coefficient_index]
        rgb = x * coefficient[0] + y * coefficient[1] + z * coefficient[2]
        rgb = 255.0 * (1.055 * pow(rgb, (1.0 / self.__gamma)) - 0.055 if rgb > 0.0031308 else 12.92 * rgb)
        return min(max(round(rgb), 0), 255)


    def main(self, colors):
        assert colors > 0, '"colors" should be a positive integer.'

        white_15x_3y_1z = self.__white_x + 15 * self.__white_y + 3 * self.__white_z
        u = 4.0 * self.__white_x / white_15x_3y_1z
        v = 9.0 * self.__white_y / white_15x_3y_1z
        Y = pow((self.__luminance + 16.0) / 116.0, 3) if (self.__luminance > self.__eps * self.__kappa) else self.__luminance / self.__kappa
        B = - 5.0 * Y
        if not (self.__max_hue - self.__min_hue) % 360:
            self.__max_hue -= 360.0 / colors
        hue_step = (self.__max_hue - self.__min_hue) / (colors - 1.0) if colors > 1 else 0.0
        for color_index in range(colors):
            hue = min(max(self.__min_hue + (color_index + 1.0 if color_index else 0.0) * hue_step, 0.0), 360.0) * pi / 180.0
            A = 1.0 / 3.0 * (52.0 * self.__luminance / (self.__chroma * cos(hue) + 13.0 * self.__luminance * u) - 1.0)
            X = (Y * (39.0 * self.__luminance / (self.__chroma * sin(hue) + 13.0 * self.__luminance * v) - 5.0) - B) / (A + 1.0 / 3.0)
            Z = X * A + B
            yield (self.__xyz2rgb(X, Y, Z, 0) / 255.0, self.__xyz2rgb(X, Y, Z, 1) / 255.0, self.__xyz2rgb(X, Y, Z, 2) / 255.0)
        return self


def __init__(parameters):
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'Plot benchmarkss assessed by Amber or CheckM.\n\npython3 {0} -i AMBER_OUTPUT|CHECKM_OUTPUT [-o OUPUT]'.format(sys.argv[0]),
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    parser.add_argument(
        '-i', '--input', type = str, required = True, metavar = 'file|directory',
        help = '''Compatible with 2 types of inputs:
-------------------------------------------------------------------------------------------
1 - TSV formatted CheckM\'s output [FILE].
    The format of the first column must be "Dataset.Program.ClusterID".
    You can define the name of each cluster file as "Dataset.Program.ClusterID.fasta",
    and put all files in the "FASTAS" folder, then run:
    "checkm lineage_wf --tab_table --extension fasta --file CHECKM.TSV FASTAS CHECKM.TEMP"
    to get the result "CHECKM.TSV".
-------------------------------------------------------------------------------------------
2 - Amber\'s output [DIRECTORY].
-------------------------------------------------------------------------------------------'''
    )
    parser.add_argument(
        '-o', '--output', type = str, required = False, metavar = 'prefix',
        help = 'Prefix of output files.'
    )
    parser.add_argument(
        '-t', '--threshold', nargs = '+', default = [0.90, 0.95], type = float, required = False, metavar = 'float',
        help = 'Thresholds of precision or 1.0 - contamination to plot (0.0 - 1.0).\nDefault: 0.90 & 0.95.'
    )
    return parser.parse_args(parameters)


def read_amber_metrics_file(input_file, precisions):
    benchmark = numpy.zeros(shape = (len(precisions), 10), dtype = numpy.int64)
    recall = numpy.array([0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95], dtype = numpy.float64)
    open_file = open(input_file, 'r')
    open_file.readline() # Remove header line. #
    # sample_id	Bin ID	Most abundant genome	Purity (bp)	Completeness (bp)	Bin size (bp) #
    for line in open_file:
        lines = line.rstrip('\n').split('\t', maxsplit = 5)
        for index, precision in enumerate(precisions):
            if float(lines[3]) >= precision:
                benchmark[index, float(lines[4]) >= recall] += 1
    open_file.close()
    return benchmark


def read_amber_directory(amber, precisions):
    unique_programs = dict()
    sample2program2benchmark = dict()
    genome = os.path.join(amber, 'genome')
    if os.access(genome, os.F_OK):
        program2benchmark = sample2program2benchmark.setdefault('Dataset', dict())
        for entry in os.scandir(genome):
            if entry.is_dir() and entry.name.lower() != 'gold standard':
                metrics = os.path.join(entry.path, 'metrics_per_bin.tsv')
                assert os.access(metrics, os.R_OK), 'Can not find "metrics_per_bin.tsv" in {0}.'.format(entry.path)
                unique_programs[entry.name] = 0
                program2benchmark[entry.name] = read_amber_metrics_file(metrics, precisions)
    else: # only for debug #
        for entry in os.scandir(amber):
            if entry.is_file() and entry.name.endswith('.amber_metrics.tsv'):
                sample, program, _ = entry.name.split('.', maxsplit = 2)
                unique_programs[program] = 0
                sample2program2benchmark.setdefault(sample, dict())[program] = read_amber_metrics_file(entry.path, precisions)
    return (sample2program2benchmark, sorted(unique_programs.keys()))


def read_checkm_file(input_file, contaminations):
    unique_programs = dict()
    sample2program2benchmark = dict()
    completeness = numpy.array([50, 55, 60, 65, 70, 75, 80, 85, 90, 95], dtype = numpy.float64)
    open_file = open(input_file, 'r')
    open_file.readline() # Remove header line. #
    # Bin Id	Marker lineage	# genomes	# markers	# marker sets	0	1	2	3	4	5+	Completeness	Contamination	Strain heterogeneity #
    # Bin Id: dataset.program.cluster #
    for line in open_file:
        lines = line.strip('\n').split('\t')
        try:
            sample, program, _ = lines[0].split('.', maxsplit = 2)
        except Exception:
            print(
                'The format of the first column must be "Dataset.Program.ClusterID".',
                'You can name each cluster file as "Dataset.Program.ClusterID.fasta"',
                'and put all files in the "FASTAS" folder, then run:',
                '"checkm lineage_wf --tab_table --extension fasta --file CHECKM.TSV FASTAS CHECKM.TEMP"',
                'to get the result "CHECKM.TSV".',
                sep = '\n'
            )
            sys.exit()
        unique_programs[program] = 0
        benchmark = sample2program2benchmark.setdefault(sample, dict()).setdefault(program, numpy.zeros(shape = (len(contaminations), 10), dtype = numpy.int64))
        for index, contamination in enumerate(contaminations):
            if float(lines[12]) <= contamination * 100.0:
                benchmark[index, float(lines[11]) >= completeness] += 1
    open_file.close()
    return (sample2program2benchmark, sorted(unique_programs.keys()))


def write_output(sample2program2benchmark, benchmark_type, thresholds, output_file):
    open_file = open(output_file, 'w')
    open_file.write('\t'.join(['Sample', 'Program', benchmark_type[0], benchmark_type[1], 'Count']) + '\n')
    print('Sample', 'Program', benchmark_type[0], benchmark_type[1] + '≥0.90', benchmark_type[1] + '≥0.80', benchmark_type[1] + '≥0.70', benchmark_type[1] + '≥0.60', benchmark_type[1] + '≥0.50', sep = '\t', flush = True)
    for sample, program2benchmark in sample2program2benchmark.items():
        for program, benchmark in program2benchmark.items():
            for threshold, row_benchmark in zip(thresholds, benchmark.astype(numpy.str_)):
                threshold = str(round(threshold, 2))
                for column_index, column_benchmark in enumerate(row_benchmark):
                    open_file.write('\t'.join([sample, program, threshold, str(round(0.05 * column_index + 0.5, 2)), column_benchmark]) + '\n')
                print('\t'.join([sample, program, threshold] + row_benchmark[[8, 6, 4, 2, 0]].tolist()), flush = True)
    open_file.close()
    return None


def define_colors(unique_programs):
    hcl = HCL()
    program2color = dict()
    for index, color in enumerate(hcl.main(len(unique_programs))):
        program2color[unique_programs[index]] = color
    return program2color


def plot(sample2program2benchmark, benchmark_index, x_label, program2color, output_file):
    x = numpy.array([0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95], dtype = numpy.float64)
    samples = len(sample2program2benchmark)
    columns = ceil(sqrt(samples))
    rows = ceil(samples / columns)
    figure, subplots = matplotlib.pyplot.subplots(nrows = rows, ncols = columns, figsize = (3.6, 3.6) if samples == 1 else (columns * 2.0, rows * 2.0), constrained_layout = True, squeeze = False)
    subplots = subplots.flatten()
    for subplot in subplots[samples : ]:
        subplot.remove()
    for subplot, (sample, program2benchmark) in zip(subplots, sorted(sample2program2benchmark.items())):
        for program, benchmark in program2benchmark.items():
            subplot.plot(x, benchmark[benchmark_index], '.-', markersize = 4, alpha = 1.0, linewidth = 0.6, label = program, color = program2color[program])
            subplot.set_title(sample, fontdict = {'fontsize': 10})
            subplot.set_xlim(0.96, 0.49)
            subplot.set_xticks((0.9, 0.8, 0.7, 0.6, 0.5))
            subplot.tick_params(labelsize = 8)
            subplot.xaxis.set_minor_locator(AutoMinorLocator(2))
            subplot.yaxis.set_minor_locator(AutoMinorLocator(2))
            subplot.grid(True, which = 'both', linestyle = '--', linewidth = 0.5)
    figure.legend(sorted(program2color.keys()), title = 'Program', bbox_to_anchor = (1, 0.5), loc = 'center left', fontsize = 8, borderaxespad = 0.0)
    figure.set_constrained_layout_pads(w_pad = 0.05, h_pad = 0.05, wspace = 0.05, hspace = 0.05)
    figure.supxlabel(x_label, fontsize = 9)
    figure.supylabel('Number of clusters', fontsize = 9)
    figure.savefig(output_file, dpi = 300, bbox_inches = 'tight')
    return None


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])
    if not parameters.output:
        parameters.output = parameters.input
    for threshold in parameters.threshold:
        assert threshold >= 0.0 and threshold <= 1.0, 'Each element of "--threshold" must be a float from 0.0 to 1.0.'

    if os.path.isdir(parameters.input):
        sample2program2benchmark, unique_programs = read_amber_directory(parameters.input, parameters.threshold)
        benchmark_type = ('Precision', 'Recall')
    else:
        parameters.threshold = [1.0 - threshold for threshold in parameters.threshold]
        sample2program2benchmark, unique_programs = read_checkm_file(parameters.input, parameters.threshold)
        benchmark_type = ('Contamination', 'Completeness')

    write_output(sample2program2benchmark, benchmark_type, parameters.threshold, parameters.output + '.tsv')

    program2color = define_colors(unique_programs)
    for benchmark_index, threshold in enumerate(parameters.threshold):
        plot(sample2program2benchmark, benchmark_index, benchmark_type[1], program2color, '{0} ({1}{2:.2f}).pdf'.format(parameters.output, benchmark_type[0], threshold))
