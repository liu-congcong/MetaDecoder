import os
from argparse import ArgumentParser, RawTextHelpFormatter

import numpy
from .detect_permission import is_readable, is_writeable
from . import metadecoder_coverage
from . import metadecoder_seed
from . import metadecoder_cluster


def __init__():
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'The MetaDecoder algorithm.',
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )

    subparser = parser.add_subparsers(title = 'program', dest = 'program', required = True, help = 'Run metadecoder coverage|seed|cluster --help for help.')

    coverage_parser = subparser.add_parser('coverage', formatter_class = RawTextHelpFormatter, help = 'Calculate sequence coverages from BAM files.\nUsage: metadecoder coverage -s sample1.bam ... sampleN.bam -o output.coverage.')
    coverage_parser.add_argument(
        '-b', '--bam', type = str, nargs = '+', required = True, metavar = 'str',
        help = 'The bam formatted files with headers "@SQ".'
    )
    coverage_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = 'str',
        help = 'The output coverage file.'
    )
    coverage_parser.add_argument(
        '--bin_size', default = 2000, type = int, required = False, metavar = 'int',
        help = 'The size (bp) of each bin.\nThe value should be a positive integer.\nDefault: 2000.'
    )
    coverage_parser.add_argument(
        '--mapq', default = 0, type = int, required = False, metavar = 'int',
        help = 'The MAPQ threshold.\nThe value should be a non-negative integer.\nDefault: 0.'
    )
    coverage_parser.add_argument(
        '--aligned', default = 0.95, type = float, required = False, metavar = 'float',
        help = 'The minimum aligned ratio calculted from CIGAR.\nThe value should be from 0 to 1.\nDefault: 0.95.'
    )
    coverage_parser.add_argument(
        '--threads', default = 10, type = int, required = False, metavar = 'int',
        help = 'The number of threads.\nThe value should be a positive integer.\nDefault: 10.'
    )

    seed_parser = subparser.add_parser('seed', formatter_class = RawTextHelpFormatter, help = 'Map the marker genes (HMMs) to sequences.\nUsage: metadecoder seed -f input.fasta -o output.seed.')
    seed_parser.add_argument(
        '-f', '--fasta', default = '', type = str, required = True, metavar = 'str',
        help = 'The fasta formatted assembly file.'
    )
    seed_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = 'str',
        help = 'The output seed file.'
    )
    seed_parser.add_argument(
        '--threads', default = 20, type = int, required = False, metavar = 'int',
        help = 'The number of threads.\nThe value should be a positive integer.\nDefault: 20.'
    )
    seed_parser.add_argument(
        '--coverage', default = 0.5, type = float, required = False, metavar = 'float',
        help = 'The min coverage of domain.\nThe value should be from 0 to 1.\nDefault: 0.5.'
    )
    seed_parser.add_argument(
        '--accuracy', default = 0.6, type = float, required = False, metavar = 'float',
        help = 'The min accuracy.\nThe value should be from 0 to 1.\nDefault: 0.6.',
    )

    cluster_parser = subparser.add_parser('cluster', formatter_class = RawTextHelpFormatter, help = 'Running the MetaDecoder to cluster sequences.\nUsage: metadecoder cluster -f input.fasta -c input.coverage -s input.seed -o output.metadecoder')
    cluster_parser.add_argument(
        '-f', '--fasta', type = str, required = True, metavar = 'str',
        help = 'The fasta formatted assembly file.'
    )
    cluster_parser.add_argument(
        '-c', '--coverage', type = str, nargs = '+', required = True, metavar = 'str',
        help = 'The coverage files generated by "metadecoder coverage".'
    )
    cluster_parser.add_argument(
        '-s', '--seed', type = str, required = True, metavar = 'str',
        help = 'The seed file generated by "metadecoder seed".'
    )
    cluster_parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = 'str',
        help = 'The prefix of output files.'
    )
    cluster_parser.add_argument(
        '--output_unclustered_sequences', default = False, required = False, action = 'store_true',
        help = 'Generate "*.unclustered.fasta" file for all unclustered sequences.\nDefault: False.'
    )
    cluster_parser.add_argument(
        '--clustering_probability', default = 0.0, type = float, required = False, metavar = 'float',
        help = 'The minimum clustering probability.\nThe value should be from 0.0 to 1.0\nDefault: 0.0.'
    )
    cluster_parser.add_argument(
        '--weight', default = 1.0, type = float, required = False, metavar = 'float',
        help = 'The scale factor of the weight in regularization term.\nThe coverage probabilistic model uses "--weight" / "number of sequencing samples" as the weight of priors.\nThe value should be a non-negative float.\nDefault: 1.'
    )
    cluster_parser.add_argument(
        '--kmer', default = 4, type = int, required = False, metavar = 'int',
        help = 'The length of kmer.\nThe value should be a positive integer.\nDefault: 4.'
    )
    cluster_parser.add_argument(
        '--min_cluster_size', default = 200000, type = int, required = False, metavar = 'int',
        help = 'The minimum size of a cluster to output.\nThe value should be a positive integer.\nDefault: 200000.'
    )
    cluster_parser.add_argument(
        '--min_dpgmm_size', default = 500000, type = int, required = False, metavar = 'int',
        help = 'The minimum size of a DPGMM cluster.\nDPGMM stops shrinking the number of clusters, until the size of the smallest cluster is greater than or equal to this value.\nThe value should be a positive integer.\nDefault: 500000.'
    )
    cluster_parser.add_argument(
        '--max_dpgmm_distance', default = 0.04, type = float, required = False, metavar = 'float',
        help = 'Sequences in a DPGMM cluster will be removed, if the average Euclidean distance of all paired kmer frequencies is greater than this value.\nWe recommend that you do not change this value.\nThe value should be a positive float.\nDefault: 0.04.'
    )
    cluster_parser.add_argument(
        '--min_sequence_length', default = 2500, type = int, required = False, metavar = 'int',
        help = 'Sequences with the length being greater than or equal to this value can be involved in the MetaDecoder algorithm.\nThe value should be a integer greater than 2000.\nDefault: 2500.'
    )
    cluster_parser.add_argument(
        '--outlier', default = 0.1, type = float, required = False, metavar = 'float',
        help = 'The proportion of outliers that will be removed by the isolation forest.\nThe value should be from 0 to 0.5.\nDefault: 0.1.'
    )
    cluster_parser.add_argument(
        '--sampling_length1', default = (1500, 2500), nargs = 2, type = int, required = False, metavar = 'int',
        help = 'The range of the length of sub-seeds sampled from each seed used to train the classifier.\nThe value should be two positive integers from low to high.\nDefault: 1500 2500.'
    )
    cluster_parser.add_argument(
        '--sampling_number1', default = 20, type = int, required = False, metavar = 'int',
        help = 'The number of sub-seeds sampled from each seed used to train the classifier.\nThe value should be a positive integer.\nDefault: 20.'
    )
    cluster_parser.add_argument(
        '--sampling_length2', default = (3000, 4000), type = int, nargs = 2, required = False, metavar = 'int',
        help = 'The range of the length of sub-seeds sampled from each extended seed used to train the classifier.\nThe value should be two positive integers from low to high.\nDefault: 3000 4000.'
    )
    cluster_parser.add_argument(
        '--sampling_number2', default = 300, type = int, required = False, metavar = 'int',
        help = 'The number of sub-seeds sampled from each extended seed used to train the classifier.\nThe value should be a positive integer.\nDefault: 300.'
    )
    cluster_parser.add_argument(
        '--random_number', default = 0, type = int, required = False, metavar = 'int',
        help = 'The random number generator seeded by the given integer\nDefault: 0.'
    )
    cluster_parser.add_argument(
        '--disable_gpu', default = False, required = False, action = 'store_true',
        help = 'Disable the GPU if it does not have enough memory to accommodate the data.\nDefault: False.'
    )
    cluster_parser.add_argument(
        '--no_clusters', default = False, required = False, action = 'store_true',
        help = 'Do not output fasta format cluster files.\nDefault: False.'
    )

    parser.add_argument(
        '-v', '--version', action = 'version', version = '%(prog)s 1.2.1'
    )
    return parser.parse_args()


def parse_parameters(parameters):

    if parameters.program == 'coverage':
        for index, file in enumerate(parameters.bam):
            is_readable(file)
            parameters.bam[index] = os.path.abspath(file)
        is_writeable(parameters.output)
        parameters.output = os.path.abspath(parameters.output)
        assert parameters.threads > 0, 'The "threads" should be a positive integer.'
        assert parameters.bin_size >= 1000, 'The "bin_size" should be greater than or equal to 1000.'
        assert parameters.mapq >= 0, 'The "mapq" should be a non-negative integer.'
        assert parameters.aligned >= 0 and parameters.aligned <= 1, 'The "aligned" should be from 0 to 1.'

    elif parameters.program == 'seed':
        # all files are readable.
        parameters.fasta = os.path.abspath(parameters.fasta)
        is_readable(parameters.fasta)
        # the output file is writeable.
        is_writeable(parameters.output)
        parameters.output = os.path.abspath(parameters.output)
        assert parameters.threads >= 1, 'The "threads" should be a positive integer.'

    else:
        #all files are readable.
        is_readable(parameters.fasta)
        parameters.fasta = os.path.abspath(parameters.fasta)
        is_readable(parameters.seed)
        parameters.seed = os.path.abspath(parameters.seed)
        for index, file in enumerate(parameters.coverage):
            is_readable(file)
            parameters.coverage[index] = os.path.abspath(file)
        #the output file is writeable.
        is_writeable(parameters.output)
        parameters.output = os.path.abspath(parameters.output)
        #all parameters are valid
        assert parameters.clustering_probability >= 0.0 and parameters.clustering_probability <= 1.0, 'The "clustering_probability" should be from 0 to 1.'
        parameters.clustering_probability = numpy.log(parameters.clustering_probability + numpy.finfo(numpy.float64).tiny)
        assert parameters.weight >= 0, 'The "weight" should be a non-negative float.'
        assert parameters.kmer > 0, 'The "kmer" should be a positive integer.'
        assert parameters.min_cluster_size > 0, 'The "min_cluster_size" should be a positive integer.'
        assert parameters.max_dpgmm_distance > 0, 'The "min_cluster_size" should be a positive float.'
        assert parameters.min_dpgmm_size > 0, 'The "min_dpgmm_size" should be a positive integer.'
        assert parameters.min_sequence_length >= 2000, 'The "min_sequence_length" should be a integer greater than or equal to 2000.'
        assert parameters.outlier >= 0.0 and parameters.outlier <= 0.5, 'The "outlier" should be from 0 to 0.5.'
        assert (parameters.sampling_length1[0] > 0) and (parameters.sampling_length1[0] <= parameters.sampling_length1[1]), 'The "sampling_length1" should be two positive integers (low and high) with 0 < low <= high.'
        assert parameters.sampling_number1 > 0, 'The "sampling_number1" should be a positive integer.'
        assert (parameters.sampling_length2[0] > 0) and (parameters.sampling_length2[0] <= parameters.sampling_length2[1]), 'The "sampling_length2" should be two positive integers (low and high) with 0 < low <= high.'
        assert parameters.sampling_number2 > 0, 'The "sampling_number2" should be a positive integer.'
    return None


def main():
    parameters = __init__()

    # Check and ensure that all parameters are valid. #
    parse_parameters(parameters)

    if parameters.program == 'coverage':
        metadecoder_coverage.main(parameters)

    elif parameters.program == 'seed':
        metadecoder_seed.main(parameters)

    else:
        metadecoder_cluster.main(parameters)
