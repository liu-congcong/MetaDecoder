import os
from ctypes import c_double, c_longlong
from datetime import datetime
from math import ceil
from multiprocessing import JoinableQueue, Lock, Process, sharedctypes

import numpy
from metadecoder.coverage_model import GMM
from metadecoder.fasta_utility import read_fasta_file
from metadecoder.isolation_forest import isolation_forest
from metadecoder.kmer_frequency_model import generate_kmer_frequency, kmer_to_index, run_svm, sample_kmer_frequency
from metadecoder.seed_selection import generate_seed, select_seed
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from threadpoolctl import threadpool_limits


def read_coverage_file(file, sequence_id2sequence, sequences):
    sequence2bin_coverage = dict()
    open_file = open(file, 'r')
    # header: sequence_id BinIndex BinSize Coverage1 Coverage2 ... #
    coverages = len(open_file.readline().rstrip('\n').split('\t')) - 3
    for line in open_file:
        lines = line.rstrip('\n').split('\t')
        if lines[0] in sequence_id2sequence:
            sequence2bin_coverage.setdefault(sequence_id2sequence[lines[0]], list()).append([float(coverage_) + 1 / float(lines[2]) for coverage_ in lines[3:]])
    open_file.close()

    bin_coverage = list()
    coverage = numpy.empty(shape = (sequences, coverages), dtype = numpy.float64)
    for sequence in range(sequences):
        bin_coverage_ = sequence2bin_coverage.get(sequence, [[1e-5 for coverage_ in range(coverages)]])
        if len(bin_coverage_) > 1:
            del bin_coverage_[-1]
        bin_coverage.append(numpy.log(bin_coverage_, dtype = numpy.float64))
        coverage[sequence] = numpy.mean(bin_coverage[sequence], axis = 0)
    return (coverage, bin_coverage)


def read_seed_file(file, sequence_id2sequence):
    '''
    Load seed file.
    '''
    sequence2markers = dict()
    marker2sequences = dict()
    clusters2weight = dict()
    open_file = open(file, 'r')
    for line in open_file:
        lines = line.rstrip('\n').split('\t')
        sequences = [sequence_id2sequence[sequence_id] for sequence_id in lines[1 : ] if sequence_id in sequence_id2sequence]
        if sequences:
            marker2sequences[lines[0]] = sequences
            clusters2weight[len(sequences)] = clusters2weight.get(len(sequences), 0) + 1
    open_file.close()
    clusters = max((weight, clusters) for clusters, weight in clusters2weight.items())[1]
    if clusters == 1: # Extreme case. #
        clusters = max(clusters2weight.keys())
        if clusters > 1:
            for marker, sequences in list(marker2sequences.items()):
                if len(sequences) == 1:
                    del marker2sequences[marker]
    for marker, sequences in marker2sequences.items():
        if len(sequences) <= clusters * 2:
            for sequence in set(sequences):
                sequence2markers.setdefault(sequence, list()).append(marker)
    return sequence2markers


'''
# Unpublished function #
def read_cluster_file(file, sequence_ids):
    # Load cluster file. #
    sequence_ids_list = list()
    open_file = open(file, 'r')
    for line in open_file:
        sequence_ids_ = set(line.rstrip('\n').split('\t')) & sequence_ids
        if sequence_ids_:
            sequence_ids_list.append(sorted(sequence_ids_))
    open_file.close()
    return sequence_ids_list
'''


def cluster_sequences(sequences, probabilities, min_probability):
    cluster_indices = numpy.argmax(probabilities, axis = 1)
    probabilities_ = probabilities[numpy.arange(sequences.shape[0]), cluster_indices]
    mask = probabilities_ >= min_probability
    return ([sequences[(cluster_indices == cluster_index) & mask] for cluster_index in numpy.unique(cluster_indices[mask])], sequences[~mask])


def output_clusters(sequence_ids, SEQUENCES, cluster_indices, output_unclustered_sequences, file_prefix, length_per_line = 100):
    for file_index, cluster_index in enumerate(numpy.unique(cluster_indices)):
        if file_index:
            open_file = open('{0}.{1}.fasta'.format(file_prefix, file_index), 'w')
        elif output_unclustered_sequences:
            open_file = open('{0}.{1}.fasta'.format(file_prefix, 'unclustered'), 'w')
        else:
            continue
        for sequence in numpy.flatnonzero(cluster_indices == cluster_index):
            open_file.write('>' + sequence_ids[sequence] + '\n')
            SEQUENCE = SEQUENCES[sequence]
            index = 0
            while open_file.write(SEQUENCE[index : index + length_per_line]):
                open_file.write('\n')
                index += length_per_line
        open_file.close()
    return None


#debug#
def dump_kmer_frequency(file, kmer_frequency):
    # Dump the kmer frequencies to file. #
    open_file = open(file, 'w')
    # Write header to file. #
    open_file.write('\t'.join(str(index) for index in range(kmer_frequency.shape[1])) + '\n')
    for kmer_frequency_ in kmer_frequency:
        open_file.write('\t'.join(kmer_frequency_.astype(numpy.str_).tolist()) + '\n')
    open_file.close()
    return None


def load_kmer_frequency(file):
    kmer_frequency = list()
    open_file = open(file, 'r')
    open_file.readline()
    for line in open_file:
        kmer_frequency.append(line.rstrip('\n').split('\t'))
    open_file.close()
    return numpy.array(kmer_frequency, dtype = numpy.float64)


def dump_dpgmm_prediction(file, predictions):
    open_file = open(file, 'w')
    for prediction in predictions:
        open_file.write(str(prediction) + '\n')
    open_file.close()
    return None


def load_dpgmm_prediction(file):
    predictions = list()
    open_file = open(file, 'r')
    for line in open_file:
        predictions.append(line.rstrip('\n'))
    open_file.close()
    return numpy.array(predictions, dtype = numpy.int64)


def read_mapping_file(input_file):
    open_file = open(input_file, 'r')
    open_file.readline() # Remove header. #
    for line in open_file:
        lines = line.rstrip('\n').split('\t')
        yield (lines[0], lines[1])
    open_file.close()
    return None


def calculate_average_distance(x):
    xx = numpy.sum(numpy.square(x), axis = 1, keepdims = True) # (samples, 1) #
    distance_matrix = x @ x.T
    distance_matrix *= -2.0
    distance_matrix += xx
    distance_matrix += xx.T
    distance_matrix[distance_matrix < 0.0] = 0.0
    return numpy.sum(numpy.sqrt(distance_matrix)) / (x.shape[0] * (x.shape[0] - 1) + numpy.finfo(numpy.float64).eps)


def run_models(process_queue, container, offset, kmer, kmer2index, kmers, sampling_length1, sampling_number1, sampling_length2, sampling_number2, weight, min_clustering_probability, outlier, random_number):
    while True:
        sequences, clusters, seeds_list, length, kmer_frequency, coverage, bin_coverage, SEQUENCES = process_queue.get()
        if clusters > 0:
            # Run the seed selection model. #
            seeds = select_seed(
                sequences,
                seeds_list, # groups of seeds #
                clusters, # number of clusters estimated by seeds #
                SEQUENCES, # sequence #
                coverage, # coverage #
                kmer_frequency, # kmer frequence matrix #
                kmer, # k #
                kmer2index, # hash #
                kmers, # number of kinds of kmer #
                sampling_length1, # length of sequence #
                sampling_number1, # number of sampling #
                random_number # random number #
            )
            if seeds:
                # Run the kmer frequency probabilistic model. #
                x4training = sample_kmer_frequency(seeds, kmer, kmer2index, kmers, sampling_length2, sampling_number2, 1, random_number)
                y4training = numpy.repeat(numpy.arange(clusters, dtype = numpy.int64), sampling_number2)
                clustering_probability = run_svm(x4training, y4training, kmer_frequency, True, random_number)
                
                # Run the coverage probabilistic model. #
                gmm = GMM(clusters, bin_coverage, clustering_probability, weight)
                gmm.main()
                clustering_probability = gmm.log_responsibility
            else:
                gmm = GaussianMixture(n_components = clusters, covariance_type = 'full', n_init = 1, init_params = 'kmeans', random_state = random_number)
                x = numpy.concatenate((kmer_frequency, coverage), axis = 1)
                gmm.fit(x)
                clustering_probability = numpy.log(gmm.predict_proba(x) + numpy.finfo(numpy.float64).eps)

            clustered_sequences_list, unclustered_sequences = cluster_sequences(sequences, clustering_probability, min_clustering_probability)
            if len(clustered_sequences_list) > 1:
                for clustered_sequences in clustered_sequences_list:
                    container_value = numpy.min(clustered_sequences)
                    for sequence in clustered_sequences:
                        container[sequence] = container_value
                if unclustered_sequences.shape[0] > 0:
                    container_value = -offset - 1
                    for sequence in unclustered_sequences:
                        container[sequence] = container_value
            else:
                sequences_, unclustered_sequences = isolation_forest(sequences, numpy.concatenate((kmer_frequency, coverage), axis = 1), length, 1, outlier, random_number)
                if unclustered_sequences.shape[0] > 0:
                    container_value = numpy.min(sequences_)
                    for sequence in sequences_:
                        container[sequence] = container_value
                    container_value = -offset - 1
                    for sequence in unclustered_sequences:
                        container[sequence] = container_value
                else:
                    container_value = -offset - 1
                    for sequence in sequences_:
                        container[sequence] = container_value

            process_queue.task_done()
        else: # All tasks have been finished. #
            process_queue.task_done()
            break
    return None


def main(parameters):
    # Import DPGMM class. #
    if parameters.disable_gpu:
        from metadecoder.dirichlet_process_gaussian_mixture import DPGMM
    else:
        try:
            from metadecoder.dirichlet_process_gaussian_mixture_gpu import DPGMM
        except Exception:
            from metadecoder.dirichlet_process_gaussian_mixture import DPGMM

    # read fasta file, return list, list #
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Loading fasta file.', flush = True)
    sequence_ids = list()
    SEQUENCES = list()
    length = list()
    sequence_id2sequence = dict()
    sequence_index = 0
    for sequence_id, SEQUENCE in read_fasta_file(parameters.fasta):
        length_ = len(SEQUENCE)
        if length_ >= parameters.min_sequence_length:
            sequence_ids.append(sequence_id)
            SEQUENCES.append(SEQUENCE)
            length.append(length_)
            sequence_id2sequence[sequence_id] = sequence_index
            sequence_index += 1
    length = numpy.array(length)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)

    sequences = numpy.arange(len(sequence_ids), dtype = numpy.int64)

    # Read coverage file, return array. #
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Loading coverage file.', flush = True)
    coverage, bin_coverage = read_coverage_file(parameters.coverage, sequence_id2sequence, sequences.shape[0])
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)

    # Kmer hash table. #
    kmer2index, kmers = kmer_to_index(parameters.kmer)

    # Calculate kmer frequency. #
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Counting kmers of all sequences.', flush = True)
    if os.access(os.path.basename(parameters.fasta) + '.' + str(parameters.min_sequence_length) + '.metadecoder.kmers', os.R_OK):
        kmer_frequency = load_kmer_frequency(os.path.basename(parameters.fasta) + '.' + str(parameters.min_sequence_length) + '.metadecoder.kmers')
    else:
        kmer_frequency = generate_kmer_frequency(SEQUENCES, parameters.kmer, kmer2index, kmers, os.cpu_count())
        dump_kmer_frequency(os.path.basename(parameters.fasta) + '.' + str(parameters.min_sequence_length) + '.metadecoder.kmers', kmer_frequency)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)

    # Read seed file, return hash. #
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Loading seeds file.', flush = True)
    sequence2markers = read_seed_file(parameters.seed, sequence_id2sequence)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)

    total_sequences = sequences.shape[0]
    # Initialize the container. #
    container = sharedctypes.RawArray(c_longlong, total_sequences)
    processes = list()
    lock = Lock()

    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Running the DPGMM algorithm to obtain clusters.', flush = True)
    pca = PCA(n_components = 0.9, whiten = False, random_state = parameters.random_number)
    if not os.access(os.path.basename(parameters.fasta) + '.' + str(parameters.min_sequence_length) + '.metadecoder.dpgmm', os.R_OK):
        dpgmm_unit = numpy.mean(length)
        clusters = max(generate_seed(sequences, sequence2markers)[0], 10)
        if clusters <= 500:
            clusters *= 3.0
        elif clusters <= 1000:
            clusters *= 2.0
        dpgmm = DPGMM(
            int(clusters), # number of clusters #
            parameters.min_dpgmm_size / dpgmm_unit, # size of the smallest cluster #
            numpy.concatenate((pca.fit_transform(kmer_frequency), coverage), axis = 1), # input data #
            length / dpgmm_unit, # data weights #
            parameters.random_number # random seed used by multinomial fai initialization #
        )
        dpgmm.main()
        dpgmm_predictions = numpy.argmax(dpgmm.multinomial_fai, axis = 1)
        dump_dpgmm_prediction(os.path.basename(parameters.fasta) + '.' + str(parameters.min_sequence_length) + '.metadecoder.dpgmm', dpgmm_predictions)
    else:
        dpgmm_predictions = load_dpgmm_prediction(os.path.basename(parameters.fasta) + '.' + str(parameters.min_sequence_length) + '.metadecoder.dpgmm')
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)

    for dpgmm_prediction in numpy.unique(dpgmm_predictions):
        sequences_ = numpy.flatnonzero(dpgmm_predictions == dpgmm_prediction)
        if calculate_average_distance(kmer_frequency[sequences_]) <= parameters.max_dpgmm_distance:
            container_value = numpy.min(sequences_)
        else:
            container_value = -total_sequences - 1
        for sequence in sequences_:
            container[sequence] = container_value

    # Start all processes. #
    process_queue = JoinableQueue(int(os.cpu_count() * 10))
    with threadpool_limits(limits = 1):
        for process_index in range(os.cpu_count()):
            processes.append(
                Process(
                    target = run_models,
                    args = [
                        process_queue,
                        container,
                        total_sequences,
                        parameters.kmer,
                        kmer2index,
                        kmers,
                        parameters.sampling_length1,
                        parameters.sampling_number1,
                        parameters.sampling_length2,
                        parameters.sampling_number2,
                        parameters.weight / coverage.shape[1],
                        parameters.clustering_probability,
                        parameters.outlier,
                        parameters.random_number
                    ]
                )
            )
            processes[-1].start()


    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Running the kmer frequency and coverage models for clustering.', flush = True)
    while True:
        container_array = numpy.asarray(container, dtype = numpy.int64)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', '{0:.2%} of the total sequences have been processed.'.format(numpy.sum(container_array < 0) / total_sequences), flush = True)
        if numpy.any(container_array >= 0):
            for container_value in numpy.unique(container_array[container_array >= 0]):
                sequences = numpy.flatnonzero(container_array == container_value)
                clusters, seeds_list = generate_seed(sequences, sequence2markers)
                if clusters > 1:
                    process_queue.put(
                        [
                            sequences,
                            clusters,
                            seeds_list,
                            length[sequences],
                            kmer_frequency[sequences],
                            coverage[sequences],
                            [bin_coverage[sequence] for sequence in sequences], # list of coverage matrix #
                            [SEQUENCES[sequence] for sequence in sequences]
                        ]
                    )
                else: # clusters <= 1 #
                    if clusters == 1 and numpy.sum(length[sequences]) >= parameters.min_cluster_size:
                        container_value = numpy.min(sequences) - total_sequences
                    else:
                        container_value = -total_sequences - 1
                    for sequence in sequences:
                        container[sequence] = container_value
        else:
            break
        process_queue.join()
    for process_index in range(os.cpu_count()):
        process_queue.put([None, 0, None, None, None, None, None, None])
    for process in processes:
        process.join()
    process_queue.join()
    processes.clear()

    output_clusters(sequence_ids, SEQUENCES, numpy.asarray(container, dtype = numpy.int64), parameters.output_unclustered_sequences, parameters.output)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Finished.', flush = True)