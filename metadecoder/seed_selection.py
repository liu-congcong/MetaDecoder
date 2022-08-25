from math import ceil
from operator import itemgetter

import numpy
from sklearn.cluster import SpectralClustering

from .kmer_frequency_model import run_svm, sample_kmer_frequency


def generate_seed(sequences, sequence2markers):
    '''
    Obtain sequences mapped with marker genes.
    Parameters:
        sequences: the list of sequences.
        sequence2markers: the hash of sequence - markers pairs.
    Return:
        (clusters, sequences_list)
    '''
    sequences_list = list()
    marker2sequences = dict()
    clusters2sequences_list = dict()
    clusters2weight = dict()
    for sequence in sequences:
        if sequence in sequence2markers:
            for marker in sequence2markers[sequence]:
                marker2sequences.setdefault(marker, list()).append(sequence)
    if marker2sequences:
        for marker, sequences_ in marker2sequences.items():
            sequences_ = tuple(sorted(sequences_)) # sequences_ must be a sorted tuple, not an array. #
            clusters = len(sequences_)
            clusters2sequences_list.setdefault(clusters, list()).append(sequences_)
            clusters2weight[clusters] = clusters2weight.get(clusters, 0) + 1
        lower_clusters, weight = max(clusters2weight.items(), key = itemgetter(1))
        clusters = max(clusters_ for clusters_, weight_ in clusters2weight.items() if weight_ >= weight * 0.5)
        for clusters_, sequences_list_ in clusters2sequences_list.items():
            if (clusters_ >= max(lower_clusters, 2)) and (clusters_ <= ceil(clusters * 2)):
                sequences_list.extend(sequences_list_)
        if clusters == 1 and len(sequences_list) > 5:
            clusters = 2
    else:
        clusters = 0
    return (clusters, sequences_list)


def select_seed(sequences, seeds_list, clusters, SEQUENCES, coverage, kmer_frequency, kmer, kmer2index, kmers, sampling_length, sampling_number, random_number):
    seeds = list()
    unique_seeds = numpy.unique([seed for seeds in seeds_list for seed in seeds])
    seed_mappings = numpy.flatnonzero(numpy.isin(sequences, unique_seeds, assume_unique = True))

    if unique_seeds.shape[0] > clusters:
        seeds_list = [numpy.array(seeds) for seeds in sorted(set(seeds_list))]
        seed_indices_list = [numpy.flatnonzero(numpy.isin(unique_seeds, seeds, assume_unique = True)) for seeds in seeds_list]

        # Weights of all seeds not in seeds_list are set to zeros. #
        seeds_list_length = len(seeds_list)
        scores = numpy.zeros(shape = (unique_seeds.shape[0], unique_seeds.shape[0]), dtype = numpy.float64)

        # create training data #
        x4training = numpy.concatenate(
            (
                sample_kmer_frequency(
                    [SEQUENCES[seed_mapping] for seed_mapping in seed_mappings],
                    kmer,
                    kmer2index,
                    kmers,
                    sampling_length,
                    sampling_number,
                    1,
                    random_number
                ),
                numpy.repeat(coverage[seed_mappings], sampling_number, axis = 0)
            ),
            axis = 1
        )
        y4training = numpy.repeat(numpy.arange(unique_seeds.shape[0]), sampling_number)

        # create evaluating data #
        x4evaluating = numpy.concatenate((kmer_frequency[seed_mappings], coverage[seed_mappings]), axis = 1)

        for seed_indices in seed_indices_list:
            training_indices = numpy.isin(y4training, seed_indices, assume_unique = False)
            svm_predictions = run_svm(x4training[training_indices], y4training[training_indices], x4evaluating, False, random_number)
            valid_count = 0
            for seed_indices_ in seed_indices_list:
                if numpy.unique(svm_predictions[seed_indices_]).shape[0] == min(seed_indices.shape[0], seed_indices_.shape[0]):
                    valid_count += 1
            if valid_count / seeds_list_length >= 0.5: # This group of seeds is valid. #
                for seed_indices_ in seed_indices_list:
                    svm_predictions_ = svm_predictions[seed_indices_]
                    for svm_prediction in numpy.unique(svm_predictions_):
                        scores[svm_prediction, seed_indices_[svm_predictions_ == svm_prediction]] += 1.0 / numpy.sum(svm_predictions_ == svm_prediction)
        scores += scores.T
        scores_mask = numpy.flatnonzero(numpy.sum(scores, axis = 0))
        if scores_mask.shape[0] > clusters:
            scores = scores[scores_mask][ : , scores_mask]
            scores += 1e-6
            # Spectral clustering may fail due to an extremely rare vulnerability. #
            try:
                spectral_clustering = SpectralClustering(n_clusters = clusters, random_state = random_number, affinity = 'precomputed')
                spectral_predictions = spectral_clustering.fit_predict(scores)
                for spectral_prediction in range(clusters):
                    seeds.append(''.join(SEQUENCES[seed_mappings[sequence_index]] for sequence_index in scores_mask[spectral_predictions == spectral_prediction]))
            except:
                pass
    return seeds