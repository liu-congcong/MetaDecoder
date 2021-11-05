from ctypes import c_longlong
from math import ceil, floor
from multiprocessing import JoinableQueue, Process, sharedctypes

import numpy
from sklearn.svm import SVC


def kmer_to_index(k):
    '''
    Parameters:
        k: the length of k-mers.
    Return:
        (kmer2index, #kmers)
    '''
    complement = str.maketrans('ACGT', 'TGCA')
    acgt = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    index = 0
    kmer2index = dict()
    kmer2index_ = dict()
    for i in range(4 ** k):
        kmer = ''.join(acgt[(i >> 2 * j) & 3] for j in range(k))
        rev_com_kmer = kmer.translate(complement)[ : : -1]
        if kmer > rev_com_kmer:
            kmer_ = rev_com_kmer + kmer
        else:
            kmer_ = kmer + rev_com_kmer
        if kmer_ not in kmer2index_:
            kmer2index_[kmer_] = index
            index += 1
        kmer2index[kmer] = kmer2index_[kmer_]
        kmer2index[rev_com_kmer] = kmer2index_[kmer_]
    return (kmer2index, len(kmer2index_))  # kmer2index, number of kmers #


def generate_kmer_frequency_worker(process_queue, container, k, kmer2index):
    '''
    Parameters:
        process_queue: [sequence, container_offset].
        container: an array to store the result.
        k: the length of k-mers.
        kmer2index: the hash of kmer - index pairs.
    Return:
        None
    '''
    while True:
        sequence, container_offset = process_queue.get()
        if sequence != None:
            for base_index in range(len(sequence) - k + 1):
                kmer = sequence[base_index : base_index + k]
                if kmer in kmer2index:
                    container_index = container_offset + kmer2index[kmer]
                    container[container_index] += 1
            process_queue.task_done()
        else:
            process_queue.task_done()
            break
    return None


def generate_kmer_frequency(sequences, k, kmer2index, kmers, threads):
    '''
    Parameters:
        sequences: the list of sequences.
        k: the length of k-mers.
        kmer2index: the hash of kmer - index pairs.
        kmers: the number of kmers.
        threads: threads.
    Return:
        an array of frequencies of kmers of all sequences.
    '''
    sequences_ = len(sequences)
    container = sharedctypes.RawArray(c_longlong, kmers * sequences_)
    processes = list()
    process_queue = JoinableQueue(threads * 10)
    # Start all threads. #
    for process in range(threads):
        processes.append(
            Process(
                target = generate_kmer_frequency_worker,
                args = (process_queue, container, k, kmer2index)
            )
        )
        processes[-1].start()
    # Put all sequences to queue. #
    for sequence_index, sequence in enumerate(sequences):
        process_queue.put([sequence, kmers * sequence_index])
    for process in range(threads):
        process_queue.put([None, None])
    for process in processes:
        process.join()
    process_queue.join()
    processes.clear()

    kmers_ = numpy.array(container, dtype = numpy.float64).reshape(sequences_, kmers)
    kmers_ += numpy.finfo(numpy.float64).eps
    kmers_ /= numpy.sum(kmers_, axis = 1, keepdims = True)
    return kmers_


def sample_kmer_frequency(sequences, k, kmer2index, kmers, sampling_length, sampling_number, threads, random_number):
    '''
    Parameters:
        sequences: the list of sequences.
        k: the length of k-mers.
        kmer2index: the hash of kmer - index pairs.
        kmers: the number of kmers.
        sampling_length: the length of sampling sequences.
        sampling_number: the number of sampling sequences.
        threads: threads.
        random_number: the random number.
    Return:
        an array of frequencies of kmers of all sampling sequences.
    '''
    random_number_generator = numpy.random.default_rng(random_number)
    sampling_length_ = random_number_generator.integers(sampling_length[0], high = sampling_length[1], size = sampling_number, dtype = numpy.int64, endpoint = True)
    sampling_sequences = list()
    for sequence in sequences:
        sequence_length = len(sequence)
        if sequence_length < sampling_length[1] * sampling_number:
            copies = ceil(sampling_length[1] * sampling_number / sequence_length)
            sequence *= copies
            sequence_length *= copies
        step = floor(sequence_length / sampling_number)
        for index in range(sampling_number):
            sampling_sequences.append(sequence[index * step : index * step + sampling_length_[index]])
    return generate_kmer_frequency(sampling_sequences, k, kmer2index, kmers, threads)


def run_svm(x4training, y4training, x4predicting, probability, random_number):
    '''
    Parameters:
        x4training: the training data.
        y4training: the training data.
        x4predicting: the unknown data.
        probability: True | False.
        random_number: the random number.
    Return:
        an array of predicted labels (probability = False) or probabilities (probability = True) of unknown data.
    '''
    # initialize multi-class svm #
    svm = SVC(
        C = 1,
        kernel = 'rbf',
        gamma = 1 / (x4training.shape[1]),
        probability = probability,
        cache_size = 4096,
        decision_function_shape = 'ovo',
        random_state = random_number
    )
    # standardization #
    x4training_mean = numpy.mean(x4training, axis = 0)
    x4training_std = numpy.std(x4training, axis = 0)
    x4training_std[x4training_std < numpy.finfo(numpy.float64).eps] = numpy.finfo(numpy.float64).eps
    # fit model #
    svm.fit((x4training - x4training_mean) / x4training_std, y4training, sample_weight = None)
    # make predictions #
    if probability:
        predictions = svm.predict_log_proba((x4predicting - x4training_mean) / x4training_std)
    else:
        predictions = svm.predict((x4predicting - x4training_mean) / x4training_std)
    return predictions


def run_svm_worker(container, offset, x4training_list, y4training_list, x4predicting, random_number):
    '''
    Parameters:
        container: an array to store the result.
        offset: the start index of the container to save the result produced by a worker.
        x4training_list: the list of training data.
        y4training_list: the list of training data.
        x4predicting: the unknown data.
        random_number: the random number.
    Return:
        None
    '''
    for x4training, y4training in zip(x4training_list, y4training_list):
        container[offset : offset + x4predicting.shape[0]] = run_svm(x4training, y4training, x4predicting, False, random_number)
        offset += x4predicting.shape[0]
    return None
