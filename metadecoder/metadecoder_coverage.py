from ctypes import c_int64
from datetime import datetime
from math import ceil
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value

import numpy
from threadpoolctl import threadpool_limits

from .bam import getUngappedRegions, indexBam, readIndices
from .plot import plotBar


def workerProcess(queue, files, mapq, identity, binSize, output, n, N):
    '''
    queue: queue
    files: list
    mapq: int
    identity: float
    binSize: int
    output: string
    n: int
    N: Value(c_int64, 0)
    '''
    m = len(files)
    while True:
        sequence, length, fileOffsets, dataOffsets, dataSizes = queue.get()
        if sequence is None:
            break
        x = numpy.zeros(shape = (ceil(length / binSize), m), dtype = numpy.float64)
        binSizes = numpy.full(shape = x.shape[0], fill_value = binSize, dtype = numpy.int64)
        binSizes[-1] = length % binSize
        if binSizes[-1] == 0:
            binSizes[-1] = binSize
        for i, (file, fileOffset, dataOffset, dataSize) in enumerate(zip(files, fileOffsets, dataOffsets, dataSizes)):
            for position, regions in getUngappedRegions(file, fileOffset, dataOffset, dataSize, mapq, identity):
                binIndex = ceil(position / binSize)
                for regionStart, regionEnd in regions:
                    while True:
                        if regionEnd <= binIndex * binSize: # binStart <= regionStart <= regionEnd <= binEnd #
                            x[binIndex - 1, i] += regionEnd - regionStart + 1
                            break
                        elif regionStart <= binIndex * binSize: # binStart <= regionStart <= binEnd < regionEnd #
                            x[binIndex - 1, i] += binIndex * binSize - regionStart + 1
                            regionStart = binIndex * binSize + 1
                            binIndex += 1
                        else: # binStart <= binEnd < regionStart <= regionEnd #
                            binIndex += 1
        x /= binSizes[ : , None]
        N.acquire()
        openFile = open(output, 'a')
        for i, (binSizeI, xi) in enumerate(zip(binSizes, x), start = 1):
            xi = '\t'.join(xi.astype(numpy.str_).tolist())
            openFile.write(f'{sequence}\t{i}\t{binSizeI}\t{xi}\n')
        openFile.close()
        N.value += 1
        plotBar(N.value / n)
        N.release()
    return None


def createProcesses(files, n, mapq, identity, binSize, threads, output):
    N = Value(c_int64, 0)
    queue = Queue()
    processes = list()
    with threadpool_limits(limits = 1):
        for i in range(threads):
            processes.append(Process(target = workerProcess, args = (queue, files, mapq, identity, binSize, output, n, N)))
            processes[-1].start()
    return (queue, processes, N)


def freeProcesses(queue, processes):
    for process in processes:
        queue.put((None, None, None, None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def main(parameters):

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Indexing all bam files.', flush = True)
    indexBam(parameters.bam, parameters.threads)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Loading all index files.', flush = True)
    sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndices(parameters.bam, parameters.threads)

    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Computing coverage.', flush = True)
    openFile = open(parameters.output, 'w')
    openFile.write('\t'.join(['sequence id', 'bin index', 'bin size'] + ['coverage' + str(coverage_index + 1) for coverage_index in range(len(parameters.bam))]) + '\n')
    openFile.close()
    queue, processes, N = createProcesses(parameters.bam, len(sequences), parameters.mapq, parameters.aligned, parameters.bin_size, parameters.threads, parameters.output)
    for sequence, length, fileOffset, dataOffset, dataSize in zip(sequences, lengths, fileOffsets, dataOffsets, dataSizes):
        queue.put((sequence, length, fileOffset, dataOffset, dataSize))
    freeProcesses(queue, processes)
    print(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} -> Finished.', flush = True)
