import os
from datetime import datetime
from math import ceil
from multiprocessing import Pool
from re import compile

import numpy
from metadecoder.make_file import make_file
from metadecoder.sam_utility import generate_block, get_sequence_information, read_sam_file


def get_coverage(input_file, block_start, block_end, output_file, bin_size, mapq, aligned_threshold):
    template = compile('(\d+)(\D)')
    sequence2bin_coverages = dict()
    regions = list()
    for _, sequence_id, position, cigar in read_sam_file(input_file, block_start, block_end, mapq):
        bin_index = ceil(position / bin_size)
        # Determine the valid read. #
        aligned_length = 0  # Number of matches and mismatches. #
        total_length = 0
        for cigar_ in template.finditer(cigar):
            operations = int(cigar_.group(1))
            operation = cigar_.group(2)
            if operation in '=MX': # Consume reference. #
                aligned_length += operations
                regions.append((position, position + operations - 1))
                position += operations
            elif operation in 'DN': # Consume reference. #
                position += operations
            total_length += operations
        if aligned_length / total_length >= aligned_threshold: # Valid read. #
            for region_start, region_end in regions:
                while True:
                    if region_end <= bin_index * bin_size: # bin_start <= region_start <= region_end <= bin_end #
                        sequence2bin_coverages[(sequence_id, bin_index)] = sequence2bin_coverages.get((sequence_id, bin_index), 0) + region_end - region_start + 1
                        break
                    elif region_start <= bin_index * bin_size: # bin_start <= region_start <= bin_end < region_end #
                        sequence2bin_coverages[(sequence_id, bin_index)] = sequence2bin_coverages.get((sequence_id, bin_index), 0) + bin_index * bin_size - region_start + 1
                        region_start = bin_index * bin_size + 1
                        bin_index += 1
                    else: # bin_start <= bin_end < region_start <= region_end #
                        bin_index += 1
        regions.clear()
    open_file = open(output_file, 'w')
    for (sequence, bin_index), coverage in sequence2bin_coverages.items():
        open_file.write('\t'.join([sequence, str(bin_index - 1), str(coverage)]) + '\n')
    open_file.close()
    return None


def main(parameters):
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Loading sam files.', flush = True)
    file2temp_files = dict()
    process_pool = Pool(os.cpu_count())
    for file in parameters.sam:
        for block_start, block_end in generate_block(file, parameters.threads):
            temp_file = make_file()
            file2temp_files.setdefault(file, list()).append(temp_file)
            process_pool.apply_async(
                get_coverage,
                (
                    file,
                    block_start,
                    block_end,
                    temp_file,
                    parameters.bin_size,
                    parameters.mapq,
                    parameters.aligned
                )
            )
    process_pool.close()
    process_pool.join()
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Writing to file.', flush = True)
    features = len(parameters.sam)
    sequence2bin_coverages = dict()
    sequence2bin_sizes = dict()
    for sequence, sequence_length in get_sequence_information(parameters.sam[0]):
        sequence2bin_coverages[sequence] = numpy.zeros(shape = (ceil(sequence_length / parameters.bin_size), features), dtype = numpy.int64)
        sequence2bin_sizes[sequence] = numpy.array(
            [
                min(parameters.bin_size * (bin_index + 1), sequence_length) - parameters.bin_size * bin_index for bin_index in range(ceil(sequence_length / parameters.bin_size))
            ],
            dtype = numpy.int64
        )
    for coverage_index, temp_files in enumerate(file2temp_files.values()):
        for temp_file in temp_files:
            open_file = open(temp_file, 'r')
            for line in open_file:
                lines = line.rstrip('\n').split('\t')
                sequence2bin_coverages[lines[0]][int(lines[1]), coverage_index] += float(lines[2])
            open_file.close()
            os.remove(temp_file)
    open_file = open(parameters.output, 'w')
    open_file.write('\t'.join(['sequence id', 'bin index', 'bin size'] + ['coverage' + str(coverage_index + 1) for coverage_index in range(features)]) + '\n')
    for sequence, BinCoverages in sequence2bin_coverages.items():
        for bin_index, (bin_size, bin_coverage) in enumerate(zip(sequence2bin_sizes[sequence], BinCoverages), start = 1):
            open_file.write('\t'.join([sequence, str(bin_index), str(bin_size)] + (bin_coverage / bin_size).astype(numpy.str_).tolist()) + '\n')
    open_file.close()
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Finished.', flush = True)
