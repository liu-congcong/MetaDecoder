import gzip
import os
from math import ceil
from uuid import uuid4


def read_fasta_file(input_fasta):
    '''
    Parameters:
        input_fasta: the path to the (compressed) fasta file.
    Return:
        a generator (sequence_id, sequence)
    '''
    container = list()
    open_file = open(input_fasta, 'rb')
    magic_code = open_file.read(2)
    open_file.close()

    if magic_code == b'\x1f\x8b':
        open_file = gzip.open(input_fasta, mode = 'rt')
    else:
        open_file = open(input_fasta, mode = 'rt')
    for line in open_file:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if container:
                yield (sequence_id, ''.join(container))
            sequence_id = line.split(' ', maxsplit = 1)[0][1 : ]
            container.clear()
        else:
            container.append(line)
    if container:
        yield (sequence_id, ''.join(container))
    open_file.close()
    return None


def split_fasta(input_fasta, output_files):
    '''
    Split a fasta into small files.
    Parameters:
        input_fasta: the path to the fasta file.
        output_files: the number of output files.
    Return:
        a generator of path of each output file.
    '''
    total_size = os.path.getsize(input_fasta)
    block_size = ceil(total_size / output_files) # block_size <= total_size #
    file_position = 0
    file_position_ = 0
    open4r = open(input_fasta, 'rb')
    while file_position < total_size:
        line = open4r.readline()
        file_position += len(line)
        if line.startswith(b'>'):
            file_position -= len(line)
            if file_position > 0:
                open4r.seek(file_position_, os.SEEK_SET)
                output_file = uuid4().hex
                open4w = open(output_file, 'wb')
                while file_position_ < file_position:
                    file_position_ += open4w.write(open4r.read(min(10485760, file_position - file_position_)))
                open4w.close()
                yield output_file
                # file_position_ will be equal to file_position, open4r.tell() will be equal to file_position #
            file_position = open4r.seek(min(file_position + block_size, total_size), os.SEEK_SET)
    open4r.seek(file_position_, os.SEEK_SET)
    output_file = uuid4().hex
    open4w = open(output_file, 'wb')
    while file_position_ < file_position:
        file_position_ += open4w.write(open4r.read(min(10485760, file_position - file_position_)))
    open4w.close()
    yield output_file
    open4r.close()
    return None
