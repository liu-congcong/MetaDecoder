import os
from math import ceil


def read_sequence_header(input_sam):
    '''
    Parameters:
        input_sam: the input sam file (with @SQ headers).
    Return:
        a generator of @SQ line.
    '''
    open4r = open(input_sam, 'r')
    for line in open4r:
        line = line.rstrip('\n')
        if line.startswith('@'):
            if line.startswith('@SQ'):
                yield line
        else: # consume all headers #
            break
    open4r.close()
    return None


def get_sequence_information(input_sam):
    '''
    Extract the information of all sequences from a sam file.
    Parameters:
        input_sam: the input sam file (with @SQ headers).
    Return:
        a generator of (sequence_id, sequence_length)
    '''
    for line in read_sequence_header(input_sam):
        lines = line.rstrip('\n').split('\t')
        for tag_value in lines:
            if tag_value.startswith('SN'):
                sequence_id = tag_value[3:]
            elif tag_value.startswith('LN'):
                sequence_length = int(tag_value[3:])
        yield (sequence_id, sequence_length)
    return None


def generate_block(input_sam, blocks):
    '''
    Split a sam file into some parts.
    Parameters:
        input_sam: the input sam file (with @SQ headers).
        blocks: the number of blocks.
    Return:
        a generator of start and end positions of a block of a sam file.
    '''
    open4r = open(input_sam, 'rb')
    while True:
        line = open4r.readline()
        if line:
            if not line.startswith(b'@'):
                block_start = open4r.tell()-len(line)
                break
        else: # EOF #
            break
    file_size = os.path.getsize(input_sam)
    block_size = ceil((file_size - block_start) / blocks)
    while block_start < file_size:
        open4r.seek(block_size, 1)
        open4r.readline()
        block_end = open4r.tell()
        yield(block_start, min(block_end, file_size))
        block_start = block_end
    open4r.close()
    return None


def read_sam_file(input_sam, block_start, block_end, mapq = 0):
    '''
    Parameters:
        input_sam: the input sam file (with @SQ headers).
        block_start: the start position of a block.
        block_end: the end position of a block.
        mapq: MAPQ.
    Return:
        (read_id, ref_id, pos, cigar)
    '''
    open4r = open(input_sam, 'rb')
    open4r.seek(block_start, 0)
    while block_start < block_end:
        line = open4r.readline()
        if not line:
            break
        lines = line.rstrip(b'\n').split(b'\t')
        if (~ int(lines[1]) & 4) and (int(lines[4]) >= mapq):
            #read id, reference sequence, position, cigar
            yield (lines[0].decode('ascii'), lines[2].decode('ascii'), int(lines[3]), lines[5].decode('ascii'))
        block_start += len(line)
    open4r.close()
    return None