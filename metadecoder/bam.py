import ctypes
import os
from gzip import GzipFile
from hashlib import md5
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import RawArray, Value
from struct import Struct

import numpy

from .c import findCPath
from .plot import plotBar


class BGZF(ctypes.Structure):
    _fields_ = [
        ('filePointer', ctypes.c_void_p),
        ('_1', ctypes.POINTER(ctypes.c_uint8)),
        ('_2', ctypes.POINTER(ctypes.c_uint8)),
        ('x', ctypes.POINTER(ctypes.c_uint8)),
        ('_3', ctypes.c_uint64),
        ('_4', ctypes.c_uint64),
        ('fileOffset', ctypes.c_uint64),
        ('_5', ctypes.c_uint64),
        ('dataSize', ctypes.c_uint64),
        ('dataOffset', ctypes.c_uint64),
    ]

BGZFpointer = ctypes.POINTER(BGZF)
bgzf = ctypes.cdll.LoadLibrary(findCPath('bgzf'))
bgzf.bgzfOpen.argtypes = [ctypes.c_char_p, ctypes.c_uint64]
bgzf.bgzfOpen.restype = BGZFpointer
bgzf.bgzfRead.argtypes = [BGZFpointer, ctypes.c_uint64]
bgzf.bgzfRead.restype = ctypes.c_int
bgzf.bgzfClose.argtypes = [BGZFpointer]
bgzf.bgzfClose.restype = ctypes.c_int

struct1B = Struct('<B') # uint8 #
struct1c1i = Struct('<ci')
struct1H = Struct('<H') # uint16 #
struct1i = Struct('<i')
struct1I = Struct('<I') # uint32 #
struct2i = Struct('<2i')
struct1Q = Struct('<Q') # uint64 #
struct3Q = Struct('<QQQ')
struct1i3Q = Struct('<iQQQ')
struct2i2B2x2H1i = Struct('<2i2B2x2Hi') # 20 bytes #
struct1i2B2x2H1i = Struct('<i2B2x2Hi') # 16 bytes #
struct3c = Struct('<3c')
structHash = {
    b'A': (0, Struct('<c'), 1), b'B': (2, None, None)     ,
    b'c': (0, Struct('<b'), 1), b'C': (0, Struct('<B'), 1),
    b'f': (0, Struct('<f'), 4), b'H': (1, None, None)     ,
    b'i': (0, struct1i, 4)    , b'I': (0, struct1I, 4)    ,
    b's': (0, Struct('<h'), 2), b'S': (0, Struct('<H'), 2),
    b'Z': (1, None, None)     ,
}


def isBGZFBam(file):
    openFile = open(file, 'rb')
    magic = openFile.read(2)
    openFile.seek(-28, 2)
    eof = openFile.read(28)
    openFile.close()
    assert magic == b'\x1f\x8b', f'\"{file}\" is not a BGZF compressed bam file.'
    assert eof == b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00', f'\"{file}\" is an incomplete bam file.'
    return None


def readAlignment(file, fileOffset, dataOffset, dataSize):
    openFile = open(file, 'rb')
    openFile.seek(fileOffset, 0)
    gzipFile = GzipFile(fileobj = openFile, mode = 'rb')
    gzipFile.read(dataOffset)
    while dataSize:
        n = struct1i.unpack(gzipFile.read(4))[0]
        alignment = gzipFile.read(n)
        yield (alignment, n)
        dataSize -= 4 + n
    gzipFile.close()
    openFile.close()
    return None


def readHeader(file):
    sequences = list()
    lengths = list()
    bgzfPointer = bgzf.bgzfOpen(file.encode('utf-8'), 0)
    bgzf.bgzfRead(bgzfPointer, 8)
    x = ctypes.string_at(bgzfPointer.contents.x, 8)
    assert x[ : 4] == b'BAM\1', f'\"{file}\" is not a valid bam file.'
    headerLength = struct1i.unpack_from(x, 4)[0]
    bgzf.bgzfRead(bgzfPointer, headerLength + 4)
    x = ctypes.string_at(bgzfPointer.contents.x, headerLength + 4)
    hdLine = x[ : -4].split(b'\n', maxsplit = 1)[0]
    assert b'SO:coordinate' in hdLine, f'\"{file}\" is not a sorted bam file.' # @ lines #
    n = struct1i.unpack_from(x, headerLength)[0]
    assert n > 0, f'\"{file}\" contains no sequences.'
    marker = md5()
    for i in range(n): # all sequences and lengths #
        bgzf.bgzfRead(bgzfPointer, 4)
        x = ctypes.string_at(bgzfPointer.contents.x, 4)
        sequenceLength = struct1i.unpack(x)[0]
        bgzf.bgzfRead(bgzfPointer, sequenceLength + 4)
        x = ctypes.string_at(bgzfPointer.contents.x, sequenceLength + 4)
        marker.update(x[ : sequenceLength - 1])
        sequences.append(x[ : sequenceLength - 1])
        marker.update(x[sequenceLength : ])
        lengths.append(x[sequenceLength : ]) # <i #
    fileOffset = bgzfPointer.contents.fileOffset
    dataOffset = bgzfPointer.contents.dataOffset
    bgzf.bgzfClose(bgzfPointer)
    return (fileOffset, dataOffset, marker, sequences, lengths)


def createIndex(marker, sequences, lengths, fileOffsets, dataOffsets, dataSizes, file):
    '''
    sequences: [bytes , ... , bytes]
    lengths: [bytes , ... , bytes]
    '''
    sequences.append(b'') # umapped sequence #
    lengths.append(b'\x00\x00\x00\x00') # umapped sequence length #
    gzipFile = GzipFile(filename = file, mode = 'wb', compresslevel = 9)
    gzipFile.write(marker.digest())
    gzipFile.write(struct1Q.pack(len(sequences))) # n = #sequences + 1 #
    for (sequence, length, fileOffset, dataOffset, dataSize) in zip(sequences, lengths, fileOffsets, dataOffsets, dataSizes):
        gzipFile.write(struct1i.pack(len(sequence)))
        gzipFile.write(sequence)
        gzipFile.write(length)
        gzipFile.write(struct3Q.pack(fileOffset, dataOffset, dataSize))
    gzipFile.close()
    return None


def indexProcess(queue, n, N):
    while True:
        file = queue.get()
        if file is None:
            break
        if not os.access(f'{file}.index', os.R_OK):
            headerFileOffset, headerDataOffset, marker, sequences, lengths = readHeader(file)
            nSequences = len(sequences)
            fileOffsets = numpy.zeros(shape = nSequences + 1, dtype = numpy.int64)
            dataOffsets = numpy.zeros(shape = nSequences + 1, dtype = numpy.int64)
            dataSizes = numpy.zeros(shape = nSequences + 1, dtype = numpy.int64)

            dataSize = 0
            reference_ = -2
            fileOffset = headerFileOffset
            fileOffset_ = headerFileOffset
            dataOffset = headerDataOffset
            dataOffset_ = headerDataOffset
            bgzfPointer = bgzf.bgzfOpen(file.encode('utf-8'), headerFileOffset)
            bgzf.bgzfRead(bgzfPointer, headerDataOffset) # skip header #
            while not bgzf.bgzfRead(bgzfPointer, 8):
                x = ctypes.string_at(bgzfPointer.contents.x, 8)
                alignmentSize, reference = struct2i.unpack(x) # reference: 0 ... n, -1 #
                bgzf.bgzfRead(bgzfPointer, alignmentSize - 4)
                if reference != reference_:
                    if reference_ >= -1:
                        fileOffsets[reference_] = fileOffset_
                        dataOffsets[reference_] = dataOffset_
                        dataSizes[reference_] = dataSize
                    reference_ = reference
                    fileOffset_ = fileOffset
                    dataOffset_ = dataOffset
                    dataSize = 0
                fileOffset = bgzfPointer.contents.fileOffset
                dataOffset = bgzfPointer.contents.dataOffset
                dataSize += 4 + alignmentSize
            if reference_ >= -1:
                fileOffsets[reference_] = fileOffset_
                dataOffsets[reference_] = dataOffset_
                dataSizes[reference_] = dataSize
            bgzf.bgzfClose(bgzfPointer)
            createIndex(marker, sequences, lengths, fileOffsets, dataOffsets, dataSizes, f'{file}.index')
        n.acquire()
        n.value += 1
        plotBar(n.value / N)
        n.release()
    return None


def indexBam(files, threads):
    '''
    files: sorted bam files.
    threads: int
    '''
    N = len(files)
    n = Value(ctypes.c_int64, 0)
    queue = Queue()
    processes = list()
    for i in range(threads):
        processes.append(Process(target = indexProcess, args = (queue, n, N)))
        processes[-1].start()
    for i in files:
        queue.put(i)
    for process in processes:
        queue.put(None)
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return None


def readIndex(file):
    sequences = list()
    lengths = list()
    fileOffsets = list()
    dataOffsets = list()
    dataSizes = list()
    gzipFile = GzipFile(filename = file, mode = 'rb', compresslevel = 9)
    marker = gzipFile.read(16)
    for i in range(struct1Q.unpack(gzipFile.read(8))[0]):
        sequences.append(gzipFile.read(struct1i.unpack(gzipFile.read(4))[0]).decode('ascii'))
        length, fileOffset, dataOffset, dataSize = struct1i3Q.unpack(gzipFile.read(28))
        lengths.append(length)
        fileOffsets.append(fileOffset)
        dataOffsets.append(dataOffset)
        dataSizes.append(dataSize)
    gzipFile.close()
    return (marker, sequences, lengths, fileOffsets, dataOffsets, dataSizes)


def readIndexProcess(queue, x, n, m, marker):
    '''
    n: #sequences
    m: #files
    '''
    x = numpy.ndarray(shape = (3, n, m), dtype = numpy.int64, buffer = x)
    while True:
        i, file = queue.get()
        if i is None:
            break
        markerI, sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndex(file)
        assert markerI == marker, 'All bam files should have the same header.'
        x[0, : , i] = fileOffsets[ : -1]
        x[1, : , i] = dataOffsets[ : -1]
        x[2, : , i] = dataSizes[ : -1]
    return None


def readIndices(files, threads):
    m = len(files)
    marker, sequences, lengths, fileOffsets, dataOffsets, dataSizes = readIndex(f'{files[0]}.index')
    n = len(sequences) - 1
    X = RawArray(ctypes.c_int64, 3 * n * m)
    x = numpy.ndarray(shape = (3, n, m), dtype = numpy.int64, buffer = X)
    x[0, : , 0] = fileOffsets[ : -1]
    x[1, : , 0] = dataOffsets[ : -1]
    x[2, : , 0] = dataSizes[ : -1]
    queue = Queue()
    processes = list()
    for i in range(threads):
        processes.append(Process(target = readIndexProcess, args = (queue, X, n, m, marker)))
        processes[-1].start()
    for i, file in enumerate(files[1 : ], start = 1):
        queue.put((i, f'{file}.index'))
    for process in processes:
        queue.put((None, None))
    queue.close()
    queue.join_thread()
    for process in processes:
        process.join()
        process.close()
    return (sequences[ : -1], lengths[ : -1], x[0], x[1], x[2])


def getUngappedRegions(file, fileOffset, dataOffset, dataSize, minMAPQ, minIdentity):
    regions = list()
    for alignment, alignmentSize in readAlignment(file, fileOffset, dataOffset, dataSize):
        Position, readIDLength, mapq, nCigars, flag, readLength = struct1i2B2x2H1i.unpack_from(alignment, offset = 4)
        position = Position + 1
        if (flag & 0x104 == 0) and (mapq >= minMAPQ) and nCigars:
            offset = 32 + readIDLength # 32 + readIDLength #
            cigarOffset = 4 * nCigars
            alignedLength = 0
            totalLength = 0
            for cigar in struct1I.iter_unpack(alignment[offset : offset + cigarOffset]):
                operations = cigar[0] >> 4
                operation = cigar[0] & 0x0F # MIDNSHP=X #
                if operation in {0, 7, 8}: # Consume reference. #
                    alignedLength += operations
                    regions.append((position, position + operations - 1))
                    position += operations
                elif operation in {2, 3}: # Consume reference. #
                    position += operations
                totalLength += operations
            if alignedLength >= minIdentity * totalLength:
                yield (Position + 1, regions)
        regions.clear()
    return None
