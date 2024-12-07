import os
from struct import Struct
from zlib import MAX_WBITS, decompress


class BGZF:
    def __init__(self, file):
        self.file = file
        self.fileSize = os.path.getsize(file)
        self.temp = list()
        self.buffer8 = bytearray(8)
        self.buffer18 = bytearray(18)
        self.struct1I = Struct('<I') # uint32 #
        self.struct1H4x1H  = Struct('<H4xH') # uint16 4bits uint16 # #
        return None

    def __readBGZF__(self):
        '''
        ID1: uint8, ID2: uint8, CM: uint8, FLG: uint8, MTIME: uint32, XFL: uint8, OS: uint8, [10 bits]
        XLEN: uint16
        SI1: uint8, SI2: uint8, SLEN: uint16, [4 bits]
        BSIZE: uint16,
        xxx: XLEN-6 bytes, [XLEN-6 bits]
        CDATA: uint8[BSIZE-XLEN-19]
        CRC32: uint32, ISIZE: uint32
        '''
        self.stream.readinto(self.buffer18)
        xlen, bsize = self.struct1H4x1H.unpack_from(self.buffer18, 10)
        self.stream.seek(xlen - 6, 1)
        data = decompress(self.stream.read(bsize - xlen - 19), wbits = - MAX_WBITS, bufsize = 65536)
        self.stream.readinto(self.buffer8)
        isize = self.struct1I.unpack_from(self.buffer8, 4)[0]
        return (bsize + 1, data, isize)

    def open(self, fileOffset = 0):
        self.data = b''
        self.dataSize = 0
        self.dataOffset = 0
        self.__dataOffset = 0
        self.fileOffset = fileOffset
        self.__fileOffset = fileOffset
        self.stream = open(self.file, 'rb')
        self.stream.seek(fileOffset, 0)
        return self

    def read(self, n):
        if n <= self.dataSize:
            self.__dataOffset += n
            self.dataOffset += n
        else:
            if self.dataSize:
                self.temp.append(self.data[- self.dataSize : ])
            while self.dataSize < n:
                if self.__fileOffset == self.fileSize:
                    n = self.dataSize
                    dataSize = 0
                    break
                self.fileOffset = self.__fileOffset
                fileSize, data, dataSize = self.__readBGZF__()
                self.__fileOffset += fileSize
                self.temp.append(data)
                self.dataSize += dataSize
            self.dataOffset = dataSize - self.dataSize + n
            self.data = b''.join(self.temp)
            self.temp.clear()
            self.__dataOffset = n
        self.dataSize -= n
        return (self.data[self.__dataOffset - n : self.__dataOffset])

    def close(self):
        self.stream.close()
        return self
