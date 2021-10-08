from tempfile import mkstemp
import os


def make_file(prefix = None, suffix = None):
    '''
    Parameters:
        prefix: the prefix of the temp file [None].
        suffix: the suffix of the temp file [None].
    Return:
        the path to the temp file.
    '''
    file_descriptor, file = mkstemp(prefix = prefix, suffix = suffix, dir = os.getcwd())
    os.close(file_descriptor)
    return file