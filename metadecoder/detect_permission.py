import os


def is_readable(path):
    '''
    Detect whether the path is readable.
    Parameters:
        path: the path to the file.
    Return:
        None
    '''
    assert os.access(path, os.R_OK), 'The path "{0}" is not readable.'.format(os.path.abspath(path))
    return None


def is_writeable(path):
    '''
    Detect whether the path is writeable.
    Parameters:
        path: the path to the file.
    Return:
        None
    '''
    try:
        open4w = open(path, 'wb')
        open4w.close()
        os.remove(path)
    except Exception:
        assert False, 'The path "{0}" is not writeable.'.format(os.path.abspath(path))
    return None