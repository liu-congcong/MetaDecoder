import os
import platform

hash = {
    ('Windows', 'x86'): '.dll',
    ('Windows', 'AMD64'): '.dll',
    ('Linux', 'x86_64'): '.so',
    ('Darwin', 'arm64'): '.dylib'
}

def findCPath(x):
    extension = hash.get((platform.system(), platform.machine()))
    assert extension is not None, f'Unsupported platform: {platform.system()} {platform.machine()}.'
    return os.path.join(os.path.dirname(__file__), x + extension)
