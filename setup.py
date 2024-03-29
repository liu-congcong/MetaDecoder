#!/usr/bin/env python3
import os
from stat import S_IRWXU, S_IRGRP, S_IXGRP, S_IROTH, S_IXOTH
from setuptools import setup, find_packages
from setuptools.command.install import install


class set_permissions(install):
    def run(self): # Rewrite run function. #
        install.run(self) # Raw function. #
        '''
        distutils/command/install.py
            class install(Command)
                def get_outputs(self)
        '''
        for file_path in self.get_outputs():
            if os.path.basename(file_path) in ('fraggenescan', 'hmmsearch'):
                os.chmod(file_path, S_IRWXU + S_IRGRP + S_IXGRP + S_IROTH + S_IXOTH) # 0o755 #


def main():
    setup(
        name = 'metadecoder',
        version = '1.0.19',
        url = 'https://github.com/liu-congcong/metadecoder/',
        author = 'Liucongcong',
        author_email = 'congcong_liu@icloud.com',
        license = 'GPLv3',
        description = 'An algorithm for clustering metagenomic sequences.',
        install_requires = [
            'numpy',
            'scipy',
            'scikit-learn',
            'threadpoolctl'
        ],
        scripts = [
            'bin/metadecoder',
        ],
        packages = find_packages(),
        package_data = {
            '': ['LICENSE', 'markers.hmm', 'fraggenescan', 'hmmsearch', 'train/*'],
        },
        cmdclass = {
            'install': set_permissions
        }
    )


if __name__ == '__main__':
    main()
