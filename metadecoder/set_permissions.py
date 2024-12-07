import os
from stat import S_IRWXU, S_IRGRP, S_IXGRP, S_IROTH, S_IXOTH
from setuptools.command.install import install

class set_permissions(install):
    def run(self): # Rewrite run function. #
        super().run()
        for file_path in self.get_outputs():
            if os.path.basename(file_path) in ('fraggenescan', 'hmmsearch'):
                os.chmod(file_path, S_IRWXU + S_IRGRP + S_IXGRP + S_IROTH + S_IXOTH) # 0o755 #
