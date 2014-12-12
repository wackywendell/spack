"""
Spack
*****

A 2D or 3D sphere packing analysis package
------------------------------------------

`Code available`_ on Github.

`Full documentation`_ not yet available at Read the Docs. 

.. _Code available: https://github.com/wackywendell/spack

.. _Full documentation: https://spack.readthedocs.org

Description
-----------

Analysis of 2D or 3D spherical packings.

"""

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.sdist import sdist as _sdist
import os.path
import versioneer

versioneer.VCS = 'git'
versioneer.versionfile_source = 'spack/_version.py'
versioneer.versionfile_build = 'spack/_version.py'
versioneer.tag_prefix = 'v' # tags are like 1.2.0
versioneer.parentdir_prefix = 'spack-v' # dirname like 'myproject-1.2.0'
# create the extension and add it to the python distribution
setup(  name='spack', 
        author="Wendell Smith",
        author_email="wackywendell@gmail.com",
        description = ("A module for analyzing packings of 2D and 3D spheres"),
        license = "BSD",
        keywords = "spheres spherical packing jammed jamming",
        url = "https://spack.readthedocs.org",
        long_description=__doc__.strip(),
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering :: Mathematics",
            "Topic :: Scientific/Engineering :: Physics",
            "Intended Audience :: Science/Research",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: BSD License",
        ],
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        )

