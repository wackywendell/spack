"""
spack
*****

`Code available`_ on Github, and `documentation`_ at Read the Docs. 

.. _Code available: https://github.com/wackywendell/spack

.. _documentation: https://spack.readthedocs.org

A 2D or 3D sphere packing analysis package
------------------------------------------

This package exists to enable fast, simple, and easy analysis of packings of spheres (3D) and
disks (2D). 
"""

from setuptools import setup
from codecs import open # To use a consistent encoding
from os import path
import os
del os.link
import versioneer

here = path.abspath(path.dirname(__file__))
# Get the long description from the relevant file
versioneer.VCS = 'git'
versioneer.versionfile_source = 'spack/_version.py'
versioneer.versionfile_build = 'spack/_version.py'
versioneer.tag_prefix = 'v' # tags are like 1.2.0
versioneer.parentdir_prefix = 'spack-v' # dirname like 'myproject-1.2.0'


setup(  name='spack', 
        version=versioneer.get_version(),
        
        description = ("A module for analyzing packings of 2D and 3D spheres"),
        long_description=__doc__.strip(),
        
        url = "https://spack.readthedocs.org",
        
        author="Wendell Smith",
        author_email="wackywendell@gmail.com",
        
        license = "BSD",
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
        
        keywords = "spheres spherical packing jammed jamming",
        packages = ['spack'],
        
        cmdclass=versioneer.get_cmdclass(),
        install_requires = [
            'numpy',
        ],
        extras_require = dict(
            voro=['tess'],
            scene=['vapory', 'matplotlib'],
            trees=['pyparm'],
        ),
        
        # If there are data files included in your packages that need to be
        # installed, specify them here. If using Python 2.6 or less, then these
        # have to be included in MANIFEST.in as well.
        #package_data={
        #   'sample': ['package_data.dat'],
        #},
        # Although 'package_data' is the preferred approach, in some case you may
        # need to place data files outside of your packages.
        # see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
        # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
        #data_files=[('my_data', ['data/data_file'])],
        # To provide executable scripts, use entry points in preference to the
        # "scripts" keyword. Entry points provide cross-platform support and allow
        # pip to create the appropriate form of executable for the target platform.
        #entry_points={
        #'console_scripts': [
        #   'sample=sample:main',
        #],
        
)

