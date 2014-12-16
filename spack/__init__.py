"""
spack: a package for spherical packing analysis
"""
from .packing import Packing, rand_disk, rand_sphere

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
