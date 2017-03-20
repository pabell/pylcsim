from __future__ import print_function
from __future__ import absolute_import

"""
A python package to simulate X-ray lightcurves
from coherent signals and power spectrum models.
"""

__author__ = 'Riccardo Campana'
__version__ = '0.4.0'

from .simulation import Simulation
from .lcsinusoid import lcsinusoid
from .lcpsd import lcpsd
from .utils import poisson_randomization
from .utils import psd
from .utils import rebin
from .utils import logrebin
from .utils import saveFITSLC
from .utils import saveFITSPSD

__all__ = ['Simulation','lcsinusoid', 'lcpsd', 'poisson_randomization', 'psd', 'rebin', 'logrebin', 'saveFITSLC', 'saveFITSPSD']


