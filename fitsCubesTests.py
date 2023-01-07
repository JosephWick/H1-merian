# fitsCubesTests.py
# Joseph Wick

import sys
sys.path.insert(0, '/home/jw1624/miniconda3/pkgs')

import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.utils import data

import matplotlib.pyplot as plt

from spectral_cube import SpectralCube

from astroquery.esasky import ESASky
from astroquery.utils import TableList
from astropy.wcs import WCS
from reproject import reproject_interp
