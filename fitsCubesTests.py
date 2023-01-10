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

# let's load some data
dir = '/data/REPOSITORY/romulus_dwarf_zooms/r597.romulus25.3072g1HsbBH/'
f = 'r597.alfalfa_hi.x.fits'

hi_data = fits.open(dir+f)
print(hi_data[0].header['CUNIT3'])
cube = SpectralCube.read(hi_data)
