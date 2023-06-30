import gizmo_analysis as gizmo
import utilities as ut

import pynbody

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

import scipy.stats as ss

import pandas as pd
import numpy as np

import glob

import sys

##
# RC prams
matplotlib.rc('font', family='serif')
rcParams['axes.linewidth'] = 3
rcParams['xtick.major.width'] = 3
rcParams['xtick.major.size'] = 12
rcParams['ytick.major.width'] = 3
rcParams['ytick.major.size'] = 12

rcParams['font.family'] = 'serif'

rcParams['text.usetex'] = False

rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.transparent'] = False

##
# command line args
if len(sys.argv) != 2:
    print('Usage: python3 makeCoringFig.py [galaxy idx]')
    sys.exit()

idx = int(sys.argv[1])

datapath = '/home/jmwick/H1-merian/fire/csvs'
currentGals = glob.glob(datapath+'/*.txt')
currentGals.sort()
gal = currentGals[idx]

##
# Do figure
gal = currentGals[1]

zees = [0,1,2,3]
snapshots = [600, 534, 172, 120]
simdir = '/data/groups/leauthaud/yluo42/FIRE/'

pmin = '0.01 kpc'
pmax = '15.0 kpc'

fig = plt.figure(facecolor='w', figsize=(10,12))

# do the needful
for snap in snapshots:
    sgiz = pynbody.load(simdir+gal[-14:-4]+'/output/snapshot_'+str(snap)+'.hdf5')
    sgiz.physical_units()
    hgiz = sgiz.halos()[1]

    # center
    cen_pot = pynbody.analysis.halo.center(hgiz, mode='pot', retcen=True)
    sgiz['pos'] -= cen_pot

    pynbody.analysis.angmom.faceon(hgiz)

    # profile dark matter density
    pdDM = pynbody.analysis.profile.Profile(sgiz.d, rmin=pmin, rmax=pmax, type='lin')

    plt.plot(pdDM['rbins'], np.log10(pdDM['density']))

# prettify
plt.title(gal[-14:-4], fontsize=18)

plt.xlabel('log Radius [kpc]', fontsize=14)
plt.ylabel(r'log $\rho \; (M_\odot/kpc^3)$', fontsize=14)

plt.xlim([0,1.2])

plt.savefig('/home/jmwick/H1-merian/fire/figures/coring/'+gal[-14:-4]+'.png')
