# makeEinastoFit.py
# Joseph Wick
# Plots DM profile and Einasto fit with alpha

import pynbody
import numpy as np

import scipy.optimize as opt

import glob

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams

import sys
sys.path.insert(0, '../util/')
from util import util

##

##
# Function Definitions

# Einastro profile function for fitting
def einasto(r, A, alpha, c):
    return c*np.exp(-A*(r**alpha))
#

# Returns: radial bins, DM density per bin
def getDMProfile(sCDM):
    pynbody.analysis.angmom.faceon(hCDM)

    # profile range
    pmin = '0.0001 kpc'
    #pmax = pynbody.analysis.halo.virial_radius(h1) / 4
    pmax = '50 kpc'

    # rotation curve
    pdCDM = pynbody.analysis.profile.Profile(hCDM.d, rmin=pmin, rmax=pmax, type='lin', nbins=500)

    rbins = pdCDM['rbins']
    dmdensity = pdCDM['density']

    return rbins, dmdensity
#

# fits to einasto profile and returns einasto fit density
def getEinastoProfile(rbins, dmdensity):
    A, alpha, c = opt.curve_fit(einasto, rbins, dmdensity,
                                maxfev=1000000,
                                p0=[1, 0.5, max(dmdensity)])[0]

    return einasto(rbins, A, alpha, c), alpha
#

# sets rc params
def setPltParams():
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
#

##
# Main

# get galaxies
galIDs = util.getGalaxies()[0]

for gal in galIDs:
    cdmPath = util.getfilepath(gal)[0]
    cdmFile = cdmPath + '/r'+str(gal)+'.romulus25.3072g1HsbBH.004096'

    sCDM = pynbody.load(cdmFile)
    sCDM.physical_units()
    hCDM = sCDM.halos()[1]

    cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
    sCDM['pos'] -= cen_pot

    pynbody.analysis.angmom.faceon(hCDM)

    # get rbins and DMdensity
    rbins, dmdensity = getDMProfile(sCDM)
    # fit to einasto profile
    einastoP, alpha = getEinastoProfile(rbins, dmdensity)

    # do figure
    setPltParams()

    fig = plt.figure(figsize=(8,8), facecolor='w')

    plt.plot(rbins, dmdensity, linewidth=3, c='tab:blue')
    plt.plot(rbins, einastoP, linewidth=2, linestyle='--', c='tab:green')

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel('Radius', fontfamily='serif', fontsize=20)
    plt.ylabel('DM Density', fontfamily='serif', fontsize=20)
    plt.title('r'+str(gal), weight='bold', fontfamily='serif', fontsize=24)

    plt.text(s = r'$\alpha=$'+str(alpha)[0:4],
             x = 40, y = max(dmdensity))

    plt.savefig('/home/jw1624/H1-merian/figures/einastoFit/r'+str(gal)+'.png')
