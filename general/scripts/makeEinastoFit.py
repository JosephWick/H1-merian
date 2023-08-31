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
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util_os import util_os

##
# Function Definitions

# Einastro profile function for fitting
def einasto(r, A, alpha, c, d):
    return c*np.exp(-A*(r**alpha)) + d
#

# power law
def powerlaw(r, alpha, c):
    return c*(r**alpha)

# Returns: radial bins, DM density per bin
def getDMProfile(sCDM):
    pynbody.analysis.angmom.faceon(hCDM)

    # profile range
    pmin = '0.0001 kpc'
    #pmax = pynbody.analysis.halo.virial_radius(h1) / 4
    pmax = '50 kpc'

    # rotation curve
    pdCDM = pynbody.analysis.profile.Profile(hCDM.d, rmin=pmin, rmax=pmax,
                type='lin', nbins=500)

    rbins = pdCDM['rbins']
    dmdensity = pdCDM['density']

    return rbins, dmdensity
#

# fits to einasto profile and returns einasto fit density
def getEinastoProfile(rbins, dmdensity):
    A, alpha, c, d = opt.curve_fit(einasto, rbins, dmdensity,
                                maxfev=1000000,
                                p0=[1, 0.5, max(dmdensity), min(dmdensity)])[0]

    return einasto(rbins, A, alpha, c, d), alpha
#

def fitPowerlaw(rbins, dmdensity):
    alpha, c = opt.curve_fit(powerlaw, rbins, dmdensity,
                                maxfev=10000,
                                p0 = [1,max(dmdensity)])[0]

    return powerlaw(rbins, alpha, c), alpha

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
galIDs = util_os.getGalaxies()[0]

for gal in galIDs:
    cdmFile = util_os.getfilepath_cdm(gal, 0)

    sCDM = pynbody.load(cdmFile)
    sCDM.physical_units()
    hCDM = sCDM.halos()[1]

    cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
    sCDM['pos'] -= cen_pot

    pynbody.analysis.angmom.faceon(hCDM)

    # get rbins and DMdensity
    rbins, dmdensity = getDMProfile(sCDM)

    # get slope of dmdensity
    sidx = 10

    # fit to einasto profile
    einastoP, alpha = getEinastoProfile(rbins[:sidx], dmdensity[:sidx])
    pl, alphap= fitPowerlaw(rbins[:sidx], dmdensity[:sidx])

    # do figure
    setPltParams()

    fig = plt.figure(figsize=(8,8), facecolor='w')

    plt.plot(rbins, dmdensity, linewidth=5, c='tab:blue')
    #plt.plot(rbins[:sidx], pl[:sidx], linewidth=3, linestyle='--', c='tab:green')

    plt.xscale('log')
    plt.yscale('log')

    #plt.xlabel('Radius', fontfamily='serif', fontsize=20)
    #plt.ylabel('DM Density', fontfamily='serif', fontsize=20)
    plt.title('r'+str(gal), weight='bold', fontfamily='serif', fontsize=42)

    plt.text(s = r'$\alpha=$'+str(alphap)[0:5],
             x = min(rbins), y = min(dmdensity),
             fontfamily='monospace', fontsize=30, weight='bold')

    plt.savefig('/home/jw1624/H1-merian/general/figures/slopeFit/r'+str(gal)+'.png')

    print(str(gal)+': done')
