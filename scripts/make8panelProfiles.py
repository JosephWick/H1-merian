import sys
sys.path.insert(0, '../util/')
from util import util

import tangos
import pynbody

import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.utils import data
from astropy.wcs import WCS

import matplotlib.pyplot as plt

from spectral_cube import SpectralCube

from astroquery.esasky import ESASky
from astroquery.utils import TableList
from astropy.wcs import WCS
from reproject import reproject_interp

import glob

def eightPanelProfiles(hID, withSIDM=False):
    # get data paths for specified halo
    cdmPath, sidmPath = util.getfilepath(hID)

    cdmFile = cdmPath + '/r'+str(hID)+'.romulus25.3072g1HsbBH.004096'

    fig, axs = plt.subplots(2,4, figsize=(16,8), facecolor='white')

    sCDM = pynbody.load(cdmFile)
    sCDM.physical_units()
    hCDM = sCDM.halos()[1]
    sCDM.physical_units()
    hCDM = sCDM.halos()[1]

    # center based on potential
    cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
    sCDM['pos'] -= cen_pot

    # profile range
    pmin = '0.01 kpc'
    #pmax = pynbody.analysis.halo.virial_radius(h1) / 4
    pmax = '50 kpc'

    # rotation curve
    pdCDM = pynbody.analysis.profile.Profile(hCDM.d, rmin=pmin, rmax=pmax, type='lin')
    pgCDM = pynbody.analysis.profile.Profile(hCDM.g, rmin=pmin, rmax=pmax, type='lin')
    psCDM = pynbody.analysis.profile.Profile(hCDM.s, rmin=pmin, rmax=pmax, type='lin')
    pCDM  = pynbody.analysis.profile.Profile(hCDM,   rmin=pmin, rmax=pmax, type='lin')

    pdCDM2 = pynbody.analysis.profile.v_circ(pdCDM)

    cdmC = 'firebrick'
    sidmC = 'royalblue'
    lw = 3

    axs[0,0].plot(pdCDM['rbins'], pdCDM['density'], c=cdmC, linewidth=lw)
    axs[0,1].plot(pgCDM['rbins'], pgCDM['density'], c=cdmC, linewidth=lw)
    axs[1,0].plot(psCDM['rbins'], psCDM['density'], c=cdmC, linewidth=lw)
    axs[1,1].plot(pCDM['rbins'], pCDM['density'], c=cdmC, linewidth=lw)

    axs[0,0].set_xscale('log')
    axs[0,0].set_yscale('log')
    axs[0,1].set_xscale('log')
    axs[0,1].set_yscale('log')
    axs[1,0].set_xscale('log')
    axs[1,0].set_yscale('log')
    axs[1,1].set_xscale('log')
    axs[1,1].set_yscale('log')

    plt.suptitle('Profiles for Galaxy '+str(hID), fontsize=20)

    axs[0,0].set_title('Density of Dark Matter', fontsize=16)
    axs[0,1].set_title('Density of Gas', fontsize=16)
    axs[1,0].set_title('Density of Stars', fontsize=16)
    axs[1,1].set_title('Density of All', fontsize=16)

    axs[0,2].plot(pdCDM['rbins'], pdCDM['v_circ'], c=cdmC, linewidth=lw)
    axs[0,3].plot(pgCDM['rbins'], pgCDM['v_circ'], c=cdmC, linewidth=lw)
    axs[1,2].plot(psCDM['rbins'], psCDM['v_circ'], c=cdmC, linewidth=lw)
    axs[1,3].plot(pCDM['rbins'], pCDM['v_circ'], c=cdmC, linewidth=lw)

    axs[0,2].set_xscale('log')
    axs[0,2].set_yscale('log')
    axs[0,3].set_xscale('log')
    axs[0,3].set_yscale('log')
    axs[1,2].set_xscale('log')
    axs[1,2].set_yscale('log')
    axs[1,3].set_xscale('log')
    axs[1,3].set_yscale('log')

    axs[0,2].set_title('Rotation of Dark Matter', fontsize=16)
    axs[0,3].set_title('Rotation of Gas', fontsize=16)
    axs[1,2].set_title('Rotation of Stars', fontsize=16)
    axs[1,3].set_title('Rotation of All', fontsize=16)

    if withSIDM:
        sidmFile = sidmPath + '/r'+str(hID)+'.romulus25cvdXsec.3072g1HsbBH.004096'

        sSIDM = pynbody.load(sidmFile)
        sSIDM.physical_units()
        hSIDM = sSIDM.halos()[1]

        cen_pot = pynbody.analysis.halo.center(hSIDM, mode='pot', retcen=True)
        sSIDM['pos'] -= cen_pot

        pdSIDM = pynbody.analysis.profile.Profile(hSIDM.d, rmin=pmin, rmax=pmax, type='lin')
        pgSIDM = pynbody.analysis.profile.Profile(hSIDM.g, rmin=pmin, rmax=pmax, type='lin')
        psSIDM = pynbody.analysis.profile.Profile(hSIDM.s, rmin=pmin, rmax=pmax, type='lin')
        pSIDM  = pynbody.analysis.profile.Profile(hSIDM,   rmin=pmin, rmax=pmax, type='lin')

        axs[0,0].plot(pdSIDM['rbins'], pdSIDM['density'], c=sidmC, linewidth=lw)
        axs[0,1].plot(pgSIDM['rbins'], pgSIDM['density'], c=sidmC, linewidth=lw)
        axs[1,0].plot(psSIDM['rbins'], psSIDM['density'], c=sidmC, linewidth=lw)
        axs[1,1].plot(pSIDM['rbins'], pSIDM['density'], c=sidmC, linewidth=lw)

        pdSIDM2 = pynbody.analysis.profile.v_circ(pdSIDM)

        axs[0,2].plot(pdSIDM['rbins'], pdSIDM['v_circ'], c=sidmC, linewidth=lw)
        axs[0,3].plot(pgSIDM['rbins'], pgSIDM['v_circ'], c=sidmC, linewidth=lw)
        axs[1,2].plot(psSIDM['rbins'], psSIDM['v_circ'], c=sidmC, linewidth=lw)
        axs[1,3].plot(pSIDM['rbins'], pSIDM['v_circ'], c=sidmC, linewidth=lw)

        axs[0,0].legend(['CDM', 'SIDM'])

    fig.tight_layout()

    if withSIDM:
        plt.savefig('../figures/DenRotProfiles/r'+str(hID)+'_8panel_2.png')
    else:
        plt.savefig('../figures/DenRotProfiles/r'+str(hID)+'_8panel.png')
    # end 

# get haloIDs
cdmHalos, sidmHalos = util.getGalaxies()

print('Making 8 Panel Profiles')
# make figs 
for g in cdmHalos:
    print(' halo '+str(g)+'...', end='')
    if g in sidmHalos:
        eightPanelProfiles(g, withSIDM=True)
    eightPanelProfiles(g)
    print('done')


