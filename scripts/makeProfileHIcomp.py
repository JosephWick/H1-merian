import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

import pynbody
import numpy as np
import matplotlib.pyplot as plt
import glob

def profileHI(hID, withSIDM=False):
    # four panel figure, 1x4
    # rot curve; gas surface density; gas dispersion; HI

    cdmPath, sidmPath, _ = util.getfilepath(hID)
    cdmFile = cdmPath + '/r'+str(hID)+'.romulus25.3072g1HsbBH.004096'

    fig, axs = plt.subplots(1,4, figsize=(16,4), facecolor='white')
    cdmC = 'firebrick'
    sidmC = 'royalblue'
    adiabatC = 'g'
    lw = 3


    sCDM = pynbody.load(cdmFile)
    sCDM.physical_units()
    hCDM = sCDM.halos()[1]

    # center based on potential
    cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
    sCDM['pos'] -= cen_pot

    # gas dispersion
    vdisp = hCDM.g['v_disp']

    pynbody.analysis.angmom.faceon(hCDM)

    # profile range; based on El-Bhadry 2018 fig A1
    pmin = '0.01 kpc'
    pmax = '15 kpc'

    # rotation curve (global?)
    pCDM = pynbody.analysis.profile.Profile(hCDM, rmin=pmin, rmax=pmax, type='lin')

    # gas surface density
    sigma = pCDM['mass']/(12.57*pCDM['rbins']*pCDM['rbins'])

    # plot first three panels
    axs[0].plot(pCDM['rbins'], pCDM['v_circ'], c=cdmC, linewidth=lw)
    axs[0].set_title('rotation curve')
    axs[0].set_xlabel('radius [kpc]')
    axs[0].set_ylabel('velocity [km/s]')
    axs[0].set_xlim([0,15])

    axs[1].plot(pCDM['rbins'], sigma, c=cdmC, linewidth=lw)
    axs[1].set_title('gas surface density')
    axs[1].set_xlabel('radius [kpc]')
    axs[1].set_ylabel(r'$\Sigma$ (r) [$M_\odot$ kpc$^{-2}$]')
    axs[1].set_xlim([0,15])

    axs[2].plot([0,15], [np.median(vdisp),np.median(vdisp)], c=cdmC, linewidth=lw)
    axs[2].set_title('median gas dispersion')
    axs[2].set_xlabel('radius [kpc]')
    axs[2].set_ylabel('velocity [km/s]')
    axs[2].set_ylim([0,15])

    # do HI
    # using 'y' orientation for all gals, as x and z sometimes are poor
    dpath = '/home/jw1624/H1-merian/h1lines/'
    fcdm = glob.glob(dpath+'r'+str(hID)+'_cdm*')
    fcdm.sort()

    cdmx = pd.read_csv(fcdm[1], sep='\s+', header=None)[0]
    cdmy = pd.read_csv(fcdm[1], sep='\s+', header=None)[1]

    axs[3].plot(cdmx, cdmy, c=cdmC, linewidth=lw)
    axs[3].set_title('HI profile')
    axs[3].set_xlabel('velocity [km/s]')
    axs[3].set_ylabel(r'mass [M$_\odot$]')



    fig.tight_layout()
    plt.savefig('/home/jw1624/H1-merian/figures/toyModel/r'+str(hID)+'_4panel.png')
    # end

# get galaxies
cdmHalos,sidmHalos,_ = util.getGalaxies()

for g in cdmHalos:
    profileHI(g, withSIDM=False)
    break
