import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

import pynbody
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

# flux density
def getFluxDensity(hi, dv, D):
    S21 = []
    for m in hi:
        s = (1/(2.36e5)) * (m/(1)) * (1/dv)*((1/D)**2)
        S21.append(s)
    return np.array(S21)

def profileHI(hID, withSIDM=False):
    # four panel figure, 1x4
    # rot curve; gas surface density; gas dispersion; HI

    # flux density params
    D = 70 #Mpc
    dv= 11.2 #km/s

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
    s_r = sCDM.g['r']
    s_vd = sCDM.g['v_disp']

    s_r_int = np.array(s_r, dtype=int)
    sigmaOfR = np.zeros(max(s_r_int))
    numPerR = np.zeros(max(s_r_int))

    for idx,r in enumerate(s_r_int):
        sigmaOfR[r-1] += s_vd[r-1]
        numPerR[r-1] += 1

    sigmaOfR = sigmaOfR/numPerR

    rxaxisFiltered = sigmaOfR[sigmaOfR>0]
    sigmaFiltered = sigmaOfR[sigmaOfR>0]
    vDispMed = np.median(sigmaFiltered[0:15000])

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

    axs[2].plot([0,15], [vDispMed,vDispMed], c=cdmC, linewidth=lw)
    axs[2].set_title('median gas dispersion')
    axs[2].set_xlabel('radius [kpc]')
    axs[2].set_ylabel('velocity [km/s]')
    axs[2].set_xlim([0,15])
    axs[2].set_ylim([0,100])

    # do HI
    # using 'y' orientation for all gals, as x and z sometimes are poor
    dpath = '/home/jw1624/H1-merian/h1lines/'
    fcdm = glob.glob(dpath+'r'+str(hID)+'_cdm*')
    fcdm.sort()

    cdmx = pd.read_csv(fcdm[1], sep='\s+', header=None)[0]
    cdmy = pd.read_csv(fcdm[1], sep='\s+', header=None)[1]
    cdmy = getFluxDensity(cdmy, dv,D)

    axs[3].plot(cdmx, cdmy, c=cdmC, linewidth=lw)
    axs[3].set_title('HI profile')
    axs[3].set_xlabel('velocity [km/s]')
    axs[3].set_ylabel('flux density [Jy]')

    if withSIDM:
        sidmFile = sidmPath + '/r'+str(hID)+'.romulus25cvdXsec.3072g1HsbBH.004096'

        sSIDM = pynbody.load(sidmFile)
        sSIDM.physical_units()
        hSIDM = sSIDM.halos()[1]

        # center based on potential
        cen_pot = pynbody.analysis.halo.center(hSIDM, mode='pot', retcen=True)
        sSIDM['pos'] -= cen_pot

        # gas dispersion
        s_r = sSIDM.g['r']
        s_vd = sSIDM.g['v_disp']

        s_r_int = np.array(s_r, dtype=int)
        sigmaOfR = np.zeros(max(s_r_int))
        numPerR = np.zeros(max(s_r_int))

        for idx,r in enumerate(s_r_int):
            sigmaOfR[r-1] += s_vd[r-1]
            numPerR[r-1] += 1

        sigmaOfR = sigmaOfR/numPerR

        rxaxisFiltered = sigmaOfR[sigmaOfR>0]
        sigmaFiltered = sigmaOfR[sigmaOfR>0]
        vDispMed = np.median(sigmaFiltered[0:15000])

        pynbody.analysis.angmom.faceon(hSIDM)

        # rotation curve (global?)
        pSIDM = pynbody.analysis.profile.Profile(hSIDM, rmin=pmin, rmax=pmax, type='lin')

        # gas surface density
        sigma = pSIDM['mass']/(12.57*pSIDM['rbins']*pSIDM['rbins'])

        # plot first three panels
        axs[0].plot(pSIDM['rbins'], pSIDM['v_circ'], c=sidmC, linewidth=lw)
        axs[1].plot(pSIDM['rbins'], sigma, c=sidmC, linewidth=lw)
        axs[2].plot([0,15], [vDispMed,vDispMed], c=sidmC, linewidth=lw)

        # do HI
        # using 'y' orientation for all gals, as x and z sometimes are poor
        dpath = '/home/jw1624/H1-merian/h1lines/'
        fsidm = glob.glob(dpath+'r'+str(hID)+'_sidm*')
        fsidm.sort()

        sidmx = pd.read_csv(fsidm[1], sep='\s+', header=None)[0]
        sidmy = pd.read_csv(fsidm[1], sep='\s+', header=None)[1]
        sidmy = getFluxDensity(sidmy, dv,D)

        axs[3].plot(sidmx, sidmy, c=sidmC, linewidth=lw)

    fig.tight_layout()
    if withSIDM:
        plt.savefig('/home/jw1624/H1-merian/figures/toyModel/r'+str(hID)+'_4panel_2.png')
    else:
        plt.savefig('/home/jw1624/H1-merian/figures/toyModel/r'+str(hID)+'_4panel.png')
    # end

# get galaxies
cdmHalos,sidmHalos,_ = util.getGalaxies()

for g in cdmHalos:
    if g in sidmHalos:
        profileHI(g, withSIDM=True)
    profileHI(g, withSIDM=False)
    print('halo '+str(g)+'... done')
