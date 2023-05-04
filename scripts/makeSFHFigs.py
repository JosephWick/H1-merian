import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util')
from util import util

import glob
import pynbody
import numpy as np
import matplotlib.pyplot as plt

def sfhFigs(hID, withSIDM=False):
    fig,axs = plt.subplots(1,2, figsize=(12,5))

    cdmC = 'firebrick'
    sidmC= 'royalblue'

    cdmPath, sidmPath, _ = util.getfilepath(hID)
    cdmFile = cdmPath + '/r'+str(hID)+'.romulus25.3072g1HsbBH.004096'

    sCDM = pynbody.load(cdmFile)
    sCDM.physical_units()

    sfhcdm1 = pynbody.plot.stars.sfh(sCDM, subplot=axs[0], color=cdmC)
    sfhcdm2 = pynbody.plot.stars.sfh(sCDM, subplot=axs[1], massform=False, color=cdmC)

    if withSIDM:
        sidmFile = sidmPath + '/r'+str(hID)+'.romulus25cvdXsec.3072g1HsbBH.004096'

        sSIDM = pynbody.load(sidmFile)
        sSIDM.physical_units()

        sfhsidm1 = pynbody.plot.stars.sfh(sSIDM, subplot=axs[0], color=sidmC)
        sfhsidm2 = pynbody.plot.stars.sfh(sSIDM, subplot=axs[1], massform=False, color=sidmC)

    # prettify
    axs[0].set_title('SFH with Original Star mass')
    axs[1].set_title('SFH with Final (Observed) Star Mass')

    axs[0].set_xlabel('Time [Gyr]')
    axs[0].set_ylabel(r'SFR $[M_{\odot} yr^{-1}]$')

    axs[1].set_xlabel('Time [Gyr]')
    axs[1].set_ylabel(r'SFR $[M_{\odot} yr^{-1}]$')

    axs[0].legend(['CDM', 'SIDM'])

    plt.tight_layout()

    if withSIDM:
        plt.savefig('/home/jw1624/H1-merian/figures/sfh/r'+str(hID)+'_SFH_2.png')
    else:
        plt.savefig('/home/jw1624/H1-merian/figures/sfh/r'+str(hID)+'_SFH.png')
#

cdmHalos, sidmHalos, adiabaticHalos = util.getGalaxies()

for g in cdmHalos:
    sfhFigs(g)
    if g in sidmHalos:
        sfhFigs(g, withSIDM=True)
    print('halo '+str(g)+' ... done')
