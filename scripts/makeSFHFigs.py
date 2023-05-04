import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util')
from util import util

import glob
import pynbody
import numpy as np
import matplotlib.pyplot as plt

def sfhFigs(hID, withSIDM=False):
    fig = plt.figure()

    cdmPath, sidmPath, _ = util.getfilepath(hID)
    cdmFile = cdmPath + '/r'+str(hID)+'.romulus25.3072g1HsbBH.004096'

    sCDM = pynbody.load(cdmFile)
    sCDM.physical_units()

    sfhcdm1 = pynbody.plot.stars.sfh(sCDM)

    if withSIDM:
        cdmFile = sidmPath + '/r'+str(hID)+'.romulus25cvdXsec.3072g1HsbBH.004096'

        sSIDM = pynbody.load(sidmFile)
        sSIDM.physical_units()

        sfhsidm1 = pynbody.plot.stars.sfh(sSIDM)

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
