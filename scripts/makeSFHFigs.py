import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util')
from util import util

import glob
import pynbody
import numpy as np
import matplotlib.pyplot as plt

def sfhFigs(hID, withSIDM=False):
    cdmPath, sidmPath, _ = util.getfilepath(hID)

    cdmFile = cdmPath + '/r'+str(hID)+'.romulus25.3072g1HsbBH.004096'

    sCDM = pynbody.load(cdmFile)
    sCDM.physical_units()
    hCDM = sCDM.halos()[1]

    sfh1 = pynbody.plot.stars(sCDM)
    plt.savefig('/home/jw1624/H1-merian/figures/sfh/r'+str(hID)+'_SFH.png')

#

cdmHalos, sidmHalos, adiabaticHalos = util.getGalaxies()

for g in cdmalos:
    sfhFigs(g)
    if g in sidmHalos:
        sfhFigs(g, withSIDM=True)
