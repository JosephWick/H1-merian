# getAxisRatios.py
# Joseph Wick
#
# finds axis ratios ba and ca for specified galaxy

import pynbody
import glob

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util_os import util_os
from util_galaxies import util_galaxies

##

gals = util_os.getGalaxies()[0]

for gal in gals:

    tsfile = util_os.getfilepath_cdm(gal, 0)

    sCDM = pynbody.load(tsfile)
    sCDM.physical_units()

    hCDM = sCDM.halos()[1]
    cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
    sCDM['pos'] -= cen_pot

    gasFrac = sum(hCDM.g['mass'])/(sum(hCDM.s['mass']) + sum(hCDM.g['mass']))

    print('r'+str(gal)+': '+str(gasFrac)[:6])
