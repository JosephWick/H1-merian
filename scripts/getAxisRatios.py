# getAxisRatios.py
# Joseph Wick
#
# finds axis ratios ba and ca for specified galaxy

import pynbody
import glob

import sys
sys.path.insert(0, '../util/')
from util_os import util_os
from util_galaxies import util_galaxies

## command line arguments
# usage: python3 getAxisRatios galIdx tsIdx radius
if len(sys.argv) != 4:
    print('Usage: python3 getAxisRatios [galaxy] [galaxy radius]')
    sys.exit()

gal = int(sys.argv[1])
tsidx = int(sys.argv[2])
R = float(sys.argv[3])

tsfile = util_os.getfilepath_cdm(gal, 0)

sCDM = pynbody.load(tsfile)
sCDM.physical_units()

hCDM = sCDM.halos(write_fpos=False)[1]
cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
sCDM['pos'] -= cen_pot

rbin, ba, ca, angle, es = pynbody.analysis.halo.halo_shape(hCDM, rout=1.1*R, N=1, rin = 1*R)

print('r'+str(gal))
print('radius: ' + str(rbin[0]))
print('ba: ' + str(ba[0]))
print('ca: ' + str(ca[0]))
