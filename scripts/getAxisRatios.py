# getAxisRatios.py
# Joseph Wick
#
# finds axis ratios ba and ca for specified galaxy

import pynbody
import glob

import sys
sys.path.insert(0, '../util/')
from util import util

## command line arguments
# usage: python3 getAxisRatios galIdx tsIdx radius
if len(sys.argv) != 4:
    print('Usage: python3 getAxisRatios [galaxy idx] [ts idx] [galaxy radius]')
    sys.exit()

gal = util.getGalaxies()[0][int(sys.argv[1])]
tsidx = int(sys.argv[2])
R = float(sys.argv[3])

baseDir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'
galDir = baseDir+ '/r' + str(gal)+'.romulus25.3072g1HsbBH'

timesteps = glob.glob(galDir+'/r*.romulus25.3072g1HsbBH.0*')
timesteps.sort(reverse=True)

tsnum = timesteps[tsidx][-4:]
tsfile = glob.glob(timesteps[tsidx]+'/r*.00'+tsnum)[0]

sCDM = pynbody.load(tsfile)
sCDM.physical_units()

hCDM = -1
try:
    hCDM = sCDM.halos(write_fpos=False)[1]
except:
    # center manually if missing halo; taken from pynbody source code
    #print('HNF for halo ' + str(gal) + ', timestep '+str(tstepnumber))
    i = sCDM['phi'].argmin()
    cen_pot = sCDM['pos'][i].copy()
    sCDM['pos'] -= cen_pot
else:
    hCDM = sCDM.halos(write_fpos=False)[1]
    cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
    sCDM['pos'] -= cen_pot

rbin, ba, ca, angle, es = pynbody.analysis.halo.halo_shape(sCDM, rout=1.1*R, N=1, rin = 1*R)

print(gal)
print(rbin)
print(ba)
print(ca)
