# makeHalfmassImages.py
# Joseph Wick
#
# Sanity check for half mass radii  of a single galaxy

import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

import pynbody

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

def makeHalfmassImg(gal, ts, hmr, width=20):
    print('Timestep: '+str(ts)+'...', end='')

    # get sim file
    baseDir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'
    galDir = baseDir+ '/r' + str(gal)+'.romulus25.3072g1HsbBH'

    f = glob.glob(galDir+'/r*.romulus25.3072g1HsbBH.00'+str(ts)+'/r*.00'+str(ts))[0]

    s = pynbody.load(f)
    s.physical_units()

    # center
    hCDM = -1
    try:
        hCDM = s.halos(write_fpos=False)[1]
    except:
        # center manually if missing halo; taken from pynbody source code
        mtot = s.s['mass'].sum()
        cen = np.sum(s.s['mass'] * s.s['pos'].transpose(), axis=1) / mtot
        cen.units = s.s['pos'].units
        s['pos'] -= cen
    else:
        hCDM = s.halos(write_fpos=False)[1]
        cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
        s['pos'] -= cen_pot

        pynbody.analysis.angmom.faceon(hCDM)

    # do image
    fig, axs = plt.subplots(1,3, figsize=(18,6), facecolor='w')

    c1 = plt.Circle((0,0), hmr, edgecolor='g', linewidth=3, fill=False)
    c2 = plt.Circle((0,0), hmr, edgecolor='g', linewidth=3, fill=False)
    c3 = plt.Circle((0,0), hmr, edgecolor='g', linewidth=3, fill=False)

    axs[0].scatter(s.s['pos'][:,0], s.s['pos'][:,1], s=1)
    axs[1].scatter(s.s['pos'][:,0], s.s['pos'][:,2], s=1)
    axs[2].scatter(s.s['pos'][:,1], s.s['pos'][:,2], s=1)

    axs[0].add_patch(c1)
    axs[1].add_patch(c2)
    axs[2].add_patch(c3)

    axs[1].set_title('r'+str(gal)+': '+str(ts), fontsize=20, fontfamily='serif')
    axs[2].set_title('hmr')

    lim = 25
    for ax in axs:
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])

    plt.tight_layout()

    plt.savefig('/home/jw1624/H1-merian/figures/breathingModes/hmrCheck/r'+str(gal)+'/'+str(ts)+'.png')
    print('done')
# end func
##

##
# main
##

# handle command line argument
if len(sys.argv) != 2 and len(sys.argv) != 3:
    print('Usage: python3 makeHalfmassImages.py [galaxy idx]')
    sys.exit
idx = int(sys.argv[1])

# get galaxies and iterate
galIDs = util.getGalaxies()[0]
gal = galIDs[idx]

# load in list of timesteps and radii
csvf = '/home/jw1624/H1-merian/csvs/breathingModes_clean/r'+str(gal)+'_qtys.txt'
df = pd.read_csv(csvf)

timesteps = df['timestep']
radii = df['R_halfmass']

if len(sys.argv) == 2:
    for i,ts in enumerate(timesteps):
        if i<25:
            makeHalfmassImg(gal, ts, radii[i])
else:
    tsidx = int(sys.argv[2])
    makeHalfmassImg(gal, timesteps[tsidx], radii[tsidx])
