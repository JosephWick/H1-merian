# makeHalfmassImages.py
# Joseph Wick
#
# Sanity check for half mass radii

import pandas as pd
import glob
import matplotlib.pyplot as plt

import pynbody

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

def makeHalfmassImg(gal, ts, hmr, width=20):
    # get sim file
    baseDir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'
    galDir = baseDir+ '/r' + str(gal)+'.romulus25.3072g1HsbBH'

    f = glob.glob(galDir+'/r*.romulus25.3072g1HsbBH.00'+str(ts)+'/r*.00'+str(ts))[0]

    s = pynbody.load(f)
    s.physical_units()

    h = s.halos()

    # center
    hCDM = -1
    try:
        hCDM = s.halos(write_fpos=False)[1]
    except:
        # center manually if missing halo; taken from pynbody source code
        #print('HNF for halo ' + str(gal) + ', timestep '+str(tstepnumber))
        i = s['phi'].argmin()
        cen_pot = s['pos'][i].copy()
        s['pos'] -= cen_pot
    else:
        hCDM = s.halos(write_fpos=False)[1]
        cen_pot = pynbody.analysis.halo.center(hCDM, mode='pot', retcen=True)
        s['pos'] -= cen_pot

        pynbody.analysis.angmom.faceon(hCDM)

    # do image
    fig, ax = plt.subplots(1,1, figsize=(8,8), facecolor='w')
    im = pynbody.plot.stars.render(s, width=str(width)+' kpc', plot=False, ret_im=True)

    pxlwid = len(im[0])
    pxlKpcRatio = pxlwid/width

    plt.title('r'+str(gal), fontsize=20, fontfamily='serif')

    circle = plt.Circle((pxlwid/2,pxlwid/2), hmr*pxlKpcRatio, edgecolor='g', linewidth=3, fill=False)
    ax.imshow(im)
    ax.add_patch(circle)
    plt.tight_layout()

    plt.savefig('/home/jw1624/H1-merian/figures/breathingModes/hmrCheck/r'+str(gal)+'/'+str(ts)+'.png')

# end func
##

##
# main
##

# get galaxies and iterate
galIDs = util.getGalaxies()[0]

for gal in galIDs:
    # load in list of timesteps and radii
    csvf = '/home/jw1624/H1-merian/csvs/breathingModes_clean/r'+str(gal)+'_qtys.txt'
    df = pd.read_csv(csvf)

    timesteps = df['timestep']
    radii = df['R_halfmass']

    for i,ts in enumerate(timesteps):
        if i<20:
            makeHalfmassImg(gal, ts, radii[i])
