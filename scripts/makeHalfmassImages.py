# makeHalfmassImages.py
# Joseph Wick
#
# Sanity check for half mass radii

import pandas as pd
import glob
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

def makeHalfmassImg(gal, ts, r, width=20):
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

        pynbody.analysis.angmom.faceon(hCDM)

    # do image
    im = pynbody.plot.stars.render(s, width=str(width)+' kpc', plot=False, ret_im=True)

    pxlwid = len(im[0])
    pxlKpcRatio = pxlwid/width

    plt.title('r'+str(hID), fontsize=20, fontfamily='serif')
    # find radius and plot as circle
    hmr = halfMassRadius_bisect(s, width, 0.01)

    circle = plt.Circle((pxlwid/2,pxlwid/2), hmr*pxlKpcRatio, edgecolor='g', linewidth=3, fill=False)
    ax.imshow(im)
    ax.add_patch(circle)
    plt.tight_layout()

    plt.savefig('/home/jw1624/H1-merian/figures/hmrCheck/r'+str(gal)+'/'+str(ts)+'.png')

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
        makeHalfmassImg(gal, ts, radii[i])
        break
    break
