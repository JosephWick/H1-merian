# makeHalfmassImages.py
# Joseph Wick
#
# Sanity check for half mass radii

import pandas as pd
import glob
import pyplot as plt

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

def makeHalfmassImg(gal, ts, r, width=20):
    # get sim file
    baseDir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'
    galDir = baseDir+ '/r' + str(hID)+'.romulus25.3072g1HsbBH'

    f = glob.glob(galDir+'/r*.romulus25.3072g1HsbBH.00'+str(ts)+'/r*.00'+str(ts))
    print(f)

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
