# parseFire.py
# Joseph Wick
# Parses a single Fire sim to make csv

import gizmo_analysis as gizmo
import utilities as ut

import pandas as pd
import numpy as np

import glob
import sys
sys.path.insert(0, '/home/jmwick/H1-merian/util/')
from util import util

##
# METHODS
##

#
# positions are 3D; must already be centered
def compute_Rhalfmass(positions, masses, startR, incBy, acc):
    radii = np.linalg.norm(correctedpos, ord=2, axis=1)
    mTot = np.sum(masses)

    # start at startR and increment
    r = startR - incBy
    hm = mTot
    while(hm < (0.50-acc)*mTot or hm > (0.5+acc)*mTot):
        r+=incBy
        hm = np.sum(masses[radii < r])

        # check bounds
        if (r>max(radii)):
            return -1

    return r

def compute_Rhalfmass_bisect(positions, masses, outerR, acc, maxiter=100000):
    innerLim = 0.0
    outerLim = outerR

    pRadii = np.linalg.norm(positions, ord=2, axis=1)
    mTot = sum(masses)

    r = outerR/2
    hm = sum(masses[pRadii < r])
    n=0
    while(hm < (0.50-acc)*mTot or hm > (0.5+acc)*mTot):
        if hm > 0.5*mTot: # too big, decrease r
            outerLim = r
            r = innerLim + (outerLim-innerLim)/2
        elif hm < 0.5*mTot: # too small, increase r
            innerLim = r
            r = innerLim + (outerLim-innerLim)/2
        hm = sum(masses[pRadii < r])

        n += 1
        if n>maxiter:
            return -1
            break

    return r

#
def makeFireCSV(gal):
  '''
  makeFireCSV()
    inputs
    ------
    gal : string - filepath to sims
  '''

  # set up output file
  outfile = '/home/jmwick/breathingModes/csvs/'+str(gal[-10:])+'.txt'

  # get redshifts
  zfile = gal+'/snapshot_times.txt'
  simtimes = pd.read_csv(zfile, sep='\s+', skiprows=6, header=None)

  # keys for snapshot_times
  st_tskey = 0
  st_zkey = 2
  st_tkey = 3
  st_lbkey = 4

  timesteps = np.flip(simtimes[st_tskey].to_numpy())
  redshifts = np.flip(simtimes[st_zkey].to_numpy())
  snaptimes = np.flip(simtimes[st_tkey].to_numpy())
  lookbacks = np.flip(simtimes[st_lbkey].to_numpy())

  # write header
  fout = open(outfile, 'w')
  fout.write('galaxyID,tstep,t,tlookback,z,')
  fout.write('M_star,R_halfmass,')
  fout.write('sigma_allstars,sigma_allstars_wtd,')
  fout.write('sigma_youngstar,sigma_youngstar_wtd,')
  fout.write('sigma_allgas,sigma_allgas_wtd,')
  fout.write('sigma_coldgas,sigma_coldgas_wtd,')
  fout.write('SFR_10,SFR_100,sSFR_10,sSFR_100\n')
  fout.close()

  for i,z in enumerate(redshifts):
    # load particles for timestep
    particles = gizmo.io.Read.read_snapshots('all', 'redshift', z, gal)

    # quantities used for deriving other quantities
    starmasses = particles['star'].prop('mass')
    starages = particles['star'].prop('age')
    starpos = particles['star'].prop('position')
    starvel = particles['star'].prop('host.velocity')

    agemask = starages<0.01

    # quantities for csv
    Mstar = np.sum(starmasses)

    # center of mass
    posDM = particles['dark'].prop('position')
    massDM= particles['dark'].prop('mass')
    com = np.sum(posDM * massDM[:, None], axis=0) / np.sum(massDM)

    posDM -= com
    starpos -= com

    rHM = compute_Rhalfmass_bisect(starpos-com,starmasses, 20000, 0.01)

    # velocity dispersions
    vel_youngstars = starvel[agemask]
    mass_youngstars = starmasses[agemask]

    sigma_star = util.compute_vdisp_std(starvel, starmasses, starvel)
    sigma_star_wtd = util.compute_vdisp_wtd(starvel, starmasses, starvel,
                        starmasses)

    sigma_youngstar = util.compute_vdisp_std(starvel, starmasses, vel_youngstars)
    sigma_youngstar_wtd = util.compute_vdisp_wtd(starvel, starmasses,
                            vel_youngstars, mass_youngstars)

    pos_allgas = particles['gas'].prop('position')
    pos_allgas -= com
    vel_allgas = particles['gas'].prop('host.velocity')
    mass_allgas = particles['gas'].prop('mass')

    sigma_allgas = util.compute_vdisp_std(vel_allgas, mass_allgas, vel_allgas)
    sigma_allgas = util.compute_vdisp_wtd(vel_allgas, mass_allgas, vel_allgas,
                    mass_allgas)

    sfr10 = np.sum(starmasses[starages<0.01])
    sfr100= np.sum(starmasses[starages<0.10])

    ssfr10 = 0
    if sfr10>0: ssfr10 = np.log10(sfr10/(Mstar*1e7))

    ssfr100 = 0
    if sfr100>0: ssfr100 = np.log10(sfr100/(Mstar*1e8))

    # write to csv
    fout = open(outfile, 'a')
    fout.write(gal[-10:]+','+str(timesteps[i])+','+str(snaptimes[i])+','
                +str(lookbacks[i])+',')
    fout.write(str(redshifts[i])+',')
    fout.write(str(Mstar)+','+str(rHM)+',')
    fout.write(str(sigma_star)+','+str(sigma_youngstar)+',')
    fout.write(str(sfr10)+','+str(sfr100)+','+str(ssfr10)+','+str(ssfr100)+'\n')
    fout.close()


##
# MAIN
##

# handle arguments
if len(sys.argv) != 2:
    print('Usage: python3 parseFire.py [galaxy idx]')
    sys.exit()
idx = int(sys.argv[1])

# get sim files
dpath = '/data/groups/leauthaud/yluo42/FIRE/'
simfiles = glob.glob(dpath+'/m*')
gal = simfiles[idx]

# do the needful
makeFireCSV(gal)
