# parseFire.py
# Joseph Wick
# Parses a single Fire sim to make csv

import gizmo_analysis as gizmo
import utilities as ut

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

import glob
import sys
sys.path.insert(0, '/home/jmwick/H1-merian/util/')
from util_os import util_os
from util_galaxies import util_galaxies

##
# METHODS
##

#
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
  fout.write('sigma_allstars_global,sigma_allstars_los,')
  fout.write('sigma_youngstar_global,sigma_youngstar_los,')
  fout.write('sigma_allgas_global,sigma_allgas_los,')
  fout.write('sigma_coldgas_global,sigma_coldgas_los,')
  fout.write('alpha,')
  fout.write('SFR_10,SFR_100,sSFR_10,sSFR_100\n')
  fout.close()

  for i,z in enumerate(redshifts):
    # load particles for timestep
    particles = gizmo.io.Read.read_snapshots('all', 'redshift', z, gal)

    # quantities used for deriving other quantities
    mass_allstars = particles['star'].prop('mass')
    age_allstars = particles['star'].prop('age')
    pos_allstars = particles['star'].prop('position')
    vel_allstars = particles['star'].prop('host.velocity')

    # quantities for csv
    Mstar = np.sum(mass_allstars)

    # center of mass
    pos_DM = particles['dark'].prop('position')
    mass_DM= particles['dark'].prop('mass')
    com = np.sum(pos_DM * mass_DM[:, None], axis=0) / np.sum(mass_DM)

    # make a cut for all particles within 2000 kpc of dark matter CoM
    posmask = np.linalg.norm(pos_allstars-com, axis=1)<2000
    mass_allstars = mass_allstars[posmask]
    age_allstars = age_allstars[posmask]
    pos_allstars = pos_allstars[posmask]
    vel_allstars = vel_allstars[posmask]

    rHM = compute_Rhalfmass_bisect(pos_allstars-com,mass_allstars, 10000, 0.01)

    # velocity dispersions
    agemask = age_allstars<0.01
    pos_youngstars = pos_allstars[agemask]
    vel_youngstars = vel_allstars[agemask]
    mass_youngstars = mass_allstars[agemask]

    sigma_star_global = util_galaxies.compute_vdisp_global(vel_allstars, mass_allstars,
                            vel_allstars, mass_allstars)
    sigma_star_los = util_galaxies.compute_vdisp_los(vel_allstars, mass_allstars,
                            pos_allstars-com, vel_allstars, rHM, mass_allstars)

    sigma_youngstar_global = util_galaxies.compute_vdisp_global(vel_allstars, mass_allstars,
                            vel_youngstars, mass_youngstars)
    sigma_youngstar_los = util_galaxies.compute_vdisp_los(vel_allstars, mass_allstars,
                            pos_youngstars-com, vel_youngstars, rHM, mass_youngstars)

    pos_allgas = particles['gas'].prop('position')
    vel_allgas = particles['gas'].prop('host.velocity')
    mass_allgas = particles['gas'].prop('mass')

    # make that 2000kpc cut again
    posmask = np.linalg.norm(pos_allgas-com)<2000
    mass_allgas = mass_allgas[posmask]
    pos_allgas = pos_allgas[posmask]
    vel_allgas = vel_allgas[posmask]

    #sigma_allgas_global = util_galaxies.compute_vdisp_global(vel_allgas, mass_allgas,
    #                        vel_allgas, mass_allgas)
    #sigma_allgas_los = util_galaxies.compute_vdisp_los(vel_allgas, mass_allgas,
    #                        pos_allgas-com, vel_allgas, rHM, mass_allgas)

    sfr10 = np.sum(mass_allstars[age_allstars<0.01])
    sfr100= np.sum(mass_allstars[age_allstars<0.10])

    ssfr10 = 0
    if sfr10>0: ssfr10 = np.log10(sfr10/(Mstar*1e7))

    ssfr100 = 0
    if sfr100>0: ssfr100 = np.log10(sfr100/(Mstar*1e8))

    # can we fit alpha

    # first center DM positions
    pos_DM -= com

    # we want to fit DM density within inner 2 kpc
    # get radii
    radii_DM = np.linalg.norm(pos_DM, axis=1)
    # TODO

    # QA scatter the dm, stars
    QA=False
    if QA:
        fig, axs = plt.subplots(1,3, figsize=(15,5), facecolor='w')

        pos_allstars-=com

        axs[0].scatter(pos_DM[:,0], pos_DM[:,1], c='grey', alpha=0.01, s=1)
        axs[0].scatter(pos_allstars[:,0], pos_allstars[:,1], c='magenta', s=1)

        axs[1].scatter(pos_DM[:,0], pos_DM[:,2], c='grey', alpha=0.01, s=1)
        axs[1].scatter(pos_allstars[:,2], pos_allstars[:,1], c='magenta', s=1)

        axs[2].scatter(pos_DM[:,1], pos_DM[:,2], c='grey', alpha=0.01, s=1)
        axs[2].scatter(pos_allstars[:,1], pos_allstars[:,2], c='magenta', s=1)

        fpath = '/home/jmwick/H1-merian/QA/parse/fire/m11d/'+str(i)+'.png'
        plt.savefig(fpath)
        plt.close()

    # write to csv
    fout = open(outfile, 'a')
    fout.write(gal[-10:]+','+str(timesteps[i])+','+str(snaptimes[i])+','
                +str(lookbacks[i])+',')
    fout.write(str(redshifts[i])+',')
    fout.write(str(Mstar)+','+str(rHM)+',')
    fout.write(str(sigma_star_global)+','+str(sigma_star_los)+',')
    fout.write(str(sigma_youngstar_global)+','+str(sigma_youngstar_los)+',')
    fout.write(str(sigma_allgas_global)+','+str(sigma_allgas_los)+',')
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
