import pynbody
import glob

import numpy as np

import os
import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util


# input: a simulation
# output: half mass radius of that simulation
def halfMassRadius(sim, startR, incBy, acc):
    pRadii = np.array(sim.s['r'].in_units('kpc'))
    pMass = np.array(sim.s['mass'])
    mTot = sum(pMass)

    # start at starR kpc and increment
    r = startR-incBy
    hm = mTot
    while(hm < (0.50-acc)*mTot or hm > (0.5+acc)*mTot):
        r += incBy
        hm = sum(pMass[pRadii < r])

        if (r>max(pRadii)):
            return -1

    return r

def halfMassRadius_bisect(sim, outerR, acc, maxiter=100000):
    innerLim = 0.0
    outerLim = outerR

    pRadii = np.array(sim.s['r'].in_units('kpc'))
    pMass = np.array(sim.s['mass'])
    mTot = sum(pMass)

    r = outerR/2
    hm = sum(pMass[pRadii < r])
    n=0
    while(hm < (0.5-acc)*mTot or hm > (0.5+acc)*mTot):
        if hm > 0.5*mTot: # too big, decrease r
            outerLim = r
            r = innerLim + (outerLim-innerLim)/2
        elif hm < 0.5*mTot: # too small increase r
            innerLim = r
            r = innerLim + (outerLim-innerLim)/2
        hm = sum(pMass[pRadii < r])

        n += 1
        if n>maxiter:
            return -1
            break

    return r

# Yifei's code for mass weighted 3D velocity distribution
def compute_vdisp(velocities, masses):
    center_of_mass_velocity = np.sum(velocities * masses[:, None], axis=0) / np.sum(masses)
    velocity_differences = velocities - center_of_mass_velocity
    velocity_differences_squared = np.sum(velocity_differences**2, axis=1)
    mass_weighted_dispersion = np.sqrt(np.sum(velocity_differences_squared * masses) / np.sum(masses))
    return mass_weighted_dispersion

#
def makeGalQtyCSV(gal, startTS=0):
    baseDir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'
    galDir = baseDir+ '/r' + str(gal)+'.romulus25.3072g1HsbBH'

    timesteps = glob.glob(galDir+'/r*.romulus25.3072g1HsbBH.0*')
    timesteps.sort(reverse=True)

    # setup file for exporting values
    outfile = '/home/jw1624/H1-merian/csvs/breathingModes/r'+str(gal)+'_qtys.txt'
    if startTS == 0:
        fout = open(outfile,'w')
        fout.write('galaxyID,timestep,t,z,')
        fout.write('M_star,R_halflight_s,R_halflight_c,R_halfmass,')
        fout.write('sigma_youngstar,sigma_younstar_wtd,')
        fout.write('sigma_coldgas,sigma_coldgas_wtd,')
        fout.write('log_sigma_pred_10,log_sigma_pred_100,')
        fout.write('SFR_10,SFR_100,sSFR_10,sSFR_100\n')
        fout.close()

    # iterate through each timestep
    for timestep in timesteps[startTS:]:
        fout = open(outfile, 'a')
        tstepnumber = timestep[-6:]

        # I'm assuming that if there's another subfolder, the sim is inside it

        # try to find sim file in current folder
        simFile = timestep+'/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tstepnumber
        # check if there's another folder
        a=glob.glob(timestep+'/*')
        if len(a)>0:
            # find sim in folder
            simFile = timestep+'/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tstepnumber
            #simFile+= '/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tstepnumber
            #print(glob.glob(simFile))

        #print(glob.glob(simFile))
        # handle sim file not existing
        if len(glob.glob(simFile)) == 0:
            print('FNF for halo ' + str(gal) + ', timestep '+tstepnumber)
            continue

        # check for addtional folder
        if os.path.isdir(simFile):
            simFile = simFile+'/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tstepnumber
        # check if simFile exists
        if len(glob.glob(simFile))==0: continue


        # open simfile
        sCDM = pynbody.load(simFile)

        sCDM.physical_units()

        vdispg = sCDM.g['v_disp']

        # SFR (do first otherwise we sometimes get an error)
        SFR_10  = sum(sCDM.s['mass'][sCDM.s['age'].in_units('Myr')<10])
        SFR_100 = sum(sCDM.s['mass'][sCDM.s['age'].in_units('Myr')<100])

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

        # get age of universe
        uage = pynbody.analysis.cosmology.age(sCDM)
        # get redshift
        stepZ = pynbody.analysis.cosmology.redshift(sCDM, uage)

        # Mass
        mStar = sum(sCDM.s['mass'])

        # Sizes
        rVir = -1
        #rVir = pynbody.analysis.halo.virial_radius(sCDM)
        rHL = pynbody.analysis.luminosity.half_light_r(sCDM).in_units('kpc')
        rHL_c=pynbody.analysis.luminosity.half_light_r(sCDM, cylindrical=True).in_units('kpc')

        rHM = halfMassRadius_bisect(sCDM, 20000, 0.01)

        # sometimes the conversion to kpc doesn't work; let's do it manually if needed
        if rHL > 1000: rHL = rHL/1000
        if rHL_c > 1000: rHL_c = rHL_c/1000
        if rHM > 1000: rHM = rHM/1000

        # sSFR
        sSFR_10 = np.log10(SFR_10/(mStar*1e7))
        sSFR_100 = np.log10(SFR_100/(mStar*1e8))

        # velocity dispersion
        vel_allstars = sCDM.s['vel']
        mass_allstars = sCDM.s['mass']

        agemask = sCDM.s['age'].in_units('Myr')<10
        vel_youngstars = vel_allstars[agemask]
        mass_youngstars = mass_youngstars[agemask]

        vdisp_youngstar_uwtd = util.compute_vdisp_std(vel_allstars, mass_allstars, vel_youngstars)
        vdisp_youngstar_wtd = util.compute_vdisp_wtd(vel_allstars, mass_allstars, vel_youngstars, mass_youngstars)

        # cold gas vdisp
        cgmask = sCDM.g['temp']<1000
        coldgasVel = sCDM.g['vel'][cgmask]
        coldgasMass= sCDM.g['mass'][cgmask]
        vdisp_coldgas_uwtd = compute_vdisp_std(vel_allstars, mass_allstars, coldgasVel)
        vdisp_coldgas_wtd = compute_vdisp_wtd(vel_allstars, mass_allstars, coldgasVel, coldgasMass)

        # line of sight sigma pred from Hirtenstein et al 2019 eqn (1)
        log_sigma_pred_10 = 0.1006*sSFR_10 + 0.3892*np.log10(mStar) + 0.0126*np.log10(mStar)*sSFR_10
        log_sigma_pred_100= 0.1006*sSFR_100+ 0.3892*np.log10(mStar) + 0.0126*np.log10(mStar)*sSFR_100

        # write to file
        fout.write(str(gal)+','+str(tstepnumber)+','+str(uage)+','+str(stepZ)+',')
        fout.write(str(mStar)+','+str(rHL)+','+str(rHL_c)+','+str(rHM)+',')
        fout.write(str(vdisp_youngstar_uwtd)+','+str(vdisp_youngstar_wtd)+',')
        fout.write(str(vdisp_coldgas_uwtd)+','+str(vdisp_coldgas_wtd)+',')
        fout.write(str(log_sigma_pred_10)+','+str(log_sigma_pred_100)+',')
        fout.write(str(SFR_10)+','+str(SFR_100)+','+str(sSFR_10)+','+str(sSFR_100)+'\n')

        fout.close()
        print('r'+str(gal)+' '+str(tstepnumber)+' done')
##

if len(sys.argv) != 3:
    print('Usage: python3 parseSimFiles.py [galaxy idx] [ts idx]')
    sys.exit()

idx = int(sys.argv[1])
currentGals = util.getGalaxies()[0]
gal = currentGals[idx]

ts_idx = int(sys.argv[2])

makeGalQtyCSV(gal, ts_idx)
