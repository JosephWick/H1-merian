import pynbody
import glob

import numpy as np

import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util


# input: a simulation
# output: half mass radius of that simulation
def halfMassRadius(sim, startR, incBy, acc):
    pRadii = np.array(sim.s['r'])
    pMass = np.array(sim.s['mass'])
    mTot = sum(pMass)

    # start at 5 kpc and increment
    r = startR-incBy
    hm = mTot
    while(hm < (0.50-acc)*mTot or hm > (0.5+acc)*mTot):
        r += incBy
        hm = sum(pMass[pRadii < r])

        if (r>max(pRadii)):
            return -1

    return r

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
        fout.write('sigma_gas,sigma_star,sigma_youngstar,')
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
        a=glob.glob(timestep+'/*/*')
        if len(a)>0:
            # find sim in folder
            simFile = timestep+'/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tstepnumber
            simFile+= '/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tstepnumber
            #print(glob.glob(simFile))

        #print(glob.glob(simFile))
        # handle sim file not existing
        if len(glob.glob(simFile)) == 0:
            print('FNF for halo ' + str(gal) + ', timestep '+tstepnumber)
            continue

        # open simfile
        sCDM = pynbody.load(simFile)

        sCDM.physical_units()

        # SFR and sSFR (do first otherwise we sometimes get an error)
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
        rVir = pynbody.analysis.halo.virial_radius(sCDM)
        rHL = pynbody.analysis.luminosity.half_light_r(sCDM).in_units('kpc')
        rHL_c=pynbody.analysis.luminosity.half_light_r(sCDM, cylindrical=True).in_units('kpc')

        rHM = halfMassRadius(sCDM, 0, 1, 0.1)
        rHM = halfMassRadius(sCDM,rHM*0.85, 0.05, 0.01)

        # sSFR
        sSFR_10 = np.log10(SFR_10/(mStar*1e7))
        sSFR_100 = np.log10(SFR_100/(mStar*1e8))

        # velocity dispersion
        sCDM.properties['boxsize'] = 3.0e4

        # sigma HII
        #hiimask = sCDM.g['HII']>max(sCDM.g['HII'])*0.95
        #sigma_gas = np.median(sCDM.g['v_disp'][hiimask])

        sigma_star= np.median(sCDM.s['v_disp'])
        sigma_youngstar = np.median(sCDM.s['v_disp'][sCDM.s['age'].in_units('Myr') < 10])
        # line of sight sigma pred from Hirtenstein et al 2019 eqn (1)
        log_sigma_pred_10 = 0.1006*sSFR_10 + 0.3892*np.log10(mStar) + 0.0126*np.log10(mStar)*sSFR_10
        log_sigma_pred_100= 0.1006*sSFR_100+ 0.3892*np.log10(mStar) + 0.0126*np.log10(mStar)*sSFR_100

        # write to file
        fout.write(str(gal)+','+str(tstepnumber)+','+str(uage)+','+str(stepZ)+',')
        #fout.write(str(mStar)+','+str(rVir)+','+str(rHL)+','+str(rHM)+',')
        fout.write(str(mStar)+','+str(rVir)+','+str(rHL)+','+str(rHL_c)+','+str(rHM)+',')
        fout.write(str(sigma_gas)+','+str(sigma_star)+','+str(sigma_youngstar)+',')
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
