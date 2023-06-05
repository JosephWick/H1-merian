import pynbody
import glob

import numpy as np

import sys
sys.path.insert(0, '../util/')
from util import util

def makeGalQtyCSV(gal):
    baseDir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'
    galDir = baseDir+ '/r' + str(gal)+'.romulus25.3072g1HsbBH'

    timesteps = glob.glob(galDir+'/r*.romulus25.3072g1HsbBH.0*')
    timesteps.sort()

    # setup file for exporting values
    outfile = '/home/jw1624/H1-merian/csvs/breathingModes/r'+str(gal)+'_qtys.txt'
    fout = open(outfile,'w')
    fout.write('galaxyID,timestep,t,z,')
    fout.write('M_star,R_vir,R_halflight,R_halfmass,')
    fout.write('sigma_gas,sigma_star,')
    fout.write('SFR_10,SFR_100,sSFR_10,sSFR_100\n')

    # iterate through each timestep
    for timestep in timesteps:
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
        #rVir = pynbody.analysis.halo.virial_radius(sCDM)
        rHL = pynbody.analysis.luminosity.half_light_r(sCDM)

        rHM = halfMassRadius(sCDM, 0, 1, 0.1)
        rHM = halfMassRadius(sCDM,rHM*0.85, 0.05, 0.01)

        # velocity dispersion
        sCDM.properties['boxsize'] = 3.0e4
        sigma_gas = np.median(sCDM.g['v_disp'])
        sigma_star= np.median(sCDM.s['v_disp'])

        # sSFR
        sSFR_10 = SFR_10/mStar
        sSFR_100 = SFR_100/mStar

        # write to file
        fout.write(str(gal)+','+str(tstepnumber)+','+str(uage)+','+str(stepZ)+',')
        #fout.write(str(mStar)+','+str(rVir)+','+str(rHL)+','+str(rHM)+',')
        fout.write(str(mStar)+','+str(rHL)+','+str(rHM)+',')
        fout.write(str(sigma_gas)+','+str(sigma_star)+',')
        fout.write(str(SFR_10)+','+str(SFR_100)+','+str(sSFR_10)+','+str(sSFR_100)+'\n')

        print('r'+str(gal)+' '+str(tstepnumber)+' done')

    fout.close()
##

if len(sys.argv) != 2:
    print('Usage: python3 parseSimFiles.py [galaxy idx]')
    sys.exit()

idx = sys.argv[1]
currentGals = util.getGalaxies()[0]
gal = currentGals[idx]

makeGalQtyCSV(gal)
