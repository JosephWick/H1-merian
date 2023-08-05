import pynbody
import glob

import numpy as np

import scipy.optimize as opt

import os
import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util import util

##

# power law
def powerlaw(r, alpha, c):
    return c*(r**alpha)

def halfMassRadius_bisect(positions, masses, outerR, acc, maxiter=100000):
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
def makeGalQtyCSV(gal):
    baseDir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'
    galDir = baseDir+ '/r' + str(gal)+'.romulus25.3072g1HsbBH'

    timesteps = glob.glob(galDir+'/r*.romulus25.3072g1HsbBH.0*')
    timesteps.sort(reverse=True)

    # setup file for exporting values
    outfile = '/home/jw1624/H1-merian/csvs/breathingModes/r'+str(gal)+'_qtys.txt'
    fout = open(outfile,'w')
    fout.write('galaxyID,timestep,t,z,')
    fout.write('M_star,R_halflight_s,R_halflight_c,R_halfmass,')
    fout.write('sigma_allstars,sigma_allstars_wtd,')
    fout.write('sigma_youngstar,sigma_youngstar_wtd,')
    fout.write('sigma_allgas,sigma_allgas_wtd,')
    fout.write('sigma_coldgas,sigma_coldgas_wtd,')
    fout.write('log_sigma_pred_10,log_sigma_pred_100,')
    fout.write('alpha,')
    fout.write('SFR_10,SFR_100,sSFR_10,sSFR_100\n')
    fout.close()

    # get initial center and hmr
    simfileprev = timesteps[0]+'/r'+str(gal)+'.romulus25.3072g1HsbBH.004096'
    a=glob.glob(timesteps[0]+'/*')
    if len(a)>0:
        # find sim in folder
        simfileprev = timesteps[0]+'/r'+str(gal)+'.romulus25.3072g1HsbBH.004096'

    # open simfile
    sPrev = pynbody.load(simfileprev)
    sPrev.physical_units()

    hZero = sPrev.halos()[1]

    # center based on potential
    cen = pynbody.analysis.halo.center(hZero, mode='pot', retcen=True)
    sPrev['pos'] -= cen

    hmrPrev = halfMassRadius_bisect(hZero.s['pos'], hZero.s['mass'], 1000, 0.01)

    # iterate through each timestep
    for timestep in timesteps:
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

        # make cut based on particles within rfac*hmrPrev of prev timestep parsed
        rfac = 30
        dmfac = 1.5
        starmask = np.linalg.norm(sCDM.s['pos'].in_units('kpc') - cen,axis=1)<rfac*hmrPrev
        gasmask = np.linalg.norm(sCDM.g['pos'], axis=1)<=rfac*hmrPrev
        # larger radius for dm
        darkmask = np.linalg.norm(sCDM.d['pos'], axis=1)<=dmfac*rfac*hmrPrev

        # center by mass
        mtot = sCDM.s['mass'][starmask].sum()
        cen = np.sum(sCDM.s['mass'][starmask] * sCDM.s['pos'][starmask].transpose(),
                     axis=1) / mtot
        cen.units = sCDM.d['pos'].units
        sCDM['pos'] -= cen

        # SFR
        SFR_10  = sum(sCDM.s['mass'][starmask][sCDM.s['age'][starmask].in_units('Myr')<10])

        SFR_100 = sum(sCDM.s['mass'][starmask][sCDM.s['age'][starmask].in_units('Myr')<100])

        # get age of universe
        uage = pynbody.analysis.cosmology.age(sCDM)
        # get redshift
        stepZ = pynbody.analysis.cosmology.redshift(sCDM, uage)

        # Mass
        mStar = sum(sCDM.s['mass'][starmask])

        # Sizes
        rVir = -1
        # NOTE: these are prolly wrong until we have halo catalog
        rHL = pynbody.analysis.luminosity.half_light_r(sCDM).in_units('kpc')
        rHL_c=pynbody.analysis.luminosity.half_light_r(sCDM,
                cylindrical=True).in_units('kpc')

        # this is accurate to our star cut
        rHM = halfMassRadius_bisect(sCDM.s['pos'][starmask], sCDM.s['mass'][starmask],
                                    20000, 0.01)

        # sometimes the conversion to kpc doesn't work?
        #if rHL > 1000: rHL = rHL/1000
        #if rHL_c > 1000: rHL_c = rHL_c/1000
        #if rHM > 1000: rHM = rHM/1000

        # sSFR
        sSFR_10 = np.log10(SFR_10/(mStar*1e7))
        sSFR_100 = np.log10(SFR_100/(mStar*1e8))

        # velocity dispersion
        vel_allstars = sCDM.s['vel'][starmask]
        mass_allstars = sCDM.s['mass'][starmask]

        agemask = sCDM.s['age'][starmask].in_units('Myr')<10
        vel_youngstars = vel_allstars[agemask]
        mass_youngstars = mass_allstars[agemask]

        vdisp_allstars_uwtd = util.compute_vdisp_std(vel_allstars,
                                mass_allstars, vel_allstars)
        vdisp_allstars_wtd = util.compute_vdisp_wtd(vel_allstars,
                                mass_allstars, vel_allstars, mass_allstars)

        vdisp_youngstar_uwtd = util.compute_vdisp_std(vel_allstars,
                                mass_allstars, vel_youngstars)
        vdisp_youngstar_wtd = util.compute_vdisp_wtd(vel_allstars,
                                mass_allstars, vel_youngstars, mass_youngstars)

        # cold gas vdisp
        vel_allgas = sCDM.g['vel'][gasmask]
        mass_allgas = sCDM.g['mass'][gasmask]

        cgmask = sCDM.g['temp'][gasmask]<1000
        vel_coldgas = vel_allgas[cgmask]
        mass_coldgas= mass_allgas[cgmask]

        vdisp_allgas_uwtd = util.compute_vdisp_std(vel_allgas, mass_allgas,
                            vel_allgas)
        vdisp_allgas_wtd = util.compute_vdisp_wtd(vel_allgas, mass_allgas,
                            vel_allgas, mass_allgas)

        vdisp_coldgas_uwtd = util.compute_vdisp_std(vel_allgas, mass_allgas,
                                vel_coldgas)
        vdisp_coldgas_wtd = util.compute_vdisp_wtd(vel_allgas, mass_allgas,
                                vel_coldgas, mass_coldgas)

        # line of sight sigma pred from Hirtenstein et al 2019 eqn (1)
        log_sigma_pred_10 = 0.1006*sSFR_10 + 0.3892*np.log10(mStar)
        log_sigma_pred_10+= 0.0126*np.log10(mStar)*sSFR_10

        log_sigma_pred_100 = 0.1006*sSFR_100+ 0.3892*np.log10(mStar)
        log_sigma_pred_100+= .0126*np.log10(mStar)*sSFR_100

        # alpha from power fit
        # profile range
        pmin = '0.0001 kpc'
        #pmax = pynbody.analysis.halo.virial_radius(h1) / 4
        pmax = '50 kpc'

        # rotation curve
        pdCDM = pynbody.analysis.profile.Profile(sCDM.d[darkmask], rmin=pmin, rmax=pmax,
            type='lin', nbins=500)
        rbins = pdCDM['rbins']
        dmdensity = pdCDM['density']

        fitidx = 10
        alpha, c = opt.curve_fit(powerlaw, rbins[:fitidx], dmdensity[fitidx],
                                    maxfev=10000,
                                    p0 = [1,max(dmdensity)])[0]

        # write to file
        fout.write(str(gal)+','+str(tstepnumber)+','+str(uage)+','+str(stepZ)+',')
        fout.write(str(mStar)+','+str(rHL)+','+str(rHL_c)+','+str(rHM)+',')
        fout.write(str(vdisp_allstars_uwtd)+','+str(vdisp_allstars_wtd)+',')
        fout.write(str(vdisp_youngstar_uwtd)+','+str(vdisp_youngstar_wtd)+',')
        fout.write(str(vdisp_allgas_uwtd)+','+str(vdisp_allstars_wtd)+',')
        fout.write(str(vdisp_coldgas_uwtd)+','+str(vdisp_coldgas_wtd)+',')
        fout.write(str(log_sigma_pred_10)+','+str(log_sigma_pred_100)+',')
        fout.write(str(alpha)+',')
        fout.write(str(SFR_10)+','+str(SFR_100)+','+str(sSFR_10)+
                    ','+str(sSFR_100)+'\n')

        fout.close()

        # update relevant quantities for next ts
        hmrPrev = rHM

        print('r'+str(gal)+' '+str(tstepnumber)+' done')
##

if len(sys.argv) != 2:
    print('Usage: python3 parseSimFiles.py [galaxy idx]')
    sys.exit()

idx = int(sys.argv[1])
currentGals = util.getGalaxies()[0]
gal = currentGals[idx]

makeGalQtyCSV(gal)
