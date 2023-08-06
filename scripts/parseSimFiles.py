import pynbody
import glob

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt

import os
import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util_os import util_os

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

# makeQAfig()
# makes QA figure that depicts starmask, dmmask and new Rhm
def makeQAfig(pos_allstars, pos_allDM, center, Rhm, Rstar, Rdm, hw, outdir):
    fig, axs = plt.subplots(1,3, figsize=(18,6))

    # make the figure
    axes = [0,1,2]
    idxX = [0,0,1]
    idxY = [1,2,2]
    for i,ax in enumerate(axes):
        # scatter the star and DM particles
        axs[ax].scatter(pos_allstars[:,idxX[i]], pos_allstars[:,idxY[i]], s=1)
        axs[ax].scatter(pos_allDM[:,idxX[i]], pos_allDM[:,idxY[i]], s=1, alpha = 0.01)

        # draw a circle around the star and DM particles we've selected
        c1 = plt.Circle((center[idxX[i]],center[idxY[i]]), Rstar, edgecolor='g',
                            linewidth=1, fill=False)
        c2 = plt.Circle((center[idxX[i]],center[idxY[i]]), Rdm, edgecolor='cyan',
                            linewidth=1, fill=False)
        c3 = plt.Circle((center[idxX[i]],center[idxY[i]]), Rhm, edgecolor='magenta',
                            linewidth=1, fill=False)
        axs[ax].add_patch(c1)
        axs[ax].add_patch(c2)
        axs[ax].add_patch(c3)

        # set view
        axs[ax].set_xlim([center[idxX[i]]-hw,center[idxX[i]]+hw])
        axs[ax].set_ylim([center[idxY[i]]-hw,center[idxY[i]]+hw])

        # label
        fs = 16
        axs[ax].set_xlabel('x', fontsize=fs)
        axs[ax].set_ylabel('y', fontsize=fs)

    # title and clean up
    plt.suptitle('r'+str(gal), weight='bold', fontsize=20)
    plt.tight_layout()

    # save
    plt.savefig(outdir)
    plt.close()

# does the parsing and creation of csvs
def makeGalQtyCSV(gal, doQA=False):
    numTS = util_os.getNumTimesteps(gal)

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

    # define QA directory
    QAdir = '/home/jw1624/H1-merian/QA/parse/r'+str(gal)

    simfileprev = util_os.getfilepath_cdm(gal,0)

    # open simfile
    sPrev = pynbody.load(simfileprev)
    sPrev.physical_units()

    hZero = sPrev.halos()[1]

    # center based on potential
    cen = pynbody.analysis.halo.center(hZero, mode='pot', retcen=True)

    hmrPrev = halfMassRadius_bisect(hZero.s['pos']-cen, hZero.s['mass'], 1000, 0.01)

    # iterate through each timestep
    for timestep in range(numTS):
        fout = open(outfile, 'a')

        # open simfile
        simFile = util_os.getfilepath_cdm(gal, timestep)
        sCDM = pynbody.load(simFile)
        sCDM.physical_units()

        # center by mass
        mtot = sCDM.s['mass'].sum()
        cen = np.sum(sCDM.s['mass'] * sCDM.s['pos'].transpose(),
                     axis=1) / mtot
        cen.units = sCDM.s['pos'].units

        # make cut based on particles within rfac*hmrPrev of current center
        rfac = 50
        dmfac = 2.0
        starmask = np.linalg.norm(sCDM.s['pos'] - cen, axis=1)<rfac*hmrPrev
        gasmask = np.linalg.norm(sCDM.g['pos'] - cen, axis=1)<=rfac*hmrPrev
        # larger radius for dm
        darkmask = np.linalg.norm(sCDM.d['pos'] - cen, axis=1)<=dmfac*rfac*hmrPrev

        # Mass
        mStar = sum(sCDM.s['mass'][starmask])

        # SFR
        SFR_10  = sum(sCDM.s['mass'][starmask][sCDM.s['age'][starmask].in_units('Myr')<10])
        SFR_100 = sum(sCDM.s['mass'][starmask][sCDM.s['age'][starmask].in_units('Myr')<100])

        # get age of universe
        uage = pynbody.analysis.cosmology.age(sCDM)
        # get redshift
        stepZ = pynbody.analysis.cosmology.redshift(sCDM, uage)

        # Sizes
        rVir = -1
        # NOTE: these are prolly wrong until we have halo catalog
        rHL = pynbody.analysis.luminosity.half_light_r(sCDM).in_units('kpc')
        rHL_c=pynbody.analysis.luminosity.half_light_r(sCDM,
                cylindrical=True).in_units('kpc')

        # this is accurate to our star cut
        rHM = halfMassRadius_bisect(sCDM.s['pos'][starmask]-cen, sCDM.s['mass'][starmask],
                                    20000, 0.01)

        # sSFR
        sSFR_10 = np.log10(SFR_10/(mStar*1e7))
        sSFR_100 = np.log10(SFR_100/(mStar*1e8))

        # velocity dispersion
        vel_allstars = sCDM.s['vel'][starmask]
        mass_allstars = sCDM.s['mass'][starmask]

        agemask = sCDM.s['age'][starmask].in_units('Myr')<10
        vel_youngstars = vel_allstars[agemask]
        mass_youngstars = mass_allstars[agemask]

        vdisp_allstars_uwtd = util_os.compute_vdisp_std(vel_allstars,
                                mass_allstars, vel_allstars)
        vdisp_allstars_wtd = util_os.compute_vdisp_wtd(vel_allstars,
                                mass_allstars, vel_allstars, mass_allstars)

        vdisp_youngstar_uwtd = util_os.compute_vdisp_std(vel_allstars,
                                mass_allstars, vel_youngstars)
        vdisp_youngstar_wtd = util_os.compute_vdisp_wtd(vel_allstars,
                                mass_allstars, vel_youngstars, mass_youngstars)

        # cold gas vdisp
        vel_allgas = sCDM.g['vel'][gasmask]
        mass_allgas = sCDM.g['mass'][gasmask]

        cgmask = sCDM.g['temp'][gasmask]<1000
        vel_coldgas = vel_allgas[cgmask]
        mass_coldgas= mass_allgas[cgmask]

        vdisp_allgas_uwtd = util_os.compute_vdisp_std(vel_allgas, mass_allgas,
                            vel_allgas)
        vdisp_allgas_wtd = util_os.compute_vdisp_wtd(vel_allgas, mass_allgas,
                            vel_allgas, mass_allgas)

        vdisp_coldgas_uwtd = util_os.compute_vdisp_std(vel_allgas, mass_allgas,
                                vel_coldgas)
        vdisp_coldgas_wtd = util_os.compute_vdisp_wtd(vel_allgas, mass_allgas,
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

        # create QA figur if desired
        if doQA:
            figfout = QAdir+'/'+str(tstepnumber)+'.png'
            hw = 600
            makeQAfig(sCDM.s['pos'], sCDM.d['pos'], cen, rHM,
                        hmrPrev*rfac, hmrPrev*rfac*dmfac, hw,
                        figfout)


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

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print('Usage: python3 parseSimFiles.py [galaxy id] [optional QA flag]')
    print('     [galaxy id] : integer id of zoomed galaxy; e.g. 431')
    print('     [optional QA flag] : 0 for no QA images, 1 for QA. Defaults to 0.')
    sys.exit()

gal = int(sys.argv[1])
currentGals = util_os.getGalaxies()[0]
if gal not in currentGals:
    print('Error: Galaxy ID not known.')
    sys.exit()

doQA = False
if len(sys.argv) == 3 and sys.argv[2] == '1':
    doQA = True

makeGalQtyCSV(gal, doQA)
