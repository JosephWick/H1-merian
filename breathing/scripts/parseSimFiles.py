import pynbody
import glob

import numpy as np

import scipy.optimize as opt

import matplotlib.pyplot as plt

import os
import sys
sys.path.insert(0, '/home/jw1624/H1-merian/util/')
from util_os import util_os
from util_galaxies import util_galaxies

# Convention: all coordinates are UNCENTERED in this file
# may be centered in a method or centered as arguments to a method

##

# power law
def powerlaw(r, alpha, c):
    return c*(r**alpha)

# makeQAfig_mask()
# makes QA figure that depicts starmask
def makeQAfig_mask(pos_allstars, pos_allDM, haloMask, center, Rhm, Rdm, hw, outdir):
    fig, axs = plt.subplots(1,3, figsize=(18,6))

    starRads = np.linalg.norm(pos_allstars-center, axis=1)
    starmask = starRads<=Rdm

    # make the figure
    xdims = [0,0,1]
    ydims = [1,2,2]
    for i,ax in enumerate(axs):
        # all DM
        ax.scatter(pos_allDM[:,xdims[i]], pos_allDM[:,ydims[i]], s=1,
            alpha = 0.01, c='grey')
        # halo DM
        ax.scatter(pos_allDM[:,xdims[i]][haloMask], pos_allDM[:,ydims[i]][haloMask], s=1,
            c='k')

        # all stars
        ax.scatter(pos_allstars[:,xdims[i]], pos_allstars[:,ydims[i]], s=1)

        # draw a circle around the star particles we've selected
        c1 = plt.Circle((center[xdims[i]],center[ydims[i]]), Rdm, edgecolor='cyan',
                            linewidth=1, fill=False)
        # circle around stellar half mass
        c2 = plt.Circle((center[xdims[i]],center[ydims[i]]), Rhm, edgecolor='r',
                            linewidth=1, fill=False)

        ax.add_patch(c1)
        ax.add_patch(c2)

        # set view
        ax.set_xlim([center[xdims[i]]-hw,center[xdims[i]]+hw])
        ax.set_ylim([center[ydims[i]]-hw,center[ydims[i]]+hw])

        # label
        fs = 16
        ax.set_xlabel('x', fontsize=fs)
        ax.set_ylabel('y', fontsize=fs)

    # title and clean up
    plt.suptitle('r'+str(gal), weight='bold', fontsize=20)
    plt.tight_layout()

    # save
    plt.savefig(outdir)
    plt.close()

def makeQAfig_profile(pd, outdir):
    plt.plot(pd['rbins'], pd['density'])

    plt.xscale('log')
    plt.yscale('log')

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
    fout.write('sigma_allstars_global,sigma_allstars_los,')
    fout.write('sigma_youngstar_global,sigma_youngstar_los,')
    fout.write('sigma_allgas_global,sigma_allgas_los,')
    fout.write('sigma_coldgas_global,sigma_coldgas_los,')
    fout.write('sigma_hotgas_global,sigma_hotgas_los,')
    fout.write('sigma_gasNearYS_global,sigma_gasNearYS_los,')
    fout.write('alpha,')
    fout.write('SFR_10,SFR_100,sSFR_10,sSFR_100\n')
    fout.close()

    # define QA directory
    QAdir = '/home/jw1624/H1-merian/QA/parse/romulus'

    simfileprev = util_os.getfilepath_cdm(gal,0)

    # open simfile
    sZero = pynbody.load(simfileprev)
    sZero.physical_units()

    hZero = sZero.halos()[1]

    # center based on potential
    cen = pynbody.analysis.halo.center(hZero, mode='pot', retcen=True)

    hmrPrev = util_galaxies.compute_massRadius(hZero.s['pos']-cen, hZero.s['mass'],
        1000, 0.01)

    # make DM mask based on halo from z=0 timestep
    xmask = np.nonzero(np.in1d(sZero.d['pos'][:,0], hZero.d['pos'][:,0]))
    ymask = np.nonzero(np.in1d(sZero.d['pos'][:,1], hZero.d['pos'][:,1]))
    zmask = np.nonzero(np.in1d(sZero.d['pos'][:,2], hZero.d['pos'][:,2]))

    DMmask = np.intersect1d(xmask, ymask)
    DMmask = np.intersect1d(DMmask, zmask)

    # iterate through each timestep
    for timestep in range(numTS):
        if gal == 431:
            if timestep in [1, 2, 3, 7, 11, 16, 18, 20, 21, 25, 31, 33]:
                continue
        fout = open(outfile, 'a')

        # open simfile
        simFile = util_os.getfilepath_cdm(gal, timestep)
        tstepnumber = simFile[-6:]

        # check if its a directoy :)
        if os.path.isdir(simFile):
            simFile = glob.glob(simFile+'/r'+str(gal)+'*.'+tstepnumber)
            if len(simFile) == 0:
                continue
            simFile = simFile[0]

        # double check simfile exists
        if len(glob.glob(simFile)) == 0:
            print('FNF: r'+str(gal)+' '+tstepnumber)
            continue

        sCDM = pynbody.load(simFile)
        sCDM.physical_units()

        hCDM = sCDM.d[DMmask]

        # center by mass of DM halo
        mtot = sCDM.d['mass'][DMmask].sum()
        cen = np.sum(sCDM.d['mass'][DMmask] * sCDM.d['pos'][DMmask].transpose(),
                     axis=1) / mtot
        cen.units = sCDM.d['pos'].units

        # make cut based on radius of masked DM particles
        rfac = 3
        hmrDM = util_galaxies.compute_massRadius(sCDM.d['pos'][DMmask],
            sCDM.d['mass'][DMmask], 1000, 0.1, frac=0.50, maxiter=1000)

        starmask = np.linalg.norm(sCDM.s['pos'] - cen, axis=1)<=rfac*hmrDM
        gasmask = np.linalg.norm(sCDM.g['pos'] - cen, axis=1)<=rfac*hmrDM

        # Mass
        mStar = np.sum(sCDM.s['mass'][starmask])

        # SFR
        SFR_10  = sum(sCDM.s['mass'][starmask][sCDM.s['age'][starmask].in_units('Myr')<10])
        SFR_100 = sum(sCDM.s['mass'][starmask][sCDM.s['age'][starmask].in_units('Myr')<100])

        # get age of universe
        uage = pynbody.analysis.cosmology.age(sCDM)
        # get redshift
        stepZ = pynbody.analysis.cosmology.redshift(sCDM, uage)

        # Sizes
        rVir = -1
        # NOTE: these are questionable until we have halo catalog
        rHL = pynbody.analysis.luminosity.half_light_r(sCDM).in_units('kpc')
        rHL_c=pynbody.analysis.luminosity.half_light_r(sCDM,
                cylindrical=True).in_units('kpc')

        # this is accurate to our star cut
        starcen = np.sum(sCDM.s['mass'][starmask] * sCDM.s['pos'][starmask].transpose(),
                         axis=1) / mStar
        rHM = util_galaxies.compute_massRadius(sCDM.s['pos'][starmask]-starcen,
            sCDM.s['mass'][starmask], 20000, 0.01)

        # sSFR
        sSFR_10 = np.log10(SFR_10/(mStar*1e7))
        sSFR_100 = np.log10(SFR_100/(mStar*1e8))

        # velocity dispersion
        pos_youngstars = sCDM.s['pos'][starmask]
        vel_allstars = sCDM.s['vel'][starmask]
        mass_allstars = sCDM.s['mass'][starmask]

        agemask = sCDM.s['age'][starmask].in_units('Myr')<10
        pos_youngstars = pos_youngstars[agemask]
        vel_youngstars = vel_allstars[agemask]
        mass_youngstars = mass_allstars[agemask]

        vdisp_allstars_global = util_galaxies.compute_vdisp_global(vel_allstars,
                                mass_allstars, vel_allstars, mass_allstars)
        vdisp_allstars_los = util_galaxies.compute_vdisp_los(vel_allstars, mass_allstars,
                            pos_youngstars-starcen, vel_youngstars, rHM, mass_youngstars)

        vdisp_youngstar_global = util_galaxies.compute_vdisp_global(vel_allstars,
                                mass_allstars, vel_youngstars, mass_youngstars)
        vdisp_youngstar_los = util_galaxies.compute_vdisp_los(vel_allstars, mass_allstars,
                            pos_youngstars-starcen, vel_youngstars, rHM, mass_youngstars)

        # cold gas vdisp
        pos_allgas = sCDM.g['pos'][gasmask]
        vel_allgas = sCDM.g['vel'][gasmask]
        mass_allgas = sCDM.g['mass'][gasmask]
        temp_allgas = sCDM.g['temp'][gasmask]

        cgmask = temp_allgas<1000
        pos_coldgas = pos_allgas[cgmask]
        vel_coldgas = vel_allgas[cgmask]
        mass_coldgas= mass_allgas[cgmask]

        vdisp_allgas_global = util_galaxies.compute_vdisp_global(vel_allgas, mass_allgas,
                            vel_allgas, mass_allgas)
        vdisp_allgas_los = util_galaxies.compute_vdisp_los(vel_allgas, mass_allgas,
                                pos_allgas-starcen, vel_allgas, rHM, mass_allgas)

        vdisp_coldgas_global = util_galaxies.compute_vdisp_global(vel_allgas, mass_allgas,
                                vel_coldgas, mass_coldgas)
        vdisp_coldgas_los = util_galaxies.compute_vdisp_los(vel_allgas, mass_allgas,
                                pos_coldgas-starcen, vel_coldgas, rHM, mass_coldgas)

        # hot gas vdisp
        hgmask = temp_allgas>10000
        pos_hotgas = pos_allgas[hgmask]
        vel_hotgas = vel_allgas[hgmask]
        mass_hotgas = mass_allgas[hgmask]

        vdisp_hotgas_global = util_galaxies.compute_vdisp_global(vel_allgas, mass_allgas,
                                vel_hotgas, mass_hotgas)
        vdisp_hotgas_los = util_galaxies.compute_vdisp_los(vel_allgas, mass_allgas,
                                pos_coldgas-starcen, vel_coldgas, rHM, mass_coldgas)

        # vdisp of gas near young stars (this is slow)
        indexes = []
        for x,pos in enumerate(pos_youngstars[:]):
            mask = np.linalg.norm(pos_allgas-pos,axis=1)<0.1

            for i in np.where(mask==True)[0]:
                indexes.append(i)
        indexes = np.unique(indexes)

        gaspos_sel = pos_allgas[indexes]
        gasmass_sel = (mass_allgas*sCDM.g['hydrogen'][gasmask])[indexes]
        gasvel_sel = vel_allgas[indexes]
        gastemp_sel = temp_allgas[indexes]

        vdisp_gasNearYS_los = util_galaxies.compute_vdisp_los(vel_allgas, mass_allgas,
                gaspos_sel-starcen, gasvel_sel, rHM, gasmass_sel)

        hgmask = gastemp_sel>10000
        gaspos_selhot = gaspos_sel[hgmask]
        gasmass_selhot = gasmass_sel[hgmask]
        gasvel_selhot = gasvel_sel[hgmask]

        vdisp_hotgasNearYS_global = util_galaxies.compute_vdisp_global(vel_allgas,
                mass_allgas, gasvel_selhot, gasmass_selhot)
        vdisp_hotgasNearYS_los = util_galaxies.compute_vdisp_los(vel_allgas, mass_allgas,
                gaspos_selhot-starcen, gasvel_selhot, rHM, gasmass_selhot)

        # alpha from power fit
        # profile range
        pmin = '0.01 kpc'
        #pmax = pynbody.analysis.halo.virial_radius(h1) / 4
        pmax = '50 kpc'

        # rotation curve
        # now we must center
        hCDM.d['pos']-=cen
        pynbody.analysis.angmom.faceon(hCDM)
        pdCDM = pynbody.analysis.profile.Profile(hCDM.d, rmin=pmin, rmax=pmax,
            type='lin', nbins=200)
        rbins = pdCDM['rbins']
        dmdensity = pdCDM['density']

        fitidx = 8
        alpha, c = opt.curve_fit(powerlaw, rbins[:fitidx], dmdensity[:fitidx],
                                    maxfev=1000,
                                    p0 = [-1,max(dmdensity)])[0]

        # create QA figure if desired
        if doQA:
            figfout_mask = QAdir+'/mask/r'+str(gal)+'/'+str(tstepnumber)+'.png'
            hw = 300
            #makeQAfig_mask(sCDM.s['pos'], sCDM.d['pos'], DMmask, cen, rHM,
            #            rfac*hmrDM, hw, figfout_mask)
            figfout_prof = QAdir+'/prof/r'+str(gal)+'/'+str(tstepnumber)+'.png'
            makeQAfig_profile(pdCDM, figfout_prof)


        # write to file
        fout.write(str(gal)+','+str(tstepnumber)+','+str(uage)+','+str(stepZ)+',')
        fout.write(str(mStar)+','+str(rHL)+','+str(rHL_c)+','+str(rHM)+',')
        fout.write(str(vdisp_allstars_global)+','+str(vdisp_allstars_los)+',')
        fout.write(str(vdisp_youngstar_global)+','+str(vdisp_youngstar_los)+',')
        fout.write(str(vdisp_allgas_global)+','+str(vdisp_allstars_los)+',')
        fout.write(str(vdisp_coldgas_global)+','+str(vdisp_coldgas_los)+',')
        fout.write(str(vdisp_hotgas_global)+','+str(vdisp_hotgas_los)+',')
        fout.write(str(vdisp_hotgasNearYS_global)+','+str(vdisp_gasNearYS_los)+',')
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