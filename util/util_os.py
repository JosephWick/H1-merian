# util.py
# Joseph Wick
# H1-merian
#
# utility functions for navigating files system
#

import numpy as np
import pandas as pd
import glob

class util_os:

    # getfilepath()
    # returns filepath for zoom sim of gal, at a given timestep
    # tsidx refers to the index, not actual timestep number. The most
    # recent time step is indexed as 0
    @staticmethod
    def getfilepath(gal, tsidx):
        '''
        gets file path for data of halo h in [CDM, SIDM] order
        '''

        basedir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'

        # cdm
        galdir = basedir+'/r'+str(gal)+'.romulus25.3072g1HsbBH'

        timesteps = glob.glob(galdir+'/r*.romulus25.3072g1HsbBH.0*')
        timesteps.sort(reverse=True)
        timestep = timesteps[tsidx]
        tsnum = timestep[-6:]

        CDMsimFile = timestep+'/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tsnum

        # sidm
        galdir = basedir+'/r'+str(gal)+'.romulus25cvdXsec.3072g1HsbBH'

        timesteps = glob.glob(galdir+'/r*.romulus25cvdXsec.3072g1HsbBH')
        timesteps.sort(reverse=True)
        timestep = timesteps[tsidx]
        tsnum = timestep[-6:]

        SIDMsimFile = timestep+'/r'+str(gal)+'.romulus25cvdXsec.3072g1HsbBH.'+tsnum

        # return
        return CDMsimFile, SIDMsimFile
    # end getfilepath()

    # getfilepath_adiabatic()
    # returns z=0 adiabatic snapshot for given galaxy
    def getfilepath_adiabatic(gal):
        '''
        gets files path for adiabatic z=0 snapshot of galaxy gal
        '''

        basedir = '/data/REPOSITORY/e11gals/romulus_dwarf_zooms'

        adiaSimFile = basedir+'r'+str(gal)+'.romulus25.3072g1HsbBH/adiabatic'
        adiaSimFile+= '/r'+str(gal)+'/romulus25.3072g1HsbBH.004096'

    # getGalaxies()
    # returns halo ids for [CDM, SIDM] galaxies
    @staticmethod
    def getGalaxies():
        '''
        returns halo ids foor present [CDM,SIDM,adiabatic] galaxies
        '''

        df = pd.read_csv('/home/jw1624/H1-merian/util/currentGals.tsv', sep='\t')
        cdm = df[df['CDM']==1]['HaloID'].to_numpy().flatten()
        sidm = df[df['SIDM']==1]['HaloID'].to_numpy().flatten()
        adiabatic = df[df['adiabatic']==1]['HaloID'].to_numpy().flatten()

        return [cdm,sidm, adiabatic]
    # end getGalaxies()

    # compute_vdisp_std()
    # compute velocity dispersion using standard deviation
    @staticmethod
    def compute_vdisp_std(vel_all, mass_all, vel_t):
        '''
        compute_vdisp_std

        Computes unweighted velocity dispersion using standard deviation.

        Parameters
        -----
        vel_all : array like
            velocities of all stars. Used for center of mass velocity

        mass_all: array like
            masses of all stars. Used for center of mass velocity

        vel_t : array_like
            velocities of target stars we want to calculate v disp for

        Returns
        -----
        v_disp : float
            unweighted standard deviation of difference from systemic velocity
            for stars in vel_t and mass_t
        '''

        # check that both `_all` arrays have the same lenth and both `_t` arrays do too
        if len(vel_all)!=len(mass_all):
            print('Error: lengths of `vel_all` and `mass_all` must be identical')
            sys.exit()
        if len(vel_all)!=len(mass_all):
            print('Error: lengths of `vel_t` and `mass_t` must be identical')
            sys.exit()

        # this is COM velocity; systemic velocity
        v_CoM = np.sum(vel_all * mass_all[:, None], axis=0) / np.sum(mass_all)

        # difference from systemic velocity
        v_diffs = vel_t - v_CoM
        # take magnitude
        v_diffs_tot = np.linalg.norm(v_diffs, axis=1)

        return np.std(v_diffs_tot)
    # end compute_vdisp_std

    # compute_vdisp_wtd
    # compute mass weighted velocity dispersion
    @staticmethod
    def compute_vdisp_wtd(vel_all, mass_all, vel_t, mass_t):
        '''
        compute_vdisp_wtd

        Computes unweighted velocity dispersion using standard deviation.

        Parameters
        -----
        vel_all : array like
            velocities of all stars. Used for center of mass velocity

        mass_all: array like
            masses of all stars. Used for center of mass velocity

        vel_t : array_like
            velocities of target stars we want to calculate v disp for

        mass_t : array_like
            masses of target stars we want to calculate vdisp for

        Returns
        -----
        v_disp : float
            mass weighted standard deviation of difference from systemic velocity
            for stars in vel_t and mass_t
        '''

        # check that both `_all` arrays have the same length
        if len(vel_all)!=len(mass_all):
            print('Error: lengths of `vel_all` and `mass_all` must be identical')
            sys.exit()
        # check that both `_t` arrays have the same length
        if len(vel_all)!=len(mass_all):
            print('Error: lengths of `vel_t` and `mass_t` must be identical')
            sys.exit()

        # center of mass velocity / systemic velocity
        v_CoM = np.sum(vel_all * mass_all[:, None], axis=0) / np.sum(mass_all)

        # difference from systemic velocity
        v_diffs = vel_t - v_CoM
        # take magnitude
        v_diffs_tot = np.linalg.norm(v_diffs, axis=1)

        # mass weighted difference from mean vdiff
        vdisp = np.sum( ((v_diffs_tot-np.mean(v_diffs_tot))**2) * mass_t)
        # divide by total mass
        vdisp = vdisp/np.sum(mass_t)
        # sqrt
        vdisp = np.sqrt(vdisp)

        return float(vdisp)
