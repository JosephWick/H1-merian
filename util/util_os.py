# util_os.py
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
    # returns filepath for zoom cdm sim of gal, at a given timestep
    # Will throw an error if gal is not in getGalaxies()
    #
    # tsidx refers to the index, not actual timestep number.
    # The most recent time step is indexed as 0
    @staticmethod
    def getfilepath_cdm(gal, tsidx):
        '''
        getfilepath_cdm()

        Gets file path for data of halo h in [CDM, SIDM] order

        Parameters
        ----------
        gal : integer
            Galaxy ID from romulus

        tsidx : integer
            Index of desired timestep. Note that tsidx=0 refers to the latest
            timestep.

        Returns
        -------
        CDMsimFile : string
            filepath to the requested sim file ran under CDM
        '''

        # check against getGalaxies()
        if gal not in util_os.getGalaxies()[0]:
            raise Exception('Specified galaxy is not in list of available galaxies.')

        basedir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'

        # cdm
        galdir = basedir+'/r'+str(gal)+'.romulus25.3072g1HsbBH'

        timesteps = glob.glob(galdir+'/r*.romulus25.3072g1HsbBH.0*')
        timesteps.sort(reverse=True)
        timestep = timesteps[tsidx]
        tsnum = timestep[-6:]

        CDMsimFile = timestep+'/r'+str(gal)+'.romulus25.3072g1HsbBH.'+tsnum

        # return
        return CDMsimFile

    # getfilepath_sidm()
    # Returns filepath for sidm zoom sim of gal, at a given timestep.
    # Will throw an error if specified gal is not in getGalaxies().
    #
    # tsidx refers to the index, not actual timestep number.
    # The most recent time step is indexed as 0
    @staticmethod
    def getfilepath_sidm(gal, tsidx):
        '''
        getfilepath_cdm

        Gets file path for data of halo h in [CDM, SIDM] order

        Parameters
        ----------
        gal : integer
            Galaxy ID from romulus

        tsidx : integer
            Index of desired timestep. Note that tsidx=0 refers to the latest
            timestep.
        '''

        # check against getGalaxies()
        if gal not in util_os.getGalaxies()[1]:
            raise Exception('Specified galaxy is not in list of available galaxies.')

        basedir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'

        # cdm
        galdir = basedir+'/r'+str(gal)+'.romulus25cvdXsec.3072g1HsbBH'

        timesteps = glob.glob(galdir+'/r*.romulus25cvdXsec.3072g1HsbBH.0*')
        timesteps.sort(reverse=True)
        timestep = timesteps[tsidx]
        tsnum = timestep[-6:]

        CDMsimFile = timestep+'/r'+str(gal)+'.romulus25cvdXsec.3072g1HsbBH.'+tsnum

        # return
        return CDMsimFile

    # getfilepath_adiabatic()
    # returns z=0 adiabatic snapshot for given galaxy
    @staticmethod
    def getfilepath_adiabatic(gal):
        '''
        getfilepath_adiabatic()

        Returns filepath for adiabatic z=0 sim for specified galaxy.

        Parameters
        ----------
        gal : integer
            galaxy to retrieve adiabatic filepath for

        Returns
        -------
        adiaSimFile : string
            filepath to adiabatic z=0 simulation file
        '''
        # check against getGalaxies()
        if gal not in util_os.getGalaxies()[1]:
            raise Exception('Specified galaxy is not in list of available galaxies.')

        basedir = '/data/REPOSITORY/e11gals/romulus_dwarf_zooms'

        adiaSimFile = basedir+'/r'+str(gal)+'.romulus25.3072g1HsbBH/adiabatic'
        adiaSimFile+= '/r'+str(gal)+'.romulus25.3072g1HsbBH.004096'

        return adiaSimFile

    # getNumTimesteps()
    # returns number of timesteps for galaxy gal
    @staticmethod
    def getNumTimesteps(gal):
        '''
        getNumTimesteps()

        Returns number of timesteps for galaxy gal

        Parameters
        ----------
        gal : integer
            galaxy to retrieve timesteps for

        Returns
        -------
        num : integer
            number of timesteps present for galaxy gal
        '''

        basedir = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms'

        # cdm
        galdir = basedir+'/r'+str(gal)+'.romulus25.3072g1HsbBH'

        timesteps = glob.glob(galdir+'/r*.romulus25.3072g1HsbBH.0*')
        return len(timesteps)

    # getGalaxies()
    # returns halo ids for [CDM, SIDM, adiabatic] galaxies
    @staticmethod
    def getGalaxies():
        '''
        getGalaxies()

        Returns galaxy ids for present [CDM,SIDM,adiabatic] galaxies

        Parameters
        ----------
        None

        Returns
        -------
        [CDM, SIDM, Adiabatic] : array-like
            Arrays of galaxy IDs for which simulations of a given type are
            present. Information taken from `currentGals.tsv`
        '''

        df = pd.read_csv('/home/jw1624/H1-merian/util/currentGals.tsv', sep='\t')
        cdm = df[df['CDM']==1]['HaloID'].to_numpy().flatten()
        sidm = df[df['SIDM']==1]['HaloID'].to_numpy().flatten()
        adiabatic = df[df['adiabatic']==1]['HaloID'].to_numpy().flatten()

        return [cdm,sidm, adiabatic]
    # end getGalaxies()
