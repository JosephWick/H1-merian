# util.py
# Joseph Wick
# H1-merian
#
# utility functions for H1 merian
#

import numpy as np
import pandas as pd

class util:

    # getfilepath()
    # returns filepath for zoom of halo h
    @staticmethod
    def getfilepath(h):
        '''
        gets file path for data of halo h in [CDM, SIDM] order
        '''

        cdm = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r'+str(h)+'.romulus25.3072g1HsbBH/r'+str(h)
        cdm += '.romulus25.3072g1HsbBH.004096'

        sidm = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r'+str(h)+'.romulus25cvdXsec.3072g1HsbBH/r'+str(h)
        sidm +='.romulus25cvdXsec.3072g1HsbBH.004096'

        adiabatic = '/data/REPOSITORY/e11Gals/romulus_dwarf_zooms/r'+str(h)+'.romulus25.3072g1HsbBH/adiabatic'

        return [cdm,sidm, adiabatic]
    # end getfilepath()

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
