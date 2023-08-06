# util_galaxies.py
# Joseph Wick
# H1-merian
#
# utility functions for galaxy analysis
#

import numpy as np

class util_galaxies:

    # compute_vdisp_std()
    # compute velocity dispersion using standard deviation
    @staticmethod
    def compute_vdisp_std(vel_all, mass_all, vel_t):
        '''
        compute_vdisp_std()

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
        vel_all : array_like
            velocities of all stars. Used for center of mass velocity

        mass_all: array_like
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
