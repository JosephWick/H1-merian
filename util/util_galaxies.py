# util_galaxies.py
# Joseph Wick
# H1-merian
#
# utility functions for galaxy analysis
#

import numpy as np

class util_galaxies:

    # compute_CoM
    # computes center of mass given arrays of positions and masses
    @staticmethod
    def compute_CoM(pos, mass):
        '''
        compute_CoM()

        Computes center of mass given positions and masses of particles. Used
        when there is not a halo catalog for a given sim file, otherwise it is
        better to use Pynbody centering (momentum)

        Parameters
        ----------
        pos : array_like
            positions of particles

        mass : array_like
            masses of particles

        Returns
        -------
        cen : array_like
            Center of mass

        '''
        pos_np = np.array(pos)
        mass_np = np.array(mass)

        mtot = mass_np.sum()
        cen = np.sum(mass_np * pos_np.transpose(), axis=1) / mtot

        return cen


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

    # compute_massRadius()
    # computes radius within percentage of mass lie
    def compute_massRadius(positions, masses, outerR, acc, frac=0.5, maxiter=100000):
        innerLim = 0.0
        outerLim = outerR

        # center by mass
        mtot = masses.sum()
        cen = np.sum(masses * positions.transpose(),
             axis=1) / mtot
        cenpos = positions - cen

        # take radii 
        pRadii = np.linalg.norm(cenpos, ord=2, axis=1)
        mTot = sum(masses)

        r = outerR/2
        hm = sum(masses[pRadii < r])
        n=0
        while(hm < (frac-acc)*mTot or hm > (frac+acc)*mTot):
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
