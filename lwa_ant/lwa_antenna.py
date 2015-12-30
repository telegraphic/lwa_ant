import numpy as np
import pylab as plt
import h5py

import healpy as hp

from scipy.interpolate import interp1d, RectBivariateSpline

from grid_utils import ungrid, meshgrid_pairs, grid_xyz
from spherical_coords import thetaphi_to_azel

from pkg_resources import resource_filename

H5_DATA = resource_filename("lwa_ant", "antenna_data.h5")

class LwaBeamPattern(object):
    """ Beam pattern class for the LWA antenna """
    def __init__(self):
        self.h5 = h5py.File(H5_DATA)
        self.data   = self.h5["data"][:]
        self.theta  = self.h5["theta"][:]
        self.phi    = self.h5["phi"][:]
        self.freqs  = self.h5["frequency"][:]
        
        self.freq_interpolator = interp1d(self.freqs, self.data, axis=0, kind='cubic')
        
        self.generated_beam_data  = None
        self.generated_beam_freq  = None
        self.generated_beam_theta = None
        self.generated_beam_phi   = None
    
    def _check_generated(self):
        """ Function to check if beam is already generated"""
        if self.generated_beam_data is None:
            raise RuntimeError("Beam pattern not generated yet. Run generate() first.")
    
    def generate(self, freq, n_theta=None, n_phi=None, linear=False):
        """ Generate a beam for the LWA antenna at a given frequency 
        
        Parameters
        ----------
        freq: float
            Frequency in MHz for beam pattern
        n_theta: int
            Number of points to generate on theta axis (0 to 360 deg). Optional.
        n_phi: int
            Number of points to generate on phi axis (-90 to 90 deg). Optional.
        linear: bool
            Output beam pattern in linear units. Defaults to False (dB).
        
        Notes
        -----
        This uses EM simulations done in 4NEC2, every 5MHz, 1 degree resolution,
        from 20 MHz to 90 MHz. File used:
        `LWA_DUAL_0-90DEG-DRIVE_300x300x1cmGRID_Symbols.nec`
        
        Interpolation (cubic) is done over frequency, to generate a 1 degree 
        resolution pattern at the new frequency. If more points are required 
        (i.e. higher res) then interpolation is done once again to achieve this,
        this time using RectBivariateSpline().
        """
        beam = self.freq_interpolator(freq)
        
        if n_theta is not None and n_phi is not None:
            theta_new = np.linspace(-90, 90, n_theta)
            phi_new   = np.linspace(0, 360, n_phi)
            
            # Do fit in linear space
            beam_lin = 10**(beam / 10.0)
            b_interp = RectBivariateSpline(self.phi, self.theta, beam_lin)
            beam     = 10 * np.log10(b_interp(phi_new, theta_new))
        else:
            theta_new, phi_new = self.theta, self.phi
            
        self.generated_beam_data  = beam
        self.generated_beam_theta = theta_new
        self.generated_beam_phi   = phi_new
        self.generated_beam_freq  = freq
        
        return beam
    
    def to_azel(self):
        """ Convert theta-phi beam to azimuth-elevation 
        
        Returns
        -------
        beam_azel: np.array
            Beam pattern on azimuth-elevation grid.
        
        """
        self._check_generated()
        theta, phi = np.deg2rad(self.generated_beam_theta), np.deg2rad(self.generated_beam_phi)
        tp = meshgrid_pairs(theta, phi)
        tpu = ungrid(theta, phi, self.generated_beam_data)

        az, el = thetaphi_to_azel(tpu[:, 0], tpu[:, 1])
        beam_lin = 10**(tpu[:, 2] / 10.0)
        azelb = np.column_stack((az, el, beam_lin))
        
        beam_azel = grid_xyz(azelb, 361, 361)
        return beam_azel
    
    def view(self, show=True, view='healpix'):
        """ View generated beam, defaulting to view healpix orthographic """     
        
        if view in ('healpix', 'orthographic', 'orthview'):
            self.view_healpix(show=show)
        elif view in ('thetaphi', 'theta_phi', 'antcoords'):
            self.view_thetaphi(show=show)
        else:
            self.view_azel(show=show)
    
    def view_thetaphi(self, show=True):
        """ View generated beam, in theta-phi coordinates """
        
        self._check_generated()
        
        plt.figure(figsize=(8,4))
        plt.subplot(121)
        tmin, tmax = self.generated_beam_theta[0], self.generated_beam_theta[-1]
        pmin, pmax = self.generated_beam_phi[0], self.generated_beam_phi[-1]
        plt.imshow(10**(self.generated_beam_data / 10), extent=(tmin, tmax, pmin, pmax), aspect='auto')
        plt.xlabel("Theta [deg]")
        plt.ylabel("Phi [deg]")
        #plt.colorbar(orientation='horizontal')
        
        plt.subplot(122)
        beam_slice = self.generated_beam_data[self.generated_beam_data.shape[0]/2]
        print self.generated_beam_phi.shape, beam_slice.shape
        plt.plot(self.generated_beam_theta, beam_slice, c='#333333')
        plt.xlabel("Theta [deg]")
        plt.ylabel("Normalized gain [dB]")
        plt.xlim(-91, 91)
        plt.ylim(-30, 3)
        plt.minorticks_on()
        
        plt.tight_layout()
        if show:
            plt.show()

    def view_azel(self, show=True):
        """ View generated beam, in azimuth-elevation coordinates """
        
        self._check_generated()
        beam = self.to_azel()
        
        plt.figure(figsize=(8,4))
        tmin, tmax = -180, 180
        pmin, pmax = -90, 90
        plt.imshow(beam, extent=(tmin, tmax, pmin, pmax), aspect='auto')
        plt.xlabel("Azimuth [deg]")
        plt.ylabel("Elevation [deg]")
        plt.colorbar()
        #plt.colorbar(orientation='horizontal')
        if show:
            plt.show()
    
    def to_healpix(self, n_side=32):
        """ Convert beam pattern to a healpix 
        
        Returns
        -------
        hmap: np.array
            Numpy array representing healpix map. Array is in observer frame,
            i.e. zenith aligned with equatorial coordinate system.
        """
        beam_azel = self.to_azel()
        
        el = np.linspace(0, np.pi, beam_azel.shape[1])
        az = np.linspace(-np.pi, np.pi, beam_azel.shape[0])[:, None]
        
        # Generate beam pattern with n_side = 32 (for gridding reasons)
        hmap = hp.ma(np.zeros(hp.nside2npix(32)))
        pix = hp.ang2pix(32, el, az)
        hmap[pix] = beam_azel
        
        # Upsample if required
        if n_side != 32:
            hmap = hp.ud_grade(hmap, nside_out=n_side)

        # Apply mask
        n_pix = hp.nside2npix(n_side)
        theta, phi = hp.pix2ang(n_side, np.arange(n_pix))
        mask1 = phi + np.pi/2 >  2 * np.pi 
        mask2 = phi < np.pi / 2
        hmap.mask = np.invert(np.logical_or(mask1, mask2))
                
        return hmap
    
    def view_healpix(self, n_side=32, show=True):
        """ View the beam pattern as a healpix map
        
        Shows orthographic view of beam pattern, using healpy plotting routines.
        """
        hmap = self.to_healpix(n_side)
        
        hp.orthview(hmap)
        if show:
            plt.show()
        
            
if __name__ == "__main__":
    lwa = LwaBeamPattern()
    lwa.generate(60)
    lwa.view_healpix(n_side=128)
    
        