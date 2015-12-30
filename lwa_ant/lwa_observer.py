from lwa_antenna import *
import ephem
from astropy import units as u

from pygsm import GSMObserver

ZEN_AZ_RAD  = 0           # Zentih azimuth in radians
ZEN_ALT_RAD = np.pi / 2   # Zenith altitude in radians

class LwaObserver(GSMObserver):
    """ Observer class for LWA. 
    
    Adds functionality from pyephem and pygsm, so that 
    * Observed skies can be generated from Global Sky Model
    * Beam patterns can be shifted into galactic frame
    
    """
    def __init__(self, site='OVRO'):

        super(LwaObserver, self).__init__()

        # Owen's Valley
        if site in ('LWA1', 'LWANM', 'LWA-NM'):
            self.site = 'LWA-NM'
            (latitude, longitude, elevation) = ('34.070', '-107.628', 2133.6)
        elif site in ('LWA-OVRO', 'OVRO', 'LEDA-OVRO', 'LEDAOVRO'):
            self.site = 'LWA-OVRO'
            (latitude, longitude, elevation) = ('37.2', '-118.2', 1222)
        elif site in ('NORTH-POLE', 'NORTH POLE', 'NP'):
            self.site = 'NORTH-POLE'
            (latitude, longitude, elevation) = ('90', '0', 0)
        elif site in ('SOUTH-POLE', 'SOUTH POLE', 'SP'):
            self.site = 'SOUTH-POLE'
            (latitude, longitude, elevation) = ('-90', '0', 0)
        elif site in ('EQUATOR', 'EQ'):
            self.site = 'EQUATOR'
            (latitude, longitude, elevation) = ('0', '0', 0)

        self.lon = longitude
        self.lat = latitude
        self.elev = elevation
        
        self.beam = LwaBeamPattern()

    def ra_dec_zenith(self, date=None):
        """ Compure the RA and DEC of zenith

        Notes
        -----
        This does not just return LST and latitude, it uses pyephem
        to compute the RA and DEC with precession and nutation
        corrections.

        Parameters
        ----------
        date: datetime.datetime
            date-time of observation (optional)

        Returns
        -------
        (ra, dec) of zenith, as astropy.unit.Quantities
            ra is in hourangle, and dec is in degrees.
        """
        if date:
            self.date = date

        zen_ra_rad, zen_dec_rad = self.radec_of(ZEN_AZ_RAD, ZEN_ALT_RAD)
        zen_ra = zen_ra_rad * u.rad
        zen_dec = zen_dec_rad * u.rad

        return zen_ra.to('hourangle'), zen_dec.to('deg')

    def altaz2radec(self, alt, az):
        """ Convert alt, az to right ascension and declination

        Parameters
        ----------
        alt: int or np.array
            Altitude values. Can be 1D array
        az:  int or np.array
            Azimuth values. Can be 1D array

        Returns
        -------
        ra_dec: np.array
            Returns array of RA values and array of DEC values.
            Shape is (n_pointings, 2).

        Notes
        -----
        This is just a vectorized version of radec_of,
        it doesn't offer much of a speed boost but is convenient.
        """
        el = alt    #elevation == altitude
        radec_vfunc = np.vectorize(self.radec_of)
        ra, dec = radec_vfunc(az, el)
        return np.column_stack((ra, dec))
    
    def view_beam(self, show=True):
        self.beam.view(show=show)
    
    def view_gsm(self, show=True):
        self.view(show=show)

if __name__ == "__main__":
    leda = LwaObserver()
    
    leda.generate(50)
    leda.beam.generate(50)
    
    leda.view_gsm(show=True)
    leda.view_beam(show=True)
    

