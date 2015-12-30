"""
# spherical_coords.py

Convert between (theta, phi) and (azimuth, elevation) coordinate systems.

Author: Danny Price
License: MIT

## Overview

Both (theta, phi) and (azimuth, elevation) are spherical coordinate systems. 
The first is used commonly in antenna pattern measurements, while the second 
is used extensively in ground-based astronomy. Conversion between the two 
essentially requires moving the location of the pole from one axis to another.

In (theta, phi) coordinates, phi is the angle from the y-axis toward the z-axis, 
as measured from the yz plane. Phi runs [0, 360) degrees. Theta runs [0, 180).

In (azimuth, elevation) coordinates, azimuth is the angle from the x-axis
toward  y-ax, as measured from the xy plane. Azimuth runs [-180, 180) degrees. 
Elevation runs [-90 to 90) degrees.

A more complete explanation of the coordinate systems can be found 
in [1] and [2]. Coordinate transforms, as coded up, are based on the
formulas presented in [1].

## SUPER IMPORTANT NOTE TO READ BEFORE USE

When writing some unit tests, I noticed that a round-trip conversion does
not return the identical angle pair that you started with. My first thought
was that this is due to the interval bounds on inverse functions (arcsin etc),
but some things seem a little odd still. So I'm releasing this as a gist, so
people can decide whether it's a useful start, and hopefully someone will 
point to the problem and fix it one day.

## References
[1] phitheta2azel MATLAB documentation, Mathworks. 
    url: http://www.mathworks.com/help/phased/ref/phitheta2azel.html
[2] COORDINATE SYSTEM PLOTTING FOR ANTENNA MEASUREMENTS, Gregory F. Masters 
    and Stuart F. Gregson. Tech. Report, Nearfield Systems, 2007. 
    url: http://ww2.nearfield.com/amta/AMTA07-0092-GFM_SFG.pdf
"""

import numpy as np

def azel_to_thetaphi(az, el):
    """ Az-El to Theta-Phi conversion.
  
    Args:
        az (float or np.array): Azimuth angle, in radians
        el (float or np.array): Elevation angle, in radians
  
    Returns:
      (theta, phi): Tuple of corresponding (theta, phi) angles, in radians
    """
    
    cos_theta = np.cos(el) * np.cos(az)
    tan_phi   = np.tan(el) / np.sin(az)
    theta     = np.arccos(cos_theta)
    phi       = np.arctan2(np.tan(el), np.sin(az))
    phi = (phi + 2 * np.pi) % (2 * np.pi)
    
    return theta, phi
  

def thetaphi_to_azel(theta, phi):
    """ Az-El to Theta-Phi conversion.
  
    Args:
        theta (float or np.array): Theta angle, in radians
        phi (float or np.array): Phi angle, in radians
  
    Returns:
      (az, el): Tuple of corresponding (azimuth, elevation) angles, in radians
    """
    sin_el = np.sin(phi) * np.sin(theta)
    tan_az = np.cos(phi) * np.tan(theta)
    el = np.arcsin(sin_el)
    az = np.arctan(tan_az)
      
    return az, el

def thetaphi_to_uv(theta, phi):
    """ Theta-phi to U-V direction cosines
  
    Args:
        theta (float or np.array): Theta angle, in radians
        phi (float or np.array): Phi angle, in radians
  
    Returns:
      (u, v): Tuple of corresponding (u, v) direction cosines
    """
    u = np.sin(theta) * np.cos(phi)
    v = np.sin(theta) * np.sin(phi)
    w = np.cos(theta)
    
    return u, v

def azel_to_uv(az, el):
    """ Azimuth-elevation to U-V direction cosines
  
    Args:
        az (float or np.array): azimuth angle, in radians
        el (float or np.array): elevation angle, in radians
  
    Returns:
      (u, v): Tuple of corresponding (u, v) direction cosines
    """    
    u = np.cos(el) * np.sin(az)
    v = np.sin(el)
    w = np.cos(az) * np.cos(el)
    
    return u, v

def uv_to_thetaphi(u, v):
    """ U-V direction cosines to theta-phi coordinates
  
    Args:
        u (float or np.array): U direction cosine
        v (float or np.array): V direction cosine
  
    Returns:
      (theta, phi): Tuple of corresponding (theta, phi) angles, in radians
    """ 
    theta = np.arcsin(np.sqrt(u**2 + v**2))
    phi = np.arctan2(u, v)
    
    phi = (phi + 2 * np.pi) % (2 * np.pi)
    
    return theta, phi
    
def uv_to_azel(u, v):
    """ U-V direction cosines to azimuth-elevation coordinates
  
    Args:
        u (float or np.array): U direction cosine
        v (float or np.array): V direction cosine
  
    Returns:
      (az, el): Tuple of corresponding (azimuth, elevation) angles, in radians
    """ 
    az = np.arctan2(u, np.sqrt(1 - u**2 - v**2))
    el = np.arcsin(v)
    
    return az, el

def plot_compare(a, b, ap, bp, title='', show=False, nf=True):
    """ Helper function for debugging, plotting before & after """
    if nf:
        plt.figure()
    plt.subplot(211)
    plt.plot(a, ap, '.', label=title)
    plt.legend()
    plt.subplot(212)
    plt.plot(b, bp, '.', label=title) 
    if show:
        plt.show()

if __name__ == "__main__":
    import pylab as plt
    
    # Generate vectors covering intervals of theta & phi
    theta = np.linspace(0, np.pi,   101)
    phi   = np.linspace(0, 2*np.pi, 101)
    
    # Convert to other coordinates and back
    u, v = thetaphi_to_uv(theta, phi)
    theta_p, phi_p = uv_to_thetaphi(u, v)
    plot_compare(theta, phi, theta_p, phi_p, 
                 title='T-P -> U-V -> T-P')

    az, el = thetaphi_to_azel(theta, phi)
    theta_p, phi_p = azel_to_thetaphi(az, el)
    plot_compare(theta, phi, theta_p, phi_p, 
                 title='T-P -> A-E -> T-P', nf=False)
    
    u, v = thetaphi_to_uv(theta, phi)
    az, el = uv_to_azel(u, v)
    theta_p, phi_p = azel_to_thetaphi(az, el)
    plot_compare(theta, phi, theta_p, phi_p, 
                 title='T-P -> U-V -> A-E -> T-P', nf=False)

    u, v = thetaphi_to_uv(theta, phi)
    az, el = uv_to_azel(u, v)
    u, v = azel_to_uv(az, el)
    theta_p, phi_p = uv_to_thetaphi(u, v)
    plot_compare(theta, phi, theta_p, phi_p, 
                 title='T-P -> U-V -> A-E -> U-V -> T-P', nf=False)    
    
    # Convert to Az-el, two methods
    u, v = thetaphi_to_uv(theta, phi)
    az, el = uv_to_azel(u, v)
    plot_compare(theta, phi, az, el, 
                 title='T-P -> U-V -> A-E')    

    az, el = thetaphi_to_azel(theta, phi)
    plot_compare(theta, phi, az, el, 
                 title='T-P -> A-E', nf=False)    
    

    # Generate vectors covering intervals of azimuth & elevation
    az = np.linspace(-np.pi, np.pi,   101)
    el = np.linspace(-np.pi/2, np.pi/2, 101)
    
    # Convert to UV and back, compare
    u, v = azel_to_uv(az, el)
    az_p, el_p = uv_to_azel(u, v)
    plot_compare(az, el, az_p, el_p, 
                 title='A-E -> U-V -> A-E')
                 
    theta, phi = azel_to_thetaphi(az, el)
    az_p, el_p = thetaphi_to_azel(theta, phi)
    plot_compare(az, el, az_p, el_p, 
                 title='A-E -> T-P -> A-E', nf=False)

    u, v = azel_to_uv(az, el)
    theta, phi = uv_to_thetaphi(u, v)
    az_p, el_p = thetaphi_to_azel(theta, phi)
    plot_compare(az, el, az_p, el_p, 
                 title='A-E -> U-V -> T-P -> A-E', nf=False)

    u, v = azel_to_uv(az, el)
    theta, phi = uv_to_thetaphi(u, v)
    u, v = thetaphi_to_uv(theta, phi)
    az_p, el_p = uv_to_azel(u, v)
    plot_compare(az, el, az_p, el_p, 
                 title='A-E -> U-V -> T-P -> U-V -> A-E', nf=False)

    # Convert to Theta-phi, two methods
    u, v = thetaphi_to_uv(az, el)
    theta, phi = uv_to_thetaphi(u, v)
    plot_compare(az, el, theta, phi,
                 title='A-E -> U-V -> T-P')    

    theta, phi = azel_to_thetaphi(az, el)
    plot_compare(az, el, theta, phi,
                 title='A-E -> T-P', nf=False)   

    plt.show()


  