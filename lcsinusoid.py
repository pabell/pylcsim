"""
.. module:: lcsinusoid
   :synopsis: Generate coherent signals.

.. moduleauthor:: Riccardo Campana <campana@iasfbo.inaf.it>
"""
import numpy as np

def lcsinusoid(dt=1., nbins=65536, mean=0., freq=None, nha=1, amp=None, phi=None, verbose=False):
    """
    Generate coherent signals as a sequence of sinusoids (if len(freq) > 1) or
    of the fundamental frequency plus nha-1 harmonics (if len(freq) == 1),
    each with normalized pulsed fraction amp[i].

    Kwargs:
        dt:     time resolution of the lightcurve to be simulated (default: 1.0).

        nbins:  Number of bins of the simulated lightcurve (default: 65536).

        freq:   if float: frequency of the fundamental harmonic, if array: frequencies of sinusoids.

        nha:    number of harmonics (>1) 

        amp:    array with nha/nfreq elements; pulsed fraction for each frequency

        phi:    array with nha/nfreq elements; phases (in degrees!) for each frequency

    Returns:
        time:   time array

        rate:   array of count rates

    History:
        v1: Initial python implementation, from the IDL procedure lcharmonics.pro v0.0.3 
        by I. Donnarumma & R. Campana. Riccardo Campana, 2014. 
        
    """
    # Sanity checks
    if not amp:
        raise Exception('ERROR: Need at least one pulsed fraction!')
    
    try:
        nfreq = len(freq)
    except TypeError:
        nfreq = 1
            
    if nfreq == 1:
        # Sum of harmonics
        if amp:
            assert len(amp) == nha, 'ERROR: Need as many pulsed fraction as harmonics!'
            if not phi:
                phi = np.zeros(nha)
    else:
        # Sum of sinusoids
        if amp:
            assert len(amp) == nfreq, 'ERROR: Need as many pulsed fraction as frequencies!'
            if not phi:
                phi = np.zeros(nfreq)
      
    if verbose:     
        print("mean", mean)            

    # Array of time bins
    time = dt*np.arange(nbins)
    rate = np.zeros(nbins)
    
    phir = np.radians(phi)

    if nfreq == 1:
        # Cycle on harmonics
        for n in range(nha):
            rate += mean * amp[n] * np.sin(2.*np.pi*(n+1)*freq*time + phir[n]) 
    else:
        # Cycle on sinusoids        
        for i in range(nfreq):
            rate += mean * amp[i] * np.sin(2.*np.pi*freq[i]*time + phir[i]) 

    rate += mean
    
    return time, rate
    