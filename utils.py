from __future__ import division
from __future__ import print_function
from builtins import range
import numpy as np
import datetime
import astropy.io.fits as pyfits

def poisson_randomization(rate, dt=1., bkg=0., seed=None):
    """
    Poisson randomization of the rate array.
      
    Args:
        rate:   input array of count rates (in cts/s).
        
    Kwargs:    
        dt:     time resolution of the lightcurve to be simulated.
        
        bkg:    Mean count rate of the simulated lightcurve (default: 0.0).
        
        seed:   Seed for the random number generator.
    
    Returns:
        newrate: Array of Poisson randomized count rates of length n.
    
    History:	 
        v1:   Initial python implementation. Riccardo Campana, 2014.

    """
    
    # Set the seed for the random generator
    if seed:
        np.random.seed(seed)
    
    n = len(rate)
    # Arrays of total counts per bin
    # If l is less than 0 insert only background...
    l = (np.clip(rate, 0., np.max(rate)) + bkg)*dt 
    newrate = np.zeros(n)
    
    for i in range(n):
    	newrate[i] = np.random.poisson(l[i])/float(dt)

    return newrate



def psd(time, rate, norm='leahy'):
    """
    Returns power spectral density from a (real) time series
    with Leahy normalization.
    
    Args:
        time:   array of times (evenly binned).
        
        rate:   array of rate in counts/s.
        
    Kwargs:    
        norm:   Normalization (only Leahy for now).
        
    Returns:
        f:      array of frequencies.
        
        p:      power spectrum.
        
    History:
        v1:   Initial python implementation. Riccardo Campana, 2014.

    """
    
    # Sanity checks
    assert len(time) == len(rate), 'ERROR: Time and rate arrays with different dimensions'
    # Check that time array is evenly spaced (also taking into account rounding errors)
    assert np.std(np.diff(time)) <= 1e-12, 'ERROR: Time array is not evenly spaced'
    
    dt = time[1]-time[0]
    
    # Total number of photons
    counts = rate*dt
    n_phot = np.sum(counts)
    print("Total number of photons: {}".format(n_phot))
    
    # Number of points
    n = len(counts)
    
    # Fourier transform
    f = np.fft.rfftfreq(n, d=dt)
    a = np.fft.rfft(counts)
    
    if norm == 'leahy':
        psd = 2./n_phot * np.abs(a)**2
    
    return f[1:], psd[1:]


def rebin(x, y, factor, mode='rate', verbose=False):
    """
    Linearly rebins the (x, y) 1-D arrays by factor,
    i.e. with len(x)/factor new bins, 
    taking for each new element the mean of the corresponding x elements,
    and the sum (if mode=='counts') or the mean (if mode=='rate') of the corresponding y elements.
    If the new number of bins is not a factor of the old one, the array is cropped.
      
    Args:
        x:   x array (same length as y)
        
        y:   y array (same length as x)
        
        factor:   rebinning factor
    
    Kwargs:
        mode: 'counts' or 'rate'; returns sum or mean of y-array elements
        
    Returns:
        xreb:   rebinned x array
        
        yreb:   rebinned y array
        
    History:
        v2:    Switched to rebinning factor instead of number of new bins. Riccardo Campana, 2014.
        v1:    Initial python implementation. Riccardo Campana, 2014.

    """
    # Sanity checks
    assert len(x) == len(y), 'ERROR: x and y must have the same length!'
    assert (mode == 'counts' or mode == 'rate'), 'ERROR: keyword mode should be counts or rate'
    
    oldbins = len(x)
    cropping_limit = oldbins - (oldbins % int(factor))
    
    if verbose:
        newbins = int(cropping_limit/float(factor))
        print("Rebinning a {} long array by a factor {} with {} new bins".format(oldbins, factor, newbins))
    
    cropped_x = x[:cropping_limit]
    cropped_y = y[:cropping_limit]
    
    xreb = np.mean( np.concatenate([[cropped_x[i::factor] for i in range(factor)] ]), axis=0)
    
    if mode == 'counts':
        yreb = np.sum( np.concatenate([[cropped_y[i::factor] for i in range(factor)] ]), axis=0)
    elif mode == 'rate':
        yreb = np.mean( np.concatenate([[cropped_y[i::factor] for i in range(factor)] ]), axis=0)  
         
    return xreb, yreb


def logrebin(x, y, factor, mode='rate', verbose=False):
    """
    Logarithmically rebins the (x, y) 1-D arrays by a constant logarithmic bin log(1+1/factor),
    i.e. each new bin has a width (1+1/factor) greater than the preceding;
    taking for each new element the logarithmic mean of the corresponding x elements,
    and the sum (if mode=='counts') or the mean (if mode=='rate') of the corresponding y elements.
      
    Args:
        x:   x array (same length as y)
        
        y:   y array (same length as x)
        
        factor:   logarithmic rebinning factor
    
    Kwargs:
        mode: 'counts' or 'rate'; returns sum or mean of y-array elements
        
    Returns:
        xreb:   rebinned x array
        
        yreb:   rebinned y array
        
    History:
        v2:    Switched to logarithmic rebinning factor. Riccardo Campana, 2014.     
        v1:    Initial python implementation. Riccardo Campana, 2014.

    """
    # Sanity checks
    assert len(x) == len(y), 'ERROR: x and y must have the same length!'
    assert (mode == 'counts' or mode == 'rate'), 'ERROR: keyword mode should be counts or rate'
    
    xreb = []
    newx = x[0]
    while newx <= x[-1]:
        xreb.append(newx)
        newx = newx * (1 + 1./factor)
    
    
    digitized = np.digitize(x, xreb)
    newbins = len(xreb)
    
    if mode == 'counts':
        yreb = [y[digitized == i].sum() for i in range(newbins)]
    elif mode == 'rate':
        import warnings; warnings.filterwarnings('ignore') # To avoid RuntimeWarning if empty slice
        yreb = [y[digitized == i].mean() for i in range(newbins)]
                
    return np.asarray(xreb), np.asarray(yreb)
    
    
def saveFITSLC(outfilename, time, rate, clobber=True):
    """
    Produce an output FITS file containing a lightcurve.

    Args:
        outfilename: Name of the output FITS file

        time: Array of times

        rate: Array of count rates

    Kwargs:
        clobber: if True, overwrites existing files with same name

    Returns:
        none

    History:
        v2:   OGIP-compliance (OGIP 93-003). Riccardo Campana, 2014.
        
        v1:   Initial python implementation. Riccardo Campana, 2014.
        
    """
    # Sanity check
    assert len(time) == len(rate), 'ERROR: Time and rate should have the same length!'
    
    timedelta = time[1]-time[0]
    
    # Dummy start date for observation
    start_date = datetime.datetime(2015, 1, 1)
    stop_date  = start_date + datetime.timedelta(seconds=time[-1])
    
    
    col1 = pyfits.Column(name='TIME', unit='s', format='D', array=time)
    col2 = pyfits.Column(name='RATE', unit='count/s', format='D', array=rate)
    
    cols = pyfits.ColDefs([col1, col2])
    tbhdu = pyfits.BinTableHDU.from_columns(cols)
    
    tbhdu.header['EXTNAME']    = ('RATE', 'Name of this binary table extension') 
    tbhdu.header['TELESCOP']   = ('PYLCSIM', 'Mission or telescope name') 
    tbhdu.header['INSTRUME']   = ('PYLCSIM', 'Instrument name') 
    tbhdu.header['ORIGIN']     = ('PYLCSIM', 'Who produced this file') 
    tbhdu.header['TIMVERSN']   = ('OGIP/93-003', 'OGIP memo describing the convention used') 
    tbhdu.header['AUTHOR']     = ('PYLCSIM', 'Name of the program that produced this file') 
    
    tbhdu.header['RA']         = (0, 'Source right ascension in degrees') 
    tbhdu.header['DEC']        = (0, 'Source declination in degrees') 
    
    tbhdu.header['DATE-OBS']   = (start_date.strftime("%d/%m/%y"), 'Date of observation start') 
    tbhdu.header['TIME-OBS']   = (start_date.strftime("%H:%M:%S.%f"), 'Time of observation start') 
    tbhdu.header['DATE-END']   = (stop_date.strftime("%d/%m/%y"), 'Date of observation end') 
    tbhdu.header['TIME-END']   = (stop_date.strftime("%H:%M:%S.%f"), 'Time of observation end') 
    
    tbhdu.header['TSTART']     = (0., 'Start time')
    tbhdu.header['TSTOP']      = (time[-1], 'Stop time')
    tbhdu.header['TIMEZERO']   = (0., 'Zero time used to calculate the n-th event')
    tbhdu.header['TIMESYS']    = ('2015.0', 'System used to define time') 
    tbhdu.header['TIMEUNIT']   = ('s', 'Unit for TSTART, TSTOP, TIMEZERO') 
    tbhdu.header['CLOCKCOR']   = ('YES', 'If time corrected to UT') 
    tbhdu.header['MJDREF']     = (57023.0, 'MJD for reference time') 
    
    tbhdu.header['TIMEDEL']    = (timedelta, 'Source declination in degrees') 
    
    prihdu = pyfits.PrimaryHDU()
    thdulist = pyfits.HDUList([prihdu, tbhdu])
    thdulist.writeto(outfilename, overwrite=clobber)
    

def saveFITSPSD(outfilename, freq, psd, clobber=True):
    """
    Produce an output FITS file containing a power spectrum.

    Args:
        outfilename: Name of the output FITS file

        freq: Array of frequencies

        psd: Array of power spectrum

    Kwargs:
        clobber: if True, overwrites existing files with same name

    Returns:
        none

    History:
        v1:   Initial python implementation. Riccardo Campana, 2014.
        
    """
    # Sanity check
    assert len(freq) == len(psd), 'ERROR: Frequencies and PSD should have the same length!'

    col1 = pyfits.Column(name='Frequency', unit='Hz', format='E', array=freq)
    col2 = pyfits.Column(name='PSD', format='E', array=psd)

    cols = pyfits.ColDefs([col1, col2])
    tbhdu = pyfits.BinTableHDU.from_columns(cols)
    tbhdu.header['EXTNAME'] = ('PSD', 'Name of this binary table extension') 
    prihdu = pyfits.PrimaryHDU()
    thdulist = pyfits.HDUList([prihdu, tbhdu])
    thdulist.writeto(outfilename, overwrite=clobber)
