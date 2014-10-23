import numpy as np
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
        v0.1:   Riccardo Campana, 2014.
        Initial python implementation.
    """
    
    # Set the seed for the random generator
    if seed:
        np.random.seed(seed)
    
    n = len(rate)
    # Arrays of total counts per bin
    # If l is less than 0 insert only background...
    l = (np.clip(rate, 0., np.max(rate)) + bkg)*dt 
    newrate = np.zeros(n)
    
    for i in xrange(n):
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
        v0.1:   Riccardo Campana, 2014.
        Initial python implementation.
    """
    
    # Sanity checks
    assert len(time) == len(rate), 'ERROR: Time and rate arrays with different dimensions'
    # Check that time array is evenly spaced (also taking into account rounding errors)
    assert np.std(np.diff(time)) <= 1e-12, 'ERROR: Time array is not evenly spaced'
    
    dt = time[1]-time[0]
    
    # Total number of photons
    counts = rate*dt
    n_phot = np.sum(counts)
    print "Total number of photons:", n_phot
    
    # Number of points
    n = len(counts)
    
    # Fourier transform
    f = np.fft.rfftfreq(n, d=dt)
    a = np.fft.rfft(counts)
    
    if norm == 'leahy':
        psd = 2./n_phot * np.abs(a)**2
    
    return f, psd


def rebin(x, y, nbins, mode='rate'):
    """
    Linearly rebins the (x, y) 1-D arrays using #nbins new bins, 
    taking for each new element the mean of the corresponding x elements,
    and the sum (if mode=='counts') or the mean (if mode=='rate') of the corresponding y elements.
      
    Args:
        x:   x array (same length as y)
        
        y:   y array (same length as x)
        
        nbins:   Number of new bins. Should be a factor of len(x) == len(y).
    
    Kwargs:
        mode: 'counts' or 'rate'; returns sum or mean of y-array elements
        
    Returns:
        xreb:   rebinned x array
        
        yreb:   rebinned y array
        
    History:
        v0.1:   Riccardo Campana, 2014.
        Initial python implementation.
    """
    # Sanity checks
    assert len(x) == len(y), 'ERROR: x and y must have the same length!'
    assert len(x) % nbins == 0,  'ERROR: nbins should be a factor of len(x) = len(y)'
    assert (mode == 'counts' or mode == 'rate'), 'ERROR: keyword mode should be counts or rate'
    
    
    factor = len(x)/nbins
    
    xreb = x.reshape((nbins, factor)).mean(axis=1)
    
    if mode == 'counts':
        yreb = y.reshape((nbins, factor)).sum(axis=1)
    elif mode == 'rate':
        yreb = y.reshape((nbins, factor)).mean(axis=1)
                
    return xreb, yreb


def logrebin(x, y, nbins, mode='rate'):
    """
    Logarithmically rebins the (x, y) 1-D arrays using #nbins new bins, 
    taking for each new element the logarithmic mean of the corresponding x elements,
    and the sum (if mode=='counts') or the mean (if mode=='rate') of the corresponding y elements.
      
    Args:
        x:   x array (same length as y)
        
        y:   y array (same length as x)
        
        nbins:   Number of new bins.
    
    Kwargs:
        mode: 'counts' or 'rate'; returns sum or mean of y-array elements
        
    Returns:
        xreb:   rebinned x array
        
        yreb:   rebinned y array
        
    History:
        v0.1:   Riccardo Campana, 2014.
        Initial python implementation.
    """
    # Sanity checks
    assert len(x) == len(y), 'ERROR: x and y must have the same length!'
    #assert len(x) % nbins == 0,  'ERROR: nbins should be a factor of len(x) = len(y)'
    assert (mode == 'counts' or mode == 'rate'), 'ERROR: keyword mode should be counts or rate'
    
    
    
    newbins = np.logspace(np.log10(np.min(x)+1e-12), np.log10(np.max(x)), nbins)
    
    xreb = newbins[:-1] + np.diff(newbins)/2
    
    digitized = np.digitize(x, newbins)
    
    if mode == 'counts':
        yreb = [y[digitized == i].sum() for i in xrange(1, len(newbins))]
    elif mode == 'rate':
        yreb = [y[digitized == i].mean() for i in xrange(1, len(newbins))]
                
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
        v0.1:   Riccardo Campana, 2014
        Initial python implementation.
    """
    # Sanity check
    assert len(time) == len(rate), 'ERROR: Time and rate should have the same length!'
    
    col1 = pyfits.Column(name='Time', unit='s', format='E', array=time)
    col2 = pyfits.Column(name='Rate', unit='counts/s', format='E', array=rate)
    
    cols = pyfits.ColDefs([col1, col2])
    tbhdu = pyfits.BinTableHDU.from_columns(cols)
    tbhdu.header['EXTNAME'] = ('LC', 'Name of this binary table extension') 
    prihdu = pyfits.PrimaryHDU()
    thdulist = pyfits.HDUList([prihdu, tbhdu])
    thdulist.writeto(outfilename, clobber=clobber)
    

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
        v0.1:   Riccardo Campana, 2014
        Initial python implementation.
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
    thdulist.writeto(outfilename, clobber=clobber)
