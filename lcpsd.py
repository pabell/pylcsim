import numpy as np
from .psd_models import *

def lcpsd(dt=1., nbins=65536, mean=0., rms=1., seed=None, models=None, phase_shift=None):
    """
    Simulate a light-curve with a general power spectrum shape.
    For the underlying algorithm see: 
    J. Timmer & M. Koenig, "On generating power law noise", 
    A&A, 300, 707-710 (1995).
    
    Kwargs:
        dt:     time resolution of the lightcurve to be simulated
        
        nbins:  Number of bins of the simulated lightcurve (default:65536).
        Should be power of two for optimum performance (FFT...)
                
        mean:   Mean count rate of the simulated lightcurve (default: 0.).
        
        rms:    Total fractional RMS of the simulated lightcurve.
        
        seed:   Seed for the random number generator.
        
        models: List of tuples, each containing:
        a. Name of a function returning the desired PSD shape,
        b. Parameters of the model (argument to model).
        Total model is the sum of these tuples.
                
        phase_shift: Constant phase shift (in degrees) to the FFT.
        
    Returns:
        time:   time array.
        
        rate:   array of count rates.
    
    History:	 
        v0.1:   Riccardo Campana, 2014. 
        Initial python implementation. 
        Based on AITLIB IDL procedure timmerlc.pro
    """


    # Raise exception if a model is not explicitly given
    if not models:
        raise Exception('ERROR: Need at least one model to continue!') 
    
    # Set the seed for the random generator
    if seed:
        np.random.seed(seed)
        
    
    # Frequencies at which the PSD is to be computed 
    # (only positive frequencies, since the signal is real)
    # simfreq = ( np.arange(0, nbins/2)+1. ) / float(dt*nt)
    simfreq = np.fft.rfftfreq(nbins, d=dt)[1:]
    
    simpsd = np.zeros(len(simfreq))
    # Compute PSD from models
    for single_model in models:
        model  = single_model[0] 
        params = single_model[1]
        simpsd += eval(model + "(simfreq, params)")

    # print "len(simfreq)", len(simfreq)
    # print "len(simpsd)", len(simpsd)
    
    
    fac = np.sqrt(simpsd/2.)
 
    print nbins
    if phase_shift:
        ph_sh_rad = np.radians(phase_shift)
        pos_real_i = np.random.normal(size=nt/2)*fac
        pos_imag_i = np.random.normal(size=nt/2)*fac
        pos_real =  pos_real_i * np.cos(ph_sh_rad) - pos_imag_i * np.sin(ph_sh_rad)
        pos_imag =  pos_real_i * np.sin(ph_sh_rad) + pos_imag_i * np.cos(ph_sh_rad)
    else:
        pos_real = np.random.normal(size=nbins/2)*fac
        pos_imag = np.random.normal(size=nbins/2)*fac

    # print " n_elements(pos_real)", len(pos_real)

    pos_freq_transform = pos_real + 1j * pos_imag

    # Simulate light curve from its Fourier transform
    arg  = np.concatenate(([mean], pos_freq_transform))

    # print "len(arg)", len(arg)
    
    # Inverse Fourier transform
    rate = np.fft.irfft(arg)
    

    print len(rate)
    # Array of time bins
    time = dt*np.arange(nbins)

    # Average and standard deviation
    avg = np.mean(rate)
    std = np.std(rate)
    # print "Generated curve avg ", avg, " std ", std

    # Rescaling of the light curve
    # 1. to zero mean and normalized to one standard deviation
    # 2. to "mean" mean and std. dev. = mean*rms
    # (where rms is the fractional rms)
    rate = (rate-avg)/std * mean * rms + mean

    return time, rate

