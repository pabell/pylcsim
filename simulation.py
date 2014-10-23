import numpy as np
from lcpsd import *
from lcsinusoid import *
from utils import *


class Simulation(object):
    """
    Main simulation class.
    
    History:
        v0.1:   Riccardo Campana, 2014. 
        Initial python implementation. 
    """
    def __init__(self, kind='psd'):
        # Sanity check
        assert kind in ['psd', 'coherent'], \
               "ERROR! Kind must be 'psd' or 'coherent'. See documentation for details."
        # Initialize with empty members
        self.kind = kind
        self.models = [] 
        self.time = []
        self.rate = []
        self.freq = []
        self.psd  = []
        # self.binsToKill = 0
        
    def addModel(self, modelName, modelParams):
        """
        Append simulation model to model dictionary
        """
        # Sanity check
        assert self.kind == 'psd', 'ERROR! You can add models only if simulation kind is PSD'
        self.models.append((modelName, modelParams))
        
    def info(self):
        """
        Prints simulation informations
        """
        modelnames = [x[0] for x in self.models]
        print "Simulation info"
        print "Kind: ", self.kind
        print "Model: ", "+".join(modelnames)
        
    def run(self, dt, nbins_old, mean, rms=None, freq=None, nha=None, amp=None, phi=None, verbose=False):
        """
        Run the simulation
        """
        # TODO: add sanity checks!
        
        # Rounding to the nearest power of 2 for the FFT
        # (greater than or equal to the original value)
        lg2  = np.log2(nbins_old)
        nbins = long(2**(np.ceil(lg2)))
        if verbose:
            print "Original nbins:\t\t",   nbins_old
            print "Rounded nbins:\t\t",    nbins
            print "Power of 2:\t\t",       int(np.ceil(lg2))
            print "Original exp. time (s):\t",  nbins_old*dt
            print "Rounded exp. time (s):\t",  nbins*dt
        if nbins > 2**24:
            nbins = 2**24
            if verbose:
                print "Maximum number of bins (2^24) exceeded. Reset to 2^24."
        # else:
        #     self.binsToKill = nbins-nbins_old  
        
        # Put control on input params
        if self.kind == 'psd':
            self.time, self.rate = lcpsd(dt=dt, nbins=nbins, mean=mean, rms=0.5, models=self.models)
        
        elif self.kind == 'coherent':
            self.time, self.rate = lcsinusoid(dt=dt, nbins=nbins, mean=mean, freq=freq, amp=amp, phi=phi, nha=nha)
        
        # self.time = self.time[:-self.binsToKill]
        # self.rate = self.rate[:-self.binsToKill]
        
    def poissonRandomize(self, dt, bkg):
        self.rate = poisson_randomization(self.rate, dt=dt, bkg=bkg)
    
    def getLightCurve(self):
        """
        Get lightcurve as time and rate arrays
        """
        return np.asarray(self.time), np.asarray(self.rate)
    
    def getPowerSpectrum(self):
        """
        Get power spectrum as frequency and power arrays
        """
        self.freq, self.psd = psd(self.time, self.rate)
        return np.asarray(self.freq), np.asarray(self.psd)

