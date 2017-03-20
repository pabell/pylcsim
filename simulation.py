from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import object
import numpy as np
from astropy.modeling import models
from .lcpsd import *
from .lcsinusoid import *
from .utils import *


class Simulation(object):
    """
    Main simulation class.

    History:
        v1.3:   Added getPSDModel method. Riccardo Campana, 2014.
        
        v1.2:   Added reset method. Riccardo Campana, 2014.
        
        v1.1:   Bugfix. Riccardo Campana, 2014.
        
        v1:     Initial python implementation. Riccardo Campana, 2014. 
        
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
        self.binsToKill = 0
        
        
    def reset(self):
        """
        Reset the simulation, emptying models array and other members
        """
        self.models = [] 
        self.time = []
        self.rate = []
        self.freq = []
        self.psd  = []
        self.binsToKill = 0
        
        
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
        modelnames = []
        try:
            basestring
        except NameError:
            basestring = str
        for x in self.models:
            if isinstance(x[0], basestring):
                modelnames.append(x[0])
            else:
                modelnames.append(x[0].__name__)     
        print("Simulation info")
        print("Kind: ".format(self.kind))
        print("Model: {:s}".format("+".join(modelnames)))
        
        
    def run(self, dt, nbins_old, mean, rms=None, freq=None, nha=None, amp=None, phi=None, verbose=False):
        """
        Run the simulation
        """
        # TODO: add sanity checks!
        
        # Rounding to the nearest power of 2 for the FFT
        # (greater than or equal to the original value)
        lg2  = np.log2(nbins_old)
        nbins = int(2**(np.ceil(lg2)))
        if verbose:                                          
            print(("Original nbins:\t\t",   nbins_old         ))
            print(("Rounded nbins:\t\t",    nbins             ))
            print(("Power of 2:\t\t",       int(np.ceil(lg2)) ))
            print(("Original exp. time (s):\t",  nbins_old*dt ))
            print(("Rounded exp. time (s):\t",  nbins*dt      ))
        if nbins > 2**24:
            nbins = 2**24
            if verbose:
                print("Maximum number of bins (2^24) exceeded. Reset to 2^24.")
        else:
            self.binsToKill = int(nbins-nbins_old)
        
        # TODO: Put control on input params
        if self.kind == 'psd':
            self.time, self.rate = lcpsd(dt=dt, nbins=nbins, mean=mean, rms=rms, models=self.models)
        
        elif self.kind == 'coherent':
            self.time, self.rate = lcsinusoid(dt=dt, nbins=nbins, mean=mean, freq=freq, amp=amp, phi=phi, nha=nha)
        
        self.time = self.time[:-self.binsToKill]
        self.rate = self.rate[:-self.binsToKill]
        
        
    def poissonRandomize(self, dt, bkg):
        """
        Add Poissonian noise to lightcurve (and background, if present)
        """
        self.rate = poisson_randomization(self.rate, dt=dt, bkg=bkg)
        
    
    def getPSDModel(self, dt, nbins, freq=1000):
        """
        Get PSD model
        Returns a tuple with: 
        frequency array, total model, array with single components
        """
        # Sanity check
        assert self.kind == 'psd', 'ERROR! You can get models only if simulation kind is PSD'
        try:
            basestring
        except NameError:
            basestring = str
                
        # Get frequency array
        f_min = 1/(dt*nbins)
        f_max = 0.5/dt
        f = np.linspace(f_min, f_max, freq)
        p_tot = np.zeros_like(f)
        p_components = []
        for single_model in self.models:
            model  = single_model[0] 
            params = single_model[1]
            # If it is a built-in model:
            if isinstance(model, basestring):
                p_c = eval(model + "(f, params)")
            # If it is an astropy.modeling.models object:
            elif isinstance(model, type(astropy.modeling.models.Gaussian1D)):
                assert type(params) == dict, 'ERROR: params should be a dict, for astropy.modeling models'
                modelingfunc = model(**params)
                p_c = modelingfunc(f)   
            # If it is an user-defined model (i.e. a function object):
            else:
                p_c = model(f, params)
            p_components.append(p_c)
            p_tot += p_c
        
        return f, p_tot, p_components        
            
    
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

