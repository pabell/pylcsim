Documentation for pyLCSIM
**************************

============
Introduction
============
pyLCSIM is python package to simulate X-ray lightcurves from coherent signals and power spectrum models.

Coherent signals can be simulated as a sum of one or more sinusoids, each with its frequency, pulsed fraction or phase shift; or as a series of harmonics of a fundamental frequency.

Power spectra can be simulated from a model of the power spectrum density (PSD), using as a template one or more library functions. The user can also define his/her custom models.
Models are additive.

Warning: the current release (0.1) is HIGHLY EXPERIMENTAL!
Use at your risk...


============
Installation
============

TBC


============
Examples
============
Let's begin with a PSD model simulation.

We import the usual packages::

	import matplotlib.pyplot as plt
	import numpy as np
	import pyLCSIM

We assume that our source has a rate of 300000 counts/s. Moreover, we have a 3000 counts/s background rate, a 50 s exposure time and our observation has a time resolution of 10 ms.
We want to simulate a QPO at a frequency of 10 Hz, superimposed to a continuum modelled as a smoothly-varying broken power law.::

	rate_src    = 300000.0
	rate_bkg    = 3000.0
	t_exp       = 50.0
	dt          = 0.01
	nbins = long(t_exp/dt)

Then, the simulation follows as::

	# Instantiate a simulation object
	sim = pyLCSIM.Simulation()

	# Add two PSD models: a smooth broken power law and a Lorentzian representing a QPO.
	sim.addModel('smoothbknpo', [1., 2, 2, 10])
	sim.addModel('lorentzian', [10., 1., 100., 10])

	# Run the simulation
	sim.run(dt, nbins, rate_src, rms=0.01)

	# Add Poisson noise to the light curve
	sim.poissonRandomize(dt, rate_bkg)

	# Get lightcurve and power spectrum as 1-D arrays
	time, rate = sim.getLightCurve()
	f, psd = sim.getPowerSpectrum()

Done! We can view the results as::

	# Plot the lightcurve and power spectrum
	fig0 = plt.figure()
	plt.plot(time, rate)

	fig1 = plt.figure()
	plt.loglog(f, psd, drawstyle='steps-mid', color='black')

	# Save FITS files with lightcurve and spectrum
	pyLCSIM.saveFITSLC("myLC.fits", time, rate)
	pyLCSIM.saveFITSPSD("myPSD.fits", f, psd)

	plt.show()
	
The following is an example using a sum of sinusoids::

	import matplotlib.pyplot as plt
	import numpy as np
	import pyLCSIM

	rate_src    = 300000.0
	rate_bkg    = 3000.0
	t_exp       = 1.0
	dt          = 0.0001
	nbins = long(t_exp/dt)

	print nbins

	# Instantiate a simulation object, this time as coherent
	sim = pyLCSIM.Simulation(kind='coherent')


	# Run the simulation:
	# four sinusoidal frequencies: 340, 550, 883, 1032 Hz
	# with pulsed fractions 1%, 5%, 7% and 15% respectively
	# The third frequency has a 35 degree phase shift with respect to the others
	sim.run(dt, nbins, rate_src, freq=[340, 550, 883, 1032], amp=[0.01, 0.05, 0.07, 0.15], phi=[0., 0, 35., 0.])

	# Add Poisson noise to the light curve
	sim.poissonRandomize(dt, rate_bkg)

	# Get lightcurve and power spectrum as 1-D arrays
	time, rate = sim.getLightCurve()
	f, psd = sim.getPowerSpectrum()

	# Plot the lightcurve and power spectrum
	fig0 = plt.figure()
	plt.plot(time, rate)

	fig1 = plt.figure()
	plt.loglog(f, psd, drawstyle='steps-mid', color='black')

	# Save FITS files with lightcurve and spectrum
	pyLCSIM.saveFITSLC("myLC.fits", time, rate)
	pyLCSIM.saveFITSPSD("myPSD.fits", f, psd)

	plt.show()


Finally, an example with a fundamental frequency and two harmonics::

	import matplotlib.pyplot as plt
	import numpy as np
	import pyLCSIM

	rate_src    = 300000.0
	rate_bkg    = 3000.0
	t_exp       = 1.0
	dt          = 0.0001
	phase_shift = 0.
	nbins = long(t_exp/dt)

	print nbins

	# Instantiate a simulation object, this time as coherent
	sim = pyLCSIM.Simulation(kind='coherent')


	# Run the simulation:
	# Fundamental at 500 Hz, 3 harmonics (500, 1000, 1500 Hz)
	# with pulsed fractions 10%, 5% and 15% respectively
	sim.run(dt, nbins, rate_src, freq=500, nha=3, amp=[0.1, 0.05, 0.15])

	# Add Poisson noise to the light curve
	sim.poissonRandomize(dt, rate_bkg)

	# Get lightcurve and power spectrum as 1-D arrays
	time, rate = sim.getLightCurve()
	f, psd = sim.getPowerSpectrum()

	# Plot the lightcurve and power spectrum
	fig0 = plt.figure()
	plt.plot(time, rate)

	fig1 = plt.figure()
	plt.loglog(f, psd, drawstyle='steps-mid', color='black')

	# Save FITS files with lightcurve and spectrum
	pyLCSIM.saveFITSLC("myLC.fits", time, rate)
	pyLCSIM.saveFITSPSD("myPSD.fits", f, psd)

	plt.show()



===========
Main module
===========

.. automodule:: pyLCSIM
   :members:

=====================
Submodule: psd_models
=====================

.. automodule:: pyLCSIM.psd_models
   :members:
