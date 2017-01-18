from __future__ import division

def lorentzian(x, p):
    """
    Generalized Lorentzian function.

    (**WARNING**: for n != 2 the function is no more normalized!)

    Args:
        x:      (non-zero) frequencies.

        p[0]:   x0 = peak central frequency.

        p[1]:   gamma = FWHM of the peak.

        p[2]:   value of the peak at x = x0.

        p[3]:   n = power coefficient.

        The quality factor is given by Q = x0/gamma.

    Returns:
        f:  psd model.

    History:
        v1:   Initial python implementation. Riccardo Campana, 2014.

    """
    assert p[3] > 0., "The power coefficient should be greater than zero."
    f = p[2] * (p[1]/2.)**p[3] * 1./( abs(x - p[0])**p[3] + (p[1]/2.)**p[3] ) 
    return f
