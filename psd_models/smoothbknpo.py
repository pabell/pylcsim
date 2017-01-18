from __future__ import division

def smoothbknpo(x, p):
    """
    Smooth broken power law.

    Args:
        x:      (non-zero) frequencies.

        p[0]:   Normalization.

        p[1]:   power law index for f--> 0.

        p[2]:   power law index for f--> oo.

        p[3]:   break frequency.

    Returns:
        f:  psd model.

    History:
        v1:   Initial python implementation. Riccardo Campana, 2014.

    """
    f = p[0]*x**(-p[1])/(1.+(x/p[3])**2.)**(-(p[1]-p[2])/2.)
    return f
	