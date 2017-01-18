from __future__ import division

def smoothbknpo_n(x, p):
    """
    Generalised smooth broken power law.
    Note: for p[4] == 2 the function is the same as 'smoothbknpo'

    Args:
        x:      (non-zero) frequencies.

        p[0]:   Normalization.

        p[1]:   power law index for f--> 0.

        p[2]:   power law index for f--> oo.

        p[3]:   break frequency.
    
        p[4]:   transition parameter (the higher the faster the transition).
    
    Returns:
        f:  psd model.

    History:
        v1:   Initial python implementation. Riccardo Campana, 2016.

    """
    f = p[0] * x**(-p[1]) * (1.+(x/p[3])**p[4])**((p[1]-p[2])/p[4])
    return f
	