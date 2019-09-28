#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   This module provides a function to calculate the known expansion coefficients
   of the energy density of the Principal Chiral model.

   References
   ----------
   Paolo Rossi and Ettore Vicari,
   "Two-dimensional SU(N)xSU(N) chiral models on the lattice",
   Phys. Rev. D 49(3), pp. 1621--1628, 1994
   doi : 10.1103/PhysRevD.49.1621
   arxiv : hep-lat/9307014v2
"""


#---------------------------------------------------------------------------
# Known PCM coefficients (\todo add ref!)
#---------------------------------------------------------------------------

# Numerical constants
Q1 = 0.0958876
Q2 = -0.0670

def a1(N):
    """
       A helper function for coef(i,N)
    """
    return (N**2-2.)/(32.*N**2)

def a2(N):
    """
       A helper function for coef(i,N)
    """
    res = (3.*N**4-14.*N**2+20.)/(768.*N**4.)
    res += (N**4.-4.*N**2+12.)*Q1/(64.*N**4) + \
           (N**4-8.*N**2+24.)*Q2/(64.*N**4)
    return res


def coef(i, N):
    """
       Returns the expansion coefficient of order i (in \beta^{1/2}) of the
       energy density in the PCM for rank N.

       Note that the coefficients are only known up to order 6!
    """
    if i%2 != 0 or float(i) == 0.0:
        return 0.0
    elif i > 6:
        return 1./0.
    else:
        const = (N**2-1.)/(8.*N**2*1.)
        if i < 4:
            return const
        elif i < 6:
            return const*a1(N)
        else:
            return const*a2(N)

### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
