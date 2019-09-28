#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   Script to extract asymptotic renormalon behaviour from Principal Chiral model
   data for different rank N.

   In the limit of n \to \infty the ratio r_n := E_n/(E_(n-1)*n) should approach
   a fixed value, independent of N.  We can therefore (for large n) fit a
   constant to the values r_n(N) for different N.

   As a cross check we also compute the mean over N of the r_n(N) (with proper
   propagation of uncertainties).
"""

import glob
import re
import numpy as np
from scipy.optimize import curve_fit
from optparse import OptionParser


#-------------------------------------------------------------------------------
# Define command line options
#-------------------------------------------------------------------------------
parser = OptionParser("fit.py [options]")
parser.add_option("-d", "--dir", dest="odir", default="./",
                  help="Working directory. Default is './'")


#-------------------------------------------------------------------------------
# Define "ratio" and error propagation
#-------------------------------------------------------------------------------
def ratio(xx, yy):
    """
       Returns the ratio xx/yy
    """
    return xx/yy

def eratio(xx, yy, dxx, dyy):
    """
       Returns the uncertainty in the ratio xx/yy.
       (Gaussian uncertainty propagation)
    """
    res = dxx**2*1./(yy)**2 + dyy**2*(xx*1./yy**2)**2
    return np.sqrt(res)

def rat_err(xx, yy, dxx, dyy, nn):
    """
        Returns the ratio xx/(yy*nn) and the corresponding uncertainty.
        (Gaussian uncertainty propagation)
    """
    rat_ = ratio(xx, yy*nn)
    err_ = eratio(xx, nn*yy, dxx, nn*dyy)
    return rat_, err_

#-------------------------------------------------------------------------------
# Define "fit function"
#-------------------------------------------------------------------------------
def const(xx, a):
    """
       Returns the constant a for any xx.
    """
    return a

#-------------------------------------------------------------------------------
# Regex to extract parameters from file names
#-------------------------------------------------------------------------------

regex = r"\S*Ratios_N(?P<N>[0-9]{1,2})_Vfit_PCM_ord(?P<ord>[0-9]{2})_" \
       "eps(?P<eps>[0-9]\.[0-9]*[eE][+-][0-9]*)(?P<suf>\S*)\.dat"

regexc = re.compile(regex)




if __name__ == "__main__":
    #---------------------------------------------------------------------------
    # Parse options
    #---------------------------------------------------------------------------
    (options, args) = parser.parse_args()

    my_dir = options.odir.rstrip("/")

    glob_name = ""

    if my_dir:
        glob_name = my_dir+"/"+"*.dat"
    else:
        glob_name = "*.dat"

    EXIT_SUCCESS = 0

    #---------------------------------------------------------------------------
    # Get data files
    #---------------------------------------------------------------------------
    print "Processing data files ..."

    resdict = {}
    for fname in glob.glob(glob_name):
        m = regexc.match(str(fname.rstrip()))
        d = {}

        if not m:
            continue

        d = m.groupdict()

        N = int(d["N"])
        eps = float(d["eps"])
        order = int(d["ord"])
        suf = str(d["suf"])


        # Load data
        data = np.loadtxt(fname, unpack=True)

        ratios, errors = rat_err(data[1], data[2], data[3], data[4], data[0]/2.)

        out = zip(ratios, errors)

        print "\t", fname

        for n in data[0]/2:

            n = int(n)-2
            #print n
            try:
                resdict[(suf, eps, n+1)].append(out[n])# [fname, out[n]])

            except KeyError:
                resdict[(suf, eps, n+1)] = [out[n]]#[fname, out[n]]]



    #tst =resdict[("_b0", 0.000, 3)]
    #print tst, len(tst)



outdict = {}
for key in resdict:

    suf, eps, n = key
    entries = resdict[key]
    #print "*"*60
    #print  entries
    rat, err = zip(*entries)
    rat = np.array(rat)
    err = np.array(err)
    #dummy x for fit
    x = np.array(range(len(rat)))


    # print rat
    # print err
    popt = 0
    pcov = 0
    if np.any(err):
        popt, pcov = curve_fit(const, x, rat, sigma=err, absolute_sigma=True)
    else:
        popt, pcov = curve_fit(const, x, rat)
    r = popt[0]
    e = np.sqrt(np.diag(pcov))[0]
    m = np.mean(rat)
    s = np.sqrt(np.sum(err*err))/len(err)

    try:
        outdict[(suf, eps)].append([n+1, r, e, m, s])
    except KeyError:
        outdict[(suf, eps)] = [[n+1, r, e, m, s]]



head = "# r_n \t Fit \t FitError \t Mean \t StdErr\n"
for key in outdict:
    suf, eps = key
    ofname = "CombinedRatio_eps{:6.4e}{:s}.dat".format(eps, suf)


    export = sorted(outdict[key], key=lambda a: a[0])
    export = np.array(export)
    np.savetxt(ofname,
               export,
               fmt="%02d\t%8.6e\t%8.6e\t%8.6e\t%8.6e",
               delimiter="\t",
               header=head)




### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
