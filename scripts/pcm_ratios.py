#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
   This is a script to parse NSPT data files for the PCM model and to generate
   new output files suitable for plotting the ratios of the expansion
   coefficients with gnuplot.
"""

import glob
import re
from optparse import OptionParser
import numpy as np


#-------------------------------------------------------------------------------
# Define command line options
#-------------------------------------------------------------------------------
parser = OptionParser("PCM_ratios.py [options]")
parser.add_option("-d", "--dir", dest="odir", default="",
                  help="Working directory. Default is './'")


#-------------------------------------------------------------------------------
# Regex to extract parameters from file names
#-------------------------------------------------------------------------------

regex = r"\N(?P<N>[0-9]{1,2}).(?P<LT>[0-9]{1,2})x(?P<LS>[0-9]{1,2})\S*_PCM_\S*"\
       "ord(?P<ord>[0-9]{2})\Seps(?P<eps>[0-9]\.[0-9]*[eE][+-][0-9]*).\S*\.stat"

regex_fit = r"\N(?P<N>[0-9]{1,2})_PCM\S*_ord(?P<ord>[0-9]{2})\S*"\
           "eps(?P<eps>[0-9]\.[0-9]*[eE][+-][0-9]*)(?P<suf>\S*)\.dat"

regexc = re.compile(regex)
regex_fitc = re.compile(regex_fit)


if __name__ == "__main__":
    #---------------------------------------------------------------------------
    # Parse options
    #---------------------------------------------------------------------------
    (options, args) = parser.parse_args()

    my_dir = options.odir.rstrip("/")

    glob_name = ""

    if my_dir:
        glob_name = my_dir+"/"+"*.*at"
    else:
        glob_name = "*.*at"

    EXIT_SUCCESS = 0


    #---------------------------------------------------------------------------
    # Get data files
    #---------------------------------------------------------------------------
    print "Parsing data files ..."

    # Fixed V and eps
    for fname in glob.glob(glob_name):
        m = regexc.match(str(fname.rstrip()))
        d = {}


        if not m:
            continue

        #print fname
        d = m.groupdict()

        N = int(d["N"])
        LS = int(d["LS"])
        LT = int(d["LT"])
        L = (LT, LS)
        eps = float(d["eps"])
        order = int(d["ord"])

        data = np.loadtxt(fname)

        data = np.array(data[3:order+1:2])

        data[:, 0] -= 1
        rdata = np.array(data[1:])
        data = data[:-1]

        # print data[:, 1]
        # print rdata[:, 1]

        #try:
        export = np.array([rdata[:, 0],
                           rdata[:, 1],
                           data[:, 1],
                           rdata[:, 2]*np.sqrt(rdata[:, 6]*2.),
                           data[:, 2]*np.sqrt(data[:, 6]*2.)])
        #except:


        head = " Original data in {:s}\n".format(fname)
        head += " n\tc_n\tc_(n-1)\terror c_n\terror c_(n-1)"
        export_name = "Ratios_N{:02d}.{:02d}x{:02d}_PCM_ord{:02d}"\
                      "_eps{:10.4e}.dat".format(N, LT, LS, order, eps)

        #print export_name
        np.savetxt(export_name,
                   export.transpose(),
                   fmt="%8.6e",
                   delimiter="\t",
                   header=head)



    # Fixed eps, V \to \infty fits
    for fname in glob.glob(glob_name):
        m = regex_fitc.match(str(fname.rstrip()))
        d = {}
        #print fname

        if not m:
            continue
        #print fname

        d = m.groupdict()

        N = int(d["N"])
        eps = float(d["eps"])
        order = int(d["ord"])
        suffix = str(d["suf"])

        data = np.loadtxt(fname)

        data = data.T
        export = np.array([data[0][1:],
                           data[1][1:],
                           data[1][:-1],
                           data[2][1:],
                           data[2][:-1]])

#        export = np.array(data[:, 1:-1],

        head = " Original data in {:s}\n".format(fname)
        head += " n\tc_n\tc_(n-1)\terror c_n\terror c_(n-1)"
        export_name = "Ratios_N{:02d}_Vfit_PCM_ord{:02d}_eps{:10.4e}{:s}"\
                      ".dat".format(N, order, eps, suffix)

        np.savetxt(export_name,
                   export.transpose(),
                   fmt="%8.6e",
                   delimiter="\t",
                   header=head)


    print "... done!"
    exit(EXIT_SUCCESS)


### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
