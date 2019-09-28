#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This is a script to parse NSPT data files for the PCM model.

    The script looks for files containing NSPT data and collects results for
    fixed lattice volume and fixed stochastic time step in two separate output
    files.  If several time steps epsilon are available for a given volume an
    extrapolation epsilon to zero is also performed. See comments for
    implementation details.
"""

import glob
import re
import numpy as np
import warnings
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from optparse import OptionParser
from math import sqrt, ceil


#-------------------------------------------------------------------------------
# Define command line options
#-------------------------------------------------------------------------------
parser = OptionParser("PCM_analyser.py [options]")
parser.add_option("-d", "--dir", dest="odir", default="",
                  help="Working directory. Default is './'")


#-------------------------------------------------------------------------------
# Regex to extract parameters from file names
#-------------------------------------------------------------------------------

regex = r"\S*N(?P<N>[0-9]{1,2})\.(?P<LT>[0-9]{1,2})x(?P<LS>[0-9]{1,2})\S*_PCM_"\
  "\S*ord(?P<ord>[0-9]{2})\Seps(?P<eps>[0-9]\.[0-9]*[eE][+-][0-9]*)\.\S*\.stat"

regexc = re.compile(regex)

#-------------------------------------------------------------------------------
# Define fit function
#-------------------------------------------------------------------------------

def fitfunc_e2(epsilon, a, b):
    """
       Fit function for epsilon dependence of NSPT results.

       By construction the Runge-Kutta used in the calculations the
       part linear in epsilon is cancelled.

       Parameters
       ----------
       epsilon : float
           NSPT time step epsilon
       a : float
           Fit parameter for quadratic term
       b : float
           Fit parameter for constant term

       Returns
       -------
       a*epsilon**2+b : float

    """
    return a*epsilon**2 + b

# def fitfunc_e3(eps, a, b, c):
#     """
#        Fit function for epsilon dependence of NSPT results.

#        By construction the Runge-Kutta used in the calculations the
#        part linear in epsilon is cancelled.

#        Parameters
#        ----------
#        eps : float
#            NSPT time step epsilon
#        a : float
#            Fit parameter for cubic term
#        b : float
#            Fit parameter for quadratic term
#        c : float
#            Fit parameter for constant term

#        Returns
#        -------
#        a*eps**3 + b*eps**2 +c : float
#    """

if __name__ == "__main__":
    #---------------------------------------------------------------------------
    # Parse options
    #---------------------------------------------------------------------------
    (options, args) = parser.parse_args()

    my_dir = options.odir.rstrip("/")

    glob_name = ""

    if my_dir:
        glob_name = my_dir+"/"+"*.stat"
    else:
        glob_name = "*.stat"

    EXIT_SUCCESS = 0


    #---------------------------------------------------------------------------
    # Get data files
    #---------------------------------------------------------------------------
    print "Parsing data files ..."

    para_dict = {}
    vol_dict = {}

    #print glob_name
    for fname in glob.glob(glob_name):
        #print fname
        m = regexc.match(str(fname.rstrip()))
        d = {}


        if not m:
            continue

        d = m.groupdict()

        N = int(d["N"])
        LS = int(d["LS"])
        LT = int(d["LT"])
        L = (LT, LS)
        eps = float(d["eps"])
        order = int(d["ord"])

        #-----------------------------------------------------------------------
        # Build dict for files with fixed volume
        #-----------------------------------------------------------------------
        if N not in para_dict:
            para_dict[N] = {}

        if L not in para_dict[N]:
            para_dict[N][L] = {}

        if order not in para_dict[N][L]:
            para_dict[N][L][order] = []

        if eps not in para_dict[N][L][order]:
            para_dict[N][L][order].append([eps, fname])


        #-----------------------------------------------------------------------
        # Build dict for files with fixed epsilon
        #-----------------------------------------------------------------------

        if N not in vol_dict:
            vol_dict[N] = {}

        if eps not in vol_dict[N]:
            vol_dict[N][eps] = {}

        if order not in vol_dict[N][eps]:
            vol_dict[N][eps][order] = []

        if L not in vol_dict[N][eps][order]:
            vol_dict[N][eps][order].append([L, fname])

        #print para_dict

        #print "-"*90
        #print fname
        #print d


    print "... done!\n"

    #---------------------------------------------------------------------------
    # Writing files for fixed volume
    #---------------------------------------------------------------------------
    print "Writing files for fixed volume ... "
    for Nkey in sorted(para_dict.keys()):
        for Ls in  sorted(para_dict[Nkey].keys()):
            for Ord in  para_dict[Nkey][Ls]:
                # For fixed N, Ls, Lt and order sort by eps
                lst = sorted(para_dict[Nkey][Ls][Ord], key=lambda ls: ls[0])

                outfile_name = "N{:02d}.{:02d}x{:02d}_PCM_ord{:02d}"\
                               "_all_eps.stat".format(Nkey, Ls[0], Ls[1], Ord)

                if my_dir:
                    outfile_name = my_dir+"/"+outfile_name

                # Sort data for different eps by order of coefficient
                collected = {}
                for fobj in lst:
                    #print fobj[1]
                    for line in open(fobj[1], "r"):
                        if line[0] == "#":
                            continue
                        line = line.split()

                        oo = int(line[0])
                        if oo < 1 or oo > Ord:
                            continue

                        if oo not in collected:
                            collected[oo] = []

                        #print line
                        if float(line[6]) >= 0:
                            tau = sqrt(2.*float(line[6]))
                        else:
                            tau = sqrt(2.*float(line[5]))
                            print "\tWARNING: Encountered "\
                                "negative auto-correlation time estimate"
                            print "\t\t file: {:s} coef: {:02d}".format(fobj[1],
                                                                        oo)


                        collected[oo].append([float(fobj[0]),
                                              float(line[1]),
                                              float(line[2])*tau])



                print "\t", outfile_name
                outfile = open(outfile_name, "w")

                header = "# File generated from data in\n"
                for f in [x[1] for x in lst]:
                    header += "#\t{:s}\n".format(f)

                header += "#\n# Cols: eps\tcoef\tcoef err "\
                          "(inc. autocorrelation)\t(syst. error)\n"

                outfile.write(header)

                for key in collected:
                    res = []
                    for entry in collected[key]:
                        outfile.write("{: 8.6e}\t{: 8.6e}\t{: 8.6e}"\
                                      "\t 0.0 \n".format(*entry))
                        res.append(entry)
                    res = np.array(res).T

                    # Check for inconsistencies in input.
                    if 1 == key and len(res[0]) != len(np.unique(res[0])):
                        print "WARNNING: Epsilon values are not unique"\
                              " for {:s}".format(outfile_name)

                    # Fit linear (in eps^2) dependence of the coefficients
                    if len(res[0]) > 1:
                        outfile.write("# Linear (in \epsilon^2) fit "\
                                      "extrapolation to \epsilon=0\n")
                        # Check that all sigmas/erros are finite ...
                        if np.all(res[2]):

                            # Treat optimisation warnings from curve_fit
                            # as errors. We want to be thorough here.
                            with warnings.catch_warnings():
                                warnings.simplefilter("error", OptimizeWarning)
                                try:
                                    popt, pcov = curve_fit(fitfunc_e2,
                                                           res[0],
                                                           res[1],
                                                           sigma=res[2],
                                                           absolute_sigma=True)

                                    perr = np.sqrt(np.diag(pcov))

                                except OptimizeWarning:
                                    print "Error: "\
                                          "curve_fit throws OptimizeWarning"
                                    exit(-1)

                        # ... otherwise fit without errors
                        else:
                            # This should only happen in the leading orders
                            if key > 2:
                                print "\tWARNING: Strange errors for"\
                                      " order {:2d}".format(key)

                            # Fit with no errors
                            popt, pcov = curve_fit(fitfunc_e2, res[0], res[1])
                            perr = np.sqrt(np.diag(pcov))



                       # Estimate systematic error:
                       # Ratio between the result for (median epsilon [the
                       # smaller one if number ob epsilons is even]) over (fit
                       # result) should be 1 in the limit where (median epsilon)
                       # -> 0. We take the difference from 1 as a measure for
                       # the syst. error (following Bali et a. from PRD 87,
                       # 094517 (2013))

                        pn = res[1][int(ceil(len(res[1])/2.)-1)]
                        dn = 1
                        if pn != 0:
                            dn = abs(1 - popt[1]/pn)
                        outfile.write("{: 8.6e}\t{: 8.6e}\t{: 8.6e}"\
                                      "\t{: 8.6e}\n".format(0,
                                                            popt[1],
                                                            perr[1],
                                                            abs(pn*dn)))

                    outfile.write("\n\n")

                outfile.close()
    print "... done!\n"
    #---------------------------------------------------------------------------
    # Writing files for fixed eps
    #---------------------------------------------------------------------------

    print "Writing files for fixed eps ... "
    for Nkey in sorted(vol_dict.keys()):
        for epskey in sorted(vol_dict[Nkey].keys()):
            for Ord in  vol_dict[Nkey][epskey]:
                # For fixed N, eps and order sort by  V=Ls*Lt
                lst = sorted(vol_dict[Nkey][epskey][Ord],
                             key=lambda ls: ls[0][0]*ls[0][1])

                outfile_name = "N{:02d}.eps{:10.4e}_PCM_ord{:02d}"\
                             "_all_V.stat".format(Nkey,
                                                  float(epskey),
                                                  Ord)

                if my_dir:
                    outfile_name = my_dir+"/"+outfile_name

                print "\t", outfile_name

                # Sort data for different Vol by order of coefficient
                collected = {}
                for fobj in lst:
                    for line in open(fobj[1], "r"):
                        if line[0] == "#":
                            continue
                        line = line.split()

                        oo = int(line[0])
                        if oo < 1 or oo > Ord:
                            continue

                        if oo not in collected:
                            collected[oo] = []

                        if float(line[6]) >= 0:
                            tau = sqrt(2.*float(line[6]))
                        else:
                            tau = sqrt(2.*float(line[5]))
                            print "\tWARNING: Encountered negative "\
                                  "autocorrelation time estimate"
                            print "\t\t file: {:s} coef: {:02d}".format(fobj[1],
                                                                        oo)

                        collected[oo].append([fobj[0][0],
                                              fobj[0][1],
                                              float(line[1]),
                                              float(line[2])*tau])

                outfile = open(outfile_name, "w")

                header = "# File generated from data in\n"
                for f in [x[1] for x in lst]:
                    header += "#\t{:s}\n".format(f)

                header += "#\n# Cols: Lt\tLs\tcoef"\
                          "\tcoef err (inc. autocorrelation)\n"

                outfile.write(header)

                for key in collected:
                    for entry in collected[key]:
                        outfile.write("{: 3d}\t{: 3d}\t{: 8.6e}"\
                                      "\t{: 8.6e}\n".format(*entry))

                    outfile.write("\n\n")

                outfile.close()

    print "... done!"
    #print vol_dict[12]
    exit(EXIT_SUCCESS)

### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
