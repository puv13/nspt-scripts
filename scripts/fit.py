#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
 This script simultaneously fits the volume dependence of all expansion
 coefficients.  Data are read from files fitting a user defined regex.
"""

import glob
import re
import numpy as np
import pcm_coeff as pcm
import fit_function_gen as FFG
import os.path
import time
from scipy.optimize import curve_fit
from optparse import OptionParser
from simultaneous_fit import sim_fit


#-------------------------------------------------------------------------------
# Define command line options
#-------------------------------------------------------------------------------
parser = OptionParser("fit.py [options]")
parser.add_option("-d", "--dir", dest="odir", default="./",
                  help="Working directory. Default is './'")

parser.add_option("-u", "--update", dest="fit_update_flag",
                  action="store_false",
                  default=True,
                  help="If the update flag is set, fits are only performed if "
                  "existing fit files are older than the input file.")

parser.add_option("-z", "--zero-epsilon-only", dest="zero_epsilon_flag",
                  action="store_true",
                  default=False,
                  help="If the zero epsilon flag is set, only data files "
                  "with epsilon = 0.0 are considered.")


#-------------------------------------------------------------------------------
# Regex to extract parameters from file names
#-------------------------------------------------------------------------------

regex = r"\S*N(?P<N>[0-9]{1,2})\.eps(?P<eps>[0-9]\.[0-9]*[eE][+-][0-9]*)\S*" \
       "_PCM_\S*ord(?P<ord>[0-9]{2})_all_V(?P<suf>\S*)\.stat"

regexc = re.compile(regex)

if __name__ == "__main__":
    #---------------------------------------------------------------------------
    # Parse options
    #---------------------------------------------------------------------------
    (options, args) = parser.parse_args()

    my_dir = options.odir.rstrip("/")
    update_flag = options.fit_update_flag
    zero_epsilon_flag = options.zero_epsilon_flag 

    glob_name = ""

    if my_dir:
        glob_name = my_dir+"/"+"*.stat"
    else:
        glob_name = "*.stat"

    #print update_flag
    if not update_flag:
        print "Existing fits will not be updated if younger than the data file."

    EXIT_SUCCESS = 0


    #---------------------------------------------------------------------------
    # Output helper function
    #---------------------------------------------------------------------------
    def write_to_file(ofname, known, coeffs, c_err, f_params, fp_err, chi_sq,
                      c_params=[], cp_err=[], bb1=False, quiet=False):
        """ Helper function to write the fitted coefficients and
            fit parameters to a file.

            Parameters
            ----------
            ofname : string
                Name of output file.
            known : int
                Number of known coefficients.
            coeffs : array_type
                Array holding the expansion coefficients.
            c_err : array_type
                Array holding the errors (from the fit) of the coefficients.
            f_params : array_type
                Array holding the fit function parameters f_i.
            fp_err : array_type
                Array holding the errors (from the fit) of the parameters f_i
            c_params : array_type
                Array holding the fit function parameters c_i.
            cp_err : array_type
                Array holding the errors (from the fit) of the parameters c_i
            bb1 : bool
                Flag to mark if the fit was done with beta function coefficient
                beta_1 set to zero (FALSE) or to its physical value (TRUE).
            quiet : bool
                Flag to toggle logging output of the function.

        """

        outfile = open(ofname, "w")
        if not quiet:
            print "\t\tWriting file {:s} ...".format(ofname)
        #outfile.write("# Data taken from '{:s}'\n".format(fname))
        outfile.write("# {:d} exact coefficients were used in the "\
                      "fits '\n".format(known-1))
        outfile.write("# chi2dof:{:4.3e}\n".format(chi_sq))
        outfile.write("# Columns: \n# n c_n  fit error c_n\n")

        # Write coefficients and errors to file
        for i in range(len(coeffs)):
            outfile.write("{:02d}\t{: 6.5e}\t{: 6.5e}\n".format(2*(i+1),
                                                                coeffs[i],
                                                                c_err[i]))

        # Write fit coefficients fn and errors to file
        outfile.write("\n\n# Fit coefficients f_n  (note that n here "
                      "corresponds to n/2 for the coefficients above, i.e., "
                      "we treat the expansion as an expansion in \\beta^-1, "
                      "not \\beta^{-1/2} .)\n")
        for i in range(len(f_params)):
            outfile.write("#{:02d}\t{: 6.5e}\t{: 6.5e}"\
                          "\n".format(i, f_params[i], fp_err[i]))

        if len(c_params) > 0:
            # Write fit coefficients cn and errors to file
            outfile.write("\n\n# Fit coefficients c_n  (note that n here "
                          "corresponds to n/2 for the coefficients above, "
                          "i.e., we treat the expansion as an "
                          "expansion in \\beta^-1, not \\beta^{-1/2} .)\n")
            for i in range(len(c_params)):
                outfile.write("#{:02d}\t{: 6.5e}\t{: 6.5e}"\
                              "\n".format(i, c_params[i], cp_err[i]))


        # Write fit parameters and functions in a gnuplot parsable  format
        gpl_string = "\n\n\n#gnuplot: b0=1/(8.*pi); "
        bb = "b0"
        if bb1:
            gpl_string += "b1=1./(128.*pi**2); "
            bb = "b1"

        for ff in range(len(f_params)+1):
            c_flag = (len(c_params) > 0)
            if ff < len(f_params):
                gpl_string += "f_{:s}_{:d} = {:10.9e}; ".format(bb,
                                                                ff,
                                                                f_params[ff])
            if c_flag and ff < len(c_params):
                gpl_string += "c_{:s}_{:d} = {:10.9e}; ".format(bb,
                                                                ff,
                                                                c_params[ff])


            if ff > 0:
                F_str = "F_{:s}_{:d}(x) = {: 10.9e} + (".format(bb,
                                                                ff,
                                                                coeffs[ff-1])
                F_str += FFG.epc_string(ff-1, bb1, c_flag)
                #print ff, F_str
                gpl_string += F_str+")/x**2; "


        outfile.write(gpl_string)


        outfile.close()
        if not quiet:
            print "\t\t...done!"


    #---------------------------------------------------------------------------
    # Fit function parameters
    #---------------------------------------------------------------------------
    b0 = 1./(8.*np.pi)
    b1 = 1./(128.*np.pi**2)

    #---------------------------------------------------------------------------
    # Get data files
    #---------------------------------------------------------------------------
    print "Processing data files ..."

    precomputed = 0
    for fname in glob.glob(glob_name):
        m = regexc.match(str(fname.rstrip()))
        d = {}


        if not m:
            continue

        print "\t", fname

        fname_age = 0
        if not update_flag:
            fname_age = os.path.getmtime(fname)


        d = m.groupdict()
        #print d
        order = int(d['ord'])
        N = int(d['N'])
        eps = float(d['eps'])
        suf = str(d['suf'])
        #print suf

        if zero_epsilon_flag and eps > 0.0:
            continue

        of_name_b0 = "N{:02d}_PCM_ord{:02d}_Vdep_fits_eps{:5.4e}{:s}"\
                     "_b0.dat".format(int(d['N']),
                                      int(d['ord']),
                                      float(d['eps']),
                                      suf)
        of_name_b0_c = "N{:02d}_PCM_ord{:02d}_Vdep_fits_eps{:5.4e}{:s}"\
                       "_b0_c.dat".format(int(d['N']),
                                          int(d['ord']),
                                          float(d['eps']),
                                          suf)

        of_name_b1 = "N{:02d}_PCM_ord{:02d}_Vdep_fits_eps{:5.4e}{:s}"\
                     "_b1.dat".format(int(d['N']),
                                      int(d['ord']),
                                      float(d['eps']),
                                      suf)

        of_name_b1_c = "N{:02d}_PCM_ord{:02d}_Vdep_fits_eps{:5.4e}{:s}"\
                       "_b1_c.dat".format(int(d['N']),
                                          int(d['ord']),
                                          float(d['eps']),
                                          suf)

        of_name_b0 = my_dir+"/"+of_name_b0
        of_name_b0_c = my_dir+"/"+of_name_b0_c
        of_name_b1 = my_dir+"/"+of_name_b1
        of_name_b1_c = my_dir+"/"+of_name_b1_c

        of_files = [of_name_b0, of_name_b1, of_name_b0_c, of_name_b1_c]

        fit_age = fname_age*2
        if not update_flag:
            #print "Checking mtimes"
            for aux in of_files:
                try:
                    fit_age = min(fit_age, os.path.getmtime(aux))
                # if file is not found
                except OSError:
                    fit_age = 0
                    break

        if fit_age > fname_age and not update_flag:
            print "\t\tFits are younger than file for '{:s}'".format(fname)
            print "\t\tSkipping fit as requested."
            continue


        # Load data
        data = np.loadtxt(fname, unpack=True, dtype=float)

        num_vs = len(data[0])/order
        x = data[0][:num_vs]
        if num_vs < 2:
            print "\tONLY DATA FOR ONE VOLUME IN  FILE '"+fname+"'!!!"
            print "\tWhat were you thinking? Skipping this one ...\n"
            continue

        # Pre-compute expansion coefficients for fit function
        if precomputed < order:
            print "\t\tPre-computing expansion coefficients. "\
                  "This might take a while ..."
            print "\t\t... done!"
            FFG.expansion_coefficient(order/2, True)
            precomputed = order


        # Get rid of "zero" coefficients
        y = np.array([data[2][2*i*num_vs:(2*i+1)*num_vs]
                      for i in range(1, order/2+1)]).flatten()
        y_err = np.array([data[3][2*i*num_vs:(2*i+1)*num_vs]
                          for i in range(1, order/2+1)]).flatten()

        # ----------------------------------------------------------------------
        # Generate fit function for simultaneous fitting
        # ----------------------------------------------------------------------
        # First coefficient is identically 0 for PCM (by construction)
        known_coeffs = 1
        if eps == 0.0:
            # In the epsilon -> 0 limit the first four coefficients are known
            known_coeffs = 4

        FF_b0 = map(lambda i: lambda xx, *myargs:
                    FFG.fit_func_b0(xx,
                                    N,
                                    i+1,
                                    order,
                                    b0,
                                    known_coeffs,
                                    *myargs),
                    range(len(y)/num_vs))

        FF_b1 = map(lambda i: lambda xx, *myargs:
                    FFG.fit_func_b1(xx,
                                    N,
                                    i+1,
                                    order,
                                    b0,
                                    b1,
                                    known_coeffs,
                                    *myargs),
                    range(len(y)/num_vs))

        FF_b0_c = map(lambda i: lambda xx, *myargs:
                      FFG.fit_func_b0_c(xx,
                                        N,
                                        i+1,
                                        order,
                                        b0,
                                        known_coeffs,
                                        *myargs),
                      range(len(y)/num_vs))

        FF_b1_c = map(lambda i: lambda xx, *myargs:
                      FFG.fit_func_b1_c(xx,
                                        N,
                                        i+1,
                                        order,
                                        b0,
                                        b1,
                                        known_coeffs,
                                        *myargs),
                      range(len(y)/num_vs))


        try:
            num_params = len(y)/num_vs-(known_coeffs-1)+order/2
            popt_b0, pcov_b0, chi2_b0 = sim_fit(FF_b0,
                                                x,
                                                y,
                                                y_err,
                                                num_params,
                                                absolute_sigma=True)
            perr_b0 = np.sqrt(np.diag(pcov_b0))

        except TypeError as te:
            print "\t\tWARNING: TypeError in b0 fit"
            print "\t\t\t{:s}".format(str(te))
            print "\t\t"+fname
            continue

        except RuntimeError as rte:
            print "\t\tWARNING: RuntimeError in b0 fit"
            print "\t\t\t{:s}".format(str(rte))
            print "\t\t"+fname
            continue

        try:

            pp0 = np.zeros(len(y)/num_vs-(known_coeffs-1)+2*(order/2)-2)
            #pp0 = np.ones(len(y)/num_vs-(known_coeffs-1)+2*(order/2)-2)
            pp0[:len(y)/num_vs-(known_coeffs-1)+order/2] = popt_b0
            #print "\t", len(pp0),  len(y), num_vs
            popt_b0_c, pcov_b0_c, chi2_b0_c = sim_fit(FF_b0_c,
                                                      x,
                                                      y,
                                                      y_err,
                                                      pp0,
                                                      absolute_sigma=True,
                                                      method="trf",
                                                      max_nfev=600*len(pp0))
            perr_b0_c = np.sqrt(np.diag(pcov_b0_c))

        except TypeError as te:
            print "\t\tWARNING: TypeError in b0_c fit"
            print "\t\t\t{:s}".format(str(te))
            print "\t\t" + fname
            continue

        except RuntimeError as rte:
            print "\t\tWARNING: RuntimeError in b0_c fit"
            print "\t\t\t{:s}".format(str(rte))
            print "\t\t" + fname
            continue

        try:

            num_params = len(y)/num_vs-(known_coeffs-1)+order/2
            popt_b1, pcov_b1, chi2_b1 = sim_fit(FF_b1, x, y, y_err,
                                                num_params,
                                                absolute_sigma=True)
            perr_b1 = np.sqrt(np.diag(pcov_b1))

        except TypeError as te:
            print "\t\tWARNING: TypeError in b1 fit"
            print "\t\t\t{:s}".format(str(te))
            print "\t\t"+fname
            continue

        except RuntimeError as rte:
            print "\t\tWARNING: RuntimeError in b1 fit"
            print "\t\t\t{:s}".format(str(rte))
            print "\t\t"+fname
            continue

        try:

            pp0 = np.zeros(len(y)/num_vs-(known_coeffs-1)+2*(order/2)-2)
            #pp0 = np.ones(len(y)/num_vs-(known_coeffs-1)+2*(order/2))
            pp0[:len(y)/num_vs-(known_coeffs-1)+order/2] = popt_b1
            popt_b1_c, pcov_b1_c, chi2_b1_c = sim_fit(FF_b1_c,
                                                      x,
                                                      y,
                                                      y_err,
                                                      pp0,
                                                      absolute_sigma=True,
                                                      method="trf",
                                                      max_nfev=600*len(pp0))
            perr_b1_c = np.sqrt(np.diag(pcov_b1_c))

        except TypeError as te:
            print "\t\tWARNING: TypeError in b1_c fit"
            print "\t\t\t{:s}".format(str(te))
            print "\t\t"+fname
            continue

        except RuntimeError as rte:
            print "\t\tWARNING: RuntimeError in b1_c fit"
            print "\t\t\t{:s}".format(str(rte))
            print "\t\t"+fname
            continue



        # Coefficients and Errors
        cfs_b0 = np.zeros(len(y)/num_vs)
        cfs_b0_err = cfs_b0.copy()
        cfs_b0_c = np.zeros(len(y)/num_vs)
        cfs_b0_err_c = cfs_b0_c.copy()
        cfs_b1 = np.zeros(len(y)/num_vs)
        cfs_b1_err = cfs_b1.copy()
        cfs_b1_c = np.zeros(len(y)/num_vs)
        cfs_b1_err_c = cfs_b1_c.copy()

        # Set exact coefficients if used
        for i in range(known_coeffs-1):
            cfs_b0[i] = pcm.coef(2*(i+1), N)
            cfs_b1[i] = pcm.coef(2*(i+1), N)
            cfs_b0_c[i] = pcm.coef(2*(i+1), N)
            cfs_b1_c[i] = pcm.coef(2*(i+1), N)


        coeffs_idx = order/2-known_coeffs+1
        cfs_b0[known_coeffs-1:] = popt_b0[:coeffs_idx]
        cfs_b1[known_coeffs-1:] = popt_b1[:coeffs_idx]
        cfs_b0_c[known_coeffs-1:] = popt_b0_c[:coeffs_idx]
        cfs_b1_c[known_coeffs-1:] = popt_b1_c[:coeffs_idx]

        cfs_b0_err[known_coeffs-1:] = perr_b0[:coeffs_idx]
        cfs_b1_err[known_coeffs-1:] = perr_b1[:coeffs_idx]
        cfs_b0_err_c[known_coeffs-1:] = perr_b0_c[:coeffs_idx]
        cfs_b1_err_c[known_coeffs-1:] = perr_b1_c[:coeffs_idx]


        # f_i and Errors
        fs_b0 = popt_b0[coeffs_idx:]
        fs_b0_err = perr_b0[coeffs_idx:]

        fs_b1 = popt_b1[coeffs_idx:]
        fs_b1_err = perr_b1[coeffs_idx:]

        fs_b0_c = popt_b0_c[coeffs_idx:2*(order/2)-known_coeffs+1]
        fs_b0_err_c = perr_b0_c[coeffs_idx:2*(order/2)-known_coeffs+1]

        fs_b1_c = popt_b1_c[coeffs_idx:2*(order/2)-known_coeffs+1]
        fs_b1_err_c = perr_b1_c[coeffs_idx:2*(order/2)-known_coeffs+1]


        # # c_i and Errors
        cc_b0 = popt_b0_c[2*(order/2)-known_coeffs+1:]
        cc_b0_err = perr_b0_c[2*(order/2)-known_coeffs+1:]

        cc_b1 = popt_b1_c[2*(order/2)-known_coeffs+1:]
        cc_b1_err = perr_b1_c[2*(order/2)-known_coeffs+1:]


        # Append c[ORD/2-2], which we set to zero by hand
        cc_b0 = np.append(cc_b0, 0.0)
        cc_b0_err = np.append(cc_b0_err, 0.0)

        cc_b1 = np.append(cc_b1, 0.0)
        cc_b1_err = np.append(cc_b1_err, 0.0)


        # Write output files
        write_to_file(of_name_b0, known_coeffs,
                      cfs_b0, cfs_b0_err, fs_b0, fs_b0_err, chi2_b0)

        write_to_file(of_name_b1, known_coeffs,
                      cfs_b1, cfs_b1_err, fs_b1, fs_b1_err, chi2_b1, bb1=True)

        write_to_file(of_name_b0_c, known_coeffs,
                      cfs_b0_c, cfs_b0_err_c, fs_b0_c, fs_b0_err_c, chi2_b0_c,
                      cc_b0, cc_b0_err)

        write_to_file(of_name_b1_c, known_coeffs,
                      cfs_b1_c, cfs_b1_err_c, fs_b1_c, fs_b1_err_c, chi2_b1_c,
                      cc_b1, cc_b1_err, bb1=True)



    print "... done!"
    exit(EXIT_SUCCESS)
    # print popt
    #  print popt_tst
    # print popt_tst_b1



### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
