#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This is a script to extract the epsilon to 0 results from NSPT data files
    for the PCM model.
"""

import glob
import re
import numpy as np
from optparse import OptionParser


#-------------------------------------------------------------------------------
# Define command line options
#-------------------------------------------------------------------------------
parser = OptionParser("Make_eps0.00_allV.py [options]")
parser.add_option("-d", "--dir", dest="odir", default="",
                  help="Working directory. Default is './'")



#-------------------------------------------------------------------------------
# Regex to extract parameters from file names
#-------------------------------------------------------------------------------

#"N03.08x08_PCM_ord09_all_eps.stat"
regex = r"\S*N(?P<N>[0-9]{1,2}).(?P<LT>[0-9]{1,2})x(?P<LS>[0-9]{1,2})"\
        "\S*_PCM_\S*ord(?P<ord>[0-9]{2})\S*all_eps.stat"
regexc = re.compile(regex)



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

    #para_dict = {}
    vol_dict = {}

    #print glob_name
    for fname in glob.glob(glob_name):
        #print fname
        m = regexc.match(str(fname.rstrip()))
        d = {}


        if not m:
            continue

        d = m.groupdict()
        # print fname
        # print d

        N = int(d["N"])
        LS = int(d["LS"])
        LT = int(d["LT"])
        eps = 0.0
        L = (LT, LS)
        order = int(d["ord"])

        #-----------------------------------------------------------------------
        # Build dict for files
        #-----------------------------------------------------------------------

        if N not in vol_dict:
            vol_dict[N] = {}


        if order not in vol_dict[N]:
            vol_dict[N][order] = []

        if L not in vol_dict[N][order]:
            vol_dict[N][order].append([L, fname])


        #print vol_dict

    #---------------------------------------------------------------------------
    # Writing files
    #---------------------------------------------------------------------------

    print "Writing files for fixed eps=0.0 ... "
    for Nkey in sorted(vol_dict.keys()):
        for Ord in  vol_dict[Nkey]:
            # For fixed N, eps and order sort by  V = Ls*Lt
            lst = sorted(vol_dict[Nkey][Ord],
                         key=lambda ls: ls[0][0]*ls[0][1])

            outfile_name = "N{:02d}.eps{:10.4e}_PCM_ord{:02d}"\
                           "_all_V.stat".format(Nkey, eps, Ord)
            outfile_sys_name = "N{:02d}.eps{:10.4e}_PCM_ord{:02d}"\
                               "_all_V_sys.stat".format(Nkey, eps, Ord)

            #print outfile_name
            if my_dir:
                outfile_name = my_dir+"/"+outfile_name
                outfile_sys_name = my_dir+"/"+outfile_sys_name

            print "\t", outfile_name
            print "\t", outfile_sys_name

            # Sort data for different Vol by order of coefficient
            collected = {}
            for fobj in lst:
                #print "Fobj: \t", fobj
                raw_data = np.loadtxt(fobj[1], dtype=float)
                # Filter for eps = 0.0
                raw_data = raw_data[np.where(raw_data[:, 0] == 0.0)].T

                # Combination of sys and stat errors added in quadrature
                tau = np.sqrt(raw_data[2]**2+raw_data[3]**2)
                # Stat error
                stat = raw_data[2]
                # Get rid of infinite errors for first two coeffs
                tau[tau == np.inf] = 0.0
                stat[stat == np.inf] = 0.0


                for i in range(len(raw_data[1])):
                    if i not in collected:
                        collected[i] = []

                    collected[i].append([fobj[0][0], fobj[0][1],
                                         float(raw_data[1][i]),
                                         float(tau[i]),
                                         float(stat[i])])


                outfile_sys = open(outfile_sys_name, "w")
                outfile = open(outfile_name, "w")

                header = "# Data for fits epsilon to 0\n"\
                         "# File generated from data in\n"
                for f in [x[1] for x in lst]:
                    header += "#\t{:s}\n".format(f)

                header_sys = header + "#\n# Cols: "\
                             "Lt\tLs\tcoef\tcoef err (inc. autocorr.+syst)\n"
                header += "#\n# Cols: "\
                          "Lt\tLs\tcoef\tcoef err (inc. autocorr., w/o syst)\n"

                outfile_sys.write(header_sys)
                outfile.write(header)

                for key in collected:
                    for entry in collected[key]:
                        output_sys = entry[:-1]
                        output = entry[:-2]
                        output.append(entry[-1])

                        outfile_sys.write("{: 3d}\t{: 3d}\t{: 8.6e}"
                                          "\t{: 8.6e}\n".format(*output_sys))
                        outfile.write("{: 3d}\t{: 3d}\t{: 8.6e}"
                                      "\t{: 8.6e}\n".format(*output))

                    outfile_sys.write("\n\n")
                    outfile.write("\n\n")

                outfile_sys.close()
                outfile.close()

    print "... done!"

### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
