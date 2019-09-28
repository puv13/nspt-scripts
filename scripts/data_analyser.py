#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is a script to do basic statistical analysis on data files.

It is assumed that the data are given by white-space separated columns of
numbers. Lines starting with a '#' character are treated as comments. By default
the mean over the columns and the corresponding standard error of the mean are
computed.

If blocking is switched on the blocking mean and standard error are also
computed. Moreover a blocking and a 'sliding window' estimate of the integrated
autocorrelation time and are calculated.
"""


import glob
from optparse import OptionParser
from timeit import default_timer as timer
import numpy as np
import os.path



#-------------------------------------------------------------------------------
#Define command line options
#-------------------------------------------------------------------------------

parser = OptionParser("data_analyser.py [options]")

parser.add_option("-d", "--dir", dest="odir", default="./",
                  help="Working directory. Default is './'")

parser.add_option("-f", "--file", dest="file", default="",
                  help="Name of input file. Empty by default")

parser.add_option("-s", "--suffix", dest="suffix", default="dat",
                  help="Suffix of data files. Default is 'dat'.")

parser.add_option("-b", "--blocking", action="store_true", dest="blocking",
                  default=False,
                  help="Switch on blocking for data. Default is off")

parser.add_option("-D", "--debug", action="store_true", dest="debug",
                  default=False, help="Switch for debug mode. Default is off")

parser.add_option("-u", "--update-only", action="store_true", dest="update",
                  default=False,
                  help="Switch for update mode. Default is off."
                  "If set the output file is only updated if it "
                  "is older than the input file.")


parser.add_option("-L", "--skip-lines", dest="skip_lines", default=0,
                  help="Skip the lines in the input file up to the line "
                  "number given. Note that commented/empty lines are not "
                  "counted. By default no lines are skipped.")


parser.add_option("-K", "--discard-lines", dest="disc_lines", default=0,
                  help="Descard the lines in the input file after the line "
                  "number given. Note that commented/empty lines are not "
                  "counted. By default no lines are discarded.")


#-------------------------------------------------------------------------------
# Helper functions
#-------------------------------------------------------------------------------

def blcking(var, max_iter=4):
    """ Blocking analysis of the data in var.

        The array var is repeatedly blocked with a block size of two. The mean
        and its standard error are computed for the blocked arrays.

        Parameters
        ---------- 
        var : array_type 
            1D-Array of data. Any type that works with numpy.array()
        max_iter : int
            Maximal number of blocking iterations. Default is 4.

        Returns
        -------
        means : numpy array
            N-th entry is the mean for blocking with block size 2^N. 
            (Starting with N=0)
        errs : numpy array 
            Numpy array containing the standard errors for the entries in means.
    """
    means = []
    errs = []
    vvar = np.array(var)
    for dummy in range(0, max_iter):
        s = vvar.size
        if s < 2:
            break
        mm = np.mean(vvar)
        ee = np.std(vvar, ddof=1)/np.sqrt(s)
        means.append(mm)
        errs.append(ee)
        #print "%10d\t%6.4f\t%6.4e" % (s,mm,ee)
        vvar = block(vvar)

    #print means, errs
    return means, errs

def block(x_array):
    """ Blocking with block size 2

        For an array x_array of length L an array y of size ceil(L/2) is
        returned. The entries of y are the means of two consecutive values of x.

        Parameters
        ---------- 
        x_array : array_type 
            1D-Array of data. Any data type that works with numpy.array()

        Returns
        -------  
        res : numpy array 
            1D-Array. Entry res[i] = 0.5*(x[2*i]+x[2*i+1]). If the number of
            entries in x is odd res[-1] = x[-1].
    """
    xx = np.array(x_array)
    s = xx.size
    l = xx[0::2]
    r = xx[1::2]
    bs = int(np.ceil(s/2.))
    res = np.zeros(bs)
    for i in range(0, s/2):
        res[i] = 0.5*(l[i]+r[i])

        if bs > s/2:
            res[bs-1] = xx[-1]

    return res

def check_mono(x_array):
    """ Check the monotony of the elements in the array x_array.
        
        Parameters
        ---------- 
        x_array : array_type 
            1D-Array of data. Any data type that works with numpy.array()

        Returns
        ------- 
        idx : int 
            The index up to which the size of consecutive entries in x show
            (strinct) monotonic growth.
    """
    xx = np.array(x_array)
    s = xx.size
    idx = 0

    for i in xrange(1, s):
        if xx[idx] >= xx[i]:
            return idx
        idx = idx+1

    return s-1


def auto_cor(var, max_iter=4):
    """ Computes an estimate for the autocorrelation time of the entries in var.

        This function repeatedly applies blocking with block size 2 to the array
        var. Every time the mean and standard error of the blocked array is
        computed. The idea is that the blocked standard error increases with
        increasing block size until the entries of the blocked array are
        uncorrelated. The autocorrelation time is estimated by comparing the
        standard error for block size 1 with the largest standard error
        encountered during blocking.

        Parameters
        ---------- 
        var : array_type 
            1D-Array of data. Any type that works with numpy.array()
        max_iter : int 
            Maximal number of blocking iterations. Default is 4.

        Returns
        ------- 
        ac : float  
            Estimate of the autocorrelation time
        mean : float  
            Blocking estimate for the mean of var
        error : float 
            Standard error of mean.
    """
    means, errors = blcking(var, max_iter)

    idx = check_mono(errors)

    ac = 0
    if 0 != errors[idx]:
        ac = (errors[idx]/errors[0])**2*0.5

    return means[idx], errors[idx], ac

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == "__main__":

    #---------------------------------------------------------------------------
    # Parse options
    #---------------------------------------------------------------------------

    (options, args) = parser.parse_args()

    my_dir = options.odir.rstrip("/")
    my_file = options.file.rstrip()
    my_suffix = options.suffix.rstrip()
    BLOCK_FLAG = options.blocking
    DEBUG_FLAG = options.debug
    SKIP_LINES = int(options.skip_lines)
    DISC_LINES = int(options.disc_lines)
    UPDATE_FLAG = options.update


    if DISC_LINES > 0 and SKIP_LINES > 0 and DISC_LINES <= SKIP_LINES:
        print "You are skipping everything up to line {:d} and discarding "\
              "everything after line {:d}.\nNothing to do here!".format(
                  SKIP_LINES, DISC_LINES)
        exit()

    #---------------------------------------------------------------------------
    # Main Loop
    #---------------------------------------------------------------------------

    if DEBUG_FLAG:
        print "Debug Mode On"

    data_files = []
    if not my_file:
        data_files = glob.glob("*."+str(my_suffix))
    else:
        data_files.append(my_file)



    START = timer()
    print "\nWorking on ..."
    for dfile in data_files:
        print "\t->   ", dfile
        basename, extension = os.path.splitext(dfile)

        #-----------------------------------------------------------------------
        # Open file and collect data
        #-----------------------------------------------------------------------

        # Extension of output file
        xtn = ""
        if SKIP_LINES > 0:
            # Indicate number of skipped lines in output file
            xtn = "_S{:d}".format(SKIP_LINES)
        if DISC_LINES > 0:
            # Indicate number of discarded lines in output file
            xtn += "_D{:d}".format(DISC_LINES)

        xtn += ".stat"
        ofile = basename+xtn


        # If requested: Only update if input file is newer than output file.
        if UPDATE_FLAG and os.path.isfile(ofile):
            in_time = os.path.getctime(dfile)
            out_time = os.path.getctime(ofile)

            if out_time > in_time:
                print "Outfile '{:s}' is newer than input file".format(ofile)
                continue

        df = open(dfile, "r")
        data = []
        first = True
        frist_len = 0
        lnumber = 0
        for line in df:
            if line[0] == "#" or (not line.strip()):
                #print line
                continue
            lnumber += 1
            #print lnumber

            if lnumber <= SKIP_LINES:
                #print SKIP_LINES, lnumber
                continue

            if lnumber >= DISC_LINES and DISC_LINES > 0:
                # Discard the rest
                break

            lstr = line.strip().split()
            if first:
                first_len = len(lstr)
                first = False

            # Sanity check: We expect every line to have the same number of cols
            if first_len != len(lstr):
                print "\t\tWarning: Something seems to be wrong with the"\
                      "formating of {:s} in line {: d}".format(dfile, lnumber)
                print "\t\tExpected line colcount: {: d}\t"\
                      "This line: {: d}".format(first_len, len(lstr))
                continue

            data.append(lstr)

            # if DEBUG_FLAG:
            #     print data[-1]


        if SKIP_LINES > lnumber:
            print "You skipped all the lines in the data file!"
            exit()


        try:
            dar = np.array(data)
            dar = np.asfarray(dar, float)
        except ValueError:
            print "\t\tWarning: Something seems to be" \
                  "wrong with the formating of ", dfile
            continue


        # Go on if data is empty
        if dar.size < 1:
            print "\t\tWarning: No data in  ", dfile
            continue

        mean = dar.mean(axis=0)
        std = dar.std(axis=0, ddof=1)/np.sqrt(dar.shape[0])
        #print dar.size

        #-----------------------------------------------------------------------
        # Additional calculations if blocking is requested
        #-----------------------------------------------------------------------

        berrs = []
        bmeans = []
        bacors = []
        iacors = []

        if BLOCK_FLAG:
            for b in xrange(mean.size):
                #print dar[:,b]
                #print auto_cor(dar[:,b])
                m, e, a = auto_cor(dar[:, b], 12)
                berrs.append(e)
                bmeans.append(m)
                bacors.append(a)

                # Compute an estimate for the integrated autocorrelation time
                # using the sliding window ansatz of:
                #
                #    Madras, N. & Sokal, A.D.
                #    J Stat Phys (1988)50:109
                #    https://doi.org/10.1007/BF01022990
                #    (Appendix C)
                #
                #(with FFT to calculate the convolution)


                x = np.asarray(dar[:, b])
                x = x - x.mean()
                N = len(x)
                # Zero padding to the next power of two larger than 2N-1
                # to speed up FFT and avoid problems due to periodic convolution.
                # See, e.g.,  https://dsp.stackexchange.com/a/2187
                # (Answer by user 'Jean-louis Durrieu')
                # for a brief explanation.
                Npad = int(2**np.ceil(max(np.log2(N*2-1), 1)))

                # Compute convolution via FFT
                x_fft = np.fft.fft(x, Npad)
                result = np.real(np.fft.ifft(x_fft * np.conjugate(x_fft), Npad))
                # Truncate inverse FFT to desired output
                result = result[:N]

                # Do not rescale by 0
                if not 0 == result[0]:
                    result /= result[0]

                SUCCESS_tau = False
                for M in xrange(1, N):
                    tau = 0.5 + np.sum(result[1:M])
                    if M > 5.*tau:
                        SUCCESS_tau = True
                        iacors.append(tau)
                        if DEBUG_FLAG:
                            print "{: 6d}\t{: 4.2e}\t{: 4.2e}".format(M,
                                                                      result[M],
                                                                      tau)
                        break

                # Something has gone terribly wrong ...
                if not SUCCESS_tau:
                    # Set tau negative to signal problems
                    iacors.append(-1.)
                    print "Warning: Could not find an estimate for" \
                          "tau for column {:d}".format(b)

                #print tau, a
                #print iacors


        #-----------------------------------------------------------------------
        # Write output file
        #-----------------------------------------------------------------------
        header = "#Lines analysed: {:d}\n".format(dar.shape[0])
        header += "#Col\t{:^12s}\t{:^12s}".format("Mean", "StdError")
        if BLOCK_FLAG:
            header += "\t{:^12s}\t{:^12s}\t{:^12s}\t{:^12s}".format(
                "BlockingMean", "BlckStdError", "BlockingAC", "EstimatedAC")
        header += "\n"

        #print ofile

        of = open(ofile, "w")
        of.write(header)
        for m in xrange(mean.size):
            of.write("{:03d}\t{: 8.6e}\t{: 8.6e}".format(m, mean[m], std[m]))
            if BLOCK_FLAG:
                of.write("\t{: 8.6e}\t{: 8.6e}\t{: 8.6e}"
                         "\t{: 8.6e}".format(bmeans[m],
                                             berrs[m],
                                             bacors[m],
                                             iacors[m]))
            of.write("\n")

        of.close()
    END = timer()
    print "This took {: 8.6e} seconds".format(END-START)





### Local Variables:
### fill-column: 80
### eval: (auto-fill-mode)
### End:
