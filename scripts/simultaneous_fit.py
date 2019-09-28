#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This module defines a function for simultaneous fits to several
    data sets.

    The fit function can be the same for all data sets or a different
    function for every data set. The important point is that all the
    functions have to depend on the same set of fit parameters.
"""


import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def sim_fit(F, x, Y, Y_err=None, init_params=2, **kwargs):
    """
        Perform a simultaneous fit to several data sets.

        Parameters
        ----------

        F : function or list of functions
            A function or a list of functions to be fitted
        x : array_type
            Values for the independent fit variable
        Y : array_type
            List of data sets to be fitted. The number of data points
            has to be an integer multiple of the number of x values.
        Y_err : array_type
            Errors for the data points. This parameter is obtional.
            Default value is 'None'.
        init_params : int
            Number of fit parameters. Default is 2.
        kwargs : dict
            Key-word arguments that are passed on to the scipy curve_fit
            function.

        Returns
        -------
        popt : array_type
            The optimal fit parameters.
        pcov : array_type
            The covariance matrix for the optimal fit parameters.
       chi2dof : float
            Chi^2 per degree of freedom.

        See Also
        --------
        - scipy.optimize.curve_fit

    """

    # --------------------------------------------------------------------------
    # Convert input to (1D) numpy arrays
    # --------------------------------------------------------------------------
    x = np.array(x)
    Y = np.array(Y).flatten()

    # --------------------------------------------------------------------------
    #Sanity checks
    # --------------------------------------------------------------------------
    xl = len(x)
    Yl = len(Y)

    if Y_err is not None:
        Y_err = np.array(Y_err).flatten()
        Yel = len(Y_err)

        # Todo: Proper error handling
        if Yl != Yel:
            print "\tError: Y and Y_err have different lenght."
            return

    if Yl%xl != 0:
        print "\tError: x and Y dimensions do not match."
        return

    # --------------------------------------------------------------------------
    # Ugly check if we are dealing with one function or a list of
    # functions
    # --------------------------------------------------------------------------
    if type([]) == type(F):
        if len(F) != Yl/xl:
            print "\tError: Mismatched length of F"
            return


    #print xl, Yl

    # --------------------------------------------------------------------------
    # Larger X for simultaneous fit
    # --------------------------------------------------------------------------
    X = []

    for i in xrange(Yl/xl):
        X.append(x.copy())

    X = np.array(X).flatten()
    #print X


    # --------------------------------------------------------------------------
    # Define fit function
    # It is possible to have different fit functions for different
    # data sets. The same parameters, however, should be shared by
    # all fit functions.
    # --------------------------------------------------------------------------
    def fit_F(xx, *args):
        """
           Fit function used for the simultaneous fits.
        """
        if type([]) == type(F):
            yy = []
            for i in xrange(Yl/xl):
                yy.append(F[i](xx[xl*i:xl*(i+1)], *args))

            return np.array(yy).flatten()

        else:
            return F(xx, *args)

    # --------------------------------------------------------------------------
    # Fit
    # --------------------------------------------------------------------------

    # Determine number of params
    try:
        dummy = len(init_params)
        init_params = tuple(init_params)

    except TypeError:
        init_params = tuple(np.zeros(init_params))


    # Make sure the user does not use init_params AND the 'p0' named argument
    my_args = {i:kwargs[i] for i in kwargs if i != 'p0'}
    if "p0" in kwargs.keys():
        print "\tWARNING: Named argument 'p0' not permitted for this function."\
              " Use init_params instead."





    popt, pcov = curve_fit(fit_F, X, Y, sigma=Y_err, p0=init_params, **my_args)
    chi2 = np.sum(((fit_F(X, *popt)-Y)/Y_err)**2)
    chi2dof = chi2/(len(Y)-len(popt))

    #print (len(Y)-len(popt))

    return  popt, pcov, chi2dof


if __name__ == "__main__":

    def test_func_one(x, a, b):
        """
           Returns 2*b*exp(-a*x).
        """
        return 2.*b*np.exp(-a*x)


    def test_func_two(x, a, b):
        """
           Returns b*exp(-a*x).
        """
        return b*np.exp(-a*x)

    def test_func_three(x, a, b):
        """
           Returns 0.5*b*exp(-a*x).
        """
        return 0.5*b*np.exp(-a*x)


    x = np.linspace(0, 2, 12)

    # 'True' parameter values
    a = .75
    b = 2.

    yo = test_func_one(x, a, b)
    yt = test_func_two(x, a, b)
    yd = test_func_three(x, a, b)

    #np.random.seed(0)
    #np.random.normal(0, .2, 1000000)


    # Add some noise
    rands = np.random.normal(0, .2, len(yo))
    yoe = yo*(1. + rands)
    yoee = np.abs(yoe-yo)
    rands = np.random.normal(0, .2, len(yo))
    yte = yt*(1. + rands)
    ytee = np.abs(yte-yt)
    rands = np.random.normal(0, .2, len(yo))
    yde = yd*(1. + rands)
    ydee = np.abs(yde-yd)


    # Extract original parameters from noisy data
    oopt, ocov = curve_fit(test_func_one, x, yoe, sigma=yoee,
                           absolute_sigma=True)
    oerr = np.sqrt(np.diag(ocov))
    a1 = oopt[0]; a1e = oerr[0]
    b1 = oopt[1]; b1e = oerr[1]

    topt, tcov = curve_fit(test_func_two, x, yte, sigma=ytee,
                           absolute_sigma=True)
    terr = np.sqrt(np.diag(tcov))
    a2 = topt[0]; a2e = terr[0]
    b2 = topt[1]; b2e = terr[1]

    dopt, dcov = curve_fit(test_func_three, x, yde, sigma=ydee,
                           absolute_sigma=True)
    derr = np.sqrt(np.diag(dcov))
    a3 = dopt[0]; a3e = derr[0]
    b3 = dopt[1]; b3e = derr[1]


    # Extract original parameters by a simultaneous fit to all data sets
    com_opt, com_cov, chi = sim_fit([test_func_one, test_func_two,
                                     test_func_three],
                                    x,
                                    [yoe, yte, yde],
                                    Y_err=[yoee, ytee, ydee],
                                    absolute_sigma=True)

    #print com_cov
    com_err = np.sqrt(np.diag(com_cov))

    ac = com_opt[0]; ace = com_err[0]
    bc = com_opt[1]; bce = com_err[1]

    af = (a1+a2+a3)/3.; afe = np.sqrt(a1e**2+a2e**2+a3e**2)/3.
    bf = (b1+b2+b3)/3.; bfe = np.sqrt(b1e**2+b2e**2+a3e**2)/3.

    ofit = test_func_one(x, a1, b1)
    tfit = test_func_two(x, a2, b2)
    dfit = test_func_three(x, a3, b3)

    cofit = test_func_one(x, ac, bc)
    ctfit = test_func_two(x, ac, bc)
    cdfit = test_func_three(x, ac, bc)

    xmin = -0.1
    xmax = 2.1

    ymin = 0
    ymax = 4.5

    plt.subplot(2, 1, 1)
    axes = plt.gca()
    axes.set_xlim([xmin, xmax])
    axes.set_ylim([ymin, ymax])
    plt.plot(x, yo, "k", x, ofit, "r.",
             x, tfit, "r.", x, dfit, "r.", x, yt, "k", x, yd, "k")
    plt.errorbar(x, yoe, yoee, fmt="bo")
    plt.errorbar(x, yte, ytee, fmt="go")
    plt.errorbar(x, yde, ydee, fmt="mo")
    plt.text(0.75, 4.00, r'$a1={:6.4e}\pm{:6.4e} \quad b1={:6.4e}'\
                          '\pm{:6.4e}$'.format(a1, a1e, b1, b1e))
    plt.text(0.75, 3.85, r'$a2={:6.4e}\pm{:6.4e} \quad b2={:6.4e}'\
                          '\pm{:6.4e}$'.format(a2, a2e, b2, b2e))
    plt.text(0.75, 3.70, r'$a3={:6.4e}\pm{:6.4e} \quad b3={:6.4e}'\
                          '\pm{:6.4e}$'.format(a3, a3e, b3, b3e))
    plt.text(0.75, 3.45, r'$af={:6.4e}\pm{:6.4e} \quad bf={:6.4e}'\
                          '\pm{:6.4e}$'.format(af, afe, bf, bfe))
    plt.title('Simultaneous fit test')
    plt.ylabel('y')

    plt.subplot(2, 1, 2)
    axes = plt.gca()
    axes.set_xlim([xmin, xmax])
    axes.set_ylim([ymin, ymax])
    plt.errorbar(x, yoe, yoee, fmt="bo")
    plt.errorbar(x, yte, ytee, fmt="go")
    plt.errorbar(x, yde, ydee, fmt="mo")
    plt.plot(x, cofit, "r.", x, ctfit, "r.", x, cdfit, "r.", x, yt, "k",
             x, yo, "k", x, yd, "k")
    plt.text(0.75, 3.75, r'$as={:6.4e}\pm{:6.4e} \quad bs={:6.4e}'\
                          '\pm{:6.4e}$'.format(ac, ace, bc, bce))
    plt.xlabel('x')
    plt.ylabel('y')

    plt.show()
    plt.close()



### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
