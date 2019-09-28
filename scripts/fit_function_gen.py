#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This module defines functions to generate the fit function for the volume
    dependence of the NSPT calculations of the PCM model.
"""


#-------------------------------------------------------------------------------
# Testing sympy expansions and stuff
#-------------------------------------------------------------------------------
import numpy as np
import sympy as sp
import pcm_coeff as pcm
import pickle
import re

#-------------------------------------------------------------------------------
# Housekeeping
#-------------------------------------------------------------------------------
# Switch debugging output of function on/off
DEBUG_FLAG = 0


#-------------------------------------------------------------------------------
# Globally used Symbols and Functions
#-------------------------------------------------------------------------------
x, y, f, c, b0, b1, N = sp.symbols('x y f c b0 b1 N', real=True)
beta = sp.Function('beta')
alpha = sp.Function('alpha')
# Beta function of the theory
beta = -(b0*alpha(x)**2+b1*alpha(x)**3)

def dalpha(n):
    """ Returns the n-th derivative of the coupling alpha

        Parameters
        ----------
        n : int
            Order of the derivative

        Returns
        -------
        dalpha : sympy_expression
            The n-th derivative of the coupling alpha.

        Notes
        -----
        This function uses a sloppy way to get a kind of static variable in
        python, see, e.g., https://stackoverflow.com/questions/279561 and the
        answer https://stackoverflow.com/a/279586 by user 'Claudiu', in order to
        reuse the computed derivatives.
        Note that 'dalpha.pre_computed' is accesible outside of the function
        dalpha(n) and in general this is a very bad idea and should not be
        done. I warned you!
    """

    # Check if the desired derivative is known  ...
    try:
        return dalpha.pre_computed[n]

    # ... and compute it if necessary.
    except IndexError:

        while len(dalpha.pre_computed) <= n:
            dalpha.pre_computed.append(sp.simplify(
                sp.diff(dalpha.pre_computed[-1], x).subs(sp.diff(alpha(x), x),
                                                         beta)))

        return dalpha.pre_computed[n]
dalpha.pre_computed = [alpha(x)]




def rec_power(expr, n):
    """ Helper function to recursively compute powers
        of the sympy expression expr.

        In this function some care is taken to keep the size of the sympy
        expressions involving expansions (with 'order terms' sympy.O) as small
        as possible.  This is achieved by expanding intermediate results and
        getting rid of high order terms that would not contribute to the final
        result anyway.


        Parameters
        ----------
        expr : sympy_expression
            An arbitrary sympy expression to be raised to the power n
        n : int

        Returns
        -------
        expr**n : sympy_expression


        See Also
        --------
        - sympy.O
        - sympy.expand
    """

    # Sanitise input
    #---------------------------------------------------------------------------
    if not int(n) == n:
        print "\t Warning: Non-integer n was cast to "\
              "integer in function 'rec_power'"
        n = int(n)

    # Take care of negative powers
    if n < 0:
        expr = 1./expr
        n = -n

    # Evaluate power
    #---------------------------------------------------------------------------
    if 0 == n:
        return 1
    elif 1 == n:
        return expr
    # even powers
    elif 0 == n%2:
        aux = rec_power(expr, n/2)
        try:
            res = aux*aux
            res = res.expand()
            return res
        except AttributeError:
            return aux*aux
    # specialisation for powers that are multiples of  3
    elif 0 == n%3:
        aux = rec_power(expr, n/3)
        try:
            res = (aux*aux).expand()*aux
            res = res.expand()
            return res
        except AttributeError:
            return aux*aux*aux
    # general case
    else:
        aux = rec_power(expr, n-1)
        try:
            aux = aux.expand()
        except:
            pass

        return aux*expr


def expansion_slow(maxorder=4):
    """ This is the naive (and very slow) implementation of the
        ("hard mode" part of the) expansion.

        This function is intended to be a debugging tool.
        Do not use this in production code.

        Parameters
        ----------
        maxorder : int
            Highest order of \alpha taken into account in the expansion.

        See also
        --------
        - expansion
    """
    res = 0
    f = sp.IndexedBase('f')

    kres = 0
    for k in range(maxorder):
        kres += dalpha(k)/sp.factorial(k)*((-sp.log(N))**k)

    for i in xrange(maxorder):
        res += sp.Indexed(f, i)*(kres**(i+1))

    return res


def expansion(maxorder=4):
    """
        A function to compute the expansion of finite volume
        correction parameters f_i(N) in terms of f_i, c_i,  N, and the beta
        function of the theory.

        Parameters
        ----------
        maxorder : int
            Highest order of \alpha taken into account in the expansion.

        References
        ----------
        For details see, e.g.,
        'Bali et al.,  PRD  89, 054505 (2014)'
    """
    res = 0
    #_ii = sp.symbols("_ii")
    f = sp.IndexedBase('f')
    kres = 0

    # "hard mode" expansion (coefficients f[i] )
    #---------------------------------------------------------------------------
    for k in range(maxorder):
        kres += dalpha(k)/sp.factorial(k)*((-sp.log(N))**k)

    # Use a simpler symbol for alpha once we have the expansion
    kres = kres.subs(alpha(x), y).expand()
    #kres = kres.removeO()

    for i in xrange(maxorder):

        order_cutoff = maxorder+1
        if i > 0:
            order_cutoff = maxorder + 1 - i

        if DEBUG_FLAG:
            print "\tterm:{:02d}/{:02d}".format(i+1, maxorder)

        # Remove powers that, when kres is raised to power i, will be
        # of order larger than maxorder
        kres_to_order = kres + sp.O(y**order_cutoff)

        # Make aux an expansion up to order (maxorder)
        aux = kres_to_order.removeO()
        # The previous line is necessary!
        # (spympy.O(x**n)+sympy.O(x**m) =  spympy.O(x**min(n, m)))
        aux = aux + sp.O(y**(maxorder+1))

        aux = rec_power(aux, (i+1))
        # Make aux an expansion up to order (maxorder) again
        # spympy.O(x**n)*sympy.O(x**m) = sympy.O(x**(n+m))
        aux = aux + sp.O(y**(maxorder+1))
        aux = aux.expand(y)

        res += sp.Indexed(f, i)*(aux)

    res = res.expand()

    # "soft mode" expansion (coefficients c[i])
    #---------------------------------------------------------------------------
    hard_exp = 1
    for i in xrange(maxorder-1):
        hard_exp += sp.Indexed(c, i)*y**(i+1)

    res = res*hard_exp
    res = res.expand()

    return res


expansion_coefficient_list = []
expansion_coeff_brine = "_FFG_coefficient_list_precomp"
def expansion_coefficient(n, precom=True):
    """ Compute the n-th coefficient for the expansion in alpha.

        Note that the first entry, i.e. expansion_coefficinet(0), is the
        factor in front of alpha**1 and not the constant term.

        In general the term alpha**m comes with the coefficient
        expansion_coefficient(m+1).

        Parameters
        ----------
        n : int
            Coefficient index
        precom : bool
            If precom is True coefficients are read from a file if available
            and newly computed coefficients are written to a file.
            Default is True.

        Returns
        -------
        The n-th coefficient for the expansion in alpha. See above for the
        index convention.
    """
    global expansion_coefficient_list
    # Sanitise input
    if n < 0 or int(n) != n:
        print "\tWarning: Parameter error in 'expansion_coefficient'."
        print "\t         Parameter has to be a non-negative integer. "
        return 1.


    # Check if the desired coefficient is known  ...
    try:
        return expansion_coefficient_list[n]

    # ... and compute it if necessary.
    except IndexError:

        # Check if coeff is pre-computed
        if precom:
            pcl_list = []
            try:
                pcl_list = pickle.load(open(expansion_coeff_brine, "rb"))
            except:
                pass


            if len(pcl_list) > n:
                expansion_coefficient_list = pcl_list
                return expansion_coefficient_list[n]


        res = expansion(n+1).expand()

        # Fill lookup table
        for j in range(len(expansion_coefficient_list), n+1):
            expansion_coefficient_list.append(res.coeff(y, j+1).collect(
                sp.log(N).removeO()))

        # Check if pre-computed list should be pickled
        if precom:
            pickle.dump(expansion_coefficient_list,
                        open(expansion_coeff_brine, "wb"))

        return expansion_coefficient_list[-1]


epc_strings = {}
def epc_string(n, bb1=False, c_coeffs=False, ToFloat=True):
    """ Return a string of expansion coefficient n for use in output *.dat
        files.

        Numerical expressions are evaluated as far as possible and N
        is replaced by x (for plotting in gnuplot). Moreover f[j] is replaced
        by f_b0_j or f_b1_j, depending on the value of bb1.
        Results are cached for later use.

        Parameters
        ----------
        n : int
            Desired coefficient order. Keep in mind the index convention of
            the expansion_coefficient function!

        bb1 : bool
            Flag to toggle inclusion/exclusion of beta_1 terms.
            Default is False.

        c_coeffs : bool
            Flag to toggle inclusion/exclusion of c_i coefficients.
            Default is False.

        ToFloat : bool
            If True numerical expressions (like ratios) are evaluates as far
            as possible. Default is True.

        Returns
        -------
        The expansion coefficient of order n converted to a string.

    """

    try:

        ep_str = epc_strings[(n, bb1, c_coeffs, ToFloat)]

        return ep_str

    except KeyError:

        ep_n = expansion_coefficient(n)
        ep_n = ep_n.subs(N, x)

        if ToFloat:
            ep_n = sp.N(ep_n)
            #print ep_n

        if not bb1:
            ep_n = ep_n.subs(b1, 0)

        if not c_coeffs:
            coeff_c_subs_dict = {sp.Indexed(c, i):0 for i in range(n)}
            ep_n = ep_n.subs(coeff_c_subs_dict)


        ep_n = str(ep_n)


        # Update global dict
        if bb1:
            ep_n = re.sub(r'f\[(\d*)\]', r'f_b1_\1', ep_n)
            ep_n = re.sub(r'c\[(\d*)\]', r'c_b1_\1', ep_n)

        else:
            ep_n = re.sub(r'f\[(\d*)\]', r'f_b0_\1', ep_n)
            ep_n = re.sub(r'c\[(\d*)\]', r'c_b0_\1', ep_n)

        epc_strings[(n, bb1, c_coeffs, ToFloat)] = ep_n


        return ep_n




logsum_b0_dict = {}
def fit_func_b0(xx, NN, i, ORD, bb0, known, *args):
    """ 
        Fit function describing the finite volume effects for a given order for a
        symmetric lattice with V=L*L. 

        This function is specialised for the case where the c_i coefficients and
        the higher order \beta-function coefficients \beta_j with j>0 are set to
        zero. 

        WARNING!!!
        ----------
        Note that the coeffs in the PCM calculation are indexed differently
        (c_n corresponds to order \alpha^n , not \alpha^(n+1)


        Parameters
        ---------- 
        xx : float 
            The linear lattice extend L (V=L*L). 
        NN : int 
            The rank N of the SU(N) matrices of the Principal Chiral model. 
        i : int 
            The expansion order. 
        ORD : int 
            The highest expansion order considered in the calculation. This is
            has to be known to calculate the total number of free parameters in 
            the fits. 
        bb0 : float 
            The value of the \beta_0 coefficient of the \beta-function of the
            Principal Chiral model. 
        known : int  
            The number of known infinite volume expansion coefficients used as
            input for the fits.  
        args : array_type 
            An array that holds all the fit coefficients

        Returns
        ------- 
        The fit function for given order i, rank N and lattice size xx. 
          
    """ 

    # Check if logsum already computed ...
    try:
        logsum = logsum_b0_dict[i]

    # ... and compute if unknown
    except KeyError:
        #print expansion_coefficient(i)
        # Symbols for lambdify
        lgN, f = sp.symbols("lgN, f")
        coeff_c_subs_dict = {sp.Indexed(c, i):0 for i in range(ORD+1)}
        # Correct for different indexing in expansion_coefficient
        log_exp = expansion_coefficient(i-1)
        log_exp = log_exp.subs({b1:0, sp.log(N):lgN}).subs(coeff_c_subs_dict)
        #print log_exp
        #Generate logsum function of given order i
        logsum = sp.lambdify((b0, lgN, f), log_exp, modules="numpy")

        # Update the global dict
        logsum_b0_dict[i] = logsum

    numfitcoef = int(np.floor(ORD/2.)) +1 - known

    xx = np.array(xx)
    xx_log = np.log(xx)


    if i < known:
        return pcm.coef(2*i, NN) + logsum(bb0, xx_log, args[numfitcoef:])/xx**2
    else:
        return args[int(i)-known] + logsum(bb0, xx_log, args[numfitcoef:])/xx**2

logsum_b1_dict = {}
def fit_func_b1(xx, NN, i, ORD, bb0, bb1, known, *args):
    """ 
        Fit function describing the finite volume effects for a given order for a
        symmetric lattice with V=L*L. 

        This function is specialised for the case where the c_i coefficients are
        set to zero and the \beta-function coefficients \beta_0 and \beta_1 are
        set to their physical values. All \beta_j with j>1 are set to zero. 
        
        WARNING!!!
        ----------
        Note that the coeffs in the PCM calculation are indexed differently
        (c_n corresponds to order \alpha^n , not \alpha^(n+1)


        Parameters
        ---------- 
        xx : float 
            The linear lattice extend L (V=L*L). 
        NN : int 
            The rank N of the SU(N) matrices of the Principal Chiral model. 
        i : int 
            The expansion order. 
        ORD : int 
            The highest expansion order considered in the calculation. This is
            has to be known to calculate the total number of free parameters in 
            the fits. 
        bb0 : float 
            The value of the \beta_0 coefficient of the \beta-function of the
            Principal Chiral model. 
        bb1 : float 
            The value of the \beta_1 coefficient of the \beta-function of the
            Principal Chiral model. 
        known : int  
            The number of known infinite volume expansion coefficients used as
            input for the fits.  
        args : array_type 
            An array that holds all the fit coefficients

        Returns
        ------- 
        The fit function for given order i, rank N and lattice size xx.           
    """ 
    # Check if logsum already computed ...
    try:
        logsum = logsum_b1_dict[i]

    # ... and compute if unknown
    except KeyError:
        #print expansion_coefficient(i)
        # Symbols for lambdify
        lgN, f = sp.symbols("lgN, f")
        coeff_c_subs_dict = {sp.Indexed(c, i):0 for i in range(ORD+1)}
        # Correct for different indexing in expansion_coefficient
        log_exp = expansion_coefficient(i-1)
        log_exp = log_exp.subs({sp.log(N):lgN}).subs(coeff_c_subs_dict)
        #Generate logsum function of given order i
        logsum = sp.lambdify((b0, b1, lgN, f), log_exp, modules="numpy")
        logsum_b1_dict[i] = logsum

    numfitcoef = int(np.floor(ORD/2.)) +1 - known

    xx = np.array(xx)
    xx_log = np.log(xx)


    if i < known:
        return pcm.coef(2*i, NN) + logsum(bb0, bb1, xx_log,
                                          args[numfitcoef:])/xx**2
    else:
        return args[int(i)-known] + logsum(bb0, bb1, xx_log,
                                           args[numfitcoef:])/xx**2


# Note that the coeffs in the PCM calculation are indexed differently
# (c_n corresponds to order \alpha^n , not \alpha^(n+1)
logsum_b0_c_dict = {}
def fit_func_b0_c(xx, NN, i, ORD, bb0, known, *args):
    """ 
        Fit function describing the finite volume effects for a given order for a
        symmetric lattice with V=L*L. 

        This function is specialised for the case where the c_i coefficients are
        fit parameters and the higher order \beta-function coefficients \beta_j
        with j>0 are set to zero. 

        WARNING!!!
        ----------
        Note that the coeffs in the PCM calculation are indexed differently
        (c_n corresponds to order \alpha^n , not \alpha^(n+1)


        Parameters
        ---------- 
        xx : float 
            The linear lattice extend L (V=L*L). 
        NN : int 
            The rank N of the SU(N) matrices of the Principal Chiral model. 
        i : int 
            The expansion order. 
        ORD : int 
            The highest expansion order considered in the calculation. This is
            has to be known to calculate the total number of free parameters in 
            the fits. 
        bb0 : float 
            The value of the \beta_0 coefficient of the \beta-function of the
            Principal Chiral model. 
        known : int  
            The number of known infinite volume expansion coefficients used as
            input for the fits.  
        args : array_type 
            An array that holds all the fit coefficients

        Returns
        ------- 
        The fit function for given order i, rank N and lattice size xx. 
          
    """ 

    # Check if logsum already computed ...
    try:
        logsum = logsum_b0_c_dict[(i, ORD)]

    # ... and compute if unknown
    except KeyError:
        #print expansion_coefficient(i)
        # Symbols for lambdify
        global c
        lgN, f = sp.symbols("lgN, f")
        # c[ORD/2-2] and f[ORD/2-1] are not distinguishable in fits
        log_exp = expansion_coefficient(i-1)
        log_exp = log_exp.subs({b1:0, sp.log(N):lgN,
                                sp.Indexed(c, int(np.floor(ORD/2.))-2):0})
        #print log_exp
        #Generate logsum function of given order i
        logsum = sp.lambdify((b0, lgN, f, c), log_exp, modules="numpy")

        # Update the global dict
        logsum_b0_c_dict[(i, ORD)] = logsum

    numfitcoef = int(np.floor(ORD/2.)) + 1 - known
    idx_c_coef = int(np.floor(ORD/2.)) + numfitcoef

    xx = np.array(xx)
    xx_log = np.log(xx)


    if i < known:
        return pcm.coef(2*i, NN) + logsum(bb0, xx_log,
                                          args[numfitcoef:idx_c_coef],
                                          args[idx_c_coef:])/xx**2
    else:
        return args[int(i)-known] + logsum(bb0, xx_log,
                                           args[numfitcoef:idx_c_coef],
                                           args[idx_c_coef:])/xx**2


# Note that the coeffs in the PCM calculation are indexed differently
# (c_n corresponds to order \alpha^n , not \alpha^(n+1)
logsum_b1_c_dict = {}
def fit_func_b1_c(xx, NN, i, ORD, bb0, bb1, known, *args):
    """ 
        Fit function describing the finite volume effects for a given order for a
        symmetric lattice with V=L*L. 

        This function is specialised for the case where the c_i coefficients are
        fit parameters and the \beta-function coefficients \beta_0 and \beta_1 are
        set to their physical values. All \beta_j with j>1 are set to zero. 

        WARNING!!!
        ----------
        Note that the coeffs in the PCM calculation are indexed differently
        (c_n corresponds to order \alpha^n , not \alpha^(n+1)

        
        Parameters
        ---------- 
        xx : float 
            The linear lattice extend L (V=L*L). 
        NN : int 
            The rank N of the SU(N) matrices of the Principal Chiral model. 
        i : int 
            The expansion order. 
        ORD : int 
            The highest expansion order considered in the calculation. This is
            has to be known to calculate the total number of free parameters in 
            the fits. 
        bb0 : float 
            The value of the \beta_0 coefficient of the \beta-function of the
            Principal Chiral model. 
        bb1 : float 
            The value of the \beta_1 coefficient of the \beta-function of the
            Principal Chiral model. 
        known : int  
            The number of known infinite volume expansion coefficients used as
            input for the fits.  
        args : array_type 
            An array that holds all the fit coefficients

        Returns
        ------- 
        The fit function for given order i, rank N and lattice size xx. 
          
    """ 
    # Check if logsum already computed ...
    try:
        logsum = logsum_b1_c_dict[(i, ORD)]

    # ... and compute if unknown
    except KeyError:
        #print expansion_coefficient(i)
        # Symbols for lambdify
        lgN, f = sp.symbols("lgN, f")
        global c
        # c[ORD/2-2] and f[ORD/2-1] are not distinguishable in fits
        log_exp = expansion_coefficient(i-1)
        log_exp = log_exp.subs({sp.log(N):lgN,
                                sp.Indexed(c, int(np.floor(ORD/2.))-2):0})
        #Generate logsum function of given order i
        logsum = sp.lambdify((b0, b1, lgN, f, c), log_exp, modules="numpy")

        # Update the global dict
        logsum_b1_c_dict[(i, ORD)] = logsum

    numfitcoef = int(np.floor(ORD/2.)) + 1 - known
    idx_c_coef = int(np.floor(ORD/2.)) + numfitcoef

    xx = np.array(xx)
    xx_log = np.log(xx)


    if i < known:
        return pcm.coef(2*i, NN) + logsum(bb0, bb1, xx_log,
                                          args[numfitcoef:idx_c_coef],
                                          args[idx_c_coef:])/xx**2
    else:
        return args[int(i)-known] + logsum(bb0, bb1, xx_log,
                                           args[numfitcoef:idx_c_coef],
                                           args[idx_c_coef:])/xx**2



# def print_Fs(order):
#     """ Print the expansion coefficients F_i
#     """
#     pr_exp = sp.expand(expansion(order))
#     pr_exp = sp.collect(pr_exp, y)

#     #print pr_exp
#     #print pr_exp.coeff(y, 2)µ
#     #print sp.collect(pr_exp.coeff(y, 2), sp.log(N))
#     for k in range(order):
#         print "f_{:d}: \t ".format(k), \
#               (sp.collect(pr_exp.coeff(y, k+1), sp.log(N))), "\n"



if __name__ == "__main__":
    import time
    def debug_exp(n):
        """ Debug Stuff """
        res2 = 0
        start = time.clock()
        res2 = expansion(n).removeO()
        res2 = res2.expand()
        end = time.clock()
        print end-start
        return res2


    def print_it(res2, m):
        """ More Debug Stuff """
        print res2.coeff(y, m).coeff(sp.log(N), m-3).subs(b1, 0),\
              res2.coeff(y, m).coeff(sp.log(N), 2).subs(b1, 0), (m-1)*(m-2)/2

    def ese(n):
        """ Even More Debug Stuff """
        n = n-1
        rese = sp.Sum((b0*sp.log(N))**y*sp.Indexed('f', n-y)*sp.binomial(n, y),
                      (y, 0, n)).doit()
        return rese


    def check(n):
        """
           A simple consistency check for the expansion
        """

        coeff_c_subs_dict = {sp.Indexed(c, i):0 for i in range(n+1)}
        dd = debug_exp(n).subs({b1:0}).subs(coeff_c_subs_dict).removeO()


        for i in range(0, n):
            j = i+1
            #print ese(j)
            #print dd.coeff(y, j)
            print "\torder {:2d} -> {:s}".format(j,
                                                 str(ese(j) == dd.coeff(y, j)))


    print "done ..."
    check(4)
    #fit_func_b0(1, 3, 4, 4, 2, 0, [1, 2, 3, 4, 5, 6])

    #exit(0)

### Local Variables:
### mode: python
### fill-column: 80
### eval: (auto-fill-mode)
### End:
