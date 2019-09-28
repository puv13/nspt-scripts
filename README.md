# nspt-scripts
A collection of Python scripts to analyse data from Numerical Stochastic Perturbation Theory (NSPT) simulations.

## Typical use case
In the current state the scripts are meant to process data from NSPT simulations of the
Principal Chiral Model. A generalisation to any other model or theory should be 
straightforward.

The scripts will typically be called in the following order:

1. **data_analyser.py**: Performs a statistical analysis of the NSPT data and computes the
mean of the observables over the Langevin trajectory.

2. **pcm_analyser.py**: Collects the data for different lattice volumes and stochastic
time steps. Performs the extrapolation to vanishing time step size at fixed lattice
volumes if possible.

3. **make_eps0.00_all_v.py**: Collects the data for vanishing time step and all available
lattice volumes in one file.

4. **fit.py**: Performs an infinite volume extrapolation by fitting to the known
functional form of the finite volume effects and writes the result to a file.

5. **pcm_ratios.py**: Generates gnuplot compatible files to plot ratios of expansion
coefficients.

6. **renormalon_fit.py**: Asymptotically the ratios of expansion coefficients should not
depend on the rank N. This script performs a fit to a constant for ratios at fixed
expansion order as a function of the rank N.


### Code for extrapolation to infinite lattice size

- **fit_function_gen.py**: Code to generate the fit functions for the infinite volume
    extrapolation. Relies on SymPy to perform the series expansions.

- **simultaneous_fit.py**: Code to perform a simultaneous fit to several data sets. Fit
    functions for different data sets need not be the same, but the fit parameters have to
    be.


## Contributors 
The source code and data files in this repository originated from a research project. The
[contributors file](Contributors.md) list all persons involved in that project,
not just the authors of the source code.

## Licensing 
Unless otherwise stated all the files in this repository are licensed under the [MIT
license](License.md).  An exception are the data files in the *data/PCM* sub-folder, which
are licensed under a [Creative Commons ](data/PCM/LICENSE.md) license.


