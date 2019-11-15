========================================================
	Colin F. Turley [cft114 at psu dot edu]
	This file is part of the VERITAS-IC40 analysis
	published in arxiv:1608.08983

	Copyright (c) by
	Colin F. Turley
=========================================================

Requirements: numpy, scipy, pyGPs

These scripts were developed to perform a search for neutrinos
from six northern blazars using IceCube 40 string data, and
VERITAS neutrino data.

VERITAS TeV data files are located in the data directory, as is the
neutrino data from IC40.

The first script to run is kdetest.py, which will calculate the flux
value we use for the meanfunction value in the Gaussian process
interpolation.

Second, the lightcurves are interpolated using the interpcurves.py
script. This will save the interpolation results in the
itms directory.

Times of interest can then be found for mrk 421 using
mrk421searchingprogram.py, and for the other blazars using
curvescombined.py

Next, the neutrino expectations can be calculated using
nulimitest.py. The scrambled search is carried out using searching.py
and the unscrambled search by using realsearch.py

The expected neutrino flux can then be found hsing nuflux.py

These scripts should let anyone reproduce the analysis found in the
paper. Understand that these scripts are some of my first efforts at
coding a full analysis and therefore are very distinctly suboptimal
and almost certainly not the best written compared to my later work.
