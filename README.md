# LOMR
Life Of a Meandering River
This repository was created by Tobias Hasse (tobiack@udel.edu) to host source code related to my thesis research.

This code is a Fork of Jon Schwenk's Life Of a Meander Bend (LOMB) source code published in conjunction with his fascinating article: 
Schwenk, Jon, Stefano Lanzoni, and Efi Foufoula‐Georgiou. "The life of a meander bend: Connecting shape and dynamics via analysis of a numerical model." Journal of Geophysical Research: Earth Surface 120.4 (2015): 690-710. 
- Open access article is available here: https://doi.org/10.1002/2014JF003252 
- Readme and code repositories are available as part of the supplementary information here:
- https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2014JF003252&file=jgrf20377-sup-0001-Readme.rtf
- https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2014JF003252&file=jgrf20377-sup-0002-CodesS1.zip 

The code in this repository is required for (my) Tobias Hasse's dissertation: Hasse, Tobias Raphael. "Storage Time Dynamics of Meandering River Floodplain Sediments: A Modeling Study." PhD diss., University of Delaware, 2021.

The dissertation includes the code as appendices and that code runs to completion but is dependant on two MATLAB toolboxes:
- The Signal Processin Toolbox, and
- The Statistics and Machine Learning Toolbox.

# This repository contains additional code to eliminate those dependencies:


# Other changes from the original LOMB to this repository LOMR include
- Adjustments to algorithms here and there to boost computational performance.
- Writing portions of the code in c and compiling it for MATLAB, particularly the code which calculates the hydrodynamic flowfield of the meandering river
- Truncating the algorithm to calculate the flowfield


This repo is specifically set up to be called by functions in my repositories MRDAST & MRLATS and will generate output files which are part of the analysis in the repositories MRDAST & MRVIST, and together will duplicate the bulk of the dissertation work referenced above.

A package of code files will be added to this repository at a later time
