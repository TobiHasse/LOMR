# LOMR
Life Of a Meandering River
This repository was created by Tobias Hasse (tobiack@udel.edu) to host source code related to my thesis research.

### Function and use
- LOMR will generate the planform evolution and history of a meandring river, output as a MATLAB struct of river (X,Y) coordinate pairs for the evolution timeseries
- LOMR output can be used for further analysis by other repositories, MRDAST, MRLATS, MRVIST studying:
- Meandering River Dynamics and Storage Time
- Meandering River Length and Time Scales
- Meandering River Variability in Storage Time

This code is a Fork of Jon Schwenk's Life Of a Meander Bend (LOMB) source code published in conjunction with his fascinating article: 
Schwenk, Jon, Stefano Lanzoni, and Efi Foufoula‐Georgiou. "The life of a meander bend: Connecting shape and dynamics via analysis of a numerical model." Journal of Geophysical Research: Earth Surface 120.4 (2015): 690-710. 
- Open access article is available here: https://doi.org/10.1002/2014JF003252 
- Readme and code repositories are available as part of the supplementary information here:
- https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2014JF003252&file=jgrf20377-sup-0001-Readme.rtf
- https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2014JF003252&file=jgrf20377-sup-0002-CodesS1.zip 

The code in this repository is required for my dissertation: Hasse, Tobias Raphael. "Storage Time Dynamics of Meandering River Floodplain Sediments: A Modeling Study." PhD diss., University of Delaware, 2021.

The dissertation includes the code as appendices and that code runs to completion but is dependant on two MATLAB toolboxes:
- The Signal Processin Toolbox, and
- The Statistics and Machine Learning Toolbox.

### This repository contains additional code to eliminate those dependencies:
- TRH_rm_cutoffs.m is a slightly slower workaround to use without the Statistics and Machine Learning Toolbox
- sgolayfilt.m & sgoaly_octave.m are files to workaround the absence of the Signal Processing Toolbox.  If this toolbox is available, see sgolayfilt.m and edit to call sgolay.m
- additionally this repository updates calls to the function nanmean which is only available in the Statistics and Machine Learning Toolbox to call mean(_____,'omitnan') with the flag: 'omitnan' which can be used to edit other function calls e.g. nanstd() etc.

### Other changes from the original LOMB to this repository LOMR include
- Adjustments to algorithms here and there to boost computational performance.
- Writing portions of the code in c and compiling it for MATLAB, particularly the code which calculates the hydrodynamic flowfield of the meandering river
- Truncating the algorithm to calculate the flowfield

### Additional changes which are recommended:
- There is more commentary in the code and in my dissertation on why these changes are appropriate
- These changes will be made in a fork of this repo.  This repo is designed to replicate the dissertation work and so these improvements have not been made...yet
- in migration_model_TRH_CH2.m line 260, change the function call for computing ub, the downstream velocity variable, to use the constant hypothetical bank-full channel depth Do rather than the ever-changing channel depth D_ch(t).  
- in the flowfield_TRH.m file, fix the A = alphaa + 1, it should be -1 and should also be handled in a higher level (and more visible) function

This repo is specifically set up to be called by functions in my repositories MRDAST & MRLATS and will generate output files which are part of the analysis in the repositories MRDAST & MRVIST, and together will duplicate* the bulk of the dissertation work referenced above.  *It may be impossible to duplicate the simulation output without duplicating the machine which was a 32(?) bit MAC.  This is because the code written in c and compiled became a mexmaci64 file which is specific to the software and hardware set up.  Comparisons of the simulations show remarkable similarity over the first few thousand simulation years, before chaotic divergence.

A package of code files will be added to this repository at a later time: NOW COMPLETED!
