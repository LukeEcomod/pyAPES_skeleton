# README

pyAPES_skeleton branch 'mxl' couples pyAPES to convective mixed boundary layer slab model.

\mxl

Simple implementation of Mixed Boundary Layer (MXL) model based on:

Janssens & Pozzer, 2015. GeoSci Model Dev. 8, 453 - 471.
Vilà-Guerau de Arellano et al. 2015. Atmospheric Boundary Layer: Integrating Air Chemistry and Land Interactions. Cambridge University Press, New York, 2015, 265 pp
Stull, 1998; Siqueira et al. 2009. J. Hydrometeorol. 10, 96-112.
   
See 'mxl_demo.py' and 'mxl_demo_bowen_ratio.py'

Stand-alone from other packages

Samuli Launiainen 9/2018 

mxl_apes.py: couples mxl-model with solution of plant canopy-asl interactions. Currently assumes soil surface state and fluxes = constant & surface layer height = 0.1 x mxl-height

\mxl_canopy

Contains canopy model, simplified planttype-code and 1st order closure models for canopy-surface layer momentum & scalar budgets. 
Leaf energy balance is solved iteratively with solution of scalar profiles throughout the asl.
Note that dict asl_profs give momentum and scalar profiles from ground to asl top but these are not saved in outputs yet. Netcdf-writing needs to be updated.

Todo:

* implement intial states for MXL growth in morning
* what happens when H < 0? Now no collapse of MXL.
* Test mxl-growth against ABL height from ceilometer & soundings when forced by measured surface fluxes.








