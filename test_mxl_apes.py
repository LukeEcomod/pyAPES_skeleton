# -*- coding: utf-8 -*-
"""
**** Runs mxl_apes ****

Created on Fri Apr  5 20:38:51 2019

@author: slauniai
"""
import numpy as np
from matplotlib import pyplot as plt
from mxl_apes import driver
#from pyAPES_utilities.plotting import plot_timeseries_xr
#from tools.iotools import read_results


# %% Run model

# yksi ajo
#outputfile=driver(create_ncf=True)


results, model = driver(create_ncf=False, parametersets={})

#run1 = driver(create_ncf=False, parametersets={})
#%%
t = model.forcing.index.hour.values + model.forcing.index.minute.values / 60
res = results[0]
fix, ax = plt.subplots(3,2)


ax[0,0].plot(t, res['canopy_SH'], label='H')
ax[0,0].plot(t, res['canopy_LE'], label='LE')
#rn = res['canopy_SWnet'] + res['canopy_LWnet']
#ax[0,0].plot(t, rn, label='Rn')
ax[0,0].legend()
ax[0,0].set_ylabel('surface flux [Wm-2]')

ax[0,1].plot(t, res['mxl_h'], label='mxl h');
ax[0,1].set_ylabel('mxl h [m]')

ax[1,0].plot(t, res['mxl_theta'] - 273.15, 'r-', label='asl top')
ax[1,0].plot(t, res['canopy_temperature'])
ax[1,0].legend()
ax[1,0].set_ylabel('Ta [degC]')

ax[1,1].plot(t, res['mxl_vpd'])
ax[1,1].set_ylabel('VPD [kPa]')

ax[2,0].plot(t, res['mxl_ca'], 'r-', label='asl top')
ax[2,0].plot(t, res['canopy_co2'])
ax[2,0].legend()
ax[2,0].set_ylabel('CO2 [ppm]')
