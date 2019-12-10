# -*- coding: utf-8 -*-
"""
SCRIPT FOR TESTING MIXED-LAYER MODEL USING ONE DAYTIME DATA FROM HYYTIÄLÄ

Created on Thu Sep 20 14:55:56 2018

@author: slauniai
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# --- import constants from module mxl
from mxl import CP_AIR_MASS, MAIR_DRY, MH2O, NT, R
# --- import model class and utility functions from module mxl
from mxl import MXLmodel, air_density, read_forcing

#import mxl # mixed-layer model

EPS = np.finfo(float).eps  # machine epsilon

# ---- parameters & initial conditions for mxl model

wdir = r'c:\repositories\pyAPES-MXL'
os.chdir(wdir)
print('---- working dir: ' + os.getcwd())

print('---- reading forcing ---')
# --- read forcing for testing mxl growth
ffile = r'forcing\forc_2010_d184_186.dat'

fday = '2010-07-03'
lday = '2010-07-03'

# read forcing into pd.dataframe
dat, tvec = read_forcing(ffile)
dat.index = tvec

# select one day and conditions when H > 0
forc = dat.loc['2010-07-03' : '2010-07-03'][['H','E', 'NEE', 'P', 'Ta','U','ust']]
# ix = forc['H'] < 0
#forc['H'][ix] = 0.0

# this selects periods when H >0 from forcing.
forc = forc[forc['H'] > 0]

# plot figure
plt.figure()
plt.plot(forc['H']); plt.ylabel('H (Wm-2)')

forc['P'] = 1e2*forc['P']  # Pa
forc['Ta'] += NT  # K

# -- convert units of surface fluxes: 
# originally in H Wm-2, E mmol m-2s-1, NEE umolm-2s-1, P hPa, Ta degC, u, ust ms-1

F_h = forc['H'].values / (1.2*CP_AIR_MASS)  # K ms-1, use constant air density
MAIR_MOLAR = (forc['P'] / (R*forc['Ta'])).values
F_h2o = 1e-3*forc['E'].values / MAIR_MOLAR *MH2O / MAIR_DRY # mmol m-2s-1 --> kg/kg ms-1

F_co2 = forc['NEE'].values / MAIR_MOLAR  # ppm ms-1, <0 is sink

ustar = forc['ust'].values # m s-1

# --- initialize mxl model    

mxlpara = {'dt': 1800.0, # s 
           'f': 1e-4,  # s-1
           'beta': 0.2, # closure constant
           'divU': 0.0, # large-scale subsidence due horizontal wind divergence s-1
           'ctr': {'Wind': True}
            }

ini = {'h': 200.,           # m
       'theta': 288.0,      # K
       'q': 8.0e-3,         # kg kg-1
       'ca': 422.0,         # ppm
       'theta_jump': 1.0,   # K
       'gamma_theta': 6e-3, # K m-1
       'q_jump': -1.0e-3,   # kg kg-1
       'gamma_q': -1.45e-6, # kg kg-1 m-1
       'ca_jump': -40.0,    # ppm
       'gamma_ca': 0.0,     # ppm m-1
       'u': 5.0,            # m s-1
       'u_jump': 8.0,       # m s-1, geostrophic wind is u_jump + u
       'gamma_u': 0.0,      # s-1
       'Psurf': 101.3       # kPa
      }

print('---- creating MXL-model----')
# --- Create model instance
run1 = MXLmodel(ini, mxlpara)
## print run1.__dict_

print('---- running MXL-model----')
nsteps = len(forc['H']) # len(F_h)
tt = 30.*np.arange(0,nsteps) # time vector, min

# initialize results dictionany, fill with NaN's
res = {'h': np.ones(nsteps)*np.NaN, 'theta': np.ones(nsteps)*np.NaN, 
       'q':np.ones(nsteps)*np.NaN, 'ca': np.ones(nsteps)*np.NaN,
       'h_lcl': np.ones(nsteps)*np.NaN, 'vpd': np.ones(nsteps)*np.NaN,
       'U': np.ones(nsteps)*np.NaN, 'u': np.ones(nsteps)*np.NaN
      }

# run model for nsteps
for k in range(nsteps):
    run1.run_timestep(F_h[k], F_h2o[k], F_co2[k], ustar[k])

    res['h'][k] = run1.h
    res['theta'][k] = run1.theta
    res['q'][k] = run1.q
    res['ca'][k] = run1.ca
    res['h_lcl'][k] = run1.h_lcl
    res['vpd'][k] = run1.vpd
    res['U'][k] = run1.U
    res['u'][k] = run1.u

print('---- making graphs----')

plt.figure(1)
plt.subplot(221); plt.plot(tt, forc['H'], 'r'); plt.ylabel('H (Wm-2)')
plt.subplot(222); plt.plot(tt, forc['E'], 'b'); plt.ylabel('E (mmol m-2)')
plt.subplot(223); plt.plot(tt, forc['NEE'], 'g'); plt.ylabel('NEE (umol m-2 s-1)')
plt.subplot(224); plt.plot(tt, forc['ust'], 'k'); plt.ylabel('ustar (m s-1)')

plt.savefig('forc.png')

plt.figure(2)
plt.subplot(321); plt.plot(tt, res['h'], label='mxl'); plt.plot(tt, res['h_lcl'], label='lcl'); 
plt.ylabel('mxl and lcl height (m)'); plt.legend(fontsize=8)
plt.subplot(322); plt.plot(tt, res['theta']); plt.ylabel('Theta (K)')
plt.subplot(323); plt.plot(tt, res['q']); plt.ylabel('q (kg/kg)')
plt.subplot(324); plt.plot(tt, res['vpd']); plt.ylabel('vpd (kPa)')
plt.subplot(325); plt.plot(tt, res['ca']); plt.ylabel('ca (ppm)'); plt.xlabel('time (min)')
plt.subplot(326); plt.plot(tt, res['U'], label='U'); plt.plot(tt, res['u'], label='u horiz')
plt.ylabel('velocity (ms-1)'); plt.xlabel('time (min)'); plt.legend(fontsize=8)
plt.savefig('mxl-test.png')

print('---- done ! ----')