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
LMOLAR = 44100.0  # J/mol latent heat of vaporization at 20 deg C

# ---- parameters & initial conditions for mxl model

wdir = 'c:\\Repositories\\mxl\\'
os.chdir(wdir)
print('---- working dir: ' + os.getcwd())

print('---- reading forcing ---')
# --- read forcing for testing mxl growth
ffile = 'forc_2010_d184_186.dat'

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

# -- convert units of surface fluxes; first adjust H and LE with bowen ratios
MAIR_MOLAR = (forc['P'] / (R*forc['Ta'])).values

H = forc['H'].values
LE = 1e-3*forc['E'].values * LMOLAR
AE = H + LE

beta1 = 0.3
le1 = AE / (beta1 + 1)
h1 = AE - le1
F_h2o1 = le1 / (beta1 + 1) / LMOLAR / MAIR_MOLAR *MH2O / MAIR_DRY
F_h1 = h1 / (1.2*CP_AIR_MASS)

beta2 = 1.35
le2 = AE / (beta2 + 1)
h2 = AE - le2
F_h2o2 = le2 / LMOLAR / MAIR_MOLAR *MH2O / MAIR_DRY
F_h2 = h2 / (1.2*CP_AIR_MASS)

# originally in H Wm-2, E mmol m-2s-1, NEE umolm-2s-1, P hPa, Ta degC, u, ust ms-1

#F_h = forc['H'].values / (1.2*CP_AIR_MASS)  # K ms-1, use constant air density
#MAIR_MOLAR = (forc['P'] / (R*forc['Ta'])).values
#F_h2o = 1e-3*forc['E'].values / MAIR_MOLAR *MH2O / MAIR_DRY # mmol m-2s-1 --> kg/kg ms-1
#
F_co2 = forc['NEE'].values / MAIR_MOLAR  # ppm ms-1, <0 is sink
#
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
run2 = MXLmodel(ini, mxlpara)
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

res2 = {'h': np.ones(nsteps)*np.NaN, 'theta': np.ones(nsteps)*np.NaN, 
       'q':np.ones(nsteps)*np.NaN, 'ca': np.ones(nsteps)*np.NaN,
       'h_lcl': np.ones(nsteps)*np.NaN, 'vpd': np.ones(nsteps)*np.NaN,
       'U': np.ones(nsteps)*np.NaN, 'u': np.ones(nsteps)*np.NaN
      }

# run model for nsteps
for k in range(nsteps):
    # run1
    run1.run_timestep(F_h1[k], F_h2o1[k], F_co2[k], ustar[k])

    res['h'][k] = run1.h
    res['theta'][k] = run1.theta
    res['q'][k] = run1.q
    res['ca'][k] = run1.ca
    res['h_lcl'][k] = run1.h_lcl
    res['vpd'][k] = run1.vpd
    res['U'][k] = run1.U
    res['u'][k] = run1.u

    run2.run_timestep(F_h2[k], F_h2o2[k], F_co2[k], ustar[k])

    res2['h'][k] = run2.h
    res2['theta'][k] = run2.theta
    res2['q'][k] = run2.q
    res2['ca'][k] = run2.ca
    res2['h_lcl'][k] = run2.h_lcl
    res2['vpd'][k] = run2.vpd
    res2['U'][k] = run2.U
    res2['u'][k] = run2.u
print('---- making graphs----')

plt.figure(1)
plt.subplot(221); plt.plot(tt, h1, tt, h2); plt.ylabel('H (Wm-2)')
plt.subplot(222); plt.plot(tt, le1, tt, le2); plt.ylabel('LE (Wm-2)')
plt.subplot(223); plt.plot(tt, forc['NEE'], 'g'); plt.ylabel('NEE (umol m-2 s-1)')
plt.subplot(224); plt.plot(tt, forc['ust'], 'k'); plt.ylabel('ustar (m s-1)')

plt.savefig('forc.png')

plt.figure(2)
plt.subplot(321); plt.plot(tt, res['h'], label=r'$\beta$ = 0.3'); #plt.plot(tt, res['h_lcl'], label='lcl'); 
plt.ylabel('mxl and lcl height (m)'); plt.legend(fontsize=8)
plt.subplot(322); plt.plot(tt, res['theta']); plt.ylabel('Theta (K)')
plt.subplot(323); plt.plot(tt, res['q']); plt.ylabel('q (kg/kg)')
plt.subplot(324); plt.plot(tt, res['vpd']); plt.ylabel('vpd (kPa)')
plt.subplot(325); plt.plot(tt, res['ca']); plt.ylabel('ca (ppm)'); plt.xlabel('time (min)')
plt.subplot(326); plt.plot(tt, res['U'], label='U'); # plt.plot(tt, res['u'], label='u horiz')
plt.ylabel('velocity (ms-1)'); plt.xlabel('time (min)'); plt.legend(fontsize=8)
# beta=1.3
plt.subplot(321); plt.plot(tt, res2['h'], '-', label=r'$\beta$ = 1.35'); #plt.plot(tt, res2['h_lcl'], label='lcl'); 
plt.ylabel('mxl and lcl height (m)'); plt.legend(fontsize=8)
plt.subplot(322); plt.plot(tt, res2['theta'], '-'); plt.ylabel('Theta (K)')
plt.subplot(323); plt.plot(tt, res2['q'], '-'); plt.ylabel('q (kg/kg)')
plt.subplot(324); plt.plot(tt, res2['vpd'], '-'); plt.ylabel('vpd (kPa)')
plt.subplot(325); plt.plot(tt, res2['ca'], '-'); plt.ylabel('ca (ppm)'); plt.xlabel('time (min)')
plt.subplot(326); plt.plot(tt, res2['U'], '-', label='U'); # plt.plot(tt, res2['u'], label='u horiz')
plt.ylabel('velocity (ms-1)'); plt.xlabel('time (min)'); plt.legend(fontsize=8)

plt.subplot(321)
plt.plot(tt, res['h_lcl'], '--', label=r'$\beta$ = 0.3'); 
plt.plot(tt, res2['h_lcl'], '--', label=r'$\beta$ = 1.35') 
plt.legend(fontsize=6)
plt.title(r'$\beta$ = 0.3 and 1.35, respectively')
plt.savefig('mxl-variable-bowen.png')

print('---- done ! ----')