# -*- coding: utf-8 -*-
"""

MXL-model parameters and initial state

Created on Wed Mar 27 19:04:55 2019

@author: slauniai
"""

# ---  mxl model parameters and initial state    

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