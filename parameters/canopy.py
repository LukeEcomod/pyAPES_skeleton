# -*- coding: utf-8 -*-
"""
CANOPY MODEL PARAMETERS for pyAPES
imports forestfloor parameters from parameters.forestfloor
"""

import numpy as np
from tools.utilities import lad_weibul #, lad_constant

from .forestfloor import forestfloor

# site location
loc = {'lat': 61.51,  # latitude
       'lon': 24.0  # longitude
       }

# grid
grid = {'zmax': 25.0,  # heigth of grid from ground surface [m]
        'Nlayers': 101  # number of layers in grid [-]
        }

# --- control flags (True/False) ---
ctr = {'Eflow': True,  # ensemble flow
       'WMA': False,  #True,  #  well-mixed assumption
       'Ebal': True,  #False,  #  computes leaf temperature by solving energy balance
       'WaterStress': 'Rew',  #'PsiL',  # Rew or PsiL or None
       'seasonal_LAI': True,  # account for seasonal LAI dynamics
       'pheno_cycle': True  # account for phenological cycle
       }

# --- micrometeo ---
micromet = {'zos': 0.01,  # forest floor roughness length [m]  -- not used?
            'dPdx': 0.01,  # horizontal pressure gradient
            'Cd': 0.15,  # drag coefficient
            'Utop': 5.0,  # ensemble U/ustar
            'Ubot': 0.0,  # lower boundary
            'Sc': {'T': 2.0, 'H2O': 2.0, 'CO2': 2.0}  # Schmidt numbers
            }

# --- radiation ---
radiation = {'clump': 0.7,  # clumping index [-]
             'leaf_angle': 1.0,  # leaf-angle distribution [-]
             'Par_alb': 0.12,  # shoot Par-albedo [-]
             'Nir_alb': 0.55,  # shoot NIR-albedo [-]
             'leaf_emi': 0.98
             }

# --- interception ---
interception = {'wmax': 0.2,  # maximum interception storage capacity for rain [kg m-2 per unit of LAI]  - Watanabe & Mizunani coniferous trees
                'wmaxsnow': 0.8,  # maximum interception storage capacity for snow [kg m-2 per unit of LAI] - about 4 * wmax (Koivusalo & Kokkonen 2002)
                'w_ini': 0.0,  # initial canopy storage [kg m-2]
                'Tmin': 0.0,  # temperature below which all is snow [degC]
                'Tmax': 2.0,  # temperature above which all is water [degC]- Koivusalo & Kokkonen 2002
                'leaf_orientation': 0.5, # leaf orientation factor for randomdly oriented leaves
                }

# --- define two planttypes ---

z = np.linspace(0, grid['zmax'], grid['Nlayers'])  # grid [m] above ground

pt1 = { 'name': 'plant1',
        'LAImax': 3.0, # maximum annual LAI m2m-2
        'lad': lad_weibul(z, LAI=1.0, h=15.0, hb=3.0, species='pine'),  # leaf-area density m2m-3
        # cycle of photosynthetic activity
        'phenop': {
            'Xo': 0.0,
            'fmin': 0.1,
            'Tbase': -4.67,  # Kolari 2007
            'tau': 8.33,  # Kolari 2007
            'smax': 18.0  # Kolari 2014
            },
        # cycle of LAI
        'laip': {
            'lai_min': 0.8,
            'lai_ini': None,
            'DDsum0': 0.0,
            'Tbase': 5.0,
            'ddo': 45.0,
            'ddmat': 250.0,
            'sdl': 12.0,
            'sdur': 30.0
            },
        # A-gs model
        'photop': {
            'Vcmax': 45.0,
            'Jmax': 85.0,  # 1.97*Vcmax (Kattge and Knorr, 2007)
            'Rd': 0.9,  # 0.023*Vcmax
            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                'Vcmax': [72., 200., 649.],
                'Jmax': [50., 200., 646.],
                'Rd': [33.0]
                },
            'alpha': 0.2,   # quantum efficiency parameter -
            'theta': 0.7,   # curvature parameter
            'g1': 2.1,      # stomatal slope kPa^(0.5)
            'g0': 5.0e-3,   # residual conductance mol m-2 s-1
            'kn': 0.5,      # nitrogen attenuation coefficient -
            'beta': 0.95,   # co-limitation parameter -
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.02,     # leaf length scale m
            },
        # root zone
        'rootp': {
            'root_depth': 0.5, # rooting depth [m]
            'beta': 0.943, # root distribution shape [-]
            'RAI_LAI_multiplier': 2.0, # fine-root to leaf-area ratio [-]
            'fine_radius': 2.0e-3, # [m]
            'root_cond': 5.0e8, # [s]
            }
        }


pt2 = { 'name': 'plant2',
        'LAImax': 1.0, # maximum annual LAI m2m-2
        'lad': lad_weibul(z, LAI=1.0, h=5.0, hb=0.0, species='pine'),  # leaf-area density m2m-3
        # cycle of photosynthetic activity
        'phenop': {
            'Xo': 0.0,
            'fmin': 0.1,
            'Tbase': -4.67,  # Kolari 2007
            'tau': 8.33,  # Kolari 2007
            'smax': 18.0  # Kolari 2014
            },
        # cycle of LAI
        'laip': {
            'lai_min': 0.8,
            'lai_ini': None,
            'DDsum0': 0.0,
            'Tbase': 5.0,
            'ddo': 45.0,
            'ddmat': 250.0,
            'sdl': 12.0,
            'sdur': 30.0
            },
        # A-gs model
        'photop': {
            'Vcmax': 45.0,
            'Jmax': 85.0,  # 1.97*Vcmax (Kattge and Knorr, 2007)
            'Rd': 0.9,  # 0.023*Vcmax
            'tresp': { # temperature response parameters (Kattge and Knorr, 2007)
                'Vcmax': [72., 200., 649.],
                'Jmax': [50., 200., 646.],
                'Rd': [33.0]
                },
            'alpha': 0.2,   # quantum efficiency parameter -
            'theta': 0.7,   # curvature parameter
            'g1': 2.1,      # stomatal slope kPa^(0.5)
            'g0': 5.0e-3,   # residual conductance mol m-2 s-1
            'kn': 0.5,      # nitrogen attenuation coefficient -
            'beta': 0.95,   # co-limitation parameter -
            'drp': [0.39, 0.83, 0.31, 3.0] # Rew-based drought response
            },
        'leafp': {
            'lt': 0.02,     # leaf length scale m
            },
        # root zone
        'rootp': {
            'root_depth': 0.5, # rooting depth [m]
            'beta': 0.943, # root distribution shape [-]
            'RAI_LAI_multiplier': 2.0, # fine-root to leaf-area ratio [-]
            'fine_radius': 2.0e-3, # [m]
            'root_cond': 5.0e8, # [s]
            }
        }

# --- canopy-model parameter dictionary

cpara = {'loc': loc,
         'ctr': ctr,
         'grid': grid,
         'radiation': radiation,
         'micromet': micromet,
         'interception': interception,
         'planttypes': {'plant1': pt1, 'plant2': pt2},
         'forestfloor': forestfloor
         }
