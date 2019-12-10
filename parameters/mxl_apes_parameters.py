# -*- coding: utf-8 -*-
"""
Defines all parameters for running mxl-apes!

Created on Fri Apr  5 20:38:51 2019

@author: slauniai
"""
import numpy as np
from tools.utilities import lad_weibul #, lad_constant

gpara = {
    # 'pyAPES_path': r'c:\repositories\pyAPES_Samuli',
    'dt0' : 1800.0,  # timestep in forcing data file [s]
    'dt': 1800.0, # model timestep [s]
    'start_time' : "2014-07-03 08:00",  # start time of simulation [yyyy-mm-dd]
    'end_time' : "2014-07-03 18:00",  #"2018-01-01",  # end time of simulation [yyyy-mm-dd]
    'forc_filename' : "forcing/FIHy_mxl_forcing_2014.dat",  # 'Hyde_forcing_1997_2016.csv', # forcing data file*
    'results_directory':'results/',
    'variables': [# variable name, description [units], (dimensions)
                  
                  # copy of forcing variables
                  ['forcing_air_temperature', 'above canopy air temperature [degC]', ('date', 'simulation')],
                  ['forcing_precipitation', 'precipitation [m s-1]', ('date', 'simulation')],
                  ['forcing_pressure', 'ambient pressure [Pa]', ('date', 'simulation')],
                  ['forcing_h2o','H2O concentration [mol mol-1]', ('date', 'simulation')],
                  ['forcing_co2','CO2 concentration [ppm]', ('date', 'simulation')],
                  ['forcing_wind_speed','wind speed [m s-1]', ('date', 'simulation')],
                  ['forcing_friction_velocity','friction velocity [m s-1]', ('date', 'simulation')],
                  ['forcing_par','downward par (direct + diffuse) [W m-2]', ('date', 'simulation')],
                  ['forcing_nir','downward nir (direct + diffuse) [W m-2]', ('date', 'simulation')],
                  ['forcing_lw_in','downward lw [W m-2]', ('date', 'simulation')],                        
    
                  # canopy state and model control statistics
                  ['canopy_LAI','canopy LAI [m2 m-2]', ('date', 'simulation')],                      
                  ['canopy_lad','leaf area density [m3 m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_z', 'canopy model grid node elevations [m]', ('canopy')],
                  ['canopy_planttypes', 'canopy planttype names', ('planttype')],
                  ['canopy_WMA_assumption','WMA assumed (1=True, 0=False)', ('date', 'simulation')],
                  ['canopy_IterWMA', 'number of iterations [-]', ('date', 'simulation')],
                  ['canopy_energy_closure', 'energy closure in canopy [W m-2]', ('date', 'simulation')],
                  ['canopy_fr_source', 'Frsource in canopy [W m-2]', ('date', 'simulation')],  #error related to isothermal long-wave balance
                  
                  # micromet profiles and canopy-layer average leaf temperatures
                  ['canopy_h2o','H2O concentration [mol mol-1]', ('date', 'simulation', 'canopy')],
                  ['canopy_co2','CO2 concentration [ppm]', ('date', 'simulation', 'canopy')],
                  ['canopy_temperature','air temperature [degC]', ('date', 'simulation', 'canopy')],
                  ['canopy_wind_speed','canopy wind speed [m s-1]', ('date', 'simulation', 'canopy')],
                  ['canopy_friction_velocity','canopy friction velocity [m s-1]', ('date', 'simulation', 'canopy')],
                  ['canopy_Tleaf', 'leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                  ['canopy_Tleaf_wet', 'wet leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                  ['canopy_Tleaf_sl', 'sunlit leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                  ['canopy_Tleaf_sh', 'shaded leaf temperature [degC]', ('date', 'simulation', 'canopy')],
                  
                  # radiation 
                  ['canopy_sunlit_fraction','fraction of sunlit leafs [-]', ('date', 'simulation', 'canopy')],
                  ['canopy_SWnet', 'net shortwave radiation balance at canopy top [W m-2]', ('date', 'simulation')],
                  ['canopy_LWnet', 'net longwave radiation balance at canopy top [W m-2]', ('date', 'simulation')],
                  ['canopy_Rnet', 'net radiation balance at canopy top [W m-2]', ('date', 'simulation')],
                  # leaf scale, per m-2 leaf
                  ['canopy_leaf_net_LW', 'net leaf longwave radiation [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_leaf_net_SW', 'net leaf shortwave radiation [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_par_absorbed_sunlit', 'absorbed PAR of sunlit leaves [W m-2]', ('date', 'simulation', 'canopy')],               
                  ['canopy_par_absorbed_shaded', 'absorbed PAR of shaded leaves [W m-2]', ('date', 'simulation', 'canopy')],                     
                  ['canopy_nir_absorbed_sunlit', 'absorbed NIR of sunlit leaves [W m-2]', ('date', 'simulation', 'canopy')],               
                  ['canopy_nir_absorbed_shaded', 'absorbed NIR of shaded leaves [W m-2]', ('date', 'simulation', 'canopy')],
                  # vertical profiles, per m-2 ground
                  ['canopy_par_down', 'downward PAR [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_par_up', 'upward PAR [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_nir_down', 'downward NIR [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_nir_up', 'upward NIR [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_LW_down', 'downward LW [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_LW_up', 'upward LW [W m-2]', ('date', 'simulation', 'canopy')],
                  
                  # interception sub-model results
    #              ['canopy_interception', 'canopy interception [m s-1]', ('date', 'simulation')],
    #              ['canopy_interception_storage', 'canopy interception storage [m]', ('date', 'simulation')],
    #              ['canopy_evaporation', 'evaporation from interception storage [m s-1]', ('date', 'simulation')],
    #              ['canopy_condensation', 'condensation to canopy interception storage [m s-1]', ('date', 'simulation')],
    #              ['canopy_condensation_drip', 'condensation to canopy that drips [m s-1]', ('date', 'simulation')],
    #              ['canopy_transpiration','transpiration [m s-1]', ('date', 'simulation')],
    #              ['canopy_throughfall', 'throughfall to moss or snow [m s-1]', ('date', 'simulation')],
    #              ['canopy_evaporation_ml', 'evaporation from interception storage, profile (condensation incl.) [m s-1]', ('date', 'simulation', 'canopy')],
    #              ['canopy_throughfall_ml', 'throughfall within canopy, profile [m s-1]', ('date', 'simulation', 'canopy')],
    #              ['canopy_condensation_drip_ml', 'condensation drip within canopy, profile [m s-1]', ('date', 'simulation', 'canopy')],
    #              ['canopy_water_closure', 'interception model mass balance error [m s-1]', ('date', 'simulation')],
    
                  # ecosystem-level fluxes (at highest gridpoint, per m-2 ground)
                  ['canopy_SH', 'sensible heat flux [W m-2]', ('date', 'simulation')],
                  ['canopy_LE', 'latent heat flux [W m-2]', ('date', 'simulation')],
                  ['canopy_NEE', 'net ecosystem exchange [umol m-2 s-1]', ('date', 'simulation')],
                  ['canopy_GPP', 'ecosystem gross primary production [umol m-2 s-1]', ('date', 'simulation')],
                  ['canopy_Reco', 'ecosystem respiration [umol m-2 s-1]', ('date', 'simulation')],
                  ['canopy__transpiration', 'transpiration of all planttypes [m s-1]', ('date', 'simulation')],
    
                  # flux profiles within canopy
                  ['canopy_co2_flux', 'co2 flux [umol m-2 s-1]', ('date', 'simulation', 'canopy')],
                  ['canopy_latent_heat_flux', 'latent heat flux [W m-2]', ('date', 'simulation', 'canopy')],
                  ['canopy_sensible_heat_flux', 'sensible heat flux [W m-2]', ('date', 'simulation', 'canopy')],                      
                  
                  # root sink profile: Kersti - this should be taken out from Soil? Now len(rootsink) != len(Soil.z)
                  #['canopy_root_sink', 'root water uptake profile [m s-1]', ('date', 'simulation', 'soil')], 
                  
                  # planttype -specific outputs: lists of length 'planttype'
    #              ['pt_total_gpp', 'gross-primary productivity [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
    #              ['pt_total_dark_respiration', 'dark (or leaf + wood?) respiration [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
    #              ['pt_total_transpiration', 'transpiration [ms-1]', ('date', 'simulation', 'planttype')],                     
    #              ['pt_total_stomatal_conductance_h2o', 'stomatal conductance for H2O [mol m-2 s-1]', ('date', 'simulation', 'planttype')],   
    #              ['pt_total_boundary_conductance_h2o', 'leaf boundary layer conductance for H2O [mol m-2 s-1]', ('date', 'simulation', 'planttype')],
    #              ['pt_root_water_potential', 'root water potential [m?]', ('date', 'simulation', 'planttype')], # CHECK UNITS!!!            
    
                  # vertical profiles: lists of length 'planttype'; layers where lad == 0 are set to np.NaN
    #              ['pt_leaf_temperature', 'leaf temperature mean [degC]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_leaf_temperature_sunlit', 'leaf temperature, sunlit leaves [degC]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_leaf_temperature_shaded', 'leaf temperature, shaded leaves [degC]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_net_co2_sunlit', 'net co2 uptake, sunlit leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_net_co2_shaded', 'net co2 uptake, shaded leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_dark_respiration_sunlit', 'dark respiration, sunlit leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_dark_respiration_shaded', 'dark respiration, shaded leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_transpiration_sunlit', 'transpiration, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_transpiration_shaded', 'transpiration, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_latent_heat_sunlit', 'latent heat flux, sunlit leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_latent_heat_shaded', 'latent heat flux, shaded leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_sensible_heat_sunlit', 'sensible heat flux, sunlit leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_sensible_heat_shaded', 'sensible heat flux, shaded leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_stomatal_conductance_h2o_sunlit', 'stomatal conductance for H2O, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_stomatal_conductance_h2o_shaded', 'stomatal conductance for H2O, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_boundary_conductance_h2o_sunlit', 'boundary-layer conductance for H2O, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_boundary_conductance_h2o_shaded', 'boundary-layer conductance for H2O, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_leaf_internal_co2_sunlit', 'leaf internal CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_leaf_internal_co2_shaded', 'leaf internal CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_leaf_surface_co2_sunlit', 'leaf surface CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
    #              ['pt_leaf_surface_co2_shaded', 'leaf surface CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
                      
                  # mxl outputs
                  # variables and their units
                  ['mxl_h', 'mxl height [m]', ('date', 'simulation')],
                  ['mxl_h_lcl', 'lcl height [m]', ('date', 'simulation')],
                  ['mxl_T_lcl', 'lcl temperature [K]', ('date', 'simulation')],
                  ['mxl_theta', 'mxl potential temperature [K]', ('date', 'simulation')],
                  ['mxl_q', 'mxl specific humidity [kg/kg]', ('date', 'simulation')],
                  ['mxl_thetav', 'mxl virtual potential temperature [K]', ('date', 'simulation')],
                  ['mxl_ca', ' mxl CO2 mixing ratio [ppm]', ('date', 'simulation')],
                  ['mxl_Ws', 'subsidene velocity [ms-1]', ('date', 'simulation')],
                  ['mxl_wstar', 'convective velocity scale [ms-1]', ('date', 'simulation')],
                  ['mxl_sigmaw', 'turbulent velocity scale [ms-1]', ('date', 'simulation')],
                  ['mxl_u', ' mxl horizontal wind speed [ms-1]', ('date', 'simulation')],
                  ['mxl_U', 'mxl wind speed [ms-1]', ('date', 'simulation')],
                  ['mxl_vpd', 'surface vapor pressure deficit [kPa]', ('date', 'simulation')],
                  ['mxl_rh', 'surface relative humidity [-]', ('date', 'simulation')],
                  ['mxl_Psurf', 'surface pressure [kPa]', ('date', 'simulation')],
                
                  # entrainment zone
                  ['mxl_Pe', 'entrainment zone pressure [kPa]', ('date', 'simulation')],
                  ['mxl_Te', 'entrainment zone temperature [K]', ('date', 'simulation')],
                  ['mxl_theta_jump', 'potential temperature jump [K]', ('date', 'simulation')],
                  ['mxl_q_jump', 'specific humidity jump [kg/kg]', ('date', 'simulation')],
                  ['mxl_ca_jump', 'CO2 jump [ppm]', ('date', 'simulation')],
                  ['mxl_thetav_jump', 'virtual potential temperature jump [K]', ('date', 'simulation')],
                
                  # surface forcing to mxl
                  ['mxl_wthetas', 'surface kinematic heat flux [Kms-1]', ('date', 'simulation')],
                  ['mxl_wqs', 'surface kinematic moisture flux [kg kg-1 ms-1]', ('date', 'simulation')],
                  ['mxl_wcs', 'surface kinematic CO2 flux [ppm ms-1]', ('date', 'simulation')],
                  ['mxl_ust', 'surface friction velocity [ppm ms-1]', ('date', 'simulation')],
                  ]
        }

# ---  mxl model parameters and initial state    

mxlpara = {'dt': 1800.0, # s 
           'f': 1e-4,  # s-1
           'beta': 0.2, # closure constant
           'divU': 0.0, # large-scale subsidence due horizontal wind divergence s-1
           'ctr': {'Wind': True}
            }

mxl_ini = {'h': 100.,           # m
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

""" --- compile canopy model parameters ---"""

# site location
loc = {'lat': 61.51,  # latitude
       'lon': 24.0  # longitude
       }

# grid
grid = {'zmax': 30.0,  # heigth of grid from ground surface [m]
        'Nlayers': 100  # number of layers in grid [-]
        }

# --- control flags (True/False) ---
ctr = {'Eflow': True,  # ensemble flow
       'WMA': False, # well-mixed assumption
       'Ebal': True,  # computes leaf temperature by solving energy balance
       'WaterStress': None, #'Rew',  # Rew or PsiL or None
       'seasonal_LAI': False,  # account for seasonal LAI dynamics
       'pheno_cycle': False  # account for phenological cycle
       }

# --- micrometeo ---
micromet = {'zos': 0.01,  # forest floor roughness length [m]  -- not used?
            'dPdx': 0.00, #0.01,  # horizontal pressure gradient
            'Cd': 0.15,  # drag coefficient
            'Utop': 5.0,  # U
            'Ubot': 0.0,  # m/s, no-slip
            'Sc': {'T': 1.0, 'H2O': 1.0, 'CO2': 1.0}  # Schmidt numbers
            }

# --- radiation ---
radiation = {'clump': 0.7,  # clumping index [-]
             'leaf_angle': 1.0,  # leaf-angle distribution [-]
             'Par_alb': 0.12,  # shoot Par-albedo [-]
             'Nir_alb': 0.55,  # shoot NIR-albedo [-]
             'leaf_emi': 0.98
             }

# --- interception ---
interception = {'wmax': 0.2e-03,  # maximum interception storage capacity for rain [m per unit of LAI]  - Watanabe & Mizunani coniferous trees
                'wmaxsnow': 1.6e-03,  # maximum interception storage capacity for snow [m per unit of LAI]
                'w_ini': 0.0,  # initial canopy storage [m]
                'Tmin': 0.0,  # temperature below which all is snow [degC]
                'Tmax': 1.0,  # temperature above which all is water [degC]
                'leaf_orientation': 0.5 # leaf orientation factor for randomdly oriented leaves
                }

# create planttype; now LAI and phenology cycles are omitted
z = np.linspace(0, grid['zmax'], grid['Nlayers'])  # grid [m] above ground

pt1 = { 'name': 'tall_tree',                                         
        'LAImax': 3.0, # maximum annual LAI m2m-2 
        'lad': lad_weibul(z, LAI=1.0, h=15.0, hb=3.0, species='pine'),  # leaf-area density m2m-3
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
        }

# --- ground surface (for testing)
ground = {'soildepth': 10.0,
          'temperature': 10.0, 
          'emissivity': 0.98,
          'albedo': {'PAR':0.05, 'NIR': 0.2},
          # fluxes for testing
          'sensible_heat': 0.0,
          'latent_heat': 0.0, 
          'net_co2': 0.0
         }
        
cpara = {'loc': loc,
         'ctr': ctr,
         'grid': grid,
         'radiation': radiation,
         'micromet': micromet,
         'interception': interception,
         'planttypes': {pt1['name']: pt1},
         'ground': ground
         }



logging_configuration = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
                'default': {'format': '%(asctime)s %(levelname)s %(name)s %(message)s'},
                'model': {'format': '%(levelname)s %(name)s %(funcName)s %(message)s'},
                },
        'handlers': {
                'console': {
                        'class' : 'logging.StreamHandler',
                        'formatter': 'model',
                        'level': 'INFO'  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        },
                'file': {
                        'class': 'logging.FileHandler',
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'formatter': 'model',
                        'filename': 'pyAPES.log',
                        'mode': 'w',  # a == append, w == overwrite
                        },
                },
        'loggers': {
                'pyAPES': {
                        'handlers': ['file', 'console'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
                'canopy':{
                        'handlers': ['file', 'console'],
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
                'soil':{
                        'handlers': ['file', 'console'],
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
                },
        }

parallel_logging_configuration = {
        'version': 1,
        'formatters': {
                'default': {
                    'class': 'logging.Formatter',
                    'format': '%(asctime)s %(levelname)s %(name)s %(message)s'},
                'model': {
                    'class': 'logging.Formatter',
                    'format': '%(process)d %(levelname)s %(name)s %(funcName)s %(message)s'},
                },
        'handlers': {
                'console': {
                        'class' : 'logging.StreamHandler',
                        'formatter': 'model',
                        'level': 'INFO'  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        },
                'pyAPES_file': {
                        'class': 'logging.FileHandler',
                        'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'formatter': 'model',
                        'filename': 'pyAPES.log',
                        'mode': 'w',  # a == append, w == overwrite
                        },
                'parallelAPES_file': {
                        'class': 'logging.FileHandler',
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'formatter': 'default',
                        'filename': 'parallelAPES.log',
                        'mode': 'w',  # a == append, w == overwrite
                        },
        },
        'loggers': {
                'pyAPES': {
                        #'handlers': ['file'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagate': True,
                        },
        #        'canopy':{
        #                #'handlers': ['file'],
        #                'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
        #                },
        #        'soil':{
        #                #'handlers': ['file'],
        #                'level': 'DEBUG',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
        #                },
                },
        'root': {
                'level': 'DEBUG',
                'handlers': ['console', 'parallelAPES_file']
                }
        }