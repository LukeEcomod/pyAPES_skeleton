# -*- coding: utf-8 -*-
"""
GENERAL PARAMETERS FOR RUNNING pyAPES
"""

gpara = {
        'dt' : 1800.0,  # timestep in forcing data file [s]
        'start_time' : "2005-06-01",  # start time of simulation [yyyy-mm-dd]
        'end_time' : "2005-06-05",  # end time of simulation [yyyy-mm-dd]
        'forc_filename' : "Hyde_forcing_1997_2016.csv",  # forcing data file*
        'results_directory':'results/pyAPES/',

        # output variables listed here: note variable group and variable names case-sensitive.
        # comment unnecessary
        
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
                      ['canopy_interception', 'canopy interception [m s-1]', ('date', 'simulation')],
                      ['canopy_interception_storage', 'canopy interception storage [m]', ('date', 'simulation')],
                      ['canopy_evaporation', 'evaporation from interception storage [m s-1]', ('date', 'simulation')],
                      ['canopy_condensation', 'condensation to canopy interception storage [m s-1]', ('date', 'simulation')],
                      ['canopy_condensation_drip', 'condensation to canopy that drips [m s-1]', ('date', 'simulation')],
                      ['canopy_transpiration','transpiration [m s-1]', ('date', 'simulation')],
                      ['canopy_throughfall', 'throughfall to moss or snow [m s-1]', ('date', 'simulation')],
                      ['canopy_evaporation_ml', 'evaporation from interception storage, profile (condensation incl.) [m s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_throughfall_ml', 'throughfall within canopy, profile [m s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_condensation_drip_ml', 'condensation drip within canopy, profile [m s-1]', ('date', 'simulation', 'canopy')],
                      ['canopy_water_closure', 'interception model mass balance error [m s-1]', ('date', 'simulation')],

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
                      ['pt_total_gpp', 'gross-primary productivity [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['pt_total_dark_respiration', 'dark (or leaf + wood?) respiration [umol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['pt_total_transpiration', 'transpiration [ms-1]', ('date', 'simulation', 'planttype')],                     
                      ['pt_total_stomatal_conductance_h2o', 'stomatal conductance for H2O [mol m-2 s-1]', ('date', 'simulation', 'planttype')],   
                      ['pt_total_boundary_conductance_h2o', 'leaf boundary layer conductance for H2O [mol m-2 s-1]', ('date', 'simulation', 'planttype')],
                      ['pt_root_water_potential', 'root water potential [m?]', ('date', 'simulation', 'planttype')], # CHECK UNITS!!!            

                      # vertical profiles: lists of length 'planttype'; layers where lad == 0 are set to np.NaN
                      ['pt_leaf_temperature', 'leaf temperature mean [degC]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_leaf_temperature_sunlit', 'leaf temperature, sunlit leaves [degC]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_leaf_temperature_shaded', 'leaf temperature, shaded leaves [degC]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_net_co2_sunlit', 'net co2 uptake, sunlit leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_net_co2_shaded', 'net co2 uptake, shaded leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_dark_respiration_sunlit', 'dark respiration, sunlit leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_dark_respiration_shaded', 'dark respiration, shaded leaves [umol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_transpiration_sunlit', 'transpiration, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_transpiration_shaded', 'transpiration, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_latent_heat_sunlit', 'latent heat flux, sunlit leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_latent_heat_shaded', 'latent heat flux, shaded leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_sensible_heat_sunlit', 'sensible heat flux, sunlit leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_sensible_heat_shaded', 'sensible heat flux, shaded leaves [W m-2 (leaf)]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_stomatal_conductance_h2o_sunlit', 'stomatal conductance for H2O, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_stomatal_conductance_h2o_shaded', 'stomatal conductance for H2O, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_boundary_conductance_h2o_sunlit', 'boundary-layer conductance for H2O, sunlit leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_boundary_conductance_h2o_shaded', 'boundary-layer conductance for H2O, shaded leaves [mol m-2 (leaf) s-1]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_leaf_internal_co2_sunlit', 'leaf internal CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_leaf_internal_co2_shaded', 'leaf internal CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_leaf_surface_co2_sunlit', 'leaf surface CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
                      ['pt_leaf_surface_co2_shaded', 'leaf surface CO2 mixing ratio, shaded leaves [ppm]', ('date', 'simulation', 'planttype', 'canopy')],
                      
                      # soil model state and fluxes
                      ['soil_z', 'soil model grid node elevations [m]', ('soil')],
                      ['soil_water_potential','soil water potential [m]', ('date', 'simulation', 'soil')],
                      ['soil_pond_storage', 'pond storage [m]', ('date', 'simulation')],
                      ['soil_ground_water_level', 'ground water level [m]', ('date', 'simulation')],
                      ['soil_infiltration', 'infiltration [m s-1]', ('date', 'simulation')],
                      ['soil_surface_runoff', 'surface runoff [m s-1]', ('date', 'simulation')],
                      ['soil_evaporation', 'evaporation from soil surface [m s-1]', ('date', 'simulation')],
                      ['soil_transpiration', 'transpiration from soil [m s-1]', ('date', 'simulation')],
                      ['soil_drainage', 'subsurface drainage [m s-1]', ('date', 'simulation')],
                      ['soil_temperature', 'soil temperature [degC]', ('date', 'simulation', 'soil')],
                      ['soil_volumetric_water_content', 'soil water content [m3/m3]', ('date', 'simulation', 'soil')],
                      ['soil_volumetric_ice_content', 'soil ice content [m3/m3]', ('date', 'simulation', 'soil')],
                      ['soil_heat_flux', 'soil heat flux [W m-2]', ('date', 'simulation', 'soil')],
                      ['soil_thermal_conductivity', 'thermal conductivity [W m-1 K-1]', ('date', 'simulation', 'soil')],
                      ['soil_water_closure', 'soil water balance error [m s-1]', ('date', 'simulation')],
                      ['soil_energy_closure', 'soil heat balance error [W m-2]', ('date', 'simulation')],
 
                      # forest floor outputs
                      ['ffloor_potential_infiltration', 'potential infiltration to soil [m s-1]', ('date', 'simulation')],
                      ['ffloor_snow_water_equivalent', 'snow water equivalent [m]', ('date', 'simulation')],
                      ['ffloor_ground_heat', 'ground heat flux (forest floor) [W m-2]', ('date', 'simulation')],
                      ['ffloor_sensible_heat', 'sensible heat flux (forest floor) [W m-2]', ('date', 'simulation')],
                      ['ffloor_latent_heat', 'latent heat flux (forest floor) [W m-2]', ('date', 'simulation')],
                      ['ffloor_snow_water_closure', "water balance error (snowcover) [m s-1]", ('date', 'simulation')],
                      ['ffloor_bryo_water_closure', "water balance error (bryophytes) [m s-1]", ('date', 'simulation')],
                      ['ffloor_bryo_energy_closure', "energy balance error (bryophytes) [W m-2]", ('date', 'simulation')],
                      ['ffloor_soil_energy_closure', "energy balance error (soil) [W m-2]", ('date', 'simulation')],
                      ['ffloor_bryo_carbon_pool', 'carbon pool (bryophyte) [kg C m-2]', ('date', 'simulation')],
                      ['ffloor_bryo_photosynthesis', 'photosynthesis rate (bryophyte) [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_bryo_respiration', 'respiration rate (bryophyte) [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_litter_respiration', 'respiration rate (litter) [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_respiration', 'respiration rate [umol m-2(ground) s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation', 'evaporation (forest floor) [m s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation_bryo', 'evaporation (bryophytes) [m s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation_litter', 'evaporation (litter) [m s-1]', ('date', 'simulation')],
                      ['ffloor_evaporation_soil', 'evaporation (soil) [m s-1]', ('date', 'simulation')],
                      ['ffloor_temperature', 'temperature (forest floor) [degC]', ('date', 'simulation')],
                      ['ffloor_litter_temperature', 'temperature (litter) [degC]', ('date', 'simulation')],
                      ['ffloor_bryo_temperature', 'temperature (bryophyte) [degC]', ('date', 'simulation')],
                      ['ffloor_soil_temperature', 'temperature (soil) [degC]', ('date', 'simulation')],
                      ['ffloor_bryo_water_storage', 'water storage (bryophytes) [kg m-2]', ('date', 'simulation')],
                      ['ffloor_litter_water_storage', 'water storage (litter) [kg m-2]', ('date', 'simulation')],
                      ['ffloor_capillar_rise', 'capillary rise to bryophyte layer [m s-1]', ('date', 'simulation')],
                      ]}

# --- logger configuration. Antti / Kersti: add option to define logger output file name?
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

# for parallel simulations
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
                        'level': 'WARNING',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
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
                'canopy':{
                        #'handlers': ['file'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagete': True,
                        },
                'soil':{
                        #'handlers': ['file'],
                        'level': 'INFO',  # CRITICAL, ERROR, WARNING, INFO, DEBUG
                        'propagete': True,
                },
        },
        'root': {
                'level': 'INFO',
                'handlers': ['console', 'parallelAPES_file']
        }
    }
