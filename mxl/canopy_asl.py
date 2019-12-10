#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: canopy
    :synopsis: pyAPES-MXL-model component
.. moduleauthor:: Kersti Haahti & Samuli Launiainen
LAST EDIT 4.4.19 / SL
"""


import logging
# from copy import deepcopy
import numpy as np

from .mxl import MOLAR_MASS_H2O, EPS, LATENT_HEAT, STEFAN_BOLTZMANN, DEG_TO_KELVIN

# changes to absolute imports: these are all inside canopy-module
from canopy.radiation import Radiation
from .canopy_asl_flow import Micromet
from .interception import Interception
from .planttype.planttype import PlantType

#from radiation import Radiation
#from canopy_asl_flow import Micromet
#from interception import Interception
#from planttype.planttype import PlantType
logger = logging.getLogger(__name__)


class CanopyModel(object):
    r""" Main model for Topsoil - Canopy -component of pyAPES-MXL
    """

    def __init__(self, cpara):
        r""" Initializesobject and submodel objects using given parameters.

        Args:
            cpara (dict):
                'ctr' (dict): switches and specifications for computation
                    'Eflow' (bool): ensemble flow assumption (False solves U_normed on daily timestep)
                    'WMA' (bool): well-mixed assumption (False solves H2O, CO2, T)
                    'Ebal' (bool): True solves energy balance
                    'WaterStress' (str): account for water stress in planttypes via 'Rew', 'PsiL' or None omits
                    'seasonal_LAI' (bool): account for seasonal LAI dynamics
                    'pheno_cycle' (bool): account for phenological cycle
                'grid' (dict):
                    'zmax': heigth of grid from ground surface [m]
                    'Nlayers': number of layers in grid [-]
                'radiation' (dict): radiation model parameters
                'micromet' (dict): micromet model parameters --- ROUGHNESS HEIGHT SHOULD from ffloor?
                'interception' (dict): interception and snow model parameters
                'planttypes' (list):
                    i. (dict): properties of planttype i
#                'forestfloor': forestfloor parameters
#                    'bryophytes'(list):
#                        i. (dict): properties of bryotype i
#                    'baresoil' (dict): baresoil parameters
#                    'snowpack' (dict): smow model parameters
#                    'initial_conditions': initial conditions for forest floor
            dz_soil (array): thickness of soilprofile layers, needed for rootzone

        Returns:
            self (object):
                .location (dict):
                    'lat': latitude [deg]
                    'lon': longitude [deg]
                .z (array): canopy model nodes, height from soil surface (= 0.0) [m]
                .dz (float): thickness of canopy layers [m]
                .ones (array): dummy of len(z)
                .Switch_Eflow (bool): ensemble flow assumption
                .Switch_WMA (bool): well-mixed assumption (H2O, CO2, T)
                .Switch_Ebal (bool): solve energy balance
                .LAI (array): total leaf area index [m2 m-2]
                .lad (array): total leaf area density [m2 m-3]
                .hc (float): canopy heigth [m]
                .rad (array): normalized total fine root density distribution [-]
                .planttypes (list):
                    i. (object): planttype object i
                .radiation (object): radiation model (SW, LW)
                .micromet (object): micromet model (U, H2O, CO2, T)
                .interception (object): interception model
                .ground_surface (object): ground surface (replaces forest floor)
                #.forestfloor (object): forest floor object (bryotype/baresoil/snow)
        """

        # --- grid ---
        self.z = np.linspace(0, cpara['grid']['zmax'], cpara['grid']['Nlayers'])  # grid [m] above ground
        self.dz = self.z[1] - self.z[0]  # gridsize [m]
        self.ones = np.ones(len(self.z))  # dummy

        # --- switches ---
        self.Switch_Eflow = cpara['ctr']['Eflow']
        self.Switch_WMA = cpara['ctr']['WMA']
        self.Switch_Ebal = cpara['ctr']['Ebal']

        logger.info('Eflow: %s, WMA: %s, Ebal: %s',
                    self.Switch_Eflow,
                    self.Switch_WMA,
                    self.Switch_Ebal)

        # --- Plant types (with phenological models) ---
        dz_soil = cpara['ground']['soildepth']
        ptypes = []
        for pt in cpara['planttypes']:        
            ptypes.append(PlantType(self.z, cpara['planttypes'][pt], dz_soil, ctr=cpara['ctr'], loc=cpara['loc']))
        self.planttypes = ptypes
    
        # --- stand characteristics ---
        # total leaf area index [m2 m-2]
        self.LAI = sum([pt.LAI for pt in self.planttypes])
        # total leaf area density [m2 m-3]
        self.lad = sum([pt.lad for pt in self.planttypes])
        
         # layerwise mean leaf characteristic dimension [m] for interception model
        self.leaf_length = sum([pt.leafp['lt'] * pt.lad for pt in self.planttypes]) / (self.lad + EPS)
        
        # below does not work if root zone depth is different among planttypes
        # root area density
        rad = sum([pt.Roots.rad for pt in self.planttypes])  # total fine root density [m2 m-3]
        self.rad = rad / sum(rad)  # normalized total fine root density distribution [-]
        
#        # workaround for multi-species in multi-layer soil
#        rad = np.zeros(np.shape(dz_soil))
#
#        for pt in self.planttypes:
#            ix = np.where(pt.Roots.rad > 0)[0]
#            rad[ix] += pt.Roots.rad[ix]
#        rad = rad[rad > 0]
#        self.rad = rad / sum(rad)  # normalized total fine root density distribution [-]
        
        # canopy height [m]
        if len(np.where(self.lad > 0)[0]) > 0:
            f = np.where(self.lad > 0)[0][-1]
            self.hc = self.z[f].copy()
        else:
            self.hc = 0.0

        # --- radiation, micromet, interception, and forestfloor instances
        self.radiation = Radiation(cpara['radiation'], self.Switch_Ebal)

        self.micromet = Micromet(self.z, self.lad, self.hc, cpara['flow'])

        self.interception = Interception(cpara['interception'], self.lad * self.dz)
        
        # --- ground surface for testing ---
        self.ground_surface = cpara['ground']
        
        #self.forestfloor = ForestFloor(cpara['forestfloor'])

    def run_daily(self, doy, Ta, PsiL=0.0, Rew=1.0):
        r""" Computatations at daily timestep.
        Updates planttypes and total canopy leaf area index and phenological state.
        Recomputes normalize flow statistics with new leaf area density profile.

        Args:
            doy (float): day of year [days]
            Ta (float): mean daily air temperature [degC]
            PsiL (float): leaf water potential [MPa] --- CHECK??
            Rew (float): relatively extractable water (-)
        """

        """ update physiology and leaf area of planttypes and canopy"""
        for pt in self.planttypes:
            pt.update_daily(doy, Ta, PsiL=PsiL, Rew=Rew)  # updates pt properties

        # total leaf area index [m2 m-2]
        self.LAI = sum([pt.LAI for pt in self.planttypes])
        # total leaf area density [m2 m-3]
        self.lad = sum([pt.lad for pt in self.planttypes])
         # layerwise mean leaf characteristic dimension [m]
        self.leaf_length = sum([pt.leafp['lt'] * pt.lad for pt in self.planttypes]) / (self.lad + EPS)

        #""" normalized flow statistics in canopy with new lad """
        #if self.Switch_Eflow and self.planttypes[0].Switch_lai:
        #    self.micromet.normalized_flow_stats(self.z, self.lad, self.hc)

    def run_timestep(self, dt, z_asl, forcing, parameters):
        r""" Calculates one timestep and updates state of CanopyModel object.

        Args:
            dt: timestep [s]
            z_asl: surface layer height [m]
            forcing (dataframe): meteorological and soil forcing data  !! NOT UP TO DATE
                'precipitation': precipitation rate [m s-1]
                'dir_par': direct fotosynthetically active radiation [W m-2]
                'dif_par': diffuse fotosynthetically active radiation [W m-2]
                'dir_nir': direct near infrared radiation [W m-2]
                'dif_nir': diffuse near infrare active radiation [W m-2]
                'lw_in': Downwelling long wave radiation [W m-2]
                'air_temperature': air temperature at z_asl [\ :math:`^{\circ}`\ C]
                'co2': CO2 mixing ratio at z_asl [ppm]
                'h2o': ambient H2O mixing ratio at z_asl [mol mol-1]
                'wind_speed': mean wind speed at z_asl [m s-1]
                'air_pressure': pressure [Pa]
                'zenith_angle': solar zenith angle [rad]
                'soil_temperature': [\ :math:`^{\circ}`\ C] properties of first soil node
                'soil_water_potential': [m] properties of first soil node
                'soil_volumetric_water': [m m\ :sup:`-3`\ ] properties of first soil node
            'parameters':
                'date'
                'thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ] properties of first soil node
                'hydraulic_conductivity': [m s\ :sup:`-1`\ ] properties of first soil node
                'depth': [m] properties of first soil node

        Returns:
            fluxes (dict)
            states (dict)
        """
        logger = logging.getLogger(__name__)

        """ --- Canopy - ASL flow statistics --- """
        # new grid
        zf = np.arange(0, z_asl + self.dz, self.dz)

        # computes Km for whole profile and returns U, ustar for canopy model grid
        U, ustar = self.micromet.flow_stats(zf, self.lad, self.hc, forcing['wind_speed'])
        #print('ust0', ustar[-1], 'Uo', U[-1], 'Km', self.micromet.Km[-1])
        del zf

        """ --- SW profiles within canopy --- """
        #ff_albedo = self.forestfloor.shortwave_albedo()

        radiation_profiles = {}

        radiation_params = {
            'ff_albedo': self.ground_surface['albedo'],
            'LAIz': self.lad * self.dz,
        }

        # --- PAR ---
        radiation_params.update({
            'radiation_type': 'par',
        })

        radiation_profiles['par'] = self.radiation.shortwave_profiles(
            forcing=forcing,
            parameters=radiation_params
        )

        f_sl = radiation_profiles['par']['sunlit']['fraction']
        sunlit_fraction = radiation_profiles['par']['sunlit']['fraction']

        if self.Switch_Ebal:
            # --- NIR ---
            radiation_params['radiation_type'] = 'nir'

            radiation_profiles['nir'] = self.radiation.shortwave_profiles(
                forcing=forcing,
                parameters=radiation_params
            )
    
            # absorbed radiation by leafs [W m-2(leaf)]
            radiation_profiles['sw_absorbed'] = (
                radiation_profiles['par']['sunlit']['absorbed'] * sunlit_fraction
                + radiation_profiles['nir']['sunlit']['absorbed'] * sunlit_fraction
                + radiation_profiles['par']['shaded']['absorbed'] * (1. - sunlit_fraction)
                + radiation_profiles['nir']['shaded']['absorbed'] * (1. - sunlit_fraction)
            )
 
        """ --- start iterative solution of H2O, CO2, T, Tleaf and Tsurf --- """

        max_err = 0.01  # maximum relative error
        max_iter = 25  # maximum iterations
        gam = 0.2  # weight for new value in iterations
        #err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
        err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 0.
        Switch_WMA = self.Switch_WMA

        # initialize state variables within canopy layer
        T, H2O, CO2, Tleaf = self._restore(forcing)
        sources = {
            'h2o': None,  # [mol m-3 s-1]
            'co2': None,  # [umol m-3 s-1]
            'sensible_heat': None,  # [W m-3]
            'latent_heat': None,  # [W m-3]
            'fr': None  # [W m-3]
        }

        iter_no = 0
        while (err_t > max_err or err_h2o > max_err or
               err_co2 > max_err or err_Tl > max_err or
               err_Ts > max_err) and iter_no <= max_iter:

            iter_no += 1
            Tleaf_prev = Tleaf.copy()
            #Tsurf_prev = self.ground_surface.temperature

            if self.Switch_Ebal:
                
                # ---  LW profiles within canopy ---
                # assumes pre-set state of surface below plant canopy
                lw_up = self.ground_surface['emissivity'] * STEFAN_BOLTZMANN * \
                    np.power((self.ground_surface['temperature'] + DEG_TO_KELVIN), 4.0)
                
                lw_forcing = {
                    'lw_in': forcing['lw_in'],
                    'lw_up': lw_up,
                    'leaf_temperature': Tleaf_prev,
                }

                #print('lwup',lw_forcing['lw_up'], 'lwin',lw_forcing['lw_in'])
                lw_params = {
                    'LAIz': self.lad * self.dz,
                    'ff_emissivity': self.ground_surface['emissivity']
                }

                radiation_profiles['lw'] = self.radiation.longwave_profiles(
                    forcing=lw_forcing,
                    parameters=lw_params
                )

            # --- heat, h2o and co2 source terms
            for key in sources.keys():
                sources[key] = 0.0 * self.ones

            # --- wet leaf water and energy balance ---
            interception_forcing = {
                'h2o': H2O,
                'wind_speed': U,
                'air_temperature': T,
                'air_pressure': forcing['air_pressure'],
                'leaf_temperature': Tleaf_prev,
                'precipitation': forcing['precipitation'],
            }

            if self.Switch_Ebal:
                interception_forcing.update({
                    'sw_absorbed': radiation_profiles['sw_absorbed'],
                    'lw_radiative_conductance': radiation_profiles['lw']['radiative_conductance'],
                    'net_lw_leaf': radiation_profiles['lw']['net_leaf'],
                })

            interception_params = {
                'LAIz': self.lad * self.dz,
                'leaf_length': self.leaf_length
            }

            interception_controls = {
                'energy_balance': self.Switch_Ebal,
                'logger_info': 'date: {} iteration: {}'.format(
                    parameters['date'],
                    iter_no
                )
            }

            if self.Switch_Ebal:
                interception_forcing.update({
                    'sw_absorbed': radiation_profiles['sw_absorbed'],
                    'lw_radiative_conductance': radiation_profiles['lw']['radiative_conductance'],
                    'net_lw_leaf': radiation_profiles['lw']['net_leaf'],
                })

            wetleaf_fluxes = self.interception.run(
                dt=dt,
                forcing=interception_forcing,
                parameters=interception_params,
                controls=interception_controls
            )

            # dry leaf fraction
            df = self.interception.df

            # update source terms
            for key in wetleaf_fluxes['sources'].keys():
                sources[key] += wetleaf_fluxes['sources'][key] / self.dz

            # canopy layer leaf temperature
            Tleaf = self.interception.Tl_wet * (1 - df) * self.lad

            # --- dry leaf gas-exchange ---
            pt_stats = []
            for pt in self.planttypes:

                forcing_pt = {
                    'h2o': H2O,
                    'co2': CO2,
                    'air_temperature': T,
                    'air_pressure': forcing['air_pressure'],
                    'wind_speed': U,
                    'par': radiation_profiles['par'],
                    'leaf_temperature': Tleaf_prev,
                    'nir': radiation_profiles['nir'],
                    'lw': radiation_profiles['lw']
                    }

                parameters_pt = {
                    'dry_leaf_fraction': self.interception.df,
                    'sunlit_fraction': sunlit_fraction
                }

                controls_pt = {
                    'energy_balance': self.Switch_Ebal,
                    'logger_info': 'date: {} iteration: {}'.format(
                        parameters['date'],
                        iter_no
                    )
                }

                if self.Switch_Ebal:
                    forcing_pt.update({
                        'nir': radiation_profiles['nir'],
                        'lw': radiation_profiles['lw'],
                    })

                # --- sunlit and shaded leaves
                pt_stats_i, pt_sources = pt.run(
                    forcing=forcing_pt,
                    parameters=parameters_pt,
                    controls=controls_pt
                )

                # update source terms
                # Dictionary IS modiefied in a loop. Wrapped in a list.
                for key in pt_sources.keys():
                    sources[key] += pt_sources[key]

                # append results
                pt_stats.append(pt_stats_i)

                # canopy leaf temperature
                Tleaf += pt_stats_i['Tleaf'] * df * pt.lad

            # canopy leaf temperature as weighted average
            Tleaf = Tleaf / (self.lad + EPS)

            err_Tl = max(abs(Tleaf - Tleaf_prev))

            
            """ --- solve ground surface --- """
            
            H_gr = 0.0 # Wm-2
            LE_gr = 0.0 # Wm-2
            Fc_gr = 0.0 # umolm-s2s-1
            
            fluxes_ffloor = {'T': self.ground_surface['temperature'], 'H': H_gr, 'LE': LE_gr, 'Fc': Fc_gr}
            
            """  --- solve scalar profiles (H2O, CO2, T) in Canopy and ASL --- """
            
            if Switch_WMA is False:
                # to recognize oscillation
                if iter_no > 1:
                    T_prev2 = T_prev.copy()
                T_prev = T.copy()
                
                # whole profiles in _asl_profs (dict), those of canopy layer in H2O, CO2, T
                H2O, CO2, T, err_h2o, err_co2, err_t, asl_profs = self.micromet.scalar_profiles(
                        gam, H2O, CO2, T, forcing['air_pressure'],
                        source=sources,
                        lbc={'H2O': LE_gr / LATENT_HEAT,
                             'CO2': Fc_gr,
                             'T': H_gr},
                        ubc={'H2O': forcing['h2o'], # to match concentrations in growing asl
                             'CO2': forcing['co2'],
                             'T': forcing['air_temperature']},
                        Ebal=self.Switch_Ebal)
                
                # to recognize oscillation
                if iter_no > 5 and np.mean((T_prev - T)**2) > np.mean((T_prev2 - T)**2):
                    T = 0.5*(T_prev + T) 
                    gam = max(gam / 2, 0.25)

                if (iter_no == max_iter or any(np.isnan(T)) or
                    any(np.isnan(H2O)) or any(np.isnan(CO2))):

                    if (any(np.isnan(T)) or any(np.isnan(H2O)) or any(np.isnan(CO2))):
                        logger.debug('%s Solution of profiles blowing up, T nan %s, H2O nan %s, CO2 nan %s',
                                         parameters['date'],
                                         any(np.isnan(T)), any(np.isnan(H2O)), any(np.isnan(CO2)))
                    elif max(err_t, err_h2o, err_co2, err_Tl, err_Ts) < 0.05:
                        if max(err_t, err_h2o, err_co2, err_Tl, err_Ts) > 0.01:
                            logger.debug('%s Maximum iterations reached but error tolerable < 0.05',
                                         parameters['date'])
                        break
                    
                    Switch_WMA = True  # if no convergence, re-compute with WMA -assumption

                    logger.debug('%s Switched to WMA assumption: err_T %.4f, err_H2O %.4f, err_CO2 %.4f, err_Tl %.4f, err_Ts %.4f',
                                 parameters['date'],
                                 err_t, err_h2o, err_co2, err_Tl, err_Ts)
                    #logger.debug('%s Max err_T at %d', parameters['date'], np.where(abs(T - T_prev) == err_t)[0])
                    
                    # reset values
                    iter_no = 0
                    err_t, err_h2o, err_co2, err_Tl, err_Ts = 999., 999., 999., 999., 999.
                    T, H2O, CO2, Tleaf = self._restore(forcing)
            else:
                err_h2o, err_co2, err_t = 0.0, 0.0, 0.0

        """ --- update state variables --- """
        self.interception.update()
        
        #self.forestfloor.update()
        asl_profs.update({'U': self.micromet.U, 'ust': self.micromet.ust})
        
        """ ---  integrate to ecosystem fluxes (per m-2 ground) --- """

        flux_co2 = (np.cumsum(sources['co2']) * self.dz
                    + Fc_gr)  # [umol m-2 s-1]
        flux_latent_heat = (np.cumsum(sources['latent_heat']) * self.dz
                            + LE_gr)  # [W m-2]
        flux_sensible_heat = (np.cumsum(sources['sensible_heat']) * self.dz
                              + H_gr)  # [W m-2]

        # net ecosystem exchange [umol m-2 s-1]
        NEE = flux_co2[-1]
        # ecosystem respiration [umol m-2 s-1]
        Reco = sum([pt_st['dark_respiration'] for pt_st in pt_stats]) +  Fc_gr#fluxes_ffloor['respiration']
        # ecosystem GPP [umol m-2 s-1]
        GPP = - NEE + Reco
        # stand transpiration [m s-1]
        Tr = sum([pt_st['transpiration'] * MOLAR_MASS_H2O * 1e-3 for pt_st in pt_stats])

        if self.Switch_Ebal:
            # energy closure of canopy  -- THIS IS EQUAL TO frsource (the error caused by linearizing sigma*ef*T^4)
            energy_closure =  sum((radiation_profiles['sw_absorbed'] +
                                   radiation_profiles['lw']['net_leaf']) * self.lad * self.dz) - (  # absorbed radiation
                              sum(sources['sensible_heat'] * self.dz)  # sensible heat
                              + sum(sources['latent_heat'] * self.dz))  # latent heat

        """ --- RESULTS FROM TIMESTEP --- """

#        fluxes_ffloor.update({
#                'potential_infiltration': fluxes_ffloor['potential_infiltration'],
#                'evaporation_bryo': fluxes_ffloor['bryo_evaporation'] * MOLAR_MASS_H2O * 1e-3,  # [m s-1]
#                'evaporation_litter': fluxes_ffloor['litter_evaporation'] * MOLAR_MASS_H2O * 1e-3,  # [m s-1]
#                'evaporation_soil': fluxes_ffloor['soil_evaporation'] * MOLAR_MASS_H2O * 1e-3,  # [m s-1]
#                'evaporation': fluxes_ffloor['evaporation'] * MOLAR_MASS_H2O * 1e-3  # [m s-1]
#                })

        
        # return state and fluxes in dictionary
        state_canopy = {
                'interception_storage': sum(self.interception.W),
                'LAI': self.LAI,
                'lad': self.lad,
                'sunlit_fraction': f_sl,
                'phenostate': sum([pt.LAI * pt.pheno_state for pt in self.planttypes])/self.LAI,
                'IterWMA': iter_no
                }

        fluxes_canopy = {
                'wind_speed': U,
                'friction_velocity': ustar,
                'throughfall': wetleaf_fluxes['throughfall'],
                'interception': wetleaf_fluxes['interception'],
                'evaporation': wetleaf_fluxes['evaporation'],
                'condensation': wetleaf_fluxes['condensation'],
                'condensation_drip': wetleaf_fluxes['condensation_drip'],
                'evaporation_ml': wetleaf_fluxes['evaporation_ml'],
                'throughfall_ml': wetleaf_fluxes['throughfall_ml'],
                'condensation_drip_ml': wetleaf_fluxes['condensation_drip_ml'],
                'transpiration': Tr,
                'SH': flux_sensible_heat[-1],
                'NEE': NEE,
                'GPP': GPP,
                'TER': Reco,
                'LE': flux_latent_heat[-1],
                'co2_flux': flux_co2,  # [umol m-2 s-1]
                'latent_heat_flux': flux_latent_heat,  # [W m-2]
                'pt_transpiration': np.array([pt_st['transpiration'] * MOLAR_MASS_H2O * 1e-3 for pt_st in pt_stats]),
                'pt_gpp': np.array([pt_st['net_co2'] + pt_st['dark_respiration'] for pt_st in pt_stats]),
                'pt_dark_respiration': np.array([pt_st['dark_respiration'] for pt_st in pt_stats]),
                'pt_stomatal_conductance_h2o':  np.array([pt_st['stomatal_conductance'] for pt_st in pt_stats]),
                'pt_boundary_conductance_h2o':  np.array([pt_st['boundary_conductance'] for pt_st in pt_stats]),
                'pt_leaf_internal_co2':  np.array([pt_st['leaf_internal_co2'] for pt_st in pt_stats]),
                'pt_leaf_surface_co2':  np.array([pt_st['leaf_surface_co2'] for pt_st in pt_stats]),
                'water_closure': wetleaf_fluxes['water_closure'],
                }

        if self.Switch_WMA is False:
            state_canopy.update({'h2o': H2O,
                          'co2': CO2,
                          'temperature': T,
                          'WMA_assumption': 1.0*Switch_WMA})

        if self.Switch_Ebal:
            # layer - averaged leaf temperatures are averaged over plant-types
            Tleaf_sl = np.where(self.lad > 0.0,
                                sum([pt_st['Tleaf_sl'] for pt_st in pt_stats]) / (self.lad + EPS),
                                np.nan)
            Tleaf_sh = np.where(self.lad > 0.0,
                                sum([pt_st['Tleaf_sh'] for pt_st in pt_stats]) / (self.lad + EPS),
                                np.nan)
            Tleaf_wet = np.where(self.lad > 0.0,
                                 self.interception.Tl_wet,
                                 np.nan)
            SWnet = (radiation_profiles['nir']['down'][-1] - radiation_profiles['nir']['up'][-1] +
                     radiation_profiles['par']['down'][-1] - radiation_profiles['par']['up'][-1])
            LWnet = (radiation_profiles['lw']['down'][-1] - radiation_profiles['lw']['up'][-1])

            state_canopy.update({
                    'Tleaf_wet': Tleaf_wet,
                    'Tleaf_sl': Tleaf_sl,
                    'Tleaf_sh': Tleaf_sh,
                    'Tleaf': np.where(self.lad > 0.0, Tleaf, np.nan)
                    })

            fluxes_canopy.update({
                    'leaf_SW_absorbed': radiation_profiles['sw_absorbed'],
                    'leaf_net_LW': radiation_profiles['lw']['net_leaf'],
                    'sensible_heat_flux': flux_sensible_heat,  # [W m-2]
                    'energy_closure': energy_closure,
                    'SWnet': SWnet,
                    'LWnet': LWnet,
                    'PARdn': radiation_profiles['par']['down'],
                    'PARup': radiation_profiles['par']['up'],
                    'NIRdn': radiation_profiles['nir']['down'],
                    'NIRup': radiation_profiles['nir']['up'],
                    'LWdn':  radiation_profiles['lw']['down'],
                    'LWup':  radiation_profiles['lw']['up'],
                    'fr_source': sum(sources['fr'] * self.dz)})

        return fluxes_canopy, state_canopy, fluxes_ffloor, asl_profs

    def _restore(self, forcing):
        """ initialize state variables """

        T = self.ones * ([forcing['air_temperature']])
        H2O = self.ones * ([forcing['h2o']])
        CO2 = self.ones * ([forcing['co2']])
        Tleaf = T.copy() * self.lad / (self.lad + EPS)
        #self.forestfloor.restore()
        self.interception.Tl_wet = T.copy()
        for pt in self.planttypes:
            pt.Tl_sh = T.copy()
            pt.Tl_sl = T.copy()

        return T, H2O, CO2, Tleaf

# EOF
