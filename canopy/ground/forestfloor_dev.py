#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: forestfloor
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho

Describes forest floor consisting of bryophytes, baresoil and snowpack.

canopy.forestfloor is interface that handles inputs and outputs from 
Bryophyte, Baresoil, Litter and Snowpack -classes that can co-exist at
forest floor. Allows many Bryophyte -types.

Note:
    migrated to python3
    - absolute imports
    - dictionary.items() in for-loops

Note: roughness length z0 is [1/15 - 1/30] * x where x is height of element.
        given at the moment as a parameter, but for mosses this is problematic.
        There should be species specific 'effective roughness' (and also
        something that takes into account roughness of foresfloor)

DEVELOPMENT VERSION OF SAMULI 9.1.2020
    
Created on Tue Mar 13 11:59:23 2018
"""

from canopy.constants import EPS

from .organiclayer import OrganicLayer
from .snowpack import DegreeDaySnow
from .carbon import SoilRespiration

import logging
logger = logging.getLogger(__name__)

class ForestFloor(object):
    r"""Describes forest floor consisting of organic layers and possible snowpack.
    """

    def __init__(self, para, respiration_profile=None):
        r""" Initializes forestfloor object 
        Args
            para (dict):
             org_types (dict):
                'name' (str) 
                'layer_type': 'bryophyte' or 'litter'
                'coverage':  [-]
                'height': [m]
                'roughness_height': [m]
                #'leaf_area_index': [m\ :sup:`2` m :sup:`-2`\ ]
                #'specific_leaf_area': [m\ :sup:`3` m :sup:`-3`\ ]
                'dry_mass': [kg m\ :sup:`-2`]
                'bulk_density': [kg m\ :sup:`-3`]
                'max_water_content': [g g\ :sup:`-1`\ ]
                'min_water_content': [g g\ :sup:`-1`\ ]
                'water_retention' (dict):
                    'theta_s': saturated water content [m\ :sup:`3` m :sup:`-3`\ ]
                    'theta_r': residual water content [m\ :sup:`3` m :sup:`-3`\ ]
                    'alpha': air entry suction [cm\ :sup:`-1`]
                    'n': pore size distribution [-]
                    'saturated conductivity': [m s\ :sup:`-1`]
                    'pore connectivity': (l) [-]
                'porosity': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'photosynthesis' (list):
                    0. Amax [\ :math:`\mu`\ mol m\ :sup:`-1`\ :sub:`leaf` s\ :sup:`-1`]
                    1. b [mol mol\ :sup:`-1`]
                'respiration' (list):
                    0. Q10 [-]
                    1. R10 [\ :math:`\mu`\ mol m\ :sup:`-1`\ :sub:`leaf` s\ :sup:`-1`]
                'optical_properties' (dict):
                    'albedo_par': [-] photosynthetically active radiation (PAR)
                    'albedo_nir': [-] near infrared radiation (NIR)
                    'emissivity': [-]
                'initial_conditions' (dict)
            snowpack (dict):
                'kmelt' Melting coefficient [m degC-1 s-1] (=2.0 mm/C/d)
                'kfreeze': Freezing  coefficient [m degC-1 s-1] (=0.5 mm/C/d)
                'retention': max fraction of liquid water in snow [-]
                'Tmelt': temperature when melting starts [degC]
                'optical_properties':
                        'emissivity': 
                        'albedo_PAR':
                        'albedo_NIR':
                'swe_initial': [kg m-2]
            soil_respiration (dict):
                    'r10': base respiration rate [umolm-2s-1]
                    'q10': temperature sensitivity [-]
                    'moisture_coeff' (list): moisture response parameters 
        Returns:
            self (object)
        """
        
        # -- forest floor tiled surface of organic layers. snowpack can overly ground
        self.snowpack = DegreeDaySnow(para['snowpack'])

        self.soilrespiration = SoilRespiration(para['soil_respiration'], weights=respiration_profile)

        # Organiclayers; includes both bryophytes and litter. Append to list and
        # Compute area-weighted forest floor temperature and optical properties
        gtypes = []
        gtnames = list(para['bottom_layer_types'].keys())
        gtnames.sort()

        f_organic = 0.0
        ff_temperature = 0.0
        ff_water_storage = 0.0
        albedo = {'PAR': 0.0, 'NIR': 0.0}
        emissivity = 0.0

        for gt in gtnames:
            if para['bottom_layer_types'][gt]['coverage'] > 0:
                print(gt)
                # case coverage > 0
                gtypes.append(OrganicLayer(para['bottom_layer_types'][gt]))
                
                f_organic += gtypes[-1].coverage
                ff_temperature += gtypes[-1].coverage * gtypes[-1].temperature
                ff_water_storage += gtypes[-1].coverage * gtypes[-1].water_storage
                albedo['PAR'] += gtypes[-1].coverage * gtypes[-1].albedo['PAR']
                albedo['NIR'] += gtypes[-1].coverage * gtypes[-1].albedo['NIR']
                emissivity += gtypes[-1].coverage * gtypes[-1].emissivity
        
        if abs(f_organic - 1.0) > EPS:
            raise ValueError("The sum of organic type coverage "
                             + "should be one! Now %.2f" % f_organic)

        self.groundtypes = gtypes
        
        self.temperature = ff_temperature
        self.water_storage = ff_water_storage
        self.albedo = albedo
        self.emissivity = emissivity
        
        if self.snowpack.swe > 0:
            self.temperature = self.snowpack.temperature
            self.albedo = self.snowpack.optical_properties['albedo']
            self.emissivity = self.snowpack.optical_properties['emissivity']

    def update(self):
        """ Updates forestfloor states
        """
        self.snowpack.update()

        ff_temperature = 0.0
        ff_water_storage = 0.0
        albedo = {'PAR': 0.0, 'NIR': 0.0}
        emissivity = 0.0        

        for gt in self.groundtypes:
            gt.update_state()
            
            ff_temperature += gt.coverage * gt.temperature
            ff_water_storage += gt.coverage * gt.water_storage
            albedo['PAR'] += gt.coverage * gt.albedo['PAR']
            albedo['NIR'] += gt.coverage * gt.albedo['NIR']
            emissivity += gt.coverage * gt.emissivity
        
        if self.snowpack.swe > 0:
            self.temperature = self.snowpack.temperature
            self.albedo = self.snowpack.optical_properties['albedo']
            self.emissivity = self.snowpack.optical_properties['emissivity']
        else:
            self.temperature = ff_temperature
            self.water_storage = ff_water_storage
            self.albedo = albedo
            self.emissivity = emissivity

    def run(self, dt, forcing, parameters, controls):
        r"""Water and energy balance at the forestfloor; handles 'tiled' surfaces
        at forest floor and aggregates average fluxes and state.

        Args:
            dt: timestep [s]
            forcing (dict): states of microclimate
                'precipitation_rain': [kg m-2 s\ :sup:`-1`\ ]
                'precipitation_snow': [kg m-2 s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ] if energy_balance is True
                'lw_dn': [W m\ :sup:`-2`\ ] if energy_balance is True
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
                'soil_volumetric_water': [m\ :sup:`3`\  m\`-3`\ ]
                'soil_volumetric_air': [m\ :sup:`3`\  m\`-3`\ ]
                'soil_pond_storage': [kg m-2]
            parameters (dict):
                'soil_thermal_conductivity': [] if energy_balance is True
                'soil_hydraulic_conductivity': []
                'depth': [m] first soil calculation node
                'reference_height': [m] first canopy calculation node
            controls (dict):
                'energy_balance': boolean
                'logger_info': str
        Returns:
            fluxes (dict): forestfloor aggregated fluxes
                'net_radiation' [W m-2]
                'sensible_heat'
                'latent_heat'
                'ground_heat'
                'energy_closure'
                
                'evaporation': tiles + soil below [kg m-2 s-1]
                'soil_evaporation': from soil [kg m-2 s-1]
                'potential_infiltration'
                'capillar_rise'
                'pond_recharge'
                'water_closure'
                
                'co2_flux' [umolm m-2 (ground) s-1]
                'photosynthesis'
                'respiration'
                'soil_respiration' - from belowground
                
            state (dict): forestfloor aggregated state
                'temperature' [degC]
                'water_storage' [kg m-2]
                'snow_water_equivalent' [kg m-2]
            
            gt_outputs (dict of lists): groundtype fluxes
                'net_radiation' [W m\ :sup:`-2`\ ]
                'latent_heat' [W m\ :sup:`-2`\ ]
                'sensible_heat' [W m\ :sup:`-2`\ ]
                'ground_heat' [W m\ :sup:`-2`\ ] (negative towards soil)
                'heat_advection' [W m\ :sup:`-2`\ ]
                'water_closure' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'energy_closure' [W m\ :sup:`-2`\ ]
                'evaporation' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'interception' [kg m\ :sup:`-2`\ s\ :sup`-1`\]            
                'pond_recharge' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'capillary_rise' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                'throughfall' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
                
                'temperature': [\ :math:`^{\circ}`\ C]
                'volumetric_water': [m\ :sup:`3` m\ :sup:`-3`\ ]
                'water_potential': [m]
                'water_content': [g g\ :sup:`-1`\ ]
                'water_storage': [kg m\ :sup:`-2`\ ]
                'hydraulic_conductivity': [m s\ :sup:`-1`\]
                'thermal_conductivity': [W m-1 K-1]
        """
        # initialize fluxes and states
        fluxes = {
            'net_radiation': 0.0, # [W m-2]
            'sensible_heat': 0.0, # [W m-2]
            'latent_heat': 0.0, # [W m-2]
            'ground_heat': 0.0, # [W m-2]
            'energy_closure': 0.0, # [W m-2]

            'evaporation': 0.0,  # [kg m-2 s-1]
            'soil_evaporation': 0.0, # [kg m-2 s-1]            
            'throughfall': 0.0,  # [kg m-2 s-1]
            'capillary_rise': 0.0,  # [kg m-2 s-1]
            'pond_recharge': 0.0, # [kg m-2 s-1]
            'water_closure': 0.0, # [kg m-2 s-1]
            
            'net_co2': 0.0, # [umol m-2(ground) s-1] 
            'photosynthesis': 0.0,  # [umol m-2(ground) s-1]
            'respiration': 0.0,  # [umol m-2(ground) s-1]
            'soil_respiration': 0.0,  # [umol m-2(ground) s-1]
        }

        state = {
            'temperature': 0.0,  # [degC]
            'water_storage': 0.0, # [kg m-2]
            'snow_water_equivalent': 0.0, # [kg m-2]
            # not needed as we take optical properties from previous dt
            #'albedo': None,
            #'emissivity': None
         }
        
        # ground-type specific fluxes and state for output: list of dicts
        gt_results = []

        # --- Soil respiration
        fluxes['soil_respiration'] = self.soilrespiration.respiration(
                                        forcing['soil_temperature'],
                                        forcing['soil_volumetric_water'],
                                        forcing['soil_volumetric_air'])
        
        fluxes['respiration'] += fluxes['soil_respiration']
        fluxes['net_co2'] += fluxes['soil_respiration']

        # --- Snow: now simple degree-day model ---
        snow_forcing = {
            'precipitation_rain': forcing['precipitation_rain'],
            'precipitation_snow': forcing['precipitation_snow'],
            'air_temperature': forcing['air_temperature'],
        }

        fluxes_snow, states_snow = self.snowpack.run(dt=dt, forcing=snow_forcing)

        if self.snowpack.swe > 0:  # snow on the ground
            
            fluxes['throughfall'] += fluxes_snow['potential_infiltration']
            
            # some groundheat flux to keep soil temperatures reasonable
            fluxes['ground_heat'] += (
                parameters['soil_thermal_conductivity']
                / abs(parameters['soil_depth'])
                * (min(forcing['air_temperature'],0.0) - forcing['soil_temperature'][0])
            )
            
            state['snow_water_equivalent'] = states_snow['snow_water_equivalent']
            state['temperature'] = states_snow['temperature']
            
        else:
            # solve organic layer tiles and compute effective forestfloor fluxes
            # and state
            org_forcing = forcing.copy()
            del org_forcing['precipitation_rain'], org_forcing['precipitation_snow']
            
            org_forcing.update(
                    {'precipitation': fluxes_snow['potential_infiltration'],
                    'soil_temperature': forcing['soil_temperature'][0]}
                    )
            
            for gt in self.groundtypes:
                gt_flx, gt_state = gt.run(dt, org_forcing, parameters, controls)
                                
                # effective forest floor fluxes and state
                for key in fluxes.keys():
                    if key in gt_flx.keys():
                        fluxes[key] += gt.coverage * gt_flx[key]
                
                state['temperature'] += gt.coverage * gt_state['temperature']
                state['water_storage'] += gt.coverage * gt_state['water_storage']
                
                # merge dicts and append to gt_results
                gt_flx.update(gt_state)
                gt_results.append(gt_flx)
                del gt_flx, gt_state
                
            fluxes['evaporation'] += fluxes['soil_evaporation']
        
        # groundtype specific results (fluxes & state): convert list of dicts to dict of lists

        gt_outputs = {}
        for k,v in gt_results[0].items():
            gt_outputs[k] = [x[k] for x in gt_results]
            
        return fluxes, state, gt_outputs

# EOF
