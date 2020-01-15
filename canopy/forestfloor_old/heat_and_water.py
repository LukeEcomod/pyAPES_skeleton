#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:23:08 2018

Note:
    migrated to python3
    - absolute import

@author: ajkieloaho
"""

import numpy as np
from scipy.integrate import odeint

from canopy.constants import EPS, MOLAR_MASS_H2O, SPECIFIC_HEAT_H2O
from canopy.constants import SPECIFIC_HEAT_ORGANIC_MATTER, LATENT_HEAT
from canopy.constants import DEG_TO_KELVIN, STEFAN_BOLTZMANN
from canopy.constants import SPECIFIC_HEAT_AIR, WATER_DENSITY, GAS_CONSTANT
from canopy.constants import MOLECULAR_DIFFUSIVITY_CO2, MOLECULAR_DIFFUSIVITY_H2O
from canopy.constants import THERMAL_DIFFUSIVITY_AIR, GRAVITY
from canopy.constants import AIR_VISCOSITY, AIR_DENSITY
from .odesolver import solver_array, ForwardEuler_array
#from canopy.constants import *

import logging
logger = logging.getLogger(__name__)

class heat_and_water_balance:
    """ 
    Defines heat and water balance equation to be solved
    with ode-solvers
    """

    def __init__(self, properties, forcing, parameters):
        r""" Inititalises a heat and water balance calculation object

        Args:
            properties (dict): characteristics of Bryophyte instance

            forcing (dict):
                'throughfall': [kg m\ :sup'-2' s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ]
                'lw_dn': [W m\ :sup:`-2`\ ]
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'precipitation_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [Pa]
            parameters (dict):
                'soil_depth': [m]
                'soil_hydraulic_conductivity'
                'soil_thermal_conductivity'
                'reference_height'
        """

        self.dudt = np.zeros(12)

        self.forcing = forcing
        
        # -- bryotype properties
        self.properties = properties
        
        self.zm = 0.5 * properties['height']        
        
        # water storage limits [kg m-2 (bryophyte)]
        self.min_storage = (
                properties['min_water_content'] * properties['dry_mass'])

        self.max_storage = (
                properties['max_water_content'] * properties['dry_mass'])
        
        # -- extract parameters; do not change during object life
        self.zs = abs(parameters['soil_depth'])
        self.soil_water_potential = forcing['soil_water_potential']
        self.soil_thermal_conductivity = parameters['soil_thermal_conductivity']
        self.soil_hydraulic_conductivity = parameters['soil_thermal_conductivity']
        
        self.zref = parameters['reference_height']

    def __call__(self, y, dt):
        """
        Callable function to solve coupled water and heat budget of bryophyte layer
        Args:
            y (array): initial values: returns derivatives 
                0. temperature [degC] --> computes (du/dt)
                1. water content [kg m-2] --> (du/dt)
                2. pond_recharge
                3. capillary_rise
                4. interception
                5. evaporation/condensation
                6. radiative_heat
                7. sensible_heat
                8. latent_heat 
                9. conducted_heat
                10 advected_heat
                11 energy_closure
            dt (float): time step
        Returns:
            derivatives (du/dt) of y.
        
        Units:  temperature [K s-1]
                water content and water fluxes [kg m-2 s-1 = mm s-1]
                energy fluxes [J m-2 s-1 = W m-2]
        """

        if dt == 0.0:
            dt = dt + EPS

        # --- compute water balance ---
        # nomenclature: 
        #   water_content [g g-1], water_storage [kg m-2 s-1 == mm]
        #   volumetric water content [-], water potential [m]
        #   all fluxes [kg m-2(bryophyte) s-1]

        # initial state
        
        # [kg m-2] or [mm]
        water_storage = min(self.max_storage, y[1])
        water_storage = max(self.min_storage, water_storage)

        # [g g-1]
        water_content = (water_storage / self.properties['dry_mass'])

        # [m m-3]
        volumetric_water = (water_content / WATER_DENSITY
                            * self.properties['bulk_density'])

        # [m]
        water_potential = convert_hydraulic_parameters(
                volumetric_water,
                self.properties['water_retention'],
                'volumetric_water')

        # --- constraints for recharge & evaporation
        max_recharge = max(self.max_storage - water_storage, 0.0)
        max_recharge = min(self.max_storage, max_recharge)

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_evaporation_rate = (y[1] - (self.min_storage + EPS)) / dt

        if np.isinf(max_evaporation_rate) or max_evaporation_rate < 0.0:
            max_evaporation_rate = 0.0

        # [kg m-2 s-1] or [mm s-1]
        max_condensation_rate = -((self.max_storage - EPS) - y[1]) / dt

        if np.isinf(max_condensation_rate) or max_condensation_rate > 0.0:
            max_condensation_rate = 0.0

        # --- compute actual evaporation/condensation rate [kg m-2 s-1] ---

        # boundary layer conductances for H2O, heat and CO2  [mol m-2 s-1]
        dT = y[0] - self.forcing['air_temperature']


        conductance_to_air = surface_atm_conductance(wind_speed=self.forcing['wind_speed'],
                                                     height=self.zref,
                                                     friction_velocity=self.forcing['friction_velocity'],
                                                     dT=dT)

        # Relative conductance is from Williams and Flanagan (1996), Oecologia. 
        relative_conductance = min(1.0,
                                   (0.1285 * y[1]
                                    / self.properties['dry_mass'] - 0.1285)) * 0.5

        # [mol m-2 s-1]
        conductance_to_air_h2o = (
            conductance_to_air['h2o'] * relative_conductance)

        # [kg m-2 s-1]
        evaporation_demand = (conductance_to_air_h2o
                            * (611.0 / self.forcing['air_pressure']
                                * np.exp(17.502 * y[0] / (y[0] + 240.0))
                                - self.forcing['h2o'])
                            * MOLAR_MASS_H2O)

        evaporation_rate = min(evaporation_demand, max_evaporation_rate)
        evaporation_rate = max(evaporation_rate, max_condensation_rate)

        max_recharge_rate = max(max_recharge_rate + evaporation_rate, 0.0)
        max_recharge = max_recharge_rate * dt

        # -- recharge from rainfall interception and/or from pond storage

        interception = min(max_recharge - EPS, self.forcing['throughfall'] * dt)
        
        # [kg m-2 s-1] or [mm s-1]
        interception_rate = interception / dt

        # [kg m-2] or [mm]
        max_recharge = max(max_recharge - interception, 0.0)
        
        pond_recharge = min(max_recharge - EPS, self.forcing['soil_pond_storage'])
        
        # [kg m-2 s-1] or [mm s-1]
        pond_recharge_rate = pond_recharge / dt

        # [kg m-2 s-1] or [mm s-1]
        max_recharge_rate = max(max_recharge - pond_recharge, 0.0) / dt
        
        # --- compute capillary rise from soil [ kg m-2 s-1 = mm s-1]

        # estimate soil-moss hydraulic conductivity assuming two resistors in series
        # this is same Kersti was using for heat trasport
        
        Km = hydraulic_conductivity(water_potential, self.properties['water_retention'])
        Ks = self.soil_hydraulic_conductivity
        
        # conductance of layer [s-1]
        g_moss = Km / self.zm
        g_soil = Ks / self.zs
        
        # [m s-1]
        Kh = (g_moss * g_soil / (g_moss + g_soil)) * (self.zm + self.zs)
        
        capillary_rise = WATER_DENSITY * max(0.0, - Kh * ((water_potential - self.forcing['soil_water_potential']) 
                                                          / (self.zm + self.zs) + 1.0))

        # [kg m-2 s-1] or [mm s-1]
        capillary_rise = min(capillary_rise, max_recharge_rate)

        # calculate mass balance of water

        # [kg m-2 s-1] or [mm s-1]
        dy_water = (
                interception_rate
                + pond_recharge_rate
                + capillary_rise
                - evaporation_rate
                )
        
        #--- compute energy balance
        
        # radiation balance # [J m-2 s-1] or [W m-2]
        
        # [-]
        albedo = bryophyte_shortwave_albedo(water_content,
                                            self.properties)
        emissivity = self.properties['optical_properties']['emissivity']
        
        # [J m-2 s-1] or [W m-2]
        net_shortwave_radiation = ( 
                self.forcing['par'] * (1.0 - albedo['PAR'])
                + self.forcing['nir'] * (1.0 - albedo['NIR'])
                )
        
    
        net_longwave_radiation = emissivity * (self.forcing['lw_dn']
                                - STEFAN_BOLTZMANN * (y[0] + DEG_TO_KELVIN)**4.0)


        # [J m-2 s-1] or [W m-2]
        net_radiation = net_shortwave_radiation + net_longwave_radiation

        # [J m-2 s-1] or [W m-2]
        sensible_heat_flux = (SPECIFIC_HEAT_AIR
                              * conductance_to_air['heat']
                              * (y[0] - self.forcing['air_temperature']))

        # [J m-2 s-1] or [W m-2]
        latent_heat_flux = LATENT_HEAT / MOLAR_MASS_H2O * (evaporation_rate + EPS)

        # heat conduction between moss and soil
        # [W m-2 K-1]
        moss_thermal_conductivity = thermal_conductivity(volumetric_water)

        # thermal conductance [W m-2 K-1]; assume the layers act as two resistors in series
        g_moss = moss_thermal_conductivity / self.zm
        g_soil = self.soil_thermal_conductivity / self.zs
        
        thermal_conductance = (g_moss * g_soil) / (g_moss + g_soil)
        
        # [J m-2 s-1] or [W m-2]
        heat_conduction = thermal_conductance *(y[0] - self.forcing['soil_temperature'])

        # heat lost or gained with liquid water removing/entering 

        # [J m-2 s-1] or [W m-2]
        heat_advection = SPECIFIC_HEAT_H2O * (
                        interception_rate * self.forcing['air_temperature']
                        + capillary_rise * self.forcing['soil_temperature']
                        + pond_recharge_rate * self.forcing['soil_temperature']
                        )

        # heat capacities

        # [J K-1]
        heat_capacity_old = (SPECIFIC_HEAT_ORGANIC_MATTER
                         * self.properties['dry_mass']
                         + SPECIFIC_HEAT_H2O * y[1])

        heat_capacity_new = (SPECIFIC_HEAT_ORGANIC_MATTER
                             * self.properties['dry_mass']
                             + SPECIFIC_HEAT_H2O * (y[1] + dy_water * dt))

        # calculate new temperature from heat balance

        heat_fluxes = (
                net_radiation
                - sensible_heat_flux
                - latent_heat_flux
                - heat_conduction
                + heat_advection
                )

        new_temperature = (heat_fluxes * dt + heat_capacity_old * y[0]) / heat_capacity_new

        # [K m-2 s-1]
        self.dudt[0] = (new_temperature - y[0]) / dt
        # [kg m-2 s-1] or [mm s-1]
        self.dudt[1] = dy_water

        # water fluxes
        # [kg m-2 s-1] or [mm s-1]
        self.dudt[2] = pond_recharge_rate
        self.dudt[3] = capillary_rise
        self.dudt[4] = interception_rate
        self.dudt[5] = evaporation_rate

        # energy fluxes
        # [J m-2 s-1] or [W m-2]
        self.dudt[6] = net_radiation
        self.dudt[7] = sensible_heat_flux
        self.dudt[8] = latent_heat_flux
        self.dudt[9] = heat_conduction
        self.dudt[10] = heat_advection
        self.dudt[11] = heat_fluxes

        return self.dudt

def heat_and_water_exchange(properties,
                            temperature,
                            water_content,
                            dt,
                            forcing,
                            parameters,
                            solver='forward_euler',
                            nsteps=20,
                            logger_info=''):
    r""" Solves moss or litter layer water and energy balance. Interface for
        ode-solvers.

    Args:
        properties (dict): characteristics of Bryophyte instance
        moss_temperature: [\ :math:`^{\circ}`\ C]
        moss_water_content: [g g\ :sup:`-1`\ ]
        dt: [s], timestep
        forcing (dict):
            'throughfall': [kg m\ :sup'-2' s\ :sup:`-1`\ ]
            'par': [W m\ :sup:`-2`\ ]
            'nir': [W m\ :sup:`-2`\ ]
            'lw_dn': [W m\ :sup:`-2`\ ]
            'h2o': [mol mol\ :sup:`-1`\ ]
            'air_temperature': [\ :math:`^{\circ}`\ C]
            'precipitation_temperature': [\ :math:`^{\circ}`\ C]
            'air_pressure': [Pa]
            #'soil_depth': [m]
            'soil_temperature': [\ :math:`^{\circ}`\ C]
            'soil_water_potential': [Pa]
        parameters (dict):
            'soil_depth': [m]
            'soil_hydraulic_conductivity'
            'soil_thermal_conductivity'
            'reference_height' [m]
        solver (str):
            'forward_euler': default
            'odeint': scipy.odeint
        nsteps (int): number of steps in odesolver
        logger_info: str

    Returns:
        fluxes (dict):
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
            'capillar_rise' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
            'throughfall' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
            
        states (dict):
            'temperature': [\ :math:`^{\circ}`\ C]
            'volumetric_water': [m\ :sup:`3` m\ :sup:`-3`\ ]
            'water_potential': [m]
            'water_content': [g g\ :sup:`-1`\ ]
            'water_storage': [kg m\ :sup:`-2`\ ]
            'hydraulic_conductivity': [m s\ :sup:`-1`\]
            'thermal_conductivity': [Wm-1K-1]
    """
    #print(properties.keys(), forcing.keys(), parameters.keys())
    
    solver = solver.lower()

    # initializing state equation object for ode-solver    
    state_equation = heat_and_water_balance(properties, forcing, parameters)

    # initial values for ode-solver
    initial_condition = np.zeros(12)
    initial_condition[0] = temperature
    initial_condition[1] = water_content * properties['dry_mass']

    # solver output timesteps
    tt = np.linspace(0.0, dt, nsteps)

    if solver == 'forward_euler':

        new_states, time_points = solver_array(
                func=state_equation,
                init=initial_condition,
                time=tt,
                method=ForwardEuler_array)
        
    elif solver == 'odeint':

        new_states = odeint(
                func=state_equation,
                y0=initial_condition,
                t=tt,
                full_output=False,
                mxstep=3600)

    del tt

    if np.any(np.isnan(new_states)):
        logger.debug(logger_info +' T: {}, W: {}, evap: {}, h2o air: {}'.format(
                    new_states[0][-1],
                    new_states[1][-1],
                    new_states[5][-1],
                    forcing['h2o']
                    ))

    # new states at every solver timestep
    # [deg C]
    temperature = new_states[:, 0]
    # [kg m-2]
    water_storage = new_states[:, 1]

    # [J m-2]: heat storage = heat capacity * temperature
    heat_storage = ((SPECIFIC_HEAT_ORGANIC_MATTER * properties['dry_mass']
                     + SPECIFIC_HEAT_H2O * water_storage) * temperature)


    # heat storage change during dt [W m-2] or [J m-2 s-1]
    heat_storage_change = (heat_storage[-1] - heat_storage[0]) / dt

    # water storage change [kg m-2 s-1] or [mm s-1]
    water_storage_change = (water_storage[-1] - water_storage[0]) / dt
    
    # now we get only last timepoint values
    temperature = temperature[-1]
    water_storage = water_storage[-1]
    
    # energy fluxes & energy closure
    net_radiation = new_states[-1, 6] / dt
    sensible_heat_flux = new_states[-1, 7] / dt
    latent_heat_flux = new_states[-1, 8] / dt
    ground_heat_flux = new_states[-1, 9] / dt
    heat_advection = new_states[-1, 10] / dt
    
    # [W m-2] or [J m-2 s-1]
    heat_fluxes = (net_radiation
                   - sensible_heat_flux
                   - latent_heat_flux
                   - ground_heat_flux
                   + heat_advection)

    energy_closure = -heat_storage_change + heat_fluxes

    # water fluxes and water balance closure
    # [kg m-2 s-1] or [mm s-1]
    pond_recharge_rate = new_states[-1, 2] / dt
    capillary_rise = new_states[-1, 3] / dt
    interception_rate = new_states[-1, 4] / dt
    evaporation_rate = new_states[-1, 5] / dt

    # [kg m-2 s-1] or [mm s-1]
    water_fluxes = (pond_recharge_rate
                    + capillary_rise
                    - evaporation_rate
                    + interception_rate)

    water_closure = -water_storage_change + water_fluxes


    # [kg m-2 s-1)] of [mm s-1]
    throughfall_rate = (forcing['throughfall'] - interception_rate)

    # --- State variables of bryophyte layer
    # [g g-1]
    water_content = (water_storage / properties['dry_mass'])

    # [m3 m-3]
    volumetric_water = (water_content / WATER_DENSITY
                        * properties['bulk_density'])
    # [m]
    matrix_potential = convert_hydraulic_parameters(volumetric_water,
                                                    properties['water_retention'],
                                                    'volumetric_water')
    
    # [m s-1]
    Kliq = hydraulic_conductivity(matrix_potential, properties['water_retention'])

    # [W m-1 K-1]
    Lambda = thermal_conductivity(volumetric_water)

    # return fluxes and state variables
    fluxes = {
        'net_radiation': net_radiation,  # [W m-2]
        'latent_heat': latent_heat_flux,  # [W m-2]
        'sensible_heat': sensible_heat_flux,  # [W m-2]
        'ground_heat': ground_heat_flux,  # [W m-2]
        'heat_advection': heat_advection,  # [W m-2]
        'water_closure': water_closure,  # [mm s-1]
        'energy_closure': energy_closure,  # [W m-2]
        'interception': interception_rate,  # [mm s-1]
        'evaporation': evaporation_rate,  # [mm s-1]
        'pond_recharge': pond_recharge_rate,  # [mm s-1]
        'capillar_rise': capillary_rise,  # [mm s-1]
        'throughfall': throughfall_rate,  # [mm s-1]
        }

    states = {
        'volumetric_water': volumetric_water,  # [m3 m-3]
        'water_potential': matrix_potential,  # [m]
        'water_content': water_content,  # [g g-1]
        'water_storage': water_storage,  # [kg m-2] or [mm]
        'temperature': temperature,  # [degC]
        'hydraulic_conductivity': Kliq,  # [m s-1]
        'thermal_conductivity': Lambda,  # [W m-1 K-1]
        }

    return fluxes, states

def water_exchange(dt,
                   water_storage,
                   properties,
                   forcing,
                   parameters):
    """
    Args:
        properties (dict): characteristics of Bryophyte instance
        moss_temperature: [\ :math:`^{\circ}`\ C]
        moss_water_content: [g g\ :sup:`-1`\ ]
        dt: [s], timestep
        forcing (dict):
            'throughfall': [kg m\ :sup'-2' s\ :sup:`-1`\ ]
            'par': [W m\ :sup:`-2`\ ]
            'nir': [W m\ :sup:`-2`\ ]
            'lw_dn': [W m\ :sup:`-2`\ ]
            'h2o': [mol mol\ :sup:`-1`\ ]
            'air_temperature': [\ :math:`^{\circ}`\ C]
            'precipitation_temperature': [\ :math:`^{\circ}`\ C]
            'air_pressure': [Pa]
            #'soil_depth': [m]
            'soil_temperature': [\ :math:`^{\circ}`\ C]
            'soil_water_potential': [Pa]
        parameters (dict):
            'soil_depth': [m]
            'soil_hydraulic_conductivity'
            'soil_thermal_conductivity'
            'reference_height' [m]
        solver (str):
            'forward_euler': default
            'odeint': scipy.odeint
        nsteps (int): number of steps in odesolver
        logger_info: str

    Returns:
        fluxes (dict):
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
            'capillar_rise' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
            'throughfall' [kg m\ :sup:`-2`\ s\ :sup`-1`\]
            
        states (dict):
            'temperature': [\ :math:`^{\circ}`\ C]
            'volumetric_water': [m\ :sup:`3` m\ :sup:`-3`\ ]
            'water_potential': [m]
            'water_content': [g g\ :sup:`-1`\ ]
            'water_storage': [kg m\ :sup:`-2`\ ]
            'hydraulic_conductivity': [m s\ :sup:`-1`\]
            'thermal_conductivity': [Wm-1K-1]
    """
    if dt == 0.0:
        dt = dt + EPS
        
    # initial state
    # moss is assumed to be at air temperature
    temperature = forcing['air_temperature']

    # [g g-1]
    water_content = (water_storage / properties['dry_mass'])

    # [m m-3]
    volumetric_water = (water_content / WATER_DENSITY
                        * properties['bulk_density'])

    # [m]
    water_potential = convert_hydraulic_parameters(
        volumetric_water,
        properties['water_retention'],
        'volumetric_water'
    )

    # [kg m-2] or [mm]
    max_storage = properties['max_water_content'] * properties['dry_mass']
    min_storage = properties['min_water_content'] * properties['dry_mass']

    water_storage = min(max_storage, water_storage)
    water_storage = max(min_storage, water_storage)

    max_recharge = max(max_storage - water_storage, 0.0)
    max_recharge = min(max_storage, max_recharge)

    # [kg m-2 s-1] or [mm s-1]
    max_recharge_rate = max_recharge / dt

    #--- evaporation / condensation from layer [kg m-2 s-1] or [mm s-1]
    max_evaporation_rate = (water_storage - (min_storage + EPS)) / dt

    if np.isinf(max_evaporation_rate) or max_evaporation_rate < 0.0:
        max_evaporation_rate = 0.0

    # [kg m-2 s-1] or [mm s-1]
    max_condensation_rate = -((max_storage - EPS) - water_storage) / dt

    if np.isinf(max_condensation_rate) or max_condensation_rate > 0.0:
        max_condensation_rate = 0.0

    # boundary layer conductances for H2O, heat and CO2
    # [mol m-2 s-1]

    conductance_to_air = surface_atm_conductance(wind_speed=forcing['wind_speed'],
                                             height=parameters['reference_height'],
                                             friction_velocity=forcing['friction_velocity'],
                                             dT=0.0)

    # water vapor conductance from moss to air
    # Relative conductance for H2O is from Williams and Flanagan (1996), Oecologia.
    relative_conductance = min(1.0,
                               (0.1285 * water_storage / properties['dry_mass'] - 0.1285)
                              )

    # [mol m-2 s-1]
    conductance_to_air_h2o = (
        conductance_to_air['h2o'] * relative_conductance
    )

    # [kg m-2 s-1]
    evaporation_rate = (
        conductance_to_air_h2o
        * (611.0 / forcing['air_pressure']
           * np.exp(17.502 * temperature / (forcing['air_temperature'] + 240.0))
           - forcing['h2o'])
        * MOLAR_MASS_H2O
    )

    evaporation_rate = min(evaporation_rate, max_evaporation_rate)
    evaporation_rate = max(evaporation_rate, max_condensation_rate)

    max_recharge_rate = max(max_recharge_rate + evaporation_rate, 0.0)
    max_recharge = max_recharge_rate * dt

    #--- interception of rainfall and recharge from ponding water
    # [mm] or [kg m-2]
    interception = min(max_recharge - EPS, forcing['throughfall'] * dt)

    # [kg m-2 s-1] or [mm s-1]
    interception_rate = interception / dt

    # [kg m-2 s-1] or [mm s-1]
    max_recharge = max(
        max_recharge - interception,
        0.0
    )
    pond_recharge = min(max_recharge - EPS, forcing['soil_pond_storage'])

    # [kg m-2 s-1] or [mm s-1]
    pond_recharge_rate = pond_recharge / dt

    # [kg m-2 s-1] or [mm s-1]
    max_recharge_rate = max(
        max_recharge - pond_recharge,
        0.0
    ) / dt


    #--- Capillary rise from underlying soil [kg m-2 s-1] or [mm s-1]

    # hydraulic conductivity from soil to moss [m s-1]
    Km = hydraulic_conductivity(water_potential, properties['water_retention'])
    Ks = parameters['soil_hydraulic_conductivity']
    
    zm = properties['height']
    zs = abs(parameters['soil_depth'])
   
    # conductance of layer [s-1]
    g_moss = Km / zm
    g_soil = Ks / zs
    
    # [m s-1]
    Kh = (g_moss * g_soil / (g_moss + g_soil)) * (zm + zs)
    
    #[kg m-2 s-1]
    capillary_rise = WATER_DENSITY * max(0.0, - Kh * ((water_potential - forcing['soil_water_potential']) 
                                            / (zm + zs) + 1.0))  
    
    capillary_rise = min(capillary_rise, max_recharge_rate)

    #-- calculate water mass balance and values after dt

    # [kg m-2 s-1] or [mm s-1]
    dy_water = (
        interception_rate
        + pond_recharge_rate
        + capillary_rise
        - evaporation_rate
    )

    # [kg m-2] or [mm]
    water_storage += dy_water * dt

    # [g g-1]
    water_content = water_storage / properties['dry_mass']

    # [m3 m-3]
    volumetric_water = (water_content / WATER_DENSITY * properties['bulk_density']
    )

    water_potential = convert_hydraulic_parameters(
        volumetric_water,
        properties['water_retention'],
        'volumetric_water'
    )

    Kliq = hydraulic_conductivity(water_potential, properties['water_retention'])

    # --- Heat exchange dummy variables ---

    # heat conduction between moss and soil
    # [W m-2 K-1]
    moss_thermal_conductivity = thermal_conductivity(volumetric_water)

    # thermal conductance [W m-2 K-1]; assume the layers act as two resistors in series
    g_moss = moss_thermal_conductivity / zm
    g_soil = parameters['soil_thermal_conductivity'] / zs
    
    thermal_conductance = (g_moss * g_soil) / (g_moss + g_soil)
    
    # [J m-2 s-1] or [W m-2]
    ground_heat = thermal_conductance *(forcing['air_temperature'] - forcing['soil_temperature'])

    latent_heat = (LATENT_HEAT * (evaporation_rate + EPS) / MOLAR_MASS_H2O)

    sensible_heat = 0.0 # moss at air temperature
    #    sensible_heat = (
    #        SPECIFIC_HEAT_AIR
    #        * conductance_to_air['heat']
    #        * (forcing['air_temperature'] - temperature)
    #    )

    fluxes = {
        'evaporation': evaporation_rate,  # [mm s-1]
        'capillar_rise': capillary_rise,  # [mm s-1]
        'pond_recharge': pond_recharge_rate,  # [mm s-1]
        'throughfall': forcing['throughfall'] - interception_rate,  # [mm s-1]
        'ground_heat': ground_heat,  # [W m-2]
        'latent_heat': latent_heat,  # [W m-2]
        'sensible_heat': sensible_heat,  # [W m-2]
    }

    states = {
        'volumetric_water': volumetric_water,  # [m3 m-3]
        'water_potential': water_potential,  # [m]
        'water_content': water_content,  # [g g-1]
        'water_storage': water_storage,  # [kg m-2] or [mm]
        'hydraulic_conductivity': Kliq,  # [m s-1]
        'thermal_conductivity': moss_thermal_conductivity,  # [W m-1 K-1]
        'temperature': temperature  # [degC]
    }

    return fluxes, states

def saturation_vapor_pressure(temperature):
    r""" Calculates saturation vapor pressure over water surface.

    Args:
        temperature: [\ :math:`^{\circ}`\ C]

    Returns:
        float: saturation vapor pressure in [Pa]
    """

    # [Pa]
    return 611.0 * np.exp((17.502 * temperature) / (temperature + 240.97))


def vapor_conductance_porous_media(volumetric_water,
                                   porosity,
                                   depth,
                                   temperature,
                                   ambient_pressure=101300.0):
    r""" Estimates molecular conductance of CO\ :sub:`2` and H\ :sub:`2`\ O in
    porous media.

    Following asumptions are made:
        1. all gas exchange is due to molecular diffusion
        2. relative diffusvity in soil is to that of air (D/D\ :sub:`o`\ ) as
           described by Millington and Quirk (1961)

    Args:
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        porosity: [-]
        depth: Transport distance [m]
        ambient_pressure: [Pa]

    Returns:
        dictionary:
            molecular conductances for
                * 'CO2': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'H2O': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
    """

    # [K]
    temperature = temperature + DEG_TO_KELVIN

    # [m3/m3], air filled porosity
    afp = np.maximum(0.0, porosity - volumetric_water)
    # afp = max(0, porosity - volumetric_water)

    # [mol m-3], air molar density
    cair = ambient_pressure / (GAS_CONSTANT * temperature)

    # D/Do, diffusivity in porous media relative to that in free air,
    # Millington and Quirk (1961)
    relative_diffusivity = (np.power((temperature / DEG_TO_KELVIN), 1.75)
                            * np.power(afp, 10.0/3.0) / porosity**2)

    # [mol/(m2 s)], through moss depth: units: mol m-3 * m2 s-1 * m-1
    conductance_h2o = cair * MOLECULAR_DIFFUSIVITY_H2O * relative_diffusivity / depth
    conductance_co2 = cair * MOLECULAR_DIFFUSIVITY_CO2 * relative_diffusivity / depth

    return {
        'co2': conductance_co2,
        'h2o': conductance_h2o
        }


def thermal_conductivity(volumetric_water, method='donnel'):
    r""" Estimates thermal conductivity (km) of bryophyte layer.

    By default organic matter heat conductivity is calculated by using equation
    by O'Donnel et al. (2009, Soil Sci. 174).

    Args:
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        flag (optional):
            optional methods are:
                * 'campbell' (Campbell et al. (1985))
                * 'constant' (0.25)
                * 'lauren' (Lauren (1999, table 4))

    Returns:
        float: heat conductivity in [W m\ :sup:`-1` K\ :sup:`-1`\ ]
    """

    method = method.lower()

    heat_conductivity = None

    if method == 'donnel':  # O'Donnell
        heat_conductivity = np.minimum(
            0.6,
            3.0e-2 + 5.0e-1 * volumetric_water)

        heat_conductivity = np.maximum(4.0e-2, heat_conductivity)

    elif method == 'campbell':
        heat_conductivity = (
            0.4 + 0.5
            * volumetric_water(0.4 - 0.06)
            * np.exp(-(1. * volumetric_water)**4))

    elif method == 'constant':
        heat_conductivity = 0.25

    elif method == 'lauren':
        heat_conductivity = -0.004 + 0.609 * volumetric_water

    # [W/(m K)]
    return heat_conductivity


def soil_boundary_layer_conductance(u, z, zo, Ta, dT, P=101300.):
    """
    Computes soil surface boundary layer conductance (mol m-2 s-1)
    assuming velocity profile logarithmic between z and z0.
    INPUT: u - mean velocity (m/s)
           z - height of mean velocity u (m)
           zo - soil surface roughness length for momentum (m)
           Ta - ambient temperature (degC)
           dT - soil surface-air temperature difference (degC)
           P - pressure(Pa)
    OUTPUT: boundary-layer conductances (mol m-2 s-1)
        gb_h - heat (mol m-2 s-1)
        gb_c- CO2 (mol m-2 s-1)
        gb_v - H2O (mol m-2 s-1)
    Based on Daamen & Simmons (1996). Note: gb decreases both in
    unstable and stable conditions compared to near-neutral;
    nonfeasible?
    Samuli Launiainen, 18.3.2014
    to python by Kersti
    """

    u = np.maximum(u, EPS)

    rho_air = 44.6*(P / 101300.0)*(273.15 / (Ta + 273.13))  # molar density of air [mol/m3]

    delta = 5.0 * GRAVITY * z * dT / ((Ta + 273.15) * u**2)
    if delta > 0:
        d = -0.75
    else:
        d = -2
    rb = (np.log(z/zo))**2 / (0.4**2*u)*(1 + delta)**d

    gb_h = rho_air * 1 / rb
    gb_v = MOLECULAR_DIFFUSIVITY_H2O / THERMAL_DIFFUSIVITY_AIR * gb_h
    gb_c = MOLECULAR_DIFFUSIVITY_CO2 / THERMAL_DIFFUSIVITY_AIR * gb_h

    return gb_h, gb_c, gb_v

def moss_atm_conductance_old(wind_speed, roughness_height):
    r""" Estimates boundary layer conductance of bryophyte canopy.

    Simple version without free convection

    Wind speed should represent vertical wind speed at ca. 20 cm above moss
    canopy (parametrsization derived from wind-tunnel studies). Roughness
    lengths scale is equal to the 'characteristic vertical height of
    individual bryophyte shoots' (typically order of 3-10 mm).

    Estimated CO\ :sub:`2` and heat conductances are related to water vapor
    by the ratio of molecular diffusivities
    (Campbell and Norman, 1998, eq. 7.29-7.33).

    References:
        Rice et al., 2001.
            Significance of variation in bryophyte canopy structure.
            Amer. J. Bot. 88:1568-1576.
        Rice, 2006.
            Towards an integrated undestanding of Bryophyte performance:
            the dimensions of space and time.
            Lindbergia 31:41-53.

    Args:
        wind_speed: [m s\ :sup:`-1`\ ]
        roughness_height: [m]

    Returns:
        dictionary:
            boundary layer conductances for
                * 'co2': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'h2o': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'heat': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
    NOTE:
        Add here approximation for free convection;
        in Campbell & Norman sect 7:

        g_free['heat'] = 0.05 * [(Tsur - Tair) / d]**0.25,
        where d [m] is characteristic dimension.
        If assume moss stems as cylinders, d ~stem height.
        g_free['h2o'] = 1.09 x g_free['heat'],
        g_free['co2'] = 0.75 x g_free['co2']

        The pain in the ass is the characteristic dimension which is hard
        to define... for range of d and Tsur - Tair we get values that
        may indicate uncertainty range... in any case these are comparable
        to conductances due to forced convection (see fig)
    """

    schmidt_number_h2o = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_H2O
    reynolds_number = wind_speed * roughness_height / AIR_VISCOSITY

    # Rice et al. (2001) eq. 1
    conductance_h2o = (
        AIR_DENSITY
        * np.power(10, -3.18)
        * np.power(reynolds_number, 1.61)
        * MOLECULAR_DIFFUSIVITY_H2O
        / roughness_height
        * np.power(schmidt_number_h2o, 1.0/3.0))

    # [mol m-2 s-1], differ from H2O by ratio of molecular diffusivities
    conductance_co2 = (MOLECULAR_DIFFUSIVITY_CO2
                       / MOLECULAR_DIFFUSIVITY_H2O
                       * conductance_h2o)

    conductance_heat = (THERMAL_DIFFUSIVITY_AIR
                        / MOLECULAR_DIFFUSIVITY_H2O
                        * conductance_h2o)

    return {
        'co2': conductance_co2,
        'h2o': conductance_h2o,
        'heat': conductance_heat
        }


def moss_atm_conductance(wind_speed, roughness_height, dT=0.0, atten_factor=0.25):
    r""" Estimates boundary layer conductance of bryophyte canopy for paralell
    forced and free convection.

    Wind speed should represent vertical wind speed at ca. 20 cm above moss
    canopy (parametrsization derived from wind-tunnel studies). Roughness
    lengths scale is equal to the 'characteristic vertical height of
    individual bryophyte shoots' (typically order of 3-10 mm).

    Estimated CO\ :sub:`2` and heat conductances are related to water vapor
    by the ratio of molecular diffusivities
    (Campbell and Norman, 1998, eq. 7.29-7.33).

    References:
        Rice et al., 2001.
            Significance of variation in bryophyte canopy structure.
            Amer. J. Bot. 88:1568-1576.
        Rice, 2006.
            Towards an integrated undestanding of Bryophyte performance:
            the dimensions of space and time.
            Lindbergia 31:41-53.
        Schuepp, 1980.
            Observations on the use of analytical and numerical models for the
            description of transfer to porous surface vegetation such as
            lichen.
            Boundary-Layer Meteorol. 29: 59-73.
        Kondo & Ishida, 1997

    Args:
        wind_speed: [m s\ :sup:`-1`\ ]
        roughness_height: [m]
        deltaT: [degC], moss-air temperature difference
        atten_factor: [-] dimensionless attenuation factor for continuous moss carpets

    Returns:
        dictionary:
            boundary layer conductances for
                * 'co2': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'h2o': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
                * 'heat': [mol m\ :sup:`-2` s\ :sup:`-1`\ ]

    """

    Sc_v = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_H2O
    Sc_c = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_CO2
    Pr = AIR_VISCOSITY / THERMAL_DIFFUSIVITY_AIR
    Re = wind_speed * roughness_height / AIR_VISCOSITY


    # Rice et al. (2001) eq. 1: ShSc**0.33 = CRe**n, where C=6.6e-4 and n=1.61.
    # however, exponent n for individual species is <1.53 so use median values of model
    # fitted to individual species.
    C = 0.0067
    n = 1.27

    Sh_v = atten_factor * C*Re**n * Sc_v**0.33 # Sherwood numbner for H2O

    conductance_h2o = Sh_v * MOLECULAR_DIFFUSIVITY_H2O / roughness_height # ms-1

    # free convection as parallell pathway, based on Condo and Ishida, 1997.
    b = 2.2e-3 #ms-1K-1 b=1.1e-3 for smooth, 3.3e-3 for rough surface
    dT = np.maximum(dT, 0.0)
    gfree = Sc_v / Pr * b * dT**0.33  # mol m-2 s-1

    # [mol m-2 s-1]
    conductance_h2o = (conductance_h2o + gfree) * AIR_DENSITY


    # [mol m-2 s-1], differ from H2O by ratio of molecular diffusivities
    conductance_co2 = Sc_c / Sc_v * conductance_h2o

    conductance_heat = Pr / Sc_v * conductance_h2o

    return {
        'co2': conductance_co2,
        'h2o': conductance_h2o,
        'heat': conductance_heat
        }

def surface_atm_conductance(wind_speed, height, friction_velocity=None, dT=0.0, zom=0.01, b=1.1e-3):
    """
    Soil surface - atmosphere transfer conductance for scalars. Two paralell
    mechanisms: forced and free convection
    Args:
        wind_speed - wind speed (m/s) at zref
        height - reference height (m). Log-profile assumed below zref.
        zom - roughness height for momentum (m), ~0.1 x canopy height
        ustar - friction velocity (m/s) at log-regime. if ustar not given,
                it is computed from Uo, zref and zom
        b - parameter for free convection. b=1.1e-3 ... 3.3e-3 from smooth...rough surface
    Returns:
        conductances for CO2, H2O and heat (mol m-2 s-1), dict
    References:
        Schuepp and White, 1975:Transfer Processes in Vegetation by Electrochemical Analog,
        Boundary-Layer Meteorol. 8, 335-358.
        Schuepp (1977): Turbulent transfer at the ground: on verification of
        a simple predictive model. Boundary-Layer Meteorol., 171-186
        Kondo & Ishida, 1997: Sensible Heat Flux from the Earth’s Surface under
        Natural Convective Conditions. J. Atm. Sci.
    """

    Sc_v = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_H2O
    Sc_c = AIR_VISCOSITY / MOLECULAR_DIFFUSIVITY_CO2
    Pr = AIR_VISCOSITY / THERMAL_DIFFUSIVITY_AIR
    kv = 0.4  # von Karman constant (-)
    d = 0.0 # displacement height

    if friction_velocity == None:
        friction_velocity = wind_speed * kv / np.log((height - d) / zom)

    delta = MOLECULAR_DIFFUSIVITY_H2O / (kv*friction_velocity + EPS)

    gb_h = (kv * friction_velocity) / (Pr - np.log(delta / height))
    gb_v = (kv * friction_velocity) / (Sc_v - np.log(delta / height))
    gb_c = (kv * friction_velocity) / (Sc_c - np.log(delta / height))

    # free convection as parallel pathway, based on Condo and Ishida, 1997.
    #b = 1.1e-3 #ms-1K-1 b=1.1e-3 for smooth, 3.3e-3 for rough surface
    dT = np.maximum(dT, 0.0)

    gf_h = b * dT**0.33  # ms-1

    # mol m-2 s-1
    gb_h = (gb_h + gf_h) * AIR_DENSITY
    gb_v = (gb_v + Sc_v / Pr * gf_h) * AIR_DENSITY
    gb_c = (gb_c + Sc_c / Pr * gf_h) * AIR_DENSITY

#    plt.figure()
#    plt.plot(friction_velocity, gb_v, '-')
    return {'co2': gb_c, 'h2o': gb_v, 'heat': gb_h}

def evaporation_through_moss(properties,
                             volumetric_water,
                             moss_temperature,
                             forcing,
                             parameters):
    r""" Estimates soil evaporation rate through bryophyte layer.

    Evaporation in bryophyte layer is limited either by atmospheric demand
    and transport or soil supply of water.

    Water vapor flow from soil to air must overcome two resistances
    in series: 1. molecular diffusion through porous living moss, 2. molecular
    diffusion from moss canopy to 1st caluclation node in the atomosphere. The
    2nd resistance is assumed to be equal to conductance from wet moss canopy.

    Method does not use state variables of BryoType instance due to iteration.

    Args:
        properties (dict): characteristics of BryoType object
        volumetric_water: [m\ :sup:`3` m\ :sup:`-3`\ ]
        air_temperature: [\ :math:`^{\circ}`\ C]
        atm_partial_pressure_h2o: [Pa]
        wind_speed: [m s\ :sup:`-1`\ ]
        ambient_pressure: [Pa]
        soil_temperature: [\ :math:`^{\circ}`\ C]
            from 1st calculation node
        soil_hydraulic_head: [m]
            from 1st calculation node
        soil_hydraulic_conductivity: [m s\ :sup:`-1`\ ]
            from 1st calculation node
        soil_depth: [m]
            as negative value

    Returns:
        float:
            evaporation in [mol m\ :sup:`-2` s\ :sup:`-1`\ ]
    """
    # [mol mol-1] -> [Pa]
    h2o = forcing['h2o'] * forcing['air_pressure']

    # [mol/(m2 s)]
    # "from soil surface to moss canopy height"
    # change air_temperature to bryo_temperature
    moss_conductance = vapor_conductance_porous_media(
        volumetric_water,
        properties['porosity'],
        properties['height'],
        moss_temperature,
        forcing['air_pressure'])

    # [mol/(m2 s)]
    # "from moss canopy to the atmosphere"
    temp_difference = moss_temperature - forcing['air_temperature']
#    atm_conductance = moss_atm_conductance(forcing['wind_speed'],
#                                           properties['roughness_height'],
#                                           dT=temp_difference)
    atm_conductance = surface_atm_conductance(wind_speed=forcing['wind_speed'],
                                              height=parameters['reference_height'],
                                              friction_velocity=forcing['friction_velocity'],
                                              dT=temp_difference)

    # [mol/(m2 s)], two resistors in series
    conductance_h2o = (
        moss_conductance['h2o'] * atm_conductance['h2o']
        / (moss_conductance['h2o'] + atm_conductance['h2o']))

    # Assuming soil is saturated (rh = 1),calculate the maximum evaporation rate

    # [mol/(m2 s)]
    # atmospheric evaporative demand
    evaporative_demand = (conductance_h2o
                          * (saturation_vapor_pressure(forcing['soil_temperature']) - h2o)
                          / forcing['air_pressure'])

    # [-, fraction]
    relative_humidity = min(1.0, h2o
                         / saturation_vapor_pressure(forcing['air_temperature']))

    # [m], in equilibrium with atmospheric relative humidity
    atm_hydraulic_head = (
        GAS_CONSTANT
        * (DEG_TO_KELVIN + forcing['air_temperature'])
        * np.log(relative_humidity)
        / (MOLAR_MASS_H2O * GRAVITY))

    # [mol/(m2 s)]: 1e3 is water density kgm-3: 1e-3 kgm-3 / (kg/mol) x m/s = mol m-2 s-1
    evaporative_supply = max(0.0,
        1e3 / MOLAR_MASS_H2O * parameters['soil_hydraulic_conductivity']
        * ((atm_hydraulic_head - forcing['soil_water_potential']) / parameters['soil_depth'] - 1.0))

    # [mol/(m2 s)]
    # return min(evaporative_demand, evaporative_supply)
    return {'evaporative_demand': evaporative_demand,
            'evaporative_supply': evaporative_supply,
            'soil_evaporation': min(evaporative_demand, evaporative_supply)}


def effective_saturation(value, water_retention, input_unit):
    r""" Calculates effective saturation from volumetric water
    content (theta, [m3 m-3]) or from water potential (psi, [m])

    Args:
        value (float):
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is VOLUMETRIC_WATER
            * [m] if input is WATER_POTENTIAL
        water_retention (dict):
            if input_unit == 'water_potential'
                'alpha'
                'n'
            if input_unit == VOLUMETRIC_WATER
                'theta_r'
                'theta_s'
        input_unit (str): VOLUMETRIC_WATER or WATER_POTENTIAL
    Returns:
        float: [-]
    """

    options = set(['water_potential', 'volumetric_water'])

    if isinstance(input_unit, str):
        input_unit = input_unit.lower()
    else:
        raise ValueError("Input unit has to be string")

    if input_unit not in options:
        raise ValueError("Input unit options are: ".format(*options))


    if input_unit == 'water_potential':
        n = water_retention['n']
        m = 1.0 - np.divide(1.0, n)
        psi = 100.0 * np.minimum(value, 0.0)

        eff_sat = (1.0 + (water_retention['alpha'] * abs(psi)) ** n) ** -m

    elif input_unit == 'volumetric_water':

        theta_r = water_retention['theta_r']
        theta_s = water_retention['theta_s']

        theta = np.minimum(value, theta_s)
        theta = np.maximum(theta, theta_r)

        eff_sat = np.divide(((theta - theta_r) + EPS), (theta_s - theta_r) + EPS)

    return eff_sat


def convert_hydraulic_parameters(value, water_retention, input_unit):
    r""" Converts between volumetric water content and water potential.

    Note:
        In case of heterogenous poresystem, linear interpolation is used
        to convert effective saturation to water potential. Water retention
        curve with 500 points with logarithmic spacing in which water potential
        is ranging from -1e-6 to -1e2. In case of effective saturation is out
        of interpolation bounds, water potential is set to 0.0 m in lower bound
        or -1e2 m in upper bound.

    For conversion, VanGenuchten-model for water retention is used with
    following parameters:
            - saturated water content (theta_s, [cm\ :sup:`3` cm\ :sup:`-3`\ ])
            - residual water content (theta_r, [cm\ :sup:`3` cm\ :sup:`-3`\ ])
            - air entry suction (alpha, [cm\ :sup:`-1`\ ])
            - pore size distribution (n, [-])

    Volumetric water content is restricted to be between residual water
    content and saturated water content.

    Water potential is restricted to be between -10000 m and 0.0 m.

    Args:
        value:
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input is VOLUMETRIC_WATER
            * [m] if input is WATER_POTENTIAL
        water_retention (dict): water retension parameters of BryoType object
        input_unit:
                * VOLUMETRIC_WATER
                * WATER_POTENTIAL
    Returns:
        float:
            * [m\ :sup:`3` m\ :sup:`-3`\ ] if input was WATER_POTENTIAL
            * [m] if input was VOLUMETRIC_WATER
    """

    if isinstance(input_unit, str):
        input_unit.lower()

    theta_s = water_retention['theta_s']  # sat. water content [m3m-3]
    theta_r = water_retention['theta_r']  # residual water content [m3m-3]

    if input_unit == 'water_potential':
        eff_sat = effective_saturation(value,
                                       water_retention,
                                       input_unit)

        return eff_sat * (theta_s - theta_r) + theta_r

    elif input_unit == 'volumetric_water':

        alpha = water_retention['alpha']  # micropore air entry suction [cm-1]
        n = water_retention['n']  # micropore shape parameter [-]

        eff_sat = effective_saturation(value, water_retention, input_unit)
        inv_effs = 1.0 / eff_sat
        m = np.divide(n, n - 1.0)

        # [m]
        psi = 1e-2 * 1.0 / -alpha * (inv_effs ** m - 1.0) ** (1.0 / n)

        if isinstance(psi, list):
            psi = np.array(list)

        if isinstance(psi, np.ndarray):
            psi[psi <= -1000.0] = -1000.0
            #psi = np.where(psi <= -100.0, -100.0, psi)
        else:
            if psi <= -1000.0:
                return -1000.0

        return psi


def hydraulic_conductivity(water_potential, water_retention, method=None):
    r""" Estimates hydraulic conductivity in porous media.

    Hydraulic conductivity is based on Voortman et al. (2014) and
    Voortman et al. (2015) for xerophilious mosses.

    The implementation is valid for unimodal and multimodal poresystems.
    Multimodal Mualem-VanGenuchten is used as introduced by
    Priesack and Durner (2006).

    09.08.2018 AJK:
    At the moment hydraulic conduction used only for micropore
    system (capillary pore system) and it is only needed to calculate capillary
    rise. Macropore system does not retain water and it is assumed that
    macropore system does not restrict hydraulic conduction.

    References:
        Voortman et al. (2014) Hydrological Processes, 28, 6251-6264
        Priesack and Durner (2006) Vadose Zone Journal, 5, 121-124

    Lauren:
    Estimation is based on Lauren (1999, table 4) and parametrization is valid
    from -4 to -80 kPa. Assumption is that in the dry states
    (lower than -80 kPa) there are no hydraulic conductivity.

    Args:
        water_potential: [m]
        water_retention (list):
            parameters of water retention curve
                0. saturated water content (theta_s) [m\ :sup:`3` m :sup:`-3`\ ]
                1. residual water content (theta_r) [m\ :sup:`3` m :sup:`-3`\ ]
                2. air entry suction (alpha) [cm\ :sup:`-1`]
                3. pore size distribution (n) [-]
                4. saturated conductivity (K_s) [m s\ :sup:`-1`]
                5. pore connectivity (l) [-]

    Returns:
        float: hydraulic conductivity in [m s\ :sup:`-1`\ ]
    """

    if isinstance(method, str):
        method.lower()

    saturated_conductivity = water_retention['saturated_conductivity']

    if method == 'lauren':
        # kPa (standard gravity 10 m/s2)
        saturated_conductivity = water_retention['saturated_conductivity']
        water_potential = -10.0 * water_potential

        # Assuming that in dry states there is no hydraulic conductivity
        if water_potential > 80.0:
            return 0.0

        if water_potential < 4.0:
            water_potential = 4.0

        # relative conductivity respect to 4 kPa
        conductivity = (
            np.power(10.0, 1.0/(-0.62 + 0.26 * np.log10(water_potential)))
            / np.power(10.0, 1.0/(-0.62 + 0.26 * np.log10(4)))
            )
        # [m/s]
        return conductivity * saturated_conductivity

    else:
        # Possibility to add more poresystems

        psi = 100.0 * np.minimum(water_potential, 0.0)

        micropores = {'alpha': water_retention['alpha'],
                      'n': water_retention['n']}

        poresystems = [micropores]

        coefficients = []
        denominators = []
        nominators = []

        for idx in range(len(poresystems)):
            m = 1.0 - np.divide(1.0, poresystems[idx]['n'])
            alpha = poresystems[idx]['alpha']

            sat_eff = effective_saturation(psi,
                                           poresystems[idx],
                                           'water_potential')

            coefficients.append(sat_eff)

            denom = 1.0 - np.power(sat_eff, 1.0 / m)
            denom = 1.0 - np.power(denom, m)
            denominators.append(alpha * denom)

            nominators.append(alpha)

        pore_connectivity = water_retention['pore_connectivity']

        coefficient = np.power(np.sum(coefficients, axis=0), pore_connectivity)

        denominator = np.power(np.sum(denominators, axis=0), 2.0)

        nominator = np.power(np.sum(nominators, axis=0), 2.0)

        return saturated_conductivity * coefficient * (denominator / nominator)


def emitted_longwave_radiation(temperature, properties=None):
    r""" Estimates emitted longwave radiation

    Args:
        temperature (float): [W m\ :sup:`-2`]
        properties (dict/float): properties dictionary or emissivity
    Returns:
        (float): [W m\ :sup:`-2`]
    """
    if isinstance(properties, dict):
        emissivity = properties['optical_properties']['emissivity']

    elif isinstance(properties, float):
        emissivity = properties

    else:
        emissivity = 0.98

    emitted_longwave_radiation = (
        emissivity * STEFAN_BOLTZMANN
        * np.power((temperature + DEG_TO_KELVIN), 4.0))

    return emitted_longwave_radiation


def bryophyte_shortwave_albedo(water_content, properties=None):
    r""" Bryophyte albedo for PAR and NIR regions as function of water content

    The effect of water content of spectral properties are based on studies by
    Vogelmann and Moss (1993) and Fernandes (1999)on Sphagnum cuspidatum and
    Pleurozium schreberi, respectively.

    The albedo is scaled specific reflectance for PAR (400-700 nm) or
    NIR (750-1400) regions. The scaling coefficient is common for both
    PAR and NIR and it is based on relationship between normalized
    reflectaces and hydration status. The species specific albedo is
    assumed to represent a reflectance when bryophyte is in full hydration.

    If bryophyte's properties are not given, estimation is based on generic
    fit of water content against reflectances separately for PAR and NIR.
    Fits are based on studies by Vogelmann and Moss (1993), and
    Fernandes (1999) on Sphagnum cuspidatum and Pleurozium schreberi,
    respectively.

    References:
        Vogelmann and Moss (1993)
            Remote Sensing of Environment 45:273-279.
        Fernandes (1999)
            PhD thesis entitled: 'Scale influences of surface
            parametrization on modelled boreal carbon and water budgets'

    Args:
        water_content (float): [g g\ :sup:`-2`]
        max_water_content (float): [g g\ :sup:`-2`]

    Returns:
            list: reflectances
                0. par albedo
                1. nir albedo
    """

    if properties is None:

        def reflectance(x, a, b):
            r""" Describes relationship between water content and reflectance

            Args:
                x (float): water_content
                a (float): fitting parameter
                b (float): fitting parameter

            Returns:
                percentage (float)
            """
            return np.abs(a * np.power(x, b))

        # Original reflectance in percents
        if isinstance(water_content, float):
            if water_content < 0.3:
                albedo_nir = 68.41/100.0
            else:
                albedo_nir = reflectance(water_content, 47.66, -0.3002) / 100.0
            albedo_par = reflectance(water_content, 8.84, -0.1303) / 100.0
        else:
            albedo_nir = np.empty(np.shape(water_content))
            albedo_nir = reflectance(water_content, 47.66, -0.3002) / 100.0
            albedo_nir[water_content < 0.3] = 68.41 / 100.0

            albedo_par = reflectance(water_content, 8.84, -0.1303) / 100.0

        return {'PAR': albedo_par, 'NIR': albedo_nir}

    else:
        albedo_nir = properties['optical_properties']['albedo_NIR']
        albedo_par = properties['optical_properties']['albedo_PAR']
        normalized_water_content = water_content / properties['max_water_content']

        scaling_coefficient = 1.0 + (4.5 - 1.0) / (1.00 + np.power(10, 4.53 * normalized_water_content))

        albedo_par = scaling_coefficient * albedo_par
        albedo_nir = scaling_coefficient * albedo_nir

        return {'PAR': albedo_par, 'NIR': albedo_nir}

#def capillary_rise(dt,
#                   properties,
#                   hydraulic_conductivity,
#                   water_potential,
#                   water_content,
#                   soil_hydraulic_conductivity,
#                   soil_water_potential,
#                   soil_depth):
#    r""" Estimates liquid water flux from soil to moss
#    :math:`q = -k(\\frac{\partial h}{\partial z}+1)` [mm s\ :sup:`-1`] to moss
#    as capillary rise from underlying soil. Capillary rise is further limited
#    by bryophytes available water storage capacity.
#
#    For water flow in soil, multiply capillary_rise by bryophyte
#    ground_coverage and add that to the root_sink.
#
#    Args:
#        dt: [s]
#        properties: dictionary containing characteristics of BryoType object
#            * 'height'
#            * 'max_water_content'
#            * 'dry_mass'
#        hydraulic_conductivity: [m s\ :sup:`-1`\ ]
#        water_potential: [m]
#        water_content: [g g\ :sup:`-1`\ ]
#        soil_hydraulic_conductivity: [m s\ :sup:`-1`\ ]
#            from 1st calculation node
#        soil_water_potential: [m]
#            from 1st calculation node
#        soil_depth: [m]
#
#    Returns:
#        float:
#            capillary rise in [mm s\ :sup:`-1`\ ]
#    """
#
#    # [m/s] geometric average of hydrological conductivities of
#    # soil and bryophyte
#    conductivities = hydraulic_conductivity * soil_hydraulic_conductivity
#    k = np.power(conductivities, 0.5)
#
#    # [mm/s] or [kg/(m2 s)]
#    capillary_rise = (1.0e3
#                      * max(0.0,
#                            - k * ((water_potential - soil_water_potential)
#                                   / (properties['height'] - soil_depth) + 1.0)))
#
#    capillary_rise = min(
#        capillary_rise,
#        (properties['max_water_content'] - water_content)
#        * properties['dry_mass'] / dt)
#
#    return capillary_rise
