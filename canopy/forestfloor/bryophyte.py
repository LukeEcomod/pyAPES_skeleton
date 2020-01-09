#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module: bryotype
    :synopsis: APES-model component
.. moduleauthor:: Antti-Jussi Kieloaho

Bryotype describes structural and functional properties and processes of
moss/lichen (bryophytes) species/groups at the forest bottom layer. Original
implementation is done in MatLab by Samuli Launiainen.

Created on Tue Mar 14 08:15:50 2017

TODO:
    only water_balance calculation.

References:


Note:
    migrated to python3
    - absolute imports
    - one dictionary comprhension changed

CHANGES (13.-14.7.2017 SL):
    * BryoType: added state variables 'hydraulic_conductivity',
        'thermal_conductivity', 'carbon_pool'
    * changes names of some functions ('estimate_' --> 'moss_'),
        some simplifications of scrips
    * energy_and_water_balance: returns separate dicts 'fluxes' & 'states'
    * arg 'time' --> 'dt' in all calls
    * _run_timestep(args): simple entry-point to energy_and_water_balance
        and carbon_exchange
"""

import numpy as np

from canopy.constants import WATER_DENSITY, MOLAR_MASS_H2O, MOLAR_MASS_C, LATENT_HEAT, EPS
from .heat_and_water import heat_and_water_exchange, water_exchange
from .heat_and_water import convert_hydraulic_parameters, evaporation_through_moss
from .carbon import carbon_exchange


class Bryophyte(object):
    r""" Represents bryophyte community-soil-atmosphere interactions.

    Characteristics of BryoType object are stored in 'properties' dictionary.
    These describes physical characteristics of a bryophyte layer.
    """

    # pylint: disable=too-many-instance-attributes
    # instance attributes are necessary to functioning of BryoModel

    def __init__(self, properties, initial_conditions=None):
        r""" Initialises a bryophyte object by using bryophyte's properties and
        initial states.

        Volumetric water content, relative water content is assumed to
        be equal to maximal retention capacity at field capacity.

        Leaf area index (*LAI*) is calculated as follows

        .. math::
            LAI = 1\\mathrm{e}^{3} \\frac{m_{dry} SLA}{1\\mathrm{e}^{4}}

        where :math:`m_{dry}` is dry mass and *SLA* is specific leaf area.

        Args:
            properties (dict):
                'ground_coverage':  [-]
                'height': [m]
                'roughness_height': [m]
                'leaf_area_index': [m\ :sup:`2` m :sup:`-2`\ ]
                'specific_leaf_area': [m\ :sup:`3` m :sup:`-3`\ ]
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
                    'compressability':
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

            initial_conditions (dict):
                initial_temperature: [\ :math:`^{\circ}`\ C]
                initial_water_content: [g g\ :sup:`-1`\ ]
        """

        dry_mass = properties['bulk_density'] * properties['height']
        properties['dry_mass'] = dry_mass

        residual_water_content = (properties['min_water_content']
                                  / WATER_DENSITY
                                  * properties['bulk_density'])

        field_capacity = (properties['max_water_content']
                          / WATER_DENSITY
                          * properties['bulk_density'])

        saturated_water_content = field_capacity

        if 'water_retention' not in properties:
            water_retention = {}

            water_retention['theta_r'] = residual_water_content
            water_retention['theta_s'] = saturated_water_content
            water_retention['field_capacity'] = field_capacity

            if 'saturated_conductivity' in properties:
                water_retention['saturated_conductivity'] = properties['saturated_conductivity']
            else:
                raise ValueError('Water retention parameter saturated_conductivity is missing')

            if 'alpha' in properties and 'n' in properties:
                water_retention['alpha'] = properties['alpha']
                water_retention['n'] = properties['n']
            else:
                raise ValueError('Water retention parameteres alpha and n are missing')

            if 'pore_connectivity' in properties:
                water_retention['pore_connectivity'] = properties['pore_connectivity']

            if 'compressibility' in properties:
                water_retention['compressability'] = properties['compressability']

            properties['water_retention'] = water_retention

        else:
            water_retention = properties['water_retention']
            water_retention['field_capacity'] = field_capacity
            # theta_s and theta_r are fixed to
            # min_water_content and max_water_content
            water_retention['theta_r'] = residual_water_content
            water_retention['theta_s'] = saturated_water_content


        self.properties = properties

        # set initial conditions
        if initial_conditions is not None:
            #: [:math:`^{\circ}`\ C]
            self.temperature = initial_conditions['temperature']

            #: [g g\ :sup:`-1`\ ]
            if initial_conditions['water_content'] <= properties['max_water_content']:
                self.water_content = initial_conditions['water_content']

            else:
                self.water_content = properties['max_water_content']

        else:
            #: [:math:`^{\circ}`\ C]
            self.temperature = 10.

            #: [g g\ :sup:`-1`\ ]
            self.water_content = (
                    properties['max_water_content']
                    + properties['min_water_content']) / 2.0

        self.water_storage = self.water_content * properties['dry_mass']

        #: [m\ :sup:`3` m\ :sup:`-3`\ ]
        self.volumetric_water = (
            self.water_content / WATER_DENSITY * properties['bulk_density'])

        #: [m]
        self.water_potential = convert_hydraulic_parameters(
                self.volumetric_water,
                self.properties['water_retention'],
                'volumetric_water')

        #: [kg C m-2], 'free carbon pool'
        self.carbon_pool = 0.0
        self.coverage = properties['ground_coverage']

        self.old_carbon_pool = self.carbon_pool
        self.old_water_content = self.water_content
        self.old_water_storage = self.water_storage
        self.old_volumetric_water = self.volumetric_water
        self.old_water_potential = self.water_potential
        self.old_temperature = self.temperature

    def update(self):
        """ Updates old states to states after converged iteration.
        """

        self.old_carbon_pool = self.carbon_pool
        self.old_water_content = self.water_content
        self.old_water_storage = self.water_storage
        self.old_volumetric_water = self.volumetric_water
        self.old_water_potential = self.water_potential
        self.old_temperature = self.temperature

    def restore(self):
        """ Restores new states back to states before iteration.
        """

        self.carbon_pool = self.old_carbon_pool
        self.water_content = self.old_water_content
        self.water_storage = self.old_water_storage
        self.volumetric_water = self.old_volumetric_water
        self.water_potential = self.old_water_potential
        self.temperature = self.old_temperature

    def run(self, dt, forcing, parameters, controls):
        r""" Calculates one timestep and updates states of Bryophyte instance.

        Args:
            dt: timestep [s]
            forcing (dict):
                'throughfall': [kg m\ :sup:`-2`\ s\ :sup:`-1`\ ]
                'par': [W m\ :sup:`-2`\ ]
                'nir': [W m\ :sup:`-2`\ ] if energy_balance is True
                'lw_dn': [W m\ :sup:`-2`\ ] if energy_balance is True
                'h2o': [mol mol\ :sup:`-1`\ ]
                'air_temperature': [\ :math:`^{\circ}`\ C]
                'air_pressure': [Pa]
                'soil_temperature': [\ :math:`^{\circ}`\ C]
                'soil_water_potential': [m]
            parameters (dict):
                'reference_height' [m]
                'soil_depth': [m]
                'soil_hydraulic_conductivity': [m s\ :sup:`-1`\ ]
                'soil_thermal_conductivity': [W m\ :sup:`-1`\  K\ :sup:`-1`\ ]
                    if energy_balance is True
            controls (dict):
                'energy_balance': boolean
                'solver': 'forward_euler', 'odeint'
                'nsteps' number of steps in odesolver
                'logger_info': str

        Returns:
            fluxes (dict)
            states (dict)
        """

        if controls['energy_balance']:
            # calculate moss energy and water balance
            fluxes, states = heat_and_water_exchange(
                properties=self.properties,
                temperature=self.old_temperature,
                water_content=self.old_water_content,
                dt=dt,
                forcing=forcing,
                parameters=parameters,
                solver=controls['solver'],
                nsteps=controls['nsteps'],
                logger_info=controls['logger_info']
            )

        else:
            # only water balance
            fluxes, states = water_exchange(
                dt=dt,
                water_storage=self.old_water_storage,
                properties=self.properties,
                forcing=forcing,
                parameters=parameters
            )

        # update state variables
        self.temperature = states['temperature']
        self.water_content = states['water_content']
        self.water_storage = states['water_storage']
        self.volumetric_water = states['volumetric_water']
        self.water_potential = states['water_potential']

        # solve photosynthesis and respiration

        # [umol m-2(ground) s-1]
        cflx = carbon_exchange(
            self.properties,
            self.water_content,  # old
            self.temperature,  # old
            forcing['par'])

        nee = -cflx['photosynthesis_rate'] + cflx['respiration_rate']
        fluxes.update({
            'photosynthesis_rate': cflx['photosynthesis_rate'],
            'respiration_rate': cflx['respiration_rate'],
            'nee': nee})

        # update bryophyte free carbon pool (g C m-2) of bryophyte
        self.carbon_pool = self.old_carbon_pool + 1e3 * MOLAR_MASS_C * 1e-6 * nee
        states['carbon_pool'] = self.carbon_pool

        # compute soil evaporation through moss layer

        # [mol m-2 s-1]
        soil_evaporation = evaporation_through_moss(
            properties=self.properties,
            volumetric_water=self.volumetric_water,  
            moss_temperature=self.temperature,
            forcing=forcing,
            parameters=parameters)

        # unit conversion from mol m-2 s-1 to kg m-2 s-1
        soil_evaporation = {key: value * MOLAR_MASS_H2O for key, value in soil_evaporation.items()}

        fluxes.update(soil_evaporation)

        return fluxes, states
