#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:23:28 2018

Module contains functions to calculate carbon exchange in bryophytes.

Note:
    migrated to python3
    - no changes made

@author: ajkieloaho
"""

import numpy as np
from canopy.constants import EPS

class BryophyteCarbon(object):
    def __init__(self, para, carbon_pool=0):
        if 'photosynthesis' in para:
            self.amax = para['photosynthesis']['amax']
            self.b = para['photosynthesis']['b']
            self.moisture_coeff = para['photosynthesis']['moisture_coeff']
            self.temperature_coeff = para['photosynthesis']['temperature_coeff']
        else:
            self.amax = 0.0
            self.b = EPS
        
        self.r10 = para['respiration']['r10']
        self.q10 = para['respiration']['q10']
        
        self.max_water_content = para['max_water_content'] # g g-1
        self.carbon_pool = carbon_pool
        
    def carbon_exchange(self,
                        water_content,
                        temperature,
                        incident_par):
        r""" Estimates photosynthesis and respiration rates of bryophyte layer.
    
        Photosynthesis is restricted by both tissue water content
        (dry conditions) and excess water film on leaves (diffusion limitation)
        as in Williams and Flanagan (1996). Water content
        coefficients are 3rd order polynomial fitted to the data represented by
        Williams and Flanagan (1996) and used to calculate effect of water
        content on photosynthesis.
    
        Empirical modifier of photosynthesis due to water content assumes that
        both light-limited and Rubisco-limited assimilation of carbon are
        affected similarly. This seems to apply for Pleurozium and Sphagnum
        when normalized water content is used as scaling. Assumes that
        there is always 5 percents left in photosynthetic capacity.Empirical
        modifier of photosynthesis due to temperature is based on
        late growing season presented in Fig. 2 in Frolking et al. (1996).
    
        References:
            Frolking et al. (1996)
                Global Change Biology 2:343-366
            Williams and Flanagan (1996)
                Oecologia 108:38-46
    
        Args:
            water_content (float): [g g-1]
            temperature (float): [degC]
            incident_par (float): [W m\ :sup:`-2`]
            

    
        Returns:
            dictionary:
                * 'photosynthesis_rate':
                  [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
                * 'respiration_rate':
                  [\ :math:`\mu`\ mol m\ :sup:`-2`:sub:`ground` s\ :sup:`-1`\ ]
        """
        # check inputs
        # [umol/(m2 s)]
        incident_par = np.maximum(EPS, 4.56 * incident_par)
        
        normalized_water_content = water_content / self.max_water_content
        

        # hyperbolic light response at community level [umolm-2(ground) s-1]
        light_response = self.amax * incident_par / (self.b + incident_par )
    
        # moisture and temperature responses [-]
        water_modifier = (self.moisture_coeff[3]
                          + self.moisture_coeff[2] * normalized_water_content
                          + self.moisture_coeff[1] * normalized_water_content ** 2.0
                          + self.moisture_coeff[0] * normalized_water_content ** 3.0)
    
        water_modifier = np.maximum(0.05, water_modifier)
    
        temperature_modifier = (self.temperature_coeff[3]
                                + self.temperature_coeff[2] * temperature
                                + self.temperature_coeff[1] * temperature ** 2.0
                                + self.temperature_coeff[0] * temperature ** 3.0)
    
        temperature_modifier = np.maximum(0.01, temperature_modifier)
    
        temperature_modifier = np.minimum(1.0, temperature_modifier)
    
        # [umol m-2 (leaf) s-1]
        photosynthetic_rate = light_response * water_modifier * temperature_modifier
        
        # --- respiration rate [umol m-2 (ground) s-1]
        
        """ replace with smooth function and move as parameter """
        if water_content < 7.0:
            water_modifier_respiration = (
                -0.45 + 0.4 * water_content
                - 0.0273 * water_content ** 2)
        else:
            water_modifier_respiration = (
                -0.04 * water_content + 1.38)
    
        # effect of water content is in the range 0.01 to 1.0
        water_modifier_respiration = np.maximum(0.01, np.minimum(1.0, water_modifier_respiration))
    
        # r = r10 * Q10^((T-10) / 10) [umol/(m2 s)]
        respiration_rate = (
            self.r10 * self.q10**((temperature - 10.0) / 10.0) * water_modifier_respiration
            )
    

        return {
            'photosynthesis_rate': photosynthetic_rate,
            'respiration_rate': respiration_rate,
            'net_co2_flux': -photosynthetic_rate + respiration_rate
            }



def soil_respiration(properties, Ts, Wliq, Wair):
    """ Soil respiration beneath forestfloor

    Heterotrophic and autotrophic respiration rate (CO2-flux) based on
    Pumpanen et al. (2003) Soil.Sci.Soc.Am

    Restricts respiration by soil moisuture as in
    Skopp et al. (1990), Soil.Sci.Soc.Am

    Args:
        properties (dict):
            'R10'
            'Q10'
            'poros'
            'limitpara'
        Ts - soil temperature [degC]
        Wliq - soil vol. moisture content [m3 m-3]
    Returns:
        rsoil - soil respiration rate [umol m-2 s-1]
        fm - relative modifier (Skopp et al.)
    """
    # Skopp limitparam [a,b,d,g] for two soil types
    # sp = {'Yolo':[3.83, 4.43, 1.25, 0.854],
    #       'Valentine': [1.65,6.15,0.385,1.03]}

    limitpara = properties['limitpara']
    r10 = properties['R10']
    q10 = properties['Q10']

    # unrestricted respiration rate
    base_respiration = r10 * np.power(q10, (Ts - 10.0) / 10.0)

    # moisture response (substrate diffusion, oxygen limitation)
    modifier = np.minimum(limitpara[0] * Wliq**limitpara[2],
                          limitpara[1] * Wair**limitpara[3])  # ]0...1]
    modifier = np.minimum(modifier, 1.0)

    respiration = base_respiration * modifier

    return respiration
