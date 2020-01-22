#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 12:39:09 2018

@author: ajkieloaho
"""

from pyAPES_utilities.soiltypes.organic import soil_properties, zh

ranges = {}

def get_parameters(scenario):
    # spefify as one values (same for all simulations) or tuple of length 'count'
    if scenario.upper() == 'CASE_1':
        parameters = {
                'count': 1,
                'general':{
                        'start_time' : "2005-06-01",
                        'end_time' : "2005-07-01"
                        },
                'canopy': {
                        # 'ctr': { # controls
                        #     'WMA': True,  # well-mixed assumption
                        #     'Ebal': False,  # no energy balance
                        #         },
                        },
            'soil': {
                    'grid': {
                            'zh': zh
                            },
                    'soil_properties': soil_properties,
                    'water_model': {
                            # 'type': 'Equilibrium',
                            'initial_condition':{
                                    'ground_water_level': -0.2
                                    },
                            'lower_boundary': {  # lower boundary condition (type, value, depth)
                                   'type': 'impermeable',
                                   'value': None,
                                   'depth': -2.0
                                   },
                           'drainage_equation': {  # drainage equation and drainage parameters
                                   'type': 'Hooghoudt',  #
                                   'depth': 0.8,  # drain depth [m]
                                   'spacing': 45.0,  # drain spacing [m]
                                   'width': 1.0,  # drain width [m]
                                   }
                            }
                    }
            }
        return parameters
    else:
        raise ValueError("Unknown parameterset!")

def iterate_parameters(parameters, default, count):
    """ Going through recursively senstivity nested parameter dictionary.
    Args:
        paramters (dict): nested dictionary
    """

    for key, value in parameters.items():
        if key == 'count':
            continue
        elif isinstance(value, dict):
            if key in default:
                default[key] = iterate_parameters(value, default[key], count)
            else:
                print(key + ' not in parameters')

        else:
            if isinstance(value, tuple):
                default[key] = value[count]
            else:
                default[key] = value

    return default
