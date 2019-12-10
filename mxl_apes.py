# -*- coding: utf-8 -*-
"""
.. module: pyAPES
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Model framework for soil-bottom layer-canopy-atmosphere interactions.
Based on MatLab implementation by Samuli Launiainen.

Created on Tue Oct 02 09:04:05 2018

Note:
    migrated to python3
    - print on same line
    - dict.keys(), but these are iterated after in for-each-loop

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy –
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.

To call model and read results:
    from tools.iotools import read_results
    from pyAPES import driver
    outputfile = driver(create_ncf=True, dbhfile="letto2014.txt")
    results = read_results(outputfile)
"""

from copy import deepcopy as copy
import numpy as np
import time
import logging

from tools.iotools import initialize_netcdf,  write_ncf
from mxl.utils import read_mxl_forcing
from mxl.mxl import MXLmodel
from mxl.canopy_dev import CanopyModel
from mxl.mxl import GAS_CONSTANT, MOLAR_MASS_H2O, MOLAR_MASS_AIR, DEG_TO_KELVIN, CP_AIR_MASS

#from soil.soil import Soil

#from parameters.parametersets import iterate_parameters


def driver(create_ncf=False, result_file=None, parametersets={}):
    """
    Args:
        create_ncf (bool): results saved to netCDF4 file
        result_file (str): name of result file
        parametersets (dict): parameter sets to overwrite default parameters
    """
    # --- LOGGING ---
    from parameters.mxl_apes_parameters import logging_configuration
    from logging.config import dictConfig
    dictConfig(logging_configuration)

    # --- PARAMETERS ---
    from parameters.mxl_apes_parameters  import gpara, cpara, mxlpara, mxl_ini

    if parametersets == {}:
        Nsim = 1
    else:
        Nsim = parametersets['count']

    # --- FORCING ---
    # Read forcing
    forcing = read_mxl_forcing(gpara['forc_filename'],
                               gpara['start_time'],
                               gpara['end_time'],
                               dt0=gpara['dt0'],
                               dt=gpara['dt'])


    #default_params = {
    #        'canopy': cpara,
    #        }

    #param_space = [iterate_parameters(parametersets, copy(default_params), count) for count in range(Nsim)]

    logger = logging.getLogger(__name__)

    logger.info('Simulation started. Number of simulations: {}'.format(Nsim))


    tasks = []

    for k in range(Nsim):
        #tasks.append(Model(gpara, param_space[k]['canopy'], param_space[k]['soil'], forcing, nsim=k))
        tasks.append(Model(forcing, gpara, cpara, mxlpara, mxl_ini, nsim=k))
        
    if create_ncf:
        timestr = time.strftime('%Y%m%d%H%M')
        if result_file:
            filename = result_file
        else:
            filename = timestr + '_results.nc'

        ncf, _ = initialize_netcdf(
                gpara['variables'],
                Nsim,
                tasks[k].Nsoil_nodes,
                tasks[k].Ncanopy_nodes,
                tasks[k].Nplant_types,
                forcing,
                filepath=gpara['results_directory'],
                filename=filename)

        for task in tasks:
            logger.info('Running simulation number (start time %s): %s' % (
                        time.strftime('%Y-%m-%d %H:%M'), task.Nsim))
            running_time = time.time()
            results = task.run()
            logger.info('Running time %.2f seconds' % (time.time() - running_time))
            write_ncf(nsim=task.Nsim, results=results, ncf=ncf)

            del results

        output = gpara['results_directory'] + filename
        logger.info('Ready! Results are in: ' + output)
        ncf.close()

    else:
        running_time = time.time()
        results = {task.Nsim: task.run() for task in tasks}
        output = results
        logger.info('Running time %.2f seconds' % (time.time() - running_time))
    
    return output, tasks[0]

class Model(object):
    """
    pyAPES - MXL full model
    """
    def __init__(self, forcing, gen_para, canopy_para, mxl_para, mxl_initial, nsim=0):

        self.dt = gen_para['dt']
        
        self.Nsteps = len(forcing)
        self.forcing = forcing
        self.Nsim = nsim

        self.Nsoil_nodes = canopy_para['ground']['soildepth']
        self.Ncanopy_nodes = canopy_para['grid']['Nlayers']

        # create canopy model instance
        self.canopy_model = CanopyModel(canopy_para)
        
        self.Nplant_types = len(self.canopy_model.planttypes)

        # create mxl model instance
        #inistate = mxlpara['initial_state']
        self.mxl_model = MXLmodel(mxl_initial, mxl_para)
        self.mxl_model.run_timestep(0.01, 0.0, 0.0, 0.01, out=False)
        print(self.mxl_model.__dict__)
        self.results = _initialize_results(gen_para['variables'],
                                           self.Nsteps,
                                           self.Nsoil_nodes,
                                           self.Ncanopy_nodes,
                                           self.Nplant_types,
                                           Nasl_nodes = np.fix(5000. / self.canopy_model.dz)
                                          )
        

    
    def run(self):
        """ Run atmosphere-canopy-soil--continuum model"""

        logger = logging.getLogger(__name__)
        logger.info('Running simulation {}'.format(self.Nsim))

        #print('RUNNING')
        k_steps=np.arange(0, self.Nsteps, int(self.Nsteps/10))
        for k in range(0, self.Nsteps):
            #print('step=', k)
            # progress bar
            if k in k_steps[:-1]:
                s = str(np.where(k_steps==k)[0][0]*10) + '%'
                print('{0}..'.format(s), end=' ')

            
            # values at mixed layer bottom: these change each timestep
            Psurf = 1e3*self.mxl_model.Psurf
            LW_in = self.forcing['LWin'].iloc[k]
            
            U_mxl = self.mxl_model.U
            T_mxl = self.mxl_model.theta - DEG_TO_KELVIN
            CO2_mxl = self.mxl_model.ca
            H2O_mxl = self.mxl_model.rh * self.mxl_model.esat / self.mxl_model.Psurf 

            # asl depth [m]
            z_asl = max(max(self.canopy_model.z), 0.1 * self.mxl_model.h)            
            #print('U', U_mxl, 'T', T_mxl, 'CO2', CO2_mxl, 'H2O', H2O_mxl)
                        
            #print('Uo', U_mxl, 'To', T_mxl, 'H2Oo', H2O_mxl, 'Cao', CO2_mxl, 'Ps', Psurf)
            
            """ Canopy model """
            canopy_forcing = {
                'zenith_angle': self.forcing['Zen'].iloc[k],
                'PAR': {'direct': self.forcing['dirPar'].iloc[k],
                        'diffuse': self.forcing['diffPar'].iloc[k]},
                'NIR': {'direct': self.forcing['dirNir'].iloc[k],
                        'diffuse': self.forcing['diffNir'].iloc[k]},
                'air_pressure': Psurf,
                'precipitation': 0.0, # 1e-3*self.forcing['Prec'].iloc[k], # precip in m/s = 1e-3 kg/m2
                
                # below should vary with MXL -state
                'lw_in': LW_in,
                'wind_speed': U_mxl,
                'air_temperature': T_mxl,
                'h2o': H2O_mxl,
                'co2': CO2_mxl,
                
                # soil 1st node state
                #'soil_temperature': Ts, #degC
                #'soil_volumetric_water': Ws, #m3/m3
                #'soil_water_potential': Psis # m
            }

            canopy_parameters = {
                'soil_depth': 0.1,
                'soil_hydraulic_conductivity': 1e-6,
                'soil_thermal_conductivity': 1.0,
                'date': self.forcing.index[k]
            }

            
            # run self.canopy_model
            canopy_flux, canopy_state, ffloor_flux, asl_profs = self.canopy_model.run(
                dt=self.dt,
                z_asl = z_asl,
                forcing=canopy_forcing,
                parameters=canopy_parameters
            )

            """ --- solve mxl --- """
            #unit conversions & mxl forcing
            H = canopy_flux['SH']
            LE = canopy_flux['LE']
            NEE = canopy_flux['NEE']
            ust = asl_profs['ust'][-1]
            Tasl = asl_profs['T'][-1]
            h2o = asl_profs['H2O'][-1]

            wt, wq, wc = _to_kinematic_fluxes(H, LE, NEE, Tasl, Psurf, h2o)
            
            #forcing_mxl = {'ust': ust, 'wt': wt, 'wq': wq, 'wc': wc}
                           
            # --- run mxl growth
            # for iterations: make copy of self.mxl_model and after convergence, replace back
            self.mxl_model.run_timestep(wt, wq, wc, ust, out=False)

            # --- output results
            
            canopy_state.update(canopy_flux)
            
            ffloor_state = {}
            ffloor_state.update(ffloor_flux)

            mxl_state = {'h': self.mxl_model.h,
                         'h_lcl': self.mxl_model.h_lcl,
                         'theta': self.mxl_model.theta,
                         'q': self.mxl_model.q,
                         'ca': self.mxl_model.ca,
                         'U': self.mxl_model.U,
                         'vpd': self.mxl_model.vpd,
                         'wt': wt,
                         'wq': wq,
                         'wc': wc,
                         'ust': ust
                         }
            
#            for nn in canopy_state.keys():
#                if isinstance(canopy_state[nn], np.ndarray):
#                    print(nn, len(canopy_state[nn]))

            self.results = _append_results('forcing', k, canopy_forcing, self.results)
            self.results = _append_results('canopy', k, canopy_state, self.results)
            self.results = _append_results('mxl', k, mxl_state, self.results)
            #self.results = _append_results('ffloor', k, ffloor_state, self.results)
            #self.results = _append_results('soil', k, soil_state, self.results)

        print('100%')


        return self.results

def _to_kinematic_fluxes(H, LE, Fc, T, P, h2o):
    """
    Converts surface fluxes to kinematic fluxes
    Args:
        H (Wm-2)
        LE (Wm-2)
        Fc (umolm-2s-1)
        T (degC)
        P (Pa)
        h2o (mol/mol)
    Returns:
        wt (K ms-1)
        wq (kg kg-1 ms-1)
        wc (ppm ms-1)
    """

    h2o = h2o * P # Pa
    Pdry = P - h2o  # Pa pressure of dry air

    rhoa = (Pdry*MOLAR_MASS_AIR + h2o*MOLAR_MASS_H2O) / (GAS_CONSTANT*(T + DEG_TO_KELVIN))  # kg m-3
    Lv = 1.0e6*(2.501 - 2.361e-3*T)  # J kg-1    
    Mair = P / (GAS_CONSTANT * (T + DEG_TO_KELVIN)) # m3 mol-1
    #print(np.mean(Mair), np.mean(rhoa))
    
    wt = H / (rhoa * CP_AIR_MASS) # K m s-1
    wq = LE / (rhoa * Lv) # kg/kg m s-1
    wc = Fc / Mair # umol/mol m s-1 = ppm m s-1
    
    return wt, wq, wc

def _initialize_results(variables, Nstep, Nsoil_nodes, Ncanopy_nodes, Nplant_types, Nasl_nodes):
    """
    Creates temporary results dictionary to accumulate simulation results
    SL 12.11.2019: removed if 'date' in dimensions and added option to save planttype profiles
    """

    results = {}

    for var in variables:

        var_name = var[0]
        dimensions = var[2]

        if 'canopy' in dimensions:
            if 'planttype' in dimensions:
                var_shape = [Nstep, Nplant_types, Ncanopy_nodes]            
            else:
                var_shape = [Nstep, Ncanopy_nodes]

        elif 'soil' in dimensions:
            var_shape = [Nstep, Nsoil_nodes]

        elif 'planttype' in dimensions and 'canopy' not in dimensions:
            var_shape = [Nstep, Nplant_types]

        else:
            var_shape = [Nstep]
        
        results[var_name] = np.full(var_shape, np.NAN)
        
        # append asl
        # print(var_name, var_shape, dimensions)
 
    return results


def _append_results(group, step, step_results, results):
    """
    Adds results from each simulation steps to temporary results dictionary
    """

    results_keys = results.keys()
    step_results_keys = step_results.keys()

    for key in step_results_keys:
        variable = group + '_' + key
        if variable in results_keys:
            if key == 'z' or key == 'planttypes':
                results[variable] = step_results[key]
        
            else:
                results[variable][step] = step_results[key]

    return results

#if __name__ == 'main':
#    import logging
#    from parameters.general import logging_configuration
#    from logging.config import dictConfig
#    dictConfig(logging_configuration)

#    driver(create_ncf=True, dbhfile='letto2014.txt')

