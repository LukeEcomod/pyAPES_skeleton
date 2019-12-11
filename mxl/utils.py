# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 19:15:07 2019

@author: slauniai
"""
import os
import numpy as np
import pandas as pd

VON_KARMAN = 0.41  # von Kárman constant (-)
GRAV_CONST = 9.81  # ms-2 acceleration due gravity
NT = 273.15  # 0 degC in Kelvin
NP = 101300.0  # Pa, sea level normal pressure
R = 8.314462175  # J mol-1 K-1, universal gas constant. One gets specific gas constant by R/M where M is molar mass
CP_AIR_MOLAR = 29.3  # J mol-1 K-1 molar heat capacity of air at constant pressure
CP_AIR_MASS = 1004.67  # J kg-1 K-1 heat capasity of the air at constant pressure
MAIR_DRY = 28.964e-3  # kg mol-1, molar mass of dry air
MH2O = 18.02e-3  # kg mol-1, molar mass of H2O
MCO2 = 44.01e-3  # kg mol-1, molar mass of CO2
L_MOLAR = 44100.0 # J mol-1, latent heat of vaporization at 20degC


# ---- data input

def read_mxl_forcing(ffile, start_time=None, end_time=None, dt0=1800.0, dt=1800.0, na_values='NaN'):
    """
    reads 30 min forcing data and returns dataframe
    
    cols = ['year','month', 'day', 'doy', 'hour', 'minute', 'U', 'ust', 'Ta', 'RH', 'CO2', 'H2O', 'O3',
            'Prec', 'P', 'dirPar', 'diffPar', 'dirNir', 'diffNir', 'Rnet', 'LWin', 'LWout',
            'LWnet', 'Tsh', 'Tsa','Tsc', 'Wh', 'emiatm', 'cloudfract', 'Rew', 'Psi_s', 'Zen', 'Azim', 'Daylength',
            'Ws', 'Tdaily', 'X', 'DDsum', 'H', 'LE', 'NEE', 'qc_H', 'qc_LE', 'qc_NEE']
    Returns:
        dataframe with columns ['ust', 'wt', 'wq', 'wc']]
    """
    
    dat = pd.read_csv(ffile, sep=';', header='infer')
    tvec = pd.to_datetime(dat[['year', 'month', 'day', 'hour', 'minute']])
    dat.index = tvec
    
    dat['doy'] = dat['doy'].astype(float)
    dat['Prec'] = dat['Prec'] / dt0 # mm s-1
    dat['P'] *= 1e3 # Pa
    dat['H2O'] *= 1e-3 # mol/mol
    
    # --- get period 
    if start_time == None:
        start_time = dat.index[0]
    if end_time == None:
        end_time = dat.index[-1]
    
    dat = dat[start_time:end_time]
    
    # -- convert surface fluxes to kinematic fluxes
    T = dat.Ta.values # degC
    P = dat.P.values # Pa
    h2o = dat.H2O.values * P # Pa
    Pdry = P - h2o  # Pa pressure of dry air

    rhoa = (Pdry*MAIR_DRY + h2o*MH2O) / (R*(T+NT))  # kg m-3
    Lv = 1.0e6*(2.501 - 2.361e-3*T)  # J kg-1    
    Mair = P / (R *(T+NT)) # m3 mol-1
    #print(np.mean(Mair), np.mean(rhoa))
    
    dat['wt'] = dat.H.values / (rhoa * CP_AIR_MASS) # K m s-1
    dat['wq'] = dat.LE / (rhoa * Lv) # kg/kg m s-1
    dat['wc'] = dat.NEE / Mair # umol/mol m s-1 = ppm m s-1
    
    #dat = dat[['doy','Ta', 'P', 'ust', 'wt', 'wq', 'wc']]
    
    # interpolate to dt
    if dt != dt0:
        step = '%dS' %dt
        sampler = dat.resample(step)
        dat = sampler.interpolate(method='linear')

    return dat


def initialize_netcdf(sim,
                      tindex,
                      filepath='results/',
                      filename='mxl2.nc',
                      description='Mxl results'):

    """ MXL netCDF4 format output file initialization
    Args:
        sim (int): number of simulations
        tindex (pd.datetimeindex)
        forcing: forcing data (pd.dataframe)
        filepath: path for saving results
        filename: filename
    """
    from netCDF4 import Dataset, date2num
    from datetime import datetime

    folder = os.getcwd()
    filepath = os.path.join(folder, filepath)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    ff = os.path.join(filepath, filename)
    
    tvec = [pd.to_datetime(k) for k in tindex]
    
    # variables and their units
    variables = [['h', 'mxl height [m]'],
                 ['h_lcl', 'lcl height [m]'],
                 ['T_lcl', 'lcl temperature [K]'],
                 ['theta', 'mxl potential temperature [K]'],
                 ['q', 'mxl specific humidity [kg/kg]'],
                 ['thetav', 'mxl virtual potential temperature [K]'],
                 ['ca', ' mxl CO2 mixing ratio [ppm]'],
                 ['Ws', 'subsidene velocity [ms-1]'],
                 ['wstar', 'convective velocity scale [ms-1]'],
                 ['sigmaw', 'turbulent velocity scale [ms-1]'],
                 ['u', ' mxl horizontal wind speed [ms-1]'],
                 ['U', 'mxl wind speed [ms-1]'],
                 ['vpd', 'surface vapor pressure deficit [kPa]'],
                 ['rh', 'surface relative humidity [-]'],
                 ['Psurf', 'surface pressure [kPa]'],

                 # entrainment zone
                 ['Pe', 'entrainment zone pressure [kPa]'],
                 ['Te', 'entrainment zone temperature [K]'],
                 ['theta_jump', 'potential temperature jump [K]'],
                 ['q_jump', 'specific humidity jump [kg/kg]'],
                 ['ca_jump', 'CO2 jump [ppm]'],
                 ['thetav_jump', 'virtual potential temperature jump [K]'],

                 # surface forcing to mxl
                 ['wthetas', 'surface kinematic heat flux [Kms-1]'],
                 ['wqs', 'surface kinematic moisture flux [kg kg-1 ms-1]'],
                 ['wcs', 'surface kinematic CO2 flux [ppm ms-1]']
                ]

    # create dataset and dimensions
    ncf = Dataset(ff, 'w')
    ncf.description = description
    ncf.history = 'created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    ncf.source = 'MXL v270319'

    ncf.createDimension('date', None)
    ncf.createDimension('simulation', sim)
    

    # time variable
    time = ncf.createVariable('date', 'f8', ('date'))
    time.units = 'days since 0001-01-01 00:00:00.0'
    time.calendar = 'standard'

    #tvec = [pd.to_datetime(k) for k in forcing.index]
    time[:] = date2num(tvec, units=time.units, calendar=time.calendar)

    for var in variables:

        var_name = var[0]
        var_unit = var[1]
        #var_dim = var[2]

        variable = ncf.createVariable(var_name, 'f4', ('date', 'simulation'))

        variable.units = var_unit

    return ncf, ff

def latent_heat_vaporization(T, units="mass"):
    """
    Latent heat of vaporization of water.
    Args:
        T - temperature (degC)
        units - output units, "mass" = J kg-1 , "molar"= J mol-1
    Returns:
        Lv - latent heat of vaporization in 'units'
    """
    if np.any(T > 200):
        T = T - NT  #T must be in degC
    Lv = 1.0e6*(2.501 - 2.361e-3*T)*MH2O  # J mol-1    
    if units=="mass":
        Lv = Lv / MH2O  # J kg-1
    return Lv

def air_density(T, P=101300.0, h2o=0.0, units="mass"):
    """
    Computes air density at temperature T, pressure P and vapor pressure H2O
    INPUT:
        T - air temperature (degC, scalar or array
        P - ambient pressure (Pa), scalar or array, optional
        H2O - water vapor partial pressure (Pa), scalar or array, optional (default = dry air)
        units - units to return density: "mass" (default), "molar")
    OUTPUT:
        rhoa - density of dry (default) or moist air (kg m-3 or mol m-3), scalar or array
    Samuli Launiainen 28.4.2014
    """
    if max(T) < 200:
        T = T + NT  # K

    # partial pressures of ideal gas are additive
    Pdry = P - h2o  # pressure of dry air

    if units == "mass":
        rhoa = (Pdry*MAIR_DRY + h2o*MH2O) / (R*T)  # kg m-3

    elif units == "molar":
        rho_d = Pdry / (R*T)  # dry air, mol m-3
        rho_v = h2o / (R*T)  # water vapor, mol m-3
        rhoa = rho_d + rho_v

    else:
        rhoa = np.nan
    return rhoa

#%%
#import datetime
#import netCDF4
#
#times = [datetime.datetime(2016, 10, 1) + datetime.timedelta(hours=hour)
#         for hour in range(84)]
#
## Create netCDF file
#calendar = 'standard'
#units = 'days since 1970-01-01 00:00'
#ds = netCDF4.Dataset('test.nc', 'w')
#timedim = ds.createDimension(dimname='time', size=len(times))
#
## Write timestamps to netCDF file using 32bit float
#timevar32 = ds.createVariable(varname='time32', dimensions=('time',),
#                              datatype='float32')
#timevar32[:] = netCDF4.date2num(times, units=units, calendar=calendar)
#
## Write timestamps to netCDF file using 64bit float
#timevar64 = ds.createVariable(varname='time64', dimensions=('time',),
#                              datatype='float64')
#timevar64[:] = netCDF4.date2num(times, units=units, calendar=calendar)
#
## Read timestamps from netCDF file
#times32 = netCDF4.num2date(timevar32[:], units=units, calendar=calendar)
#times64 = netCDF4.num2date(timevar64[:], units=units, calendar=calendar)
#for time, time32, time64 in zip(times, times32, times64):
#    print("original  ", time)
#    print("  32 bit  ", time32)
#    print("  64 bit  ", time64)
#    print("---")
                      