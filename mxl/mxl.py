# -*- coding: utf-8 -*-
"""
Created on Mon May 22 16:43:11 2017

@author: slauniai
"""

import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

#: machine epsilon
EPS = np.finfo(float).eps

#: [-], von Karman constant
VON_KARMAN = 0.41
#: [K], zero degrees celsius in Kelvin
DEG_TO_KELVIN = 273.15
#: [W m\ :sup:`-2` K\ :sup:`-4`\ ], Stefan-Boltzmann constant
STEFAN_BOLTZMANN = 5.6697e-8
#: [kg m\ :sup:`2` s\ :sup:`-1`\ ], standard gravity
GRAVITY = 9.81
#: [J mol\ :sup:`-1` K\ :sup:``-1], universal gas constant; One gets specific gas constant by R/M where M is molar mass
GAS_CONSTANT = 8.314
#: [kg mol\ :sup:`-1`\ ], molar mass of air
MOLAR_MASS_AIR = 29.0e-3
#: [kg mol\ :sup:`-1`\ ], molar mass of H\ :sub:`2`\ O
MOLAR_MASS_H2O = 18.015e-3
#: [kg mol\ :sup:`-1`\ ], molar mass of CO\ :sub:`2`\
MOLAR_MASS_CO2 = 44.01e-3
#: [m\ :sup:`2` s\ :sup:`-1`\ ], kinematic viscosity of air at 20\ :math:`^{\circ}`\ C
AIR_VISCOSITY = 15.1e-6
#: [m\ :sup:`2` s\ :sup:`-1`\ ], thermal diffusivity of air at 20\ :math:`^{\circ}`\ C
THERMAL_DIFFUSIVITY_AIR = 21.4e-6
#: [m\ :sup:`2` s\ :sup:`-1`\ ], molecular diffusvity of CO\ :sub:`2` at 20\ :math:`^{\circ}`\ C
MOLECULAR_DIFFUSIVITY_CO2 = 15.7e-6
#: [m\ :sup:`2` s\ :sup:`-1`\ ], molecular diffusvity of H\ :sub:`2`\ at 20\ :math:`^{\circ}`\ C
MOLECULAR_DIFFUSIVITY_H2O = 24.0e-6
#: [Pa], sea level normal pressure
NP = 101300.0
#: [J mol\ :sup:`-1` K\ :sup:`-1`\ ], heat capacity of air at constant pressure
SPECIFIC_HEAT_AIR = 29.3
#: [J kg\ :sup:`-1` K\ :sup:``-1] gravimetric heat capacity of air at constant pressure
CP_AIR_MASS = 1004.67
#: [J mol\ :sup:`-1`\ ], latent heat of vaporization at 20\ :math:`^{\circ}`\ C
LATENT_HEAT = 44100.0

class MXLmodel():
    """
    Simple implementation of Mixed Boundary Layer (MXL) model based on:
        Janssens & Pozzer, 2015. GeoSci Model Dev. 8, 453 - 471.
        Vilà-Guerau de Arellano et al. 2015. Atmospheric Boundary Layer: 
            Integrating Air Chemistry and Land Interactions. Cambridge University Press, 
            New York, 2015, 265 pp
        Stull, 1998.
        Siqueira et al. 2009. J. Hydrometeorol. 10, 96-112.
    
    Samuli Launiainen Luke, 20.9.2018
    """

    def __init__(self, ini, params):
        """
        Args:
            ini - initial conditions (dict)
            params - mxl model parameters (dict)
        Returns:
            MXLmodel instance
        """
        # --- parameters
        self.deltat = params['dt']                      # forcing timestep
        self._nsteps = int(np.ceil(self.deltat / 60.0))
        self._dt = self.deltat / self._nsteps           # internal timestep s
        self._steps = int(self.deltat / self._dt)       # number of subtimesteps
        self._beta = params['beta']                     # ratio of entrainment
                                                        # boyancy flux to surface boyancy flux
        self._divU = params['divU']                     # s-1, large-scale horizontal wind divergence
        self._f = params['f']                           # coriolis parameter, s-1
        self.ctr = params['ctr']                        # control structures for model [dict]
        
        # --- mixed layer state variables
        self.h = ini['h']                               # initial height, m
        self.h_lcl = None                               # lifting condensation level, m
        self.T_lcl = None                               # lifting condensation level T, m
        self.theta = ini['theta']                       # potential temperature, K
        self.q = ini['q']                               # specific humidity, kg /kg
        self.thetav = self.theta*(1.0 + 0.61*self.q)    # virtual potential temperature
        self.ca = ini['ca']                             # CO2 mixing ratio (ppm)

        # --- computed from mxl state
        self.rhoa = 1.2                                 # air density at surface (kg m-3)
        self.vpd = None                                 # vapor pressure deficit (kPa)
        self.rh = None                                  # vapor pressure (kPa)
        self.esat = None                                # saturation vapor pressure (kPa)
        self.Psurf = ini['Psurf']                       # surface pressure (kPa)
        self.Pe = None                                  # pressure at mixed layer top (kPa)
        self.Te = None

        # --- free troposphere lapse rates
        self.gamma_theta = ini['gamma_theta']   # lapse rate K m-1
        self.gamma_q = ini['gamma_q']           # kg/kg m-1
        self.gamma_ca = ini['gamma_ca']         # ppm m-1   
        
        # --- scalar jumps at entrainment layer of infinitesimal height
        self.theta_jump = ini['theta_jump']     # K        
        self.q_jump = ini['q_jump']             # kg/kg                
        self.ca_jump = ini['ca_jump']           # ppm
        
        self.thetav_jump = self.theta_jump + 0.61*(self.q*self.theta_jump \
                           + self.theta*self.q_jump + self.theta_jump*self.q_jump) # K

        # --- mixed layer height & scalar tendencies        
        self.h_tend = 0.0       # m s-1
        self.we = None          # entrainment velocity m s-1
        self.theta_tend = None  # K s-1
        self.q_tend = None      # kg/kg s-1
        self.ca_tend = None     # ppm s-1
        
        self.theta_jump_tend = None # K s-1
        self.q_jump_tend = None     # kg kg-1 s-1
        self.ca_jump_tend = None    # ppm s-1

        # --- velocity scales
        self.Ws = -self._divU * self.h  # large-scale subsidence velocity [m s-1]
        self.wstar = 1e-6               # convective velocity scale m s-1
        self.sigmaw = None              # turbulent velocity scale m s-1
        
        # --- horizontal wind velocity: we care only on magnitude, not components
        if self.ctr['Wind']:
            self.u = ini['u'] + EPS
            #self.v = ini['v'] + EPS
            self.U = np.sqrt(self.u**2. + self.wstar**2.)
            self.u_jump = ini['u_jump']
            #self.v_jump = ini['v_jump']
            self.gamma_u = ini['gamma_u'] # s-1
            #self.gamma_v = ini['gamma_v']
            

    def run_timestep(self, wthetas, wqs, wcs, ustar, out=False):
        """
        grows mixed layer and computes scalar state variables for one timestep
        IN:
            dt - gross timestep
            wthetas - kinematic surface heat flux (K m s-1)
            wqs - kinematic moisture flux (kg kg-1 m s-1)
            wcs - kinematic co2 flux (ppm m s-1)
        OUT:
            updated object
            if out = True, returns mixed-layer state variables
                self.theta
                self.rh
                self.ca
                self.U
        """
        # integrate to t + self.deltat
        for jj in range(self._nsteps):

            # --- compute mixed layer tendencies
            self.tendencies(wthetas, wqs, wcs)
            
            if self.ctr['Wind']:
                self.compute_windspeed(ustar)

            # --- integrate mixed layer 
            self.time_integrate()
        
        if out:
            return self.theta, self.rh, self.ca, self.U
        
    def tendencies(self, wthetas, wqs, wcs):
        """
        Computes mxl growth rate and scalar tendencies (dx/dt) based on MXL 
        and inversion state using surface fluxes from current timestep.
        Args:
            wthetas - surface kinematic heat flux (K m s-1)
            wqs - surface moisture flux (kg kg-1 m s-1)
            wcs - surface co2 flux (ppm s-1)
        """
        self.Ws = -self._divU * self.h # subsidence velocity m s-1, <0
        self.thetav = self.theta*(1. + 0.61*self.q) #  virtual pot. temperature K     
        wthetavs = wthetas + 0.61*self.theta*wqs # K m s-1 boyancy flux at surface
        
        # convective velocity scale
        if wthetavs > 0:
            self.wstar = (GRAVITY / self.thetav * self.h*wthetavs)**(1./3.)
        else:
            self.wstar = 1e-6
        
        # --- entrainment zone  
        wthetave = -self._beta * wthetavs # entrainment boyancy flux K m s-1
        self.we = -wthetave / self.thetav_jump  # entrainment velocity m s-1
        if self.we < 0:
            self.we = 0.0

        wthetae = -self.we * self.theta_jump
        wqe = -self.we * self.q_jump  # kg/kg s-1
        wce = -self.we * self.ca_jump # ppm s-1
        
        # --- mxl growth rate and scalar tendencies
        self.h_tend = self.we + self.Ws  # mixed layer height growth (m s-1) by entrainment and large-scale subsidence
        
        self.theta_tend = (wthetas - wthetae) / self.h  #+ self.adv_theta # K s-1
        self.q_tend = (wqs - wqe) / self.h  #+ self.adv_q # kgkg-1 s-1
        self.ca_tend = (wcs - wce) / self.h #+ self.adv_c # ppm s-1
        
        self.theta_jump_tend = self.gamma_theta *self.we - self.theta_tend
        self.q_jump_tend = self.gamma_q * self.we - self.q_tend
        self.ca_jump_tend = self.gamma_ca * self.we - self.ca_tend
        
    def time_integrate(self):
        """
        integrate mixed layer state to next time level
        """
        
        self.h += self.h_tend * self._dt
        
        self.theta += self.theta_tend * self._dt
        self.q += self.q_tend * self._dt
        self.ca += self.ca_tend * self._dt
        self.thetav = self.theta*(1. + 0.61*self.q)
        
        self.theta_jump += self.theta_jump_tend * self._dt
        self.q_jump += self.q_jump_tend * self._dt
        self.ca_jump += self.ca_jump_tend * self._dt
        
        self.thetav_jump = self.theta_jump + 0.61*(self.q*self.theta_jump \
                            + self.theta*self.q_jump + self.theta_jump*self.q_jump)
        
        # --- state variables at mxl bottom
        self.esat, qs = e_sat(self.theta, self.Psurf)
        self.rh = self.q / qs
        self.vpd = (1. - self.rh)*self.esat
        
        # --- compute lifting condensation level variables
        self.lcl()
        
        # --- pressure at surface
        self.rhoa = air_density_from_elev(self.theta, 0.0, self.Psurf)
        
        # --- pressure, temperature and rh at entrainment layer
        self.Pe = pressure_from_elev(self.h, self.Psurf)
        self.Te = self.theta - GRAVITY / CP_AIR_MASS * self.h  # K
        _, qs = e_sat(self.Te, self.Pe)
        self.rhe = self.q / qs
        del qs
        
        
    def compute_windspeed(self, ustar):
        """
        wind speed in mixed layer. We don't care about components but only on
        mean horizontal wind speed
        """
        u_tend = (-ustar**2. + self.we*self.u_jump) / self.h  # m

        u_jump_tend = self.gamma_u * self.we - u_tend # m
        
        # time integrate
        self.u += u_tend * self._dt # m s-1
        if self.u < 0:
            self.u = EPS
        self.u_jump += u_jump_tend * self._dt # m s-1
        if self.u_jump <0:
            self.u_jump = EPS
        # mean velocity        
        self.U = np.sqrt(self.u**2. + self.wstar**2.)
        

    def lcl(self):
        """
        Calculates lifting condensation level and saturation point temperature and
        pressure. 
        Modified from Gaby's notes
    
        Updates:
            h_lcl - lifting condensation level, m
            T_lcl -  --"-- temperature, K
        """
        HZ = GAS_CONSTANT * self.theta / (MOLAR_MASS_AIR * GRAVITY)  # scaling height, m
        ea =  self.rh * self.esat
        r = 0.622 * ea / (self.Psurf - ea)
        
        self.T_lcl = 2840. / (3.5*np.log(self.theta) - np.log(self.Psurf*r / (0.622 + r)) - 7.108) + 55.
        
        p_lcl = self.Psurf * (self.T_lcl / self.theta)**3.5
        self.h_lcl = HZ*(-np.log(p_lcl / self.Psurf))
        
        
""" ---- utility functions --- """

def e_sat(T, Ps=101.3):
    """
    Saturation vapor pressure and saturation specific humidity
    Args:
        T (K), Ps (kPa)
    Returns
        es (kPa), qs (kg kg-1)
    """
    es = 0.611 * np.exp(17.2694 * (T - 273.16) / (T - 35.86))  # kPa
    qs = 0.622 * es / Ps # kg kg-1
    return es, qs

def air_density_from_elev(T, Elev, Ps=101.3):
    """
    Args: T [k], Elev[m], Ps [kPa]
    Returns: rho_air [kgm-3]
    """
    p = Ps * np.exp(-Elev / 8200.0)       # kPa
    rho_air = 1000 * p * MOLAR_MASS_AIR / (GAS_CONSTANT * T)     # kg m-3 
   
    return rho_air

def pressure_from_elev(Elev, Ps=101.3):
    """
    Args: Elev [m], Ps [kPa]. From sea level OR reference to local surface
    """
    p = Ps*np.exp(-Elev / 8200.0)       # kPa
    return p

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
        T = T - DEG_TO_KELVIN  #T must be in degC
    Lv = 1.0e6 * (2.501 - 2.361e-3*T) * MOLAR_MASS_H2O  # J mol-1    
    if units=="mass":
        Lv = Lv / MOLAR_MASS_H2O  # J kg-1
    return Lv

def equilibrium_height(ustar, wthetas, thetas, f=1e-4):
    """
    Steady-state equilibrium height of stable boundary layer
    Zilitinkevich, 1972. BLM 3, 141-145.
    IN:
        ustar - ms-1, mean nocturnal friction velocity
        wthetas - K ms-1, mean nocturnal kinematic heat flux
        wtheta - K, mean nocturnal potential temperature
        f - coriolis parameter (s-1)
    OUT:
        z - equilibrium height, m
        L - Obukhov length
    """
    GAMMA_C = 0.4  # similarity constant
    L = ustar**3 / (VON_KARMAN * GRAVITY * wthetas / thetas)
    z = GAMMA_C * abs(ustar*L / f)**0.5
    return z, L

#    def velocity_scales(self, ustar, out=False):
#        """
#        updates convective (wstar) and turbulent velocity (sigmaw) scales
#        """
#        if self.thetavs > 0:
#            self.wstar = (GRAV_CONST / self.thetav * self.h*self.wthetavs)**(1./3.)
#        else:
#            self.wstar = 1e-6
#            
#        # self.sigmaw = (self.wstar**3. + (self._A / self._cf)*ustar**3.)**(1./3.)     
#
#        if out:
#            return self.wstar, self.sigmaw
