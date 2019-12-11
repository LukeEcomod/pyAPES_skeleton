# -*- coding: utf-8 -*-
"""
.. module: micromet
    :synopsis: APES-model component
.. moduleauthor:: Kersti Haahti

Describes turbulent flow and scalar profiles in canopy - surface layer regime  using
1st order closure scheme.

modified version of 'micromet' to handle varying vertical grid in growing ASL.
LAST CHANGES: 04.04.2019 / SL

References:
Launiainen, S., Katul, G.G., Lauren, A. and Kolari, P., 2015. Coupling boreal
forest CO2, H2O and energy flows by a vertically structured forest canopy – 
Soil model with separate bryophyte layer. Ecological modelling, 312, pp.385-405.
"""
import numpy as np
import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)

from mxl.mxl import EPS, SPECIFIC_HEAT_AIR, VON_KARMAN, MOLECULAR_DIFFUSIVITY_CO2, \
                       MOLECULAR_DIFFUSIVITY_H2O, THERMAL_DIFFUSIVITY_AIR, \
                       AIR_VISCOSITY, MOLAR_MASS_AIR, GRAVITY

class Micromet(object):
    r""" Computes flow and scalar profiles within canopy.
    """
    def __init__(self, z, lad, hc, p):
        r""" Initializes object for computation of flow and scalar profiles within canopy.

        Args:
            z (array): canopy model nodes; equidistant; height from soil surface (= 0.0) [m]
            lad (array): leaf area density [m2 m-3]
            hc (float): canopy heigth [m]
            p (dict):
                'zos': forest floor roughness length [m]
                'dPdx': horizontal pressure gradient
                'Cd': drag coefficient
                'Utop': Upper boundary [m/s]
                'Ubot': lower boundary [m/s]
                'Sc' (dict): {'T','H2O','CO2'} Schmidt numbers
        Returns:
            self (object)
        NOTE: this assumes now that Utop and Ubot given to model are dimensionless [ms-1]
        """

        # parameters
        #self.zos = p['zos']  # ground roughness length [m]
        self.Cd = p['Cd']  # drag coefficient
        self.Utop = p['Utop']  # wind speed at top boundary [ms-1]
        self.Ubot = p['Ubot']  # wind speed at lwer boundary [ms-1]
        self.Sc = p['Sc']  # Schmidt numbers
        self.dPdx = p['dPdx']  # horizontal pressure gradient
        
        # grid       
        self.z = z
        self.dz = z[1] - z[0] # step: must be continuous
        self.N = len(z) # number of nodes belonging to canopy model.
        
        # state variables
        self.U = None # mean velocity profile
        self.Km = None # eddy diffusivity
        self.tau = None # reynolds stress profile
        self.ust = None # friction velocity profile
        
        self.d = None # zero-place displacement height (m)
        self.zo = None # roughness length for momentum (m)

        # initialize state variables
        self.tau, self.U, self.Km, self.l_mix, self.d, self.zo = closure_1_model_U(
                z, self.Cd, lad, hc, self.Utop + EPS, self.Ubot + EPS, dPdx=self.dPdx)
        
        self.ust = np.abs(self.tau)**0.5
        
    def flow_stats(self, z, lad, hc, Utop, Ubot=None):
        r""" Computes mean velocity profile, shear stress and
        eddy diffusivity within and above horizontally homogenous plant
        canopies using 1st order closure.

        Args:
            z (array): nodes, height from soil surface (= 0.0) [m]
            lad (array): leaf area density [m2 m-3]
            hc (float): canopy heigth [m]
            Utop (float): U/ustar [-]
            ust0 (float): ustar at top; give if Utop & Ubot are normalized
        """
        
        # pad lad with zeros to account for growing z
        n = len(z)
        m = len(lad)
        if n > m:
            lad = np.pad(lad, (0,n - m), mode='constant', constant_values=[0])
        del n, m
        
        if Ubot:
            self.Ubot = Ubot
            
        tau, U, Km, l_mix, _, _ = closure_1_model_U(
                z, self.Cd, lad, hc, Utop + EPS, self.Ubot, dPdx=self.dPdx, U_ini=None)

        if any(U < 0.0):
            logger.debug('Negative U !!!')

        else:
            self.z = z.copy()
            self.l_mix = l_mix.copy()
            self.U = U.copy()
            self.Km = Km.copy()
            self.tau = tau.copy()
            self.ust = np.abs(self.tau)**0.5
            
        # outputs to canopy model:
        return self.U[0:self.N], self.ust[0:self.N]
            

    def scalar_profiles(self, gam, H2O, CO2, T, P, source, lbc, ubc, Ebal):
        r""" Solves scalar profiles (H2O, CO2 and T) in canopy-asl continuum
        Input array len equals canopy model grid; extended to whole asl here!
        
        Args:
            gam (float): weight for new value in iterations
            H2O (array): water vapor mixing ratio [mol mol-1]
            CO2 (array): carbon dioxide mixing ratio [ppm]
            T (array): ambient air temperature [degC]
            P: ambient pressure [Pa]
            source (dict):
                'H2O' (array): water vapor source [mol m-3 s-1]
                'CO2' (array): carbon dioxide source [umol m-3 s-1]
                'T' (array): heat source [W m-3]
            lbc (dict):
                'H2O' (float): water vapor lower boundary [mol m-2 s-1]
                'CO2' (float): carbon dioxide lower boundary [umol m-2 s-1]
                'T' (float): heat lower boundary [W m-2]
            ubc (dict): These are to match solution to mxl-value at each timestep
                'H2O' (float): water vapor lower boundary [mol/mol]
                'CO2' (float): carbon dioxide lower boundary [ppm]
                'T' (float): heat lower boundary [T]                
        Returns:
            H2O (array): water vapor mixing ratio [mol mol-1]
            CO2 (array): carbon dioxide mixing ratio [ppm]
            T (array): ambient air temperature [degC]
            err_h2o, err_co2, err_t (floats): maximum error for each scalar
        """
        
        # SL 04.04.2019: for growing abl, we need compute profiles from
        # soil surface to top of abl. Thus, pad constant canopy model grid with zeros
        
        source0 = source.copy() # make local copy as changes here were seen in canopy_asl !!
        
        T = T.copy()
        H2O = H2O.copy()
        CO2 = CO2.copy()
        
        n = len(H2O) # nodes in canopy
        m = len(self.z) # nodes in canopy + asl
        if n < m:
            T = np.pad(T, (0, m - n), mode='linear_ramp', end_values=[ubc['T']])
            H2O = np.pad(H2O, (0, m - n), mode='linear_ramp', end_values=[ubc['H2O']])
            CO2 = np.pad(CO2, (0, m - n), mode='linear_ramp', end_values=[ubc['CO2']])
            
            for key in source:
                source0[key] = np.pad(source0[key], (0, m - n), mode='constant', constant_values=[0])
        
        # previous guess, not values of previous time step!
        H2O_prev = H2O.copy()
        CO2_prev = CO2.copy()
        T_prev = T.copy()

        # --- H2O ---
        H2O = closure_1_model_scalar(dz=self.dz,
                                     Ks=self.Km * self.Sc['H2O'],
                                     source=source0['h2o'],
                                     ubc=H2O[-1],
                                     lbc=lbc['H2O'],
                                     scalar='H2O',
                                     T=T[-1], P=P)
        # new H2O
        H2O = (1 - gam) * H2O_prev + gam * H2O
        # limit change to +/- 10%
        if all(~np.isnan(H2O)):
            H2O[H2O > H2O_prev] = np.minimum(H2O_prev[H2O > H2O_prev] * 1.1, H2O[H2O > H2O_prev])
            H2O[H2O < H2O_prev] = np.maximum(H2O_prev[H2O < H2O_prev] * 0.9, H2O[H2O < H2O_prev])
        # relative error
        err_h2o = max(abs((H2O[0:n] - H2O_prev[0:n]) / H2O_prev[0:n]))

        # --- CO2 ---
        CO2 = closure_1_model_scalar(dz=self.dz,
                                     Ks=self.Km * self.Sc['CO2'],
                                     source=source0['co2'],
                                     ubc=CO2[-1],
                                     lbc=lbc['CO2'],
                                     scalar='CO2',
                                     T=T[-1], P=P)
        # new CO2
        CO2 = (1 - gam) * CO2_prev + gam * CO2
        
        # limit change to +/- 10%
        if all(~np.isnan(CO2)):
            CO2[CO2 > CO2_prev] = np.minimum(CO2_prev[CO2 > CO2_prev] * 1.1, CO2[CO2 > CO2_prev])
            CO2[CO2 < CO2_prev] = np.maximum(CO2_prev[CO2 < CO2_prev] * 0.9, CO2[CO2 < CO2_prev])
        # relative error
        err_co2 = max(abs((CO2[0:n] - CO2_prev[0:n]) / CO2_prev[0:n]))

        if Ebal:
            # --- T ---
            T = closure_1_model_scalar(dz=self.dz,
                                       Ks=self.Km * self.Sc['T'],
                                       source=source0['sensible_heat'],
                                       ubc=T[-1],
                                       lbc=lbc['T'],
                                       scalar='T',
                                       T=T[-1], P=P)
            # new T
            T = (1 - gam) * T_prev + gam * T
            # limit change to T_prev +/- 2degC
            if all(~np.isnan(T)):
                T[T > T_prev] = np.minimum(T_prev[T > T_prev] + 2.0, T[T > T_prev])
                T[T < T_prev] = np.maximum(T_prev[T < T_prev] - 2.0, T[T < T_prev])

            # absolut error
            err_t = max(abs(T[0:n] - T_prev[0:n]))
        else:
            err_t = 0.0
        
        # return profiles from ground to top of ASL
        asl_profs = {'T': T, 'H2O': H2O, 'CO2': CO2}
        
        # in canopy model layers
        T = T[0:n]
        H2O = H2O[0:n]
        CO2 = CO2[0:n]
        
        return H2O, CO2, T, err_h2o, err_co2, err_t, asl_profs

def closure_1_model_U(z, Cd, lad, hc, Utop, Ubot, dPdx=0.0, lbc_flux=None, U_ini=None):
    """
    Computes mean velocity profile, shear stress and eddy diffusivity
    within and above horizontally homogenous plant canopies using 1st order closure.
    Accounts for horizontal pressure gradient force dPdx, assumes neutral diabatic stability.
    Solves displacement height as centroid of drag force.
    IN:
       z - height (m), constant increments
       Cd - drag coefficient (0.1-0.3) (-)
       lad - plant area density, 1-sided (m2m-3)
       hc - canopy height (m)
       Utop - U or U /u* upper boundary
       Uhi - U or U /u* at ground (0.0 for no-slip)
       dPdx - u* -normalized horizontal pressure gradient
    OUT:
       tau - reynolds stress ([m2s-2] or[ms-1]) 
       U - mean wind speed (m/s or -) 
       Km - eddy diffusivity for momentum ([m2s-1] or [m]) 
       l_mix - mixing length (m)
       d - zero-plane displacement height (m)
       zo - roughness lenght for momentum (m)
    
    NOTE:
        inputs Utop and Ubot can be either dimensional [ms-1] or U/ust [-]
        this will change only units in outputs.
    CODE:
        Gaby Katul, Samuli Launiainen 2008-2015. Converted to Python 17.5.2017
    """
    lad = 0.5*lad  # frontal plant-area density is half of one-sided
    dz = z[1] - z[2]
    N = len(z)
    if U_ini is None:
        U = np.linspace(Ubot, Utop, N)
    else:
        U = U_ini.copy()

    nn1 = max(2, np.floor(N/20))  # window for moving average smoothing

    # --- Start iterative solution
    err = 999.9
    iter_max = 20
    eps1 = 0.5
    dPdx_m = 0.0

    iter_no = 0.0

    while err > 0.01 and iter_no < iter_max:
        iter_no += 1
        Fd = Cd*lad*U**2  # drag force
        d = sum(z*Fd) / (sum(Fd) + EPS)  # displacement height
        l_mix = mixing_length(z, hc, d)  # m

        # --- dU/dz (m-1)
        y = central_diff(U, dz)
        # y = smooth(y, nn1)

        # --- eddy diffusivity & shear stress
        Km = l_mix**2*abs(y)
        # Km = smooth (Km, nn1)
        tau = -Km * y

        # ------ Set the elements of the Tri-diagonal Matrix
        a1 = -Km
        a2 = central_diff(-Km, dz)
        a3 = Cd*lad*U

        upd = (a1 / (dz*dz) + a2 / (2*dz))  # upper diagonal
        dia = (-a1*2 / (dz*dz) + a3)  # diagonal
        lod = (a1 / (dz*dz) - a2 / (2*dz))  # subdiagonal
        rhs = np.ones(N) * dPdx  #_m ???

        # upper BC
        upd[-1] = 0.
        dia[-1] = 1.
        lod[-1] = 0.
        rhs[-1] = Utop

        if not lbc_flux:  # --- lower BC, fixed Ubot
            upd[0] = 0.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = Ubot
        else:  # --- lower BC, flux-based
            upd[0] = -1.
            dia[0] = 1.
            lod[0] = 0.
            rhs[0] = 0.  # zero-flux bc
            # how to formulate prescribed flux bc?
            # rhs[0] = lbc_flux

        # --- call tridiagonal solver
        Un = tridiag(lod, dia, upd, rhs)

        err = max(abs(Un - U))

        # --- Use successive relaxations in iterations
        U = eps1*Un + (1.0 - eps1)*U
        dPdx_m = eps1*dPdx + (1.0 - eps1)*dPdx_m  # ???
        if iter_no == iter_max:
            logger.debug('Maximum number of iterations reached: U_n = %.2f, err = %.2f',
                         np.mean(U), err)

    # ---- return values
    tau = tau # shear stress
    zo = (z[-1] - d)*np.exp(-0.4*U[-1])  # roughness length

    y = forward_diff(U, dz)
    Kmr = l_mix**2 * abs(y)  # eddy diffusivity
    Km = smooth(Kmr, nn1)
    Km[0] = Km[1]
    tau = smooth(tau, nn1)

    # --- for testing ----
    plt.figure(201)
    plt.subplot(221); plt.plot(Un, z, 'r-'); plt.title('U')
    plt.subplot(222); plt.plot(y, z, 'b-'); plt.title('dUdz')
    #plt.subplot(223); plt.plot(l_mix, z, 'r-'); plt.title('l mix')
    plt.subplot(223); plt.plot(tau, z, 'r-'); plt.title('tau')
    plt.subplot(224); plt.plot(Km, z, 'r-', Kmr, z, 'b-'); plt.title('Km')    

    return tau, U, Km, l_mix, d, zo

def closure_1_model_scalar(dz, Ks, source, ubc, lbc, scalar,
                           T=20.0, P=101300.0, lbc_dirchlet=False):
    r""" Solves stationary scalar profiles in 1-D grid
    Args:
        dz (float): grid size (m)
        Ks (array): eddy diffusivity [m2 s-1]
        source (array): sink/source term
            CO2 [umolm-3s-1], H2O [mol m-3 s-1], T [Wm-3]
        ubc (float):upper boundary condition
            value: CO2 [ppm], H2O [mol mol-1], T [degC]
        lbc (float): lower boundary condition, flux or value:
            flux: CO2 [umol m-2 s-1], H2O [mol m-2 s-1], T [Wm-2].
            value: CO2 [ppm], H2O [mol mol-1], T [degC]
        scalar (str): 'CO2', 'H2O', 'T'
        T (float/array?): air temperature [degC]
        P (float/array?): pressure [Pa]
        lbc_dirchlet - True for Dirchlet lower boundary
    OUT:
        x (array): mixing ratio profile
            CO2 [ppm], H2O [mol mol-1], T [degC]
    References:
        Juang, J.-Y., Katul, G.G., Siqueira, M.B., Stoy, P.C., McCarthy, H.R., 2008. 
        Investigating a hierarchy of Eulerian closure models for scalar transfer inside 
        forested canopies. Boundary-Layer Meteorology 128, 1–32.
    Code:
        Gaby Katul & Samuli Launiainen, 2009 - 2017
        Kersti: code condensed, discretization simplified (?)
    Note:
        assumes constant dz and Dirchlet upper boundary condition
    """

    dz = float(dz)
    N = len(Ks)
    rho_a = P / (287.05 * (T + 273.15))  # [kg m-3], air density

    CF = rho_a / MOLAR_MASS_AIR  # [mol m-3], molar conc. of air

    Ks = spatial_average(Ks, method='arithmetic')  # length N+1

    if scalar.upper() == 'CO2':  # [umol] -> [mol]
        ubc = 1e-6 * ubc
        source = 1e-6 * source
        lbc = 1e-6 * lbc

    if scalar.upper() == 'T':
        CF = CF * SPECIFIC_HEAT_AIR  # [J m-3 K-1], volumetric heat capacity of air

    # --- Set elements of tridiagonal matrix ---
    a = np.zeros(N)  # sub diagonal
    b = np.zeros(N)  # diagonal
    g = np.zeros(N)  # super diag
    f = np.zeros(N)  # rhs

    # intermediate nodes
    a[1:-1] = Ks[1:-2]
    b[1:-1] = - (Ks[1:-2] + Ks[2:-1])
    g[1:-1] = Ks[2:-1]
    f[1:-1] = - source[1:-1] / CF * dz**2

    # uppermost node, Dirchlet boundary
    a[-1] = 0.0
    b[-1] = 1.0
    g[-1] = 0.0
    f[-1] = ubc

    # lowermost node
    if not lbc_dirchlet:  # flux-based
        a[0] = 0.0
        b[0] = 1.
        g[0] = -1.
        f[0] = (lbc / CF)*dz / (Ks[1] + EPS)

    else:  #  fixed concentration/temperature
        a[0] = 0.0
        b[0] = 1.
        g[0] = 0.0
        f[0] = lbc

    x = tridiag(a, b, g, f)

    if scalar.upper() == 'CO2':  # [mol] -> [umol]
        x = 1e6*x

    return x

def mixing_length(z, h, d, l_min=None):
    """
    computes mixing length: linear above the canopy, constant within and
    decreases linearly close the ground (below z< alpha*h/VON_KARMAN)
    IN:
        z - computation grid, m
        h - canopy height, m
        d - displacement height, m
        l_min - intercept at ground, m
    OUT:
        lmix - mixing length, m
    """
    dz = z[1] - z[0]

    if not l_min:
        l_min = dz / 2.0
    
    alpha = (h - d)*VON_KARMAN / (h + EPS)
    I_F = np.sign(z - h) + 1.0
    l_mix = alpha*h*(1 - I_F / 2) + (I_F / 2) * (VON_KARMAN*(z - d))
    
    sc = (alpha*h) / VON_KARMAN
    ix = np.where(z < sc)
    l_mix[ix] = VON_KARMAN*(z[ix] + dz / 2)
    l_mix = l_mix + l_min

    return l_mix

def leaf_boundary_layer_conductance(u, d, Ta, dT, P=101300.):
    """
    Computes 2-sided leaf boundary layer conductance assuming mixed forced and free
    convection form two parallel pathways for transport through leaf boundary layer.
    INPUT: u - mean velocity (m/s)
           d - characteristic dimension of the leaf (m)
           Ta - ambient temperature (degC)
           dT - leaf-air temperature difference (degC)
           P - pressure(Pa)
    OUTPUT: boundary-layer conductances (mol m-2 s-1)
        gb_h - heat (mol m-2 s-1)
        gb_c- CO2 (mol m-2 s-1)
        gb_v - H2O (mol m-2 s-1)
        r - ratio of free/forced convection
    Reference: Campbell, S.C., and J.M. Norman (1998),
    An introduction to Environmental Biophysics, Springer, 2nd edition, Ch. 7
    Gaby Katul & Samuli Launiainen
    Note: the factor of 1.4 is adopted for outdoor environment, see Campbell and Norman, 1998, p. 89, 101.
    """

    u = np.maximum(u, EPS)

    # print('U', u, 'd', d, 'Ta', Ta, 'P', P)
    factor1 = 1.4*2  # forced conv. both sides, 1.4 is correction for turbulent flow
    factor2 = 1.5  # free conv.; 0.5 comes from cooler surface up or warmer down

    # -- Adjust diffusivity, viscosity, and air density to pressure/temp.
    t_adj = (101300.0 / P)*((Ta + 273.15) / 293.16)**1.75
    Da_v = MOLECULAR_DIFFUSIVITY_H2O*t_adj
    Da_c = MOLECULAR_DIFFUSIVITY_CO2*t_adj
    Da_T = THERMAL_DIFFUSIVITY_AIR*t_adj
    va = AIR_VISCOSITY*t_adj
    rho_air = 44.6*(P / 101300.0)*(273.15 / (Ta + 273.13))  # [mol/m3]

    # ----- Compute the leaf-level dimensionless groups
    Re = u*d / va  # Reynolds number
    Sc_v = va / Da_v  # Schmid numbers for water
    Sc_c = va / Da_c  # Schmid numbers for CO2
    Pr = va / Da_T  # Prandtl number
    Gr = GRAVITY*(d**3)*abs(dT) / (Ta + 273.15) / (va**2)  # Grashoff number

    # ----- aerodynamic conductance for "forced convection"
    gb_T = (0.664*rho_air*Da_T*Re**0.5*(Pr)**0.33) / d  # [mol/m2/s]
    gb_c=(0.664*rho_air*Da_c*Re**0.5*(Sc_c)**0.33) / d  # [mol/m2/s]
    gb_v=(0.664*rho_air*Da_v*Re**0.5*(Sc_v)**0.33) / d  # [mol/m2/s]

    # ----- Compute the aerodynamic conductance for "free convection"
    gbf_T = (0.54*rho_air*Da_T*(Gr*Pr)**0.25) / d  # [mol/m2/s]
    gbf_c = 0.75*gbf_T  # [mol/m2/s]
    gbf_v = 1.09*gbf_T  # [mol/m2/s]

    # --- aerodynamic conductance: "forced convection"+"free convection"
    gb_h = factor1*gb_T + factor2*gbf_T
    gb_c = factor1*gb_c + factor2*gbf_c
    gb_v = factor1*gb_v + factor2*gbf_v
    # gb_o3=factor1*gb_o3+factor2*gbf_o3

    #r = Gr / (Re**2)  # ratio of free/forced convection

    return gb_h, gb_c, gb_v#, r


def e_sat(T):
    """
    Computes saturation vapor pressure (Pa), slope of vapor pressure curve
    [Pa K-1]  and psychrometric constant [Pa K-1]
    IN:
        T - air temperature (degC)
    OUT:
        esa - saturation vapor pressure in Pa
        s - slope of saturation vapor pressure curve (Pa K-1)
    SOURCE:
        Campbell & Norman, 1998. Introduction to Environmental Biophysics. (p.41)
    """

    esa = 611.0 * np.exp((17.502 * T) / (T + 240.97))  # Pa
    s = 17.502 * 240.97 * esa / ((240.97 + T)**2)

    return esa, s


""" --- numerical methods --- """

def forward_diff(y, dx):
    """
    computes gradient dy/dx using forward difference
    assumes dx is constatn
    """
    N = len(y)
    dy = np.ones(N) * np.NaN
    dy[0:-1] = np.diff(y)
    dy[-1] = dy[-2]
    return dy / dx

def central_diff(y, dx):
    """
    computes gradient dy/dx with central difference method
    assumes dx is constant
    """
    N = len(y)
    dydx = np.ones(N) * np.NaN
    # -- use central difference for estimating derivatives
    dydx[1:-1] = (y[2:] - y[0:-2]) / (2 * dx)
    # -- use forward difference at lower boundary
    dydx[0] = (y[1] - y[0]) / dx
    # -- use backward difference at upper boundary
    dydx[-1] = (y[-1] - y[-2]) / dx

    return dydx

def tridiag(a, b, C, D):
    """
    tridiagonal matrix algorithm
    a=subdiag, b=diag, C=superdiag, D=rhs
    """
    n = len(a)
    V = np.zeros(n)
    G = np.zeros(n)
    U = np.zeros(n)
    x = np.zeros(n)

    V[0] = b[0].copy()
    G[0] = C[0] / V[0]
    U[0] = D[0] / V[0]

    for i in range(1, n):  # nr of nodes
        V[i] = b[i] - a[i] * G[i - 1]
        U[i] = (D[i] - a[i] * U[i - 1]) / V[i]
        G[i] = C[i] / V[i]

    x[-1] = U[-1]
    inn = n - 2
    for i in range(inn, -1, -1):
        x[i] = U[i] - G[i] * x[i + 1]
    return x

def smooth(a, WSZ):
    """
    smooth a by taking WSZ point moving average.
    NOTE: even WSZ is converted to next odd number.
    """
    WSZ = int(np.ceil(WSZ) // 2 * 2 + 1)
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid') / WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2] / r
    stop = (np.cumsum(a[:-WSZ:-1])[::2] / r)[::-1]
    x = np.concatenate((start, out0, stop))
    return x

def spatial_average(y, x=None, method='arithmetic'):
    """
    Calculates spatial average of quantity y, from node points to soil compartment edges
    Args: 
        y (array): quantity to average
        x (array): grid,<0, monotonically decreasing [m]
        method (str): flag for method 'arithmetic', 'geometric','dist_weighted'
    Returns: 
        f (array): averaged y, note len(f) = len(y) + 1
    """

    N = len(y)
    f = np.empty(N+1)  # Between all nodes and at surface and bottom
    if method is 'arithmetic':
        f[1:-1] = 0.5*(y[:-1] + y[1:])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'geometric':
        f[1:-1] = np.sqrt(y[:-1] * y[1:])
        f[0] = y[0]
        f[-1] = y[-1]

    elif method is 'dist_weighted':                                             # En ymmärrä, ei taida olla käyttössä
        a = (x[0:-2] - x[2:])*y[:-2]*y[1:-1]
        b = y[1:-1]*(x[:-2] - x[1:-1]) + y[:-2]*(x[1:-1] - x[2:])

        f[1:-1] = a / b
        f[0] = y[0]
        f[-1] = y[-1]

    return f