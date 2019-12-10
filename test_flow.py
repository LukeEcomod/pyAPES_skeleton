# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 20:47:48 2019

@author: slauniai
"""
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

EPS = np.finfo(float)

from canopy.canopy_asl_flow import Flow, mixing_length
from tools.utilities import lad_weibul

#zmax = np.arange(20.0, 300, 50)
zmax = 20.0
dz = 0.2

z = np.arange(0, zmax + dz, dz)
hc = 15.0
lad = lad_weibul(z=z, LAI=3.0, h=hc, hb=0.4*hc, species='pine') 

p = {'zos': 1e-3, 'dPdx': 0.0, 'Cd': 0.15, 'Utop': 3.0, 'Ubot': 0.0, 
     'Sc':{'T': 1.0, 'H2O': 1.0, 'CO2': 1.0}}


# test 
model = Flow(z, lad, hc, p)
plt.figure()
res = []
res.append([model.U.copy(), model.ust.copy(), model.z.copy()])

for k in range(0, 10):
    zmax += 10.0
    z = np.arange(0,zmax+dz, dz) 

    model.flow_stats(z, lad, hc, 3.0, Ubot = 0.0)
    res.append([model.U.copy(), model.ust.copy(), model.z.copy()])

    plt.subplot(121)
    plt.plot(res[k][0], res[k][2], '-'); plt.xlabel('U'); plt.ylabel('z (m)')
    plt.subplot(122)
    plt.plot(res[k][1], res[k][2], '-', label='%.1f' %zmax); plt.xlabel('ust'); plt.ylabel('z (m)')

plt.legend()