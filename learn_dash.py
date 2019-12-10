# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:14:26 2019

@author: slauniai
"""
#import pandas as pd
import numpy as np

import dash
import dash_core_components as dcc
import dash_html_components as html
from visualization.visual import read_results, timeseries
import plotly.io as pio
pio.renderers.default = "browser"

# xarray functions
# Dataset.filter_by_attrs(self, **kwargs # filter by arguments

# test apes results in netcdf
ff = r'c:\\Repositories\pyAPES_skeleton\results\pyAPES\testresults.nc'

x = read_results(ff)

# remove simulation dimension in case only one
x = x.squeeze()

za = x['canopy_z']
zs = x['soil_z']
lad = x['canopy_lad']
LAI = x['canopy_LAI']
ptypes = x['canopy_planttypes']

x = x.drop(['canopy_z', 'soil_z', 'canopy_lad', 'canopy_planttypes', 'canopy_LAI'])

variables = x.keys()
dimensions = x.dims
coords = x.coords


# get planttype_profiles into own dataset
v = []
vv = []
for k in variables:
    if 'planttype' in x[k].dims and 'canopy' in x[k].dims and 'date' in x[k].dims:
        v.append(k) 

pt_profiles = x[v]
vv.extend(v); del v

# pt timeseries
v = []
for k in variables:
    if 'planttype' in x[k].dims and 'date' in x[k].dims and 'canopy' not in x[k].dims:
        v.append(k) 

pt_ts = x[v]
vv.extend(v); del v

# get canopy_profiles into own dataset
v = []
for k in variables:
    if  'canopy' in x[k].dims and 'date' in x[k].dims and 'planttype' not in x[k].dims:
        v.append(k) 
profiles = x[v]
vv.extend(v); del v

# soil_profiles
v = []
for k in variables:
    if  'soil' in x[k].dims and 'date' in x[k].dims:
        v.append(k) 
soil_profiles = x[v]
vv.extend(v); del v

# only timeseries remaining when all others are removed
ts = x.drop(vv)

f = timeseries(ts, 'canopy_GPP')


#%% get dash app
app = dash.Dash(__name__)
app.layout = html.Div([dcc.Graph(id='plot-graph', figure=f)])

if __name__ == '__main__':
    app.run_server(debug=True)