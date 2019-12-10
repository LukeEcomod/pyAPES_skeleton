# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 16:51:35 2019

@author: slauniai
"""
import xarray as xr
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import pandas as pd

from datetime import datetime
import numpy as np

def to_datetime(date):
    """
    Converts a numpy datetime64 object to a python datetime object 
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    """
    
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    print(timestamp)
    return datetime.utcfromtimestamp(timestamp)

def read_results(outputfiles):
    """
    Opens simulation results netcdf4 dataset in xarray
    (or multiple in list of xarrays)
    Args:
        outputfiles (str or list of str)
    Returns:
        results (xarray or list of xarrays):
    """

    if type(outputfiles) != list:
        outputfiles = [outputfiles]

    results = []
    for outputfile in outputfiles:
        fp = outputfile
        result = xr.open_dataset(fp)
        result.coords['simulation'] = result.simulation.values
        result.coords['soil'] = result.soil_z.values
        result.coords['canopy'] = result.canopy_z.values
        result.coords['planttype'] = result.canopy_planttypes.values
#        result.coords['planttype'] = ['pine','spruce','decid','shrubs']
        results.append(result)

    if len(results) == 1:
        return results[0]
    else:
        return results

# plotly figure defs
        
def timeseries(x, k, text=None):

    df = pd.Series(index=x.date.values, data=x[k].values)
    t = df.index
    txt = x[k].attrs['units']

    try:
        titletext = txt.split('[')[0]
    except:
        titletext = ' '
    try:
        units = txt.split('[')[-1]
    except:
        units = ' '
    
    # create fig
    fig = go.Figure({"data": [{"type": "scatter", "x":t, "y": x[k].values}]})
    fig.update_layout( width=1200,
                      title=titletext,
                      yaxis_title=units
                      )
    return fig

def timeseries_subplots(x, k, n):
    if type(k) == str:
        k =  [k]
    df = pd.Series(index=x.date.values, data=x[k].values)
    t = df.index
    
    subtitles = []
    [subtitles.append(x[j].attrs['units']) for j in k]
    
    fig = make_subplots(rows=n+1, cols=2, subplot_titles=tuple(subtitles))

    for j in range(0, n):
        fig.add_trace(go.Scatter(x=t, y=x[k[j]].values),
                      row=j+1, col=1
                      )
    fig.update_layout(height=n*200,
                      width=1200,
                      showlegend=False
                      #title=titletext,
                      #yaxis_title=units
                      )
    return fig

def timeseries_with_gradient(x, k, z, step):
    pairs = {'canopy_NEE': 'canopy_co2_flux',
             'canopy_LE': 'canopy_latent_heat_flux',
             'canopy_SH': 'canopy_sensible_heat_flux',
             'forcing_precipitation': 'canopy_throughfall_ml',
             'forcing_par': 'canopy_par_down'}
    
    if type(k) == str:
        k =  [k]
    n = len(k)
    
    df = pd.Series(index=x.date.values, data=x[k].values)
    t = df.index
    
    subtitles = []
    for j in k:
        subtitles.append(x[j].attrs['units'])
        if j in pairs and pairs[j] in x:
            p = pairs[j]
            subtitles.append(x[p].attrs['units'])
        else:
            subtitles.append('')
    
    fig = make_subplots(rows=n+1, cols=2, subplot_titles=tuple(subtitles))

    for j in range(0, n):
        fig.add_trace(go.Scatter(x=t, y=x[k[j]].values),
                      row=j+1, col=1
                      )
        if k[j] in pairs:
            p = pairs[k[j]]
            if p in x:
                y = x[p][step,:].values
            fig.add_trace(go.Scatter(x=y, y=z),
                          row=j+1, col=2
                          )
    fig.update_layout(height=n*200,
                      width=1200,
                      showlegend=False
                      #title=titletext,
                      #yaxis_title=units
                      )
    return fig