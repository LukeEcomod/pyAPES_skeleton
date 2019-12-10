# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 20:22:22 2019

@author: slauniai
"""
import pandas as pd
from visualization.visual import read_results, timeseries, timeseries_subplots, timeseries_with_gradient

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
#import plotly.io as pio
#pio.renderers.default = "browser"

#%% get data
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
        print(x[k].attrs)
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
cprofs = x[v]
vv.extend(v); del v

# soil_profiles
v = []
for k in variables:
    if  'soil' in x[k].dims and 'date' in x[k].dims:
        v.append(k) 
sprofs = x[v]
vv.extend(v); del v

# only timeseries remaining when all others are removed
ts = x.drop(vv)
#%%
fig = timeseries(x, 'canopy_GPP')

def get_variables(data):
    keys = ([{'label': k, 'value': k} for k in data.keys()])
    return keys

#%% 
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


app.layout = html.Div([

    # Page Header
    html.Div([
        html.H1('pyAPES results explorer')
    ]),

    # Dropdown Grid
    html.Div([
        html.Div([
            #Dropdown
            html.Div([
                html.Div('Select timeseries', className='three columns'),
                html.Div(dcc.Dropdown(id='variable-selector',
                                      options=get_variables(x),
                                      multi=True),
                         className='nine columns')
            ]),

            # Dropdown
            html.Div([
                html.Div('Select profile', className='three columns'),
                html.Div(dcc.Dropdown(id='profile-selector',
                                      options=get_variables(cprofs)),
                         className='nine columns')
            ]),

        ], className='six columns'),

#        # Empty
#        html.Div(className='six columns'),
#    ], className='twelve columns'),
#
#    # Match Results Grid
#    html.Div([
#
#        # Match Results Table
#        html.Div(
#            html.Table(id='match-results'),
#            className='six columns'
#        ),

        # Season Summary Table and Graph
        html.Div([
            # summary table
            dcc.Graph(id='timeseries-graph'),

            # graph
            #dcc.Graph(id='profile-graph')
            # style={},

        ], className='twelve columns')
    ]),
])

@app.callback(
    Output('timeseries-graph', 'figure'),
    [Input('variable-selector', 'value')])
def update_timeseries(k):
    n = len(k)
    #print(type(k), n)
    if isinstance(k, str) or n == 1:
        return timeseries(x, k[0])
    else:
        return timeseries_with_gradient(x, k, za, step=24)

#app.layout = html.Div([
#        
#    html.H1('pyAPES - results explorer'),
#    
#    html.Label('Select timeseries'),
#    dcc.Dropdown(
#        options=[
#            {'label': u'canopy_GPP', 'value': 'canopy_GPP'},
#            {'label': u'canopy_LE', 'value': 'canopy_LE'},
#            {'label': u'canopy_SH', 'value': 'canopy_SH'}
#        ],
#        value='canopy_GPP'
#    ),
#
#    html.Label('Multi-Select Dropdown'),
#    dcc.Dropdown(
#        options=[
#            {'label': 'New York City', 'value': 'NYC'},
#            {'label': u'Montréal', 'value': 'MTL'},
#            {'label': 'San Francisco', 'value': 'SF'}
#        ],
#        value=['MTL', 'SF'],
#        multi=True
#    ),
#
#    dcc.Graph(
#        id='timeseries-graph',
#        figure=fig
#    ),
##    html.Label('Radio Items'),
##    dcc.RadioItems(
##        options=[
##            {'label': 'New York City', 'value': 'NYC'},
##            {'label': u'Montréal', 'value': 'MTL'},
##            {'label': 'San Francisco', 'value': 'SF'}
##        ],
##        value='MTL'
##    ),
#
##    html.Label('Checkboxes'),
##    dcc.Checklist(
##        options=[
##            {'label': 'New York City', 'value': 'NYC'},
##            {'label': u'Montréal', 'value': 'MTL'},
##            {'label': 'San Francisco', 'value': 'SF'}
##        ],
##        value=['MTL', 'SF']
##    ),
#
#    html.Label('Text Input'),
#    dcc.Input(value='MTL', type='text'),
#
##    html.Label('Slider'),
##    dcc.Slider(
##        min=0,
##        max=9,
##        marks={i: 'Label {}'.format(i) if i == 1 else str(i) for i in range(1, 6)},
##        value=5,
##    ),
#], style={'columnCount': 3})

if __name__ == '__main__':
    app.run_server(debug=True)