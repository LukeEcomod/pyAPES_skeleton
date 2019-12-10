# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 17:51:30 2019

@author: slauniai
"""
import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import datetime

app = dash.Dash(__name__)

x = df.index
y = df[column]
print(x)
print(y)
print(isinstance(df.index[0], datetime.datetime))

trace = go.Scatter(x=x, y=y)
data = [trace]
layout = dict(title=column)
fig = dict(data=data, layout=layout)

app.layout = html.Div([
    dcc.Graph(id='plot-graph', figure=fig)
    ])


if __name__ == '__main__':
    app.run_server(debug=True)