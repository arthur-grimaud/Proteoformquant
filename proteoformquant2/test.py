import dash
from dash import html
from dash import dcc
from matplotlib.pyplot import figure
import plotly.graph_objects as go
from dash.dependencies import Input, Output
import pickle
from matplotlib.pyplot import cm
from Classes.envelope import Envelope
import numpy as np
import dash_bootstrap_components as dbc
import plotly.io as pio
import json
from Utils import misc
import plotly.express as px
import math
from dash.exceptions import PreventUpdate
from scipy import stats
from kneed import KneeLocator
from statistics import mean



minSubLen = 5
l = ["a","b","c","d","e","f","g","h","i","j"]
x = [-5,2,3,4,7,10,11,12,59,58]

subsetsL = [l]
subsetsX = [x]


print(subsetsL)
print(subsetsX)

l = [None]

print(l.any(None))



