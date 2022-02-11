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


l = [10,20,10,10]
print(len(l))
for i in range(0,len(l)-1):
    print(i)

with open('pfq_out_obj_test.pkl', 'rb') as inp:
    exp = pickle.load(inp) 


for spectra in exp.spectra.values():
     print(spectra.psms[0].proteoform.theoFrag)
     print(spectra.psms[0].annotation)


# uniqueFragments = []

# psmA =  exp.spectra[].psms[].theoFrag
# psmB =  exp.spectra[].psms[].theoFrag

# for fragType in psmA.proteoform.theoFrag.keys():
#     for fragment in psmA.proteoform.theoFrag[fragType]:
#         if psmA.proteoform.theoFrag[fragType][fragment] =! psmB.proteoform.theoFrag[fragType][fragment]
#             uniqueFragments.append(fragment)

