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

from Utils import misc
# ---------------------------------------------------------------------------- #
#                                  Data import                                 #
# ---------------------------------------------------------------------------- #

with open('pfq_out_obj.pkl', 'rb') as inp:
    exp = pickle.load(inp) 

# ---------------------------------------------------------------------------- #
#                                    Layout                                    #
# ---------------------------------------------------------------------------- #

app = dash.Dash()
template = "plotly_white"

app.layout = html.Div([

html.Button('reload', id='reload_intens'),
html.Div(children=[
    html.Label('Graph'),
    dcc.Graph(id = 'line_plot_global_intens'),
    ]
),

html.Br(),

html.Button('reload', id='reload_envelope'),

html.Br(),

html.Div(children=[
    html.Label('Elution profile (Select proteoform 1)'),
    dcc.Dropdown( id="dropdown_envelope_1"),
    dcc.Graph(id = 'envelope_plot_1')],
    style={'width': '49%', 'display': 'inline-block'}
    ),

html.Div(children=[
    html.Label('Elution profile (Select proteoform 2)'),
    dcc.Dropdown( id="dropdown_envelope_2"),
    dcc.Graph(id = 'envelope_plot_2')],
    style={'width': '49%', 'display': 'inline-block'}
    ),

html.Br(),

html.Button('reload', id='reload_all_enveloppes'),
html.Div(children=[
    html.Label('All determined elution profil envelope'),
    dcc.Graph(id = 'plot_all_enveloppes'),
    ]
),

html.Br(),

html.Button('reload', id='reload_all_enveloppes_3d'),
html.Div(children=[
    html.Label('All determined elution profil envelope'),
    dcc.Graph(id = 'plot_all_enveloppes_3d'),
    ]
)

])


# ---------------------------------------------------------------------------- #
#                                   Callbacks                                  #
# ---------------------------------------------------------------------------- #



# --------------------------- spectra chromatogram --------------------------- #

@app.callback(
    Output('line_plot_global_intens', 'figure'),
    Input('reload_intens', "id")
)
def plotAnnotMsmsAndPrecIntens(input):
    #print(exp)
    #print("calling plotAnnotMsmsAndPrecIntens")
    precIntens = [spectrum.getPrecIntens() for spectrum in exp.spectra.values()]
    annotIntens  = [spectrum.getSumIntensAnnotFrag() for spectrum in exp.spectra.values()]
    rt = [spectrum.getRt() for spectrum in exp.spectra.values()]

    fig = go.Figure()
    fig.add_scatter( x=rt, y=precIntens, mode='markers', marker=dict(size=4, color="red"), name='Precursor Intensity' )
    fig.add_scatter( x=rt, y=annotIntens, mode='markers', marker=dict(size=4, color="blue"), name='Annotated Fragment Summed' )
    fig.update_layout(template=template)
    return fig  





# --------------------------------- envelope --------------------------------- #


@app.callback(
    Output('dropdown_envelope_1', 'options'),
    Input('reload_envelope', "value")
)
def dropdownOptionsEnvelopes_1(input):
    return dropdownOptionsEnvelopes(input)

@app.callback(
    Output('dropdown_envelope_2', 'options'),
    Input('reload_envelope', "value")
)
def dropdownOptionsEnvelopes_2(input):
    return dropdownOptionsEnvelopes(input)

def dropdownOptionsEnvelopes(input):
    """Create option list for dropdown item"""
    options = [{"label": proteo[1].getModificationBrno() + "   " + str(len(proteo[1].getValidatedLinkedPsm()) ) ,"value": proteo[0]} for proteo in exp.proteoforms.items() if len(proteo[1].getValidatedLinkedPsm()) > 0 ]
    #print(options)
    return options

# ---------------------------------------------------------------------------- #

@app.callback(
    Output('envelope_plot_1', 'figure'),
    Input('dropdown_envelope_1', 'value')
)
def plotEnvelope1(proteo):
    return plotEnvelope(proteo)

@app.callback(
    Output('envelope_plot_2', 'figure'),
    Input('dropdown_envelope_2', 'value')
)
def plotEnvelope2(proteo):
    return plotEnvelope(proteo)

def plotEnvelope(proteo):

    xData = [psm.spectrum.getRt() for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()]
    yData = [psm.getPrecIntensRatio() for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()]

    fig = go.Figure()
    fig.add_scatter( x=xData, y=yData, mode='markers', marker=dict(size=10, color="black"), name='Precursor Intensity' )

    for env in exp.proteoforms[proteo].envelopes: #if envelope has been computed add line to the plot

        xDataEnv = list(range(int(min(xData)),int(max(xData)),1))

        yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
        if yDataEnvEstim[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvEstim, mode='lines', marker=dict(size=4, color="orange"), name='Estimated Parameters', line_shape='spline' )

        yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
        if yDataEnvFitted[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4, color="red"), name='Fitted Parameters', line_shape='spline' )

    titleText = "Proteoform: {0} <br>Parameters Estimated: {1} <br>KS: {3} <br>Parameters Fitted: {2} <br>KS: {4}".format(exp.proteoforms[proteo].getModificationBrno(), parametersEstim , parametersFitted, env.KsEstimated, env.KsFitted)

    fig.update_layout(title=go.layout.Title(text=titleText, font=dict(
            family="Courier New, monospace",
            size=10,
        )))
    fig.update_layout(template=template)
    return fig  

# ---------------------------------------------------------------------------- #
#                                 all envelopes                                #
# ---------------------------------------------------------------------------- #

@app.callback(
    Output('plot_all_enveloppes', 'figure'),
    Input('reload_all_enveloppes', "id")
)
def plotAllEnvelopes(input):


    rt_range = [spectrum.getRt() for spectrum in exp.spectra.values()]
    fig = go.Figure()


    colors=misc.linear_gradient("#4682B4","#FFB347",len([x for x in exp.proteoforms.values() if len(x.envelopes)>0]))


    i = 0
    for proteoform in [x for x in exp.proteoforms.values() if len(x.envelopes)>0]:
        c = colors["hex"][i]
        i += 1
        if len(proteoform.envelopes) > 0:
            for env in proteoform.envelopes: #if envelope has been computed add line to the plot

                xDataEnv = list(range(int(min(rt_range)),int(max(rt_range)),1))
                zDataEnv = [proteoform.getMzFirstPsm() for x in xDataEnv]

                yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
                if yDataEnvFitted[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4, color=c), name=proteoform.getModificationBrno(), line_shape='spline' )
                else:
                    yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
                    if yDataEnvEstim[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvEstim, mode='lines', marker=dict(size=4, color=c), name=proteoform.getModificationBrno(), line_shape='spline' )
    
    
    fig.update_layout(template=template)              
    return fig  


# ---------------------------------------------------------------------------- #
#                               all envelopes 3D                               #
# ---------------------------------------------------------------------------- #

@app.callback(
    Output('plot_all_enveloppes_3d', 'figure'),
    Input('reload_all_enveloppes_3d', "id")
)
def plotAllEnvelopes3d(input):


    rt_range = [spectrum.getRt() for spectrum in exp.spectra.values()]
    fig = go.Figure()



    colors=misc.linear_gradient("#4682B4","#FFB347",len([x for x in exp.proteoforms.values() if len(x.envelopes)>0]))


    i = 0
    for proteoform in [x for x in exp.proteoforms.values() if len(x.envelopes)>0]:
        if len(proteoform.envelopes) > 0:
            for env in proteoform.envelopes: #if envelope has been computed add line to the plot

                xDataEnv = list(range(int(min(rt_range)),int(max(rt_range)),1))
                zDataEnv = [proteoform.getMzFirstPsm() for x in xDataEnv]

                yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
                if yDataEnvFitted[0] != None: fig.add_trace(go.Scatter3d( x=xDataEnv, y=yDataEnvFitted, z=zDataEnv, mode='lines', name=proteoform.getModificationBrno() ))
                else:
                    yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
                    if yDataEnvEstim[0] != None: fig.add_trace(go.Scatter3d( x=xDataEnv, y=yDataEnvEstim, z=zDataEnv, mode='lines', name=proteoform.getModificationBrno()))
    
    fig.update_layout(template=template,height=1000)
               
    return fig  


# ---------------------------------------------------------------------------- #
#                                     start                                    #
# ---------------------------------------------------------------------------- #

app.run_server(debug=True)


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #

