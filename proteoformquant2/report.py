import dash
from dash import html
from dash import dcc
from matplotlib.pyplot import figure
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State
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
import pandas as pd


# ---------------------------------------------------------------------------- #
#                                  Data import                                 #
# ---------------------------------------------------------------------------- #

with open('pfq_out_obj_test_2b.pkl', 'rb') as inp:
    exp = pickle.load(inp) 

print(exp.getDatasetMetrics())


# ----------------------------------- test ----------------------------------- #


# ----------------------------------- test ----------------------------------- #


# ---------------------------------------------------------------------------- #
#                                    Layout                                    #
# ---------------------------------------------------------------------------- #

app = dash.Dash()
template = "plotly_white"
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
app.layout = html.Div([

    html.Div(children=[
        # dcc.Input(id="input_min_mz", type= "number", placeholder="Minimum Precursor Mz"),
        # dcc.Input(id="input_max_mz", type= "number", placeholder="Maximum Precursor Mz"),
        dcc.RangeSlider( id='range_mz', min=exp.getMzRange()[0], max=exp.getMzRange()[1], step=2, value=exp.getMzRange(), tooltip={"placement": "bottom", "always_visible": True}),
        html.Label('Chromatogram of precursor and annotated fragments intensities'),
        dcc.Graph(id = 'line_plot_global_intens'),
        dcc.Graph(id = 'line_plot_itens_comp'),
        ]
    ),

    html.Br(),

    html.Div(children=[
        html.Label('All determined elution profil envelopes'),
        dcc.Graph(id = 'plot_all_enveloppes'),
        ]
    ),

    html.Br(),

    html.Div(children=[
        html.Label('Elution profile envelope (Select proteoform 1)'),
        dcc.Dropdown( id="dropdown_envelope_1"),
        dcc.Graph(id = 'envelope_plot_1')],
        style={'width': '49%', 'display': 'inline-block'}
        ),

    html.Div(children=[
        html.Label('Elution profile envelope (Select proteoform 2)'),
        dcc.Dropdown( id="dropdown_envelope_2"),
        dcc.Graph(id = 'envelope_plot_2')],
        style={'width': '49%', 'display': 'inline-block'}
        ),

    html.Br(),


    html.Div(children=[
        dcc.RangeSlider( id='range_mz_2', min=exp.getMzRange()[0], max=exp.getMzRange()[1], step=2, value=exp.getMzRange(), tooltip={"placement": "bottom", "always_visible": True}),
        html.Label('All determined elution profil envelope with m/z dimension'),
        dcc.Graph(id = 'plot_all_enveloppes_3d'),
        ]
    ),

    html.Br(),

    html.Div(children=[
        html.Label('Proteoform ratios'),
        dcc.Graph(id = 'plot_quantification'),
        ]
    ),

    html.Br(),

    html.Div(children=[
        html.Label('Additional Plot'),
        dcc.Graph(id = 'misc_plot'),
        ]
    ),

    # ---------------------------------- pop up ---------------------------------- #

   
    dbc.Modal(
        [
            dbc.ModalHeader("Header",id = 'modalheader'),
            dbc.ModalBody("This is the content of the modal"),
            dbc.ModalFooter(
                dbc.Button("Close", id="close", className="ml-auto")
            ),
        ],
        size="xl",
        id="modal",
    )

], id="mainDiv")



# ---------------------------------------------------------------------------- #
#                                   Callbacks                                  #
# ---------------------------------------------------------------------------- #

# ----------------------------------- test ----------------------------------- #

@app.callback([Output("modal", "children"),
               Output("modal", "is_open")],
              [Input("line_plot_global_intens", "clickData"),
               Input("line_plot_itens_comp", "clickData"),
               Input("plot_all_enveloppes", "clickData"),
               Input("envelope_plot_1", "clickData"),
               Input("envelope_plot_2", "clickData"),
               Input("close", "n_clicks")],
              [State("modal", "is_open"),
               State("modal", "children")])
def spectrum_popup(v1, v2, v3, v4,v5, clicked,is_open,children):
    ctx = dash.callback_context

    if ctx.triggered[0]['prop_id'] == 'close.n_clicks':
        # you pressed the closed button, keeping the modal children as is, and 
        # close the modal itself. 
        return children, False 

    elif len(ctx.triggered) == 0:
        raise dash.exceptions.PreventUpdate

    elif "points" in ctx.triggered[0]["value"] :
        # you clicked in the graph, returning the modal children and opening it
        try:
            key = ctx.triggered[0]["value"]["points"][0]["customdata"]
            #get spectrum info string
            spectrum = exp.spectra[key]
            strInfo = []
            for psm in spectrum.psms:
                strInfo.append("Rank: {0}, Proteoform: {1}, isValidated: {2}, ratio: {3}".format(psm.rank, psm.proteoform.getModificationBrno(), psm.isValidated, psm.ratio))
                strInfo.append(html.Br())

            return [dbc.ModalHeader("Spectrum"),
                    dbc.ModalBody(
                    html.Div([
                        html.P(strInfo,style = {'fontFamily':'monospace'}),
                        dcc.Graph(figure = spectrum_plot(spectrum),id = 'plot_spectrum')      
                    ])),
                    dbc.ModalFooter(dbc.Button("Close", id="close"))
                    ], True
        except KeyError: #In case click data is not acorresponding to a spectrum
            print("Cannot display additional information for that object")

    else:
        raise dash.exceptions.PreventUpdate

def spectrum_plot(spectrum):
    m = 0
    #plot PSM Annotation

    fig = go.Figure()
    for psm in spectrum.getValidatedPsm():
        print(m)
        # get all frag type mz and intens in singles lists
        fragMz = [j for i in [fragType["mzTheo"] for fragType in psm.annotation.values()] for j in i]
        fragIntens =  [j for i in [fragType["intens"] for fragType in psm.annotation.values()] for j in i]
        fragIntens = [i+(m*10000) for i in fragIntens]
        fig.add_scatter( x=fragMz , y=fragIntens, mode='markers',marker_symbol=m, marker=dict( size=7), name = psm.getModificationBrno() )
        m += 1

    for f in range(len(spectrum.fragIntens)):

        fig.add_scatter( x=[spectrum.fragMz[f],spectrum.fragMz[f]] , y=[0, spectrum.fragIntens[f]], mode='lines', marker=dict( color = "black",size=2) )

    fig.update_layout(template=template)
    return fig

# ----------------------------------- test ----------------------------------- #


# --------------------------- spectra chromatogram --------------------------- #

@app.callback(
    Output('line_plot_global_intens', 'figure'),
    Input('range_mz', "value")
    # Input('input_max_mz', "value")
)
def plotAnnotMsmsAndPrecIntens(minMaxMz):
    #print(exp)
    #print("calling plotAnnotMsmsAndPrecIntens")
    minMz = minMaxMz[0]
    maxMz = minMaxMz[1]
    print(minMz,maxMz)
    precIntens = [spectrum.getPrecIntens() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]
    annotIntens  = [spectrum.getSumIntensAnnotFrag() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]
    spectrumKey  = [spectrum for spectrum in exp.spectra.keys() if exp.spectra[spectrum].getPrecMz() > minMz and exp.spectra[spectrum].getPrecMz() < maxMz]
    rt = [spectrum.getRt() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]

    fig = go.Figure()
    fig.add_scatter( x=rt, y=precIntens, mode='markers', marker=dict(size=4, color="red"), name='Precursor Intensity',customdata=spectrumKey )
    fig.add_scatter( x=rt, y=annotIntens, mode='markers', marker=dict(size=4, color="blue"), name='Annotated Fragment Summed',customdata=spectrumKey )
    fig.update_layout(template=template)
    return fig  


@app.callback(
    Output('line_plot_itens_comp', 'figure'),
    Input('range_mz', "value")
)
def plotAnnotMsmsVsPrecIntens(minMaxMz):
    minMz = minMaxMz[0]
    maxMz = minMaxMz[1]
    precIntens = [spectrum.getPrecIntens() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]
    #annotIntens  = [spectrum.getSumIntensAnnotFrag() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]
    sumFragIntens  = [spectrum.getSumIntensFrag() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]
    spectrumKey  = [spectrum for spectrum in exp.spectra.keys() if exp.spectra[spectrum].getPrecMz() > minMz and exp.spectra[spectrum].getPrecMz() < maxMz]
    fig = go.Figure()
    #fig.add_scatter( x=np.log(annotIntens), y=np.log(precIntens), mode='markers', marker=dict(size=4, color="red"), name='Annotated Fragment Summed',customdata=spectrumKey )
    fig.add_scatter( x=np.log(sumFragIntens), y=np.log(precIntens), mode='markers', marker=dict(size=4, color="red"), name='Fragment Sum',customdata=spectrumKey )
    fig.update_layout(template=template)
    fig.update_xaxes(title_text='log(Fragments Summed Intensity)')
    fig.update_yaxes(title_text='log(Precursor Intensity)')
    return fig


# --------------------------------- envelope --------------------------------- #


@app.callback(
    Output('dropdown_envelope_1', 'options'),
    Input('mainDiv', "id")
)
def dropdownOptionsEnvelopes_1(input):
    return dropdownOptionsEnvelopes(input)

@app.callback(
    Output('dropdown_envelope_2', 'options'),
    Input('mainDiv', "id")
)
def dropdownOptionsEnvelopes_2(input):
    return dropdownOptionsEnvelopes(input)

def dropdownOptionsEnvelopes(input):
    """Create option list for dropdown item"""
    options = [{"label": proteo[1].getModificationBrno() + "  " + str(len(proteo[1].getValidatedLinkedPsm()) ) ,
                "value": proteo[0], 
                "sort":len(proteo[1].getValidatedLinkedPsm()) } for proteo in exp.proteoforms.items() if len(proteo[1].getLinkedPsm()) > 4]
    options = sorted(options, key=lambda x: x["sort"], reverse=True)

    options = [{k: v for k, v in d.items() if k != 'sort'} for d in options]


    
    #print(options)
    return options

# ---------------------------------------------------------------------------- #

@app.callback(
    Output('envelope_plot_1', 'figure'),
    Input('dropdown_envelope_1', 'value')
)
def plotEnvelope1(proteo):
    if proteo is None:
        raise PreventUpdate
    else:
        return plotEnvelope(proteo)

@app.callback(
    Output('envelope_plot_2', 'figure'),
    [Input('dropdown_envelope_1', 'value'),
    Input('dropdown_envelope_2', 'value'),]
)
def plotEnvelope2(proteo1, proteo2):
    if proteo1 is None or proteo2 is None:
        raise PreventUpdate
    else:
        return plotEnvelopeTwo(proteo1, proteo2)


def plotEnvelope(proteo):

    fig = go.Figure()

    xData = [psm.spectrum.getRt() for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()]
    yDataSpectrum = [psm.spectrum.getPrecIntens() for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()]
    yDataPsm = [psm.getPrecIntensRatio() for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()]
    spectrumKey  = [psm.spectrum.id for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()]

    fig.add_scatter( x=xData, y=yDataSpectrum, mode='markers', marker=dict(size=10, color="grey"), name='Spectrum Intensity', customdata=spectrumKey )
    fig.add_scatter( x=xData, y=yDataPsm, mode='markers', marker=dict(size=7, color="red"), name='PSM Intensity',customdata=spectrumKey )

    #add lines between PSM and spectrum intens points
    for i in range(0,len(xData),1):
        fig.add_scatter( x=[xData[i],xData[i]] , y=[yDataSpectrum[i],yDataPsm[i]], mode='lines', marker=dict(size=2, color="grey") )

    env = exp.proteoforms[proteo].getEnvelope()
    if env != None: #if envelope has been computed add line to the plot
        
        xDataEnv = list(range(int(min(xData)),int(max(xData)),1))

        yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
        if yDataEnvEstim[0] != None: 
            fig.add_scatter( x=xDataEnv, y=yDataEnvEstim, mode='lines', marker=dict(size=4, color="orange"), name='Estimated Parameters', line_shape='spline' )

        yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
        if yDataEnvFitted[0] != None: 
            fig.add_scatter( x=xDataEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4, color="red"), name='Fitted Parameters', line_shape='spline' )

            titleText = "Proteoform: {0} <br>Parameters Estimated: {1} <br>Score: {3} <br>Parameters Fitted: {2} <br>Score: {4} ".format(exp.proteoforms[proteo].getModificationBrno(), parametersEstim , parametersFitted, env.scoreEstimated, env.scoreFitted)

        fig.update_layout(title=go.layout.Title(text=titleText, font=dict(
                family="Courier New, monospace",
                size=10,
            )))
    fig.update_layout(template=template)
    return fig  


def plotEnvelopeTwo(proteo1, proteo2):

    fig = go.Figure()

    xData1 = [psm.spectrum.getRt() for psm in exp.proteoforms[proteo1].getValidatedLinkedPsm()]
    yDataSpectrum1 = [psm.spectrum.getPrecIntens() for psm in exp.proteoforms[proteo1].getValidatedLinkedPsm()]
    yDataPsm1 = [psm.getPrecIntensRatio() for psm in exp.proteoforms[proteo1].getValidatedLinkedPsm()]
    spectrumKey1  = [psm.spectrum.id for psm in exp.proteoforms[proteo1].getValidatedLinkedPsm()]

    xData2 = [psm.spectrum.getRt() for psm in exp.proteoforms[proteo2].getValidatedLinkedPsm()]
    yDataSpectrum2 = [psm.spectrum.getPrecIntens() for psm in exp.proteoforms[proteo2].getValidatedLinkedPsm()]
    yDataPsm2 = [psm.getPrecIntensRatio() for psm in exp.proteoforms[proteo2].getValidatedLinkedPsm()]
    spectrumKey2  = [psm.spectrum.id for psm in exp.proteoforms[proteo2].getValidatedLinkedPsm()]

    fig.add_scatter( x=xData1, y=yDataSpectrum1, mode='markers', marker=dict(size=10, color="grey"), name='Spectrum Intensity', customdata= spectrumKey1)

    fig.add_scatter( x=xData2, y=yDataSpectrum2, mode='markers', marker=dict(size=10, color="grey"), name='Spectrum Intensity',customdata=spectrumKey2 )

    fig.add_scatter( x=xData1, y=yDataPsm1, mode='markers', marker=dict(size=7, color="red"), name='PSM Intensity',customdata=spectrumKey1 )

    
    fig.add_scatter( x=xData2, y=yDataPsm2, mode='markers', marker=dict(size=7, color="blue"), name='PSM Intensity',customdata=spectrumKey2 )

    #add lines between PSM and spectrum intens points
    for i in range(0,len(xData1),1):
        fig.add_scatter( x=[xData1[i],xData1[i]] , y=[yDataSpectrum1[i],yDataPsm1[i]], mode='lines', marker=dict(size=2, color="grey") )

    for i in range(0,len(xData2),1):
        fig.add_scatter( x=[xData2[i],xData2[i]] , y=[yDataSpectrum2[i],yDataPsm2[i]], mode='lines', marker=dict(size=2, color="grey") )

    env1 = exp.proteoforms[proteo1].getEnvelope()
    if env1 != None: #if envelope has been computed add line to the plot
        xDataEnv1 = list(range(int(min(xData1)),int(max(xData1)),1))
        yDataEnvFitted1, parametersFitted1 = list(env1.getEnvelopeSerie(xDataEnv1, method = "fitted"))
        if yDataEnvFitted1[0] != None: 
            fig.add_scatter( x=xDataEnv1, y=yDataEnvFitted1, mode='lines', marker=dict(size=4, color="red"), name='Fitted Parameters', line_shape='spline' )

    env2 = exp.proteoforms[proteo2].getEnvelope()
    if env2 != None: #if envelope has been computed add line to the plot
        xDataEnv2 = list(range(int(min(xData2)),int(max(xData2)),1))
        yDataEnvFitted2, parametersFitted2 = list(env2.getEnvelopeSerie(xDataEnv2, method = "fitted"))
        if yDataEnvFitted2[0] != None: 
            fig.add_scatter( x=xDataEnv2, y=yDataEnvFitted2, mode='lines', marker=dict(size=4, color="blue"), name='Fitted Parameters', line_shape='spline' )
        
            
        fig.update_layout(title=go.layout.Title( font=dict(
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
    Input('mainDiv', "id")
)
def plotAllEnvelopes(input):


    rt_range = [spectrum.getRt() for spectrum in exp.spectra.values()]
    fig = go.Figure()


    colors=misc.linear_gradient("#4682B4","#FFB347",len([x for x in exp.proteoforms.values() if x.getEnvelope() != None]))


    i = 0
    for proteoform in [x for x in exp.proteoforms.values() if x.getEnvelope() != None]:
        c = colors["hex"][i]
        i += 1
        if proteoform.getEnvelope() != None:
            env = proteoform.getEnvelope() #if envelope has been computed add line to the plot

            xDataEnv = list(range(int(min(rt_range)),int(max(rt_range)),1))
            zDataEnv = [proteoform.getMzFirstPsm() for x in xDataEnv]

            yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
            if yDataEnvFitted[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4, color=c), name=proteoform.getModificationBrno(), line_shape='spline' )
            else:
                yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
                if yDataEnvEstim[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvEstim, mode='lines', marker=dict(size=4, color=c), name=proteoform.getModificationBrno(), line_shape='spline' )

    
    precIntens = [spectrum.getPrecIntens() for spectrum in exp.spectra.values()]
    rt = [spectrum.getRt() for spectrum in exp.spectra.values()]
    spectrumKey  = [spectrum for spectrum in exp.spectra.keys()]

    fig.add_scatter( x=rt, y=precIntens, mode='markers', marker=dict(size=4, color="black"), name='Precursor Intensity',customdata=spectrumKey)

 
    fig.update_layout(template=template)              
    return fig  

# ---------------------------------------------------------------------------- #
#                               all envelopes 3D                               #
# ---------------------------------------------------------------------------- #


@app.callback(
    Output('plot_all_enveloppes_3d', 'figure'),
    Input('range_mz_2', "value")
)
def plotAllEnvelopes3d(minMaxMz):
    fig = go.Figure()

    minMz = minMaxMz[0]
    maxMz = minMaxMz[1]


    #Display unassigned spectra:
    specFilt = [ spectrum for spectrum in exp.proteoform0.linkedSpectra if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz ]

    specFiltMz = [spectrum.getPrecMz() for spectrum in specFilt]
    specFiltRt = [spectrum.getRt() for spectrum in specFilt]
    specFiltIntens = [spectrum.getPrecIntens() for spectrum in specFilt]
    specFiltKey = [spectrum.getId() for spectrum in specFilt]

    fig.add_trace(go.Scatter3d( x=specFiltRt, y=specFiltIntens, z=specFiltMz, mode='markers',  marker=dict(size=2, color="grey",opacity=0.6), name="proteoform0",customdata=specFiltKey))
              



    proteoFilt = { proteoName:proteo for (proteoName,proteo) in exp.proteoforms.items() if proteo.getMzFirstPsm() > minMz and proteo.getMzFirstPsm() < maxMz}

    rt_range = exp.getRtRange()
    




    colors=misc.linear_gradient("#4682B4","#FFB347",len([x for x in proteoFilt.values() if x.getEnvelope() != None]))


    i = 0
    for proteoform in [x for x in proteoFilt.values() if x.getEnvelope() != None]:

        env = proteoform.getEnvelope() #if envelope has been computed add line to the plot

        xDataEnv = list(range(int(min(rt_range)),int(max(rt_range)),1))
        zDataEnv = [proteoform.getMzFirstPsm() for x in xDataEnv]

        yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
        if yDataEnvFitted[0] != None: 
            fig.add_trace(go.Scatter3d( x=xDataEnv, y=yDataEnvFitted, z=zDataEnv, mode='lines',  marker=dict(color=proteoform.getColor()), name=proteoform.getModificationBrno()))
        else:
            yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
            if yDataEnvEstim[0] != None: fig.add_trace(go.Scatter3d( x=xDataEnv, y=yDataEnvEstim, z=zDataEnv, mode='lines',  marker=dict(color=proteoform.getColor()), name=proteoform.getModificationBrno()))
    
    for proteoform in [x for x in proteoFilt.values()]:
        precIntens = [psm.spectrum.getPrecIntens() for psm in proteoform.getValidatedLinkedPsm()]
        rt = [psm.spectrum.getRt() for psm in proteoform.getValidatedLinkedPsm()]
        mz  = [psm.spectrum.getPrecMz() for psm in proteoform.getValidatedLinkedPsm()]
        spectrumKey  = [psm.spectrum.getId() for psm in proteoform.getValidatedLinkedPsm()]

        fig.add_trace(go.Scatter3d( x=rt, y=precIntens, z=mz, mode='markers', marker=dict(size=2, color=proteoform.getColor()),name=proteoform.getModificationBrno() ,customdata=spectrumKey))

    


    fig.update_layout(template=template,height=1000)
               
    return fig  


# ---------------------------------------------------------------------------- #
#                               proteoform ratios                              #
# ---------------------------------------------------------------------------- #



@app.callback(
    Output("plot_quantification", 'figure'),
    Input("mainDiv", "id")
)
def plotRelativeAbundanceProteo(input):

    proteoformsBrno = [proteo.getModificationBrno() for proteo in exp.proteoforms.values() if proteo.getProteoformTotalIntens() > 0]
    proteoformsIntens = [proteo.getProteoformTotalIntens() for proteo in exp.proteoforms.values() if proteo.getProteoformTotalIntens() > 0]

    proteoformsBrno.append("proteoform0")
    exp.proteoform0.setProteoformTotalIntens()
    proteoformsIntens.append(exp.proteoform0.getProteoformTotalIntens())

    proteoformsRatio = [proteoIntens/sum(proteoformsIntens) for proteoIntens in proteoformsIntens]
    
    
    

    fig = px.bar(x=proteoformsBrno, y=proteoformsRatio)
    fig.update_layout(template=template,height=800)
    return fig


# ---------------------------------------------------------------------------- #
#                                   Mics plot                                  #
# ---------------------------------------------------------------------------- #


@app.callback(
    Output('misc_plot', 'figure'),
    Input('mainDiv', "id")
    # Input('input_max_mz', "value")
)
def envelopeParams(minMaxMz):

    envS = [proteoform.envelope.fittedParam[1] for proteoform in exp.proteoforms.values() if proteoform.envelope != None]
    envM = [proteoform.envelope.fittedParam[0] for proteoform in exp.proteoforms.values() if proteoform.envelope != None]

    fig = go.Figure()
    fig.add_scatter( x=envM, y=envS, mode='markers', marker=dict(size=4, color="red"), name='')
    
    fig.update_layout(template=template)
    return fig  


# ---------------------------------------------------------------------------- #
#                                     start                                    #
# ---------------------------------------------------------------------------- #

app.run_server(debug=True)


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #

