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
# ---------------------------------------------------------------------------- #
#                                  Data import                                 #
# ---------------------------------------------------------------------------- #

with open('pfq_out_obj_1_2.pkl', 'rb') as inp:
    exp = pickle.load(inp) 

# ---------------------------------------------------------------------------- #
#                                    Layout                                    #
# ---------------------------------------------------------------------------- #

app = dash.Dash()
template = "plotly_white"

app.layout = html.Div([


html.Div(children=[
    # dcc.Input(id="input_min_mz", type= "number", placeholder="Minimum Precursor Mz"),
    # dcc.Input(id="input_max_mz", type= "number", placeholder="Maximum Precursor Mz"),
    dcc.RangeSlider( id='range_mz', min=5300, max=5800, step=2, value=[5300, 5800], tooltip={"placement": "bottom", "always_visible": True}),
    html.Label('Chromatogram of precursor and annotated fragments intensities'),
    dcc.Graph(id = 'line_plot_global_intens'),
    dcc.Graph(id = 'line_plot_itens_comp'),
    ]
),

html.Br(),

html.Div([
            dcc.Markdown("""Click on datapoint to display spectrum's PSMs details"""),
            html.Pre(id='click-data')
        ]),

html.Br(),

html.Button('reload', id='reload_all_enveloppes'),
html.Div(children=[
    html.Label('All determined elution profil envelopes'),
    dcc.Graph(id = 'plot_all_enveloppes'),
    ]
),

html.Button('reload', id='reload_envelope'),

html.Br(),

html.Div([
            dcc.Markdown("""Click on datapoint to display spectrum's PSMs details"""),
            html.Pre(id='click-data-2')
        ]),

html.Br(),

html.Div(children=[
    html.Label('Elution profile envelope (Select proteoform 1)'),
    dcc.Dropdown( id="dropdown_envelope_1"),
    dcc.Graph(id = 'envelope_plot_1'),
    dcc.Graph(id = 'envelope_plot_score_1')],
    style={'width': '49%', 'display': 'inline-block'}
    ),

html.Div(children=[
    html.Label('Elution profile envelope (Select proteoform 2)'),
    dcc.Dropdown( id="dropdown_envelope_2"),
    dcc.Graph(id = 'envelope_plot_2')],
    style={'width': '49%', 'display': 'inline-block'}
    ),

html.Br(),

html.Button('reload', id='reload_all_enveloppes_3d'),
html.Div(children=[
    dcc.RangeSlider( id='range_mz_2', min=5300, max=5800, step=2, value=[5300, 5800], tooltip={"placement": "bottom", "always_visible": True}),
    html.Label('All determined elution profil envelope with m/z dimension'),
    dcc.Graph(id = 'plot_all_enveloppes_3d'),
    ]
),

html.Div([
            dcc.Markdown("""Click on datapoint to display spectrum's PSMs details"""),
            html.Pre(id='click-data-3')
        ]),

html.Br(),

html.Br(),

html.Div(children=[
    html.Label('Proteoform ratios'),
    dcc.Graph(id = 'plot_quantification'),
    ]
)

])


# ---------------------------------------------------------------------------- #
#                                   Callbacks                                  #
# ---------------------------------------------------------------------------- #



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
    annotIntens  = [spectrum.getSumIntensAnnotFrag() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]
    spectrumKey  = [spectrum for spectrum in exp.spectra.keys() if exp.spectra[spectrum].getPrecMz() > minMz and exp.spectra[spectrum].getPrecMz() < maxMz]
    fig = go.Figure()
    fig.add_scatter( x=np.log(annotIntens), y=np.log(precIntens), mode='markers', marker=dict(size=4, color="red"), name='Annotated Fragment Summed',customdata=spectrumKey )
    fig.update_layout(template=template)
    return fig

# ----------------- click info ------------------ #

@app.callback(
    Output('click-data', 'children'),
    Input('line_plot_global_intens', 'clickData'))
    # Input('line_plot_itens_comp', 'clickData'))
def display_click_data(clickData):
    if clickData is None:
        raise PreventUpdate
    else:
        key = clickData["points"][0]["customdata"]
        spectrum = exp.spectra[key]
        strInfo = str()
        for psm in spectrum.psms:
            strInfo = strInfo + "Rank: {0}, Proteoform: {1}, isValidated: {2}, ratio: {3} \n".format(psm.rank, psm.proteoform.getModificationBrno(), psm.isValidated, psm.ratio)
        return strInfo


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
    if proteo is None:
        raise PreventUpdate
    else:
        return plotEnvelope(proteo)

@app.callback(
    Output('envelope_plot_2', 'figure'),
    Input('dropdown_envelope_2', 'value')
)
def plotEnvelope2(proteo):
    if proteo is None:
        raise PreventUpdate
    else:
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
        if yDataEnvFitted[0] != None: 
            fig.add_scatter( x=xDataEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4, color="red"), name='Fitted Parameters', line_shape='spline' )


            

            # m, s, a, k =  parametersFitted[0], parametersFitted[1], parametersFitted[2], parametersFitted[3]

            # mean = m + (math.sqrt(2/math.pi)) * ((s*a)/math.sqrt(1+(a**2)) ) 
            # std = s * math.sqrt( (1-( (2*a**2)/((1+a**2)*math.pi) )) )
            # std=s
            # print(mean)
            # print(std)
            
            # fig.add_vline(x=mean) 
            # fig.add_vline(x=mean+std) 
            # fig.add_vline(x=mean+std*2)
            # fig.add_vline(x=mean+std*3)
            # fig.add_vline(x=mean+std*10)
            # fig.add_vline(x=mean+std*20)


    titleText = "Proteoform: {0} <br>Parameters Estimated: {1} <br>KS: {3} <br>Parameters Fitted: {2} <br>KS: {4} ".format(exp.proteoforms[proteo].getModificationBrno(), parametersEstim , parametersFitted, env.KsEstimated, env.KsFitted)

        

    fig.update_layout(title=go.layout.Title(text=titleText, font=dict(
            family="Courier New, monospace",
            size=10,
        )))
    fig.update_layout(template=template)
    return fig  

# -------------------------------- score curve ------------------------------- #

# @app.callback(
#     Output('envelope_plot_score_1', 'figure'),
#     Input('dropdown_envelope_1', 'value')
# )
# def plotScore(proteo):

#     if proteo is None:
#         raise PreventUpdate
#     else:
#         fig = go.Figure()


#         xData = np.array([psm.spectrum.getRt() for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()])
#         yData = np.array([psm.getPrecIntensRatio() for psm in exp.proteoforms[proteo].getValidatedLinkedPsm()])
#         xDataEnv = list(range(int(min(xData)),int(max(xData)),1))


#         #fig.add_scatter( x=xData, y=yData, mode='markers', marker=dict(size=10, color="black"), name='Precursor Intensity' )

#         nRemoved = []
#         splitScores = []
#         for n in range(0,len(xData)-5):
#             if len(xData[:-n]) > 5:
#                 print(n)
                
#                 xDataT = xData[:-n]
#                 yDataT = yData[:-n]

#                 print(xDataT)

#                 env = Envelope(exp.proteoforms[proteo].getValidatedLinkedPsm())
                        
#                 env.estimatedParam, env.fittedParam, env.KsEstimated, env.KsFitted = env.fitSkewNormal(xDataT, yDataT)


#                 if env.KsFitted >0.7:
#                     nRemoved.append(n)
#                     splitScores.append(env.KsFitted)

#                     print(env.fittedParam)
#                     yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
#                     #if yDataEnvFitted[0] != None: 
#                         #fig.add_scatter( x=xDataEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4), name=len(xDataT), line_shape='spline')


#         kn = KneeLocator(nRemoved, splitScores, curve='concave', direction='increasing')
#         print(kn.knee)
#         fig.add_scatter( x=nRemoved, y=splitScores, mode='lines', name='R square' )

#         return fig


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
    
    
    precIntens = [spectrum.getPrecIntens() for spectrum in exp.spectra.values()]
    rt = [spectrum.getRt() for spectrum in exp.spectra.values()]
    spectrumKey  = [spectrum for spectrum in exp.spectra.keys()]

    fig.add_scatter( x=rt, y=precIntens, mode='markers', marker=dict(size=4, color="black"), name='Precursor Intensity',customdata=spectrumKey)

 
    fig.update_layout(template=template)              
    return fig  

# ----------------- click info ------------------ #

@app.callback(
    Output('click-data-2', 'children'),
    Input('plot_all_enveloppes', 'clickData'))
def display_click_data(clickData):

    if clickData is None:
        raise PreventUpdate
    else:
        key = clickData["points"][0]["customdata"]
        spectrum = exp.spectra[key]
        strInfo = str()
        for psm in spectrum.psms:
            strInfo = strInfo + "Rank: {0}, Proteoform: {1}, isValidated: {2}, ratio: {3} \n".format(psm.rank, psm.proteoform.getModificationBrno(), psm.isValidated, psm.ratio)
        return strInfo
# ---------------------------------------------------------------------------- #
#                               all envelopes 3D                               #
# ---------------------------------------------------------------------------- #


@app.callback(
    Output('plot_all_enveloppes_3d', 'figure'),
    Input('range_mz_2', "value")
)
def plotAllEnvelopes3d(minMaxMz):


    minMz = minMaxMz[0]
    maxMz = minMaxMz[1]

    specFilt = { spectrumName:spectrum for (spectrumName,spectrum) in exp.spectra.items() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz}
    proteoFilt = { proteoName:proteo for (proteoName,proteo) in exp.proteoforms.items() if proteo.getMzFirstPsm() > minMz and proteo.getMzFirstPsm() < maxMz}

    rt_range = [spectrum.getRt() for spectrum in specFilt.values()]
    fig = go.Figure()



    colors=misc.linear_gradient("#4682B4","#FFB347",len([x for x in proteoFilt.values() if len(x.envelopes)>0]))


    i = 0
    for proteoform in [x for x in proteoFilt.values() if len(x.envelopes)>0]:
        if len(proteoform.envelopes) > 0:
            for env in proteoform.envelopes: #if envelope has been computed add line to the plot

                xDataEnv = list(range(int(min(rt_range)),int(max(rt_range)),1))
                zDataEnv = [proteoform.getMzFirstPsm() for x in xDataEnv]

                yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
                if yDataEnvFitted[0] != None: fig.add_trace(go.Scatter3d( x=xDataEnv, y=yDataEnvFitted, z=zDataEnv, mode='lines',  marker=dict(color=proteoform.getColor()), name=proteoform.getModificationBrno()))
                else:
                    yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
                    if yDataEnvEstim[0] != None: fig.add_trace(go.Scatter3d( x=xDataEnv, y=yDataEnvEstim, z=zDataEnv, mode='lines',  marker=dict(color=proteoform.getColor()), name=proteoform.getModificationBrno()))
    
    for proteoform in [x for x in proteoFilt.values()]:
        precIntens = [psm.spectrum.getPrecIntens() for psm in proteoform.getValidatedLinkedPsm()]
        rt = [psm.spectrum.getRt() for psm in proteoform.getValidatedLinkedPsm()]
        mz  = [psm.spectrum.getPrecMz() for psm in proteoform.getValidatedLinkedPsm()]
        spectrumKey  = [psm.spectrum.getId() for psm in proteoform.getValidatedLinkedPsm()]

        fig.add_trace(go.Scatter3d( x=rt, y=precIntens, z=mz, mode='markers', marker=dict(size=2, color=proteoform.getColor()),customdata=spectrumKey))

    fig.update_layout(template=template,height=1000)
               
    return fig  

# ----------------- click info ------------------ #

@app.callback(
    Output('click-data-3', 'children'),
    Input('plot_all_enveloppes_3d', 'clickData'))
def display_click_data(clickData):
    if clickData is None:
        raise PreventUpdate
    else:
        key = clickData["points"][0]["customdata"]
        spectrum = exp.spectra[key]
        strInfo = str()
        for psm in spectrum.psms:
            strInfo = strInfo + "Rank: {0}, Proteoform: {1}, isValidated: {2}, ratio: {3} \n".format(psm.rank, psm.proteoform.getModificationBrno(), psm.isValidated, psm.ratio)
        return strInfo


# ---------------------------------------------------------------------------- #
#                               proteoform ratios                              #
# ---------------------------------------------------------------------------- #





@app.callback(
    Output("plot_quantification", 'figure'),
    Input('reload_all_enveloppes_3d', "id")
)
def plotAllEnvelopes3d(input):

    proteoformsBrno = [proteo.getModificationBrno() for proteo in exp.proteoforms.values() if proteo.getProteoformTotalIntens() > 0]
    proteoformsIntens = [proteo.getProteoformTotalIntens() for proteo in exp.proteoforms.values() if proteo.getProteoformTotalIntens() > 0]
    proteoformsRatio = [proteoIntens/sum(proteoformsIntens) for proteoIntens in proteoformsIntens]
    

    fig = px.bar(x=proteoformsBrno, y=proteoformsRatio)
    fig.update_layout(template=template,height=800)
    return fig


# ---------------------------------------------------------------------------- #
#                                     start                                    #
# ---------------------------------------------------------------------------- #

app.run_server(debug=True)


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #

