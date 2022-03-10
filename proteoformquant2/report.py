

import dash
from dash import html, dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


import pickle

from Classes.elution_profile import ElutionProfile
import numpy as np
import dash_bootstrap_components as dbc
import plotly.io as pio
import json
from Utils import misc
from Utils import constant
import plotly.express as px
import math

from scipy import stats
from kneed import KneeLocator
import pandas as pd
from matplotlib.pyplot import figure, cm
import plotly.graph_objects as go

# ---------------------------------------------------------------------------- #
#                                  Data import                                 #
# ---------------------------------------------------------------------------- #

with open('pfq_out_obj_test_1a.pkl', 'rb') as inp:
    exp = pickle.load(inp) 

print(exp.get_dataset_metrics())


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
        html.Label('Proteoforms'),
        dcc.Dropdown(id="dropdown_proteoforms", multi=True),
        dcc.Graph(id = 'plot_elution_profiles')]
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
        dcc.Graph(id = 'misc_plot', style={'display': 'inline-block'}),
        dcc.Graph(id = 'misc_plot2', style={'display': 'inline-block'}),
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
               Input("plot_elution_profiles", "clickData"),
               Input("plot_all_enveloppes_3d", "clickData"),
               Input("close", "n_clicks")],
              [State("modal", "is_open"),
               State("modal", "children")])


def popup(v1, v2, v3, v4, v5, clicked,is_open,children):
    ctx = dash.callback_context

    print( ctx.triggered[0]["value"])

    if ctx.triggered[0]['prop_id'] == 'close.n_clicks':
        # you pressed the closed button, keeping the modal children as is, and 
        # close the modal itself. 
        return children, False 

    elif len(ctx.triggered) == 0:
        raise dash.exceptions.PreventUpdate

    elif ctx.triggered[0]["value"] == None:
        raise dash.exceptions.PreventUpdate

    elif "points" in ctx.triggered[0]["value"] :
        # you clicked in the graph, returning the modal children and opening it
        try:
            key = ctx.triggered[0]["value"]["points"][0]["customdata"]
            #get spectrum info string
            spectrum = exp.spectra[key]
            str_info = []
            str_info.append(spectrum.get_id())
            str_info.append(html.Br())
            str_info.append("Residuals from multiple proteo quant")
            str_info.append(spectrum.quant_residuals)
            str_info.append(html.Br())
     
            for psm in spectrum.psms:
                str_info.append("Rank: {0}, Proteoform: {1}, isValidated: {2}, ratio: {3}".format(psm.rank, psm.proteoform.get_modification_brno(), psm.isValidated, psm.ratio))
                str_info.append(html.Br())

            return [dbc.ModalHeader("Spectrum"),
                    dbc.ModalBody(
                    html.Div([
                        html.P(str_info,style = {'fontFamily':'monospace'}),
                        html.Div(children=[
                            dcc.Graph(figure = get_table(spectrum.unique_matrix_r.transpose(), spectrum.proteforms_multiquant ),id = 'table_unique', style={'display': 'inline-block'}),
                            dcc.Graph(figure = get_table(spectrum.intensity_matrix_r.transpose(), spectrum.proteforms_multiquant),id = 'table_intensity', style={'display': 'inline-block'})
                        ]),
                        dcc.Graph(figure = spectrum_plot(spectrum),id = 'plot_spectrum')      
                    ])),
                    dbc.ModalFooter(dbc.Button("Close", id="close"))
                    ], True

        except KeyError: #In case click data is not acorresponding to a spectrum try to display proteoform info
            key = ctx.triggered[0]["value"]["points"][0]["customdata"]
            proteoform = exp.proteoforms[key]
            
            str_info = []
            str_info.append(proteoform.get_elution_profile().score_fitted)
            str_info.append(html.Br())
            

            return [dbc.ModalHeader("Proteoform"),
                    dbc.ModalBody(
                    html.Div([
                        html.P(str_info,style = {'fontFamily':'monospace'}),
                        dcc.Graph(figure = ms2_chromatogram_plot(proteoform),id = 'plot_spectrum')      
                    ])),
                    dbc.ModalFooter(dbc.Button("Close", id="close"))
                    ], True

        except KeyError:
            print("Cannot display additional information for that object")

    else:
        raise dash.exceptions.PreventUpdate

def get_table(mat, proteos):

    fig = go.Figure(data=[go.Table(
    header=dict(values=proteos,
                align='left'),
    cells=dict(values=mat, 
               line_color='darkslategray',
               fill_color='lightcyan',
               align='left'))
    ])

    fig.update_layout(template=template)
    return fig

def spectrum_plot(spectrum):
    m = 0
    #plot PSM Annotation

    fig = go.Figure()
    for psm in spectrum.get_validated_psm():
        print(m)
        # get all frag type mz and intens in singles lists
        fragMz = [j for i in [fragType["mzTheo"] for fragType in psm.annotation.values()] for j in i]
        fragIntens =  [j for i in [fragType["intens"] for fragType in psm.annotation.values()] for j in i]
        fragIntens = [i+(m*10000) for i in fragIntens]
        fig.add_scatter( x=fragMz , y=fragIntens, mode='markers',marker_symbol=m, marker=dict( size=7), name = psm.get_modification_brno() )
        m += 1

    for f in range(len(spectrum.fragIntens)):

        fig.add_scatter( x=[spectrum.fragMz[f],spectrum.fragMz[f]] , y=[0, spectrum.fragIntens[f]], mode='lines', marker=dict( color = "black", size=2) )

    fig.update_layout(template=template)
    return fig

def ms2_chromatogram_plot(proteoform, top_n_frag = 40):

    #get most intense fragments at elution peak:
    mean_elution_profile_model = proteoform.get_elution_profile().get_x_at_max_y() #get RT of elution peak
    linked_psms = [psm for psm in proteoform.get_linked_psm()]
    linked_psms_mz = [psm.spectrum.get_rt() for psm in proteoform.get_linked_psm()]
    peak_psm = linked_psms[misc.find_nearest(linked_psms_mz,mean_elution_profile_model)] #psm at max elution 

    annotated_frag_code_intens = peak_psm.get_annotation_pair_format("fragCode","intens")
    frag_codes, intensities = zip(*annotated_frag_code_intens) #unzip values

    if top_n_frag < len(intensities):
        idx_top_frags = sorted(range(len(intensities)), key=lambda i: intensities[i], reverse=True)[:top_n_frag]
        print(idx_top_frags)
    else: 
        idx_top_frags = range(0,len(frag_codes)-1)
        print(idx_top_frags)

    frag_codes_top = [frag_code for idx, frag_code in enumerate(frag_codes) if idx in idx_top_frags] #get fragment codes of top n most intens fragments
    #data table of frag intensities
    intensities_top_frag = pd.DataFrame(columns = frag_codes_top)
    
    for psm in linked_psms:
        psm_annot_pair = psm.get_annotation_pair_format("fragCode","intens")
        dict_psm_annot_pair = dict(psm_annot_pair)
        intens_row = [ dict_psm_annot_pair[code] if code in dict_psm_annot_pair.keys() else 0 for code in frag_codes_top ]

        intensities_top_frag.loc[len(intensities_top_frag)] = intens_row

    rt_psms = [psm.spectrum.get_rt() for psm in linked_psms]

    print("dim dataframe")
    print(intensities_top_frag.shape)

    print("len retention time")
    print(len(rt_psms))


    fig = go.Figure(data=[go.Heatmap(z= np.transpose(np.asarray(intensities_top_frag.values)), y=frag_codes_top)]) #, 
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
    spectrum_key  = [spectrum for spectrum in exp.spectra.keys() if exp.spectra[spectrum].getPrecMz() > minMz and exp.spectra[spectrum].getPrecMz() < maxMz]
    rt = [spectrum.get_rt() for spectrum in exp.spectra.values() if spectrum.getPrecMz() > minMz and spectrum.getPrecMz() < maxMz]

    fig = go.Figure()
    fig.add_scatter( x=rt, y=precIntens, mode='markers', marker=dict(size=4, color="red"), name='Precursor Intensity',customdata=spectrum_key )
    fig.add_scatter( x=rt, y=annotIntens, mode='markers', marker=dict(size=4, color="blue"), name='Annotated Fragment Summed',customdata=spectrum_key )
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
    spectrum_key  = [spectrum for spectrum in exp.spectra.keys() if exp.spectra[spectrum].getPrecMz() > minMz and exp.spectra[spectrum].getPrecMz() < maxMz]
    fig = go.Figure()
    #fig.add_scatter( x=np.log(annotIntens), y=np.log(precIntens), mode='markers', marker=dict(size=4, color="red"), name='Annotated Fragment Summed',customdata=spectrum_key )
    fig.add_scatter( x=np.log(sumFragIntens), y=np.log(precIntens), mode='markers', marker=dict(size=4, color="red"), name='Fragment Sum',customdata=spectrum_key )
    fig.update_layout(template=template)
    fig.update_xaxes(title_text='log(Fragments Summed Intensity)')
    fig.update_yaxes(title_text='log(Precursor Intensity)')
    return fig


# ---------------------------------------------------------------------------- #
#                             All Elution Profiles                             #
# ---------------------------------------------------------------------------- #

@app.callback(
    Output('plot_all_enveloppes', 'figure'),
    Input('mainDiv', "id")
)


def plotAllEnvelopes(input):

    rt_range = [spectrum.get_rt() for spectrum in exp.spectra.values()]
    fig = go.Figure()


    colors=misc.linear_gradient("#4682B4","#FFB347",len([x for x in exp.proteoforms.values() if x.get_elution_profile() != None]))


    i = 0
    for proteoform in [x for x in exp.proteoforms.values() if x.get_elution_profile() != None]:
        c = colors["hex"][i]
        i += 1
        if proteoform.get_elution_profile() != None:
            env = proteoform.get_elution_profile() #if envelope has been computed add line to the plot

            data_yEnv = list(range(int(min(rt_range)),int(max(rt_range)),1))
            zDataEnv = [proteoform.getMzFirstPsm() for x in data_yEnv]
            proteoform_key = np.repeat(proteoform.get_modification_proforma(), [len(data_yEnv)], axis=0)

            yDataEnvFitted, parametersFitted = list(env.get_y_serie(data_yEnv, method = "fitted"))
            if yDataEnvFitted[0] != None: fig.add_scatter( x=data_yEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4, color=c), name=proteoform.get_modification_brno(), line_shape='spline', customdata=proteoform_key)
            else:
                yDataEnvEstim, parametersEstim = list(env.get_y_serie(data_yEnv, method = "estimated"))
                if yDataEnvEstim[0] != None: fig.add_scatter( x=data_yEnv, y=yDataEnvEstim, mode='lines', marker=dict(size=4, color=c), name=proteoform.get_modification_brno(), line_shape='spline', customdata=proteoform_key )

    
    precIntens = [spectrum.getPrecIntens() for spectrum in exp.spectra.values()]
    rt = [spectrum.get_rt() for spectrum in exp.spectra.values()]
    spectrum_key  = [spectrum for spectrum in exp.spectra.keys()]

    fig.add_scatter( x=rt, y=precIntens, mode='markers', marker=dict(size=4, color="black"), name='Precursor Intensity',customdata=spectrum_key)

 
    fig.update_layout(template=template)              
    return fig  



# ---------------------------------------------------------------------------- #
#                               Elution Profiles                               #
# ---------------------------------------------------------------------------- #

# --------------------------------- dropdown --------------------------------- #
@app.callback(
    Output('dropdown_proteoforms', 'options'),
    Input('mainDiv', "id")
)
def dropdown_proteoforms(input):
    """Create option list with proteoforms for dropdown"""



    options = [{"label": " ".join([proteo[1].get_modification_brno(), str(len(proteo[1].get_validated_linked_psm())), str(proteo[1].get_protein_ids()) ]),
                "value": proteo[0], 
                "sort":len(proteo[1].get_validated_linked_psm()) } for proteo in exp.proteoforms.items() if len(proteo[1].get_linked_psm()) > 4]
    options = sorted(options, key=lambda x: x["sort"], reverse=True)

    options = [{k: v for k, v in d.items() if k != 'sort'} for d in options]

    return options

# ---------------------------------------------------------------------------- #


@app.callback(
    Output('plot_elution_profiles', 'figure'),
    [Input('dropdown_proteoforms', 'value')]
)
def plot_elution_profiles(proteoforms_input):
    #avoid initial callback
    print(proteoforms_input)
    if proteoforms_input is None:
        raise PreventUpdate

    #Get plot boundaries
    x_min_max = [0,0]
    y_min_max = [0,0]
    for proteo in proteoforms_input:
        x_val = [psm.spectrum.get_rt() for psm in exp.proteoforms[proteo].get_linked_psm()]
        if min(x_val) < x_min_max[0]: x_min_max[0] = min(x_val)
        if max(x_val)+100 > x_min_max[1]: x_min_max[1] = max(x_val)
        y_val = [psm.spectrum.getPrecIntens() for psm in exp.proteoforms[proteo].get_linked_psm()]
        if min(y_val) < y_min_max[0]: y_min_max[0] = min(y_val)
        if max(y_val) > y_min_max[1]: y_min_max[1] = max(y_val)



    #Instanciate figure
    fig = go.Figure()
    cols = constant.colors
    cols_n = 0

    #Plot each proteoforms:
    for proteo in proteoforms_input:

        data_x_all = [psm.spectrum.get_rt() for psm in exp.proteoforms[proteo].get_linked_psm()]
        data_y_all = [psm.spectrum.getPrecIntens() for psm in exp.proteoforms[proteo].get_linked_psm()]
        spectrum_key = [psm.spectrum.id for psm in exp.proteoforms[proteo].get_linked_psm()]
        fig.add_scatter( x=data_x_all, y=data_y_all, mode='markers', marker=dict(size=7, color="grey",  opacity=0.5), marker_symbol="x-open", name='Spectrum Intensity unvalid', customdata= spectrum_key)

        
        data_y = [psm.spectrum.get_rt() for psm in exp.proteoforms[proteo].get_validated_linked_psm()]
        data_y_spectrum = [psm.spectrum.getPrecIntens() for psm in exp.proteoforms[proteo].get_validated_linked_psm()]
        data_y_psm = [psm.get_prec_intens_ratio() for psm in exp.proteoforms[proteo].get_validated_linked_psm()]
        spectrum_key  = [psm.spectrum.id for psm in exp.proteoforms[proteo].get_validated_linked_psm()]

        fig.add_scatter( x=data_y, y=data_y_spectrum, mode='markers', marker=dict(size=9, color="grey",  opacity=0.5), name='Spectrum Intensity', customdata= spectrum_key)
        fig.add_scatter( x=data_y, y=data_y_psm, mode='markers', marker=dict(size=8, color=cols[cols_n]), name='PSM Intensity',customdata=spectrum_key )

        #add lines between PSM and spectrum intens points
        for i in range(0,len(data_y),1):
            fig.add_scatter( x=[data_y[i],data_y[i]] , y=[data_y_spectrum[i],data_y_psm[i]], mode='lines', marker=dict(size=1, color="#c9c9c9",  opacity=0.5), line={'dash': 'dash'} )


        elution_profile = exp.proteoforms[proteo].get_elution_profile()
        
        if elution_profile != None: #if elution profile model has been computed add line to the plot
            data_x_elution_profile = list(range(int(x_min_max[0]),int(x_min_max[1]),1))
            data_y_elution_profile_fitted, params_fitted = list(elution_profile.get_y_serie(data_x_elution_profile, method = "fitted"))
            data_y_elution_profile_estimated, params_fitted = list(elution_profile.get_y_serie(data_x_elution_profile, method = "estimated"))
            
            proteoform_key = np.repeat(proteo, [len(data_x_elution_profile)], axis=0)
            if data_y_elution_profile_fitted[0] != None: 
                fig.add_scatter( x=data_x_elution_profile, y=data_y_elution_profile_fitted, mode='lines', marker=dict(size=6, color=cols[cols_n]), name='Fitted Parameters', line_shape='spline', customdata=proteoform_key )
            if data_y_elution_profile_estimated[0] != None: 
                fig.add_scatter( x=data_x_elution_profile, y=data_y_elution_profile_estimated, mode='lines', marker=dict(size=2, color=cols[cols_n]), line={'dash': 'dash'}, name='Estimated Parameters', line_shape='spline',customdata=proteoform_key )

        cols_n += 1 
    
    fig.update_layout(title=go.layout.Title( font=dict(
                family="Courier New, monospace",
                size=10,
            )))
    fig.update_layout(template=template, height=1000)

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
    specFiltRt = [spectrum.get_rt() for spectrum in specFilt]
    specFiltIntens = [spectrum.getPrecIntens() for spectrum in specFilt]
    specFiltKey = [spectrum.get_id() for spectrum in specFilt]

    fig.add_trace(go.Scatter3d( x=specFiltRt, y=specFiltIntens, z=specFiltMz, mode='markers',  marker=dict(size=2, color="grey",opacity=0.6), name="proteoform0",customdata=specFiltKey))
              



    proteoFilt = { proteoName:proteo for (proteoName,proteo) in exp.proteoforms.items() if proteo.getMzFirstPsm() > minMz and proteo.getMzFirstPsm() < maxMz}

    rt_range = exp.getRtRange()
    




    colors=misc.linear_gradient("#4682B4","#FFB347",len([x for x in proteoFilt.values() if x.get_elution_profile() != None]))


    i = 0
    for proteoform in [x for x in proteoFilt.values() if x.get_elution_profile() != None]:

        env = proteoform.get_elution_profile() #if envelope has been computed add line to the plot

        data_yEnv = list(range(int(min(rt_range)),int(max(rt_range)),1))
        zDataEnv = [proteoform.getMzFirstPsm() for x in data_yEnv]

        yDataEnvFitted, parametersFitted = list(env.get_y_serie(data_yEnv, method = "fitted"))
        if yDataEnvFitted[0] != None: 
            fig.add_trace(go.Scatter3d( x=data_yEnv, y=yDataEnvFitted, z=zDataEnv, mode='lines',  marker=dict(color=proteoform.get_color()), name=proteoform.get_modification_brno()))
        else:
            yDataEnvEstim, parametersEstim = list(env.get_y_serie(data_yEnv, method = "estimated"))
            if yDataEnvEstim[0] != None: fig.add_trace(go.Scatter3d( x=data_yEnv, y=yDataEnvEstim, z=zDataEnv, mode='lines',  marker=dict(color=proteoform.get_color()), name=proteoform.get_modification_brno()))
    
    for proteoform in [x for x in proteoFilt.values()]:
        precIntens = [psm.spectrum.getPrecIntens() for psm in proteoform.get_validated_linked_psm()]
        rt = [psm.spectrum.get_rt() for psm in proteoform.get_validated_linked_psm()]
        mz  = [psm.spectrum.getPrecMz() for psm in proteoform.get_validated_linked_psm()]
        spectrum_key  = [psm.spectrum.get_id() for psm in proteoform.get_validated_linked_psm()]

        fig.add_trace(go.Scatter3d( x=rt, y=precIntens, z=mz, mode='markers', marker=dict(size=2, color=proteoform.get_color()),name=proteoform.get_modification_brno() ,customdata=spectrum_key))

    


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

    proteoformsBrno = [proteo.get_modification_brno() for proteo in exp.proteoforms.values() if proteo.get_proteoform_total_intens()]
    proteoformsIntens = [proteo.get_proteoform_total_intens() for proteo in exp.proteoforms.values() if proteo.get_proteoform_total_intens()]

    proteoformsBrno.append("proteoform0")
    exp.proteoform0.update_proteoform_total_intens()
    proteoformsIntens.append(exp.proteoform0.get_proteoform_total_intens())

    proteoformsRatio = [proteoIntens/sum(proteoformsIntens) for proteoIntens in proteoformsIntens]
    
    print(len([p for p in proteoformsRatio if p>0.001]), " have relative abundance above 0.001")
    print(len([p for p in proteoformsRatio if p>0.0001]), " have relative abundance above 0.0001")

    fig = px.bar(x=proteoformsBrno, y=proteoformsRatio)
    fig.update_layout(template=template,height=800)
    fig.update_layout( xaxis={'categoryorder':'total descending'})

    return fig


# ---------------------------------------------------------------------------- #
#                                   Mics plot                                  #
# ---------------------------------------------------------------------------- #


@app.callback(
    Output('misc_plot', 'figure'),
    Input('mainDiv', "id")
    # Input('input_max_mz', "value")
)
def plot_envelopes_parameters(minMaxMz):

    env_s = [proteoform.get_elution_profile().get_parameters_fitted()[1] for proteoform in exp.proteoforms.values() if proteoform.get_elution_profile() != None and proteoform.get_elution_profile().is_parameters_fitted() ]
    env_m = [proteoform.get_elution_profile().get_parameters_fitted()[0] for proteoform in exp.proteoforms.values() if proteoform.get_elution_profile() != None and proteoform.get_elution_profile().is_parameters_fitted() ]
    env_k = [proteoform.get_elution_profile().get_parameters_fitted()[3] for proteoform in exp.proteoforms.values() if proteoform.get_elution_profile() != None and proteoform.get_elution_profile().is_parameters_fitted() ]
    env_a = [proteoform.get_elution_profile().get_parameters_fitted()[2] for proteoform in exp.proteoforms.values() if proteoform.get_elution_profile() != None and proteoform.get_elution_profile().is_parameters_fitted() ]

    fig = go.Figure()
    fig.add_scatter( x=env_m, y=stats.zscore(env_s), mode='markers', marker=dict(size=10, color="red"), name='S')
    fig.add_scatter( x=env_m, y=stats.zscore(env_k), mode='markers', marker=dict(size=10, color="blue"), name='K')
    fig.add_scatter( x=env_m, y=stats.zscore(env_a), mode='markers', marker=dict(size=10, color="green"), name='A')


    fig.update_layout(template=template,height=500, width=500)
    return fig  



@app.callback(
    Output('misc_plot2', 'figure'),
    Input('mainDiv', "id")
    # Input('input_max_mz', "value")
)
def plot_count_proteoform(input):

    #proteocount
    # proteoform_count = {}

    # for proteo in exp.proteoforms.values():
    #     proteoform_count[proteo.get_modification_proforma()] = len(proteo.get_linked_psm())
    
    # fig = px.bar(x=proteoform_count.keys(), y=proteoform_count.values())
    # fig.update_layout(template=template,height=500, width=500)
    # fig.update_layout( xaxis={'categoryorder':'total descending', "ticks":'', "showticklabels":False})

    #all frag mz error
    errors = []
    for proteo in exp.proteoforms.values():
        for psm in proteo.get_linked_psm():
            print(psm.rank)
            for frag_type_name, fragtype in psm.annotation.items():
                if frag_type_name == "z+2":
                    errors = errors + fragtype["mzErrorList"]

    print(len(errors))
    fig = px.histogram(x=errors)


    return fig 


# ---------------------------------------------------------------------------- #
#                                     start                                    #
# ---------------------------------------------------------------------------- #

app.run_server(debug=True)


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #

