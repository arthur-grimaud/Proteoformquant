# Import dash modules
from decimal import Overflow
import dash
from dash import html
from dash import dcc
from dash import dash_table
import dash_daq as daq
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import numpy as np
from sklearn.decomposition import PCA

# Import other modules
import base64
import datetime
import io
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
import glob
from io import StringIO
import itertools
import networkx as nx

# ---------------------------------------------------------------------------- #
#                                  VARIABLES                                   #
# ---------------------------------------------------------------------------- #

# The name of the files
file_names = []
# The content of the files
file_contents = []
# pfq objects
global msruns
msruns = []
# active run index
global run_i
run_i = 0

# Sample table colnames
spl_tbl_col = ["file", "condition", "replicate"]
dt_tbl_col = ["file", "Tot MS2", "MS2 w/ID", "Chimer MS2", "Tot proteo", "Tot val proteo"]

# ---------------------------------------------------------------------------- #
#                                     MISC                                     #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
#                                    LAYOUT                                    #
# ---------------------------------------------------------------------------- #

app = dash.Dash()
app = dash.Dash(external_stylesheets=[dbc.themes.LUX])
app.layout = html.Div(
    id="main_div",
    children=[
        # =========================================================================================================================
        html.Div(
            id="page_1",  ############# PAGE 1 ############
            style={
                "width": "21cm",
                "min-height": "29.7cm",
                "padding": "1cm",
                "margin": "1cm auto",
                "border": "1px",
                "border-radius": "5px",
                "background": "white",
                "box-shadow": "0 0 5px rgba(0, 0, 0, 0.1)",
            },
            children=[
                html.Div(
                    children=[html.H1("Proteoformquant QC report", style={"textAlign": "center"})],
                    id="title_div",
                ),
                html.P(
                    "In this app you can import .pckl files from proteoformquant to explore both the processing done by proteoformquant and some visualization on the quantification results",
                    style={"text-align": "left"},
                ),
                html.P(id="placeholder"),
                html.Hr(),
                html.H2(
                    "Data upload and preparation",
                    id="title_page1_div",
                    style={"text-align": "left"},
                ),
                html.P(
                    "Load .pkl file(s) from proteoformquant",
                    style={"text-align": "left"},
                ),
                html.Div(
                    id="upload_div",
                    style={
                        "text-align": "left",
                    },
                    children=[
                        html.Div(
                            id="upload_file_div",
                            style={"width": "85%", "display": "inline-block", "margin-right": "5px"},
                            className="card border-primary mb-3",
                            children=[
                                html.Div(
                                    "Upload Spectra file (in .mgf or .mzML)",
                                    className="card-header",
                                    style={"height": "50px"},
                                ),
                                dcc.Upload(
                                    id="upload_spectra",
                                    className="btn btn-secondary",
                                    children=html.Div(["Drag and Drop or ", html.A("Select Files")]),
                                    # Allow multiple files to be uploaded
                                    multiple=True,
                                    style={"height": "50px"},
                                ),
                                html.Div(
                                    id="output_spectra_upload",
                                    children=[
                                        html.Ul(
                                            children=["None"],
                                            style={"overflow-y": "scroll", "height": "150px"},
                                        )
                                    ],
                                ),
                            ],
                        ),
                        html.Div(
                            id="submit_data_div",
                            style={"width": "5%", "display": "inline-block", "verticalAlign": "top"},
                            children=[
                                html.Button(
                                    "Submit",
                                    id="submit_data",
                                    className="btn btn-outline-primary",
                                    n_clicks=0,
                                    style={"height": "269px"},
                                )
                            ],
                        ),
                    ],
                ),
                html.P(
                    "Once submitting valid files should appear below. The table can be modified by clicking on a cell, this allow to difine replicate and condition per file, these are required for some of the vizualization found below. "
                ),
                html.Div(
                    id="upload_table_div",
                    style={"verticalAlign": "top"},
                    className="card border-primary mb-3",
                    children=[
                        html.Div(
                            "Uploaded PFQ objects summary",
                            className="card-header",
                            style={"height": "50px"},
                        ),
                        dash_table.DataTable(
                            id="table_sample",
                            columns=([{"id": p, "name": p} for p in spl_tbl_col]),
                            data=[],
                            editable=True,
                        ),
                    ],
                ),
            ],
        ),
        # =========================================================================================================================
        html.Div(
            id="page_2",
            style={
                "width": "21cm",
                "min-height": "29.7cm",
                "padding": "1cm",
                "margin": "1cm auto",
                "border": "1px",
                "border-radius": "5px",
                "background": "white",
                "box-shadow": "0 0 5px rgba(0, 0, 0, 0.1)",
            },
            children=[
                html.H2(
                    "Dataset processing overview",
                    style={"text-align": "left"},
                ),
                html.Div(
                    children=[
                        html.P(
                            "Here add a description ",
                            style={"width": "80%", "display" "text-align": "left", "display": "inline-block"},
                        ),
                        html.Div(
                            id="gen_plot_page2_div",
                            style={"width": "5%", "display": "inline-block", "verticalAlign": "top"},
                            children=[
                                html.Button(
                                    "Generate",
                                    id="gen_plot_page2_btn",
                                    className="btn btn-outline-dark",
                                    n_clicks=0,
                                    style={},
                                )
                            ],
                        ),
                    ]
                ),
                html.Div(
                    id="summary_processing_table_div",
                    style={"verticalAlign": "top"},
                    className="card border-primary mb-3",
                    children=[
                        html.Div(
                            "Dataset Metrics",
                            className="card-header",
                            style={"height": "50px"},
                        ),
                        dash_table.DataTable(
                            id="table_dataset_metrics",
                            style_cell={"font_size": "10px"},
                            columns=([{"id": p, "name": p} for p in dt_tbl_col]),
                            data=[],
                            editable=False,
                        ),
                    ],
                ),
                html.Div(
                    id="pca_div",
                    style={"verticalAlign": "top"},
                    className="card border-primary mb-3",
                    children=[
                        html.Div(
                            "PCA on peptidoforms ??features?? ",
                            className="card-header",
                            style={"height": "50px"},
                        ),
                        dcc.Graph(id="pca_plot"),
                    ],
                ),
            ],
        ),
        # =========================================================================================================================
        html.Div(
            id="page_3",
            style={
                "width": "21cm",
                "min-height": "29.7cm",
                "padding": "1cm",
                "margin": "1cm auto",
                "border": "1px",
                "border-radius": "5px",
                "background": "white",
                "box-shadow": "0 0 5px rgba(0, 0, 0, 0.1)",
            },
            children=[
                html.Div(
                    children=[
                        html.H2(
                            "MS run",
                            style={
                                "text-align": "left",
                                "display": "inline-block",
                                "width": "50%",
                            },
                        ),
                        dcc.Dropdown(
                            id="dropdown_sample_page3",
                            className="form-label mt-4",
                            style={
                                "width": "250pt",
                                "display": "inline-block",
                                "verticalAlign": "bottom",
                                "text-align": "left",
                            },
                        ),
                    ]
                ),
                html.P(
                    "Here add a description ",
                    style={
                        "width": "50%",
                        "text-align": "left",
                        "verticalAlign": "botom",
                    },
                ),
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.H6("Proteoform grouping", style={"textAlign": "center"}),
                                dcc.Graph(id="plot_network_isobaric"),
                            ],
                            style={"border": "1px black solid"},
                        ),
                        html.Div(
                            children=[
                                html.H6("Chromatogram overview", style={"textAlign": "center"}),
                                dcc.Graph(id="plot_all_enveloppes"),
                            ],
                            style={"border": "1px black solid"},
                        ),
                    ]
                ),
            ],
        ),
        # =========================================================================================================================
        html.Div(
            id="page_4",
            style={
                "width": "29.7cm",
                "min-height": "21cm",
                "padding": "1cm",
                "margin": "1cm auto",
                "border": "1px",
                "border-radius": "5px",
                "background": "white",
                "box-shadow": "0 0 5px rgba(0, 0, 0, 0.1)",
            },
            children=[
                html.Div(
                    children=[
                        html.H6("Proteoforms' elution profiles", style={"textAlign": "center"}),
                        dcc.Dropdown(
                            id="dropdown_proteoforms", multi=True, style={"backgroundColor": "#FFFFFF"}
                        ),
                        dcc.Graph(id="plot_elution_profiles"),
                    ]
                ),
            ],
        ),
        # ---------------------------------- pop up ---------------------------------- #
        dbc.Modal(
            [
                dbc.ModalHeader("Header", id="modalheader"),
                dbc.ModalBody("This is the content of the modal"),
                dbc.ModalFooter(dbc.Button("Close", id="close", className="ml-auto")),
            ],
            size="xl",
            id="modal",
        ),
        # =========================================================================================================================
        html.Div(
            id="page_5",
            style={
                "width": "21cm",
                "min-height": "29.7cm",
                "padding": "1cm",
                "margin": "1cm auto",
                "border": "1px",
                "border-radius": "5px",
                "background": "white",
                "box-shadow": "0 0 5px rgba(0, 0, 0, 0.1)",
            },
            children=[
                html.Div(
                    children=[
                        html.H6("Elution map", style={"textAlign": "center"}),
                        html.Button(
                            "Generate",
                            id="gen_plot_page5_btn",
                            className="btn btn-outline-dark",
                            n_clicks=0,
                            style={},
                        ),
                        dcc.Graph(id="plot_elution_3d"),
                    ]
                ),
            ],
        ),
        # =========================================================================================================================
        html.Div(
            id="page_6",
            style={
                "width": "29.7cm",
                "min-height": "21cm",
                "padding": "1cm",
                "margin": "1cm auto",
                "border": "1px",
                "border-radius": "5px",
                "background": "white",
                "box-shadow": "0 0 5px rgba(0, 0, 0, 0.1)",
            },
            children=[
                html.Div(
                    children=[
                        html.H6("Quantification", style={"textAlign": "center"}),
                        dcc.Dropdown(
                            id="dropdown_sample_page6", multi=False, style={"backgroundColor": "#FFFFFF"}
                        ),
                        dcc.Graph(id="barplot_peptide_quant"),
                    ]
                ),
            ],
        ),
    ],
)


# ---------------------------------------------------------------------------- #
#                                   Callbacks                                  #
# ---------------------------------------------------------------------------- #


@app.callback(
    Output("output_spectra_upload", "children"),
    Input("upload_spectra", "filename"),
    Input("upload_spectra", "contents"),
)
def load_data(list_of_names, list_of_contents):
    if list_of_names is None:
        raise dash.exceptions.PreventUpdate

    # add file to list if not added already
    for i in range(len(list_of_names)):
        if list_of_names[i] not in file_names:

            file_names.append(list_of_names[i])
            file_contents.append(list_of_contents[i])

    return html.Ul(
        children=[html.Li(i) for i in file_names], style={"overflow-y": "scroll", "height": "150px"}
    )


@app.callback(
    Output("table_sample", "data"),
    Input("submit_data", "n_clicks"),
)
def sample_table(x):
    global msruns
    msruns = []
    for fc in file_contents:
        content_type, content_string = fc.split(",")
        decoded = base64.b64decode(content_string)
        print(io.BytesIO(decoded))
        msruns.append(pickle.load(io.BytesIO(decoded)))

    data_table = [
        {spl_tbl_col[0]: file_names[j], spl_tbl_col[1]: "na", spl_tbl_col[2]: "na"}
        for j in range(len(msruns))
    ]
    print("submitting")
    return data_table


@app.callback(
    Output("table_dataset_metrics", "data"),
    Input("gen_plot_page2_btn", "n_clicks"),
)
def dataset_metrics_table(x):

    data_table = [
        {
            dt_tbl_col[0]: file_names[j],
            dt_tbl_col[1]: len(msruns[j].spectra),
            dt_tbl_col[2]: len([s for s in msruns[j].spectra.values() if s.get_number_validated_psm() > 0]),
            dt_tbl_col[3]: len([s for s in msruns[j].spectra.values() if s.get_number_validated_psm() > 1]),
            dt_tbl_col[4]: len(msruns[j].proteoforms),
            dt_tbl_col[5]: len(
                [p for p in msruns[j].proteoforms.values() if p.get_number_validated_linked_psm() > 0]
            ),
        }
        for j in range(len(msruns))
    ]

    print("gen tbl page 2")
    print(data_table)
    return data_table


@app.callback(
    Output("pca_plot", "figure"),
    Input("gen_plot_page2_btn", "n_clicks"),
)
def dpca_sample_plot(x):

    if x is None:
        raise PreventUpdate

    df_all_quant = pd.DataFrame()

    for j in range(len(msruns)):
        proteo_proforma = [proteo.get_modification_proforma() for proteo in msruns[j].proteoforms.values()]
        proteo_intens_prec = [
            proteo.update_proteoform_total_intens(method="precursor")
            for proteo in msruns[j].proteoforms.values()
        ]
        proteo_intens_prec = [x / sum(proteo_intens_prec) for x in proteo_intens_prec]

        df = pd.DataFrame({"proforma": proteo_proforma, file_names[j]: proteo_intens_prec})

        if j == 0:
            df_all_quant = df
        else:
            df_all_quant = df_all_quant.merge(right=df, how="outer", on="proforma")

    df_all_quant = df_all_quant.replace(np.nan, 0)
    df_all_quant.index = df_all_quant["proforma"]
    df_all_quant = df_all_quant.drop("proforma", axis=1)

    df_all_quant = df_all_quant.transpose()

    matrix_quant = df_all_quant.to_numpy()

    pca = PCA(n_components=2)
    components = pca.fit_transform(matrix_quant)
    print(components)

    fig = px.scatter(components, x=0, y=1, text=file_names)
    fig.update_layout(template="plotly_white")
    return fig


@app.callback(Output("placeholder", "style"), Input("dropdown_sample_page3", "value"))
def update_active_msrun(value):
    global run_i
    run_i = value


@app.callback(
    [Output("dropdown_sample_page3", "options"), Output("dropdown_sample_page6", "options")],
    Input("table_sample", "data"),
)
def update_dropdown_page3(x):

    if x is None or x == []:
        raise PreventUpdate
        "preventing update on update dropdown page3"

    opt = []
    for i in range(len(msruns)):
        opt.append({"label": file_names[i], "value": i})

    print(opt)
    return opt, opt


@app.callback(Output("plot_network_isobaric", "figure"), Input("dropdown_sample_page3", "value"))
def plot_network_isobaric_proteoforms(value):

    if value is None:
        raise dash.exceptions.PreventUpdate

    run = msruns[value]
    G = run.isobaric_proteform_graph
    G_nodes_pos = nx.spring_layout(G)

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G_nodes_pos[edge[0]]
        x1, y1 = G_nodes_pos[edge[1]]

        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, line=dict(width=0.5, color="#888"), hoverinfo="none", mode="lines"
    )

    node_x = []
    node_y = []
    proteoform_key = []
    for node in G.nodes():
        x, y = G_nodes_pos[node]
        node_x.append(x)
        node_y.append(y)
        proteoform_key.append(node)

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        customdata=proteoform_key,
        mode="markers",
        hoverinfo="text",
        marker=dict(
            showscale=True,
            colorscale="Portland",
            color=[],
            size=10,
            colorbar=dict(thickness=15, title="Theoretical mass (DA)", xanchor="left", titleside="right"),
            line_width=2,
        ),
    )

    node_adjacencies = []
    node_text = []

    for node in G.nodes:
        delta_mods_tot = sum(
            [mod["monoisotopicMassDelta"] for mod in run.proteoforms[node].get_modification_dict()]
        )

        node_adjacencies.append(delta_mods_tot)
        node_text.append(node)

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            title="",
            titlefont_size=16,
            showlegend=False,
            hovermode="closest",
            margin=dict(b=20, l=5, r=5, t=40),
            annotations=[
                dict(
                    text="",
                    showarrow=False,
                    xref="paper",
                    yref="paper",
                    x=0.005,
                    y=-0.002,
                )
            ],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        ),
    )
    fig.update_layout(template="plotly_white")

    return fig


@app.callback(Output("plot_all_enveloppes", "figure"), Input("dropdown_sample_page3", "value"))
def plotAllEnvelopes(value):

    if value is None:
        raise PreventUpdate

    print("vaaaa")
    print(value)
    if value is None:
        raise dash.exceptions.PreventUpdate

    run = msruns[value]

    rt_range = [spectrum.get_rt() for spectrum in run.spectra.values()]
    fig = go.Figure()

    colors = misc.linear_gradient(
        "#4682B4",
        "#FFB347",
        len([x for x in run.proteoforms.values() if x.get_elution_profile() != None]),
    )

    i = 0
    for proteoform in [x for x in run.proteoforms.values() if x.get_elution_profile() != None]:
        c = colors["hex"][i]
        i += 1
        if proteoform.get_elution_profile() != None:
            env = proteoform.get_elution_profile()  # if envelope has been computed add line to the plot

            data_yEnv = list(range(int(min(rt_range)), int(max(rt_range)), 1))
            zDataEnv = [proteoform.getMzFirstPsm() for x in data_yEnv]
            proteoform_key = np.repeat(proteoform.get_modification_proforma(), [len(data_yEnv)], axis=0)

            yDataEnvFitted, parametersFitted = list(env.get_y_serie(data_yEnv, method="fitted"))
            if yDataEnvFitted[0] != None:
                fig.add_scatter(
                    x=data_yEnv,
                    y=yDataEnvFitted,
                    mode="lines",
                    marker=dict(size=4, color=c),
                    name=proteoform.get_modification_brno(),
                    line_shape="spline",
                    customdata=proteoform_key,
                )
            else:
                yDataEnvEstim, parametersEstim = list(env.get_y_serie(data_yEnv, method="estimated"))
                if yDataEnvEstim[0] != None:
                    fig.add_scatter(
                        x=data_yEnv,
                        y=yDataEnvEstim,
                        mode="lines",
                        marker=dict(size=4, color=c),
                        name=proteoform.get_modification_brno(),
                        line_shape="spline",
                        customdata=proteoform_key,
                    )

    precIntens = [spectrum.get_prec_intens() for spectrum in run.spectra.values()]
    rt = [spectrum.get_rt() for spectrum in run.spectra.values()]
    spectrum_key = [spectrum for spectrum in run.spectra.keys()]

    fig.add_scatter(
        x=rt,
        y=precIntens,
        mode="markers",
        marker=dict(size=4, color="black"),
        name="Precursor Intensity",
        customdata=spectrum_key,
    )
    fig.update_layout(template="plotly_white")

    return fig


########Popup##########


@app.callback(
    [Output("modal", "children"), Output("modal", "is_open")],
    [
        # Input("line_plot_global_intens", "clickData"),
        Input("plot_network_isobaric", "clickData"),
        Input("plot_all_enveloppes", "clickData"),
        Input("plot_elution_profiles", "clickData"),
        # Input("plot_all_enveloppes_3d", "clickData"),
        # Input("plot_envelopes_parameters", "clickData"),
        Input("close", "n_clicks"),
    ],
    [State("modal", "is_open"), State("modal", "children")],
)
def popup(v1, v2, v3, clicked, is_open, children):
    ctx = dash.callback_context

    # print( ctx.triggered[0]["value"])

    print(ctx.triggered[0]["value"])

    if ctx.triggered[0]["prop_id"] == "close.n_clicks":
        # you pressed the closed button, keeping the modal children as is, and
        # close the modal itself.
        return children, False

    elif len(ctx.triggered) == 0:
        raise dash.exceptions.PreventUpdate

    elif ctx.triggered[0]["value"] == None:
        raise dash.exceptions.PreventUpdate

    elif "points" in ctx.triggered[0]["value"]:
        # you clicked in the graph, returning the modal children and opening it
        try:
            key = ctx.triggered[0]["value"]["points"][0]["customdata"]
            # get spectrum info string
            spectrum = msruns[run_i].spectra[key]
            str_info = []
            str_info.append(spectrum.get_id())
            str_info.append(html.Br())
            str_info.append(f"Precursor mz: {spectrum.get_prec_mz()}")
            str_info.append(html.Br())
            str_info.append(f"Precursor Intensity:spectrum.get_prec_intens()")
            str_info.append(html.Br())
            str_info.append("Residuals from multiple proteo quant")
            str_info.append(spectrum.quant_residuals)
            str_info.append(html.Br())
            str_info.append("Missing determining framgent ions")
            str_info.append(str(spectrum.miss_determining_ions))
            str_info.append(html.Br())

            for psm in spectrum.get_psms():
                str_info.append(
                    "Rank: {0}, Proteoform: {1}, is_validated: {2}, ratio: {3}".format(
                        psm.rank,
                        psm.proteoform.get_modification_brno(),
                        psm.is_validated,
                        psm.ratio,
                    )
                )
                str_info.append(html.Br())

            return [
                dbc.ModalHeader("Spectrum"),
                dbc.ModalBody(
                    html.Div(
                        [
                            html.P(str_info, style={"fontFamily": "monospace"}),
                            dcc.Graph(figure=fragments_plot(spectrum), id="fragment_plot"),
                            html.Div(
                                children=[
                                    dcc.Graph(
                                        figure=get_table(
                                            spectrum.unique_matrix_r,
                                            spectrum.proteforms_multiquant,
                                        ),
                                        id="table_unique",
                                    ),
                                    dcc.Graph(
                                        figure=get_table(
                                            spectrum.intensity_matrix_r,
                                            spectrum.proteforms_multiquant,
                                        ),
                                        id="table_intensity",
                                    ),
                                ]
                            ),
                            dcc.Graph(figure=spectrum_plot(spectrum), id="plot_spectrum"),
                        ]
                    )
                ),
                dbc.ModalFooter(dbc.Button("Close", id="close")),
            ], True

        except (
            KeyError,
            TypeError,
            AttributeError,
        ):  # In case click data is not acorresponding to a spectrum try to display proteoform info
            key = ctx.triggered[0]["value"]["points"][0]["customdata"]
            proteoform = msruns[run_i].proteoforms[key]
            print(proteoform)
            str_info = []
            str_info.append(proteoform.get_modification_brno())
            str_info.append(html.Br())
            str_info.append(proteoform.get_sequence())
            str_info.append(html.Br())
            str_info.append(
                f"This proteoform has been identified in {len(proteoform.get_linked_psm())} spectra"
            )
            str_info.append(html.Br())
            str_info.append("proteoform EP fit score:  ")
            str_info.append(proteoform.get_fit_score())
            str_info.append(html.Br())
            str_info.append("proteoform EP coverage:  ")
            str_info.append(proteoform.get_coverage())
            str_info.append(html.Br())
            str_info.append("proteoform gap:  ")
            str_info.append(proteoform.score_gap)
            str_info.append(html.Br())

            return [
                dbc.ModalHeader("Proteoform"),
                dbc.ModalBody(
                    html.Div(
                        [
                            html.P(str_info, style={"fontFamily": "monospace"}),
                            dcc.Graph(
                                figure=plot_psms_rank_distribution(proteoform, msruns[run_i].max_rank),
                                id="plot_psms_rank_distribution",
                            ),
                            dcc.Graph(
                                figure=ms2_chromatogram_plot(proteoform),
                                id="plot_spectrum",
                            ),
                            html.P("Expected vs predicted", style={"fontFamily": "monospace"}),
                            dcc.Graph(
                                figure=fit_expect_predict(proteoform),
                                id="plot_predict_expect",
                            ),
                        ]
                    )
                ),
                dbc.ModalFooter(dbc.Button("Close", id="close")),
            ], True

        except KeyError:
            print("Cannot display additional information for that object")

    else:
        raise dash.exceptions.PreventUpdate


def get_table(mat, proteos):
    print(mat)

    if mat is None:
        fig = go.Figure()
        fig.update_layout(template="plotly_white")
        return fig

    if mat.any() == True:

        fig = go.Figure(
            data=[
                go.Table(
                    header=dict(values=proteos, align="left"),
                    cells=dict(
                        values=mat.transpose(),
                        line_color="darkslategray",
                        fill_color="lightcyan",
                        align="left",
                    ),
                )
            ]
        )
    else:
        fig = go.Figure()

    fig.update_layout(template="plotly_white")
    return fig


def fragments_plot(spectrum):

    spectrum.update_psms_ratio(psms_rank=[], verbose=True)

    m = 0
    # plot PSM Annotation

    fig = go.Figure()

    matrix_annot_intens = []
    proteoform_y = []

    if len(spectrum.get_psms()) == 0:
        return go.Figure()

    for psm in spectrum.get_psms():

        # for both directions

        for dir in ["n-term", "c-term"]:

            proteoform_y.append(m)
            m += 1

            full_annot_intens = [
                psm.get_intensity_at_pos(pos, dir) for pos in range(1, len(psm.proteoform.peptideSequence))
            ]

            matrix_annot_intens.append(full_annot_intens)
            aminos = list(psm.proteoform.peptideSequence)
            aminos = [str(aa + str(index + 1)) for index, aa in enumerate(aminos)]

    fig = go.Figure(
        data=[
            go.Heatmap(
                z=np.asarray(matrix_annot_intens),
                x=aminos,
                y=proteoform_y,
                xgap=4,
                ygap=4,
            )
        ]
    )
    fig["layout"]["yaxis"]["autorange"] = "reversed"
    fig.update_layout(template="plotly_white")
    return fig


def spectrum_plot(spectrum):
    m = 0
    # plot PSM Annotation

    fig = go.Figure()
    for psm in spectrum.get_validated_psm():
        # print(m)
        # get all frag type mz and intens in singles lists
        fragMz = [j for i in [fragType["mzTheo"] for fragType in psm.annotation.values()] for j in i]
        fragIntens = [j for i in [fragType["intens"] for fragType in psm.annotation.values()] for j in i]
        fragIntens = [i + (m * 1000) for i in fragIntens]
        frag_code = [j for i in [fragType["fragCode"] for fragType in psm.annotation.values()] for j in i]
        # print(frag_code)

        fig.add_scatter(
            x=fragMz,
            y=fragIntens,
            mode="markers",
            marker_symbol=m,
            text=frag_code,
            marker=dict(size=7),
            name=psm.get_modification_brno(),
        )
        fig.update_traces(textposition="top center")
        m += 1

    for f in range(len(spectrum.fragIntens)):

        fig.add_scatter(
            x=[spectrum.fragMz[f], spectrum.fragMz[f]],
            y=[0, spectrum.fragIntens[f]],
            mode="lines",
            marker=dict(color="black", size=2),
        )

    fig.update_layout(template="plotly_white")
    return fig


def plot_psms_rank_distribution(proteoform, max_rank):

    ranks = []
    counts = []

    for rank in range(1, max_rank):
        ranks.append(rank)
        counts.append(proteoform.get_number_linked_psm_Rx(rank))

    fig = go.Figure()
    fig.add_trace(go.Bar(x=ranks, y=counts, base=0, marker=dict(color=ranks, colorscale="Bluered")))
    return fig


def ms2_chromatogram_plot(proteoform, top_n_frag=40):

    # if no EP:
    if proteoform.get_elution_profile() == None:
        return go.Figure()

    # get most intense fragments at elution peak:
    mean_elution_profile_model = proteoform.get_elution_profile().get_x_at_max_y()  # get RT of elution peak
    linked_psms = [psm for psm in proteoform.get_linked_psm()]
    linked_psms_mz = [psm.spectrum.get_rt() for psm in proteoform.get_linked_psm()]
    peak_psm = linked_psms[
        misc.find_nearest(linked_psms_mz, mean_elution_profile_model)
    ]  # psm at max elution

    annotated_frag_code_intens = peak_psm.get_annotation_pair_format("fragCode", "intens")
    frag_codes, intensities = zip(*annotated_frag_code_intens)  # unzip values

    if top_n_frag < len(intensities):
        idx_top_frags = sorted(range(len(intensities)), key=lambda i: intensities[i], reverse=True)[
            :top_n_frag
        ]
        # print(idx_top_frags)
    else:
        idx_top_frags = range(0, len(frag_codes) - 1)
        # print(idx_top_frags)

    frag_codes_top = [
        frag_code for idx, frag_code in enumerate(frag_codes) if idx in idx_top_frags
    ]  # get fragment codes of top n most intens fragments
    # data table of frag intensities
    intensities_top_frag = pd.DataFrame(columns=frag_codes_top)

    for psm in linked_psms:
        psm_annot_pair = psm.get_annotation_pair_format("fragCode", "intens")
        dict_psm_annot_pair = dict(psm_annot_pair)
        intens_row = [
            dict_psm_annot_pair[code] if code in dict_psm_annot_pair.keys() else 0 for code in frag_codes_top
        ]

        intensities_top_frag.loc[len(intensities_top_frag)] = intens_row

    rt_psms = [psm.spectrum.get_rt() for psm in linked_psms]

    # print("dim dataframe")
    # print(intensities_top_frag.shape)

    # print("len retention time")
    # print(len(rt_psms))

    fig = go.Figure(
        data=[
            go.Heatmap(
                z=np.transpose(np.asarray(intensities_top_frag.values)),
                y=frag_codes_top,
            )
        ]
    )  # ,
    fig.update_layout(template="plotly_white")
    return fig


def fit_expect_predict(proteoform):

    # if no EP:
    if proteoform.get_elution_profile() == None:
        return go.Figure()

    psms = proteoform.get_validated_linked_psm()
    data_x = np.array([psm.spectrum.get_rt() for psm in psms])
    data_y = np.array([psm.get_prec_intens_ratio() for psm in psms])

    data_y_expect, data_y_predict = proteoform.envelope.scoring_return(
        proteoform.envelope.skewnormal, proteoform.envelope.param_fitted, data_x, data_y
    )

    fig = go.Figure()

    fig.add_scatter(
        x=list(range(len(data_y_expect))),
        y=data_y_predict - data_y_expect,
        mode="markers",
        marker=dict(size=7, color="black"),
        name="predicted(yaxis) vs expected(xaxis)",
    )

    fig.update_layout(template="plotly_white")
    return fig


########Proteoform selection##########


@app.callback(Output("dropdown_proteoforms", "value"), Input("plot_network_isobaric", "selectedData"))
def update_proteo_selector(value):

    if value is None:
        raise PreventUpdate

    print("VALUES")
    holder = []
    if value:
        for x in value["points"]:
            print(x)
            holder.append(x["customdata"])

    print(holder)
    return list(set(holder))


@app.callback(Output("dropdown_proteoforms", "options"), Input("dropdown_sample_page3", "value"))
def dropdown_proteoforms(value):
    """Create option list with proteoforms for dropdown"""

    if value is None:
        raise PreventUpdate

    run = msruns[value]
    options = [
        {
            "label": " ".join(
                [
                    proteo[1].get_modification_brno(),
                    str(len(proteo[1].get_validated_linked_psm())),
                    str(proteo[1].get_protein_ids()),
                ]
            ),
            "value": proteo[0],
            "sort": len(proteo[1].get_validated_linked_psm()),
        }
        for proteo in run.proteoforms.items()
        if len(proteo[1].get_linked_psm()) > 4
    ]
    options = sorted(options, key=lambda x: x["sort"], reverse=True)

    options = [{k: v for k, v in d.items() if k != "sort"} for d in options]

    return options


@app.callback(Output("plot_elution_profiles", "figure"), [Input("dropdown_proteoforms", "value")])
def plot_elution_profiles(proteoforms_input):
    # avoid initial callback
    # print(proteoforms_input)

    exp = msruns[run_i]

    if proteoforms_input is None:
        raise PreventUpdate

    # Get plot boundaries
    x_min_max = [0, 0]
    y_min_max = [0, 0]
    for proteo in proteoforms_input:
        x_val = [psm.spectrum.get_rt() for psm in exp.proteoforms[proteo].get_linked_psm()]
        if min(x_val) < x_min_max[0]:
            x_min_max[0] = min(x_val)
        if max(x_val) + 100 > x_min_max[1]:
            x_min_max[1] = max(x_val)
        y_val = [psm.spectrum.get_prec_intens() for psm in exp.proteoforms[proteo].get_linked_psm()]
        if min(y_val) < y_min_max[0]:
            y_min_max[0] = min(y_val)
        if max(y_val) > y_min_max[1]:
            y_min_max[1] = max(y_val)

    # Instanciate figure
    fig = go.Figure()
    cols = constant.colors
    cols_n = 0

    # Plot each proteoforms:
    for proteo in proteoforms_input:
        if len(exp.proteoforms[proteo].get_validated_linked_psm()) > 0:
            data_x_all = [psm.spectrum.get_rt() for psm in exp.proteoforms[proteo].get_linked_psm()]
            data_y_all = [psm.spectrum.get_prec_intens() for psm in exp.proteoforms[proteo].get_linked_psm()]
            spectrum_key = [psm.spectrum.id for psm in exp.proteoforms[proteo].get_linked_psm()]
            fig.add_scatter(
                x=data_x_all,
                y=data_y_all,
                mode="markers",
                marker=dict(size=10, color="grey", opacity=0.5),
                marker_symbol="x-open",
                name="Spectrum Intensity unvalid",
                customdata=spectrum_key,
            )

            data_y = [psm.spectrum.get_rt() for psm in exp.proteoforms[proteo].get_validated_linked_psm()]
            data_y_spectrum = [
                psm.spectrum.get_prec_intens() for psm in exp.proteoforms[proteo].get_validated_linked_psm()
            ]
            data_y_psm = [
                psm.get_prec_intens_ratio() for psm in exp.proteoforms[proteo].get_validated_linked_psm()
            ]
            spectrum_key = [psm.spectrum.id for psm in exp.proteoforms[proteo].get_validated_linked_psm()]

            fig.add_scatter(
                x=data_y,
                y=data_y_spectrum,
                mode="markers",
                marker=dict(size=15, color="grey", opacity=0.2),
                name="Spectrum Intensity",
                customdata=spectrum_key,
            )
            fig.add_scatter(
                x=data_y,
                y=data_y_psm,
                mode="markers",
                marker=dict(size=10, color=cols[cols_n], opacity=1),
                name="PSM Intensity",
                customdata=spectrum_key,
            )

            # add lines between PSM and spectrum intens points
            for i in range(0, len(data_y), 1):
                fig.add_scatter(
                    x=[data_y[i], data_y[i]],
                    y=[data_y_spectrum[i], data_y_psm[i]],
                    mode="lines",
                    marker=dict(size=2, color="#c9c9c9", opacity=0.5),
                    line={"dash": "dash"},
                )

            elution_profile = exp.proteoforms[proteo].get_elution_profile()

            if elution_profile != None:  # if elution profile model has been computed add line to the plot
                data_x_elution_profile = list(range(int(x_min_max[0]), int(x_min_max[1]), 1))
                data_y_elution_profile_fitted, params_fitted = list(
                    elution_profile.get_y_serie(data_x_elution_profile, method="fitted")
                )
                data_y_elution_profile_estimated, params_fitted = list(
                    elution_profile.get_y_serie(data_x_elution_profile, method="estimated")
                )

                proteoform_key = np.repeat(proteo, [len(data_x_elution_profile)], axis=0)
                if data_y_elution_profile_fitted[0] != None:
                    fig.add_scatter(
                        x=data_x_elution_profile,
                        y=data_y_elution_profile_fitted,
                        mode="lines",
                        marker=dict(size=10, color=cols[cols_n]),
                        name="Fitted Parameters",
                        line_shape="spline",
                        customdata=proteoform_key,
                    )
                if data_y_elution_profile_estimated[0] != None:
                    fig.add_scatter(
                        x=data_x_elution_profile,
                        y=data_y_elution_profile_estimated,
                        mode="lines",
                        marker=dict(size=3, color=cols[cols_n]),
                        line={"dash": "dash"},
                        name="Estimated Parameters",
                        line_shape="spline",
                        customdata=proteoform_key,
                    )

            cols_n += 1

    fig.update_layout(
        title=go.layout.Title(
            font=dict(
                family="Courier New, monospace",
                size=10,
            )
        )
    )
    fig.update_layout(template="plotly_white", height=1000)

    return fig


@app.callback(
    Output("plot_elution_3d", "figure"),
    Input("3dplot-toggle-switch", "value"),
)
def plotAllEnvelopes3d(is_displayed):

    if is_displayed == False:
        raise PreventUpdate

    fig = go.Figure()

    minMz = minMaxMz[0]
    maxMz = minMaxMz[1]

    print("new mz range for 3dplt: ", minMaxMz)

    # Display unassigned spectra:
    specFilt = [
        spectrum
        for spectrum in exp.spectra.values()
        if spectrum.get_prec_mz() > minMz and spectrum.get_prec_mz() < maxMz
    ]

    specFiltMz = [spectrum.get_prec_mz() for spectrum in specFilt]
    specFiltRt = [spectrum.get_rt() for spectrum in specFilt]
    specFiltIntens = [spectrum.get_prec_intens() for spectrum in specFilt]
    specFiltKey = [spectrum.get_id() for spectrum in specFilt]

    fig.add_trace(
        go.Scatter3d(
            x=specFiltRt,
            y=specFiltIntens,
            z=specFiltMz,
            mode="markers",
            marker=dict(size=2, color="grey", opacity=0.6),
            name="proteoform0",
            customdata=specFiltKey,
        )
    )

    proteoFilt = {
        proteoName: proteo
        for (proteoName, proteo) in exp.proteoforms.items()
        if proteo.getTheoPrecMz() > minMz and proteo.getTheoPrecMz() < maxMz
    }

    rt_range = exp.get_rt_range()

    colors = misc.linear_gradient(
        "#4682B4",
        "#FFB347",
        len([x for x in proteoFilt.values() if x.get_elution_profile() != None]),
    )

    i = 0
    for proteoform in [x for x in proteoFilt.values() if x.get_elution_profile() != None]:

        env = proteoform.get_elution_profile()  # if envelope has been computed add line to the plot

        data_yEnv = list(range(int(min(rt_range)), int(max(rt_range)), 1))
        zDataEnv = [proteoform.getMzFirstPsm() for x in data_yEnv]

        yDataEnvFitted, parametersFitted = list(env.get_y_serie(data_yEnv, method="fitted"))
        if yDataEnvFitted[0] != None:
            fig.add_trace(
                go.Scatter3d(
                    x=data_yEnv,
                    y=yDataEnvFitted,
                    z=zDataEnv,
                    mode="lines",
                    marker=dict(color=proteoform.get_color()),
                    name=proteoform.get_modification_brno(),
                )
            )
        else:
            yDataEnvEstim, parametersEstim = list(env.get_y_serie(data_yEnv, method="estimated"))
            if yDataEnvEstim[0] != None:
                fig.add_trace(
                    go.Scatter3d(
                        x=data_yEnv,
                        y=yDataEnvEstim,
                        z=zDataEnv,
                        mode="lines",
                        marker=dict(color=proteoform.get_color()),
                        name=proteoform.get_modification_brno(),
                    )
                )

    for proteoform in [x for x in proteoFilt.values()]:
        precIntens = [psm.spectrum.get_prec_intens() for psm in proteoform.get_validated_linked_psm()]
        rt = [psm.spectrum.get_rt() for psm in proteoform.get_validated_linked_psm()]
        mz = [psm.spectrum.get_prec_mz() for psm in proteoform.get_validated_linked_psm()]
        spectrum_key = [psm.spectrum.get_id() for psm in proteoform.get_validated_linked_psm()]

        fig.add_trace(
            go.Scatter3d(
                x=rt,
                y=precIntens,
                z=mz,
                mode="markers",
                marker=dict(size=2, color=proteoform.get_color()),
                name=proteoform.get_modification_brno(),
                customdata=spectrum_key,
            )
        )

    fig.update_layout(template=template, height=1000, scene=dict(zaxis=dict(range=[minMz, maxMz])))

    return fig


@app.callback(Output("barplot_peptide_quant", "figure"), Input("dropdown_sample_page6", "value"))
def plot_relab_peptide(value):

    if value is None:
        raise PreventUpdate

    run = msruns[value]

    proteoformsBrno = [proteo.get_modification_brno() for proteo in run.proteoforms.values()]
    proteoformsIntens_prec = [
        proteo.update_proteoform_total_intens(method="precursor") for proteo in run.proteoforms.values()
    ]

    proteoformsRatio_prec = [
        proteoIntens / sum(proteoformsIntens_prec) for proteoIntens in proteoformsIntens_prec
    ]

    for proteo in run.proteoforms.values():
        print(proteo.get_protein_ids())

    proteoforms_protein = [";".join(proteo.get_protein_ids()) for proteo in run.proteoforms.values()]

    print(proteoforms_protein)

    print(
        len([p for p in proteoformsRatio_prec if p > 0.001]),
        " have relative abundance above 0.001",
    )
    print(
        len([p for p in proteoformsRatio_prec if p > 0.0001]),
        " have relative abundance above 0.0001",
    )

    fig = px.bar(x=proteoformsBrno, y=proteoformsRatio_prec, color=proteoforms_protein, barmode="group")
    fig.update_layout(template="plotly_white", height=800)
    fig.update_layout(xaxis={"categoryorder": "total descending"})

    return fig


# ---------------------------------------------------------------------------- #
#                                      CSS                                     #
# ---------------------------------------------------------------------------- #

# app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

# ---------------------------------------------------------------------------- #
#                                     start                                    #
# ---------------------------------------------------------------------------- #

app.run_server(debug=True)


# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
