# Import dash modules
from decimal import Overflow
import dash
from dash import html
from dash import dcc
from dash import dash_table
import dash_daq as daq
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

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

# Sample table colnames
spl_tbl_col = ["file", "condition", "replicate"]
dt_tbl_col = ["file", "Total MS2", "MS2 w/ ID", "Chimeric MS2", "Total proteoform", "Total valid proteoform"]

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
                html.P(
                    "Here add a description ",
                    style={"text-align": "left"},
                ),
                html.Div(
                    id="gen_plot_page2_div",
                    style={"width": "5%", "display": "inline-block", "verticalAlign": "top"},
                    children=[
                        html.Button(
                            "Submit",
                            id="gen_plot_page2_btn",
                            className="btn btn-outline-primary",
                            n_clicks=0,
                            style={"height": "10px"},
                        )
                    ],
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
                            columns=([{"id": p, "name": p} for p in dt_tbl_col]),
                            data=[],
                            editable=False,
                        ),
                    ],
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
    print(msruns)
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
    for j in range(len(msruns)):
        print(j)
        print(len([s for s in msruns[j].spectra.values() if s.get_number_validated_psm() > 0]))
    print("gen tbl page 2")
    print(data_table)
    return data_table


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
