# Import dash modules
import dash
from dash import html, dcc
import dash_daq as daq
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

# Import other modules
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

# ---------------------------------------------------------------------------- #
#                                  Data import                                 #
# ---------------------------------------------------------------------------- #

folder = None

with open("pfq_out_obj.pkl", "rb") as inp:
    exp = pickle.load(inp)

print(exp.get_dataset_metrics())
print(exp.ident_fn)

# ----------------------------------- test ----------------------------------- #

# ----------------------------------- test ----------------------------------- #

# ---------------------------------------------------------------------------- #
#                                     MISC                                     #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
#                                    Layout                                    #
# ---------------------------------------------------------------------------- #

app = dash.Dash()
template = "plotly_white"
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div(
    id="main_div",
    children=[
        # .container class is fixed, .container.scalable is scalable
        html.Div(
            className="banner",
            children=[
                # Change App Name here
                html.Div(
                    className="container scalable",
                    children=[
                        # Change App Name here
                        html.H2(
                            id="banner-title",
                            children=[
                                html.A(
                                    "Proteoformquant 2: results analysis and quality control",
                                    href="https://github.com/arthur-grimaud/ProteoformQuant2",
                                    style={
                                        "text-decoration": "none",
                                        "color": "inherit",
                                    },
                                )
                            ],
                        ),
                    ],
                )
            ],
        ),
        html.Div(
            id="body",
            className="container scalable",
            children=[
                html.Div(
                    id="app-container",
                    # className="row",
                    children=[
                        # -------------------------------- Left column ------------------------------- #
                        html.Div(
                            # className="three columns",
                            id="left-column",
                            children=[
                                # -------------- Sample selection -------------- #
                                html.Div(
                                    style={"margin": "10px 0px"},
                                    children=[
                                        html.P(
                                            children="Select a sample:",
                                            style={"margin-left": "0px"},
                                        ),
                                        dcc.Dropdown(id="dropdown_sample"),
                                    ],
                                )
                            ],
                        ),
                        # -------------------------------- Main graph -------------------------------- #
                        html.Div(
                            id="div-graphs",
                            children=dcc.Graph(
                                id="main-graph",
                                figure=dict(layout=dict()),
                            ),
                        ),
                        # -------------------------------- Right column -------------------------------- #
                        html.Div(
                            # className="three columns",
                            id="right-column",
                            children=[
                                # -------------- Sample selection -------------- #
                                html.Div(
                                    style={"margin": "0px 10px"},
                                    children=[
                                        html.P(
                                            children="Select a sample:",
                                            style={"margin-right": "3px"},
                                        ),
                                        dcc.Dropdown(),
                                    ],
                                )
                            ],
                        ),
                    ],
                )
            ],
        ),
    ],
)

# ---------------------------------------------------------------------------- #
#                                   Callbacks                                  #
# ---------------------------------------------------------------------------- #

# --------------------------------- Dropwowns -------------------------------- #


@app.callback(Output("dropdown_sample", "options"), Input("main_div", "id"))
def dropdown_sample(input):
    """
    Returns a option list for dropdown with samples in dataset
    """
    # Get samples to read (local folder if not specified):
    if folder == None:
        path = r"./*.pkl"
        files = glob.glob(path)

    options = [{"label": file, "value": file} for file in files]

    return options


# ---------------------------------------------------------------------------- #
#                                     start                                    #
# ---------------------------------------------------------------------------- #

app.run_server(debug=True)

# ---------------------------------------------------------------------------- #
#                                     utils                                    #
# ---------------------------------------------------------------------------- #
