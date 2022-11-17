from cmath import nan
from re import A
import dash
from dash import html
from dash import dcc
from matplotlib.pyplot import figure
import plotly.graph_objects as go
from dash.dependencies import Input, Output
import pickle
from matplotlib.pyplot import cm
from Classes.elution_profile import ElutionProfile
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
from statistics import mean, stdev
import operator
from itertools import compress
import itertools
from scipy.optimize import nnls
import itertools
import scipy as sc
from scipy.stats import norm
import numpy as np

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_problem
from pymoo.algorithms.soo.nonconvex.es import ES
from pymoo.core.evaluator import Evaluator
from pymoo.core.population import Population
from pymoo.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from scipy.optimize import fsolve
import pylab
import numpy as np
from sympy import symbols, Eq, solve
from scipy.integrate import quad
import mpmath as mp
from scipy.stats import skewnorm
from scipy.stats import norm
from scipy.special import owens_t
import pandas as pd
import sys

import resource

# print(resource.getrlimit(resource.RLIMIT_STACK))
# print(sys.getrecursionlimit())

max_rec = 0x100000

# May segfault without this line. 0x100 is a guess at the size of each stack frame.
resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
sys.setrecursionlimit(max_rec)

# ----------------------------- DATA PREPARATION ----------------------------- #


filenames = [
    "save_res_2022_mix1_rep1",
    "save_res_2022_mix1_rep2",
    "save_res_2022_mix2_rep1",
    "save_res_2022_mix2_rep2",
    "save_res_2022_mix3_rep1",
    "save_res_2022_mix3_rep2",
    "save_res_2022_mix4_rep1",
    "save_res_2022_mix4_rep2",
    "save_res_2022_mix5_rep1",
    "save_res_2022_mix5_rep2",
]


for file in filenames:

    with open(f"{file}.pkl", "rb") as inp:
        run = pickle.load(inp)

    number_bins = 1000

    rt_range = run.get_rt_range()
    rt_subdiv = list(np.linspace(rt_range[0], rt_range[1], number_bins))

    elution_map_df = pd.DataFrame(columns=rt_subdiv)
    print(elution_map_df)

    proforma_l = []
    brno_l = []

    for spectrum in run.spectra.values():
        if len(spectrum.get_psms()) > 0:
            psm = spectrum.get_psms()[0]
            proforma = psm.proteoform.get_modification_proforma()
            brno = psm.proteoform.get_modification_brno()

            rt = spectrum.get_rt()

            if proforma not in proforma_l:

                df2 = pd.DataFrame([[0] * elution_map_df.shape[1]], columns=elution_map_df.columns)
                elution_map_df = elution_map_df.append(df2, ignore_index=True)
                proforma_l.append(proforma)
                brno_l.append(brno)

            i_row = proforma_l.index(proforma)
            i_col_val = min([i for i in rt_subdiv if i >= rt], key=lambda x: abs(x - rt))
            i_col = rt_subdiv.index(i_col_val)

            elution_map_df.iloc[i_row, i_col] = elution_map_df.iloc[i_row, i_col] + spectrum.get_prec_intens()

            print(i_col, " : ", i_row)

            # print(elution_map_df.loc[i_row, i_col])
            # print(
            #     elution_map_df.loc[
            #         i_row,
            #     ]
            # )

    elution_map_df.insert(loc=0, column="brno", value=brno_l)
    elution_map_df.insert(loc=0, column="proforma", value=proforma_l)

    # print(elution_map_df)
    # print(brno_l)
    elution_map_df.to_csv(f"elution_map_{file}.csv")
