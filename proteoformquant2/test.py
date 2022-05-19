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
from statistics import mean
import operator
from itertools import compress
import itertools
from scipy.optimize import nnls
import itertools

import numpy as np

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_problem
from pymoo.algorithms.soo.nonconvex.es import ES
from pymoo.core.evaluator import Evaluator
from pymoo.core.population import Population
from pymoo.optimize import minimize


problem = get_problem("zdt2")

# create initial data and set to the population object
X = np.random.random((10, problem.n_var))
print(X)
print(type(X[0]))
pop = Population.new("X", X)
print("XXXX")
print(pop.get("X"))
Evaluator().eval(problem, pop)

algorithm = ES(n_offsprings=50, rule=1.0 / 7.0, sampling=pop)

minimize(problem, algorithm, ("n_gen", 10), seed=1, verbose=True)
