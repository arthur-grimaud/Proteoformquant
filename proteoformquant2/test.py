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
from statistics import mean
import operator
from itertools import compress
import itertools
from scipy.optimize import nnls
import itertools
import scipy as sc

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

# ----------------------------- DATA PREPARATION ----------------------------- #

with open("test_save_1_2.pkl", "rb") as inp:
    run = pickle.load(inp)


# proteoform_to_keep = [
#     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",
#     "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",
#     # "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",
#     "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",
#     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Trimethyl]KPHRYRPGTVALRE",
# ]

proteoform_to_keep = [
    "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",  # K36ac
    "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",  # K27ac
    "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",  # K27me3
    "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Acetyl]PHRYRPGTVALRE",  # K37ac
    "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",  # K23me3
    "ARTKQTARKSTGGKAPRKQLATK[Acetyl]AARKSAPATGGVKKPHRYRPGTVALRE",  # K23ac
    # "ARTKQTARK[Acetyl]STGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",  # K9ac
    "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Trimethyl]KPHRYRPGTVALRE",  # K36me3
    # "ARTKQTARKSTGGKAPRK[Acetyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",  # K18ac
    # "ARTKQTARKSTGGK[Acetyl]APRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",
    # "ARTKQTARKSTGGKAPRK[Trimethyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",
    # "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Trimethyl]PHRYRPGTVALRE",  # K37me3
]


p_del_list = []
for key in run.proteoforms.keys():
    if key not in proteoform_to_keep:
        p_del_list.append(key)
        for psm in run.proteoforms[key].get_linked_psm():
            psm.exclude()

for key in p_del_list:
    del run.proteoforms[key]
print("deleted proteo")


for proteo in proteoform_to_keep:
    proteo = run.proteoforms[proteo]
    print(proteo.get_modification_brno())
    print(proteo.get_number_linked_psm_Rx(rank=1))
    print(proteo.get_number_linked_psm_Rx(rank=2))
    # print(proteo.get_linked_psm())


for spectrum in run.spectra.values():
    spectrum.psms = [i for i in spectrum.get_psms() if i != 0]

# Get list of proteoform and spectra objects in the group
run.proteoform_subset = []
run.spectra_subset = []
run.variables = []

for proforma in proteoform_to_keep:
    proteoform = run.proteoforms[proforma]
    run.proteoform_subset.append(proteoform)

    # get starting values (center of mass based on psm rank and precursor intensity)
    rt_center = proteoform.get_rt_center()
    run.variables.append(rt_center)
    run.variables.append(rt_center)

    for psm in proteoform.get_linked_psm():
        if psm.spectrum not in run.spectra_subset:
            run.spectra_subset.append(psm.spectrum)

with open("test_res_1.pkl", "wb") as outp:
    pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
# ---------------------------------- TESTING --------------------------------- #


bounds = [4450, 4700, 4550, 4700, 4760, 5000, 4600, 4700]

print(proteoform_to_keep)

for g in range(0, 30):

    bounds[6] += 8
    bounds[7] += 8
    # print(bounds)

    run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
    run.update_psms_ratio_subset(run.spectra_subset)
    run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

    score_SI = run.score_signal_explained(run.proteoform_subset, run.spectra_subset, bounds)
    score_EP = run.score_elution_profile_quality(
        run.proteoform_subset, run.spectra_subset, bounds, verbose=False
    )
    score_CH = run.score_chimeric_spectra_quality(run.proteoform_subset, run.spectra_subset, bounds)

    # print("score 3 ", score_SI, score_EP, score_CH)
    # print("score mean  ", mean([score_SI, score_EP, score_CH]))

    print(
        g,
        run.proteoforms[proteoform_to_keep[3]].get_coverage_2(),
        run.proteoforms[proteoform_to_keep[3]].get_fit_score(),
        run.get_residuals_subset(run.spectra_subset),
    )

    fig = run.plot_elution_profiles(proteoform_to_keep, count=g)
    fig.write_image("images/fig_" + f"{g:03}" + ".png")


bounds = [4450, 4700, 4550, 4700, 4760, 5000, 4781, 4924]

run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
run.update_psms_ratio_subset(run.spectra_subset)
run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

with open("test_res_1.pkl", "wb") as outp:
    pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

# ------------------------------------- s ------------------------------------ #
# ------------------------------------- s ------------------------------------ #
# # def findIntersection(fun1, fun2, x0):
# #     return fsolve(lambda x: fun1(x) - fun2(x), x0)


# def skewnormal(x, m, s, a, k):
#     u = (x - m) / s
#     return a * ((2 / s) * norm.pdf(u) * norm.cdf(k * u))


# def skewnorm_cdf(x, m, s, a, k):
#     u = (x - m) / s
#     return a * (norm.cdf(u) - 2 * owens_t(u, k))


# # def skewnormal_mp(x, m, s, a, k):
# #     return a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1))


# # def skewnormal_indefinite_integral(x, m, s, a, k):
# #     return (1 / 2) * (a**2) * mp.exp(((k * (m - s)) / s) - mp.sqrt(((x - s) ** 2) * (1 / s) + 1))


# # def skewnormal_w_param(x):
# #     return skewnormal(x, *[0, 0.5, 14, 0.01])


# # def hline(x):
# #     return x


# m = 4495.662952616016
# s = 47.16108707508988
# a = 1577.3011930167195
# k = 3.079590511861099

# x_val = list(np.linspace(0, 5000, 1000))
# y_val = [skewnormal(x, *[m, s, a, k]) for x in x_val]
# y_val_cdf = [skewnorm_cdf(x, *[m, s, a, k]) for x in x_val]


# def EQ(x1x2, m, s, a, k):

#     x1, x2 = x1x2

#     # percentage area
#     total_auc = a
#     bound_auc = skewnorm_cdf(x2, m, s, a, k) - skewnorm_cdf(x1, m, s, a, k)
#     area_95 = total_auc * 0.95 - bound_auc

#     # x1 x2 at same height
#     fx1_fx2_equal = skewnormal(x1, *[m, s, a, k]) - skewnormal(x2, *[m, s, a, k])

#     return (area_95, fx1_fx2_equal)


# res = sc.optimize.fsolve(EQ, [m, m], args=(m, s, a, k), factor=0.05)
# y_val = [skewnormal(x, *[m, s, a, k]) for x in x_val]
# plt.plot(x_val, y_val)
# plt.axvline(res[0])
# plt.axvline(res[1])


# plt.show()


# # results_min = sc.optimize.fsolve(lambda x: skewnormal(x, *[m, s, a, k]) - 10, m - (s * (1 + np.abs(k))))
# # results_min = results_min[0]
# # print(results_min)
# # results_max = sc.optimize.fsolve(lambda x: skewnormal(x, *[m, s, a, k]) - 10, m + (s * (1 + np.abs(k))))
# # results_max = results_max[0]
# # print(results_max)

# # results_min, results_max = 4505.933238795374, 4593.910061378187

# # res = np.linspace(results_min, results_max, 21)
# # print(res)
# # print(len(res))


# # print([0] * 10)


# # print(" integral solve")


# # total = mp.quad(
# #     lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [-mp.inf, mp.inf]
# # )

# # print(total)


# # def eq1(x1, x2, m, s, a, k):
# #     total_auc = mp.quad(
# #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [-mp.inf, mp.inf]
# #     )
# #     bound_auc = mp.quad(
# #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [x1, x2]
# #     )
# #     return total_auc * 0.95 - bound_auc


# # def eq2(x1, x2, m, s, a, k):
# #     return skewnormal(x1, *[m, s, a, k]) - skewnormal(x2, *[m, s, a, k])


# # def EQ(x1x2, m, s, a, k):

# #     x1, x2 = x1x2

# #     # percentage area
# #     total_auc = mp.quad(
# #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [-mp.inf, mp.inf]
# #     )
# #     bound_auc = mp.quad(
# #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [x1, x2]
# #     )
# #     area_95 = total_auc * 0.95 - bound_auc

# #     # x1 x2 at same height
# #     fx1_fx2_equal = skewnormal(x1, *[m, s, a, k]) - skewnormal(x2, *[m, s, a, k])

# #     return (area_95, fx1_fx2_equal)


# # res = sc.optimize.fsolve(EQ, [m - (s * (1 + np.abs(k))), m + (s * (1 + np.abs(k)))], args=(m, s, a, k))
# # print(res)


# # plt.plot(x_val, y_val)
# # plt.axvline(res[0])
# # plt.axvline(res[1])
# # plt.axhline(skewnormal(res[0], *[m, s, a, k]))
# # plt.axhline(skewnormal(res[1], *[m, s, a, k]))
# # plt.show()
