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

import sys

import resource

# print(resource.getrlimit(resource.RLIMIT_STACK))
# print(sys.getrecursionlimit())

max_rec = 0x100000

# May segfault without this line. 0x100 is a guess at the size of each stack frame.
resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
sys.setrecursionlimit(max_rec)

# ----------------------------- DATA PREPARATION ----------------------------- #


with open("save_res_wt_1_1_mascot.pkl", "rb") as inp:
    run = pickle.load(inp)

run.set_proteoform_isobaric_groups()

i = 0
for group in run.proteoform_isobaric_group:

    group_as_list = list(group)

    print(i, " ", run.proteoforms[group_as_list[0]].get_modification_brno())
    i += 1

group_number = 28
group = run.proteoform_isobaric_group[group_number]


print("-----------***--------------")

avg_rank_scores = []
for rank in range(10):
    scrs = []
    for spectrum in run.spectra.values():
        try:
            psm = spectrum.psms[rank]
            scrs.append(psm.__dict__["Mascot:score"])
        except (IndexError):
            pass
    avg_rank_scores.append(mean(scrs))

avg_rank_scores = [x / max(avg_rank_scores) for x in avg_rank_scores]
print(avg_rank_scores)

print("-----------***--------------")
# Variable TODO hard-coded!

g = 0
n_iterations = 40
n_last_bounds = 20
min_ep_score = 0.3
score_ep_proteos = []
all_scores = []
all_bounds = []

run.proteoform_subset = []  # Proteoform objects in the subset
run.spectra_subset = []  # Spectra objects in the subset
run.rt_boundaries = []  # retention time validation range for the proteoforms in the subset

# Define proteoform subset and spectra subset
for proforma in group:
    proteoform = run.proteoforms[proforma]
    run.proteoform_subset.append(proteoform)
    for psm in proteoform.get_linked_psm():
        if psm.spectrum not in run.spectra_subset:
            run.spectra_subset.append(psm.spectrum)

# TODO could be improved (generalize)
# Sort proteoform in subset (based on psms rank count)

weighted_rank = [
    proteoform.get_weighted_number_linked_validated_psm(max_rank=10) for proteoform in run.proteoform_subset
]
n_rank_1 = [proteoform.get_number_linked_psm_Rx(rank=1) for proteoform in run.proteoform_subset]
n_rank_2 = [proteoform.get_number_linked_psm_Rx(rank=2) for proteoform in run.proteoform_subset]
n_rank_3 = [proteoform.get_number_linked_psm_Rx(rank=3) for proteoform in run.proteoform_subset]
n_rank_4 = [proteoform.get_number_linked_psm_Rx(rank=4) for proteoform in run.proteoform_subset]
n_rank_5 = [proteoform.get_number_linked_psm_Rx(rank=5) for proteoform in run.proteoform_subset]
n_rank_6 = [proteoform.get_number_linked_psm_Rx(rank=6) for proteoform in run.proteoform_subset]
n_rank_7 = [proteoform.get_number_linked_psm_Rx(rank=7) for proteoform in run.proteoform_subset]

zipped_rank_proteo = zip(
    weighted_rank,
    n_rank_1,
    n_rank_2,
    n_rank_3,
    n_rank_4,
    n_rank_5,
    n_rank_6,
    n_rank_7,
    group,
    run.proteoform_subset,
)
zipped_proteo = sorted(zipped_rank_proteo, reverse=True)
run.proteoform_subset = [list(tuple)[-1] for tuple in zipped_proteo]

i = 0
for proteoform in run.proteoform_subset:
    print(proteoform.get_modification_brno(), ",", zipped_proteo[i][:6])

    run.rt_boundaries.append([0, 0])
    i += 1

# run.update_proteoform_subset_validation(run.proteoform_subset, run.rt_boundaries)
run.validate_psms_rank_1()
run.update_psms_ratio_subset(run.spectra_subset)
run.update_proteoforms_elution_profile_subset(run.proteoform_subset)
all_rts = sorted([spectrum.get_rt() for spectrum in run.spectra_subset])
fig = run.plot_elution_profiles(run.proteoform_subset, rt_values=all_rts, count=1)
fig.write_image("images/fig_" + "start_test" + ".png")

run.unvalidate_all_psms()
run.update_psms_ratio_subset(run.spectra_subset)

n_iterations = 20
g = 0
n_last_bounds = 10
min_ep_score = 0.01
score_ep_proteos = []
all_scores = []
all_bounds = []
rd_loc = 4
rd_scale = 1
n_proteo_excluded = 0
threshold_stop = 2

print(run.rt_boundaries)
print("ALL RETENTION TIMES IN GROUP:")
print(all_rts)

scores_proteos = np.zeros((n_iterations * 100, 15))
quant_proteos = np.zeros((n_iterations * 100, 15))

for ps in range(len(run.proteoform_subset)):
    print(ps)
    if n_proteo_excluded < threshold_stop and ps <= 5:
        try:
            rt_window = run.proteoform_subset[ps].get_rt_range_centered()
            # get closest spectra rt

            rt_window[0] = min(all_rts, key=lambda x: abs(x - rt_window[0]))
            rt_window[1] = min(all_rts, key=lambda x: abs(x - rt_window[1]))
            run.rt_boundaries[ps] = rt_window
        except (ValueError):
            print("Problem in range")

        print(f"Proteo: { run.proteoform_subset[ps].get_modification_brno()} range: {run.rt_boundaries[ps]}")
        for i in range(n_iterations):
            g += 1

            print([proteo.get_modification_brno() for proteo in run.proteoform_subset[: ps + 1]])
            for p in range(len(run.proteoform_subset[: ps + 1])):

                if run.rt_boundaries[p][0] != 0 and run.rt_boundaries[p][1]:

                    all_scores.append([])
                    all_bounds.append([])
                    all_scores[i].append(i)

                    proteo_obj = run.proteoform_subset[p]

                    # print(run.proteoform_subset[p])
                    # Minimun boundarie mutation:
                    min_spec_rt, min_ep_rt = proteo_obj.get_min_max_rt_range_shift(side="min")
                    if min_spec_rt < min_ep_rt:  # if spectra lower than modeled ep
                        modifier_min = int(np.round(np.random.normal(loc=rd_loc, scale=rd_scale, size=1)[0]))
                    else:
                        modifier_min = int(np.round(np.random.normal(loc=-rd_loc, scale=rd_scale, size=1)[0]))

                    index_rt_min_start = all_rts.index(run.rt_boundaries[p][0])

                    new_index_rt_min = index_rt_min_start + modifier_min
                    if new_index_rt_min >= 0 and new_index_rt_min < len(all_rts):
                        run.rt_boundaries[p][0] = all_rts[new_index_rt_min]

                    # Maximum boundarie mutation:
                    max_spec_rt, max_ep_rt = proteo_obj.get_min_max_rt_range_shift(side="max")
                    if max_spec_rt > max_ep_rt:  # if spectra lower than modeled ep
                        modifier_max = int(np.round(np.random.normal(loc=-rd_loc, scale=rd_scale, size=1)[0]))
                    else:
                        modifier_max = int(np.round(np.random.normal(loc=rd_loc, scale=rd_scale, size=1)[0]))

                    index_rt_max_start = all_rts.index(run.rt_boundaries[p][1])
                    new_index_rt_max = index_rt_max_start + modifier_max
                    if new_index_rt_max >= 0 and new_index_rt_max < len(all_rts):
                        run.rt_boundaries[p][1] = all_rts[new_index_rt_max]

                    # Test score
                    run.update_proteoform_subset_validation(run.proteoform_subset, run.rt_boundaries)
                    run.update_psms_ratio_subset(run.spectra_subset)
                    run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

                    scores_proteos[g][p] = mean([proteo_obj.get_fit_score(), proteo_obj.get_coverage_2()])
                    quant_proteos[g][p] = proteo_obj.update_proteoform_total_intens()
            fig = run.plot_elution_profiles(run.proteoform_subset, rt_values=all_rts, count=g)
            fig.write_image("images/fig_" + f"{group_number:03}" + "_" + f"{g:04}" + ".png")

        # Check last proteoform added:
        scores_last_proteo = scores_proteos[g - n_last_bounds : g, ps]
        scores_last_proteo = list(scores_last_proteo)

        print(scores_last_proteo)

        print(stdev(scores_last_proteo))

        if stdev(scores_last_proteo) > 0.25 or mean(scores_last_proteo) < 0.5:
            run.rt_boundaries[p] = [0, 0]
            print("proteoform has benn unvalidated ")
            n_proteo_excluded += 1

        run.update_proteoform_subset_validation(run.proteoform_subset, run.rt_boundaries)
        run.update_psms_ratio_subset(run.spectra_subset)
        run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

        fig = run.plot_elution_profiles(run.proteoform_subset, rt_values=all_rts, count=g)
        fig.write_image("images/fig_" + f"{group_number:03}" + "_" + f"{g:04}" + "_intermediate" + ".png")

    else:
        break

np.savetxt("scores_test.csv", scores_proteos, delimiter=",")
np.savetxt("quant_test.csv", quant_proteos, delimiter=",")
with open("testings.pkl", "wb") as outp:
    pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

# import itertools

# combinations = list(itertools.product([0, 1], repeat=5))

# start_rg = 3700
# end_rg = 4200


# for combination in combinations:
#     if sum(combination) == 4:
#         for i in range(len(combination)):
#             if combination[i] == 0:
#                 run.rt_boundaries[i] = [0, 0]
#             if combination[i] == 1:
#                 run.rt_boundaries[i] = [start_rg, end_rg]

#         # Initial validation
#         run.update_proteoform_subset_validation(run.proteoform_subset, run.rt_boundaries)
#         run.update_psms_ratio_subset(run.spectra_subset)
#         run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

#         print(combination)
#         print(run.get_residuals_subset(run.spectra_subset))
#         print(run.get_ratio_validated_0(run.spectra_subset))

#         all_rts = sorted([spectrum.get_rt() for spectrum in run.spectra_subset])
#         # print(all_rts)

#         comb_str = "".join(str(combination))
#         fig = run.plot_elution_profiles(run.proteoform_subset, rt_values=all_rts, count=1)
#         fig.write_image("images/fig_" + f"{comb_str}" + ".png")
#         # print(run.get_coverage_score(run.proteoform_subset))


# run.rt_boundaries[0] = [start_rg, 4000]
# run.rt_boundaries[1] = [3970, end_rg]
# run.rt_boundaries[3] = [start_rg, 4000]
# run.rt_boundaries[2] = [3970, end_rg]
# run.rt_boundaries[4] = [start_rg, end_rg]
# run.rt_boundaries[5] = [start_rg, end_rg]

# run.rt_boundaries[6] = [start_rg, end_rg]
# # Initial validation
# run.update_proteoform_subset_validation(run.proteoform_subset, run.rt_boundaries)
# run.update_psms_ratio_subset(run.spectra_subset)
# run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

# comb_str = "".join(str(combination))
# fig = run.plot_elution_profiles(run.proteoform_subset, rt_values=all_rts, count=1)
# fig.write_image("images/fig_" + "fin_test" + "_final" + ".png")

# with open("testings.pkl", "wb") as outp:
#     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

#
#
#
# # proteoform_to_keep = [
# #     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",
# #     "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",
# #     # "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",
# #     "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",
# #     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Trimethyl]KPHRYRPGTVALRE",
# # ]

# proteoform_to_keep = [
#     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",  # K36ac
#     "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",  # K27ac
#     "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",  # K27me3
#     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Acetyl]PHRYRPGTVALRE",  # K37ac
#     "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",  # K23me3
#     "ARTKQTARKSTGGKAPRKQLATK[Acetyl]AARKSAPATGGVKKPHRYRPGTVALRE",  # K23ac
#     "ARTKQTARK[Acetyl]STGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",  # K9ac
#     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Trimethyl]KPHRYRPGTVALRE",  # K36me3
#     # "ARTKQTARKSTGGKAPRK[Acetyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",  # K18ac
#     # "ARTKQTARKSTGGK[Acetyl]APRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",
#     # "ARTKQTARKSTGGKAPRK[Trimethyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",
#     # "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Trimethyl]PHRYRPGTVALRE",  # K37me3
# ]


# p_del_list = []
# for key in run.proteoforms.keys():
#     if key not in proteoform_to_keep:
#         p_del_list.append(key)
#         for psm in run.proteoforms[key].get_linked_psm():
#             psm.exclude()

# for key in p_del_list:
#     del run.proteoforms[key]
# print("deleted proteo")


# # for proteo in proteoform_to_keep:
# #     proteo = run.proteoforms[proteo]
# #     print(proteo.get_modification_brno())
# #     print(proteo.get_number_linked_psm_Rx(rank=1))
# #     print(proteo.get_number_linked_psm_Rx(rank=2))
# #     # print(proteo.get_linked_psm())


# for spectrum in run.spectra.values():
#     spectrum.psms = [i for i in spectrum.get_psms() if i != 0]


# n_rank_1 = [run.proteoforms[proteoform].get_number_linked_psm_Rx(rank=1) for proteoform in proteoform_to_keep]
# n_rank_2 = [run.proteoforms[proteoform].get_number_linked_psm_Rx(rank=2) for proteoform in proteoform_to_keep]
# n_rank_3 = [run.proteoforms[proteoform].get_number_linked_psm_Rx(rank=3) for proteoform in proteoform_to_keep]

# zipped_rank_proteo = zip(n_rank_1, n_rank_2, n_rank_3, proteoform_to_keep)
# zipped_proteo = sorted(zipped_rank_proteo, reverse=True)

# print(zipped_proteo)
# proteoform_sorted = [list(tuple)[-1] for tuple in zipped_proteo]

# print(proteoform_sorted)

# # Get list of proteoform and spectra objects in the group
# run.proteoform_subset = []
# run.spectra_subset = []
# run.variables = []

# for proforma in proteoform_sorted:
#     proteoform = run.proteoforms[proforma]
#     run.proteoform_subset.append(proteoform)

#     # get starting values (center of mass based on psm rank and precursor intensity)
#     rt_center = proteoform.get_rt_center()
#     run.variables.append(rt_center)
#     run.variables.append(rt_center)

#     for psm in proteoform.get_linked_psm():
#         if psm.spectrum not in run.spectra_subset:
#             run.spectra_subset.append(psm.spectrum)

# with open("test_res_1.pkl", "wb") as outp:
#     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

# # ----------------------------------- vars ----------------------------------- #
# g = 0
# # -------------------------------- CUSTOM OPTI ------------------------------- #


# bounds = []
# for proteo in proteoform_sorted:
#     proteo_obj = run.proteoforms[proteo]
#     if proteo_obj.get_number_linked_psm_Rx(rank=1) >= 5:
#         bounds.append(proteo_obj.get_rt_range(rank=1))
#     else:
#         bounds.append([0, 0])


# print(bounds)


# run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
# run.update_psms_ratio_subset(run.spectra_subset)
# run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

# fig = run.plot_elution_profiles(proteoform_sorted, count=g)
# fig.write_image("images/fig_" + f"{g:03}" + ".png")


# with open("test_res_1.pkl", "wb") as outp:
#     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)


# score_ep_proteos = []

# for proteo in proteoform_sorted:
#     proteo_obj = run.proteoforms[proteo]
#     score = mean([proteo_obj.get_coverage_2(), proteo_obj.get_fit_score()])
#     score_ep_proteos.append(score)

# print(score_ep_proteos)

# all_rts = sorted([spectrum.get_rt() for spectrum in run.spectra_subset])

# all_scores = []
# all_bounds = []

# for i in range(100):
#     g += 1
#     print(g)
#     all_scores.append([])
#     all_bounds.append([])
#     all_scores[i].append(i)
#     for p in range(len(proteoform_sorted)):
#         if bounds[p][0] != 0:
#             proteo_obj = run.proteoforms[proteoform_sorted[p]]
#             # print(proteoform_sorted[p])
#             # Minimun boundarie mutation:
#             min_spec_rt, min_ep_rt = proteo_obj.get_min_max_rt_range_shift(side="min")
#             if min_spec_rt < min_ep_rt:  # if spectra lower than modeled ep
#                 # modifier_min = np.random.choice([+3, +1, -1, -3], p=[0.3, 0.5, 0.15, 0.05])
#                 modifier_min = int(np.round(np.random.normal(loc=2.5, scale=1.5, size=1)[0]))
#                 print("increase lower bound")
#             else:
#                 # modifier_min = np.random.choice([-3, -1, +1, +3], p=[0.3, 0.5, 0.15, 0.05])
#                 modifier_min = int(np.round(np.random.normal(loc=-2.5, scale=1.5, size=1)[0]))
#                 print("decrease lower bound")

#             index_rt_min_start = all_rts.index(bounds[p][0])
#             try:
#                 print("modifier_min:", modifier_min)
#                 bounds[p][0] = all_rts[index_rt_min_start + modifier_min]
#             except IndexError:
#                 pass

#             # Maximum boundarie mutation:
#             max_spec_rt, max_ep_rt = proteo_obj.get_min_max_rt_range_shift(side="max")
#             if max_spec_rt > max_ep_rt:  # if spectra lower than modeled ep
#                 print("decrease upper bound")
#                 # modifier_max = np.random.choice([-3, -1, +1, +3], p=[0.3, 0.5, 0.15, 0.05])
#                 modifier_max = int(np.round(np.random.normal(loc=-2.5, scale=1.5, size=1)[0]))
#             else:
#                 print("increase upper bound")
#                 # modifier_max = np.random.choice([+3, +1, -1, -3], p=[0.3, 0.5, 0.15, 0.05])
#                 modifier_max = int(np.round(np.random.normal(loc=2.5, scale=1.5, size=1)[0]))

#             index_rt_max_start = all_rts.index(bounds[p][1])
#             try:
#                 bounds[p][1] = all_rts[index_rt_max_start + modifier_max]
#                 print("modifier_max:", modifier_max)
#             except IndexError:
#                 pass

#             # TEst score

#             run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
#             run.update_psms_ratio_subset(run.spectra_subset)
#             run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

#             new_score = mean([proteo_obj.get_coverage_2(), proteo_obj.get_fit_score()])
#             if new_score >= score_ep_proteos[p] - 0.1:  # keep if improve
#                 score_ep_proteos[p] = new_score

#             else:  # undo if worse

#                 bounds[p][0] = all_rts[index_rt_min_start]
#                 bounds[p][1] = all_rts[index_rt_max_start]

#                 run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
#                 run.update_psms_ratio_subset(run.spectra_subset)
#                 run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

#             all_scores[i].append(score_ep_proteos[p])
#             all_scores[i].append(bounds[p][0])
#             all_scores[i].append(bounds[p][1])

#             all_bounds[i].append((bounds[p][0], bounds[p][1]))

#     fig = run.plot_elution_profiles(proteoform_sorted, count=g)
#     fig.write_image("images/fig_" + f"{g:03}" + ".png")

#     # print(fig)

# for p in range(len(proteoform_sorted)):
#     if bounds[p][0] != 0:
#         bounds[p][0] = mean([b[p][0] for b in all_bounds])
#         bounds[p][1] = mean([b[p][1] for b in all_bounds])

#         print(run.proteoforms[proteoform_sorted[p]].get_modification_brno())
#         print([b[p][0] for b in all_bounds])
#         print("lower ", mean([b[p][0] for b in all_bounds]))
#         print("higher ", mean([b[p][1] for b in all_bounds]))


# run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
# run.update_psms_ratio_subset(run.spectra_subset)
# run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

# fig = run.plot_elution_profiles(proteoform_sorted, count=g)
# fig.write_image("images/fig_final" + ".png")


# import csv

# with open("scores_bounds.csv", "w") as f:
#     writer = csv.writer(f)
#     writer.writerows(all_scores)
# # ---------------------------------- TESTING --------------------------------- #

# # # bounds = [4450, 4700, 4550, 4700, 4760, 5000, 4600, 4700]

# # print(proteoform_to_keep)

# # for g in range(0, 30):

# #     bounds[6] += 8
# #     bounds[7] += 8
# #     # print(bounds)

# #     run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
# #     run.update_psms_ratio_subset(run.spectra_subset)
# #     run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

# #     score_SI = run.score_signal_explained(run.proteoform_subset, run.spectra_subset, bounds)
# #     score_EP = run.score_elution_profile_quality(
# #         run.proteoform_subset, run.spectra_subset, bounds, verbose=False
# #     )
# #     score_CH = run.score_chimeric_spectra_quality(run.proteoform_subset, run.spectra_subset, bounds)

# #     # print("score 3 ", score_SI, score_EP, score_CH)
# #     # print("score mean  ", mean([score_SI, score_EP, score_CH]))

# #     print(
# #         g,
# #         run.proteoforms[proteoform_to_keep[3]].get_coverage_2(),
# #         run.proteoforms[proteoform_to_keep[3]].get_fit_score(),
# #         run.get_residuals_subset(run.spectra_subset),
# #     )

# #     fig = run.plot_elution_profiles(proteoform_to_keep, count=g)
# #     fig.write_image("images/fig_" + f"{g:03}" + ".png")


# # bounds = [4450, 4700, 4550, 4700, 4760, 5000, 4781, 4924]

# # run.update_proteoform_subset_validation(run.proteoform_subset, bounds)
# # run.update_psms_ratio_subset(run.spectra_subset)
# # run.update_proteoforms_elution_profile_subset(run.proteoform_subset)

# # with open("test_res_1.pkl", "wb") as outp:
# #     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

# # ------------------------------------- s ------------------------------------ #
# # ------------------------------------- s ------------------------------------ #
# # # def findIntersection(fun1, fun2, x0):
# # #     return fsolve(lambda x: fun1(x) - fun2(x), x0)


# # def skewnormal(x, m, s, a, k):
# #     u = (x - m) / s
# #     return a * ((2 / s) * norm.pdf(u) * norm.cdf(k * u))


# # def skewnorm_cdf(x, m, s, a, k):
# #     u = (x - m) / s
# #     return a * (norm.cdf(u) - 2 * owens_t(u, k))


# # # def skewnormal_mp(x, m, s, a, k):
# # #     return a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1))


# # # def skewnormal_indefinite_integral(x, m, s, a, k):
# # #     return (1 / 2) * (a**2) * mp.exp(((k * (m - s)) / s) - mp.sqrt(((x - s) ** 2) * (1 / s) + 1))


# # # def skewnormal_w_param(x):
# # #     return skewnormal(x, *[0, 0.5, 14, 0.01])


# # # def hline(x):
# # #     return x


# # m = 4495.662952616016
# # s = 47.16108707508988
# # a = 1577.3011930167195
# # k = 3.079590511861099

# # x_val = list(np.linspace(0, 5000, 1000))
# # y_val = [skewnormal(x, *[m, s, a, k]) for x in x_val]
# # y_val_cdf = [skewnorm_cdf(x, *[m, s, a, k]) for x in x_val]


# # def EQ(x1x2, m, s, a, k):

# #     x1, x2 = x1x2

# #     # percentage area
# #     total_auc = a
# #     bound_auc = skewnorm_cdf(x2, m, s, a, k) - skewnorm_cdf(x1, m, s, a, k)
# #     area_95 = total_auc * 0.95 - bound_auc

# #     # x1 x2 at same height
# #     fx1_fx2_equal = skewnormal(x1, *[m, s, a, k]) - skewnormal(x2, *[m, s, a, k])

# #     return (area_95, fx1_fx2_equal)


# # res = sc.optimize.fsolve(EQ, [m, m], args=(m, s, a, k), factor=0.05)
# # y_val = [skewnormal(x, *[m, s, a, k]) for x in x_val]
# # plt.plot(x_val, y_val)
# # plt.axvline(res[0])
# # plt.axvline(res[1])


# # plt.show()


# # # results_min = sc.optimize.fsolve(lambda x: skewnormal(x, *[m, s, a, k]) - 10, m - (s * (1 + np.abs(k))))
# # # results_min = results_min[0]
# # # print(results_min)
# # # results_max = sc.optimize.fsolve(lambda x: skewnormal(x, *[m, s, a, k]) - 10, m + (s * (1 + np.abs(k))))
# # # results_max = results_max[0]
# # # print(results_max)

# # # results_min, results_max = 4505.933238795374, 4593.910061378187

# # # res = np.linspace(results_min, results_max, 21)
# # # print(res)
# # # print(len(res))


# # # print([0] * 10)


# # # print(" integral solve")


# # # total = mp.quad(
# # #     lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [-mp.inf, mp.inf]
# # # )

# # # print(total)


# # # def eq1(x1, x2, m, s, a, k):
# # #     total_auc = mp.quad(
# # #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [-mp.inf, mp.inf]
# # #     )
# # #     bound_auc = mp.quad(
# # #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [x1, x2]
# # #     )
# # #     return total_auc * 0.95 - bound_auc


# # # def eq2(x1, x2, m, s, a, k):
# # #     return skewnormal(x1, *[m, s, a, k]) - skewnormal(x2, *[m, s, a, k])


# # # def EQ(x1x2, m, s, a, k):

# # #     x1, x2 = x1x2

# # #     # percentage area
# # #     total_auc = mp.quad(
# # #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [-mp.inf, mp.inf]
# # #     )
# # #     bound_auc = mp.quad(
# # #         lambda x: a * mp.exp(k * (x - m) / s - mp.sqrt((x - m) / s * (x - m) / s + 1)), [x1, x2]
# # #     )
# # #     area_95 = total_auc * 0.95 - bound_auc

# # #     # x1 x2 at same height
# # #     fx1_fx2_equal = skewnormal(x1, *[m, s, a, k]) - skewnormal(x2, *[m, s, a, k])

# # #     return (area_95, fx1_fx2_equal)


# # # res = sc.optimize.fsolve(EQ, [m - (s * (1 + np.abs(k))), m + (s * (1 + np.abs(k)))], args=(m, s, a, k))
# # # print(res)


# # # plt.plot(x_val, y_val)
# # # plt.axvline(res[0])
# # # plt.axvline(res[1])
# # # plt.axhline(skewnormal(res[0], *[m, s, a, k]))
# # # plt.axhline(skewnormal(res[1], *[m, s, a, k]))
# # # plt.show()
