# Utility script to store parameters and generate a json file with default parameters

import json
import jsonc_parser
import os
import sys


def parameters():
    """Parameters of the program"""

    return {
        # Tolerances
        "prec_mz_tol": 1.05,  # in Da
        "frag_mz_tol": 0.015,  # in Da
        # Fragments to be considered (see available options in "Utils")
        "fragments_types": ["c", "zdot", "z+1", "z+2"],
        # "fragments_types" : ["b", "y"],
        # Filtering and thresholds:
        "max_rank": 5,  # Maximum rank to be considered (1-based)
        "fdr_threshold": 0.1,  # False-Discovery rate
        "decoy_tag": "REV_",  # Tag for decoy proteins (e.g. REV_ or DECOY_)
        "score_name": "Amanda:AmandaScore",  # Score to be used for FDR control
        "intensity_threshold": 0,
        "elution_profile_score_threshold": 0,
        "min_n_psm": 3,  # Minimal number of PSM for a peptidoform to be considered
        # Optimization
        "n_iter": 25,  # number of iteration for each peptidoform optimization
        "max_rejected_proteo": 3,  # number of peptidoform rejected before stopping optimization
        # Starting range
        "window_size_rt": 400,
        # Param of the normal distribution used for random mutation of the rt_validation_ranges
        "rd_loc": 3,  # deviation of the mean from zero
        "rd_scale": 2,  # spread of the distribution
        # Param for subset optimization:
        "min_spectra_subset": 5,
        # Param for peptidoform validation:
        "min_ep_fit": 0.6,  # minimal score for validation
        "min_ep_cov": 0.4,  # minimal score for validation
        # "min_ep_gap" : 0.7,  # mini
        # Proteform groups:
        "min_connect": 1,  # Minimal number of time two peptidoform must be co identified in order to be in the same proteoform group
        # additional param:
        "only_r1": "False",
    }


def generate_default_param_file():
    """Generate a default parameter file"""
    params = parameters()
    with open("params.jsonc", "w") as f:
        json.dump(params, f, indent=4)
    print("Default parameter file generated: params.jsonc")
