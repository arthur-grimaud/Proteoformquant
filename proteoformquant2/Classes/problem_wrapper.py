# GA
import pymoo
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.factory import get_crossover, get_mutation, get_sampling
from pymoo.optimize import minimize
from pymoo.core.problem import Problem
import numpy as np
from statistics import mean
from statistics import median
from multiprocessing.pool import ThreadPool as Pool
from pymoo.core.problem import starmap_parallelized_eval
from pymoo.core.problem import ElementwiseProblem
from sqlalchemy import Constraint


class ProblemWrapper(ElementwiseProblem):
    def __init__(self, run, **kwargs):

        self.run = run

        # Define variables boundaries (min max Retention time in subset)
        self.rt_list = [spectrum.get_rt() for spectrum in run.spectra_subset]
        min_rt = min(self.rt_list)
        max_rt = max(self.rt_list)
        print(f"minimum rt: {min_rt} maximum rt: {max_rt}")

        super().__init__(
            # 2 variables per proteoforms (upper and lower retention time limit)
            n_var=len(run.proteoform_subset) * 2,
            n_obj=1,
            n_constr=0,  # len(run.proteoform_subset),
            xl=min_rt,
            xu=max_rt,
            **kwargs,
        )

    def score_subset(self, variables):

        # Update the proteoform subset based on variables
        self.run.update_proteoform_subset_validation(self.run.proteoform_subset, variables)  # WIP
        self.run.update_psms_ratio_subset(self.run.spectra_subset)
        self.run.update_proteoforms_elution_profile_subset(self.run.proteoform_subset)

        # Compute the scores
        signal_score = self.run.score_signal_explained(
            self.run.proteoform_subset, self.run.spectra_subset, variables
        )
        ep_score = self.run.score_elution_profile_quality(
            self.run.proteoform_subset, self.run.spectra_subset, variables
        )
        chimer_score = self.run.score_chimeric_spectra_quality(
            self.run.proteoform_subset, self.run.spectra_subset, variables
        )
        # print(mean([ep_score, chimer_score]))
        # print([ep_score, chimer_score])
        return mean([signal_score, ep_score, chimer_score])

    def constraints_subset(self, variables):
        const = []
        for var in range(0, len(variables), 2):
            const.append(variables[var] - variables[var + 1])
        return const

    def _evaluate(self, variables, out, *args, **kargs):

        # print(variables)
        scores_array = self.score_subset(variables)
        # constraints_array = self.constraints_subset(variables)
        # print(scores_array)

        out["F"] = scores_array
        # out["G"] = constraints_array
