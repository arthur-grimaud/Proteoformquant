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


class ProblemWrapper(Problem):
    def __init__(self, run, **kwargs):

        self.run = run

        self.rt_list = [spectrum.get_rt() for spectrum in run.spectra_subset]

        min_rt = min(self.rt_list)
        max_rt = max(self.rt_list)

        print(f"minimum rt: {min_rt} maximum rt: {max_rt}")
        # max_rt = np.repeat(max(self.rt_list) + 1, len(run.variables))

        # print(min_rt)
        # print(max_rt)

        # runner = Pool(8)
        # func_eval = starmap_parallelized_eval

        # print("calculated n_var = ", len(run.proteoform_subset) * 2)

        super().__init__(
            n_var=len(run.proteoform_subset) * 2,
            n_obj=1,
            n_constr=0,
            xl=min_rt,
            xu=max_rt,
            # runner=runner,
            # func_eval=func_eval,
            **kwargs,
        )

    def score_subset(self, variables):
        # print(variables)

        self.run.update_psms_validation_subset_2(self.run.proteoform_subset, variables)  # WIP
        self.run.update_psms_ratio_subset(self.run.spectra_subset)
        self.run.update_proteoforms_elution_profile_subset(self.run.proteoform_subset)

        # Mean of elution profiles scores
        u_fit = mean(
            [
                proteoform.get_fit_score()
                for proteoform in self.run.proteoform_subset
                if proteoform.get_fit_score() != None
            ]
        )
        # mean of elution profile evidence
        u_ep_evid = mean(
            [
                proteoform.get_ratio_left_right()
                for proteoform in self.run.proteoform_subset
                if proteoform.get_fit_score() != None
            ]
        )
        # mean of elution profile area range valid
        u_ep_area = mean(
            [
                proteoform.get_boundaries_area_ratio()
                for proteoform in self.run.proteoform_subset
                if proteoform.get_fit_score() != None
            ]
        )
        # mean coverage intensity EP:
        # u_cover = mean(self.run.get_coverage_score(self.run.proteoform_subset))
        # Mean of residual scores
        u_residual = mean([spectra.get_residuals() for spectra in self.run.spectra_subset])
        # Mean of penalities for "gap" in Elution profiles
        u_gap = mean(
            self.run.get_gap_in_validated_psms(self.run.proteoform_subset, self.run.spectra_subset, variables)
        )
        # Mean of rank score
        u_rank = mean(self.run.get_rank_score(self.run.proteoform_subset))

        # ratio of intensities explained:
        r_intens = self.run.get_intensity_explained_ratio(self.run.spectra_subset)
        # ratio multiquant = 0
        # ratio of first rank psm unused
        r_missR1 = self.run.get_ratio_missed_rank_1(self.run.spectra_subset)
        # penality for psm quantified as zero in multiquant
        r_quant0 = self.run.get_ratio_validated_0(self.run.spectra_subset)
        # ratio of unvalidated psm before last validated (ordered by rank)
        r_gap_psm = self.run.get_gap_in_rank_psms(self.run.spectra_subset)

        # print("mean fit score:", u_fit)
        # print("mean res score:", u_residual)

        # print(
        #     ",".join(
        #         [
        #             str(1 - u_fit),
        #             str(u_residual),
        #             str(1 - u_gap),
        #             # + ((r_gap_psm) * 0),
        #             # + ((1 - u_rank) * 1),
        #             str(1 - r_intens),
        #             # + ((r_missR1) * 1),
        #             str(r_quant0),
        #         ]
        #     )
        # )

        return (
            ((1 - u_fit) * 3)
            + ((1 - u_ep_evid) * 2)
            + ((1 - u_ep_area) * 2)
            + (u_residual * 3)
            + ((1 - u_gap) * 1)
            # + ((r_gap_psm) * 0)
            + ((1 - u_rank) * 1)
            + ((1 - r_intens) * 6)
            + ((r_missR1) * 1)
            + ((r_quant0) * 1)
        )

    def _evaluate(self, designs, out, *args, **kargs):
        # print(designs[0])
        res = []
        for design in designs:
            res.append(self.score_subset(design))
            # print("Results")
        # print(res)

        out["F"] = np.array(res)
