from pickle import TRUE
from pyteomics import mass
from logging import warning
import plotly.graph_objs as go
from Classes.elution_profile import ElutionProfile
from Utils import constant
from statistics import mean
import numpy as np


class Proteoform:
    def __init__(
        self,
        peptideSequence,
        modificationBrno,
        modificationProforma,
        modificationDict={},
        protein_ids=[],
    ):

        # Proteoform Description
        self.peptideSequence: str = peptideSequence
        self.modificationDict: dict = modificationDict
        self.modificationBrno: str = modificationBrno
        self.modificationProforma: str = modificationProforma

        # Protein(s) reference
        self.protein_ids = protein_ids

        # Theoretical Fragments
        self.theoFrag = None
        self.fragment_types = None  # types of fragments computedS

        # All PSMs of that proteoform
        self.linkedPsm: list() = []

        # ElutionTime~intensity envelope
        self.envelope = None

        # color
        self.color = None

        # Summed intensity of linked psm
        self.totalIntens = 0

        # Rt range for validating psms
        self.min_bound_rt = None
        self.max_bound_rt = None

    # Getters

    def getTheoFrag(self):
        return self.theoFrag

    def getMzFirstPsm(self):
        return self.linkedPsm[0].getCalculatedMassToCharge()

    def getTheoPrecMz(self):
        return self.linkedPsm[0].getCalculatedMassToCharge()

    def getMzMean(self):
        pass

    def get_rt_center(self):
        """returns rt "center" based on the weigthed average of psm's precursor mass"""
        rts = [psm.spectrum.get_rt() for psm in self.get_linked_psm()]
        weights = [(psm.spectrum.getPrecIntens() + 1) / psm.get_rank() for psm in self.get_linked_psm()]

        numerator = sum([rts[i] * weights[i] for i in range(len(rts))])
        denominator = sum(weights) + 1

        print(round(numerator / denominator, 2))

        return round(numerator / denominator, 2)

        # if len(rts) == 0:
        #     return (
        #         self.get_linked_psm()[0].spectrum.get_rt(),
        #         self.get_linked_psm()[-1].spectrum.get_rt(),
        #     )

        # rt_range = (min(rts), max(rts))
        # return rt_range

    def get_rt_range(self, rank=1):
        """returns the RT range for linked psm of rank <= rank"""
        rts = [psm.spectrum.get_rt() for psm in self.get_linked_psm() if psm.rank <= rank]

        return [min(rts), max(rts)]

    def get_peptide_sequence(self):
        return self.peptideSequence

    def get_modification_dict(self):
        return self.modificationDict

    def get_modification_brno(self):
        return self.modificationBrno

    def get_modification_proforma(self):
        return self.modificationProforma

    def get_linked_psm(self):
        """Return a list of linked PSMs that haven't been excluded"""
        return [psm for psm in self.linkedPsm if psm.is_included()]

    def get_validated_linked_psm(self):
        """Return a curated list of linked psm whose self.is_validated = true"""
        return [psm for psm in self.linkedPsm if psm.is_validated and psm.is_included()]

    def get_number_validated_linked_psm(self):
        """Return a the number of element in curated list of linked psm whose self.is_validated = true"""
        return len([psm for psm in self.linkedPsm if psm.is_validated])

    def get_number_linked_psm_Rx(self, rank):
        """Return a the number of element in curated list of linked psm whose rank is equal to x"""
        return len([psm for psm in self.linkedPsm if psm.rank == rank])

    def get_weighted_number_linked_psm(self, max_rank):
        """Returns the sum of linked PSM weigthed by their rank"""
        weights = range(max_rank + 1, 1, -1)  # score for each psm rank

        return sum([weights[psm.get_rank() - 1] for psm in self.linkedPsm])

    def get_weighted_number_linked_validated_psm(self, max_rank):
        """Returns the sum of linked PSM weigthed by their rank"""
        weights = range(max_rank + 1, 1, -1)  # score for each psm rank

        return sum([weights[psm.get_rank() - 1] for psm in self.get_validated_linked_psm()])

    def get_proteoform_total_intens(self):
        return self.totalIntens

    def get_color(self):
        return self.color

    def get_elution_profile(self):
        return self.envelope

    def get_protein_ids(self):
        return self.protein_ids

    def get_fit_score(self):

        if self.get_number_validated_linked_psm == 0:  # If there is no validated PSM do not return the score
            return 0

        if self.get_elution_profile() == None:
            return 0

        if self.get_elution_profile().is_parameters_fitted() == False:
            return 0

        return self.get_elution_profile().score_fitted

    def get_coverage(self, y_limit, n_quantiles=10):
        """Return the percentage of elution profile quantile for which a psm is found
        n_quantiles: int, number of quantiles (intervals) to compute the score
        y_limit: intensity value to determine the RT range on which to compute the quantiles
        """

        y_limit = 100000

        if (
            self.get_number_validated_linked_psm() == 0
        ):  # If there is no validated PSM do not return the score
            return None

        if self.get_number_validated_linked_psm() < 3:  # If there is only one return worst score
            return 0

        if self.get_elution_profile() == None:
            return 0  # If no fit is found return worst score

        if self.get_elution_profile().is_parameters_fitted() == False:
            return 0  # If no fit is found return worst score

        psms_rt = [psm.spectrum.get_rt() for psm in self.get_validated_linked_psm()]

        n_quantiles = 4

        # print("n quantiles:", n_quantiles)

        range_rt = self.get_elution_profile().get_x_range_for_y(Y=y_limit)

        if range_rt[0] == range_rt[1]:
            # print("no interval")
            return 0

        intervals_rt = np.linspace(range_rt[0], range_rt[1], n_quantiles + 1)

        intervals_count = [0] * n_quantiles

        print("len count", len(intervals_count))
        print("len intervals", len(intervals_rt))
        print("range rt", range_rt)
        print("intervals rt", intervals_rt)

        for i in range(len(intervals_rt) - 1):

            n_psm_in_interval = len(
                [rt for rt in psms_rt if intervals_rt[i] <= rt and rt < intervals_rt[i + 1]]
            )

            intervals_count[i] = n_psm_in_interval

        print(intervals_count)

        return len([i for i in intervals_count if i > 0]) / n_quantiles

    def get_coverage_2(self):

        if (
            self.get_number_validated_linked_psm() == 0
        ):  # If there is no validated PSM do not return the score
            return 0

        if self.get_number_validated_linked_psm() < 3:  # Not enough PSM return worst score
            return 0

        if self.get_elution_profile() == None:
            return 0  # If no fit is found return worst score

        if self.get_elution_profile().is_parameters_fitted() == False:
            return 0  # If no fit is found return worst score

        psms_rt = [psm.spectrum.get_rt() for psm in self.get_validated_linked_psm()]
        range_rt = self.get_elution_profile().get_bounds_area()
        n_quantiles = len(psms_rt) - 1

        if range_rt[0] == range_rt[1]:
            return 0

        intervals_rt = np.linspace(range_rt[0], range_rt[1], n_quantiles + 1)

        intervals_count = [0] * n_quantiles

        for i in range(len(intervals_rt) - 1):

            n_psm_in_interval = len(
                [rt for rt in psms_rt if intervals_rt[i] <= rt and rt < intervals_rt[i + 1]]
            )

            intervals_count[i] = n_psm_in_interval

        return len([i for i in intervals_count if i > 0]) / n_quantiles

    def get_min_max_rt_range_shift(self, side):
        if (
            self.get_number_validated_linked_psm() == 0
        ):  # If there is no validated PSM do not return the score
            return [0, 0]

        if self.get_elution_profile() == None:
            return [0, 0]  # If no fit is found

        if self.get_elution_profile().is_parameters_fitted() == False:
            return [0, 0]  # If no fit is found e

        psms_rt = [psm.spectrum.get_rt() for psm in self.get_validated_linked_psm()]
        range_rt = self.get_elution_profile().get_bounds_area()

        if side == "min":
            return min(psms_rt), range_rt[0]
        if side == "max":
            return max(psms_rt), range_rt[1]

    def get_ratio_left_right(self):
        """get the ratio of psm before max elution peak and after max elution peak RT"""
        left_sum = 0
        right_sum = 0

        if self.get_elution_profile() != None:
            EP = self.get_elution_profile()
            if EP.is_parameters_fitted():
                EP_peak_rt = EP.get_x_at_max_y()

                print("Peak retention time ", EP_peak_rt)

                for psm in self.get_validated_linked_psm():

                    if psm.spectrum.get_rt() <= EP_peak_rt:
                        left_sum += psm.get_prec_intens_ratio()
                    if psm.spectrum.get_rt() >= EP_peak_rt:
                        right_sum += psm.get_prec_intens_ratio()
        else:
            return 0

        if sum([left_sum, right_sum]) == 0:
            return 0

        return min([left_sum, right_sum]) / (sum([left_sum, right_sum]) / 2)

    def get_boundaries_area_ratio(self):
        """Area under EP for interval min_boud_rt max_bound rt vs total AUC"""

        if self.get_elution_profile() != None:
            EP = self.get_elution_profile()
            if EP.is_parameters_fitted():

                auc_inter = EP.get_auc(self.min_bound_rt, self.max_bound_rt)
                auc_tot = EP.get_auc(
                    self.min_bound_rt - 500, self.max_bound_rt + 500
                )  # TODO hardcoded, change this

                if auc_tot != 0:
                    return auc_inter / auc_tot
                else:
                    return 0

        else:
            return 0

    # Setters

    def set_color(self, colorInt):
        self.color = colorInt
        return self

    def link_psm(self, psm):
        self.linkedPsm.append(psm)

    def compute_theoretical_fragments(self, ionTypes):
        """Returns and set a list of m/z of fragment ions  and informations on the type/position of each fragments for a given peptidoform/proteoform"""

        # store the fragment types looked for:
        self.fragment_types = ionTypes

        ion_formulas = constant.ion_formulas
        sequence = self.peptideSequence
        modifications = self.modificationDict

        frag_masses = {}

        # fragments masses:
        for ion_type in ionTypes:
            frag_masses_iontype = {}

            if "-" in ion_type:  # Internal Fragment
                # sum of all modification masses present in the internal fragment
                sum_mods = lambda modifications, i, j: sum(
                    [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"] <= j]
                )  # sum mods delta for internal fragm ions
                # get all sub string of the peptide sequence
                sub_sequences = [
                    (
                        sequence[i - 1 : j],
                        i,
                        j,
                        ion_type,
                        [mod["location"] for mod in modifications if i <= mod["location"] <= j],
                    )
                    for i in range(1, len(sequence))
                    for j in range(i + 1, len(sequence))
                ]
                # compute internal frag masses
                frag_masses_iontype.update(
                    {
                        ",".join(str(s) for s in seq[1:4]): round(
                            mass.fast_mass(
                                sequence=seq[0],
                                ion_type=ion_type,
                                ion_comp=ion_formulas[ion_type],
                            )
                            + sum_mods(modifications, seq[1], seq[2]),
                            4,
                        )
                        for seq in sub_sequences
                    }
                )

            else:  # Terminal Fragment
                if any(i in ion_type for i in ["a", "b", "c"]):  # Nterm
                    sum_mods = lambda modifications, i, j: sum(
                        [mod["monoisotopicMassDelta"] for mod in modifications if mod["location"] <= j]
                    )
                    sub_sequences = [
                        (
                            sequence[:j],
                            1,
                            j,
                            ion_type,
                            [mod["location"] for mod in modifications if mod["location"] <= j],
                        )
                        for j in range(2, len(sequence))
                    ]
                    frag_masses_iontype.update(
                        {
                            ",".join(str(s) for s in seq[1:4]): round(
                                mass.fast_mass(
                                    sequence=seq[0],
                                    ion_type=ion_type,
                                    ion_comp=ion_formulas[ion_type],
                                )
                                + sum_mods(modifications, seq[1], seq[2]),
                                4,
                            )
                            for seq in sub_sequences
                        }
                    )

                else:  # Cterm
                    sum_mods = lambda modifications, i, j: sum(
                        [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"]]
                    )
                    sub_sequences = [
                        (
                            sequence[i - 1 :],
                            i,
                            len(sequence),
                            ion_type,
                            [mod["location"] for mod in modifications if i <= mod["location"]],
                        )
                        for i in range(1, len(sequence) + 1)
                    ]
                    frag_masses_iontype.update(
                        {
                            ",".join(str(s) for s in seq[1:4]): round(
                                mass.fast_mass(
                                    sequence=seq[0],
                                    ion_type=ion_type,
                                    ion_comp=ion_formulas[ion_type],
                                )
                                + sum_mods(modifications, seq[1], seq[2]),
                                4,
                            )
                            for seq in sub_sequences
                        }
                    )

            frag_masses[ion_type] = frag_masses_iontype

        if self.theoFrag == None:
            self.theoFrag = frag_masses
        else:
            warning(
                "Theoretical fragments already set for proteoform: " + self.modificationBrno + "OVERWRITING !"
            )
            self.theoFrag = frag_masses

        return frag_masses

    def model_elution_profile(self, elution_profile_score_threshold):
        """instanciate an envelope object by providing the list of psm associated to that proteoform"""
        # reset elution profile:
        self.envelope = None

        if len(self.get_validated_linked_psm()) > 5:
            env = ElutionProfile()
            env.model_elution_profile(self.get_validated_linked_psm(), elution_profile_score_threshold)

            if env.score_threshold == None:  # no threshold
                self.envelope = env
            else:  # threshold (min boundary !!)
                if env.score_fitted > elution_profile_score_threshold:
                    self.envelope = env
                else:
                    # self.envelope = env # !!!
                    # print("Could not fit curve to chromatogram")
                    pass
        else:
            # print("Not enough datapoints to compute envelope")
            pass

    def update_proteoform_total_intens(self, method="precursor"):
        """Return the sum of intensities of psm of that proteoform method = precursor  or annotated (correspond to the intensity value used)"""

        self.totalIntens = 0

        if method == "AUC":
            if self.get_elution_profile() != None:
                # print(self.get_elution_profile().get_auc())
                self.totalIntens += self.get_elution_profile().get_auc()
            return None

        for psm in self.get_validated_linked_psm():
            if method == "precursor":
                self.totalIntens += psm.get_prec_intens_ratio()
            if method == "annotated":
                self.totalIntens += psm.getAnnotMsmsIntensRatio()

    def update_proteoform_psm_validation(self):
        """ """

        if self.get_elution_profile() == None:  # if no envelopes are found for that proteoform
            for psm in self.linkedPsm:
                psm.is_validated = False
                psm.ratio = 0.0
        else:
            for psm in self.get_elution_profile().psms_outliers:
                # print("removing aberant psm")
                psm.is_validated = False
                psm.ratio = 0.0
