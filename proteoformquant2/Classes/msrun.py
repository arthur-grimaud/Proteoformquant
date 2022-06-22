### Import ###
from distutils.log import warn
from pickle import TRUE
from pyteomics import mzid
from pyteomics import mgf
from pyteomics import mzml
from chart_studio.plotly import plot, iplot
import plotly.express as px
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from statistics import mean
from statistics import median
from progress.bar import Bar
from alive_progress import alive_bar
import pandas as pd
import networkx as nx
from itertools import combinations
import matplotlib.pyplot as plt
import multiprocessing
from pymoo.core.problem import starmap_parallelized_eval
from warnings import warn
import multiprocessing as mp

# Custom classes
from Classes.spectrum import Spectrum
from Classes.proteoform import Proteoform
from Classes.proteoform0 import Proteoform0
from Classes.problem_wrapper import ProblemWrapper
from pymoo.core.evaluator import Evaluator
from pymoo.core.population import Population
from pymoo.factory import get_problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.visualization.scatter import Scatter

# GA
import pymoo
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.algorithms.soo.nonconvex.es import ES
from pymoo.factory import get_crossover, get_mutation, get_sampling
from pymoo.optimize import minimize
from pymoo.core.problem import Problem

# Helpers
from Utils import misc
from pprint import pprint

# For visualization (TEMPORARY)
# from report import plot_elution_profiles
import plotly.graph_objects as go


class Msrun:
    """A class used to reprensent one MSrun with identification

    The Msrun stores all data corresponding one ms run from an identification file (in .mzid format)
    and a spectrum file (in .mgf or .mzml).
    This class contains all the methods called from the mains script
    """

    def __init__(self, run_id: str = "Default run ID", dbse: str = "comet"):

        self.run_id = run_id
        self.dbse = dbse
        self.ident_fn: str = "Not Specified"
        self.spectra_fn: str = "Not Specified"

        self.spectra: dict(Spectrum) = {}
        self.proteoforms: dict(Proteoform) = {}
        self.proteoform0 = Proteoform0()  # Obsolete ?

        self.proteoform_isobaric_group = []

        self.mod_mass_to_name = {}  # store mass to name to limit queries to unimod DB

        ### PARAMETERS (should ultimately be in a separated file) ###
        self.prec_mz_tol = 1.05  # in Da
        self.frag_mz_tol = 0.02  # in Da

        # Filtering and thresholds:
        self.max_rank = 10  # Maximum rank to be considered (1-based)
        self.fdr_threshold = 0.01
        self.intensity_threshold = 100000
        self.elution_profile_score_threshold = None

        # Proteform groups:
        self.min_connect = 5

        # Genetic Algorithm parameters:
        self.n_individuals = 5
        self.n_generation = 5
        # multiprocessing

        # Fragments to be considered (see options in "Utils")
        self.fragments_types = [
            "c",
            "zdot",
            "z+1",
            "z+2",
        ]

        # self.fragments_types = ["a","x","b","y","c","zdot","z+1","z+2","c-zdot", "c-z+1", "cdot-zdot", "c-z+1", "a-x", "n-n", "b-y"]

    def __getstate__(self):

        # this method is called when you are
        # going to pickle the class, to know what to pickle
        state = self.__dict__.copy()

        # don't pickle the parameter fun. otherwise will raise
        # AttributeError: Can't pickle local object 'Process.__init__.<locals>.<lambda>'
        # del state["pool"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    # ---------------------------------- getters --------------------------------- #

    def get_rt_range(self):
        """Returns the and max retention time for all spectra in self.spectra"""
        mn = min([spectrum.get_rt() for spectrum in self.spectra.values()])
        mx = max([spectrum.get_rt() for spectrum in self.spectra.values()])
        return (mn, mx)

    def get_mz_range(self):
        """Returns the min and max mass to charge ratio for all spectra in self.spectra"""
        mn = min([spectrum.getPrecMz() for spectrum in self.spectra.values()])
        mx = max([spectrum.getPrecMz() for spectrum in self.spectra.values()])
        return (mn, mx)

    def get_dataset_metrics(self):
        """Returns a set of metrics about the run as a dictionnary"""

        return {
            "TotalSpectra": len(self.spectra),
            "AssignedSpectra(method1)": len(
                [spectrum for spectrum in self.spectra.values() if spectrum.get_number_validated_psm() > 0]
            ),
            "AssignedSpectra(method2)": 0,
            "UnassignedSpectra": len([spectrum for spectrum in self.proteoform0.linkedSpectra]),
            "ChimericSpectra": len(
                [spectrum for spectrum in self.spectra.values() if spectrum.get_number_validated_psm() > 1]
            ),
            "TotalProteoforms": len(self.proteoforms) - 1,
            # "AvgEnvelopeScore": mean([proteoform.get_elution_profile().score_fitted for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            # "MinEnvelopeScore": min([proteoform.get_elution_profile().score_fitted for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            # "MedianEnvelopeS":   median([proteoform.get_elution_profile().param_fitted[1] for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            # "MedianEnvelopeA":   median([proteoform.get_elution_profile().param_fitted[2] for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            # "MedianEnvelopeK":   median([proteoform.get_elution_profile().param_fitted[3] for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ])
        }

    # ----------------------------------- main ----------------------------------- #

    # ------------------------- DATA PREPARATION METHODS ------------------------ #

    # Here wraping functions (read identification data and read spectra data can be added to wrap multiple file format)

    def read_mzid(self, ident_fn):
        """Read a spectra identification file in .mzIdenMl and store data
        Path of the identification is stored in self.inputFn
        Spectrum objects with identification data are stored in self.spectra
        Psm objects are stored within corresponding spectra objects
        """
        self.ident_fn = ident_fn  # Store File that has been read
        mzidObj = mzid.read(ident_fn)  # Create a pyteomics' mzid iterator

        print("---Loading identification data from: {0}---".format(self.ident_fn))
        with alive_bar(0) as bar:

            for (
                identMzid
            ) in mzidObj:  # Iterate over identificationresults in mzid and create Spectrum objects for each
                self.spectra[identMzid["spectrumID"]] = Spectrum(
                    spectrumID=identMzid["spectrumID"], identMzid=identMzid, max_rank=self.max_rank
                )
                bar()
        pass

    def read_mgf(self, spectra_fn):
        """Add informations from mgf file to spectrum objects in self.spectra"""

        self.spectra_fn = spectra_fn  # Store File name that has been read
        mgf_obj = mgf.read(spectra_fn)

        # print(mgf_obj)

        # pprint(mgf_obj.__dict__)

        print("---Loading spectrum data from: {0}---".format(self.spectra_fn))
        with alive_bar(0) as bar:

            # Spectrum title are not stored the same way depending on the mgf and mzid file, therefore :
            for specID in self.spectra:
                try:
                    specMgf = mgf_obj.get_spectrum(self.spectra[specID].spectrum_title)
                except (KeyError):
                    try:
                        specMgf = mgf_obj.get_spectrum(self.spectra[specID].spectrum_title_alt)
                    except (KeyError):
                        try:
                            specMgf = mgf_obj.get_spectrum(self.spectra[specID].spectrum_title_alt_alt)
                        except (KeyError):
                            print(
                                f" \n Spectrum with ID: {self.spectra[specID].spectrum_title} \n or {self.spectra[specID].spectrum_title_alt} \n or {self.spectra[specID].spectrum_title_alt_alt}  not found"
                            )

                self.spectra[specID].set_spec_data_mgf(specMgf)
                bar()
        pass

    def read_mzml(self, spectra_fn):
        """Add informations from mzml file to spectrum objects in self.spectra"""

        self.spectra_fn = spectra_fn  # Store File name that has been read
        mzml_obj = mzml.read(spectra_fn)

        print("---Loading spectrum data from: {0}---".format(self.spectra_fn))
        with alive_bar(0) as bar:

            for specID in self.spectra:  # Take into account how DBSEs store spectra ids

                if self.dbse == "mascot":
                    index = str(int(specID.split("=")[1]) + 1)
                if self.dbse == "comet":
                    index = specID.split("=")[1]

                try:
                    spec_mzml = mzml_obj.get_by_id(index, id_key="index")  # need to be splited
                except (KeyError):
                    print("Key error spectrum not found")

                self.spectra[specID].set_spec_data_mzml(spec_mzml)
                bar()
        pass

    def fdr_filtering(self, decoy_tag, score_name):

        """
        Compute FDR and remove psms below score threshold at FDR = self.fdr_threshold

        decoy_tag: String in protein_description indicative of decoy sequence
        score_name: name of the score reported by the database search engine (IDEA: define by DBSE used)
        (comet: 'Comet:xcorr', msamanda:'Amanda:AmandaScore')
        """

        # generate list of all psms (IDEA: could be method)
        psms_obj = []
        psms_score = []
        psms_isdecoy = []

        print("---FDR filtering---")
        with alive_bar(0) as bar:

            for spectrum in self.spectra.values():
                for psm in spectrum.get_psms():
                    psms_obj.append(psm)
                    bar()
                    try:
                        psms_score.append(psm.__dict__[score_name])
                    except KeyError:
                        print(
                            "The score '{0}' was not found and fdr fitering has been aborted".format(
                                score_name
                            )
                        )
                        return 0

                    isdecoy = False
                    # pprint(psm.__dict__)
                    # print(psm.get_accessions())
                    for protein_ref in psm.get_accessions():  # TODO not generic for all DBSE
                        if decoy_tag in protein_ref:
                            isdecoy = True

                    psms_isdecoy.append(isdecoy)

            psms_score, psms_isdecoy, psms_obj = map(
                list, zip(*sorted(zip(psms_score, psms_isdecoy, psms_obj), key=lambda x: x[0], reverse=True))
            )

            decoy_hits = 0
            target_hits = 0
            score_fdr = 0

            for i in range(len(psms_isdecoy)):
                bar()
                if psms_isdecoy[i]:
                    decoy_hits += 1
                    psms_obj[i].exclude()  # remove decoy psms for further analysis
                else:
                    target_hits += 1

                # print(decoy_hits / target_hits)

                if decoy_hits / target_hits < self.fdr_threshold:
                    score_fdr = psms_score[i]

            idx_remove = [idx for idx, element in enumerate(psms_score) if element < score_fdr]

            psms_remove = [psms_obj[index] for index in idx_remove]

            for psm in psms_remove:
                psm.exclude()

            if decoy_hits == 0:
                print("No decoys hit were found")
            if score_fdr == 0:
                print("Warning: FDR threshold score calculated as 0")

        print(
            f"{decoy_hits} decoys and {target_hits} targets hits found, FDR {self.fdr_threshold} for score: {score_fdr}"
        )

        # print(idx_remove)

    def add_proteoforms(self):
        """From spectrum objects in self.spectra instanciates proteoforms object and stores them self.proteoforms"""
        i = 0
        print("---Instanciating proteoform objects---")
        with alive_bar(0) as bar:

            for specID in self.spectra:
                i += 1

                for psm in self.spectra[specID].get_psms():

                    proforma = psm.get_modification_proforma(self.mod_mass_to_name)
                    brno = psm.get_modification_brno()
                    seq = psm.getPeptideSequence()
                    brnoSeq = brno + "-" + seq  # TODO use proforma

                    # Create list with associated protein description (try except as variable name could be different depending on the DBSE used)
                    try:
                        protein_ids = [ref["DBSequence_Ref"] for ref in psm.PeptideEvidenceRef]
                    except KeyError:
                        protein_ids = [ref["protein description"] for ref in psm.PeptideEvidenceRef]

                    if (
                        proforma not in self.proteoforms.keys()
                    ):  # if a proteoform is new create a instance of Proteoform for this proteoform
                        bar()
                        self.proteoforms[proforma] = Proteoform(
                            peptideSequence=seq,
                            modificationBrno=brno,
                            modificationProforma=proforma,
                            modificationDict=psm.Modification,
                            protein_ids=protein_ids,
                        ).set_color(i)

                    self.proteoforms[proforma].link_psm(psm)  # add link to Psms in Proteoform
                    psm.setProteoform(self.proteoforms[proforma])  # add link to Proteoform in Psm

    def match_fragments(self, internal=False):
        """
        For each proteoform object generates the set of theoretical fragments
        If mgf and identification data are provided in a spectrum object,
        get the and store the annotated fragments for each PSM"""

        print("---Generating theoretical fragments---")
        with alive_bar(0) as bar:
            for proteoID in self.proteoforms:
                bar()
                self.proteoforms[proteoID].compute_theoretical_fragments(self.fragments_types)

        print("---Matching fragments---")
        with alive_bar(0) as bar:
            for spectrumID in self.spectra:
                self.spectra[spectrumID].annotateFragPsm(frag_mz_tol=self.frag_mz_tol)
                self.spectra[spectrumID].setSumIntensAnnotFrag()
                bar()

        pass

    def scale_precursor_intensities(self, target_range=[0, 100]):
        print("---Scaling precursor intensities---")
        with alive_bar(0) as bar:
            spectra_prec_intens = [spectrum.getPrecIntens() for spectrum in self.spectra.values()]
            max_prec_intens = max(spectra_prec_intens)
            min_prec_intens = min(spectra_prec_intens)

            factor = (target_range[1] - target_range[0]) / (max_prec_intens - min_prec_intens)

            for spectrum in self.spectra.values():
                spectrum.precIntens = factor * (spectrum.precIntens + target_range[0])

                bar()

    # ------------------------- METHODS FOR QUANTIFICATION ------------------------ #

    def update_proteoforms_elution_profile(self):
        """For each proteoforms in self. proteoforms model the elution profile"""
        with Bar("Updating proteoforms envelopes", max=1) as bar:

            # TODO add a function thjat set the bounds based on the entire set of envelopes
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].model_elution_profile(self.elution_profile_score_threshold)
                bar.next()
        pass

    def update_proteoform_intens(self):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        with Bar("Updating Proteoforms Total Intens", max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].update_proteoform_total_intens()
                bar.next()
        pass

    def update_psm_validation(self):
        """Update psm.is_validated"""
        with Bar("Updating PSM Validation", max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].update_proteoform_psm_validation()
        pass

    def update_unassigned_spectra(self):
        """add spectra without any validated psm to "proteoform0" in self.proteoforms"""
        for spectrumID in self.spectra:
            if self.spectra[spectrumID].get_number_validated_psm() == 0:
                self.proteoform0.linkSpectrum(self.spectra[spectrumID])

    def validate_all_psms(self):
        for spectrum in self.spectra.values():
            for psm in spectrum.get_psms():
                psm.is_validated = True

    def filter_proteform_low_count(self, min_n_psm):
        n_exclude = 0
        n_proteo = 0
        for proteoform in self.proteoforms.values():
            if len(proteoform.get_linked_psm()) < min_n_psm:
                n_proteo += 1
                for psm in proteoform.get_linked_psm():
                    psm.exclude()
                    n_exclude += 1
        print(n_proteo, " proteoforms and ", n_exclude, " PSMs have been excluded based on count")

    def update_chimeric_spectra(self, max_rank):

        """For every psm of rank < rank try to find a matching envelope, and assign to that proteoform it if it is the case"""
        for spectrum in self.spectra.values():
            spectrumMz = spectrum.getPrecMz()
            spectrumRt = spectrum.get_rt()

            for rank in range(1, max_rank):
                if len(spectrum.get_psms()) >= rank:
                    psm = spectrum.get_psms()[rank - 1]
                    psmProforma = psm.get_modification_proforma()

                    altProteo = self.proteoforms[
                        psmProforma
                    ]  # TODO might need a try except for proteo only in second rank

                    if altProteo.get_elution_profile() != None:
                        if altProteo.get_elution_profile().get_y(spectrumRt) > self.intensity_threshold:

                            # print("spectrum at RT {0} psm: {1} matches proteoform {2}".format(spectrumRt, psmProforma, altProteo.get_modification_brno()))
                            psm.is_validated = True
                            # altProteo.link_psm(psm)

            spectrum.update_psms_ratio()

    def update_psms_ratio_subset(self, spectra_subset):
        for spectrum in spectra_subset:
            spectrum.update_psms_ratio(0)

    # ----------------------- TESTING GROUP QUANTIFICATION ----------------------- #

    def set_proteoform_isobaric_groups(self):
        """From self.proteoforms define self.proteoform_isobaric_group where isobaric proteoform are grouped"""
        # print("In: set proteoform isobaric group")

        print("---Subsetting proteoform---")
        with alive_bar(0) as bar:
            all_proteoform_pairs = []  # pairs of proteoform founds together in the psms of a spectra
            all_proteoform_pairs_count = []

            for spectrum in self.spectra.values():
                bar()

                # print(spectrum.get_id())
                # print([psm.get_modification_brno() for psm in spectrum.get_psms()])

                proforma_psms = [psm.proteoform.get_modification_proforma() for psm in spectrum.get_psms()]
                proforma_psms = np.unique(proforma_psms)
                proforma_combinations = [
                    sorted((a, b)) for idx, a in enumerate(proforma_psms) for b in proforma_psms[idx:]
                ]

                for c in proforma_combinations:
                    if c not in all_proteoform_pairs:
                        all_proteoform_pairs.append(c)
                        all_proteoform_pairs_count.append(1)
                    else:
                        all_proteoform_pairs_count[all_proteoform_pairs.index(c)] += 1

            # keep pairs found at least self.min_connect times
            all_proteoform_pairs_filt = [
                all_proteoform_pairs[i]
                for i in range(len(all_proteoform_pairs))
                if all_proteoform_pairs_count[i] >= self.min_connect
            ]
            # print(all_proteoform_pairs_filt)

            # find groups (connected graphs) in the network defined from the edgelist "all-proteoform-pairs"
            G = nx.from_edgelist(all_proteoform_pairs_filt)

            self.isobaric_proteform_graph = G
            # # plor graph
            G.remove_edges_from(nx.selfloop_edges(G))
            nx.draw(G, node_size=10)  # , with_labels=True)
            plt.show()

            l = list(nx.connected_components(G))

            # print("connected : ", l)

            # for g in l:
            #     theo_masses = [self.proteoforms[p].getTheoPrecMz() for p in g]
            #     print("all theo mz")
            #     print(theo_masses)
            #     print("max mz")
            #     print(max(theo_masses))
            #     print("min mz")
            #     print(min(theo_masses))

            self.proteoform_isobaric_group = l

    # ------------------------ Proteoform subset operation ----------------------- #

    def update_psms_validation_proteoform_subset(self, spectra_subset, proteoforms_proforma):
        for spectrum in spectra_subset:
            for psm in spectrum.get_psms():
                if psm.proteoform.get_modification_proforma() in proteoforms_proforma:
                    psm.is_validated = True
                else:
                    psm.is_validated = False

    def update_psms_validation_subset(self, psm_list, validation_list):
        """Given a psm object list and a binary list of psm to validate validates the psms:"""
        for i in range(len(validation_list)):
            if validation_list[i] == 1:
                psm_list[i].is_validated = True
            else:
                psm_list[i].is_validated = False

    def update_psms_ratio_subset(self, spectra_subset):
        for spectrum in spectra_subset:
            spectrum.update_psms_ratio(verbose=False)

    def update_proteoforms_elution_profile_subset(self, proteoform_subset):
        """For each proteoforms in self. proteoforms model the elution profile"""

        # IDEA add a function that set the bounds based on the entire set of envelopes

        multi = False  # run with mp

        processes = []

        for proteoform in proteoform_subset:

            if multi:
                p = mp.Process(
                    target=proteoform.model_elution_profile, args=(self.elution_profile_score_threshold,)
                )
                print(p)
                processes.append(p)
            else:
                proteoform.model_elution_profile(self.elution_profile_score_threshold)

        if multi:
            print(processes)
            [x.start() for x in processes]

    # ---------------------------------- SCORING --------------------------------- #

    # -------------------------- SCORE: signal explained ------------------------- #

    def score_signal_explained(self, proteoform_subset, spectra_subset, boundaries):
        # Get individual scores for the subset
        r_explained = self.get_intensity_explained_ratio(spectra_subset)
        # Group scores and define weights (TODO hard-coded weights)
        scores = [r_explained]
        weights = [1]
        # Compute weighted average
        overall_score = sum([scores[i] * weights[i] for i in range(len(scores))]) / sum(weights)

        return overall_score

    def get_intensity_explained_ratio(self, spectra_subset):
        i_tot = 0
        i_explained = 0

        for spectrum in spectra_subset:
            i_tot += spectrum.getSumIntensFrag()
            # print(spectrum.getSumIntensFrag(), " ... ", spectrum.getSumIntensAnnotFrag())
            i_explained += spectrum.getSumIntensAnnotFrag()

        return 1 - (i_explained / i_tot)

    # ---------------------- SCORE: Elution profile quality ---------------------- #

    def score_elution_profile_quality(self, proteoform_subset, spectra_subset, boundaries, verbose=False):

        # Get individual scores for the subset
        u_correlation = self.get_correlation_score(proteoform_subset)
        # u_gap_spectra_valid = self.get_gap_in_validated_spectra(proteoform_subset, spectra_subset, boundaries)
        u_coverage = self.get_coverage_score(proteoform_subset)
        # u_balance = self.get_balance_evidence_score(proteoform_subset)

        # Group scores and define weights (TODO hard-coded weights)
        scores = [u_correlation, u_coverage]
        weights = [1, 1]

        # prints:
        if verbose:
            # print("u_correlation, u_gap_spectra_valid, u_coverage")
            print(scores)
        # Compute weighted average
        overall_score = sum([scores[i] * weights[i] for i in range(len(scores))]) / sum(weights)

        return overall_score

    def get_correlation_score(self, proteoform_subset):
        scores_correlation = [
            proteoform.get_fit_score()
            for proteoform in proteoform_subset
            if proteoform.get_fit_score() != None
        ]

        if len(scores_correlation) != 0:
            return 1 - mean(scores_correlation)
        else:
            return 1

    def get_gap_in_validated_spectra(self, proteoform_subset, spectra_subset, boundaries):
        ratios_missed = []  # ratio of "missed psm" in the range for each proteoforms
        p = 0
        for b in range(0, len(boundaries), 2):

            proteoform = proteoform_subset[p]
            p += 1
            # print(proteoform)
            v_proteo, v_subset = 0, 0

            # print("n linked psms: ", len(proteoform.get_linked_psm()))
            for psm in proteoform.get_linked_psm():
                if psm.spectrum.get_rt() > min(
                    boundaries[b], boundaries[b + 1]
                ) and psm.spectrum.get_rt() < max(boundaries[b], boundaries[b + 1]):
                    v_proteo += 1
            for spectrum in spectra_subset:
                if spectrum.get_rt() > min(boundaries[b], boundaries[b + 1]) and spectrum.get_rt() < max(
                    boundaries[b], boundaries[b + 1]
                ):
                    v_subset += 1

            if v_subset != 0:
                ratios_missed.append(v_proteo / v_subset)
            else:
                pass  # if no psm do not include
                # ratios_missed.append(1)  # if no validated psm add perfect score

        if len(ratios_missed) != 0:
            return 1 - mean(ratios_missed)
        else:
            return 1

    def get_balance_evidence_score(self, proteoform_subset):

        ratios_l_r = [
            proteoform.get_ratio_left_right()
            for proteoform in proteoform_subset
            if proteoform.get_fit_score() != None
        ]

        return 1 - mean(ratios_l_r)

    def get_area_under_bound_score(self, proteoform_subset):

        ratios_area = [
            proteoform.get_boundaries_area_ratio()
            for proteoform in proteoform_subset
            if proteoform.get_fit_score() != None
        ]

        return 1 - mean(ratios_area)

    def get_coverage_score(self, proteoform_subset):

        coverage = [proteoform.get_coverage_2() for proteoform in proteoform_subset]
        # print(coverage)
        coverage = [i for i in coverage if i]  # remove None values
        if len(coverage) == 0:
            return 1
        return 1 - mean(coverage)

    # ---------------------- SCORE: Chimeric Spectra Quality --------------------- #

    def score_chimeric_spectra_quality(self, proteoform_subset, spectra_subset, boundaries):
        # Get individual scores for the subset
        u_residuals = self.get_residuals_subset(spectra_subset)
        # u_miss_r1 = self.get_ratio_missed_rank_1(spectra_subset)
        u_gap_psm_valid = self.get_gap_in_rank_psms(spectra_subset)
        # u_psm_ratio_0 = self.get_ratio_validated_0(spectra_subset)

        # Group scores and define weights (TODO hard-coded weights)
        scores = [u_residuals, u_gap_psm_valid]
        weights = [1, 1]
        # Compute weighted average
        overall_score = sum([scores[i] * weights[i] for i in range(len(scores))]) / sum(weights)

        return overall_score

    def get_residuals_subset(self, spectra_subset):
        # TODO lower penality for low residuals

        residuals = [
            spectra.get_residuals(threshold=0.0)
            for spectra in spectra_subset
            if spectra.get_residuals(threshold=0.0) != None
        ]

        if len(residuals) == 0:
            return 1

        # print(residuals)

        return mean(residuals)

    def get_ratio_missed_rank_1(self, spectra_subset):
        s_unvalidated = 0
        s_tot = 0

        for spectrum in spectra_subset:
            if spectrum.get_psms()[0].is_validated == False:
                s_unvalidated += 1
            s_tot += 1

        return s_unvalidated / s_tot

    def get_gap_in_rank_psms(self, spectra_subset):
        p_0 = 0
        p_tot = 0

        for spectrum in spectra_subset:

            val_list = [psm.is_validated for psm in spectrum.get_psms()]
            val_index = [i for i, x in enumerate(val_list) if x]

            if len(val_index) != 0:
                last_val_index = val_index[-1]
                for p in range(int(last_val_index)):
                    if val_list[p] == False:
                        p_0 += 1
                    p_tot += 1

        if p_tot == 0:
            return 1
        return p_0 / p_tot

    def get_ratio_validated_0(self, spectra_subset):
        p_0 = 0
        p_tot = 0

        for spectrum in spectra_subset:
            for psm in spectrum.get_psms():
                if psm.is_validated == True and psm.ratio == 0:

                    p_0 += 1

                p_tot += 1

        return p_0 / p_tot

    # Others

    def get_rank_score(self, proteoform_subset):
        rank_scores = []
        for proteoform in proteoform_subset:

            if proteoform.get_number_validated_linked_psm() != 0:
                rank_scores.append(
                    proteoform.get_weighted_number_linked_validated_psm(self.max_rank)
                    / (self.max_rank * proteoform.get_number_validated_linked_psm())
                )
            else:
                rank_scores.append(1)

        # print(rank_scores)
        return mean(rank_scores)

    # ------------------------------- OPTIMIZATION ------------------------------- #

    def update_proteoform_subset_validation(self, proteoform_subset, boundaries):

        p = 0
        for b in range(len(boundaries)):

            proteoform = proteoform_subset[p]
            p += 1
            v, uv = 0, 0

            # min_idx = min(boundaries[b], boundaries[b + 1])
            # max_idx = max(boundaries[b], boundaries[b + 1])
            # print(boundaries[b])
            proteoform.min_bound_rt = boundaries[b][0]
            proteoform.max_bound_rt = boundaries[b][1]

            for psm in proteoform.get_linked_psm():
                if psm.spectrum.get_rt() >= boundaries[b][0] and psm.spectrum.get_rt() <= boundaries[b][1]:
                    psm.is_validated = True
                    v += 1
                else:
                    psm.is_validated = False
                    uv += 1

    def test_proteoform_subsets_scoring(self):
        group = [
            "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",
            "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",  ### K27me3
            "ARTKQTARKSTGGKAPRKQLATK[Acetyl]AARKSAPATGGVKKPHRYRPGTVALRE",
            "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Trimethyl]KPHRYRPGTVALRE",  ###K36me3
            "ARTKQTARKSTGGKAPRK[Acetyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",
            "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",  ###K36ac
            "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",  ###K27ac
            "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Acetyl]PHRYRPGTVALRE",
            "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Trimethyl]PHRYRPGTVALRE",
        ]

        # Get list of proteoform and spectra objects in the group
        self.proteoform_subset = []
        self.spectra_subset = []
        self.variables = []

        for proforma in group:
            proteoform = self.proteoforms[proforma]
            self.proteoform_subset.append(proteoform)

            # get starting values (center of mass based on psm rank and precursor intensity)
            rt_center = proteoform.get_rt_center()
            self.variables.append(rt_center)
            self.variables.append(rt_center)

            for psm in proteoform.get_linked_psm():
                if psm.spectrum not in self.spectra_subset:
                    self.spectra_subset.append(psm.spectrum)

        bounds = [
            0,
            0,
            4763,
            5000,
            0,
            0,
            4800,
            4900,
            0,
            0,
            4400,
            4600,
            4530,
            4750,
            0,
            0,
            0,
            0,
        ]

        print(group)
        print(bounds)

        self.update_proteoform_subset_validation(self.proteoform_subset, bounds)
        self.update_psms_ratio_subset(self.spectra_subset)
        self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

        print(
            self.score_signal_explained(self.proteoform_subset, self.spectra_subset, bounds),
            self.score_elution_profile_quality(self.proteoform_subset, self.spectra_subset, bounds),
            self.score_chimeric_spectra_quality(self.proteoform_subset, self.spectra_subset, bounds),
        )

        print(
            mean(
                [
                    self.score_signal_explained(self.proteoform_subset, self.spectra_subset, bounds),
                    self.score_elution_profile_quality(self.proteoform_subset, self.spectra_subset, bounds),
                    self.score_chimeric_spectra_quality(self.proteoform_subset, self.spectra_subset, bounds),
                ]
            )
        )

    def optimize_proteoform_subsets(self):
        group_number = 0

        rd_loc = 5
        rd_scale = 4

        # Process proteoform groups separately
        for group in self.proteoform_isobaric_group:

            group_number += 1

            # Variable TODO hard-coded!

            g = 0
            n_iterations = 90
            score_ep_proteos = []
            all_scores = []
            all_bounds = []

            self.proteoform_subset = []  # Proteoform objects in the subset
            self.spectra_subset = []  # Spectra objects in the subset
            self.rt_boundaries = []  # retention time validation range for the proteoforms in the subset

            # Define proteoform subset and spectra subset
            for proforma in group:
                proteoform = self.proteoforms[proforma]
                self.proteoform_subset.append(proteoform)
                for psm in proteoform.get_linked_psm():
                    if psm.spectrum not in self.spectra_subset:
                        self.spectra_subset.append(psm.spectrum)

            # TODO could be improved (generalize)
            # Sort proteoform in subset (based on psms rank count)
            n_rank_1 = [proteoform.get_number_linked_psm_Rx(rank=1) for proteoform in self.proteoform_subset]
            n_rank_2 = [proteoform.get_number_linked_psm_Rx(rank=2) for proteoform in self.proteoform_subset]
            n_rank_3 = [proteoform.get_number_linked_psm_Rx(rank=3) for proteoform in self.proteoform_subset]
            n_rank_4 = [proteoform.get_number_linked_psm_Rx(rank=4) for proteoform in self.proteoform_subset]
            n_rank_5 = [proteoform.get_number_linked_psm_Rx(rank=5) for proteoform in self.proteoform_subset]

            zipped_rank_proteo = zip(
                n_rank_1, n_rank_2, n_rank_3, n_rank_4, n_rank_5, group, self.proteoform_subset
            )
            zipped_proteo = sorted(zipped_rank_proteo, reverse=True)
            self.proteoform_subset = [list(tuple)[-1] for tuple in zipped_proteo]

            # Define rt_boundaries starting values
            for proteoform in self.proteoform_subset:
                if proteoform.get_number_linked_psm_Rx(rank=1) >= 5:
                    self.rt_boundaries.append(proteoform.get_rt_range(rank=1))
                else:
                    self.rt_boundaries.append([0, 0])

            # Sorted retention times of spectra in the subset
            all_rts = sorted([spectrum.get_rt() for spectrum in self.spectra_subset])

            # Initial validation
            self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
            self.update_psms_ratio_subset(self.spectra_subset)
            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

            fig = self.plot_elution_profiles(self.proteoform_subset, rt_values=all_rts, count=g)
            fig.write_image("images/fig_" + f"{group_number:03}" + "_" + f"{g:04}" + ".png")

            # Initial scores:
            for proteo in self.proteoform_subset:
                score = mean([proteo.get_coverage_2(), proteo.get_fit_score()])
                score_ep_proteos.append(score)

            # Start optimizing
            for i in range(n_iterations):
                g += 1
                all_scores.append([])
                all_bounds.append([])
                all_scores[i].append(i)

                for p in range(len(self.proteoform_subset)):
                    if self.rt_boundaries[p][0] != 0:
                        proteo_obj = self.proteoform_subset[p]
                        # print(proteoform_sorted[p])
                        # Minimun boundarie mutation:
                        min_spec_rt, min_ep_rt = proteo_obj.get_min_max_rt_range_shift(side="min")
                        if min_spec_rt < min_ep_rt:  # if spectra lower than modeled ep
                            # modifier_min = np.random.choice([+3, +1, -1, -3], p=[0.3, 0.5, 0.15, 0.05])
                            modifier_min = int(
                                np.round(np.random.normal(loc=rd_loc, scale=rd_scale, size=1)[0])
                            )
                            # print("increase lower bound")
                        else:
                            # modifier_min = np.random.choice([-3, -1, +1, +3], p=[0.3, 0.5, 0.15, 0.05])
                            modifier_min = int(
                                np.round(np.random.normal(loc=-rd_loc, scale=rd_scale, size=1)[0])
                            )
                            # print("decrease lower bound")

                        index_rt_min_start = all_rts.index(self.rt_boundaries[p][0])
                        try:
                            # print("modifier_min:", modifier_min)
                            self.rt_boundaries[p][0] = all_rts[index_rt_min_start + modifier_min]
                        except IndexError:
                            pass

                        # Maximum boundarie mutation:
                        max_spec_rt, max_ep_rt = proteo_obj.get_min_max_rt_range_shift(side="max")
                        if max_spec_rt > max_ep_rt:  # if spectra lower than modeled ep
                            # print("decrease upper bound")
                            # modifier_max = np.random.choice([-3, -1, +1, +3], p=[0.3, 0.5, 0.15, 0.05])
                            modifier_max = int(
                                np.round(np.random.normal(loc=-rd_loc, scale=rd_scale, size=1)[0])
                            )
                        else:
                            # print("increase upper bound")
                            # modifier_max = np.random.choice([+3, +1, -1, -3], p=[0.3, 0.5, 0.15, 0.05])
                            modifier_max = int(
                                np.round(np.random.normal(loc=rd_loc, scale=rd_scale, size=1)[0])
                            )

                        index_rt_max_start = all_rts.index(self.rt_boundaries[p][1])
                        try:
                            self.rt_boundaries[p][1] = all_rts[index_rt_max_start + modifier_max]
                            # print("modifier_max:", modifier_max)
                        except IndexError:
                            pass

                        # TEst score

                        self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
                        self.update_psms_ratio_subset(self.spectra_subset)
                        self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

                        new_score = mean([proteo_obj.get_coverage_2(), proteo_obj.get_fit_score()])
                        if new_score >= score_ep_proteos[p] - 0.1:  # keep if improve
                            score_ep_proteos[p] = new_score

                        else:  # undo if worse

                            self.rt_boundaries[p][0] = all_rts[index_rt_min_start]
                            self.rt_boundaries[p][1] = all_rts[index_rt_max_start]

                            self.update_proteoform_subset_validation(
                                self.proteoform_subset, self.rt_boundaries
                            )
                            self.update_psms_ratio_subset(self.spectra_subset)
                            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

                        all_scores[i].append(score_ep_proteos[p])
                        all_scores[i].append(self.rt_boundaries[p][0])
                        all_scores[i].append(self.rt_boundaries[p][1])

                        all_bounds[i].append((self.rt_boundaries[p][0], self.rt_boundaries[p][1]))

                fig = self.plot_elution_profiles(self.proteoform_subset, rt_values=all_rts, count=g)
                fig.write_image("images/fig_" + f"{group_number:03}" + "_" + f"{g:04}" + ".png")
            for p in range(len(self.proteoform_subset)):
                if self.rt_boundaries[p][0] != 0:
                    self.rt_boundaries[p][0] = mean([b[p][0] for b in all_bounds])
                    self.rt_boundaries[p][1] = mean([b[p][1] for b in all_bounds])

                    print(self.proteoform_subset[p].get_modification_brno())
                    print([b[p][0] for b in all_bounds])
                    print("lower ", mean([b[p][0] for b in all_bounds]))
                    print("higher ", mean([b[p][1] for b in all_bounds]))

            self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
            self.update_psms_ratio_subset(self.spectra_subset)
            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

            fig = self.plot_elution_profiles(self.proteoform_subset, rt_values=all_rts, count=g)
            fig.write_image("images/fig_" + f"{group_number:03}" + "_final" + ".png")

    # ------------------------------- TEMP VIZ FUNC ------------------------------ #

    def plot_elution_profiles(self, proteoforms_input, rt_values, count=0, function_values="NA"):

        # Get plot boundaries
        x_min_max = [min(rt_values) - 100, max(rt_values) + 100]
        y_min_max = [-150, 100]

        # Instanciate figure
        fig = go.Figure()
        cols = [
            "#9b2226",
            "#005f73",
            "#ee9b00",
            "#0a9396",
            "#94d2bd",
            "#ca6702",
            "#e9d8a6",
            "#bb3e03",
            "#001219",
        ]
        cols_n = 0

        # Plot each proteoforms:
        for proteo in proteoforms_input:
            if proteo.min_bound_rt != 0:

                data_x_all = [psm.spectrum.get_rt() for psm in proteo.get_linked_psm()]
                data_y_all = [psm.spectrum.getPrecIntens() for psm in proteo.get_linked_psm()]
                fig.add_scatter(
                    x=data_x_all,
                    y=data_y_all,
                    mode="markers",
                    marker=dict(size=10, color="grey", opacity=0.5),
                    marker_symbol="x-open",
                    name="Spectrum Intensity unvalid",
                )

                data_y = [psm.spectrum.get_rt() for psm in proteo.get_validated_linked_psm()]
                data_y_spectrum = [psm.spectrum.getPrecIntens() for psm in proteo.get_validated_linked_psm()]
                data_y_psm = [psm.get_prec_intens_ratio() for psm in proteo.get_validated_linked_psm()]

                fig.add_scatter(
                    x=data_y,
                    y=data_y_spectrum,
                    mode="markers",
                    marker=dict(size=15, color="grey", opacity=0.5),
                    name="Spectrum Intensity",
                )

                fig.add_scatter(
                    x=data_y,
                    y=data_y_psm,
                    mode="markers",
                    marker=dict(size=12, color=cols[cols_n]),
                    name="PSM Intensity",
                )

                fig.add_vrect(
                    x0=proteo.min_bound_rt,
                    x1=proteo.max_bound_rt,
                    # annotation_text=proteo.get_modification_brno(),
                    # annotation_position="top left",
                    opacity=0.05,
                    fillcolor=cols[cols_n],
                )

                fig.add_trace(
                    go.Scatter(
                        x=[proteo.min_bound_rt - 10],
                        y=[(0 - (cols_n + 1) * 2) - 1],
                        text=[proteo.get_modification_brno()],
                        mode="text",
                    )
                )

                fig.add_shape(
                    type="rect",
                    x0=proteo.min_bound_rt,
                    y0=0 - (cols_n + 1) * 2,
                    x1=proteo.max_bound_rt,
                    y1=0 - (cols_n + 2) * 2,
                    line=dict(
                        color=cols[cols_n],
                        width=2,
                    ),
                    fillcolor=cols[cols_n],
                    layer="below",
                )

                # add lines between PSM and spectrum intens points
                for i in range(0, len(data_y), 1):
                    fig.add_scatter(
                        x=[data_y[i], data_y[i]],
                        y=[data_y_spectrum[i], data_y_psm[i]],
                        mode="lines",
                        marker=dict(size=2, color="#c9c9c9", opacity=0.5),
                        line={"dash": "dash"},
                    )

                elution_profile = proteo.get_elution_profile()
                # print(elution_profile)

                if elution_profile != None:  # if elution profile model has been computed add line to the plot
                    data_x_elution_profile = list(range(int(x_min_max[0]), int(x_min_max[1]), 1))
                    data_y_elution_profile_fitted, params_fitted = list(
                        elution_profile.get_y_serie(data_x_elution_profile, method="fitted")
                    )
                    data_y_elution_profile_estimated, params_fitted = list(
                        elution_profile.get_y_serie(data_x_elution_profile, method="estimated")
                    )

                    # print(data_y_elution_profile_fitted)

                    if data_y_elution_profile_fitted[0] != None:
                        fig.add_scatter(
                            x=data_x_elution_profile,
                            y=data_y_elution_profile_fitted,
                            mode="lines",
                            marker=dict(size=15, color=cols[cols_n]),
                            name="Fitted Parameters",
                            line_shape="spline",
                        )

                    if data_y_elution_profile_estimated[0] != None:
                        fig.add_scatter(
                            x=data_x_elution_profile,
                            y=data_y_elution_profile_estimated,
                            mode="lines",
                            marker=dict(size=3, color=cols[cols_n]),
                            line={"dash": "dash"},
                            name="Estimated Parameters",
                            line_shape="spline",
                        )

                cols_n += 1

        fig.update_layout(
            title=go.layout.Title(
                font=dict(
                    family="Courier New, monospace",
                    size=10,
                )
            )
        )

        fig.update_layout(template="plotly_white", height=1500, width=1500, title=f"Iteration: {count}")

        return fig

    # ---------------------------- METHODS FOR EXPORT ---------------------------- #

    def result_dataframe_pfq1_format(self):
        df = pd.DataFrame(
            columns=(
                "protein",
                "sequence",
                "brno",
                "proforma",
                "intensity",
                "linked_psm",
                "linked_psm_validated",
                "rt_peak",
                "auc",
            )
        )

        for proteo in self.proteoforms.values():

            # get Elution profile max peak rt
            if proteo.get_elution_profile() != None:
                rt_peak = proteo.get_elution_profile().get_x_at_max_y()
                auc = proteo.get_elution_profile().get_auc(rt_peak - 1000, rt_peak + 1000)  # TODO hard coded
            else:
                rt_peak = "NA"
                auc = "NA"

            df.loc[len(df)] = [
                proteo.get_protein_ids(),
                proteo.peptideSequence,
                proteo.get_modification_brno(),
                proteo.get_modification_proforma(),
                proteo.get_proteoform_total_intens(),
                len(proteo.get_linked_psm()),
                len(proteo.get_validated_linked_psm()),
                rt_peak,
                auc,
            ]

        return df

    def get_dataframe_fragment_annotation(self, max_rank=1):
        "Create a dataframe where each row is the annotation for one PSM (TAKES FIRST RANK)"

        df = pd.DataFrame(
            columns=(
                "scan",
                "sequence",
                "brno",
                "proforma",
                "prec_intens",
                "sum_intens_frag",
                "prec_in_msms",
            )
        )

        for spectrum in self.spectra.values():

            for psm in spectrum.get_psms()[:max_rank]:

                # print(psm.spectrum.get_id())
                # print(psm.proteoform.get_protein_ids())
                # add information spectrum and psm
                dict_add = {
                    "scan": psm.spectrum.get_id(),
                    "sequence": psm.proteoform.peptideSequence,
                    "brno": psm.get_modification_brno(),
                    "proforma": psm.proteoform.get_modification_proforma(),
                    "prec_intens": psm.spectrum.getPrecIntens(),
                    "sum_intens_frag": psm.spectrum.getSumIntensFrag(),
                    "prec_in_msms": psm.spectrum.get_sum_intens_frag_at_mz(
                        [psm.spectrum.getPrecMz() - 20, psm.spectrum.getPrecMz() + 20]
                    ),
                }

                # add annotation:
                # print(psm.get_annotation())
                for frag_type in psm.get_annotation().values():
                    # print(frag_type)

                    annotation_frag_type = {
                        frag_type["fragCode"][i]: frag_type["intens"][i]
                        for i in range(0, len(frag_type["fragCode"]))
                    }

                    # print(annotation_frag_type)
                    dict_add.update(annotation_frag_type)

                # print(dict_add)
                df_add = pd.DataFrame(dict_add, index=[0])

                dfs = [df, df_add]

                df = pd.concat(dfs)

        return df
