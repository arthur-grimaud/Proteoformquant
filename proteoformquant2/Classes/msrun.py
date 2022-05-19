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
        self.prec_mz_tol = 1.5  # in Da
        self.frag_mz_tol = 0.015  # in Da

        # Filtering and thresholds:
        self.max_rank = 5  # Maximum rank to be considered (1-based)
        self.fdr_threshold = 0.01
        self.intensity_threshold = 20000
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
                    spectrumID=identMzid["spectrumID"], identMzid=identMzid
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

        with Bar("loading spectra", max=1) as bar:

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
                bar.next()
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
        for spectrum in self.spectra.values():
            for psm in spectrum.get_psms():
                psms_obj.append(psm)

                try:
                    psms_score.append(psm.__dict__[score_name])
                except KeyError:
                    warn("The score '{0}' was not found and fdr fitering has been aborted".format(score_name))
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
            if psms_isdecoy[i]:
                decoy_hits += 1
                psms_obj[i].exclude()  # remove decoy psms from further analysis
            else:
                target_hits += 1
            if decoy_hits / target_hits > self.fdr_threshold:
                score_fdr = psms_score[i]
                break

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
        with Bar("Adding PSMs to Proteoform objects", max=1) as bar:

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
                        self.proteoforms[proforma] = Proteoform(
                            peptideSequence=seq,
                            modificationBrno=brno,
                            modificationProforma=proforma,
                            modificationDict=psm.Modification,
                            protein_ids=protein_ids,
                        ).set_color(i)

                    self.proteoforms[proforma].link_psm(psm)  # add link to Psms in Proteoform
                    psm.setProteoform(self.proteoforms[proforma])  # add link to Proteoform in Psm
                    bar.next()

    def match_fragments(self, internal=False):
        """
        For each proteoform object generates the set of theoretical fragments
        If mgf and identification data are provided in a spectrum object,
        get the and store the annotated fragments for each PSM"""

        with Bar("Generating theoretical fragments", max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].compute_theoretical_fragments(self.fragments_types)
                bar.next()

        with Bar("Matching fragments", max=1) as bar:
            for spectrumID in self.spectra:
                self.spectra[spectrumID].annotateFragPsm(frag_mz_tol=self.frag_mz_tol)
                self.spectra[spectrumID].setSumIntensAnnotFrag()
                bar.next()

        pass

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
        all_proteoform_pairs = []  # pairs of proteoform founds together in the psms of a spectra
        all_proteoform_pairs_count = []

        for spectrum in self.spectra.values():

            # print(spectrum.get_id())
            # print([psm.get_modification_brno() for psm in spectrum.get_psms()])

            proforma_psms = [psm.get_modification_proforma() for psm in spectrum.get_psms()]
            proforma_psms = np.unique(proforma_psms)
            proforma_combinations = [
                sorted((a, b)) for idx, a in enumerate(proforma_psms) for b in proforma_psms[idx + 1 :]
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
        G = nx.from_edgelist(all_proteoform_pairs)
        # # plor graph
        # nx.draw(G, node_size=10, with_labels=False)
        # plt.show()

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

    # ---------------------------- methods for scoring --------------------------- #

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
            spectrum.update_psms_ratio()

    def update_proteoforms_elution_profile_subset(self, proteoform_subset):
        """For each proteoforms in self. proteoforms model the elution profile"""
        # with Bar("Updating proteoforms envelopes", max=1) as bar:

        # TODO add a function that set the bounds based on the entire set of envelopes

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

        # for proteoform in proteoform_subset:
        #     proteoform.model_elution_profile(self.elution_profile_score_threshold)
        # bar.next()
        # pass

    def get_gap_in_validated_psms(self, proteoform_subset, spectra_subset, boundaries):
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
                ratios_missed.append(1)  # if no validated psm add perfect score
        return ratios_missed

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
        return rank_scores

    def get_intensity_explained_ratio(self, spectra_subset):
        i_tot = 0
        i_explained = 0

        for spectrum in spectra_subset:
            i_tot += spectrum.getSumIntensFrag()
            # print(spectrum.getSumIntensFrag(), " ... ", spectrum.getSumIntensAnnotFrag())
            i_explained += spectrum.getSumIntensAnnotFrag()

        # print(i_explained / i_tot)
        return i_explained / i_tot

    def get_ratio_missed_rank_1(self, spectra_subset):
        s_unvalidated = 0
        s_tot = 0

        for spectrum in spectra_subset:
            if spectrum.get_psms()[0].is_validated == False:
                s_unvalidated += 1
            s_tot += 1

        return s_unvalidated / s_tot

    def get_ratio_validated_0(self, spectra_subset):
        p_0 = 0
        p_tot = 0

        for spectrum in spectra_subset:
            for psm in spectrum.get_psms():
                if psm.is_validated == True and psm.ratio == 0:

                    p_0 += 1

                p_tot += 1

        return p_0 / p_tot

    def get_gap_in_rank_psms(self, spectra_subset):
        p_0 = 0
        p_tot = 0

        for spectrum in spectra_subset:
            val_list = [psm.is_validated for psm in spectrum.get_psms()]
            # print(val_list)
            val_index = [i for i, x in enumerate(val_list) if x]
            # print(val_index)
            if len(val_index) != 0:
                last_val_index = val_index[-1]
                # print(last_val_index)
                for p in range(int(last_val_index)):
                    if val_list[p] == False:
                        p_0 += 1
                    p_tot += 1
        if p_tot == 0:
            return 1
        return p_0 / p_tot

    def get_coverage_score(proteoform_subset):

        for proteoform in proteoform_subset:
            psms_rt = [psm.spectrum.get_rt()]

    # ------------------------------- OPTIMIZATION ------------------------------- #

    def update_psms_validation_subset_2(self, proteoform_subset, boundaries):

        p = 0
        for b in range(0, len(boundaries), 2):
            proteoform = proteoform_subset[p]
            p += 1
            v, uv = 0, 0

            proteoform.min_bound_rt = min(boundaries[b], boundaries[b + 1])
            proteoform.max_bound_rt = max(boundaries[b], boundaries[b + 1])

            for psm in proteoform.get_linked_psm():
                if psm.spectrum.get_rt() > min(
                    boundaries[b], boundaries[b + 1]
                ) and psm.spectrum.get_rt() < max(boundaries[b], boundaries[b + 1]):
                    psm.is_validated = True
                    v += 1
                else:
                    psm.is_validated = False
                    uv += 1

    def find_optimal_proteoform_set_2(self):

        # Process proteoform groups separately
        for group in self.proteoform_isobaric_group:

            # Get list of proteoform and spectra objects in the group
            self.proteoform_subset = []
            self.spectra_subset = []
            self.variables = []

            for proforma in group:
                proteoform = self.proteoforms[proforma]
                self.proteoform_subset.append(proteoform)

                rt_range = proteoform.get_rt_range_r1()
                print(rt_range)
                self.variables.append(mean(rt_range))
                self.variables.append(mean(rt_range))

                for psm in proteoform.get_linked_psm():
                    if psm.spectrum not in self.spectra_subset:
                        self.spectra_subset.append(psm.spectrum)

            # Generate list of variables to optimize in the form [Proteo1-lowerbound, Proteo1-upperbound, Proteo2...]

            gen = 10

            size_start_pop = 10
            offsprings = 60

            self.variables = np.asarray([self.variables] * size_start_pop)
            noise = np.random.normal(0, 50, self.variables.shape)

            X = self.variables + noise

            print("variable: ", self.variables)

            # Problem class for optimization (could be moved somewhere else)

            # gen = self.n_generation
            # pop = self.n_individuals

            # instanciate problem
            problem = ProblemWrapper(run=self)

            # init pop
            pop = Population.new("X", X)
            Evaluator().eval(problem, pop)

            algorithm = ES(n_offsprings=offsprings, pop_size=10, rule=1.0 / 5.0, sampling=pop)
            # algorithm = GA(pop_size=pop)

            # Run optimization
            res = minimize(problem, algorithm, ("n_gen", gen), verbose=True, save_history=True, seed=1)

            # To make a gif: !
            for g in range(gen):
                print(g)
                # print(res.history[g].__dict__["opt"][0].__dict__["X"])
                vars = res.history[g].__dict__["opt"][0].__dict__["X"]
                self.update_psms_validation_subset_2(self.proteoform_subset, vars)
                self.update_psms_ratio_subset(self.spectra_subset)
                self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

                fig = self.plot_elution_profiles(group)
                fig.write_image("images/fig_{0}.png".format(g))

            # Results
            print("Best solution found: %s" % res.X.astype(int))
            print("Proteoforms: %s" % group)
            print("Function value: %s" % res.F)
            print("Constraint violation: %s" % res.CV)

            self.update_psms_validation_subset_2(self.proteoform_subset, res.X.astype(int))  # WIP
            self.update_psms_ratio_subset(self.spectra_subset)
            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

            n_evals = np.array([e.evaluator.n_eval for e in res.history])
            opt = np.array([e.opt[0].F for e in res.history])

            plt.title("Convergence")
            plt.plot(n_evals, opt, "--")
            plt.show()

    # ------------------------------- TEMP VIZ FUNC ------------------------------ #

    def plot_elution_profiles(self, proteoforms_input):

        # Get plot boundaries
        x_min_max = [4400, 4950]
        y_min_max = [-2000000, 4000000]

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

            data_x_all = [psm.spectrum.get_rt() for psm in self.proteoforms[proteo].get_linked_psm()]
            data_y_all = [psm.spectrum.getPrecIntens() for psm in self.proteoforms[proteo].get_linked_psm()]
            spectrum_key = [psm.spectrum.id for psm in self.proteoforms[proteo].get_linked_psm()]
            fig.add_scatter(
                x=data_x_all,
                y=data_y_all,
                mode="markers",
                marker=dict(size=10, color="grey", opacity=0.5),
                marker_symbol="x-open",
                name="Spectrum Intensity unvalid",
                customdata=spectrum_key,
            )

            data_y = [psm.spectrum.get_rt() for psm in self.proteoforms[proteo].get_validated_linked_psm()]
            data_y_spectrum = [
                psm.spectrum.getPrecIntens() for psm in self.proteoforms[proteo].get_validated_linked_psm()
            ]
            data_y_psm = [
                psm.get_prec_intens_ratio() for psm in self.proteoforms[proteo].get_validated_linked_psm()
            ]

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
                x0=self.proteoforms[proteo].min_bound_rt,
                x1=self.proteoforms[proteo].max_bound_rt,
                annotation_text=self.proteoforms[proteo].get_modification_brno(),
                annotation_position="top left",
                opacity=0.1,
                fillcolor=cols[cols_n],
            )

            fig.add_shape(
                type="rect",
                x0=self.proteoforms[proteo].min_bound_rt,
                y0=0 - cols_n * 100000,
                x1=self.proteoforms[proteo].max_bound_rt,
                y1=0 - (cols_n + 1) * 100000,
                line=dict(
                    color=cols[cols_n],
                    width=2,
                ),
                fillcolor=cols[cols_n],
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

            elution_profile = self.proteoforms[proteo].get_elution_profile()
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

        fig.update_layout(template="plotly_white", height=1500, width=1500)

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
