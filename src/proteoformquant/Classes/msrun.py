### Import ###

# MS file parsers
import re
from os.path import splitext
from logging import warning
from pyteomics import mzid
from pyteomics import mgf
from pyteomics import mzml
from collections import defaultdict

# import ms_deisotope

# Data manipulations
from pickle import TRUE
import numpy as np
import pandas as pd

# Parallelization and multiprocessing
import multiprocessing as mp

# Misc
from statistics import mean, stdev
from tqdm import tqdm
import networkx as nx

try:  # local modules
    from proteoformquant.Classes.spectrum import Spectrum
    from proteoformquant.Classes.proteoform import Proteoform
    from proteoformquant.Utils import constant
    from proteoformquant.Utils import exception
except ImportError:
    from Classes.spectrum import Spectrum
    from Classes.proteoform import Proteoform
    from Utils import constant
    from Utils import exception


# Visualization (TEMPORARY)
import plotly.graph_objects as go
import plotly.graph_objs as go
import matplotlib.pyplot as plt


class Msrun:
    """A class containing all informations corresponding to an MSrun and the methods for processing.

    The Msrun class stores all data corresponding to an MS run, this includes the spectra and their identification.
    This class contains all the main methods called by the main script "proteoformquant.py".
    """

    def __init__(self, run_id: str = "Default run ID", params={}, params_over={}, verbose=True):
        self.run_id = run_id
        self.ident_fn: str = "Not Specified"
        self.spectra_file: str = "Not Specified"
        self.verbose: bool = verbose

        self.spectra: dict(Spectrum) = {}
        self.proteoforms: dict(Proteoform) = {}

        self.proteoform_isobaric_group = []

        self.mod_mass_to_name = {}  # store mass to name to limit queries to unimod DB

        # log
        self.n_scans_in_mgf = 0
        self.log = [
            [
                "processing_step",
                "self.n_scans_in_mgf",
                "n_spectra",
                "n_spectra_w_id",
                "n_spectra_chimeric" "n_psms",
                "n_psms_valid",
                "n_psms_valid_r1",
                "n_proteoforms",
                "n_proteoforms_w_valid_psm",
            ]
        ]

        # Load parameters
        for key, value in params.items():
            setattr(self, key, value)

        # Parameters overwrite

        for key, value in params_over.items():
            if key in self.__dict__.keys():
                try:  # convert to float if possible
                    value = float(value)
                    if value - int(value) == 0:  # convert to int decimal is zero
                        value = int(value)
                except ValueError:
                    pass
                setattr(self, key, value)
            else:
                warning(f"\n The argument {key} is not a valid parameter and will be ignored")

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
        """
        Returns the retention time range of the MS run.

        This function looks throught the retention times of all spectra in self.spectra,
        and returns the minimum and maximum rentention time values.

        Parameters
        ----------

        Returns
        -------
        tuple of (float, float)
            The minimum and maximum rentention times
        """

        mn = min([spectrum.get_rt() for spectrum in self.spectra.values()])
        mx = max([spectrum.get_rt() for spectrum in self.spectra.values()])
        return (mn, mx)

    def get_mz_range(self):
        """
        Returns the mass to charge range of precursors in the MS run.

        This function looks throught the mass to charge of spectra precursors in self.spectra,
        and returns the minimum and maximum mass to charge values

        Parameters
        ----------

        Returns
        -------
        tuple of (float, float)
            The minimum and maximum mass to charge
        """
        mn = min([spectrum.get_prec_mz() for spectrum in self.spectra.values()])
        mx = max([spectrum.get_prec_mz() for spectrum in self.spectra.values()])
        return (mn, mx)

    def get_dataset_metrics(self):
        """Returns a set of metrics about the run as a dictionnary

        Parameters
        ----------

        Returns
        -------
        dict of (str, float)
            A set of metrics for the MSrun
        """

        return {
            "TotalSpectra": len(self.spectra),
            "AssignedSpectra(method1)": len(
                [spectrum for spectrum in self.spectra.values() if spectrum.get_number_validated_psm() > 0]
            ),
            "AssignedSpectra(method2)": 0,
            # "UnassignedSpectra": len([spectrum for spectrum in self.proteoform0.linkedSpectra]),
            "ChimericSpectra": len(
                [spectrum for spectrum in self.spectra.values() if spectrum.get_number_validated_psm() > 1]
            ),
            "TotalProteoforms": len(self.proteoforms) - 1,
        }

    def add_metrics_to_log(self, processing_step="Not given"):
        """Append the metrics from self.get_dataset_metrics() to self.log.

        This function compute a set of metrics from the MS run and add them as a new line in the dataframe self.log  .

        Parameters
        ----------
        processing_step : str
            The string to be added as the rowname for the set of metrics

        Returns
        -------
        dict of (str, float)
            A set of metrics on the MSrun
        """

        n_spectra = len(self.spectra)
        n_spectra_w_id = len([s for s in self.spectra.values() if s.get_number_validated_psm() > 0])
        n_spectra_chimeric = len([s for s in self.spectra.values() if s.get_number_validated_psm() > 1])
        n_psms = sum([len(spectrum.psms) for spectrum in self.spectra.values()])
        n_psms_valid = sum([spectrum.get_number_validated_psm() for spectrum in self.spectra.values()])
        n_psms_valid_r1 = sum(
            [
                len([p for p in spectrum.psms if p.is_validated and p.get_rank() == 1])
                for spectrum in self.spectra.values()
            ]
        )
        n_proteoforms = len(self.proteoforms)
        n_proteoforms_w_valid_psm = len(
            [p for p in self.proteoforms.values() if p.get_number_validated_linked_psm() > 0]
        )

        metrics_list = [
            processing_step,
            self.n_scans_in_mgf,
            n_spectra,
            n_spectra_w_id,
            n_spectra_chimeric,
            n_psms,
            n_psms_valid,
            n_psms_valid_r1,
            n_proteoforms,
            n_proteoforms_w_valid_psm,
        ]

        self.log.append(metrics_list)
        print(f"Dataset consist of {n_spectra} spectra with {n_psms} psms for {n_proteoforms} proteoforms")

        return metrics_list

    # ----------------------------------- MAIN ---------------------------------- #

    # ------------------------- DATA PREPARATION METHODS ------------------------ #

    def read_mzid(self, ident_fn):
        """
        Read a spectra identification file in mzIdenMl and stores the data .

        This function reads the information from a spectra identification file mzIdenMl (.mzid) format and stores the information as several classes.
        Path of the identification is stored in self.inputFn
        Spectrum are stored as spectrum.Spectrum objects in self.spectra.
        Peptide spectrum matches information are stored as psm.Psm withing the Spectrum objects in self.spectra.

        Parameters
        ----------
        ident_fn : str
            The path to the .mzid file to be read.

        Returns
        -------

        """
        self.ident_fn = ident_fn  # Store File that has been read
        mzidObj = mzid.MzIdentML(ident_fn)  # Create a pyteomics' mzid iterator

        print("\n---Loading identification data from: {0}---".format(self.ident_fn))

        # Iterate over identificationresults in mzid and create Spectrum objects for each
        for identMzid in tqdm(mzidObj, disable=not self.verbose):
            self.spectra[identMzid["spectrumID"]] = Spectrum(
                spectrumID=identMzid["spectrumID"], identMzid=identMzid, max_rank=self.max_rank
            )
            # if self.verbose:

        pass

    def read_spectra(self, spectra_file):
        """Wrapper method (for self.read_mgf() and self.read_mzml()) to read in spectra from a file, and adds information to the list of spectra objects.

        Reads in spectra from a file, and adds them to the list of spectra. The type
        of file is determined by the extension. Currently supported file types are
        .mgf and .mzml.

        Args:
            spectra_file: The name of the file spectra file.
        """
        extension = splitext(spectra_file)[1]
        if extension.lower() == ".mgf":
            self.read_mgf(spectra_file)
        if extension.lower() == ".mzml":
            self.read_MZML(spectra_file)

    def read_spectra(self, spectra_file):
        extension = splitext(spectra_file)[1]
        if extension.lower() == ".mgf":
            self.read_mgf(spectra_file)
        if extension.lower() == ".mzml":
            self.read_MZML(spectra_file)

    def read_mgf(self, spectra_file):
        """
        Add informations from an mgf file to the spectrum objects.

        Parameters
        ----------
        spectra_file : str
            The path to the .mgf file to be read.

        Returns
        -------
        """
        self.spectra_file = spectra_file  # Store File name that has been read
        self.spectra_source = mgf.IndexedMGF(spectra_file)  # Create a pyteomics' mgf iterator
        self.indices = defaultdict(dict)  # Create a dictionary to store the indices of the spectra
        self._index_MGF()

        print("\n---Loading spectrum data from: {0}---".format(self.spectra_file))

        for spec_id, spec_obj in tqdm(self.spectra.items(), disable=not self.verbose):
            # Different ways to retrieve a spectrum from an mgf file:
            spec_found = False

            # print("spectrum title: ", spec_obj.spectrum_title)
            # print("spectrum id: ", spec_id)

            # 1. By full title
            if spec_found == False:
                try:
                    spec = self.spectra_source.get_spectrum(spec_obj.spectrum_title)

                    print("spectrum: ", spec, "has been found by title")
                    spec_found = True
                except KeyError:
                    pass

            # 2: By ID in IndexedMGF
            if spec_found == False and type(spec_id) is str:
                try:
                    spec = self.spectra_source.get_by_id(spec_id)
                    spec_found = True
                except KeyError:
                    pass

            # 3: By index
            if spec_found == False and type(spec_id) is int:
                try:
                    spec = self.spectra_source.get_by_index(self.indices["scan"][spec_id])
                    spec_found = True
                except KeyError:
                    pass

            # 4: recover

            if spec_found == False and type(spec_id) is str:
                try:
                    match_index = re.match(r"index(?:\s+)?=(?:\s+)?(\d+)", spec_id)
                    match_scan = re.match(r"scan(?:\s+)?=(?:\s+)?(\d+)", spec_id)

                    if not match_index is None:
                        spec_found = True
                        spec = self.spectra_source.get_by_index(int(match_index.group(1)))

                    # elif not match_scan is None:
                    #     spec = self.spectra_source.get_by_index(int(match_scan.group(1)))

                    # elif spec_id in self.indices["title"].keys():
                    #     spec = self.spectra_source.get_by_index(self.indices["title"][spec_id])

                    # elif type(spec_id) is int or spec_id.isdigit():
                    #     spec = self.spectra_source.get_by_index(self.indices["scan"][int(spec_id)])

                except KeyError:
                    pass

            if spec_found == True:
                spec_mz = spec["params"]["pepmass"][0]
                if round(spec_mz) == round(spec_obj.experimentalMassToCharge):
                    self.spectra[spec_id].set_spec_data_mgf(spec)
                else:
                    warning(
                        f"Precursor masses does not match for spectrum_id: {spec_id} ({spec_mz}) {round(spec_mz)} and psm ({spec_obj.experimentalMassToCharge})"
                    )
                    self.spectra[spec_id].set_spec_data_mgf(spec)

            else:
                raise exception.ProteoformquantError(f"Cannot find scan in mgf from spectrum_id={spec_id}")

            print("spectrum has been correctly matched to the spectrum object")

    def _index_MGF(self):
        for index in range(len(self.spectra_source)):
            params = self.spectra_source[index]["params"]
            if "title" in params.keys():
                self.indices["title"][params["title"]] = index

            if "scans" in params.keys():
                self.indices["scan"][int(params["scans"])] = index

    def read_MZML(self, spectra_file):
        self.spectra_file = spectra_file  # Store File name that has been read
        self.spectra_source = mzml.MzML(spectra_file)
        self.indices = defaultdict(dict)
        self._index_MZML()

        print("\n---Loading spectrum data from: {0}---".format(self.spectra_file))

        for spec_id, spec_obj in tqdm(self.spectra.items(), disable=not self.verbose):
            # Diffrent ways to retrievea spectrum from an mgf file:
            spec_found = False

            # 1: By ID
            if spec_found == False:
                try:
                    spec = self.spectra_source.get_by_id(spec_obj.spectrum_title)
                    spec_found = True
                except KeyError:
                    pass

            # 2: by ID from spec_id
            if spec_found == False and type(spec_id) is str:
                try:
                    spec = self.spectra_source.get_by_id(spec_id)
                    spec_found = True
                except KeyError:
                    pass

            # 2: By index
            if spec_found == False and type(spec_id) is int:
                try:
                    spec = self.spectra_source.get_by_index(self.indices["scan"][spec_id])
                    spec_found = True
                except KeyError:
                    pass

            match_index = re.match(r"index(?:\s+)?=(?:\s+)?(\d+)", spec_id)
            match_scan = re.match(r"scan(?:\s+)?=(?:\s+)?(\d+)", spec_id)

            # 4:
            if spec_found == False and type(spec_id) is str:
                try:
                    if not match_index is None:
                        spec = self.spectra_source.get_by_index(int(match_index.group(1) + 1))
                        spec_found = True
                    # elif not match_scan is None:
                    #     spec = self.spectra_source.get_by_index(int(match_scan.group(1)+1))
                    #     spec_found = True
                    # elif spec_id in self.indices["title"].keys():
                    #     spec = self.spectra_source.get_by_index(self.indices["title"][spec_id])
                    #     spec_found = True
                    # elif type(spec_id) is int or spec_id.isdigit():
                    #     spec = self.spectra_source.get_by_index(self.indices["scan"][int(spec_id)])
                    #     spec_found = True
                    else:
                        pass
                except KeyError:
                    pass

            # Check validity of spectrum
            spec_mz = spec["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0][
                "selected ion m/z"
            ]
            if spec_found == True and spec_mz == spec_obj.experimentalMassToCharge:
                self.spectra[spec_id].set_spec_data_mgf(spec)
            elif spec_found == True:
                raise exception.ProteoformquantError(
                    f"Precursor masses does not match for spectrum_id: {spec_id} ({spec_mz}) and psm ({spec_obj.experimentalMassToCharge})"
                )
            else:
                raise exception.ProteoformquantError(f"Cannot find scan in mzML from spectrum_id={spec_id}")

        #     if type(spec_id) is int:
        #         spec = self.spectra_source.get_by_index(self.indices["scan"][spec_id])
        #     elif type(spec_id) is str:
        #         # in case of whitespaces
        #         spec_id = spec_id.strip()
        #         try:
        #             try:
        #                 # the default way
        #                 spec = self.spectra_source.get_by_id(spec_id)
        #             except KeyError:
        #                 spec = self.spectra_source.get_by_id(spec_id.split("=")[1])
        #         except KeyError:
        #             # trying to recover

        #             match_index = re.match(r"index(?:\s+)?=(?:\s+)?(\d+)", spec_id)
        #             match_scan = re.match(r"scan(?:\s+)?=(?:\s+)?(\d+)", spec_id)

        #             if not match_index is None:
        #                 spec = self.spectra_source.get_by_index(int(match_index.group(1)))

        #             elif not match_scan is None:
        #                 spec = self.spectra_source.get_by_index(int(match_scan.group(1) - 1))

        #             elif type(spec_id) is int or spec_id.isdigit():
        #                 spec = self.spectra_source.get_by_index(self.indices["scan"][int(spec_id)])

        #             else:
        #                 raise KeyError(f"Cannot infer MzML scan from spectrum_id= {spec_id}")

        #     else:
        #         raise TypeError(f"Unsupported id type ({type(spec_id)}): should be int or string")

        #     self.spectra[spec_id].set_spec_data_mzml(spec)

    def _read_by_id_MZML(self, spectra_file):
        self.spectra_file = spectra_file  # Store File name that has been read
        self.spectra_source = mzml.MzML(spectra_file)
        self.indices = defaultdict(dict)
        self._index_MZML()

        print("\n---Loading spectrum data from: {0}---".format(self.spectra_file))

        for spec_id in tqdm(
            self.spectra, disable=not self.verbose
        ):  # Take into account how DBSEs store spectra ids
            spec = self._get_by_id_MZML(spec_id)
            self.spectra[spec_id].set_spec_data_mzml(spec)

    def _get_by_id_MZML(self, id_string):
        if type(id_string) is int:
            return self.spectra_source.get_by_index(self.indices["scan"][id_string])
        elif type(id_string) is str:
            # in case of whitespaces
            id_string = id_string.strip()
            try:
                # the default way
                return self.spectra_source.get_by_id(id_string)
            except KeyError:
                # trying to recover
                match_index = re.match(r"index(?:\s+)?=(?:\s+)?(\d+)", id_string)

                if not match_index is None:
                    return self.spectra_source.get_by_index(int(match_index.group(1)))

                # if id_string in self.indices['title'].keys():
                #    return self.spectra_source.get_by_index(self.indices['title'][id_string])

                if type(id_string) is int or id_string.isdigit():
                    return self.spectra_source.get_by_index(self.indices["scan"][int(id_string)])

                raise KeyError(f"Cannot infer MzML scan from spectrum_id={id_string}")

        else:
            raise TypeError(f"Unsupported id type ({type(id_string)}): should be int or string")

    def _index_MZML(self):
        for spectrum in self.spectra_source:
            spectrumID = spectrum["id"]
            index = spectrum["index"]

            scan_match = re.match(r".+scan(?:\s+)?=(?:\s+)?(\d+)", spectrumID)

            if not scan_match is None:
                self.indices["scan"][int(scan_match.group(1))] = index

            self.indices["id"][spectrumID] = index

    def add_proteoforms(self):
        """
        From spectrum objects in self.spectra instanciates proteoform.Proteoform objects and stores them self.proteoforms

        This functions looks through all Psm objects in Spectrum objects in self.spectra and instanciate a Proteoform for each unique
        proteoform in the Psms.
        """
        i = 0
        print("\n---Instanciating proteoform objects---")

        for specID in tqdm(self.spectra, disable=not self.verbose):
            i += 1

            for psm in self.spectra[specID].get_psms():
                proforma = psm.get_modification_proforma(self.mod_mass_to_name)
                brno = psm.get_modification_brno()
                seq = psm.getPeptideSequence()

                # Create list with associated protein description (try except as variable name could be different depending on the DBSE used)
                try:
                    protein_ids = [ref["accession"] for ref in psm.PeptideEvidenceRef]
                except KeyError:
                    try:
                        protein_ids = [ref["protein description"] for ref in psm.PeptideEvidenceRef]
                    except KeyError:
                        protein_ids = [ref["DBSequence_Ref"] for ref in psm.PeptideEvidenceRef]

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

    def match_fragments(self, internal=False):
        """
        Generates theoretical fragments and annotated the spectra.

        For each Proteoform objects, generates the set of theoretical fragments
        If spectrum and identification data are provided in a spectrum object,
        get and store the annotated fragments for each PSM

        Parameters
        ----------
        internal : bool
            Should internal ions be annotated (default: false, only terminal framgent).
        """

        print("\n---Generating theoretical fragments---")
        for proteoID in tqdm(self.proteoforms, disable=not self.verbose):
            self.proteoforms[proteoID].compute_theoretical_fragments(self.fragments_types)

        print("\n---Matching fragments---")
        for spectrumID in tqdm(self.spectra, disable=not self.verbose):
            self.spectra[spectrumID].annotateFragPsm(frag_mz_tol=self.frag_mz_tol)

            self.spectra[spectrumID].set_sum_intens_annot_frag()

        pass

    def scale_precursor_intensities(self, target_range=[0, 100]):
        print("\n---Scaling precursor intensities---")

        """
        Filter precursor intensity for outliers and scale the values. 

        This methods looks through all precursor intensities and detect the outliers (>4stdev from distribution), 
        spectra with outlying precursor intensity are removed. The remaining intensity are scale to the target_range. 

        Parameters
        ----------
        target_range : list of [int, int]
            Minimum and maximum intensity value to scale to
        """
        spectra_prec_intens = [spectrum.get_prec_intens() for spectrum in self.spectra.values()]

        if len(spectra_prec_intens) > 2:
            # remove extreme values:
            highest_intens = mean(spectra_prec_intens) + (4 * stdev(spectra_prec_intens))
            lowest_intens = 0

            spectra_to_remove = []
            for spectrum in self.spectra.values():
                if spectrum.get_prec_intens() > highest_intens or lowest_intens > spectrum.get_prec_intens():
                    spectra_to_remove.append(spectrum)

            print(
                f"{len(spectra_to_remove)} outlying spectra where excluded (intensity range accepted: {lowest_intens}:{highest_intens}"
            )
            for spectrum in spectra_to_remove:
                self.remove_spectrum(spectrum)

            ###

            spectra_prec_intens = [spectrum.get_prec_intens() for spectrum in self.spectra.values()]

            max_prec_intens = max(spectra_prec_intens)
            min_prec_intens = min(spectra_prec_intens)

            factor = (target_range[1] - target_range[0]) / (max_prec_intens - min_prec_intens)

            for spectrum in tqdm(self.spectra.values(), disable=not self.verbose):
                spectrum.precIntens = factor * (spectrum.get_prec_intens() + target_range[0])

    def fdr_filtering(self):
        """
        Compute the FDR and remove psms below score threshold at FDR = self.fdr_threshold

        Parameters
        ----------
        decoy_tag: str
            String in protein_description indicative of decoy sequence
        score_name: str
            Name of the score reported by the database search engine (comet: 'Comet:xcorr', msamanda:'Amanda:AmandaScore')
        """
        
        #var
        
        decoy_tag = self.decoy_tag
        score_name = self.score_name
        

        # generate list of all psms (IDEA: could be method)
        psms_obj = []
        psms_score = []
        psms_isdecoy = []

        print("\n---FDR filtering---")

        for spectrum in tqdm(self.spectra.values(), disable=not self.verbose):
            for psm in spectrum.get_psms():
                psms_obj.append(psm)
                try:
                    psms_score.append(psm.__dict__[score_name])
                except KeyError:
                    print(
                        "The score '{0}' was not found and fdr fitering has been aborted".format(score_name)
                    )
                    return 0

                isdecoy = False
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
                psms_obj[i].exclude()  # remove decoy psms for further analysis
            else:
                target_hits += 1

            if decoy_hits / target_hits < self.fdr_threshold:
                score_fdr = psms_score[i]

        idx_remove = [idx for idx, element in enumerate(psms_score) if element < score_fdr]

        psms_remove = [psms_obj[index] for index in idx_remove]

        for psm in psms_remove:
            self.remove_psm(psm)

        if decoy_hits == 0:
            print("No decoys hit were found")
        if score_fdr == 0:
            print("Warning: FDR threshold score calculated as 0")

        print(
            f"{decoy_hits} decoys and {target_hits} targets hits found, FDR {self.fdr_threshold} for score: {score_fdr}"
        )

        self.add_metrics_to_log(processing_step="fdr_filtering")

    def filter_psms_low_score(self, deviation):
        """
        Filter Psms object based on score

        Parameters
        ----------
        min_rt: int
            Minimum RT value tolerated
        max_rt: int
            Maximum RT value tolerated
        """

    def retention_time_window_filtering(self, min_rt, max_rt):
        """
        Removes spectrum objects from self.Spectra based on rentention time threshold.

        Parameters
        ----------
        min_rt: int
            Minimum RT value tolerated
        max_rt: int
            Maximum RT value tolerated
        """
        spectra_to_remove = []
        for spectrum in self.spectra.values():
            if spectrum.get_rt() > max_rt or spectrum.get_rt() < min_rt:
                spectra_to_remove.append(spectrum)

        print(f"Removed {len(spectra_to_remove)} spectra based on rentention time")

        for spectrum in spectra_to_remove:
            self.remove_spectrum(spectrum)

        self.add_metrics_to_log(processing_step="RT_window_filtering")

    def deconvolute_ms2(self):
        """
        Perform deconvolution of msms spectra and replace the peak list with deconvoluted peak list.
        """

        print("\n---Deconvoluting spectra---")

        for spectrum in tqdm(self.spectra.values(), disable=not self.verbose):
            peaks = ms_deisotope.deconvolution.utils.prepare_peaklist(
                zip(spectrum.get_frag_intens(), spectrum.get_frag_mz())
            )
            deconvoluted_peaks, targeted = ms_deisotope.deconvolute_peaks(
                peaks,
                averagine=ms_deisotope.peptide,
                scorer=ms_deisotope.MSDeconVFitter(10.0),
                verbose=True,
            )
            self.fragMz = [p.mz for p in deconvoluted_peaks.peaks]
            self.fragIntens = [p.intensity for p in deconvoluted_peaks.peaks]

        # return mzs, intensities

    # ------------------------------ Remove Methods ------------------------------ #
    # The following method are used for filtering steps of the analysis

    def remove_proteoform(self, proteoform):
        """Removes a given proteform object and its associated PSMs,
        will also remove spectrum if this spectrum doesn't have any identification

        Parameters
        ----------
        proteoform: proteoform.Proteoform
            The instance of Proteoform to be removed

        """

        psms_to_rm = proteoform.get_all_linked_psm()

        for psm in psms_to_rm:  # remove psm in spectrum
            psm.spectrum.psms.remove(psm)
            if len(psm.spectrum.psms) == 0:  # remove spectra if no psms are left
                self.spectra[psm.spectrum.get_id()]
                del psm.spectrum

        del self.proteoforms[proteoform.get_modification_proforma()]
        del proteoform

    def remove_psm(self, psm):
        spectrum = psm.spectrum
        proteoform = psm.proteoform

        spectrum.psms.remove(psm)
        if len(spectrum.psms) == 0:  # remove spectra if no psms are left
            self.spectra[spectrum.get_id()]
            del psm.spectrum

        if isinstance(proteoform, Proteoform):
            proteoform.linkedPsm.remove(psm)  # remove proteoform if no psms are left
            if len(proteoform.linkedPsm) == 0:
                self.proteoforms[proteoform.get_modification_proforma()]
                del proteoform

        del psm

    def remove_spectrum(self, spectrum):
        """Removes a given spectrum object and its associated PSMs,
        will also remove the proteoform object associated to a removed PSM
         if this proteoform doesn't have any associated PSM.

        Parameters
        ----------
        spectrum: spetrum.Spectrum
            The instance of Spectrum to be removed

        """

        psm_to_rm = []
        for psm in spectrum.psms:
            psm_to_rm.append(psm)

        for psm in psm_to_rm:
            proteoform = psm.proteoform

            if isinstance(proteoform, Proteoform):
                proteoform.linkedPsm.remove(psm)  # remove proteoform if no psms are left
                if len(proteoform.linkedPsm) == 0:
                    self.proteoforms[proteoform.get_modification_proforma()]
                    del proteoform

            del psm

        del self.spectra[spectrum.get_id()]
        del spectrum

    # ------------------------- METHODS FOR QUANTIFICATION ------------------------ #

    def update_proteoforms_elution_profile(self):
        """For each Proteoform object in self.proteoforms model the elution profile

        This method calls the method model_elution_profile() for each proteoform object in self.proteoform.
        """
        print("\n---Updating proteoforms envelopes---")

        # TODO add a function thjat set the bounds based on the entire set of envelopes
        for proteoID in tqdm(self.proteoforms, disable=not self.verbose):
            self.proteoforms[proteoID].model_elution_profile(self.elution_profile_score_threshold)

        pass

    def update_proteoform_intens(self):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        print("\n---Updating Proteoforms Total Intensity---")
        for proteoID in tqdm(self.proteoforms, disable=not self.verbose):
            self.proteoforms[proteoID].update_proteoform_total_intens()
        self.add_metrics_to_log(processing_step="Update_proteoform_intens")
        pass

    def update_psm_validation(self):
        """Updates psm.is_validated in each PSM in self.Proteoform"""
        for proteoID in tqdm(self.proteoforms, disable=not self.verbose):
            self.proteoforms[proteoID].update_proteoform_psm_validation()
        pass

    def validate_all_psms(self):
        """Updates psm.is_validated to TRUE in each PSM in self.Proteoform"""
        for spectrum in tqdm(self.spectra.values(), disable=not self.verbose):
            for psm in spectrum.get_psms():
                psm.is_validated = True

    def unvalidate_psms(self, spectra_subset=[None]):
        """Updates psm.is_validated to TRUE for PSMs in spectra_subset

        Parameters
        ----------
        spectra_subset: list of spetrum.Spectrum
            list of Spectrum in which PSM will be unvalidated
        """
        if spectra_subset[0] == None:
            spec_list = self.spectra.values()
        else:
            spec_list = spectra_subset

        for spectrum in spec_list:
            for psm in spectrum.get_psms():
                psm.is_validated = False

    def validate_psms_rank_1(self, spectra_subset=[None]):
        """Updates psm.is_validated to TRUE for PSMs or rank 1 in spectra_subset

        Parameters
        ----------
        spectra_subset: list of spetrum.Spectrum
            list of Spectrum in which PSM of rank 1 will be unvalidated
        """

        if spectra_subset[0] == None:
            spec_list = self.spectra.values()
        else:
            spec_list = spectra_subset

        for spectrum in spec_list:
            for psm in spectrum.get_psms():
                if psm.get_rank() == 1:
                    psm.is_validated = True
                if psm.get_rank() != 1:
                    psm.is_validated = False

    def validate_psms_rank_n(self, spectra_subset=[None], rank=1):
        """Updates psm.is_validated to TRUE for PSMs with higher rank than "rank"

        Parameters
        ----------
        spectra_subset: list of spetrum.Spectrum
            list of Spectrum in which PSM should be validated
        """

        if spectra_subset[0] == None:
            spec_list = self.spectra.values()
        else:
            spec_list = spectra_subset

        for spectrum in spec_list:
            for psm in spectrum.get_psms():
                if psm.get_rank() <= rank:
                    psm.is_validated = True
                if psm.get_rank() > rank:
                    psm.is_validated = False

    def filter_proteform_low_count(self):
        """Removes Proteoform object if the number of associated PSM is inferioir to min_n_psm

        Parameters
        ----------
        min_n_psm: int
            Min number of PSM for which a proteoform will be removed
        """
        print("\n---Filtering proteoforms on PSM count---")

        n_exclude = 0
        n_proteo = 0
        proteoform_to_rm = []

        for proteoform in tqdm(self.proteoforms.values(), disable=not self.verbose):
            if len(proteoform.get_linked_psm()) < self.min_n_psm:
                n_proteo += 1
                proteoform_to_rm.append(proteoform)

        for proteoform in proteoform_to_rm:
            self.remove_proteoform(proteoform)

        print(n_proteo, " proteoforms have been excluded based on count")

        self.add_metrics_to_log(processing_step="low_count_filtering")

    def update_psms_ratio_all(self):
        for spectrum in self.spectra.values():
            spectrum.update_psms_ratio()

    def update_psms_ratio_subset(self, spectra_subset):
        for spectrum in spectra_subset:
            spectrum.update_psms_ratio(0)

    # ----------------------- GROUP QUANTIFICATION ----------------------- #

    def set_proteoform_isobaric_groups(self):
        """From self.proteoforms define self.proteoform_isobaric_group where isobaric proteoform are grouped"""

        print("\n---Determining isobaric/isomeric proteoform groups---")
        all_proteoform_pairs = []  # pairs of proteoform founds together in the psms of a spectra
        all_proteoform_pairs_count = []

        for spectrum in tqdm(self.spectra.values(), disable=not self.verbose):
            proforma_psms = [psm.proteoform.get_modification_proforma() for psm in spectrum.get_psms()]
            proforma_psms = np.unique(proforma_psms)
            proforma_combinations = [
                sorted((a, b)) for idx, a in enumerate(proforma_psms) for b in proforma_psms[idx:]
            ]

            # concatenate proforma pairs to all_proteoform_pairs:
            all_proteoform_pairs = all_proteoform_pairs + proforma_combinations

        # get unique pairs and count occurences
        unique_pairs, counts = np.unique(all_proteoform_pairs, axis=0, return_counts=True)

        # keep pairs found at least self.min_connect times
        all_proteoform_pairs_filt = [
            unique_pairs[i] for i in range(len(unique_pairs)) if counts[i] >= self.min_connect
        ]

        # find groups (connected graphs) in the network defined from the edgelist "all-proteoform-pairs"
        G = nx.from_edgelist(all_proteoform_pairs_filt)

        self.isobaric_proteform_graph = G
        # # plor graph
        G.remove_edges_from(nx.selfloop_edges(G))
        # nx.draw(G, node_size=10)  # , with_labels=True)
        # plt.show()

        l = list(nx.connected_components(G))

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
            i_tot += spectrum.get_sum_intens_frag()
            i_explained += spectrum.get_sum_intens_annot_frag()

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
            v_proteo, v_subset = 0, 0

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
        coverage = [proteoform.get_coverage() for proteoform in proteoform_subset]
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

            proteoform.min_bound_rt = boundaries[b][0]
            proteoform.max_bound_rt = boundaries[b][1]

            for psm in proteoform.get_linked_psm():
                if psm.spectrum.get_rt() >= boundaries[b][0] and psm.spectrum.get_rt() <= boundaries[b][1]:
                    psm.is_validated = True
                    v += 1
                else:
                    if psm.get_rank() != 1:  # Do not unvalidate rank 1
                        psm.is_validated = False
                        uv += 1

    def optimize_proteoform_subsets(self):
        grp = 0

        print("\n---Searching for chimeric spectra---")
        for group in tqdm(self.proteoform_isobaric_group, disable=not self.verbose):
            # objects lists of subset:
            self.all_rts = []  # list of all Retention time in the spectra subset
            self.proteoform_subset = []  # Proteoform objects in the subset
            self.spectra_subset = []  # Spectra objects in the subset
            self.rt_boundaries = []  # retention time validation range for the proteoforms in the subset

            # Define proteoform subset, spectra subset, and list of RTs:
            for proforma in group:
                proteoform = self.proteoforms[proforma]
                self.proteoform_subset.append(proteoform)
                for psm in proteoform.get_linked_psm():
                    if psm.spectrum not in self.spectra_subset:
                        self.spectra_subset.append(psm.spectrum)
            self.all_rts = sorted([spectrum.get_rt() for spectrum in self.spectra_subset])

            # Sort proteoforms in subset (based on psms rank count)
            weighted_rank = [
                proteoform.get_weighted_number_linked_validated_psm(max_rank=self.max_rank)
                for proteoform in self.proteoform_subset
            ]
            zipped_rank_proteo = zip(weighted_rank, group, self.proteoform_subset)
            zipped_proteo = sorted(zipped_rank_proteo, reverse=True)
            self.proteoform_subset = [list(tuple)[-1] for tuple in zipped_proteo]

            # --------------
            # print("-------Processing group-----")
            # for p_n in range(len(self.proteoform_subset)):
            #     print(
            #         self.proteoform_subset[p_n].get_modification_proforma(),
            #         ",",
            #         zipped_proteo[p_n][0],
            #         " ",
            #         [
            #             self.proteoform_subset[p_n].get_number_linked_psm_Rx(rank=r)
            #             for r in range(1, self.max_rank)
            #         ],
            #     )

            # ---------------

            # Validate first rank:
            self.validate_psms_rank_1(spectra_subset=self.spectra_subset)
            self.update_psms_ratio_subset(self.spectra_subset)
            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

            # ---------
            # self.plot_elution_profiles(
            #     self.proteoform_subset, rt_values=self.all_rts, count=grp, plot_all=True
            # ).write_image("images/fig_" + f"{grp:03}" + "_00_0000A" + ".png")
            # --------

            # Initialize retention time boundaries from first rank EPs
            for p in range(len(self.proteoform_subset)):
                self.rt_boundaries.append([0, 0])
                self.rt_boundaries[p] = self.proteoform_subset[p].get_boundaries_of_ep()

            # Update validation from elution ranges
            self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
            self.update_psms_ratio_subset(self.spectra_subset)
            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

            # ---------
            # self.plot_elution_profiles(
            #     self.proteoform_subset, rt_values=self.all_rts, count=grp
            # ).write_image("images/fig_" + f"{grp:03}" + "_00_0000B" + ".png")
            # ---------

            # Update retention time boundaries
            for p in range(len(self.proteoform_subset)):
                self.rt_boundaries[p] = self.proteoform_subset[p].get_boundaries_of_ep()

            # Update validation from elution ranges
            self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
            self.update_psms_ratio_subset(self.spectra_subset)
            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

            # ---------
            # self.plot_elution_profiles(
            #     self.proteoform_subset, rt_values=self.all_rts, count=grp
            # ).write_image("images/fig_" + f"{grp:03}" + "_00_0000C" + ".png")
            # ---------

            # for testing:
            # for proteo in self.proteoform_subset:
            #     if proteo.get_elution_profile() != None:
            #         proteo.gap_score = proteo.get_gap_in_validated_spectra(rts=self.all_rts)

            # Unvalidate proteoform with poor elution profile scores
            for p in range(len(self.proteoform_subset)):
                proteo = self.proteoform_subset[p]
                conditions = [
                    proteo.get_fit_score() > self.min_ep_fit,
                    proteo.get_coverage() > self.min_ep_cov,
                    # proteo.get_gap_in_validated_spectra(rts=self.all_rts) > self.min_ep_gap,
                ]
                if (
                    sum(conditions) < 2
                    or proteo.get_elution_profile() == None
                    or proteo.get_elution_profile().is_modeled() == False
                ):
                    self.rt_boundaries[p] = [0, 0]
                    self.proteoform_subset[p].envelope = None

            # Update validation from elution ranges
            self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
            self.update_psms_ratio_subset(self.spectra_subset)
            # self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

            # ---------
            # self.plot_elution_profiles(
            #     self.proteoform_subset, rt_values=self.all_rts, count=grp
            # ).write_image("images/fig_" + f"{grp:03}" + "_00_0000D" + ".png")
            # ---------

            ##Stops here if not enough Spectra###
            if len(self.spectra_subset) > self.min_spectra_subset:
                # look for proteoform "hidden in higher ranks"
                n_proteo_excluded = 0
                for p in range(len(self.proteoform_subset)):
                    if n_proteo_excluded < self.max_rejected_proteo:
                        if self.rt_boundaries[p][0] == 0:  # Add new proteoform/peptidoform
                            rt_window = self.proteoform_subset[p].get_rt_range_centered(self.window_size_rt)
                            rt_window[0] = min(self.all_rts, key=lambda x: abs(x - rt_window[0]))
                            rt_window[1] = min(self.all_rts, key=lambda x: abs(x - rt_window[1]))
                            self.rt_boundaries[p] = rt_window

                            # optimization
                            self.mutate_rt_boundaries(proteo_to_mutate=p, grp=grp, iter_try=n_proteo_excluded)

                            for proteo in [
                                self.proteoform_subset[i]
                                for i in range(len(self.proteoform_subset))
                                if self.rt_boundaries[i][0] != 0
                            ]:
                                conditions = [
                                    proteo.get_fit_score() > self.min_ep_fit,
                                    proteo.get_coverage() > self.min_ep_cov,
                                    # proteo.get_gap_in_validated_spectra(rts=self.all_rts)
                                    # > self.min_ep_gap,
                                ]

                                # print(
                                #     proteo.modificationBrno,
                                #     proteo.get_fit_score(),
                                #     proteo.get_coverage(),
                                #     proteo.get_gap_in_validated_spectra(rts=self.all_rts),
                                # )

                                if (  # if score of one proteo drops below thresholds the newly added proteoform is unvalidated
                                    sum(conditions) < 2
                                    or proteo.get_elution_profile() == None
                                    or proteo.get_elution_profile().is_modeled() == False
                                ):
                                    # print(
                                    #     f"proteoform {self.proteoform_subset[p].get_modification_proforma()} UNVALIDATED, cor:{self.proteoform_subset[p].get_fit_score()}, cov:{ self.proteoform_subset[p].get_coverage()}"
                                    # )
                                    self.rt_boundaries[p] = [0, 0]
                                    self.proteoform_subset[p].envelope = None
                                    n_proteo_excluded += 1
                                    # Update validation from elution ranges
                                    self.update_proteoform_subset_validation(
                                        self.proteoform_subset, self.rt_boundaries
                                    )
                                    self.update_psms_ratio_subset(self.spectra_subset)
                                    self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

                                    break

            # Unvalidate proteoform with poor elution profile scores
            for p in range(len(self.proteoform_subset)):
                proteo = self.proteoform_subset[p]
                conditions = [
                    proteo.get_fit_score() > self.min_ep_fit,
                    proteo.get_coverage() > self.min_ep_cov,
                    # proteo.get_gap_in_validated_spectra(rts=self.all_rts) > self.min_ep_gap,
                ]
                if (
                    sum(conditions) < 2
                    or proteo.get_elution_profile() == None
                    or proteo.get_elution_profile().is_modeled() == False
                ):
                    # print("final unvalidated proteo")
                    # print(
                    #     proteo.modificationBrno,
                    #     proteo.get_fit_score(),
                    #     proteo.get_coverage(),
                    #     proteo.get_gap_in_validated_spectra(rts=self.all_rts),
                    # )
                    self.rt_boundaries[p] = [0, 0]
                    self.proteoform_subset[p].envelope = None

            # ---------
            # self.plot_elution_profiles(
            #     self.proteoform_subset, rt_values=self.all_rts, count=grp
            # ).write_image("images/fig_" + f"{grp:03}" + "_99_999" + ".png")
            # ---------
            grp += 1

        file = open(f"score_param_ep.csv", "w")
        for proteo in self.proteoforms.values():
            if proteo.get_elution_profile() != None:
                ep = proteo.get_elution_profile()
                file.write(
                    ",".join(
                        [
                            str(proteo.peptideSequence),
                            str(proteo.modificationBrno),
                            str(proteo.get_rt_center()),
                            str(proteo.get_fit_score()),
                            str(proteo.get_coverage()),
                            # str(proteo.gap_score),
                            str(ep.param_fitted[0]),
                            str(ep.param_fitted[1]),
                            str(ep.param_fitted[2]),
                            str(ep.param_fitted[3]),
                            str(proteo.get_boundaries_of_ep()[1] - proteo.get_boundaries_of_ep()[0]),
                        ]
                    )
                )
                file.write("\n")
        file.close()

    def validate_first_rank_no_id(self):
        """validate the first rank psm of all spectra without any id"""
        # file = open(f"scores_psms_all.csv", "w")
        coverage_valid = []
        for spectrum in self.spectra.values():
            for psm in spectrum.psms:
                coverage = psm.get_fragment_coverage()
                rank = psm.get_rank()
                validated = psm.is_validated
                # file.write(",".join([str(validated), str(rank), str(coverage)]))
                # file.write("\n")

                if validated:
                    coverage_valid.append(coverage)

        min_coverage = mean(coverage_valid) - stdev(coverage_valid)

        # file.close()

        for spectrum in self.spectra.values():
            if spectrum.get_number_validated_psm() == 0 and len(spectrum.get_psms()):
                if spectrum.psms[0].get_fragment_coverage() > min_coverage:
                    spectrum.psms[0].is_validated = True
                    spectrum.update_psms_ratio()

        self.add_metrics_to_log(processing_step="Validation of rank 1")

    def mutate_rt_boundaries(self, proteo_to_mutate="all", grp=0, iter_try=0):
        """
        modify the boundaries in self.rt_boundaries for the current subset of proteoform in self.proteoform_subset,
        proteo_to_mutate: if "all", mutate all boundaries for proteoform whose rt_boundaries are not (0,0), else specify the index of the roteoform to mutate
        """

        if proteo_to_mutate == "all":  # Mutate RT range of all proteoforms
            proteo_indexes = list(range(len(self.proteoform_subset)))
        elif type(proteo_to_mutate) == int:  # Mutate Only specified proteoform
            proteo_indexes = [proteo_to_mutate]
        else:
            print("error in proteoform index, cannot mutate rt boundaries")

        x = 0
        for p in proteo_indexes:
            if self.rt_boundaries[p][0] == 0 or self.rt_boundaries[p][1] == 0:
                warning("Mutating RT range that is Zero")

            cor_score_l = []
            cov_score_l = []
            rt_ranges_l = []

            for iter in range(self.n_iter):
                proteo = self.proteoform_subset[p]

                # Elution profile model range estimate:
                EP_range = proteo.get_boundaries_of_ep()

                # -------Lower RT bound mutation-------
                # Determine if increase or decrease are most likely
                if EP_range[0] == 0:
                    modifier_min = int(np.round(np.random.normal(loc=0, scale=self.rd_scale, size=1)[0]))
                elif self.rt_boundaries[p][0] < EP_range[0]:  # if spectra lower than modeled ep
                    modifier_min = int(
                        np.round(np.random.normal(loc=self.rd_loc, scale=self.rd_scale, size=1)[0])
                    )
                else:
                    modifier_min = int(
                        np.round(np.random.normal(loc=-self.rd_loc, scale=self.rd_scale, size=1)[0])
                    )

                # Apply mutation
                index_rt_min_start = min(
                    enumerate(self.all_rts), key=lambda x: abs(x[1] - self.rt_boundaries[p][0])
                )[0]
                new_index_rt_min = index_rt_min_start + modifier_min
                if new_index_rt_min >= 0 and new_index_rt_min < len(self.all_rts):
                    self.rt_boundaries[p][0] = self.all_rts[new_index_rt_min]
                    proteo.min_bound_rt = self.all_rts[new_index_rt_min]

                # -------Upper RT bound mutation-------
                # Determine if increase or decrease are most likely
                if EP_range[0] == 0:
                    modifier_max = int(np.round(np.random.normal(loc=0, scale=self.rd_scale, size=1)[0]))
                elif self.rt_boundaries[p][1] > EP_range[1]:  # if spectra lower than modeled ep
                    modifier_max = int(
                        np.round(np.random.normal(loc=-self.rd_loc, scale=self.rd_scale, size=1)[0])
                    )
                else:
                    modifier_max = int(
                        np.round(np.random.normal(loc=self.rd_loc, scale=self.rd_scale, size=1)[0])
                    )

                # Apply mutation
                index_rt_max_start = min(
                    enumerate(self.all_rts), key=lambda x: abs(x[1] - self.rt_boundaries[p][1])
                )[0]
                new_index_rt_max = index_rt_max_start + modifier_max
                if new_index_rt_max >= 0 and new_index_rt_max < len(self.all_rts):
                    self.rt_boundaries[p][1] = self.all_rts[new_index_rt_max]
                    proteo.max_bound_rt = self.all_rts[new_index_rt_max]

                # fig = self.plot_elution_profiles(self.proteoform_subset, rt_values=self.all_rts, count=iter)
                # fig.write_image(
                #     "images/fig_"
                #     + f"{grp:03}"
                #     + "_"
                #     + f"{iter_try:03}"
                #     + "_"
                #     + f"{iter+1:03}"
                #     + "_A"
                #     + ".png"
                # )

                self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
                self.update_psms_ratio_subset(self.spectra_subset)
                self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

                # fig = self.plot_elution_profiles(self.proteoform_subset, rt_values=self.all_rts, count=iter)

                # fig.write_image(
                #     "images/fig_"
                #     + f"{grp:03}"
                #     + "_"
                #     + f"{iter_try:03}"
                #     + "_"
                #     + f"{iter+1:03}"
                #     + "_B"
                #     + ".png"
                # )
                # x += 1

                # store score and corresponding RT ranges
                cor_score_l.append(proteo.get_fit_score())
                cov_score_l.append(proteo.get_coverage())
                rt_ranges_l.append((proteo.min_bound_rt, proteo.max_bound_rt))

            # Retain the RT range with the best scoring elution profile
            pass_threshold = [
                cor_score_l[i] > self.min_ep_fit and cov_score_l[i] > self.min_ep_cov
                for i in range(len(cor_score_l))
            ]
            avg_pass_score = [
                ((cor_score_l[i] + cov_score_l[i]) / 2) * pass_threshold[i] for i in range(len(cor_score_l))
            ]

            max_i = avg_pass_score.index(max(avg_pass_score))

            self.rt_boundaries[p][0] = rt_ranges_l[max_i][0]
            proteo.min_bound_rt = rt_ranges_l[max_i][0]
            self.rt_boundaries[p][1] = rt_ranges_l[max_i][1]
            proteo.max_bound_rt = rt_ranges_l[max_i][1]

            self.update_proteoform_subset_validation(self.proteoform_subset, self.rt_boundaries)
            self.update_psms_ratio_subset(self.spectra_subset)
            self.update_proteoforms_elution_profile_subset(self.proteoform_subset)

    def quantify_chimeric_spectra_no_elution_profile(self):
        """ For spectra with only rank 1 PSM validated, validate second rank PSM and quantify the chimeric spectra"""
        
        # Validate rank 2 PSMs
        for spectrum in self.spectra.values():
            if spectrum.get_number_validated_psm() == 1:
                #if there is a second rank psm
                if len(spectrum.get_psms()) > 1:
                    spectrum.get_psms()[1].is_validated = True
                    spectrum.update_psms_ratio()
                    

        self.add_metrics_to_log(processing_step="Quantification of chimeric spectra (no elution profile / Top 2 PSMs)")


    # ------------------------------- TEMP VIZ FUNC ------------------------------ #

    def plot_elution_profiles(
        self, proteoforms_input, rt_values, count=0, function_values="NA", plot_all=False
    ):
        # Get plot boundaries
        x_min_max = [min(rt_values) - 50, max(rt_values) + 50]
        y_min_max = [-150, 100]

        # Instanciate figure
        fig = go.Figure()
        cols = constant.colors
        cols_n = 0

        # Plot each proteoforms:
        for proteo in proteoforms_input:
            if (proteo.min_bound_rt != 0 and proteo.min_bound_rt is not None) or plot_all:
                data_x_all = [psm.spectrum.get_rt() for psm in proteo.get_linked_psm()]
                data_y_all = [psm.spectrum.get_prec_intens() for psm in proteo.get_linked_psm()]
                fig.add_scatter(
                    x=data_x_all,
                    y=data_y_all,
                    mode="markers",
                    marker=dict(size=10, color="grey", opacity=0.5),
                    marker_symbol="x-open",
                    name="Spectrum Intensity unvalid",
                )

                data_y = [psm.spectrum.get_rt() for psm in proteo.get_validated_linked_psm()]
                data_y_spectrum = [
                    psm.spectrum.get_prec_intens() for psm in proteo.get_validated_linked_psm()
                ]
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

                if plot_all == False:
                    min_bound_rt = proteo.min_bound_rt
                    max_bound_rt = proteo.max_bound_rt
                else:
                    min_bound_rt = x_min_max[0]
                    max_bound_rt = x_min_max[1]

                # fig.add_vrect(
                #     x0=min_bound_rt,
                #     x1=max_bound_rt,
                #     # annotation_text=proteo.get_modification_brno(),
                #     # annotation_position="top left",
                #     opacity=0.05,
                #     fillcolor=cols[cols_n],
                # )

                fig.add_trace(
                    go.Scatter(
                        x=[min_bound_rt],
                        y=[0 - (cols_n + 1) * 2],
                        text=[
                            proteo.get_modification_brno()
                            + " "
                            + str(round(proteo.get_fit_score(), 2))
                            + " "
                            + str(round(proteo.get_coverage(), 2))
                            + " "
                            + str(proteo.min_bound_rt)
                            + " "
                            + str(proteo.get_boundaries_of_ep()[0])
                            + " "
                            + str(proteo.max_bound_rt)
                            + " "
                            + str(proteo.get_boundaries_of_ep()[1])
                        ],
                        mode="text",
                        textfont_size=25,
                    )
                )

                fig.add_shape(
                    type="rect",
                    x0=min_bound_rt,
                    y0=0 - (cols_n + 1) * 2,
                    x1=max_bound_rt,
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
        # Create empty dataframe with the columns that we want
        df = pd.DataFrame(
            columns=(
                "protein",
                "sequence",
                "brno",
                "proforma",
                "proforma_full",
                "intensity",
                "intensity_r1",
                "linked_psm",
                "linked_psm_validated",
                "rt_peak",
                "auc",
                "ambiguity",
            )
        )

        # Iterate through the proteoforms
        for proteo in self.proteoforms.values():
            # Get the elution profile
            if proteo.get_elution_profile() != None:
                rt_peak = proteo.get_elution_profile().get_x_at_max_y()
                auc = proteo.get_elution_profile().get_auc()  # TODO hard coded
            else:
                rt_peak = 0
                auc = 0

            # Add the proteoform information to the dataframe
            df.loc[len(df)] = [
                # get protein list and join with ;
                ";".join(proteo.get_protein_ids()),
                proteo.peptideSequence,
                proteo.get_modification_brno(),
                "na",
                proteo.get_modification_proforma(),
                proteo.get_proteoform_total_intens(),
                proteo.get_proteoform_total_intens_r1(),
                len(proteo.get_linked_psm()),
                len(proteo.get_validated_linked_psm()),
                float(rt_peak),
                float(auc),
                int(proteo.ambiguous_spectra()),
            ]

        print(df.head())

        # Proforma without charge info
        df["proforma"] = df["proforma_full"].str.split("/").str[0]

        # merge row with identical proforma, sequence, brno and protein.
        df = df.groupby(["proforma", "sequence", "brno", "protein"]).agg(
            {
                "intensity": "sum",
                "intensity_r1": "sum",
                "linked_psm": "sum",
                "linked_psm_validated": "sum",
                "rt_peak": lambda x: ";".join(x.astype(str)),
                "auc": lambda x: ";".join(x.astype(str)),
                "ambiguity": "sum",
            }
        )

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
                    "prec_intens": psm.spectrum.get_prec_intens(),
                    "sum_intens_frag": psm.spectrum.get_sum_intens_frag(),
                    "prec_in_msms": psm.spectrum.get_sum_intens_frag_at_mz(
                        [psm.spectrum.get_prec_mz() - 20, psm.spectrum.get_prec_mz() + 20]
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

    def psm_score_dataframe(self, file_name, score_name):
        psms_df = {
            "spec": [],
            "rank": [],
            "sequence": [],
            "brno": [],
            "proforma": [],
            "score": [],
            "validated": [],
            "frag_cov": [],
        }

        for spectrum in self.spectra.values():
            for psm in spectrum.get_psms():
                psms_df["spec"].append(spectrum.get_id())
                psms_df["rank"].append(psm.get_rank())
                psms_df["sequence"].append(psm.proteoform.peptideSequence)
                psms_df["brno"].append(psm.get_modification_brno())
                psms_df["proforma"].append(psm.proteoform.get_modification_proforma())
                psms_df["score"].append(psm.__dict__[score_name])
                psms_df["validated"].append(psm.is_validated)
                psms_df["frag_cov"].append(psm.get_fragment_coverage())

        psms_df = pd.DataFrame.from_dict(psms_df, orient="index").transpose()
        return psms_df
