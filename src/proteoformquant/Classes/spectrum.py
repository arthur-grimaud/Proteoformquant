from array import array
from pickle import FALSE
from itertools import compress
from scipy.optimize import nnls
import numpy as np
from numba import jit


from yaml import warnings

try:  # local modules
    from proteoformquant.Classes.psm import Psm
except ImportError:
    from Classes.psm import Psm


class Spectrum:
    def __init__(self, spectrumID, identMzid=None, max_rank=10):
        self.id = spectrumID  # Unique ID for the spectrum

        if "spectrum title" in identMzid.keys():
            self.spectrum_title = identMzid["spectrum title"]
        elif "spectrumID" in identMzid.keys():
            self.spectrum_title = identMzid["spectrumID"]
        else:
            print("Spectrum title not found in mzid file")
            raise Exception
            

        # Add spectrum experimental mz value from mzident file:
        self.experimentalMassToCharge = identMzid["SpectrumIdentificationItem"][0]["experimentalMassToCharge"]

        self.psms: list(Psm) = []  # list of PSM for that spectrum

        if identMzid != None:
            self.set_ident_data_mzid(identMzid, max_rank)

        self.sumIntensAnnotFrag = 0
        self.quant_residuals = 0

        # MultiQuant Attr (will be removed)
        self.unique_matrix = None
        self.intensity_matrix = None
        self.unique_matrix_r = None
        self.intensity_matrix_r = None
        self.equations = None
        self.variables = None
        self.proteforms_multiquant = None
        self.intervals = None

        # Missing determining ions:
        self.miss_determining_ions = False

        pass

    # Getters:

    def get_id(self):
        return self.id

    def get_prec_mz(self):
        return self.precMz

    def get_prec_intens(self):
        if self.precIntens is not None:
            return self.precIntens
        else:
            print("Precursor intensity missing, using fragment intensity sum")
            return self.get_sum_intens_frag()

    def get_rt(self):
        return self.rt

    def get_frag_intens(self):
        return self.fragIntens

    def get_frag_mz(self):
        return self.fragMz

    def get_sum_intens_frag(self):
        return sum(self.fragIntens)

    def get_sum_intens_frag_at_mz(self, mz_range):
        # return the sum of intensity at specified mz range

        intensities = [
            pair[1]
            for pair in zip(self.fragMz, self.fragIntens)
            if pair[0] > mz_range[0] and pair[0] < mz_range[1]
        ]
        return sum(intensities)

    def get_sum_intens_annot_frag(self):
        self.set_sum_intens_annot_frag()  # could be done one time every update of ratios
        return self.sumIntensAnnotFrag

    def get_psms(self):
        return [psm for psm in self.psms if psm.is_excluded == False]

    def get_validated_psm(self):
        return [psm for psm in self.get_psms() if psm.is_validated == True and psm.is_excluded == False]

    def get_number_validated_psm(self):
        return len([psm for psm in self.get_psms() if psm.is_validated == True and psm.is_excluded == False])

    def get_residuals(self, threshold=0):
        """returns the residuals, threshold if the minimum value for residuals to be taken into acount"""
        if self.quant_residuals < threshold:
            return None
        return self.quant_residuals

    # Setters

    def set_ident_data_mzid(self, identMzid, max_rank):
        """fill PSM list from identification result dict from pyteomics"""
        for identItem in identMzid["SpectrumIdentificationItem"]:
            # Iterate over identification item and create an instance of the object psm for each
            if len(self.get_psms()) + 1 <= max_rank:
                self.psms.append(
                    Psm(rank=len(self.get_psms()) + 1, spectrum=self, identificationItem=identItem)
                )

    def set_spec_data_mgf(self, specMgf):
        "add spetrum information from a pyteomics mgf object"
        self.fragIntens: array = specMgf["intensity array"]
        self.fragMz: array = specMgf["m/z array"]
        self.precIntens: float = specMgf["params"]["pepmass"][1]
        print("------------------>", self.precIntens)
        self.precMz: float = specMgf["params"]["pepmass"][0]
        self.rt: float = specMgf["params"]["rtinseconds"]

    def set_spec_data_mzml(self, spec_mzml):
        "add spetrum information from a pyteomics mzml object"

        self.fragIntens: array = spec_mzml["intensity array"]
        self.fragMz: array = spec_mzml["m/z array"]
        self.precIntens: float = spec_mzml["total ion current"]

        self.precMz: float = spec_mzml["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0][
            "selected ion m/z"
        ]
        self.rt: float = spec_mzml["scanList"]["scan"][0]["scan start time"]

    def set_sum_intens_annot_frag(self):
        """For proteoforms in self.psms where is_validated is true, set self.sumIntensAnnotFrag as the sum of annotated peaks"""
        self.sumIntensAnnotFrag = 0
        seenPeaks = []
        for psm in self.get_psms():
            if psm.is_validated:
                ann = psm.get_annotation()

                for type in ann.values():
                    for i in range(len(type["intens"])):
                        if type["index"][i] not in seenPeaks:
                            seenPeaks.append(type["index"][i])
                            self.sumIntensAnnotFrag += type["intens"][i]

    def annotateFragPsm(self, frag_mz_tol):
        for psm in self.get_psms():
            psm.setAnnotatedFragments(frag_mz_tol=frag_mz_tol)

    def updateValidation(self):
        pass

    # --------------------------- PSMs Ratio in Spectra -------------------------- #

    def update_psms_ratio(self, psms_rank=[], verbose=False):
        """compute the ratios of psms in "psms"
        if psms_rank is empty proteoforms that are validated (is_validated == T) are considered
        else the proteforms of rank in psms_rank are considered"""

        # reset ratios:
        for psm in self.get_psms():
            psm.ratio = 0

        if len(psms_rank) == 0:
            psms = self.get_validated_psm()
        else:
            psms = [psm for idx, psm in self.get_psms() if idx in psms_rank]

        if len(psms) == 1:
            psms[0].ratio = 1
        if len(psms) > 1:
            self.__ratios_multiple(psms, verbose=verbose)
        else:
            pass

    def __ratios_multiple(self, psms, verbose):
        # Position (1 based) of modified residues across multiple psm (e.g K14acK24ac K9me3 -> [9,14,24])
        mod_pos_list = []
        for psm in psms:
            for mod in psm.proteoform.get_modification_dict():
                if mod["location"] not in mod_pos_list:
                    mod_pos_list.append(mod["location"])
        mod_pos_list = sorted(mod_pos_list)
        self.proteforms_multiquant = [psm.get_modification_brno() for psm in psms]

        # Add start and end of the sequence, sort the postions and create "interval of interest" for both n and c term frags
        mod_pos_list_n = sorted([*mod_pos_list, 1, len(psms[0].proteoform.get_peptide_sequence()) + 1])
        mod_pos_list_c = sorted([*mod_pos_list, 0, len(psms[0].proteoform.get_peptide_sequence())])
        intervals_n = [[mod_pos_list_n[p], mod_pos_list_n[p + 1] - 1] for p in range(len(mod_pos_list_n) - 1)]
        intervals_c = [[mod_pos_list_c[p] + 1, mod_pos_list_c[p + 1]] for p in range(len(mod_pos_list_c) - 1)]

        # Create matrix with modification induced mass shift (row = proteoform, col position)
        mod_mass_matrix = np.zeros((len(psms), len(mod_pos_list)))
        for row, psm in enumerate(psms):
            for mod in psm.proteoform.get_modification_dict():
                col = mod_pos_list.index(mod["location"])
                mod_mass_matrix[row, col] = mod["monoisotopicMassDelta"]

        # Compute matrix with cumulative mass shift for n and c term fragments ( row: proteoform, column:corresponding position interval)
        additive_mass_c = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1] + 1))
        for col in range(additive_mass_c.shape[1] - 1, 0, -1):
            additive_mass_c[:, col - 1] = np.add(additive_mass_c[:, col], mod_mass_matrix[:, col - 1])

        additive_mass_n = np.zeros((mod_mass_matrix.shape[0], mod_mass_matrix.shape[1] + 1))
        for col in range(1, additive_mass_n.shape[1]):
            additive_mass_n[:, col] = np.add(additive_mass_n[:, col - 1], mod_mass_matrix[:, col - 1])

        # Compute grouping matrix, give a common index to intervals with same mass across proteoforms (e.g [[54,54][54,59]] = [[1,1][1,2]])
        unique_matrix_n = self.__group_matrix(additive_mass_n)
        unique_matrix_c = self.__group_matrix(additive_mass_c)

        # reduce intervales
        red_unique_matrix_n, red_intervals_n = self.__reduce_intervals(unique_matrix_n, intervals_n)
        red_unique_matrix_c, red_intervals_c = self.__reduce_intervals(unique_matrix_c, intervals_c)

        # compute intensity matrices
        intensity_matrix_n = self.__intensity_matrix(red_unique_matrix_n, red_intervals_n, psms, "n-term")
        intensity_matrix_c = self.__intensity_matrix(red_unique_matrix_c, red_intervals_c, psms, "c-term")

        # transpose
        t_red_unique_matrix_n = red_unique_matrix_n.transpose()
        t_red_unique_matrix_c = red_unique_matrix_c.transpose()
        t_intensity_matrix_n = intensity_matrix_n.transpose()
        t_intensity_matrix_c = intensity_matrix_c.transpose()

        # Combine matrices
        unique_matrix = np.vstack((t_red_unique_matrix_c, t_red_unique_matrix_n))
        intensity_matrix = np.vstack((t_intensity_matrix_c, t_intensity_matrix_n))
        # save matrices:
        self.unique_matrix = unique_matrix
        self.intensity_matrix = intensity_matrix

        if verbose:
            print("PSMS:")
            for p in psms:
                print(p.proteoform.protein_ids)
                print(p.get_modification_proforma())
                print(p.spectrum.id)
                print(p.rank)

            print("intervals_n: \n", intervals_n)
            print("intervals_c: \n", intervals_c)

            print("additive_mass_n: \n", additive_mass_n)
            print("additive_mass_c: \n", additive_mass_c)

            print("unique_matrix_n: \n", unique_matrix_n)
            print("unique_matrix_c: \n", unique_matrix_c)

            print("red_unique_matrix_n: \n", red_unique_matrix_n)
            print("red_unique_matrix_c: \n", red_unique_matrix_c)

            print("red_intervals_n: \n", red_intervals_n)
            print("red_intervals_c: \n", red_intervals_c)

            print("red_unique_matrix_n: \n", red_unique_matrix_n)
            print("red_unique_matrix_c: \n", red_unique_matrix_c)

            print("intensity_matrix_n: \n", intensity_matrix_n)
            print("intensity_matrix_c: \n", intensity_matrix_c)

            print("unique_matrix: \n", unique_matrix)
            print("intensity_matrix: \n", intensity_matrix)

        # reduce equations that are redundant
        # unique_matrix, intensity_matrix = self.__merge_equations(unique_matrix, intensity_matrix)
        # save matrices
        self.unique_matrix_r = unique_matrix
        self.intensity_matrix_r = intensity_matrix

        # Check whether any fragment are found in pairs
        if np.all((self.intensity_matrix_r == 0)):
            for idx, psm in enumerate(psms):
                psm.ratio = 0
            if verbose:
                print("NO FRAGMENTS FOUND FOR STOCHIOMETRY DETERMINATION")
            return 0

        # get eq system and weights
        equations, variables, W = self.__equation_system(unique_matrix, intensity_matrix)

        if verbose:
            print("Weights: \n", W)

            print("unique_matrix reduced: \n", unique_matrix)
            print("intensity_matrix reduced: \n", intensity_matrix)

            print("equations: \n", equations)
            print("variables: \n", variables)

        self.equations = unique_matrix
        self.variables = intensity_matrix

        # Non negative least square solving:
        # results = nnls(equations, variables)

        # Non negative least square solving WEIGHTED:
        results = nnls(np.sqrt(W)[:, None] * equations, np.sqrt(W) * variables)

        ratios_psms = [ratio / sum(results[0]) for ratio in results[0]]  # normalized ratios

        if verbose:
            print("ratios_psms: \n", ratios_psms)

        self.quant_residuals = results[1]

        # assign rations:
        for idx, psm in enumerate(psms):
            psm.ratio = ratios_psms[idx]

    def __group_matrix(self, A):
        "From a matrix A return a matrix M with indices corresponding to identical groups column wise"
        M = np.zeros((A.shape[0], A.shape[1]))
        for col in range(A.shape[1]):
            grp_mass = []
            for row, m in enumerate(A[:, col]):
                if m not in grp_mass:
                    grp_mass.append(m)
                grp_index = grp_mass.index(m)
                M[row, col] = grp_index + 1
        return M

    # @staticmethod
    # @jit()
    def __reduce_intervals(self, unique_matrix, intervals):
        """From a unique/group matrix where columns correspond to intervals in "intervals",
        merge or delete intervals that are uninformative on the ratios of proteoforms"""
        n_col = unique_matrix.shape[1]
        n_row = unique_matrix.shape[0]

        # first pass to remove interval without any unique ions acroos proteoforms
        col_rm = []
        for col in range(unique_matrix.shape[1]):
            if (unique_matrix[:, col] == np.ones(n_row)).all():  # if no unique ions in interval
                col_rm.append(col)

        unique_matrix = np.delete(unique_matrix, col_rm, 1)
        intervals = np.delete(intervals, col_rm, 0)

        pair_found = True
        while pair_found:
            pair_found = False
            prev_col = np.zeros(n_row)
            for col in range(unique_matrix.shape[1]):
                if (unique_matrix[:, col] == prev_col).all():
                    pair_found = True

                    unique_matrix = np.delete(unique_matrix, col - 1, 1)
                    intervals[col, :] = [intervals[col - 1, 0], intervals[col, 1]]
                    intervals = np.delete(intervals, col - 1, 0)
                    break
                prev_col = unique_matrix[:, col]

        return unique_matrix, intervals.tolist()

    # @staticmethod
    # @jit(nopython=True)
    def __intensity_matrix(self, unique_matrix, intervals, psms, direction):
        intensity_matrix = np.ones((unique_matrix.shape[0], unique_matrix.shape[1]))
        for col in range(unique_matrix.shape[1]):
            for row in range(unique_matrix.shape[0]):
                fragments_at_range = psms[row].get_fragments_at_range(intervals[col], direction=direction)
                intensity_matrix[row, col] = psms[row].spectrum.get_sum_intens_fragment_list(
                    psms[row], fragments_at_range
                )

        return intensity_matrix

    # @staticmethod
    # @jit(nopython=True)
    def __merge_equations(self, A, B):
        """from matrices A and B of the same dimenssion return A_out and B_out
        where identical rows in A_out have been merged and corresponding rows in B_out have been added"""

        y, inverse = np.unique(A, axis=0, return_inverse=True)

        A_out = np.empty((0, A.shape[1]), int)
        B_out = np.empty((0, A.shape[1]), int)

        for g in np.unique(inverse):
            indices_identical = [i for i, x in enumerate(inverse) if x == g]
            A_out = np.append(
                A_out, np.array([A[indices_identical[0], :]]), axis=0
            )  # add first group row from unique matrix
            B_out = np.append(
                B_out, np.array([B[indices_identical[0], :]]), axis=0
            )  # add first group row from intensity matrix

            for i in range(1, len(indices_identical)):  # sum up intensities of identical rows
                B_out[-1, :] = np.sum([B_out[-1, :], B[indices_identical[i], :]], axis=0)

        return A_out, B_out

    # @staticmethod
    # @jit(nopython=True)
    def __equation_system(self, unique_matrix_t, intensity_matrix_t):
        var = 0
        equations = []
        variables = []
        W = []

        for row in range(unique_matrix_t.shape[0]):
            eq = unique_matrix_t[row, :]
            unique_grps = np.unique(eq).tolist()

            for u in unique_grps[:-1]:
                retaining_ratios = [v == u for v in eq]

                sum_intens_at_loc = intensity_matrix_t[row, :]

                W = np.append(W, max(intensity_matrix_t[row, :]))

                sum_intens_at_loc = np.unique(
                    sum_intens_at_loc
                ).tolist()  # TODO problem if two group have the same intensity sum (very unlikely)
                sum_intens_at_loc = sum(np.unique(sum_intens_at_loc).tolist())

                intens_proteos_at_loc = list(compress(intensity_matrix_t[row, :], retaining_ratios))
                intens_proteos_at_loc = np.unique(intens_proteos_at_loc).tolist()
                intens_proteos_at_loc = sum(intens_proteos_at_loc)

                if sum_intens_at_loc != 0:
                    var = intens_proteos_at_loc / sum_intens_at_loc
                else:
                    var = 0

                variables.append(var)
                equations.append(1 * np.array(retaining_ratios))

        # print(equations)
        # print(variables)
        # print(W)

        W = np.append(W, np.max(W))
        W = W / np.max(W)

        variables.append(1)
        equations.append(np.ones(unique_matrix_t.shape[1]))

        ### APPEND ADDITIONAL EQUATION FOR PROTEOFORM GROUPS WITH MISSING UNIQUE IONS ###

        # print("equations")
        # print(equations)

        self.miss_determining_ions = False

        is_unique = np.full(unique_matrix_t.shape[1], False)
        for row in equations:
            row_unique_indexes = []
            for x in row:
                if list(row).count(x) == 1:
                    row_unique_indexes.append(list(row).index(x))
                    self.miss_determining_ions = True

                for index in row_unique_indexes:
                    is_unique[index] = True

        # var to check whether ratio had to be added because of missing unique ions

        for i in range(len(is_unique)):
            if is_unique[i] == False:
                add_eq = np.zeros(unique_matrix_t.shape[1])
                add_eq[i] = 1
                equations.append(add_eq)
                variables.append(1 / unique_matrix_t.shape[1])
                W = np.append(W, np.mean(W))

        return equations, variables, W

    # ------------------ quantification of two proteoforms only ------------------ #

    def __ratios_pair(self, psms):
        # TODO works only for 2 validated psms

        validatedPsms = psms

        uniquePairwise = (
            []
        )  # list of list of unique ions names formatted as such [[List Unique validatedPsms[0 and 1]], [List Unique validatedPsms[1 and 2]], .... ]
        ratioPairwise = []

        if len(validatedPsms) == 1:
            validatedPsms[0].ratio = 1.0
        else:
            for p in range(0, len(validatedPsms) - 1):
                uniqueFragments = []
                for fragType in validatedPsms[p].proteoform.theoFrag.keys():
                    for fragment in validatedPsms[p].proteoform.theoFrag[fragType]:
                        if (
                            validatedPsms[p].proteoform.theoFrag[fragType][fragment]
                            != validatedPsms[p + 1].proteoform.theoFrag[fragType][fragment]
                        ):
                            uniqueFragments.append(fragment)
                uniquePairwise.append(uniqueFragments)

            for p in range(0, len(validatedPsms) - 1):
                A = self.get_sum_intens_fragment_list(validatedPsms[p], uniquePairwise[p])
                B = self.get_sum_intens_fragment_list(validatedPsms[p + 1], uniquePairwise[p])

                validatedPsms[p].ratio = A / (A + B)
                validatedPsms[p + 1].ratio = B / (A + B)

    def __get_unique_fragments(self, proteo_a, proteo_b):
        "from two proteoform object get the list of unique fragments (frag with mass shift)"
        uniqueFragments = []
        for fragType in validatedPsms[p].proteoform.theoFrag.keys():
            for fragment in validatedPsms[p].proteoform.theoFrag[fragType]:
                if (
                    validatedPsms[p].proteoform.theoFrag[fragType][fragment]
                    != validatedPsms[p + 1].proteoform.theoFrag[fragType][fragment]
                ):
                    uniqueFragments.append(fragment)

        return uniqueFragments

    def get_sum_intens_fragment_list(self, psm, fragments):
        # could be moved to psm object
        intensities = []
        for fragType in psm.annotation.values():
            for i, fragCode in enumerate(fragType["fragCode"]):
                if fragCode in fragments:
                    intensities.append(int(fragType["intens"][i]))

        if len(intensities) == 0:
            return 0
        else:
            return sum(intensities)
