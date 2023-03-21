import unimod_mapper

um = unimod_mapper.UnimodMapper()
import warnings
from importlib.metadata import version

try:  # local modules
    from proteoformquant.Utils.misc import truncate
    import spectrum_utils.spectrum as sus
    from proteoformquant.Utils.constant import ion_direction
    from proteoformquant.Utils import constant
except ImportError:
    from Utils.misc import truncate
    import spectrum_utils.spectrum as sus
    from Utils.constant import ion_direction
    from Utils import constant


class Psm:
    def __init__(self, rank, spectrum, identificationItem):
        version("unimod_mapper")

        self.Modification = []

        # PSM characteristics
        for key in identificationItem:
            setattr(self, key, identificationItem[key])

            # print(key, identificationItem[key])

        self.rank: int = rank
        self.spectrum = spectrum
        self.proteoform = None

        self.annotation = {}

        # PSM quantification
        if rank == 1:
            self.is_validated = True
            self.ratio: float = 1.0
        else:
            self.is_validated = False
            self.ratio: float = 0.0

        # PSM exculsion
        self.is_excluded = False

        # pprint.pprint(vars(self))
        pass

    # Getters

    def is_included(self):
        return not self.is_excluded

    def get_rank(self):
        return self.rank

    def get_accessions(self):
        """return protein/accesion reference list"""
        if hasattr(self, "PeptideEvidenceRef"):
            if "accession" in self.PeptideEvidenceRef[0]:
                return [ref["accession"] for ref in self.PeptideEvidenceRef]
            if "DBSequence_Ref" in self.PeptideEvidenceRef[0]:
                return [ref["DBSequence_Ref"] for ref in self.PeptideEvidenceRef]

        return []

    def getPeptideSequence(self):
        return self.PeptideSequence

    def getChargeState(self):
        return self.chargeState

    def getCalculatedMassToCharge(self):
        return self.calculatedMassToCharge

    def get_modification_brno(self):
        """Returns modifications of a peptide in brno nomenclature: K27meK34ac"""
        deltaMods = constant.delta_mod

        if self.Modification == []:
            return "Unmodified"
        mod_list = []

        for mod in self.Modification:
            m = str(truncate(mod["monoisotopicMassDelta"], 2))
            mod_type = [value for key, value in deltaMods.items() if key.lower() in [m]]

            if mod_type != []:
                mod_type = mod_type[0]
                mod_list.append(
                    str(self.PeptideSequence[int(mod["location"]) - 1]) + str(mod["location"]) + mod_type
                )
            else:
                # print("Delta mass to brno conversion failed, {0} not found".format(mod["monoisotopicMassDelta"]))
                mod_list.append(
                    str(self.PeptideSequence[int(mod["location"]) - 1]) + str(mod["location"]) + "na"
                )

        return "".join(mod_list)

    def getModificationsSu(self):
        """Returns a dictionnary of modifications in the spectrum utils format: {'position(0-based)': 'mass shift'}"""
        modDict = {}
        for mod in self.proteoform.get_modification_dict():
            modDict[int(mod["location"])] = float(mod["monoisotopicMassDelta"])  # use list comprehension ??
        return modDict

    def get_modification_proforma(self, mod_mass_to_name=None):
        """Returns the proteoform modification in the proforma format"""
        with warnings.catch_warnings():  # catch warnings from unimod
            if self.Modification == []:
                proforma = self.PeptideSequence
            else:
                sequence_list = list(self.PeptideSequence.strip(" "))

                for mod in reversed(self.Modification):
                    modMass = mod["monoisotopicMassDelta"]

                    if mod_mass_to_name is None:  # if called without mod to mass dict
                        modName = um.id_to_name(um.mass_to_ids(modMass, decimals=2)[0])[0]
                    elif (
                        modMass not in mod_mass_to_name
                    ):  # if called with mod to mass dict but mass not added
                        modName = um.id_to_name(um.mass_to_ids(modMass, decimals=2)[0])[0]
                        mod_mass_to_name[modMass] = modName  # mass found in mass
                    else:
                        modName = mod_mass_to_name[modMass]

                    sequence_list.insert(int(mod["location"]), "[{0}]".format(modName))

                proforma = "".join(sequence_list)

            return proforma + "/" + str(self.chargeState)

        # mods = um.mass_to_ids(79.96, decimals=1)

    def get_annotation(self):
        return self.annotation

    def get_fragment_coverage(self):
        tot_frag_count = 0
        frag_count = 0
        for dir in ["n-term", "c-term"]:
            fragments_in_direction = [
                self.get_intensity_at_pos(pos, dir) for pos in range(1, len(self.proteoform.peptideSequence))
            ]

            frag_count += sum(x is not None for x in fragments_in_direction)
            tot_frag_count += len(fragments_in_direction)

        return frag_count / tot_frag_count

    def get_annotation_pair_format(self, first_var, second_var):
        "return annotated fragments in the following format: [(fragcode,intensity), ... ]" ""

        frag_code_list = []
        intens_list = []

        for frag_type in self.annotation.values():
            frag_code_list = frag_code_list + frag_type[first_var]
            intens_list = intens_list + frag_type[second_var]

        return list(zip(frag_code_list, intens_list))

    def get_prec_intens_ratio(self):
        """given self.ratio return the corresponding precursor intensity fraction for that psm"""
        return self.spectrum.get_prec_intens() * self.ratio

    def getAnnotMsmsIntensRatio(self):
        """given self.ratio return the corresponding annotated fragment intensity sum fraction for that psm"""
        return self.spectrum.get_sum_intens_annot_frag() * self.ratio

    def get_fragments_at_range(self, bounds, direction):
        """Get a list of fragments names given a range of residue position (1-based) to be present in the fragments"""
        ## TODO DOES NOT WORK FOR INTERNAL IONS
        fragments = []

        theo_frags = self.proteoform.theoFrag
        # TO FINISH

        # print(bounds)

        for frag_type_name, frag_type in theo_frags.items():
            # print(ion_direction[frag_type_name])
            if ion_direction[frag_type_name] == direction:
                # print("yes")
                for frag_code in frag_type.keys():
                    frag_code_list = frag_code.split(",")
                    # print(frag_code)
                    if direction == "n-term":
                        i = 1
                    if direction == "c-term":
                        i = 0

                    if bounds[0] <= int(frag_code_list[i]) and int(frag_code_list[i]) <= bounds[1]:
                        fragments.append(frag_code)
                        # print("nterm ion" + frag_code)

        return fragments

    def get_intensity_at_pos(self, pos, direction):
        """Get the intensity of the annotated fragment at position "pos",
        where pos define the amino acid index (1 based) to the left (n-term)
        of the fragmentation site (e.g pos=3 get the intens of a c3 ion and
        a z7 ion for an peptide of 10 AA).
        if direction = both the sum of z and c ions is returned"""

        intensity = 0

        for frag_type_name, frag_type in self.annotation.items():
            direction_frag_type = ion_direction[frag_type_name]

            if direction_frag_type == direction or direction == "both":
                if direction_frag_type == "n-term":
                    try:
                        index = frag_type["pos"].index(str("1:" + str(pos)))
                        intensity += frag_type["intens"][index]
                    except ValueError:
                        pass

                if direction_frag_type == "c-term":
                    try:
                        # print(str( str(pos) + str(len(self.proteoform.peptideSequence)) ))
                        index = frag_type["pos"].index(
                            str(str(pos + 1) + ":" + str(len(self.proteoform.peptideSequence)))
                        )
                        intensity += frag_type["intens"][index]
                    except ValueError:
                        pass
        if intensity == 0:
            return None

        return intensity

    # Setters

    def exclude(self):
        "set psm.is_excluded to false, excluded psm are not returned by spectrum.get_psms()"
        self.is_excluded = True

    def setProteoform(self, proteoform):
        self.proteoform = proteoform

    def setAnnotatedFragments(self, frag_mz_tol, maxCharge=1):
        # get information to create a spectrum utils MsmsSpectrum object
        id = self.spectrum.get_id()
        fragIntens = self.spectrum.get_frag_intens()
        fragMz = self.spectrum.get_frag_mz()

        calculatedMassToCharge = self.getCalculatedMassToCharge()
        chargeState = self.getChargeState()
        peptideSequence = self.getPeptideSequence()
        modifications = self.getModificationsSu()
        mods_brno = self.get_modification_brno()

        nAnnotFrag = 0

        # annotate each fragment type individualy

        for theoFragType in self.proteoform.getTheoFrag():
            # create msmsspectrum object
            spectrumSu = sus.MsmsSpectrum(
                id,
                calculatedMassToCharge,  # might cause an issue when charge state != 0
                chargeState,
                fragMz,
                fragIntens,
                peptide=peptideSequence,
                modifications=modifications,
            )

            # annotate fragments in msmsspectrum object
            for frag in self.proteoform.getTheoFrag()[theoFragType].items():  # TODO change this format?
                try:
                    spectrumSu = spectrumSu.annotate_mz_fragment(
                        fragment_mz=float(frag[1]),
                        fragment_charge=maxCharge,
                        fragment_tol_mass=frag_mz_tol,
                        fragment_tol_mode="Da",
                        text=frag[0],
                    )
                except (
                    ValueError
                ):  # necessary as the method .annotate_mz_fragment throw an error if the mz is not found
                    pass

            if spectrumSu.annotation is not None:
                mzTheoList, intensList, posList, indexList, mzErrorList = (
                    [],
                    [],
                    [],
                    [],
                    [],
                )

                indexInPeakList = 0
                for mz, intensity, annotation in zip(
                    spectrumSu.mz, spectrumSu.intensity, spectrumSu._annotation
                ):
                    if annotation != None:  # if the fragment has been matched/annotated
                        annotationList = str(annotation.annotation).split(sep=",")
                        mzTheoList.append(mz)
                        intensList.append(intensity)
                        posList.append(str(annotationList[0]) + ":" + str(annotationList[1]))
                        indexList.append(indexInPeakList)
                        mzErrorList.append(self.spectrum.get_frag_mz()[indexInPeakList] - mz)
                        nAnnotFrag += 1

                    indexInPeakList += 1

                # print(str(theoFragType))
                # store annotation in psm object

                self.annotation[str(theoFragType)] = {
                    "mzTheo": mzTheoList,
                    "intens": intensList,
                    "pos": posList,
                    "index": indexList,
                    "mzErrorList": mzErrorList,
                    "fragCode": [
                        ",".join((pos.split(":")[0], pos.split(":")[1], theoFragType)) for pos in posList
                    ],
                }
