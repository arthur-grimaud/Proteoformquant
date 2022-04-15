from pickle import TRUE
from pyteomics import mass
from logging import warning
import plotly.graph_objs as go
from Classes.elution_profile import ElutionProfile
from Utils import constant
from statistics import mean

class Proteoform():

    def __init__(self,peptideSequence, modificationBrno, modificationProforma, modificationDict = {}, protein_ids=[]):

        #Proteoform Description
        self.peptideSequence: str = peptideSequence
        self.modificationDict: dict = modificationDict
        self.modificationBrno: str = modificationBrno
        self.modificationProforma: str = modificationProforma

        #Protein(s) reference
        self.protein_ids = protein_ids



        #Theoretical Fragments
        self.theoFrag = None
        self.fragment_types = None #types of fragments computedS

        #All PSMs of that proteoform
        self.linkedPsm : list() = []

        #ElutionTime~intensity envelope
        self.envelope = None

        #color
        self.color = None

        #Summed intensity of linked psm
        self.totalIntens = 0  



    #Getters

    def getTheoFrag(self):
        return self.theoFrag

    def getMzFirstPsm(self):
        return self.linkedPsm[0].getCalculatedMassToCharge()

    def getTheoPrecMz(self):
        return self.linkedPsm[0].getCalculatedMassToCharge()

    def getMzMean(self):
        pass

    def get_peptide_sequence(self):
        return self.peptideSequence

    def get_modification_dict(self):
        return self.modificationDict

    def get_modification_brno(self):
        return self.modificationBrno

    def get_modification_proforma(self):
        return self.modificationProforma    

    def get_linked_psm(self):
        """Return a list of linked PSMs"""
        return [psm for psm in self.linkedPsm]

    def get_validated_linked_psm(self):
        """Return a curated list of linked psm whose self.isvalidated = true"""
        return [psm for psm in self.linkedPsm if psm.isValidated]

    def get_proteoform_total_intens(self):
        return self.totalIntens

    def get_color(self):
        return self.color

    def get_elution_profile(self):
        return self.envelope

    def get_protein_ids(self):
        return self.protein_ids     

    #Setters

    def set_color(self, colorInt):
        self.color = colorInt
        return self

    def link_psm(self, psm):
        self.linkedPsm.append(psm)


    def compute_theoretical_fragments(self, ionTypes):
        """ Returns and set a list of m/z of fragment ions  and informations on the type/position of each fragments for a given peptidoform/proteoform"""

        #store the fragment types looked for:
        self.fragment_types = ionTypes 

        ion_formulas = constant.ion_formulas 
        sequence = self.peptideSequence
        modifications = self.modificationDict

        frag_masses = {}

        #fragments masses:
        for ion_type in ionTypes:
            frag_masses_iontype =  {}

            if "-" in ion_type: #Internal Fragment
                #sum of all modification masses present in the internal fragment
                sum_mods = lambda modifications, i, j : sum( [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"] <= j  ] ) #sum mods delta for internal fragm ions
                #get all sub string of the peptide sequence
                sub_sequences = [(sequence[i-1:j],i,j,ion_type, [mod["location"] for mod in modifications if i <= mod["location"] <= j  ] ) for i in range(1,len(sequence)) for j in range(i + 1, len(sequence))]
                #compute internal frag masses
                frag_masses_iontype.update({ ','.join(str(s) for s in seq[1:4]): round(mass.fast_mass(sequence=seq[0], ion_type=ion_type, ion_comp= ion_formulas[ion_type]) + sum_mods(modifications, seq[1], seq[2]),4) for seq in sub_sequences })


            else: #Terminal Fragment
                if any(i in ion_type for i in ["a","b","c"]): #Nterm
                    sum_mods = lambda modifications,i,j : sum( [mod["monoisotopicMassDelta"] for mod in modifications if  mod["location"] <= j] )
                    sub_sequences = [(sequence[:j],1,j,ion_type, [mod["location"] for mod in modifications if  mod["location"] <= j] ) for j in range(2,len(sequence))]
                    frag_masses_iontype.update(  { ','.join(str(s) for s in seq[1:4]): round(mass.fast_mass(sequence=seq[0], ion_type=ion_type, ion_comp= ion_formulas[ion_type]) + sum_mods(modifications, seq[1], seq[2]),4) for seq in sub_sequences})
                
                else: #Cterm
                    sum_mods= lambda modifications,i,j : sum( [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"] ] )
                    sub_sequences= [(sequence[i-1:],i,len(sequence),ion_type, [mod["location"] for mod in modifications if i <= mod["location"] ] ) for i in range(1,len(sequence)+1)]
                    frag_masses_iontype.update( { ','.join(str(s) for s in seq[1:4]): round(mass.fast_mass(sequence=seq[0], ion_type=ion_type, ion_comp= ion_formulas[ion_type]) + sum_mods(modifications, seq[1], seq[2]),4) for seq in sub_sequences } )

            frag_masses[ion_type] = frag_masses_iontype

        
        if self.theoFrag == None:
            self.theoFrag = frag_masses
        else:
            warning("Theoretical fragments already set for proteoform: " + self.modificationBrno + "OVERWRITING !")
            self.theoFrag = frag_masses
        
        return(frag_masses)

    def model_elution_profile(self,elution_profile_score_threshold):
        """instanciate an envelope object by providing the list of psm associated to that proteoform"""
        if len(self.get_validated_linked_psm()) > 5:
            env = ElutionProfile()
            env.model_elution_profile(self.get_validated_linked_psm(), elution_profile_score_threshold)

            if env.score_fitted > elution_profile_score_threshold:
                self.envelope = env
            else:
                #self.envelope = env # !!! 
                print("Could not fit curve to chromatogram")
        else:
            #print("Not enough datapoints to compute envelope")
            pass



    def update_proteoform_total_intens(self, method= "precursor"):
        """Return the sum of intensities of psm of that proteoform method = precursor  or annotated (correspond to the intensity value used)"""

        self.totalIntens = 0

        
        if method == "AUC":
            if self.get_elution_profile() != None:
                #print(self.get_elution_profile().get_auc())
                self.totalIntens+=self.get_elution_profile().get_auc()
            return None

        for psm in self.get_validated_linked_psm():
            if method == "precursor":
                self.totalIntens+=psm.get_prec_intens_ratio()
            if method == "annotated":
                self.totalIntens+=psm.getAnnotMsmsIntensRatio()

        


    def update_proteoform_psm_validation(self):
        """   """

        if self.get_elution_profile() == None: #if no envelopes are found for that proteoform
            for psm in self.linkedPsm:
                psm.isValidated = False
                psm.ratio = 0.0 
        else:
            for psm in self.get_elution_profile().psms_outliers:
                #print("removing aberant psm")
                psm.isValidated = False
                psm.ratio = 0.0 