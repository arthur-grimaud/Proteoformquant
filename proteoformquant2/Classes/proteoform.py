from pickle import TRUE
from pyteomics import mass
from logging import warning
import plotly.graph_objs as go
from Classes.envelope import Envelope
from Utils import constant
from statistics import mean

class Proteoform():

    def __init__(self,peptideSequence, modificationBrno, modificationProforma, modificationDict = {}):

        #Proteoform Description
        self.peptideSequence: str = peptideSequence
        self.modificationDict: dict = modificationDict
        self.modificationBrno: str = modificationBrno
        self.modificationProforma: str = modificationProforma


        #Theoretical Fragments
        self.theoFrag = None

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


    def getModificationDict(self):
        return self.modificationDict

    def getModificationBrno(self):
        return self.modificationBrno

    def getModificationProforma(self):
        return self.modificationProforma    

    def getLinkedPsm(self):
        """Return a list of linked PSMs"""
        return [psm for psm in self.linkedPsm]

    def getValidatedLinkedPsm(self):
        """Return a curated list of linked psm whose self.isvalidated = true"""
        return [psm for psm in self.linkedPsm if psm.isValidated]

    def getProteoformTotalIntens(self):
        return self.totalIntens

    def getColor(self):
        return self.color

    def getEnvelope(self):
        return self.envelope
        

    #Setters

    def setColor(self, colorInt):
        self.color = colorInt
        return self

    def linkPsm(self, psm):
        self.linkedPsm.append(psm)


    def setTheoreticalFragments(self, ionTypes):
        """ Returns and set a list of m/z of fragment ions  and informations on the type/position of each fragments for a given peptidoform/proteoform"""

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

    def computeEnvelope(self,scoreFitThreshold):
        """instanciate an envelope object by providing the list of psm associated to that proteoform"""
        if len(self.getValidatedLinkedPsm()) > 5:
            env = Envelope()
            env.computeEnvelope(self.getValidatedLinkedPsm(), scoreFitThreshold)
            print(env)
            if env.scoreFitted > scoreFitThreshold:
                self.envelope = env
            else:
                self.envelope = env # !!! 
                print("Could not fit curve to chromatogram")
        else:
            #print("Not enough datapoints to compute envelope")
            pass



    def setProteoformTotalIntens(self, method= "precursor"):
        """Return the sum of intensities of psm of that proteoform method = precursor  or annotated (correspond to the intensity value used)"""

        self.totalIntens = 0

        
        if method == "AUC":
            if self.getEnvelope() != None:
                print(self.getEnvelope().getAUC())
                self.totalIntens+=self.getEnvelope().getAUC()
            return None

        for psm in self.linkedPsm:
            if method == "precursor":
                self.totalIntens+=psm.getPrecIntensRatio()
            if method == "annotated":
                self.totalIntens+=psm.getAnnotMsmsIntensRatio()

        


    def setProteoformPsmValidation(self):
        """   """

        if self.getEnvelope() == None: #if no envelopes are found for that proteoform
            for psm in self.linkedPsm:
                psm.isValidated = False
                psm.ratio = 0.0 
        else:
            for psm in self.getEnvelope().psmsOutliers:
                print("removing aberant psm")
                psm.isValidated = False
                psm.ratio = 0.0 