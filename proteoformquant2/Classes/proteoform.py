from pyteomics import mass
from logging import warning

from Classes.envelope import Envelope
from Utils import constant

class Proteoform():

    def __init__(self,peptideSequence, modificationBrno, modificationDict = {}):

        #Proteoform Description
        self.peptideSequence: str = peptideSequence

        self.modificationDict: dict = modificationDict
        self.modificationBrno: str = modificationBrno
        self.modificationProforma: str = None

        #Theoretical Fragments
        self.theoFrag = None

        #Proteoform found in spectra
        self.linkedPsm : list() = []

        #ElutionTime~intensity envelope
        self.envelopes : list(Envelope) = []

    #Getters

    def getTheoFrag(self):
        return self.theoFrag


    def getModificationDict(self):
        return self.modificationDict


    #Setters

    def linkPsm(self, psm):
        self.linkedPsm.append(psm)


    def setTheoreticalFragments(self, ionTypes):
        """ Returns and set a list of m/z of fragment ions  and informations on the type/position of each fragments for a given peptidoform/proteoform"""

        intern_ion_formulas = constant.intern_ion_formulas 
        sequence = self.peptideSequence
        modifications = self.modificationDict


        frag_masses = {}

        #fragments masses:
        for ion_type in ionTypes:
            frag_masses_iontype =  {}

            if "I" in ion_type: #Internal Fragment
                #sum of all modification masses present in the internal fragment
                sum_mods = lambda modifications, i, j : sum( [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"] <= j  ] ) #sum mods delta for internal fragm ions
                #get all sub string of the peptide sequence
                sub_sequences = [(sequence[i-1:j],i,j,ion_type, [mod["location"] for mod in modifications if i <= mod["location"] <= j  ] ) for i in range(1,len(sequence)) for j in range(i + 1, len(sequence))]
                #compute internal frag masses
                frag_masses_iontype.update({ ','.join(str(s) for s in seq[1:4]): round(mass.fast_mass(sequence=seq[0], ion_type=ion_type, ion_comp= intern_ion_formulas[ion_type]) + sum_mods(modifications, seq[1], seq[2]),4) for seq in sub_sequences })
            
            else: #Terminal Fragment
                if any(i in ion_type for i in ["a","b","c"]): #Nterm
                    sum_mods = lambda modifications,i,j : sum( [mod["monoisotopicMassDelta"] for mod in modifications if  mod["location"] <= j] )
                    sub_sequences = [(sequence[:j],1,j,ion_type, [mod["location"] for mod in modifications if  mod["location"] <= j] ) for j in range(2,len(sequence))]
                    frag_masses_iontype.update(  { ','.join(str(s) for s in seq[1:4]): round(mass.fast_mass(sequence=seq[0], ion_type=ion_type, ion_comp= intern_ion_formulas[ion_type]) + sum_mods(modifications, seq[1], seq[2]),4) for seq in sub_sequences})
                
                else: #Cterm
                    sum_mods= lambda modifications,i,j : sum( [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"] ] )
                    sub_sequences= [(sequence[i-1:],i,len(sequence),ion_type, [mod["location"] for mod in modifications if i <= mod["location"] ] ) for i in range(1,len(sequence)+1)]
                    frag_masses_iontype.update( { ','.join(str(s) for s in seq[1:4]): round(mass.fast_mass(sequence=seq[0], ion_type=ion_type, ion_comp= intern_ion_formulas[ion_type]) + sum_mods(modifications, seq[1], seq[2]),4) for seq in sub_sequences } )

            frag_masses[ion_type] = frag_masses_iontype


        if self.theoFrag == None:
            self.theoFrag = frag_masses
        else:
            warning("Theoretical fragments already set for proteoform: " + self.modificationBrno + "OVERWRITING !")
            self.theoFrag = frag_masses

        return(frag_masses)

