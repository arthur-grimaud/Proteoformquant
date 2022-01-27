from pickle import TRUE
from pyteomics import mass
from logging import warning
import plotly.graph_objs as go
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

    def getMzFirstPsm(self):
        return self.linkedPsm[0].spectrum.getPrecMz()

    def getModificationDict(self):
        return self.modificationDict

    def getModificationBrno(self):
        return self.modificationBrno


    def getValidatedLinkedPsm(self):
        """Return a curated list of linked psm whose self.isvalidated = true"""
        return [psm for psm in self.linkedPsm if psm.isValidated]


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

    def computeEnvelope(self):
        """instanciate an envelope object by providing the list of psm associated to that proteoform"""
        if len(self.getValidatedLinkedPsm()) > 5:
            env = Envelope(self.getValidatedLinkedPsm())
            self.envelopes.append(env)
            #self.getEnvelopePlot()
        else:
            #print("Not enough datapoints to compute envelope")
            pass



    #visualization:


    def getEnvelopePlot(self):

        xData = [psm.spectrum.getRt() for psm in self.getValidatedLinkedPsm()]
        yData = [psm.getPrecIntensRatio() for psm in self.getValidatedLinkedPsm()]

        

        fig = go.Figure()
        fig.add_scatter( x=xData, y=yData, mode='markers', marker=dict(size=10, color="black"), name='Precursor Intensity' )

        for env in self.envelopes: #if envelope has been computed add line to the plot

            xDataEnv = list(range(int(min(xData)),int(max(xData)),1))

            yDataEnvEstim, parametersEstim = list(env.getEnvelopeSerie(xDataEnv, method = "estimated"))
            if yDataEnvEstim[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvEstim, mode='lines', marker=dict(size=4, color="orange"), name='Estimated Parameters', line_shape='spline' )

            yDataEnvFitted, parametersFitted = list(env.getEnvelopeSerie(xDataEnv, method = "fitted"))
            #print(yDataEnvFitted)
            if yDataEnvFitted[0] != None: fig.add_scatter( x=xDataEnv, y=yDataEnvFitted, mode='lines', marker=dict(size=4, color="red"), name='Fitted Parameters', line_shape='spline' )

        titleText = "Proteoform: {0} <br> Parameters Estimated: {1} / KS: {3} <br> Parameters Fitted: {2} / KS: {4}".format(self.modificationBrno, parametersEstim , parametersFitted, env.KsEstimated, env.KsFitted)

        fig.update_layout(title=go.layout.Title(text=titleText, font=dict(
                family="Courier New, monospace",
                size=15,
            )))

        fig.show()




