from pickle import TRUE
from pyteomics import mass
from logging import warning
import plotly.graph_objs as go
from Classes.elution_profile import ElutionProfile
from Utils import constant

class Proteoform0():

    def __init__(self):

        #Proteoform found in spectra
        self.linkedSpectra : list() = []

        self.totalIntens = 0  

    #Getters

    def get_proteoform_total_intens(self):
        return self.totalIntens

    #Setters

    def linkSpectrum(self, sprectrum):
        self.linkedSpectra.append(sprectrum)

    def update_proteoform_total_intens(self, method= "precursor"):
        """Return the sum of intensities of psm of that proteoform method = precursor  or annotated (correspond to the intensity value used)"""

        self.totalIntens = 0

        for spectrum in self.linkedSpectra:
            if method == "precursor":
                self.totalIntens+=spectrum.getPrecIntens()

    def update_proteoform_psm_validation(self):
        """   """
        pass