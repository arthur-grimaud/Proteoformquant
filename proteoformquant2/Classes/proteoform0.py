from pickle import TRUE
from pyteomics import mass
from logging import warning
import plotly.graph_objs as go
from Classes.envelope import Envelope
from Utils import constant

class Proteoform0():

    def __init__(self):

        #Proteoform found in spectra
        self.linkedSpectra : list() = []

        self.totalIntens = 0  

    #Getters

    def getProteoformTotalIntens(self):
        return self.totalIntens

    #Setters

    def linkSpectrum(self, sprectrum):
        self.linkedSpectra.append(sprectrum)

    def setProteoformTotalIntens(self, method= "precursor"):
        """Return the sum of intensities of psm of that proteoform method = precursor  or annotated (correspond to the intensity value used)"""

        self.totalIntens = 0

        for spectrum in self.linkedSpectra:
            if method == "precursor":
                self.totalIntens+=spectrum.getPrecIntensRatio()

    def setProteoformPsmValidation(self):
        """   """
        pass