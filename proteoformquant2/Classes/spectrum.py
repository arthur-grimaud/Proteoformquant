from array import array
from Classes.psm import Psm

import pprint

class Spectrum():

    def __init__(self, spectrumID, identMzid = None):

     
        self.id = spectrumID  #Unique ID for the spectrum
        self.psms : list(Psm) = []  #list of PSM for that spectrum
        
        if(identMzid != None):
            self.setIdentDataMzid(identMzid)
            
        pass
    

    #Getters:

    def getId(self):
        return self.id

    def getFragIntens(self):
        return self.fragIntens

    def getFragMz(self):
        return self.fragMz

    #Setters


    def setIdentDataMzid(self, identMzid):
        """fill PSM list from identification result dict from pyteomics"""
        for identItem in identMzid["SpectrumIdentificationItem"]: #Iterate over identification item and create an instance of the object psm for each
            self.psms.append(Psm(rank = len(self.psms)+1, spectrum=self, identificationItem = identItem))


    def setSpecDataMgf(self, specMgf):
        "add spetrum information from a pyteomics mgf object"
        self.fragIntens: array = specMgf["intensity array"]
        self.fragMz: array = specMgf["m/z array"]
        self.precIntens: float = specMgf["params"]["pepmass"][0]
        self.precMz: float = specMgf["params"]["pepmass"][1]
        self.rt: float = specMgf["params"]["rtinseconds"]


        
    def setSpecDataMzxml(self, specMgf):
        pass

    def annotateFragPsm(self):
        for psm in self.psms:
            psm.setAnnotatedFragments()
        
    def updateValidation(self):
        pass

    def updateRatio(self):
        pass

    