from array import array
from Classes.psm import Psm

import pprint

class Spectrum():

    def __init__(self, spectrumID, identMzid = None):

     
        self.id = spectrumID  #Unique ID for the spectrum
        self.psms : list(Psm) = []  #list of PSM for that spectrum
        
        if(identMzid != None):
            self.setIdentDataMzid(identMzid)

        self.sumIntensAnnotFrag = 0
            
        pass
    

    #Getters:

    def getId(self):
        return self.id

    def getPrecMz(self):
        return self.precMz
    
    def getPrecIntens(self):
        return self.precIntens

    def getRt(self):
        return self.rt

    def getFragIntens(self):
        return self.fragIntens

    def getFragMz(self):
        return self.fragMz
    
    def getSumIntensAnnotFrag(self):
        return self.sumIntensAnnotFrag


    #Setters


    def setIdentDataMzid(self, identMzid):
        """fill PSM list from identification result dict from pyteomics"""
        for identItem in identMzid["SpectrumIdentificationItem"]: #Iterate over identification item and create an instance of the object psm for each
            self.psms.append(Psm(rank = len(self.psms)+1, spectrum=self, identificationItem = identItem))


    def setSpecDataMgf(self, specMgf):
        "add spetrum information from a pyteomics mgf object"
        self.fragIntens: array = specMgf["intensity array"]
        self.fragMz: array = specMgf["m/z array"]
        self.precIntens: float = specMgf["params"]["pepmass"][1]
        self.precMz: float = specMgf["params"]["pepmass"][0]
        self.rt: float = specMgf["params"]["rtinseconds"]
        
    def setSpecDataMzxml(self, specMgf):
        pass

    def setSumIntensAnnotFrag(self):
        """For proteoforms in self.psms where isvalidated is true, set self.sumIntensAnnotFrag as the sum of annotated peaks"""
        self.sumIntensAnnotFrag = 0
        seenPeaks=[]
        for psm in self.psms:
            if psm.isValidated:
                ann = psm.getAnnotation()
                for type in ann.values():
                    for i in range(len(type["intens"])):
                        if type["index"] not in seenPeaks:
                            seenPeaks.append(i)
                            self.sumIntensAnnotFrag += type["intens"][i]
        #print("sum intensities annotated: " + str(self.sumIntensAnnotFrag))



    #other methods:

    def annotateFragPsm(self):
        for psm in self.psms:
            psm.setAnnotatedFragments()
            
    def updateValidation(self):
        pass

    def updateRatio(self):
        pass

    