### Import ###
from pyteomics import mzid 
from pyteomics import mgf
from chart_studio.plotly import plot, iplot
import plotly.express as px
import plotly.graph_objs as go
import numpy as np
import pandas as pd

#Custom classes
from Classes.spectrum import Spectrum
from Classes.proteoform import Proteoform
from Classes.proteoform0 import Proteoform0

import pprint
class Msrun():

    def __init__(self, runId:str = "Default run ID", dbse:str = "comet" ):
        
        self.runId = runId
        self.dbse = dbse
        self.identFn: str = "Not Specified" 
        self.spectraFn: str = "Not Specified"

        self.spectra: dict(Spectrum) = {} 
        self.proteoforms: dict(Proteoform) = {}
        self.proteoform0 = Proteoform0()
 

    def readMzid(self, identFn):
        """Read a spectra identification file in .mzIdenMl whose path is specfieid in self.inputFn"""
        print("start reading mzid\n")
        self.identFn = identFn #Store File that has been read
        mzidObj = mzid.read(identFn) #Create a pyteomics' mzid iterator

        for identMzid in mzidObj: #Iterate over spectra and create Spectrum object for each 
            print(identMzid["spectrumID"])
            self.spectra[identMzid["spectrumID"]] = Spectrum(spectrumID= identMzid["spectrumID"], identMzid = identMzid)
        pass

    def addMgfData(self, spectraFn):
        """Add info from mgf file to spectrum objects in self.spectra"""
        print("start reading mgf\n")
        self.spectraFn = spectraFn #Store File that has been read
        mgfObj = mgf.read(spectraFn)

        for specID in self.spectra:
            
            if self.dbse == "mascot":
                index = str(int(specID.split("=")[1]) + 1)
            if self.dbse == "comet":
                index = specID.split("=")[1]

            specMgf = mgfObj.get_spectrum(index) #need to be splited 
            self.spectra[specID].setSpecDataMgf(specMgf)
        pass

    def addMzmlData(self):
        """Add info from mzml file to spectrum objects in self.spectra add Proteform objects to self.Proteoforms"""
        pass
    
    def addProteoforms(self):
        """From spectrum objects in self.spectra add proteoforms object to self.proteoforms"""
        print("start adding proteoforms\n")
        i=0
        for specID in self.spectra:
            i+=1
            for psm in self.spectra[specID].psms:

                brno = psm.getModificationsBrno()
                seq = psm.getPeptideSequence()
                brnoSeq = brno+"-"+seq #TODO use proforma

                if brnoSeq not in self.proteoforms.keys(): #if a proteoform is new create a instance of Proteoform for this proteoform
                    self.proteoforms[brnoSeq] = Proteoform(peptideSequence=seq, modificationBrno=brno, modificationDict=psm.Modification).setColor(i)  

                self.proteoforms[brnoSeq].linkPsm(psm) #add link to Psms in Proteoform
                psm.setProteoform(self.proteoforms[brnoSeq]) #add link to Proteoform in Psm

        #Add proteoform object for unassigned proteoform
        

            # pprint.pprint(vars(self.proteoforms[brno+"-"+psm.PeptideSequence]))
            #print(self.proteoforms[brno+"-"+psm.getPeptideSequence()].linkedPsm)


    def matchFragments(self, msmsTol= 0.02, internal = False):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        print("start matching fragments\n")
        for proteoID in self.proteoforms:
            self.proteoforms[proteoID].setTheoreticalFragments(["c","zdot","c-1","z+1","z+2"])
            #print(self.proteoforms[proteoID].theoFrag)
        for spectrumID in self.spectra:
            self.spectra[spectrumID].annotateFragPsm()
            self.spectra[spectrumID].setSumIntensAnnotFrag()

        pass

    
    def updateProteoformsEnvelope(self):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        print("start updating proteoforms envelopes\n")
        for proteoID in self.proteoforms:
            self.proteoforms[proteoID].computeEnvelope()


        pass

    def updateProteoformsTotalIntens(self):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        print("start updating proteoforms envelopes\n")
        for proteoID in self.proteoforms:
            self.proteoforms[proteoID].setProteoformTotalIntens()


        pass


    def updateProteoformsValidation(self):
        """Update psm.isValdiated and add spectrum without any validated psm to proteoform[unasigned]"""
        for proteoID in self.proteoforms:
            self.proteoforms[proteoID].setProteoformPsmValidation()

        pass

    def updateUnassignedSpectra(self):
        for spectrumID in self.spectra:
            if self.spectra[spectrumID].getNumberValidatedPsm() == 0:
                self.proteoform0.linkSpectrum(self.spectra[spectrumID])



    #Visualization
    #Methods here should return 
