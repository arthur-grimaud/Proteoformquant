### Import ###
from pyteomics import mzid 
from pyteomics import mgf
from chart_studio.plotly import plot, iplot
import plotly.express as px
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from statistics import mean
from statistics import median
from progress.bar import Bar

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

        self.precMzTolerance = 1.5
        self.intensityThreshold = 200000
        self.scoreFitThreshold = 0.7


    # ---------------------------------- getters --------------------------------- #

    def getRtRange(self):
        "Get min and max Rt in self.spectra"
        mn = min([spectrum.get_rt() for spectrum in self.spectra.values()])
        mx = max([spectrum.get_rt() for spectrum in self.spectra.values()])
        return (mn, mx)

    def getMzRange(self):
        "Get min and max precursor mz in self.spectra"
        mn = min([spectrum.getPrecMz() for spectrum in self.spectra.values()])
        mx = max([spectrum.getPrecMz() for spectrum in self.spectra.values()])
        return (mn, mx)

    def getDatasetMetrics(self):
        """returns a set of metric as a dictionnary for the ms run"""

        return {
            "TotalSpectra": len(self.spectra),
            "AssignedSpectra(method1)": len([spectrum for spectrum in self.spectra.values() if spectrum.getNumberValidatedPsm() > 0]),
            "AssignedSpectra(method2)": 0,
            "UnassignedSpectra": len([spectrum for spectrum in self.proteoform0.linkedSpectra]),
            "ChimericSpectra": len([spectrum for spectrum in self.spectra.values() if spectrum.getNumberValidatedPsm() > 1]),
            "TotalProteoforms": len(self.proteoforms)-1,
            "AvgEnvelopeScore": mean([proteoform.getEnvelope().score_fitted for proteoform in self.proteoforms.values() if proteoform.getEnvelope() != None ]),
            "MinEnvelopeScore": min([proteoform.getEnvelope().score_fitted for proteoform in self.proteoforms.values() if proteoform.getEnvelope() != None ]),
            "MedianEnvelopeS":   median([proteoform.getEnvelope().param_fitted[1] for proteoform in self.proteoforms.values() if proteoform.getEnvelope() != None ]),
            "MedianEnvelopeA":   median([proteoform.getEnvelope().param_fitted[2] for proteoform in self.proteoforms.values() if proteoform.getEnvelope() != None ]),
            "MedianEnvelopeK":   median([proteoform.getEnvelope().param_fitted[3] for proteoform in self.proteoforms.values() if proteoform.getEnvelope() != None ])
        }

    # ----------------------------------- main ----------------------------------- #
 
    def read_mzid(self, identFn):
        """Read a spectra identification file in .mzIdenMl whose path is specfieid in self.inputFn"""
        self.identFn = identFn #Store File that has been read
        mzidObj = mzid.read(identFn) #Create a pyteomics' mzid iterator

        with Bar('loading identifications', max=1) as bar:

            for identMzid in mzidObj: #Iterate over spectra and create Spectrum object for each 
                self.spectra[identMzid["spectrumID"]] = Spectrum(spectrumID= identMzid["spectrumID"], identMzid = identMzid)
                bar.next()
        pass

    def add_mgf_data(self, spectraFn):
        """Add info from mgf file to spectrum objects in self.spectra"""

        self.spectraFn = spectraFn #Store File that has been read
        mgfObj = mgf.read(spectraFn)


        with Bar('loading spectra', max=1) as bar:

            for specID in self.spectra: #Take into account how DBSEs store spectra ids 
                if self.dbse == "mascot":
                    index = str(int(specID.split("=")[1]) + 1)
                if self.dbse == "comet":
                    index = specID.split("=")[1]

                specMgf = mgfObj.get_spectrum(index) #need to be splited 
                self.spectra[specID].setSpecDataMgf(specMgf)
                bar.next()
        pass

    def addMzmlData(self):
        """Add info from mzml file to spectrum objects in self.spectra add Proteform objects to self.Proteoforms"""
        pass
    
    def add_proteoforms(self):
        """From spectrum objects in self.spectra instanciate proteoforms object to self.proteoforms"""
        i=0
        with Bar('Adding PSMs to Proteoform objects', max=1) as bar:

            for specID in self.spectra:
                i+=1
                for psm in self.spectra[specID].psms:

                    proforma = psm.getModificationProforma()
                    brno = psm.getModificationBrno()
                    seq = psm.getPeptideSequence()
                    brnoSeq = brno+"-"+seq #TODO use proforma

                    if proforma not in self.proteoforms.keys(): #if a proteoform is new create a instance of Proteoform for this proteoform
                        self.proteoforms[proforma] = Proteoform(peptideSequence=seq, modificationBrno=brno, modificationProforma=proforma, modificationDict=psm.Modification).setColor(i)  

                    self.proteoforms[proforma].linkPsm(psm) #add link to Psms in Proteoform
                    psm.setProteoform(self.proteoforms[proforma]) #add link to Proteoform in Psm
                    bar.next()



    def match_fragments(self, msmsTol= 0.02, internal = False):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""

        with Bar('Generating theoretical fragments', max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].setTheoreticalFragments(["c","zdot","z+1"])
                bar.next()

        with Bar('Matching fragments', max=1) as bar:
            for spectrumID in self.spectra:
                self.spectra[spectrumID].annotateFragPsm()
                self.spectra[spectrumID].setSumIntensAnnotFrag()
                bar.next()

        pass

    
    def updateProteoformsEnvelope(self):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        with Bar('Updating proteoforms envelopes', max=1) as bar:

            #TODO add a function thjat set the bounds based on the entire set of envelopes  
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].model_elution_profile(self.scoreFitThreshold)
                bar.next()
        pass


    def updateProteoformsTotalIntens(self):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        with Bar('Updating Proteoforms Total Intens', max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].setProteoformTotalIntens()
                bar.next()
        pass

    def updatePsmValidation(self):
        """Update psm.isValidated"""
        with Bar('Updating PSM Validation', max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].setProteoformPsmValidation()
        pass

    def updateUnassignedSpectra(self):
        """add spectra without any validated psm to "proteoform0" in self.proteoforms"""
        for spectrumID in self.spectra:
            if self.spectra[spectrumID].getNumberValidatedPsm() == 0:
                self.proteoform0.linkSpectrum(self.spectra[spectrumID])


    def updateChimericSpectra(self, maxRank):

        """ For every psm of rank = rank try to find a matching envelope, and assign to that proteoform it if it is the case"""
        for spectrum in self.spectra.values():
            spectrumMz = spectrum.getPrecMz()
            spectrumRt = spectrum.get_rt()
            
            for rank in range(1,maxRank):
                if len(spectrum.psms) >= rank:
                    psm = spectrum.psms[rank-1]
                    psmProforma = psm.getModificationProforma()
                    
                    altProteo = self.proteoforms[psmProforma] #TODO might need a try except for proteo only in second rank

                    if altProteo.getEnvelope() != None: 
                        if altProteo.getEnvelope().get_y(spectrumRt) > self.intensityThreshold:
                            
                            #print("spectrum at RT {0} psm: {1} matches proteoform {2}".format(spectrumRt, psmProforma, altProteo.getModificationBrno()))
                            psm.isValidated = True 
                            altProteo.linkPsm(psm)
                            spectrum.updateRatio()
                            
                       


                #isobaricProteoforms = [proteo for proteo in self.proteoforms.values if proteo.getTheoPrecMz() > spectrumMz - self.precMzTolerance and spectrumMz + self.precMzTolerance > proteo.getTheoPrecMz() and proteo.getModificationProforma() == psmProforma]

    #Visualization
    #Methods here should return 
