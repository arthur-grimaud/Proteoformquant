### Import ###
from pyteomics import mzid 
from pyteomics import mgf
#Custom classes
from Classes.spectrum import Spectrum
from Classes.proteoform import Proteoform

import pprint
class Msrun():

    def __init__(self, runId:str = "Default run ID" ):
        
        self.runId = runId

        self.identFn: str = "Not Specified" 
        self.spectraFn: str = "Not Specified"

        self.spectra: dict(Spectrum) = {} 
        self.proteoforms: dict(Proteoform) = {}
 

    def readMzid(self, identFn):
        """Read a spectra identification file in .mzIdenMl whose path is specfieid in self.inputFn"""
        self.identFn = identFn #Store File that has been read
        mzidObj = mzid.read(identFn) #Create a pyteomics' mzid iterator

        for identMzid in mzidObj: #Iterate over spectra and create Spectrum object for each 
            print(identMzid["spectrumID"])
            self.spectra[identMzid["spectrumID"]] = Spectrum(spectrumID= identMzid["spectrumID"], identMzid = identMzid)
        pass

    def addMgfData(self, spectraFn):
        """Add info from mgf file to spectrum objects in self.spectra"""
        print("start reading mgf:\n")
        self.spectraFn = spectraFn #Store File that has been read
        mgfObj = mgf.read(spectraFn)

        for specID in self.spectra:
            specMgf = mgfObj.get_spectrum(specID.split("=")[1]) #need to be splited 
            self.spectra[specID].setSpecDataMgf(specMgf)
        pass

    def addMzmlData(self):
        """Add info from mzml file to spectrum objects in self.spectra add Proteform objects to self.Proteoforms"""
        pass
    
    def addProteoforms(self):
        """From spectrum objects in self.psectra add proteoforms object to self.proteoforms"""
        for specID in self.spectra:
            for psm in self.spectra[specID].psms:

                brno = psm.getModificationsBrno()
                seq = psm.getPeptideSequence()
                brnoSeq = brno+"-"+seq #TODO use proforma

                if brnoSeq not in self.proteoforms.keys(): #if a proteoform is new create a instance of Proteoform for this proteoform
                    self.proteoforms[brnoSeq] = Proteoform(peptideSequence=seq, modificationBrno=brno, modificationDict=psm.Modification)  

                self.proteoforms[brnoSeq].linkPsm(psm) #add link to Psms in Proteoform
                psm.setProteoform(self.proteoforms[brnoSeq]) #add link to Proteoform in Psm

            # pprint.pprint(vars(self.proteoforms[brno+"-"+psm.PeptideSequence]))
            print(self.proteoforms[brno+"-"+psm.getPeptideSequence()].linkedPsm)


    def matchFragments(self, msmsTol= 0.02, internal = False):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        
        for proteoID in self.proteoforms:
            self.proteoforms[proteoID].setTheoreticalFragments(["c","zdot"])
            #print(self.proteoforms[proteoID].theoFrag)
        for spectrumID in self.spectra:
            self.spectra[spectrumID].annotateFragPsm()

        pass