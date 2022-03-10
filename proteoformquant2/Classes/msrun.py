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
import pandas as pd

#Custom classes
from Classes.spectrum import Spectrum
from Classes.proteoform import Proteoform
from Classes.proteoform0 import Proteoform0

import pprint
class Msrun():

    def __init__(self, run_id:str = "Default run ID", dbse:str = "comet" ):
        
        self.run_id = run_id
        self.dbse = dbse
        self.ident_fn: str = "Not Specified" 
        self.spectra_fn: str = "Not Specified"

        self.spectra: dict(Spectrum) = {} 
        self.proteoforms: dict(Proteoform) = {}
        self.proteoform0 = Proteoform0()

        #PARAMETERS
        self.prec_mz_tol = 1.5   #in Da
        self.frag_mz_tol = 0.015  #in Da
        self.intensityThreshold = 200000
        self.elution_profile_score_threshold = 0.5
        self.fragments_types = ["c","zdot","z+1","z+2"]


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

    def get_dataset_metrics(self):
        """returns a set of metric as a dictionnary for the ms run"""

        return {
            "TotalSpectra": len(self.spectra),
            "AssignedSpectra(method1)": len([spectrum for spectrum in self.spectra.values() if spectrum.get_number_validated_psm() > 0]),
            "AssignedSpectra(method2)": 0,
            "UnassignedSpectra": len([spectrum for spectrum in self.proteoform0.linkedSpectra]),
            "ChimericSpectra": len([spectrum for spectrum in self.spectra.values() if spectrum.get_number_validated_psm() > 1]),
            "TotalProteoforms": len(self.proteoforms)-1,
            #"AvgEnvelopeScore": mean([proteoform.get_elution_profile().score_fitted for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            #"MinEnvelopeScore": min([proteoform.get_elution_profile().score_fitted for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            #"MedianEnvelopeS":   median([proteoform.get_elution_profile().param_fitted[1] for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            #"MedianEnvelopeA":   median([proteoform.get_elution_profile().param_fitted[2] for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ]),
            #"MedianEnvelopeK":   median([proteoform.get_elution_profile().param_fitted[3] for proteoform in self.proteoforms.values() if proteoform.get_elution_profile() != None ])
        }

    # ----------------------------------- main ----------------------------------- #
 
    def read_mzid(self, ident_fn):
        """Read a spectra identification file in .mzIdenMl whose path is specfieid in self.inputFn"""
        self.ident_fn = ident_fn #Store File that has been read
        mzidObj = mzid.read(ident_fn) #Create a pyteomics' mzid iterator

        with Bar('loading identifications', max=1) as bar:

            for identMzid in mzidObj: #Iterate over spectra and create Spectrum object for each 

                self.spectra[identMzid["spectrumID"]] = Spectrum(spectrumID= identMzid["spectrumID"], identMzid = identMzid)
                bar.next()
        pass

    def add_mgf_data(self, spectra_fn):
        """Add info from mgf file to spectrum objects in self.spectra"""

        self.spectra_fn = spectra_fn #Store File that has been read
        mgfObj = mgf.read(spectra_fn)


        with Bar('loading spectra', max=1) as bar:

            for specID in self.spectra: #Take into account how DBSEs store spectra ids 
                if self.dbse == "mascot":
                    index = str(int(specID.split("=")[1]) + 1)
                if self.dbse == "comet":
                    index = specID.split("=")[1]

                specMgf = mgfObj.get_spectrum(index) #need to be splited 
                self.spectra[specID].set_spec_data_mgf(specMgf)
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

                    proforma = psm.get_modification_proforma()
                    brno = psm.get_modification_brno()
                    seq = psm.getPeptideSequence()
                    brnoSeq = brno+"-"+seq #TODO use proforma

                    #Create list with associated protein description (try except as variable name could be different depending on the DBSE used)
                    try:
                        protein_ids = [ref["DBSequence_Ref"] for ref in psm.PeptideEvidenceRef]
                    except KeyError:
                        protein_ids = [ref['protein description'] for ref in psm.PeptideEvidenceRef]


                    if proforma not in self.proteoforms.keys(): #if a proteoform is new create a instance of Proteoform for this proteoform
                        self.proteoforms[proforma] = Proteoform(peptideSequence=seq, modificationBrno=brno, modificationProforma=proforma, modificationDict=psm.Modification, protein_ids=protein_ids).set_color(i)  

                    self.proteoforms[proforma].link_psm(psm) #add link to Psms in Proteoform
                    psm.setProteoform(self.proteoforms[proforma]) #add link to Proteoform in Psm
                    bar.next()



    def match_fragments(self, internal = False):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""

        with Bar('Generating theoretical fragments', max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].compute_theoretical_fragments(self.fragments_types)
                bar.next()

        with Bar('Matching fragments', max=1) as bar:
            for spectrumID in self.spectra:
                self.spectra[spectrumID].annotateFragPsm(frag_mz_tol = self.frag_mz_tol)
                self.spectra[spectrumID].setSumIntensAnnotFrag()
                bar.next()

        pass

    
    def update_proteoforms_elution_profile(self):
        """For each proteoforms in self. proteoforms model the elution profile"""
        with Bar('Updating proteoforms envelopes', max=1) as bar:

            #TODO add a function thjat set the bounds based on the entire set of envelopes  
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].model_elution_profile(self.elution_profile_score_threshold)
                bar.next()
        pass


    def update_proteoform_intens(self):
        """If mgf and identification data are provided in a spectrum object, get the annotated fragments for each PSM"""
        with Bar('Updating Proteoforms Total Intens', max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].update_proteoform_total_intens()
                bar.next()
        pass

    def update_psm_validation(self):
        """Update psm.isValidated"""
        with Bar('Updating PSM Validation', max=1) as bar:
            for proteoID in self.proteoforms:
                self.proteoforms[proteoID].update_proteoform_psm_validation()
        pass

    def update_unassigned_spectra(self):
        """add spectra without any validated psm to "proteoform0" in self.proteoforms"""
        for spectrumID in self.spectra:
            if self.spectra[spectrumID].get_number_validated_psm() == 0:
                self.proteoform0.linkSpectrum(self.spectra[spectrumID])


    def validate_all_psms(self):
        for spectrum in self.spectra.values():
            for psm in spectrum.psms:
                psm.isValidated = True
            


    def update_psms_ratio(self):
        for spectrum in self.spectra.values():
            spectrum.update_psms_ratio()


    def update_chimeric_spectra(self, max_rank):

        """ For every psm of rank = rank try to find a matching envelope, and assign to that proteoform it if it is the case"""
        for spectrum in self.spectra.values():
            spectrumMz = spectrum.getPrecMz()
            spectrumRt = spectrum.get_rt()
            
            for rank in range(1,max_rank):
                if len(spectrum.psms) >= rank:
                    psm = spectrum.psms[rank-1]
                    psmProforma = psm.get_modification_proforma()
                    
                    altProteo = self.proteoforms[psmProforma] #TODO might need a try except for proteo only in second rank

                    if altProteo.get_elution_profile() != None: 
                        if altProteo.get_elution_profile().get_y(spectrumRt) > self.intensityThreshold:
                            
                            #print("spectrum at RT {0} psm: {1} matches proteoform {2}".format(spectrumRt, psmProforma, altProteo.get_modification_brno()))
                            psm.isValidated = True 
                            #altProteo.link_psm(psm)
    
                   
    def result_dataframe_pfq1_format(self):
        df = pd.DataFrame(columns=('protein','sequence', 'brno', 'proforma', 'intensity', 'linked_psm', 'linked_psm_validated', 'rt_peak', 'auc'))

        for proteo in self.proteoforms.values():

            #get Elution profile max peak rt
            if proteo.get_elution_profile() != None:
                rt_peak = proteo.get_elution_profile().get_x_at_max_y()
                auc = proteo.get_elution_profile().get_auc(rt_peak-1000, rt_peak+1000) #TODO hard coded
            else:
                rt_peak = "NA"
                auc="NA"
            
            df.loc[len(df)]= [
                proteo.get_protein_ids(),
                proteo.peptideSequence, 
                proteo.get_modification_brno(), 
                proteo.get_modification_proforma(), 
                proteo.get_proteoform_total_intens(), 
                len(proteo.get_linked_psm()), 
                len(proteo.get_validated_linked_psm()), 
                rt_peak,
                auc
                ]

        return df
