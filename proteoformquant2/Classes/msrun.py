### Import ###
from pyteomics import mzid 
from pyteomics import mgf
from pyteomics import mzml
from chart_studio.plotly import plot, iplot
import plotly.express as px
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from statistics import mean
from statistics import median
from progress.bar import Bar
import pandas as pd
import networkx as nx
from itertools import combinations

#Custom classes
from Classes.spectrum import Spectrum
from Classes.proteoform import Proteoform
from Classes.proteoform0 import Proteoform0

from Utils import misc

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
        
        self.proteoform_isobaric_group = []

        #PARAMETERS
        self.prec_mz_tol = 1.5   #in Da
        self.frag_mz_tol = 0.015  #in Da
        self.intensityThreshold = 20000
        self.elution_profile_score_threshold = 0.3
        self.fragments_types = ["c","zdot","z+1","z+2","c-zdot", "c-z+1", "cdot-zdot", "c-z+1", "a-x"]


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


    # ------------------------- DATA PREPARATION METHODS ------------------------ #
 
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

        self.spectra_fn = spectra_fn #Store File name that has been read
        mgf_obj = mgf.read(spectra_fn, index_by_scans=True)


        with Bar('loading spectra', max=1) as bar:

            for specID in self.spectra: #Take into account how DBSEs store spectra ids 

            

                if self.dbse == "mascot":
                    index = str(int(specID.split("=")[1]) + 1)
                if self.dbse == "comet":
                    index = specID.split("=")[1]


                try:
                    specMgf = mgf_obj.get_spectrum(index) #need to be splited 
                except(KeyError):
                    specMgf = mgf_obj.get_spectrum

                self.spectra[specID].set_spec_data_mgf(specMgf)
                bar.next()
        pass

    def add_mzml_data(self, spectra_fn):
        """Add info from mzml file to spectrum objects in self.spectra add Proteform objects to self.Proteoforms"""

        self.spectra_fn = spectra_fn #Store File name that has been read
        mzml_obj = mzml.read(spectra_fn)
        
        with Bar('loading spectra', max=1) as bar:

            for specID in self.spectra: #Take into account how DBSEs store spectra ids 

                

                if self.dbse == "mascot":
                    index = str(int(specID.split("=")[1]) + 1)
                if self.dbse == "comet":
                    index = specID.split("=")[1]


                print(index)
                try:
                    spec_mzml = mzml_obj.get_by_id(index, id_key="index") #need to be splited 
                except(KeyError):
                    print("Key error spectrum not found")
                
                print(spec_mzml)

                self.spectra[specID].set_spec_data_mzml(spec_mzml)
                bar.next()
        pass

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


    # ------------------------- QUANTIFICATION METHODS ------------------------ #


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
    
    def filter_proteform_low_count(self, min_n_psm ):
        for proteoform in self.proteoforms.values():
            if(len(proteoform.get_validated_linked_psm()) < min_n_psm):
                for psm in proteoform.get_linked_psm():
                    psm.isValidated == False


    


    def update_chimeric_spectra(self, max_rank):

        """ For every psm of rank < rank try to find a matching envelope, and assign to that proteoform it if it is the case"""
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

    def get_dataframe_fragment_annotation(self, max_rank = 1):
        "Create a dataframe where each row is the annotation for one PSM (TAKES FIRST RANK)"


        df = pd.DataFrame(columns=('scan','sequence', 'brno', 'proforma', 'prec_intens', "sum_intens_frag","prec_in_msms"))

        for spectrum in self.spectra.values():
            
            
            for psm in spectrum.psms[:max_rank]:

                print(psm.spectrum.get_id())
                print(psm.proteoform.get_protein_ids())
                #add information spectrum and psm
                dict_add = {'scan':psm.spectrum.get_id(),
                            'sequence':psm.proteoform.peptideSequence, 
                            'brno':psm.get_modification_brno(), 
                            'proforma':psm.proteoform.get_modification_proforma(), 
                            'prec_intens':psm.spectrum.getPrecIntens(),
                            'sum_intens_frag':psm.spectrum.getSumIntensFrag(),
                            'prec_in_msms': psm.spectrum.get_sum_intens_frag_at_mz([psm.spectrum.getPrecMz()-20, psm.spectrum.getPrecMz()+20])   }
                
                
                #add annotation:
                #print(psm.get_annotation())
                for frag_type in psm.get_annotation().values():
                    #print(frag_type)

                    annotation_frag_type = {frag_type["fragCode"][i]:frag_type["intens"][i] for i in range(0,len(frag_type["fragCode"]))}

                    #print(annotation_frag_type)
                    dict_add.update(annotation_frag_type)
                
                print(dict_add)
                df_add = pd.DataFrame(dict_add, index=[0]) 

                dfs = [df, df_add]

                df = pd.concat(dfs)


        return df

    # ----------------------- TESTING GROUP QUANTIFICATION ----------------------- #

    def set_proteoform_isobaric_groups(self):
        """ From self.proteoforms define self.proteoform_isobaric_group where isobaric proteoform are grouped """

        all_proteoform_pairs = [] #pairs of proteoform founds together in the psms of a spectra

        for spectrum in self.spectra.values():
            proforma_psms =  [psm.get_modification_proforma() for psm in spectrum.psms]
            proforma_combinations = [(a, b) for idx, a in enumerate(proforma_psms) for b in proforma_psms[idx + 1:]]

            for c in proforma_combinations:
                if c not in all_proteoform_pairs:
                    all_proteoform_pairs.append(c)

        #find groups (connected graphs) in the network defined from the edgelist "all-proteoform-pairs"
        G=nx.from_edgelist(all_proteoform_pairs)
        l=list(nx.connected_components(G))
        print(l)
        self.proteoform_isobaric_group = l


    def update_psms_validation_proteoform_subset(self, spectra_subset, proteoforms_proforma):
        for spectrum in spectra_subset:
            for psm in spectrum.psms:
                if psm.proteoform.get_modification_proforma() in proteoforms_proforma:
                    psm.isValidated = True
                else:
                    psm.isValidated = False

    def update_psms_ratio_subset(self, spectra_subset):
        for spectrum in spectra_subset:
            spectrum.update_psms_ratio()

    def update_proteoforms_elution_profile_subset(self, proteoform_subset):
        """For each proteoforms in self. proteoforms model the elution profile"""
        with Bar('Updating proteoforms envelopes', max=1) as bar:

            #TODO add a function thjat set the bounds based on the entire set of envelopes  
            for proteoform in proteoform_subset:
                proteoform.model_elution_profile(self.elution_profile_score_threshold)
                bar.next()
        pass


    def find_optimal_proteoform_set(self):

        for group in self.proteoform_isobaric_group:
            print("PROTEFORM ISOBARIC GROUP:  ",  group)

            #create a list of PSMs linked count for each proteoform in the group
            group_psm_count = [len(self.proteoforms[proforma].get_linked_psm()) for proforma in group]
            sorted_group = [x for _, x in sorted(zip(group_psm_count, group), reverse=True)]

            for proteoform_proforma in sorted_group:
                print("proteoform linked psm", len(self.proteoforms[proteoform_proforma].get_linked_psm()))

               
                # p


                # for proforma_subset in misc.combinations(proteo_group):
                    
                #     #Get the set of spectra in the proteoform group:
                #     #Get the set of proteoform object in the proteform group:
                #     spectra_subset = []
                #     proteoform_subset = []
                #     for proforma in proforma_subset:
                #         proteoform = self.proteoforms[proforma]
                #         proteoform_subset.append(proteoform)
                #         for psm in proteoform.get_linked_psm():
                #             if psm.spectrum not in spectra_subset:
                #                 spectra_subset.append(psm.spectrum)

                #     print("Proteoform combination:  ",  [p.get_modification_brno() for p in proteoform_subset])

                #     self.update_psms_validation_proteoform_subset(spectra_subset, proforma_subset)
                #     self.update_psms_ratio_subset(spectra_subset)
                #     self.update_proteoforms_elution_profile_subset(proteoform_subset)


                #     residuals = [spectrum.quant_residuals for spectrum in spectra_subset if spectrum.quant_residuals !=0]
                #     if len(residuals) > 1: 
                #         print(mean(residuals))
                #     else: 
                #         print("no mean")


                
            






            
            

        






