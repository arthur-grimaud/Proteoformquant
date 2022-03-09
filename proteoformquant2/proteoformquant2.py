#!/usr/bin/env python

### Import ###
#Modules
import sys
from pyteomics import mzid 
from objbrowser import browse #dev
import pickle
#Custom Modules
from Utils import input 
#Classes
from Classes.msrun import Msrun 
import resource
import pandas as pd

def main():

    
    import sys

    print(resource.getrlimit(resource.RLIMIT_STACK))
    print(sys.getrecursionlimit())

    max_rec = 0x100000

    # May segfault without this line. 0x100 is a guess at the size of each stack frame.
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)

    progName = "ProteoformformQuant2"

    ### Input ###
    args = input.doArgs(sys.argv[1:], progName) #Parse arguments
    args = input.checkArgs(args) #Verify arguments 

    verbose = args.verbose
    indentFn = args.indentFn
    spectra_fn = args.spectra_fn
    outputFn = args.outputFn
    dbse = args.dbse

    print("Starting " +  progName )

    # --------------------------------- Analysis --------------------------------- #

    ### Read Data ###
    run = Msrun(run_id="1", dbse = dbse)
    run.read_mzid(indentFn)
    run.add_mgf_data(spectra_fn)
    
    ### Prepare Data ###
    run.add_proteoforms()
    run.match_fragments()
    


    ## Quantification ### First rank only no valid
    run.update_proteoform_intens()
    ### Report ###
    print(run.get_dataset_metrics())
    sys.setrecursionlimit(10000)
    with open('pfq_out_obj_test_1o.pkl', 'wb') as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    run.result_dataframe_pfq1_format().to_csv('pfq_out_obj_test_3o.csv')  
     

    ## Quantification ### First rank only valid elution profile 
    run.update_proteoforms_elution_profile()
    run.update_psm_validation()
    run.update_unassigned_spectra()
    run.update_proteoform_intens()
    ### Report ### 
    print(run.get_dataset_metrics())
    sys.setrecursionlimit(10000)
    with open('pfq_out_obj_test_1a.pkl', 'wb') as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    run.result_dataframe_pfq1_format().to_csv('pfq_out_obj_test_3a.csv')  



    ## Quantification ### chimeric valid elution profile
    run.update_chimeric_spectra(max_rank=100)
    run.update_psms_ratio()
    run.update_proteoforms_elution_profile()
    run.update_psm_validation()
    run.update_unassigned_spectra()
    run.update_proteoform_intens()
    ### Report ###
    print(run.get_dataset_metrics())
    sys.setrecursionlimit(10000)
    with open('pfq_out_obj_test_1b.pkl', 'wb') as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    run.result_dataframe_pfq1_format().to_csv('pfq_out_obj_test_3b.csv')




if __name__ == '__main__':
    #sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
