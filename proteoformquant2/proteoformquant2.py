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
    #run.add_mgf_data(spectra_fn)
    run.add_mgf_data(spectra_fn)
    
    ### Prepare Data ###
    run.add_proteoforms()
    run.match_fragments()

    ### Export Fragment Annotation ###
    run.get_dataframe_fragment_annotation().to_csv(outputFn + ".csv")

    
    ### Quantification ###
    

    #run.update_proteoforms_elution_profile()
    #run.update_psm_validation()
    #run.update_unassigned_spectra()
    #run.update_psms_ratio()

    


    # with open('pfq_out_obj.pkl', 'wb') as outp:
    #     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    # run.result_dataframe_pfq1_format().to_csv('pfq_out_obj.csv')








if __name__ == '__main__':
    #sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
