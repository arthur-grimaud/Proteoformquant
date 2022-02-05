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
    spectraFn = args.spectraFn
    outputFn = args.outputFn
    dbse = args.dbse

    print("Starting " +  progName )

    ### Read Data ###
    run = Msrun(runId="1", dbse = dbse)
    run.readMzid(indentFn)
    run.addMgfData(spectraFn)
    
    ### Prepare Data ###
    run.addProteoforms()
    run.matchFragments()

    ### Quantification ###
    run.updateProteoformsEnvelope()
    run.updateProteoformsValidation()
    run.updateProteoformsTotalIntens()
    run.updateUnassignedSpectra()

    ### Output ###

    ### Report ###
    sys.setrecursionlimit(10000)
    with open('pfq_out_obj_test.pkl', 'wb') as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)




    

if __name__ == '__main__':
    #sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
