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


def main():
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
    #run.getPlotPrecVsMsMs() #viz

    ### Quantification ###
    run.updateProteoformsEnvelope()


    ### Output ###

    ### Report ###

    with open('pfq_out_obj.pkl', 'wb') as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)




    

if __name__ == '__main__':
    #sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
