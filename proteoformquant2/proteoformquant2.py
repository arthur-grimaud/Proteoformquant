#!/usr/bin/env python

### Import ###

#Modules
import sys
from pyteomics import mzid 
from objbrowser import browse #dev

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

    print("Starting " +  progName )
    ### Read Data ###
    run = Msrun(runId="1")
    run.readMzid(indentFn)
    run.addMgfData(spectraFn)
    
    ### Prepare Data ###
    run.addProteoforms()
    run.matchFragments()
    ### Quantification ###


    ### Output ###



    return 

if __name__ == '__main__':
    #sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
