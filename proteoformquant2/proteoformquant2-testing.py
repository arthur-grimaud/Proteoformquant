#!/usr/bin/env python

### Import ###
# Modules
import sys
from pyteomics import mzid
from objbrowser import browse  # dev
import pickle

# Custom Modules
from Utils import input

# Classes
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
    args = input.doArgs(sys.argv[1:], progName)  # Parse arguments
    args = input.checkArgs(args)  # Verify arguments

    verbose = args.verbose
    indent_file = args.indent_file
    spectra_file = args.spectra_file
    output_dir = args.output_dir
    dbse = args.dbse

    print("Starting " + progName)

    # --------------------------------- Analysis --------------------------------- #

    ### Read Data ###
    run = Msrun(run_id="1", dbse=dbse)
    run.read_mzid(indent_file)
    # run.read_mgf(spectra_file)
    run.read_mgf(spectra_file)

    ### Prepare Data ###
    run.add_proteoforms()
    run.match_fragments()

    ### Export Fragment Annotation ###
    # run.get_dataframe_fragment_annotation().to_csv(output_dir + ".csv")

    with open("test_save_1_1.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

    with open("test_save_1_1.pkl", "rb") as inp:
        run = pickle.load(inp)

    # For testing purpose filter protroforms:
    # proteoform_to_keep = {
    #     "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATK[Acetyl]AARKSAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARK[Acetyl]STGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Trimethyl]KPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRK[Acetyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGK[Acetyl]APRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRK[Trimethyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Acetyl]PHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Trimethyl]PHRYRPGTVALRE",
    # }
    # only correct
    proteoform_to_keep = {
        "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",
    }
    p_del_list = []
    for key in run.proteoforms.keys():
        if key not in proteoform_to_keep:
            p_del_list.append(key)
            for psm in run.proteoforms[key].get_linked_psm():
                psm.spectr.get_psms()[psm.rank - 1] = 0
                del psm

    for key in run.spectra.keys():
        if len(run.spectra[key.get_psms()) == 0:
            try:
                del run.spectra[key]
                print("deleted spec")
            except KeyError:
                pass

    for key in p_del_list:
        del run.proteoforms[key]
        print("deleted proteo")

    for spectra in run.spectra.values():
        spectr.get_psms() = [i for i in spectr.get_psms() if i != 0]

    ### Quantification groups ###

    run.set_proteoform_isobaric_groups()
    run.optimize_proteoform_subsets()

    with open("test_res_1_1.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

    ### Quantification Standard ###
    # run.update_proteoforms_elution_profile()
    # run.update_psm_validation()
    # # run.update_unassigned_spectra()
    # # run.update_chimeric_spectra(max_rank=10)
    # # run.update_proteoforms_elution_profile()
    # # run.update_psm_validation()
    # # run.update_unassigned_spectra()
    # run.update_proteoform_intens()

    # with open("pfq_out_obj_1.pkl", "wb") as outp:
    #     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    # run.result_dataframe_pfq1_format().to_csv("pfq_out_obj_1.csv")


if __name__ == "__main__":
    # sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
