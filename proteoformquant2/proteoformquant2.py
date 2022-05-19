#!/usr/bin/env python
# python3 -m cProfile -o program.prof proteoformquant2.py -i Data/mix1_1pmol_rep2_renamed.mzid -s Data/mix1_1pmol_rep2_renamed.mgf -d mascot
# snakeviz program.prof
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

    # print(resource.getrlimit(resource.RLIMIT_STACK))
    # print(sys.getrecursionlimit())

    max_rec = 0x100000

    # May segfault without this line. 0x100 is a guess at the size of each stack frame.
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)

    progName = "ProteoformformQuant2"

    ### Input ###
    args = input.doArgs(sys.argv[1:], progName)  # Parse arguments
    args = input.checkArgs(args)  # Verify arguments

    verbose = args.verbose
    indentFn = args.indentFn
    spectra_fn = args.spectra_fn
    outputFn = args.outputFn
    dbse = args.dbse

    print("Starting " + progName)

    # --------------------------------- Debugging -------------------------------- #

    # import sys
    # import traceback

    # class TracePrints(object):
    #     def __init__(self):
    #         self.stdout = sys.stdout

    #     def write(self, s):
    #         self.stdout.write("Writing %r\n" % s)
    #         traceback.print_stack(file=self.stdout)

    # sys.stdout = TracePrints()
    # --------------------------------- Analysis --------------------------------- #

    # ### Read Data ###
    # run = Msrun(run_id="1", dbse=dbse)
    # run.read_mzid(indentFn)
    # # run.read_mgf(spectra_fn)
    # run.read_mgf(spectra_fn)

    # ### Prepare Data ###

    # run.fdr_filtering(decoy_tag="decoy_", score_name="Comet:xcorr")
    # run.add_proteoforms()
    # run.filter_proteform_low_count(min_n_psm=5)
    # run.match_fragments()

    # ### Export Fragment Annotation ###
    # # run.get_dataframe_fragment_annotation().to_csv(outputFn + ".csv")

    # with open("test_save_1_1.pkl", "wb") as outp:
    #     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

    with open("test_save_1_1.pkl", "rb") as inp:
        run = pickle.load(inp)

    # For testing purpose filter protroforms:
    proteoform_to_keep = {
        "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATK[Acetyl]AARKSAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARK[Acetyl]STGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Trimethyl]KPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRK[Acetyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGK[Acetyl]APRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRK[Trimethyl]QLATKAARKSAPATGGVKKPHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Acetyl]PHRYRPGTVALRE",
        "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKK[Trimethyl]PHRYRPGTVALRE",
    }
    # only correct
    # proteoform_to_keep = {
    #     "ARTKQTARKSTGGKAPRKQLATK[Trimethyl]AARKSAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARK[Trimethyl]SAPATGGVKKPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[Acetyl]KPHRYRPGTVALRE",
    #     "ARTKQTARKSTGGKAPRKQLATKAARK[Acetyl]SAPATGGVKKPHRYRPGTVALRE",
    # }
    p_del_list = []
    for key in run.proteoforms.keys():
        if key not in proteoform_to_keep:
            p_del_list.append(key)
            for psm in run.proteoforms[key].get_linked_psm():
                psm.exclude()

    for key in p_del_list:
        del run.proteoforms[key]
    print("deleted proteo")

    for spectrum in run.spectra.values():
        spectrum.psms = [i for i in spectrum.get_psms() if i != 0]

    ### Quantification groups ###

    run.set_proteoform_isobaric_groups()
    run.find_optimal_proteoform_set_2()

    with open("test_res_1_1.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

    ### Quantification Standard ###

    # run.update_proteoforms_elution_profile()
    # run.update_psm_validation()
    # run.update_unassigned_spectra()
    # run.update_chimeric_spectra(max_rank=10)
    # run.update_psm_validation()
    # run.update_proteoforms_elution_profile()
    # run.update_unassigned_spectra()
    # run.update_proteoform_intens()

    # with open("test_res_1_1.pkl", "wb") as outp:
    #     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    # run.result_dataframe_pfq1_format().to_csv("pfq_out_obj_1.csv")


if __name__ == "__main__":
    # sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
