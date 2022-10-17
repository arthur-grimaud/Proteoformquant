#!/usr/bin/env python
# python3 -m cProfile -o program.prof proteoformquant2.py -i Data/mix1_1pmol_rep2_renamed.mzid -s Data/mix1_1pmol_rep2_renamed.mgf -d mascot
# snakeviz program.prof
### Import ###
# Modules
import sys
from pyteomics import mzid
from objbrowser import browse  # dev
import pickle
import csv

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

    progName = "ProteoformQuant2"

    ### Input ###
    args = input.doArgs(sys.argv[1:], progName)  # Parse arguments
    args = input.checkArgs(args)  # Verify arguments

    verbose = args.verbose
    indent_file = args.indent_file
    spectra_file = args.spectra_file
    output_dir = args.output_dir

    # Name of output files from input identification filename
    output_prefix = indent_file.split(".")[0].split("/")[1]
    print(output_prefix)

    # --------------------------------- Debugging -------------------------------- #

    # --------------------------------- Analysis --------------------------------- #

    print("---===Starting " + progName + "===---")

    ### Read Data ###
    run = Msrun(run_id="1")
    run.read_mzid(indent_file)
    run.read_mgf(spectra_file)

    run.add_metrics_to_log(processing_step="initial")

    ### Prepare Data ###
    run.add_proteoforms()

    # run.update_proteoforms_elution_profile()
    run.update_proteoform_intens()
    ### Print Quant results ###
    f = open(f"quant_raw_{output_prefix}.csv", "w")
    for proteoform in run.proteoforms.values():
        if proteoform.update_proteoform_total_intens(method="precursor") > 0:
            f.write(
                "".join(
                    [
                        str(proteoform.get_modification_brno()),
                        ",",
                        str(proteoform.get_peptide_sequence()),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="precursor")),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="AUC")),
                        ",",
                        str(proteoform.is_ambiguous()),
                        "\n",
                    ]
                )
            )
    f.close()

    run.fdr_filtering(decoy_tag="decoy_", score_name="Amanda:AmandaScore")
    run.filter_proteform_low_count(min_n_psm=2)
    # run.retention_time_window_filtering(0, 5300)
    run.scale_precursor_intensities()
    run.match_fragments()

    #### SAVE ###
    with open(f"save_{output_prefix}.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    with open(f"save_{output_prefix}.pkl", "rb") as inp:
        run = pickle.load(inp)

    run.update_proteoforms_elution_profile()
    run.update_proteoform_intens()
    run.add_metrics_to_log(processing_step="First_rank_quantification")

    # ### SAVE ###
    with open(f"save_inter_{output_prefix}.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    with open(f"save_inter_{output_prefix}.pkl", "rb") as inp:
        run = pickle.load(inp)

    ### Print Quant results ###
    f = open(f"quant_initial_{output_prefix}.csv", "w")
    for proteoform in run.proteoforms.values():
        if proteoform.update_proteoform_total_intens(method="precursor") > 0:
            f.write(
                "".join(
                    [
                        str(proteoform.get_modification_brno()),
                        ",",
                        str(proteoform.get_peptide_sequence()),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="precursor")),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="AUC")),
                        ",",
                        str(proteoform.is_ambiguous()),
                        "\n",
                    ]
                )
            )
    f.close()

    ### Optimize EP ###
    # print(run.proteoforms)
    run.set_proteoform_isobaric_groups()

    # ### SAVE ###
    with open(f"save_res_{output_prefix}.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    with open(f"save_res_{output_prefix}.pkl", "rb") as inp:
        run = pickle.load(inp)

    run.optimize_proteoform_subsets_2()

    # ### SAVE ###
    # with open(f"save_res_{output_prefix}.pkl", "wb") as outp:
    #     pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    # with open(f"save_res_{output_prefix}.pkl", "rb") as inp:
    #     run = pickle.load(inp)

    ### Update Quatification ###
    # run.update_proteoforms_elution_profile()
    # run.update_proteoform_intens()

    ### Print Quant results ###
    run.update_proteoform_intens()
    f = open(f"quant_opti_{output_prefix}.csv", "w")
    for proteoform in run.proteoforms.values():
        if proteoform.update_proteoform_total_intens(method="precursor") > 0:
            f.write(
                "".join(
                    [
                        str(proteoform.get_modification_brno()),
                        ",",
                        str(proteoform.get_peptide_sequence()),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="precursor")),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="AUC")),
                        ",",
                        str(proteoform.is_ambiguous()),
                        "\n",
                    ]
                )
            )
    f.close()

    pd.DataFrame(run.log).to_csv(f"log_{output_prefix}.csv", ",")

    run.validate_first_rank_no_id()

    ### Print Quant results ###
    run.update_proteoform_intens()
    f = open(f"quant_opti2_{output_prefix}.csv", "w")
    for proteoform in run.proteoforms.values():
        if proteoform.update_proteoform_total_intens(method="precursor") > 0:
            f.write(
                "".join(
                    [
                        str(proteoform.get_modification_brno()),
                        ",",
                        str(proteoform.get_peptide_sequence()),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="precursor")),
                        ",",
                        str(proteoform.update_proteoform_total_intens(method="AUC")),
                        ",",
                        str(proteoform.is_ambiguous()),
                        "\n",
                    ]
                )
            )
    f.close()

    ### SAVE ###
    with open(f"save_res_{output_prefix}.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    # with open(f"save_res_{output_prefix}.pkl", "rb") as inp:
    #     run = pickle.load(inp)


if __name__ == "__main__":
    # sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
