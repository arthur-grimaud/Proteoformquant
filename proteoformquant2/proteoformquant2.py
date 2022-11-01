#!/usr/bin/env python
# python3 -m cProfile -o program.prof proteoformquant2.py -i Data/mix1_1pmol_rep2_renamed.mzid -s Data/mix1_1pmol_rep2_renamed.mgf -d mascot
# snakeviz program.prof
### Import ###
# Modules
import sys
import pickle

# Custom Modules
from Utils import input

# Classes
from Classes.msrun import Msrun
import resource
import pandas as pd
import sys
from jsonc_parser.parser import JsoncParser


def main():
    progName = "ProteoformQuant2"

    # May segfault without these lines.
    max_rec = 0x100000
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)

    # --------------------------------- Inputs -------------------------------- #

    args, unknwownargs = input.doArgs(sys.argv[1:], progName)  # Parse arguments
    args = input.checkArgs(args)  # Verify arguments

    verbose = args.verbose
    indent_file = args.indent_file
    spectra_file = args.spectra_file
    param_file = args.param_file
    output_dir = args.output_dir
    output_file = args.output_file
    # Read Parameters File:
    params = JsoncParser.parse_file(param_file)
    # read parameter overwited in cmd line arguments:
    params_over = {unknwownargs[i][1:]: unknwownargs[i + 1] for i in range(0, len(unknwownargs), 2)}
    # Name of output prefix from input identification filename
    if output_file != False:
        output_prefix = output_file

    else:
        output_prefix = indent_file.split(".")[0].split("/")[1]

    # --------------------------------- Debugging -------------------------------- #

    # --------------------------------- Analysis --------------------------------- #

    print("---===Starting " + progName + "===---")

    ### Read Data ###
    run = Msrun(run_id="1", params=params, params_over=params_over)
    run.read_mzid(indent_file)
    run.read_mgf(spectra_file)

    run.add_metrics_to_log(processing_step="initial")

    ### Prepare Data ###
    run.add_proteoforms()
    run.update_proteoform_intens()

    ### Filtering Data ###
    run.fdr_filtering(decoy_tag="decoy_", score_name="Amanda:AmandaScore")
    run.filter_proteform_low_count(min_n_psm=2)

    ### Prepare Data ###
    run.scale_precursor_intensities()
    run.match_fragments()

    ### First Rank Quantification ###
    run.update_proteoforms_elution_profile()
    run.update_proteoform_intens()
    run.add_metrics_to_log(processing_step="First_rank_quantification")

    ### Chimeric Spectra Quantification ###
    run.set_proteoform_isobaric_groups()
    run.optimize_proteoform_subsets_2()
    run.validate_first_rank_no_id()
    run.update_proteoform_intens()

    ### OUTPUT ###
    with open(f"save_res_{output_prefix}.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)
    run.result_dataframe_pfq1_format().to_csv(f"quant_{output_prefix}.csv")
    pd.DataFrame(run.log).to_csv(f"log_{output_prefix}.csv", ",")


if __name__ == "__main__":
    # sys.argv = ["programName.py","--input","test.txt","--output","tmp/test.txt"]
    main()
