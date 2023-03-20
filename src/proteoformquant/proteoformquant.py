#!/usr/bin/env python
# python3 -m cProfile -o program.prof proteoformquant.py -i Data/mix1_1pmol_rep2_renamed.mzid -s Data/mix1_1pmol_rep2_renamed.mgf -d mascot
# snakeviz program.prof
### Import ###
# Modules
import sys
import pickle

# Custom Modules
from proteoformquant.Utils import input

# Classes
from proteoformquant.Classes.msrun import Msrun
import resource
import pandas as pd
import sys
from jsonc_parser.parser import JsoncParser


def main():
    progName = "proteoformquant"

    # May segfault without these lines.
    max_rec = 0x100000
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)
    #

    print("---===Starting " + progName + "===---")

    # --------------------------------- Inputs -------------------------------- #

    args, unknownargs = input.doArgs(sys.argv[1:], progName)  # Parse arguments
    args = input.checkArgs(args)  # Verify arguments

    verbose = args.verbose
    indent_file = args.indent_file
    spectra_file = args.spectra_file
    param_file = args.param_file
    output_dir = args.output_dir
    output_file = args.output_file

    # Read Parameters File:

    

    try:
        params = JsoncParser.parse_file(param_file)
    except JsoncParser.errors.FileError:
        print("Default file params.jsonc not found.")
        print("Please provide a valid parameter file with the -p (option.")
        print("you can generate a default parameter file with the -cg option.")
        sys.exit(1)

    # read parameter overwrited in cmd line arguments:
    params_over = {unknownargs[i][1:]: unknownargs[i + 1] for i in range(0, len(unknownargs), 2)}
    # Name of output prefix from input identification filename

    if output_file != None:
        output_prefix = output_file
    else:
        output_prefix = indent_file.split(".")[0].split("/")[1]

    max_rank_validation = args.rankquant

    # --------------------------------- Debugging -------------------------------- #

    # --------------------------------- Analysis --------------------------------- #

    ### Read Data ###
    run = Msrun(run_id="1", params=params, params_over=params_over, verbose=verbose)
    ###
    run.read_mzid(indent_file)
    run.read_spectra(spectra_file)
    # run.read_mzml(spectra_file)

    run.add_metrics_to_log(processing_step="initial")

    ### Prepare Data ###
    run.add_proteoforms()
    # run.deconvolute_ms2()
    run.update_proteoform_intens()

    ### Filtering Data ###
    run.fdr_filtering(decoy_tag="decoy_", score_name="Amanda:AmandaScore")
    # run.filter_psms_low_score()
    run.filter_proteform_low_count()

    ### Prepare Data ###
    run.scale_precursor_intensities()
    run.match_fragments()

    ### First Rank Quantification ###
    run.update_proteoforms_elution_profile()
    run.update_proteoform_intens()

    if max_rank_validation == 0:
        ### Chimeric Spectra Quantification ###
        run.set_proteoform_isobaric_groups()
        run.optimize_proteoform_subsets()
        run.validate_first_rank_no_id()
        run.update_proteoform_intens()
    else:  # validate n first rank and perform chimeric spectra quantification
        run.validate_psms_rank_n(rank=max_rank_validation)
        run.update_psms_ratio_all()
        run.update_proteoform_intens()

    # --------------------------------- Output --------------------------------- #

    # with open("save_res_2022_mix2_rep1.pkl", "rb") as f:
    #     run = pickle.load(f)

    # Quantification Table
    run.result_dataframe_pfq1_format().to_csv(f"{output_dir}/quant_{output_prefix}.csv")
    pd.DataFrame(run.log).to_csv(f"{output_dir}/log_{output_prefix}.csv", ",")

    # Pickle Object for report
    with open(f"{output_dir}/obj_{output_prefix}.pkl", "wb") as outp:
        pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

    # PSMs Table
    run.psm_score_dataframe(
        file_name=f"{output_dir}/psms_{output_prefix}.csv", score_name="Amanda:AmandaScore"
    )


if __name__ == "__main__":
    main()
