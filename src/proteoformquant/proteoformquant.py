import sys
import pickle
import os
import resource
import pandas as pd
from jsonc_parser.parser import JsoncParser
from pathlib import Path

try:  # local modules
    from proteoformquant.Utils import input
    from proteoformquant.Classes.msrun import Msrun
except ImportError:
    from Utils import input
    from Classes.msrun import Msrun

import os
import pandas as pd
import pickle


def main():
    progName = "proteoformquant"

    # May segfault without these lines.
    max_rec = 0x100000
    resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
    sys.setrecursionlimit(max_rec)

    print("\n---===Starting " + progName + "===---")

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

    # read parameter overwriten in cmd line arguments:
    params_over = {
        unknownargs[i][1:]: unknownargs[i + 1] for i in range(0, len(unknownargs), 2)
    }
    # Name of output prefix from input identification filename

    if output_file != None:
        output_prefix = output_file
    else: # if output file is not provided, use the input file name as output prefix
        output_prefix = Path(indent_file).stem
        

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
    run.fdr_filtering()
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
        
    ### Chimeric spectra quantification no elution profile ###
    # Go through all spectra and perform quantification of the top 2 proteoforms if only the first rank was validated
    if max_rank_validation >= 2 or max_rank_validation == 0:
        #print("Quantify chimeric spectra no elution profile")
        run.quantify_chimeric_spectra_no_elution_profile()
        run.update_proteoform_intens()

    # --------------------------------- Output --------------------------------- #
    print("Outputting results")
    try:
        # Before writing output files, ensure the output directory exists
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)

        # print the output directory
        print("Output directory: ", output_dir_path)
        
        print("Output prefix: ", output_prefix)
        

        # Quantification Table
        quant_file_path = output_dir_path / f"quant_{output_prefix}.csv"
        
        print("Quant file path: ", quant_file_path)
        run.result_dataframe_pfq1_format().to_csv(quant_file_path)

        # Log Table
        log_file_path = output_dir_path / f"log_{output_prefix}.csv"
        pd.DataFrame(run.log).to_csv(log_file_path, ",")

        # Pickle Object for report
        obj_file_path = output_dir_path / f"obj_{output_prefix}.pkl"
        with open(obj_file_path, "wb") as outp:
            pickle.dump(run, outp, pickle.HIGHEST_PROTOCOL)

        # PSMs Table
        psm_file_path = output_dir_path / f"psms_{output_prefix}.csv"
        run.psm_score_dataframe(
            file_name=psm_file_path, score_name="Amanda:AmandaScore"
        ).to_csv(psm_file_path)

    except FileNotFoundError:
        print("Output directory does not exist.")
        # displayan error message or exiting the program
        sys.exit(1)

    print("Outputting results done")

    # terminate the program
    sys.exit(0)


if __name__ == "__main__":
    print("in __name__ == __main__")

    main()
