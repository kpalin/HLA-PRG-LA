#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Make known reference file 

Created at Tue Apr 11 09:56:33 2017 by Kimmo Palin <kpalin@merit.ltdk.helsinki.fi>
"""


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Make known reference for HLA-PRG-LA")

    parser.add_argument(
        "-f",
        "--fai",
        help="Fai index for reference to use (assume GRCh37) [default:%(default)s]",
        required=True)

    parser.add_argument(
        "-V",
        "--verbose",
        default=False,
        const=True,
        nargs="?",
        help="Be more verbose with output [and log to a file] [default:%(default)s]"
    )

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')
        if args.verbose != True:
            log_file_handler = logging.FileHandler(args.verbose)
            log_file_handler.setFormatter(
                logging.getLogger().handlers[0].formatter)
            logging.getLogger().addHandler(log_file_handler)

    return args


def process_fai_to_known(fname):
    import pandas as pd
    fai = pd.read_table(
        fname,
        delim_whitespace=True,
        header=None,
        names=[
            "contigID", "contigLength", "contigStart", "linelen", "steplen"
        ])
    known = pd.DataFrame([fai.contigID, fai.contigLength]).T
    known["ExtractCompleteContig"] = 0
    known = known.append(pd.DataFrame([["*", 0, 1]], columns=known.columns))
    known["PartialExtraction_Start"] = None
    known["PartialExtraction_Stop"] = None
    MHC_chr = known.contigID.str.contains("^(chr)?6$")
    known.loc[MHC_chr, "PartialExtraction_Start"] = 28477796
    known.loc[MHC_chr, "PartialExtraction_Stop"] = 33448355

    known.loc[known.contigID.str.startswith("gi|"),
              "ExtractCompleteContig"] = 0
    import sys
    known.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == '__main__':
    args = main()
    process_fai_to_known(args.fai)
