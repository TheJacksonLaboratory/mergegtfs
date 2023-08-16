#!/usr/bin/env python

"""
Contains initialization and argument handling code. 
"""

import argparse
import os
import re
import socket
import sys
import textwrap
import time

## local:

import intypes
import util

version = "20230816a "

description = (
    "Merges redundant transcripts from multiple GTF2.2 formatted "
    "files listed in gtf_list_file, resulting in a non-redundant union gtf."
)

epilog = (
    "Input GTF2.2 formatted files are required to have exon features "
    "with attributes that include 'gene_id' and 'transcript_id'. Output "
    "GTF2.2 formatted file includes 'exon' and 'transcript' features, both "
    "with attributes (in order) 'gene_id' and 'transcript_id'. New gene_ids "
    "will be assigned in accordance with --p_exon_overlap (min overlap for "
    "matching exons) and --p_exons_overlap (min proportion of matched exons "
    "for matching gene_ids). New transcript_ids are assigned sequentially "
    "for each gene. New gene_ids and transcript_ids both have prefix "
    "specified by --gene_prefix. It is assumed that the order of exons in "
    "the input GTF file reflects the actual order in the spliced transcript; "
    "it is also assumed that for negative strand transcripts, the last exon "
    "is listed first, and the first exon is listed last. If the order of "
    "exon records for negative strand transcripts in the input GTF file "
    "begins with the last exon and ends with the first, make sure to set "
    "'--rev_neg_exons'. If the order of exons in the input file is "
    "indeterminate and does not reflect exon order in the transcript, please "
    "set '--sort_exons'; however, this may collapse isoforms where exon "
    "order really is shuffled due to e.g. local genomic rearrangement."
)

def build_parser():

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = textwrap.dedent(f'''\
        mergegtfs.py version {version}
        --------------------------------
        {description}
        '''),
        epilog = epilog
    )

    parser.add_argument(
	    "gtf_list_file", 
	    help="Tab-delimited file with 1 unique label followed by 1 GTF2.2 file path per line"
    )

    parser.add_argument(
	    "--tol_sj", 
	    type=intypes.non_negative_int, 
	    default=0,
	    help="Tolerance (bp) for matching splice junction coordinates"
    )

    parser.add_argument(
        "--tol_tss", 
        type=intypes.non_negative_int, 
        default=0,
        help="Tolerance (bp) for matching transcript start coordinates"
    )

    parser.add_argument(
	    "--tol_tts", 
	    type=intypes.non_negative_int, 
	    default=0,
	    help="Tolerance (bp) for matching transcript end coordinates"
    )

    parser.add_argument(
	    "--p_exon_overlap",
	    type=intypes.float_proportion,
	    default=0.5,
	    help="Minimum proportion overlap between two exons needed for gene matching"
    )

    parser.add_argument(
	    "--p_exons_overlap",
	    type=intypes.float_proportion,
	    default=0.1,
	    help="Minimum proportion of overlapping exons needed for gene matching"
    )

    parser.add_argument(
        "--max_intron_length",
        type=intypes.non_negative_int,
        default=1250000,
        help="Maximum expected intron size; used for identifying fusions"
    )

    parser.add_argument(
        "--sort_exons",
        action="store_true",
        help="For each transcript, sort exons in ascending order by (begin, end)"
    )

    parser.add_argument(
        "--rev_neg_exons",
        action="store_true",
        help="For each transcript, reverse order of negative strand exons"
    )

    parser.add_argument(
	    "--gene_prefix",
	    type=intypes.non_whitespace_str,
	    default='LOC.',
	    help="Prefix for gene_ids and transcript_ids"
    )

    parser.add_argument(
	    "--output_prefix",
	    type=intypes.non_whitespace_str,
	    default='union',
	    help="Prefix for output file names"
    )

    return parser


def initialize():
 
    time_start = time.time()

    parser = build_parser()
    params = parser.parse_args()
    params.time_start = time_start

    print(
        f"version: {version}\n"
        f"begin: {util.time_stamp()}\n"
        f"hostname: {socket.gethostname()}\n"
        f"working_directory: {os.getcwd()}"
    )

    arg_names = [
        'gtf_list_file',
        'tol_sj',
        'tol_tss',
        'tol_tts',
        'p_exon_overlap',
        'p_exons_overlap',
        'rev_neg_exons',
        'sort_exons',
        'max_intron_length',
        'output_prefix',
    ]

    for arg_name in arg_names:
        print(f"{arg_name}: {getattr(params, arg_name)}")

    return params


if __name__ == "__main__":
    params = initialize()
    print(f"parse_args() yielded: {params}")
    sys.exit(0)

