#!/usr/bin/env python

"""
Types for parameter checking. 
"""

import argparse
import re

def float_proportion(val):

    try:
        val = float(val)
        if not (0 <= val <= 1):
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' not a proportion in the interval [0, 1]"
        )

    return val


def strictly_positive_int(val):

    try:
        val = int(val)
        if val < 1:
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' not a strictly positive integer"
        )    

    return val


def non_negative_int(val):

    try:
        val = int(val)
        if val < 0:
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' not a non-negative integer"
        )

    return val


def non_whitespace_str(val):

    pat = re.compile(r'\s')

    try:
        val = str(val)
        if pat.search(val):
            raise Exception()
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"val '{val}' should be a string with no whitespace"
        )

    return val
    

