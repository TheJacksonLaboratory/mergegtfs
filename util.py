#!/usr/bin/env python

import time

def elapsed(params, places=3):
    return round(time.time() - params.time_start, places)


def time_stamp():
    """
    returns human readable, sortable timestamp str
    """

    time_stamp = time.strftime(
        '%Y%m%d%H%M%S',
        time.localtime(time.time())
    )

    return time_stamp
