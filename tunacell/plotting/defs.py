#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tunacell package
============

plotting/defs.py module
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

DEFAULT_COLORS = (
    "red",
    "blue",
    "purple",
    "green",
    "yellowgreen",
    "cyan",
    "magenta",
    "indigo",
    "darkorange",
    "pink",
    "yellow",
)
colors = DEFAULT_COLORS

# plotting parameters
params = {
    "length": {"bottom": 1.0, "top": 8.0, "delta": 2.0, "unit": "$\mu$m"},
    "dot_length": {"bottom": 1e-2, "top": 1e-1, "delta": 3e-2, "unit": "$\mu$m/hr"},
    "dotlog_length": {"bottom": 0.5, "top": 2.5, "delta": 0.5, "unit": "dbs/hr"},
    "width": {"bottom": 0.5, "top": 1.5, "delta": 0.2, "unit": "$\mu$m"},
    "fluo": {"bottom": 1e5, "top": 2e6, "delta": 5e5, "unit": "A.U."},
    "dot_fluo": {"bottom": 1e2, "top": 5e4, "delta": 1e4, "unit": "A.U./hr"},
    "dotlog_fluo": {"bottom": 0.1, "top": 3, "delta": 0.5, "unit": "dbs/hr"},
    "concentration": {"bottom": 2e5, "top": 5e5, "delta": 1e5},
    "volume": {"bottom": 0.0, "top": 4.0, "delta": 1.0, "unit": "$\mu$m$^3$"},
    "area": {"bottom": 1.0, "top": 8.0, "delta": 2.0, "unit": "$\mu$m$^2$"},
    "dotlog_area": {"bottom": 0.5, "top": 2.5, "delta": 0.5, "unit": "dbs/hr"},
    "density": {"bottom": 1e5, "top": 4e5, "delta": 1e5},
    "ALratio": {"bottom": 0.1, "top": 1.5, "delta": 0.4, "unit": "$\mu$m"},
    "age": {"bottom": 0.0, "top": 1.0},
}


def get_params(obs, params, *keys):
    return [params[obs][k] for k in keys]
