#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tuna package
============

filters/cells.py module
~~~~~~~~~~~~~~~~~~~~~~

Class definitions to filter cells.
"""
from __future__ import print_function

import copy
import numpy as np

from tuna.filters.main import FilterGeneral
from tuna.base.container import Container


class FilterContainer(FilterGeneral):
    """General class to filter containers"""

    _type = 'CONTAINER'


class FilterContainerAny(FilterContainer):
    """True for any container"""

    def __init__(self):
        label = 'True for any container'
        self.label = label
        return

    def func(self, container):
        return True


class FilterContainerMetadataEquals(FilterContainer):
    """Test a given metadata value

    Parameters
    ----------
    key : str
    value : requested type (str, float, int)

    Raises
    ------
    TypeError : if test called upon non-Container object
    KeyError : when requested key is not present in container metadata
    """

    def __init__(self, key, value):
        self.key = key
        self.value = value
        label = '{}={}'.format(key, value)
        self.label = label
        return

    def func(self, container):
        if not isinstance(container, Container):
            raise TypeError('argument not a container')
        if not hasattr(container.metadata, self.key):
            raise KeyError
        svalue = getattr(container.metadata, self.key)
        if type(self.value) == str:
            if self.value == svalue:
                return True
            else:
                return False
        elif type(self.value) == int:
            if int(svalue) == self.value:
                return True
            else:
                return False
        elif type(self.value) == float:
            if float(svalue) == self.value:
                return True
            else:
                return False
        else:
            raise TypeError('type not understood')
