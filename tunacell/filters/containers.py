#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines filters for Container instances.
"""
from __future__ import print_function

from tunacell.filters.main import FilterGeneral


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
    KeyError : when requested key is not present in container metadata
    """

    def __init__(self, key, value):
        self.key = key
        self.value = value
        label = '{}={}'.format(key, value)
        self.label = label
        return

    def func(self, container):
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
