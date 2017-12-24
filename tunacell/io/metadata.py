#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Metadata is imported as pandas DataFrame. There is only one such file for a
given experiment and must contain at least one row, and 2 mandatory columns:
'label', and 'period'.
"""
import pandas as pd
import warnings


class MetadataError(Exception):
    pass


class MetadataMissingMainLabel(MetadataError):
    pass


class MissingEntry(MetadataError):
    pass


class MissingPeriod(MissingEntry):
    pass


class Metadata(object):
    """Experiment metadata
    
    Parameters
    ----------
    ddict : dict of dict
        keys are experiment and container labels, values are metadata content
        as dict
    exp : :class:`tunacell.base.experiment.Experiment` instance
    """
    
    def __init__(self, ddict, exp_label):
        if not isinstance(ddict, dict):
            raise TypeError('Input file must be a dictionary')
        # check experiment label
        if exp_label not in ddict.keys():
            raise MetadataMissingMainLabel('Missing experiment label')
        self._ddict = ddict
        # define experiment, top level metadata
        self.top_level = {}
        self._locals = {}
        for key, dic in self._ddict.items():
            if 'level' in dic and dic['level'] == 'experiment':
                self.top_level = dic
            self._locals[key] = LocalMetadata(self, key)
    
    @property
    def loc(self):
        """Kept for legacy with pandas.DataFrame.loc. Erase?"""
        return self._locals
    
    @property
    def period(self):
        return self.top_level['period']
    
    def __getitem__(self, key):
        try:
            return self.top_level[key]
        except KeyError:
            warnings.warn('Metadata key error {}, value is assigned to None'.format(key))
            return None
    
    def __str__(self):
        msg = str(self._locals[self.top_level['label']])
        return msg
    
    def __repr__(self):
        return str(self)
    

class LocalMetadata(object):
    """Container metadata"""

    def __init__(self, meta, label):
        self._top = meta
        self._dict = meta._ddict[label]
        keys = [key for key in meta.top_level.keys()]
        for key in self._dict.keys():
            if key not in keys:
                keys.append(key)
        self._keys = keys

    def __getitem__(self, key):
        try:
            return self._dict[key]
        except KeyError:
            return self._top[key]
    
    @property
    def loc(self):
        return self
    
    def __str__(self):
        msg = '{:<20s} {:<20s}\n'.format('Parameter', 'Value')
        for key in self._keys:
            if key == 'level':
                continue
            val = '{}'.format(self[key])
            msg += '{:<20s} {:<20s}\n'.format(key, val)
        return msg.strip()
    
    def __repr__(self):
        return str(self)
        

def load_from_csv(filename, sep=','):
    """Gets pandas DataFrame from metadata file

    Parameters
    ----------
    metadata_file : str
        absolute path to metadata csv file
    sep : str (default ',')
        default separator for csv files

    Returns
    -------
    pandas.DataFrame

    Note
    ----
    There must be a column named 'label' that stores either the experiment
    label, and/or container file labels.
    """
    possible_names = ['interval',
                      'time interval', 'time_interval', 'time-interval',
                      'dt',
                      'delta time', 'delta_time', 'delta-time',
                      'delta-t', 'delta_t', 'delta t',
                      'period']
    df = pd.read_csv(filename, sep=sep, index_col='label')
    boo = False
    for name in possible_names:
        if name in df.columns:
            boo = True
            break
    if not boo:
        msg = 'Period acquisition not found in metadata!'
        raise MissingPeriod(msg)
    # rename period column in case other name is used
    else:
        df = df.rename(columns={name: 'period'})
    return df


def fill_rows(df, exp_label, cont_labels):
    """Fill df rows for each container label.

    NA values in such rows are filled with experiment row metadata.

    Parameters
    ----------
    df : pandas.DataFrame
        There must be at least one row, indexed with the experiment label; and
        there must be at least one column that reports for the chosen
        acquisition period.
    exp_label : str
        experiment label
    cont_labels : list of str
        list of container labels

    Returns
    -------
    pandas.DataFrame
    """
    #  reindex with exp label and each container label as row indices
    rows = [exp_label, ] + cont_labels
    df = df.reindex(rows)  # initially non-reported containers are set NaNs
    # fill na values with experiment values
    meta = df.fillna(df.loc[exp_label])
    del df
    return meta


def get_period(meta, label):
    """Gets the period acquisition for corresponding label.

    Parameters
    ----------
    meta : pandas.DataFrame
        must have a 'period' column and label is a valid index
    label : str

    Returns
    -------
    float
    """
    return meta.loc[label, 'period']


if __name__ == '__main__':
    md = Metadata({'simutest': {'level': 'experiment', 
                                'label': 'simutest',
                                'period': 12.,
                                'medium': 'M9'},
                   'container_01': {'label': 'container_01',
                                    'species': 'e.coli'},
                   'container_02': {'label': 'container_02',
                                    'species': 'b.subtilis',
                                    'medium': 'RDM'}
                                    }, 'simutest')