#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Metadata is imported as pandas DataFrame. There is only one such file for a
given experiment and must contain at least one row, and 2 mandatory columns:
'label', and 'period'.
"""
import pandas as pd
import warnings
import yaml


class MetadataError(Exception):
    pass


class MetadataMissingMainLabel(MetadataError):
    pass


class MissingEntry(MetadataError):
    pass


class MissingPeriod(MissingEntry):
    pass


def load_metadata(path_to_yaml_file):
    """Use yaml parser and call Metadata constructor"""
    stream = open(path_to_yaml_file, 'r')
    docs = list(yaml.load_all(stream))
    stream.close()
    md = Metadata(docs)
    return md


def _update_local_dict(local_dict, new_items):
    """Update dictionary local_dict with items from new_items
    
    Emits warning when an item is tried for overiding
    
    Parameters
    ----------
    local_dict : dict
    new_items : dict
    """
    for key, value in new_items.items():
        if key in local_dict.keys():
            warnings.warn('Multiple values for one key, last value is kept.\n'
                          ' Check input file to remove duplicate entries.')
            local_dict[key] = value


class Metadata(object):
    """Experiment metadata
    
    Parameters
    ----------
    ddict : iterable of dict
        each dict is a metadata entry corresponding to either top, experiment
        level ('level': 'experiment'), or to specific containers, in which
        case container label is mandatory (e.g. 'label': 'container_01');
        this iterable should be returned by the yaml.load_all() method
    
    Attributes
    ----------
    top : :class:`LocalMetadata` instance
        metadata corresponding to top, experiment level
    locals : dict
        key: container label, value is :class:`LocalMetadata` instance
    period : float (or int)
        minimal period found in top level (when multiple periods are saved)
    """
    def __init__(self, iter_dict):
        # parse iter_dict argument and build internal dict representation
        self._iter_dict = list(iter_dict)  # store input as a list of dicts
        self._ddict = {}  # dict: container label -> metadata dict

        # build internal metadata dicts

        self._setup_internals()
        # build LocalMetadata API
        self._setup_meta()
    
    def _setup_internals(self):
        self._ddict = {}
        for dic in self._iter_dict:
            if dic is None:
                continue
            if 'level' in dic and dic['level'] in ['experiment', 'top']:
                self._ddict['top'] = dic
            else:
                labs = dic['label']
                container_labs = []
                if isinstance(labs, list):
                    container_labs = labs
                else:  # labs is only one label
                    container_labs = [labs, ]
                for lab in container_labs:
                    if lab in self._ddict.keys():
                        _update_local_dict(self._ddict[lab], dic)
                    else:
                        self._ddict[lab] = dic
        # check that top/experiment level has been found
        if 'top' not in self._ddict.keys():
            raise MetadataMissingMainLabel('Missing "top"/"experiment" level')
    
    def _setup_meta(self):
        self.top = LocalMetadata(self._ddict['top'], self._ddict['top'])
        self.locals = {}
        for key, value in self._ddict.items():
            if key != 'top':
                self.locals[key] = LocalMetadata(value, self._ddict['top'])
    
    def from_container(self, container_label):
        """Get LocalMetadata instance corresponding to container label"""
        try:
            return self.locals[container_label]
        except KeyError:
            return self.top
    
    @property
    def loc(self):
        """Kept for legacy with pandas.DataFrame.loc. Erase?"""
        return self.locals
    
    @property
    def period(self):
        """Get top level (experiment) period
        
        Note
        -----
        There might be more than acquisition periods when more than 1
        channel are used; the smallest is taken as bare reference for now
        """
        reduced = [(k, v) for k, v in self._ddict['top'].items() if 'period' in k]
        sorted_ = sorted(reduced, key=lambda x: x[1], reverse=False)
        return sorted_[0][1]
    
    def __getitem__(self, key):
        """Use parameter key to retrieve metadata in top level"""
        if key == 'period':
            return self.period
        return self.top[key]
    
    def __str__(self):
        msg = str(self.top)
        return msg
    
    def __repr__(self):
        return str(self)
    
    def to_yaml(self, stream=None):
        """Exports metadata to yaml file"""
#        locs = [self.top, ] + [dic for k, dic in self.locals.items()]
#        iter_dict = []
#        for loc in locs:
#            # print(loc)
#            iter_dict.append(loc._dict)
        out = yaml.dump_all(self._iter_dict, stream=stream, default_flow_style=False)
        if stream is None:
            print(out)
    

class LocalMetadata(object):
    """Class to search for lower level metadata (containers)
    
    Parameters
    ----------
    meta : :class:`Metadata` instance
        the instance from which the local metadata is derived; items not found
        in local dict are returned from top level (experiment) dict (when they
        exist)"""

    def __init__(self, local_dict, top_dict):
        self._top = top_dict
        self._dict = local_dict
        keys = [key for key in self._top.keys()]
        for key in self._dict.keys():
            if key not in keys:
                keys.append(key)
        self._keys = keys

    def __getitem__(self, key):
        """Use instance as a dict; look first locally; if not found go top"""
        try:
            return self._dict[key]
        except KeyError:
            return self._top[key]
    
    @property
    def loc(self):
        """Legacy attribute"""
        return self
    
    @property
    def period(self):
        """Returns period"""
        reduced = [(k, v) for k, v in self._dict.items() if 'period' in k]
        if reduced:
            sorted_ = sorted(reduced, key=lambda x: x[1], reverse=False)
            return sorted_[0][1]
        else:
            return self._top['period']
    
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
        

#def load_from_csv(filename, sep=','):
#    """Gets pandas DataFrame from metadata file
#
#    Parameters
#    ----------
#    metadata_file : str
#        absolute path to metadata csv file
#    sep : str (default ',')
#        default separator for csv files
#
#    Returns
#    -------
#    pandas.DataFrame
#
#    Note
#    ----
#    There must be a column named 'label' that stores either the experiment
#    label, and/or container file labels.
#    """
#    possible_names = ['interval',
#                      'time interval', 'time_interval', 'time-interval',
#                      'dt',
#                      'delta time', 'delta_time', 'delta-time',
#                      'delta-t', 'delta_t', 'delta t',
#                      'period']
#    df = pd.read_csv(filename, sep=sep, index_col='label')
#    boo = False
#    for name in possible_names:
#        if name in df.columns:
#            boo = True
#            break
#    if not boo:
#        msg = 'Period acquisition not found in metadata!'
#        raise MissingPeriod(msg)
#    # rename period column in case other name is used
#    else:
#        df = df.rename(columns={name: 'period'})
#    return df
#
#
#def fill_rows(df, exp_label, cont_labels):
#    """Fill df rows for each container label.
#
#    NA values in such rows are filled with experiment row metadata.
#
#    Parameters
#    ----------
#    df : pandas.DataFrame
#        There must be at least one row, indexed with the experiment label; and
#        there must be at least one column that reports for the chosen
#        acquisition period.
#    exp_label : str
#        experiment label
#    cont_labels : list of str
#        list of container labels
#
#    Returns
#    -------
#    pandas.DataFrame
#    """
#    #  reindex with exp label and each container label as row indices
#    rows = [exp_label, ] + cont_labels
#    df = df.reindex(rows)  # initially non-reported containers are set NaNs
#    # fill na values with experiment values
#    meta = df.fillna(df.loc[exp_label])
#    del df
#    return meta
#
#
#def get_period(meta, label):
#    """Gets the period acquisition for corresponding label.
#
#    Parameters
#    ----------
#    meta : pandas.DataFrame
#        must have a 'period' column and label is a valid index
#    label : str
#
#    Returns
#    -------
#    float
#    """
#    return meta.loc[label, 'period']


if __name__ == '__main__':
    md = load_metadata('/home/joachim/tmptunacell/test_yaml.yml')
#    print(md)
    output_stream = open('/home/joachim/tmptunacell/myyaml.yml', 'w')
    print(md.to_yaml(output_stream))
    output_stream.close()
    