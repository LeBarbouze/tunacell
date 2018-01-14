#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Metadata can be used to filter containers, or to use the 'period' parameter
when computing local_rates by joining mother and daughter cells (specifically,
the joining timepoint is set at the first registered timepoint in current
daughter cell, minus half the acquisition period).

Metadata can be loaded from two different types of file:
    * csv, which can be exported from spreadsheet software, but is hardly
      readable from a raw text file;
    * yaml, which is both readable and easy to fill in as a raw text file.
    
One entry must be given for the constructor to work properly: the ``level``
entry that allows to distinguish 'top' experiment-level metadata, to
lower-level container metadata.

Minimal example for a csv file:
    
    level,label,author,strain
    top,my_experiment,J. Rambeau,e.coli
    container_12,weird_container,,weirdo.weirdus

and the same for yaml file:
    
    level     : top
    label     : my_experiment
    author    : J. Rambeau
    strain    : e.coli
    # use 3 dashes to separate levels
    ---
    level     : container_12
    label     : weird_container
    strain    : weirdo.weirdus
    
These metadata files indicate that the author is the same for all containers,
and that e.coli is used in all containers but container_12 which has been
labeled 'weird_container' and contains the new species 'weirdo.weirdus'.

A hidden file will be written when data is parsed for the first time, and
will report for a certain number of 
"""
import warnings
import yaml
import csv
import os


class MetadataError(Exception):
    pass


class MetadataNotFound(MetadataError):
    pass


class MetadataMissingMainLabel(MetadataError):
    pass


class MissingEntry(MetadataError):
    pass


class MissingPeriod(MissingEntry):
    pass


def load_metadata(experiment_path):
    """Loads metadata using experiment absolute path
    
    Parameters
    ----------
    experiment_path : str
        absolute path to experiment root directory
    
    Returns
    -------
    :class:`Metadata` instance
        metadata found in experiment root folder
    
    Raises
    ------
    MetadataNotFound
        when file is not found
    """
    ls = os.listdir(experiment_path)
    path = None
    ext = None
    for bn in ls:
        fn = os.path.join(experiment_path, bn)
        if os.path.isfile(fn):
            if 'metadata' in os.path.basename(fn):
                path = fn
                _, ext = os.path.splitext(fn)
                break
    if path is None:
        raise MetadataNotFound()
    if ext in ['.yml', '.yaml']:
        md = load_from_yaml(path)
    elif ext == '.csv':
        md = load_from_csv(path, sep=',')
    elif ext == '.tsv':
        md = load_from_csv(path, sep='\t')
    return md


def load_from_yaml(path_to_yaml_file):
    """Builds metadata from yaml files
    
    Parameters
    ----------
    path_to_yaml_file : str
        absolute path to yaml metadata file
    
    Returns
    -------
    :class:`Metadata` instance
    """
    with open(path_to_yaml_file, 'r') as f:
        docs = list(yaml.load_all(f))
    md = Metadata(docs)
    return md


def load_from_csv(filename, sep=','):
    """Builds metadata from csv files

    Parameters
    ----------
    metadata_file : str
        absolute path to metadata csv file
    sep : str (default ',')
        default separator for csv files

    Returns
    -------
    :class:`Metadata` instance

    Note
    ----
    There must be a column named 'label' that stores either the experiment
    label, and/or container file labels.
    """
    list_dict = []
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            list_dict.append({k: v for k, v in row.items() if v})
    md = Metadata(list_dict)
    return md


def _update_local_dict(local_dict, new_items):
    """Update dictionary local_dict with items from new_items
    
    Emits warning when an item is tried for overiding
    
    Parameters
    ----------
    local_dict : dict
    new_items : dict
        (must have an .items() method)
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
                labs = dic['level']
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
            raise MetadataMissingMainLabel('Missing "top" experiment level')
    
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
        s = yaml.dump_all(self._iter_dict, default_flow_style=False)
        return s
#        msg = str(self.top)
#        return msg
    
    def __repr__(self):
        return str(self)
    
    def to_yaml(self, stream=None):
        """Exports metadata to yaml file"""
        out = yaml.dump_all(self._iter_dict, stream=stream, default_flow_style=False)
        if stream is None:
            print(out)
        
    def to_csv(self, stream=None):
        """Exports metadata to csv file"""
        # TODO: write this method
        pass
    

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
        """String output based on yaml.dump"""
        s = yaml.dump(self._dict, default_flow_style=False)
        return s

    def __repr__(self):
        return str(self)


def load_counts(experiment_abspath):
    """Reads hidden file to load dict of counts (for any given filterset)
    
    Parameters
    ----------
    experiment_abspath : str
    
    Returns
    -------
    dict of dicts
        primary key is the representation of a FilterSet, secondary key is one
        from 'cells', 'colonies', 'lineages', 'containers'
    """
    pass



if __name__ == '__main__':
    md = load_from_yaml('/home/joachim/tmptunacell/test_yaml.yml')
#    print(md)
    output_stream = open('/home/joachim/tmptunacell/myyaml.yml', 'w')
    print(md.to_yaml(output_stream))
    output_stream.close()
    