#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Metadata is imported as pandas DataFrame. There is only one such file for a
given experiment and must contain at least one row, and 2 mandatory columns:
'label', and 'period'.
"""
import pandas as pd


class MetadataError(Exception):
    pass


class MissingEntry(MetadataError):
    pass


class MissingPeriod(MissingEntry):
    pass


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
