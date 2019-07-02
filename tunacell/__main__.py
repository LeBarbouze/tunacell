#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse


from tunacell import __version__


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-V', '--version', action='store_true', help='Print tunacell version')
    args = argparser.parse_args()

    if args.version:
        print('tunacell version {}'.format(__version__))