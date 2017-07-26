#!/usr/bin/env python2
# -*- coding: utf-8 -*-
descr = """Analysis of time-lapse data of growing, and dividing cells

Analysis of time-lapse data from microscopy of growing micro-organisms,
including tree reconstruction, time-series visualization, computation of
statistics of dynamic, and cell-cycle variables.
"""

from setuptools import setup, find_packages

setup(name='tunacell',
      version='0.0.6',
      description='Time-lapse of unicellular organisms analyzer',
      long_description=descr,
      url='',
      author='Joachim Rambeau',
      author_email='joachim.rambeau@gmail.com',
      license='MIT',
      #packages=['tuna'],
      packages=find_packages(),  
      install_requires=['numpy',
                        'pandas',
                        'treelib',
                        'matplotlib>=2',  # new color code...
                        'tables',
                        'future'],
#      scripts=['bin/commande.py'],
      zip_safe=False)
