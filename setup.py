#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
setuptools for tuna package
"""

from setuptools import setup, find_packages

setup(name='tuna',
      version='0.0.6',
      description='Time-lapse UNicellular Analyzer',
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
