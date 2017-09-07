#!/usr/bin/env python2
# -*- coding: utf-8 -*-
descr = """Analysis of Timeseries from dividing UNicellular microorganisms

Analysis of time-lapse data from microscopy of growing micro-organisms,
including tree reconstruction, time-series visualization, computation of
statistics of dynamic, and cell-cycle variables.
"""

from setuptools import setup, find_packages

setup(name='tunacell',
      version='0.0.7',
      description='Analysis of Timeseries from dividing UNicellular microorganisms',
      long_description=descr,
      url='',
      author='Joachim Rambeau',
      author_email='joachim.rambeau@gmail.com',
      license='MIT',
      classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        #'Development Status :: 1 - Planning',
        #'Development Status :: 2 - Pre-Alpha',
        'Development Status :: 3 - Alpha',
        #'Development Status :: 4 - Beta',
        #'Development Status :: 5 - Production/Stable',
        #'Development Status :: 6 - Mature',
        #'Development Status :: 7 - Inactive',
        "License :: OSI Approved :: MIT License"
      ],
      #packages=['tuna'],
      packages=find_packages(),
      install_requires=['numpy',
                        'scipy',
                        'pandas',
                        'treelib',
                        'matplotlib>=2',  # new color code...
                        'future',
                        'dill'],
#      scripts=['bin/commande.py'],
      zip_safe=False)
