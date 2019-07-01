#!/usr/bin/env python2
# -*- coding: utf-8 -*-
descr = """Analysis of Timeseries from dividing UNicellular microorganisms

Analysis of time-lapse data from microscopy of growing micro-organisms,
including tree reconstruction, time-series visualization, computation of
statistics of dynamic, and cell-cycle variables.
"""
import codecs
import os
import re

from setuptools import setup, find_packages

# pieces below are taken from pip setup.py file to detect version
here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    # intentionally *not* adding an encoding option to open
    # see here: https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    return codecs.open(os.path.join(here, *parts), 'r').read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")
## end of pip setup extract


setup(name='tunacell',
      version=find_version('tunacell','__init__.py'),
      description='Analysis of Timeseries from dividing UNicellular microorganisms',
      long_description=descr,
      url='https://github.com/LeBarbouze/tunacell',
      author='Joachim Rambeau',
      author_email='joachim.rambeau@gmail.com',
      license='MIT',
      classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        #'Development Status :: 1 - Planning',
        #'Development Status :: 2 - Pre-Alpha',
        #'Development Status :: 3 - Alpha',
        'Development Status :: 4 - Beta',
        #'Development Status :: 5 - Production/Stable',
        #'Development Status :: 6 - Mature',
        #'Development Status :: 7 - Inactive',
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
      ],
      keywords='microscopy time-lapse bacteria cell-division analysis statistics dynamics cross-correlations',
      #packages=['tuna'],
      packages=find_packages(),
      install_requires=['numpy',
                        'scipy',
                        'pandas',
                        'treelib',
                        'matplotlib>=2',  # new color code...
                        'future',
                        'dill',  # serialization
                        'PyYAML',  # yaml parser
                        'tqdm',  # status bar
                        'tabulate',  # tables
#                        'pathlib2',  # backward compatibility for python2.7
                       ],
      scripts=['bin/tunasimu'],
      zip_safe=False)
