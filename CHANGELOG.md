# Changelog
All notables change to this project will be documented in this file

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.1] - 2019-07-03
### Added
- __main__ module and its entry point to print out tunacell version `tunacell -V`, or `python -m tunacell -V`
  which simplifies the cumbersome `python -c "from tunacell import __version__; print(__version__)"`

### Fixed
- Parser.info_samples() when no sample in Parser instance
- hidden files under the `containers` folder are skipped (not retained as a container)

### Changed
- bin/tunasimu script is moved as a module: tunacell/simu/run.py and its new entry point is properly declared in setup.py
  (call to executable `tunasimu` still valid)

### Removed
- Parser.iter_cells, Parser.iter_colonies, Parser.iter_containers methods (were moved to Experiment instances)

## [0.2.0] - 2019-07-01
### Added
- supersegger module to read directly from supersegger output
- Travis CI support
- CHANGELOG
- documentation source, hosted by readthedocs at [tunacell.readthedocs.io](https://tunacell.readthedocs.io/en/latest/)

### Changed
- Use analysis module instead of everything in text module
- Metadata yaml and csv syntax
- cell identifiers are kept to their original type (not translated to str)

### Fixed
- Root cells can point to non-zero label
  
## [0.1.0] - 2018-01-29
### Added
- tunacell API
- 10-minute tutorial script
- analysis example scripts

### Deprecated
- Parser.iter_cells, Parser.iter_colonies, Parser.iter_containers methods are now attached to Experiment instances
