# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [4.0.2] - 2017-08-21
### Added
- load: max_size parameter added that filters out any families with family sizes larger than the value
- lambda: warning added when calculated posterior probability for a family is 0
- Use Intel compiler and MKL library for calculations if they are available

### Changed
- Fixes to allow compilation with gcc 4.9
- Errormodel command no longer fails if model has fewer lines than the expected maximum family size


## [4.0.1] - 2017-05-09
### Added
- Configure command to allow more precise tuning for a user's system
- report: "json" flag added to output data in that format
- load: Files with "CSV" extension will now be interpreted as comma-separated values
- C++11 support is now required for compilation

### Changed
- Buffer overflow error when setting random family sizes fixed


## [4.0.0] - 2017-03-15
### Added
- License file and this changelog. With this release a new website and Google group for questions are available.
- Command line history now available if compiled with the USE_READLINE option
- report: html flag added and the html format modified to be more useful
- A -version flag was added to the exe so version can be detected in one command

### Changed
- Manual was updated to reflect website and Google group
- Makefile was modified for easier compilation under Ubuntu
- lhtest: Fixed crashing bug if directory not specified
- genfamily: Fixed crashing bug if directory not in correct format

