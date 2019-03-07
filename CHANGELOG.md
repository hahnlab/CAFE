# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [4.2.1] - 2019-03-15
### Changed
- Fixed a memory leak in situations where branch lengths were almost identical to each other
- Improved error messages when there is an error in the family file
- Fixed bug where the tree was printed twice in the report

## [4.2] - 2018-07-01
### Added
- tree: a -i parameter is now supported to load a tree from a file.
- A wider variety of compiler/operating system combinations are now tested before releasing
- lambda: a -score parameter has been added which logs the calculated 
log-likelihood for the given set of lambdas.
- The PGI compiler will now be autodetected and used if available.

### Changed
- Fixed bug where running the optimizer could change the initial values passed to the optimizer.
(This did not seem to have an effect on the final values but made the optimizer less predictable)
- load: Fixed occasional crash if the log file specified with the -l parameter is invalid 
- report: Fixed occasional crash if multiple threads were specified in the load command
- lambda: Fixed crash if there is a species in the family file that does not appear in the tree 
- tree: A warning is now emitted if the tree is not ultrametric (to within 0.01% of the max 
root-to-leaf length)

## [4.1] - 2017-10-27
### Added
- seed: new command that sets the random seed so commands with randomness can be replicated

### Changed
- report: fixed bug where viterbi values were showing incorrect data
- load: fixed bug in handling trees with more than about 100 species

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

