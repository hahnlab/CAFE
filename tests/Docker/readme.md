# Testing in Various Environments

Dockerfiles for various environments to be supported by CAFE

## Getting Started

The Dockerfiles represent the supported compilation platforms officially supported by CAFE: 
* Ubuntu 14.04 (trusty) with Clang 3.4
* Ubuntu 16.04 (xenial) with PGI Community Edition 17.10 
* Ubuntu 17.10 (artful) with GCC 7.2,
* Debian 8 (jessie) with gcc 5.4.


### Running the tests

Build the dockerfiles. Make sure the cafe_test.sh script is available when the image is created.

For each image, create an instance and run the cafe_test.sh script. The latest CAFE code should be 
downloaded, compiled in the environment, and tested.

## Authors

* **Ben Fulton** - *Initial work* - (https://github.com/benfulton)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

Adapted by Indiana University from Educational Community License, Version 2.0 (ECL-2.0), available at http://opensource.org/licenses/ECL-2.0.

See the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

