#!/bin/bash

set -e 

git clone https://github.com/hahnlab/CAFE.git

cd CAFE

autoconf
./configure
make prep
make test
test/runtests

