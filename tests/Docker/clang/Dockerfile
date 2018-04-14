
FROM        ubuntu:trusty
MAINTAINER  Robin Sommer <robin@icir.org>

# Setup environment.
ENV PATH /opt/llvm/bin:$PATH

# Default command on startup.
CMD bash

# Setup packages.
RUN apt-get update && apt-get -y install clang git build-essential automake cpputest

ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++

COPY cafe_test.sh .

