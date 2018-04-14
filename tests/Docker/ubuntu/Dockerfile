# Set the base image
FROM ubuntu:17.10

# Dockerfile author / maintainer 
MAINTAINER Name befulton@iu.edu

RUN apt-get update && apt-get -y upgrade && apt-get -y install git build-essential automake cpputest

COPY cafe_test.sh .

