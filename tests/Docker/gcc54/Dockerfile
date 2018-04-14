# Set the base image
FROM gcc:5.4

# Dockerfile author / maintainer 
MAINTAINER Name befulton@iu.edu

RUN apt-get update && apt-get -y upgrade && apt-get -y install git automake cpputest

COPY cafe_test.sh .

