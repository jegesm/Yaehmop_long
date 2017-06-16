#FROM ubuntu:16.04
FROM debian:8.8

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y mc \
                       wget \
                       vim \
                       git \
                       gfortran \
                       libblas-dev \
                       liblapack-dev \
                       libf2c2-dev \
                       make \
                       ssh 

ADD tightbind /opt/tightbind
WORKDIR /opt/tightbind/
RUN make 
RUN cp eht_parms.dat /usr/local/lib/eht_parms.dat
                   

WORKDIR /opt

