FROM rootproject/root-ubuntu16

MAINTAINER dcraik@cern.ch

USER root

# need kerberos to access eos
ADD krb5.conf /etc/krb5.conf
RUN apt-get -y update \
    && apt-get -y install krb5-user

RUN apt-get -y update \
    && apt-get -y install wget

# adding the analysis code to the container
WORKDIR /code
ADD . /code

# setup EvtGen
#RUN chmod +x setupEvtGen.sh
#RUN export ROOTSYS=/usr/local/ && printenv && ./setupEvtGen.sh

# build RapidSim
RUN export ROOTSYS=/usr/local/ && mkdir -p build && cd build && cmake ../ && make

# build RapidSim with EvtGen
#RUN export ROOTSYS=/usr/local/ && export LD_LIBRARY_PATH=/code/EvtGen/external/HepMC/lib:/code/EvtGen/external/pythia8186/lib:/code/EvtGen/external/PHOTOS/lib:/code/EvtGen/external/TAUOLA/lib:/code/EvtGen/evtgen/lib:$LD_LIBRARY_PATH && export PYTHIA8DATA=/code/EvtGen/external/pythia8186/xmldoc && export EVTGEN_ROOT=/code/EvtGen/evtgen/ && printenv && mkdir -p build-evtgen && cd build-evtgen && cmake ../ && make

RUN chmod +x run.sh
CMD ["./run.sh"]
