FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git \
	cmake \
	sudo \
	python3-pip \
	wget \
	unzip \
	vim \
	openjdk-8-jdk \
	ant \
	ca-certificates-java && \
	apt-get clean

RUN	pip3 install configobj \
	numpy \
	scipy \
	biopython \
	pytest

# Fix certificate issues
RUN update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

RUN mkdir -p tral_repository
WORKDIR /tral_repository 
RUN git clone https://github.com/acg-team/tral.git

WORKDIR /tral_repository/tral/easy_setup
RUN rm configTRAL_path.cfg
COPY ./configTRAL_path.cfg /tral_repository/tral/easy_setup

COPY ./.tral /

RUN chmod +x *.sh && \
	sudo ./setupTRAL.sh setup
# RUN sudo ./install_ext_software.sh # install all ext. software in one layer
RUN ./install_ext_software/alf.sh
RUN ./install_ext_software/castor.sh
RUN ./install_ext_software/hhrepid.sh
RUN ./install_ext_software/hmmer.sh
RUN ./install_ext_software/mafft.sh
RUN ./install_ext_software/phobos.sh
RUN ./install_ext_software/phyml.sh
RUN ./install_ext_software/treks.sh
RUN ./install_ext_software/trf.sh
# RUN ./install_ext_software/trust.sh # Currently not supported
RUN ./install_ext_software/xstream.sh

WORKDIR /tral_repository/tral

### HINTS usage:
# create Github.com Token: https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line
# $ cat ~/TOKEN.txt | docker login https://docker.pkg.github.com -u USERNAME --password-stdin
# sudo docker docker pull docker.pkg.github.com/acg-team/tral/tral_docker:latest
# sudo docker run -ti docker.pkg.github.com/acg-team/tral/tral_docker

### HINTS development:
# Read: https://help.github.com/en/packages/publishing-and-managing-packages/about-github-packages
# cat ~/TOKEN.txt | docker login https://docker.pkg.github.com -u USERNAME --password-stdin
# cd ~/tral_repository/tral/docker/
# sudo docker build -t docker.pkg.github.com/acg-team/tral/tral_docker:latest -f ~/tral_repository/tral/docker/Dockerfile --no-cache 
# sudo docker push docker.pkg.github.com/acg-team/tral/tral_docker:latest
