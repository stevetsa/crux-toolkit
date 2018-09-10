#### 
# Original Repo - https://github.com/crux-toolkit/crux-toolkit
# http://crux.ms/
# The Crux mass spectrometry analysis toolkit is an open source project that aims to provide users with a cross-platform suite of 
# analysis tools for interpreting protein mass spectrometry data.
####

FROM ubuntu:16.04
MAINTAINER Steve Tsang <mylagimail2004@yahoo.com>
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 git-all \ 
 cmake \
 wget

RUN apt-get update -y && \
 apt-get upgrade -y && \
 apt-get dist-upgrade -y && \
 apt-get install build-essential software-properties-common -y && \
 add-apt-repository ppa:ubuntu-toolchain-r/test -y && \
 apt-get update -y && \
 apt-get install gcc-7 g++-7 -y && \
 update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7 && \
 update-alternatives --config gcc
RUN apt-get install -y subversion 

WORKDIR /opt/
RUN git clone https://github.com/stevetsa/crux-toolkit.git
WORKDIR /opt/crux-toolkit

#cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=<install directory> .
#where <install directory> is the location where you wish to install the crux programs. If you want crux to go in $HOME/bin, then <install directory> would be $HOME. The installation process will automatically put the programs in a directory called bin.
#Note, this configuration will build with optimizations turned on, and will not include debug symbols. To build with optimizations turned off, and debug symbols included, use the command:

RUN cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX:PATH=/usr/local .
RUN make
RUN make install

#RUN cp /opt/fastq_screen_v0.12.0/fastq_screen /usr/local/bin
#RUN cp /opt/bowtie2-2.3.4.1-linux-x86_64/bowtie* /usr/local/bin

