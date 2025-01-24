FROM ubuntu:24.04
MAINTAINER JTB 
RUN apt-get update && apt-get install wget cmake gcc -y
RUN wget https://sourceforge.net/projects/arma/files/armadillo-10.6.0.tar.xz
RUN apt-get install gzip   xz-utils -y
RUN tar -xf armadillo-10.6.0.tar.xz
RUN apt-get install gcc libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev -y
RUN apt-get install build-essential -y
RUN apt-get install git libopenblas-dev liblapack-dev libboost-all-dev -y
RUN cd /armadillo-10.6.0 && ./configure && make && make install 
#https://xuchang0201.github.io/mtPGS/documentation/02_installation.html
RUN git clone https://github.com/xuchang0201/mtPGS.git 
RUN cd mtPGS/src && make
RUN apt-get install python3-pip python3-dev -y
RUN ln -s /usr/bin/python3 /usr/bin/python
#RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py
RUN pip3 install scipy h5py numpy --break-system-packages
RUN git clone https://github.com/getian107/PRScs.git  && cd PRScs
RUN git clone https://github.com/getian107/PRScsx.git  && cd PRScsx
ENV PATH="$PATH:/PRScs/"
ENV PATH="$PATH:/PRScsx/"
#
### Install PRSice
RUN git clone https://github.com/choishingwan/PRSice.git && cd PRSice && mkdir build && cd build && cmake ../ && make
RUN apt-get install r-base-dev -y
RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_github('borangao/GECKO')"
RUN R -e "install.packages('ggplot2')"
ENV PATH="$PATH:/mtPGS/src/"

