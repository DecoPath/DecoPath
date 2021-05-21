FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

MAINTAINER Sarah Mubeen "sarah.mubeen@scai.fraunhofer.de"

RUN apt-get update
RUN apt-get -y upgrade

# Install essentials
RUN apt-get install -y --no-install-recommends build-essential software-properties-common libffi6 libffi-dev libmagic-dev

# Install Python 3.7
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update
RUN apt-get install -y --no-install-recommends python3.7 python3.7-dev
RUN apt-get install -y --no-install-recommends python3-pip python3-setuptools libcurl4-openssl-dev libssl-dev libxml2-dev

# Install R 4.0
RUN apt-get install -y --no-install-recommends ca-certificates software-properties-common gnupg2 gnupg1
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80/ --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt-get update
RUN apt-get install -y r-base

# Upgrade pip
RUN python3.7 -m pip install --upgrade pip

# R configuration
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('DESeq2')"

# Copy python requirements
COPY ./requirements.txt /opt/decopath/
WORKDIR /opt/decopath

# Install python requirements
RUN python3.7 -m pip install -r requirements.txt

# Add --user python modules to PATH
ENV PATH="/home/decopath/.local/bin:$PATH"

# Website Port
EXPOSE 8000

# Setup File Structure
COPY . /opt/decopath